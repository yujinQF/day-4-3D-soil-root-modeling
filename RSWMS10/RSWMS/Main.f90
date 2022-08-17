!> \file Main.f90
!! \brief Main program for RSWMS

!> Program RSWMS_MAIN
Program RSWMS_MAIN
  Use iso_c_binding
  Use typedef
  Use ParamData, Only: lvtk,lOutPartrace,iter_tot,i_noConv
  Use tmctrl, Only: tmax,kout
  Use SolData, Only: theta,conc,theta_old,ConO,Kode,con,hold,cap,hTemp,hnew
  Use CumData, Only: WatIn, SolIn,Q,WBalR
  Use DoussanMat, Only: savelast,nplant
  Use RootData, Only: lDou
  Use sparsematrix
  Use RhizoData, Only: thetaNonEq, thetaNonEqOld,thetaTot, rhizoModel, lRhizo
  Use GridData, Only: nPt
  Use Output, Only: OutFEM, OutVTK, FlxOut, PartraceOut, OutDou, OutDouVTK

#ifdef WITH_PARTRACE
  Use GridData, Only: continu,nx,ny,nz,nPt,nElm,sink_cube
  Use SolData, Only: Vx,Vy,Vz,Vx_old,Vy_old,Vz_old,SoilSoluteConcentration,SoilSoluteMass,SoilSoluteMass_new,Uptake,Cum_uptake
  Use TempData, Only: tstart
  Use ParamData, Only: lPartUp,lSalinity
  Use ParticlesInRoot
  Use SoluteMod, Only: Veloc
  Use SolutePartrace
  Use Environmental, Only: DegradSoil,CalcFieldCapacity,CalcDepthFactor

#endif
  Use RSWMSMod, Only: RSWMS, RSWMS_INI, timestep
  Use MPIUtils
  
  Implicit None

  Integer(ap) :: jj,ipl
  Real(dp) ::t, dt
  Logical :: master=.true.

#ifdef WITH_PARTRACE
  Integer(ap):: icount, periodic=0
  Real(dp) :: timePt,dtPt,dummy_soil, dummy_root, dummy_diff
  Real(dp), Pointer, Dimension(:) :: sinkE => null()
  !Real(dp), Pointer, Dimension(:) :: fac_deg => null()
  Integer :: isavelast
  ipl=1
  Call MPIInit()
  Call MPIRank(myrank)
  mpi=.True.
  master=(myrank==0)
#endif

  If(master) Call RSWMS_INI(t,dt)

#ifdef WITH_PARTRACE
  Call BroadcastReal(t)
  Call BroadcastReal(dt)
  Call BroadcastReal(tMax)
  Call BroadcastReal(tStart)
  Call initPartrace()
  If(master) Then 
     If(lPartUp) Then
	 Call IniSolute(t)		 
	 If (l_degrad) Then
           Call CalcFieldCapacity
           Call CalcDepthFactor
        End If
     End If
     Allocate (SoilSoluteConcentration(nElm))
     Allocate (sinkE(nElm))
     Allocate (SoilSoluteMass(nElm))
     Allocate (SoilSoluteMass_new(nElm))
     Allocate (Uptake(1))
     Allocate (Cum_uptake(1))
     Uptake=0._dp
     Cum_uptake=0._dp
  End If
  Call BroadcastInt(nElm)
  If (continu) Then
  Call BroadcastInt((nx+1)*(ny+1)*nz)
  else
  Call BroadcastInt(nPt)
  End if
  If(continu) periodic=1
  Call BroadcastInt(periodic)
#endif

  If(master) Then
     ! Initiate Rhizo if rhizodynamics is considered
     If (lRhizo .And. RhizoModel .Eq. 3) Then
        Do jj=1,nPt
           thetaNonEqOld(jj) = thetaTot(jj)
        Enddo
     Endif
  Endif


!> setup partrace field before running rswms to get the correct amount of solute mass
#ifdef WITH_PARTRACE
  If(master) Then
        Call Veloc
  Endif
  If(master) Then
        !> Setup partrace field
        Select Case(uptakeorder)
        Case(1) !exclusion
           sinkE = 0._dp
        Case(2) !passiv uptake, advection only
           sinkE = sink_cube
        Case(3) !active uptake, Michaelis Menten
           sinkE = sink_cube
        Case(4) !passiv uptake, incl. diffusive uptake
           sinkE = fact
        Case DEFAULT
           Call stop_program('Sorry, uptakeorder in Input_Partrace_etc. has to be set to a value between 1 and 4. Program terminated.')
        End Select   
  Endif
 
  Call setpartracefields(nPt,Vx,Vy,Vz,theta_old,sinkE,fac_deg,nElm,periodic)
  Call setConcentration(nElm, SoilSoluteConcentration,SoilSoluteMass)
#endif
  Timeloop: Do While (.Not.(Abs(t-tMax).Le.0.001_dp*dt))

     If(master) Then
        Call timestep(t,dt)
        If(master) Write (*,'(/,a,1pe12.5,a,1pe11.3,a,i6,a,1pe11.3)') ' t = ',t,' dt = ',dt,' iter = ',iter_tot, ' RelE = ',WBalR
     Endif

#ifdef WITH_PARTRACE     
     ! First call veloc subroutine to get velocities at the beginning of the time step, and save them as vx_old, vy_old and vz_old to pass onto partrace
     If(master) Then
        Call Veloc
        Vx_old=Vx
        Vy_old=Vy
        Vz_old=Vz
     Endif
#endif
	 
     If(master) Then
        Call RSWMS(t,dt)
        ! Updating thetaNonEqOld
        If (lRhizo .And. RhizoModel .Eq. 3) Then
           Do jj=1,nPt
              thetaNonEqOld(jj) = thetaNonEq(jj)
           Enddo
        Endif
     Endif
	 
#ifdef WITH_PARTRACE
     Call BroadcastReal(t)
     Call BroadcastReal(dt)

     If(master) Then
        Call Veloc
        dtPt = dt!*timestepfactor !time step size particle trackers
        timePt = t-tstart
     Endif
     Call BroadcastReal(dtPt)
     Call BroadcastReal(timePt)

     If(master) Then
           dummy_soil =  Sum(SoilSoluteMass)  ! solute mass in soil before calling partrace
           !> Calculate solute uptake using partrace and SoluteRoot
      If(lPartUp) Then
              If(l_degrad) Call DegradSoil(t)	
              SoilSoluteUptake = 0._dp
              !Do ipl=1,nplant
              Call SoluteRoot(icount,dt,t,ipl)
              !End do
              dummy_root = Sum(SoilSoluteUptake)  ! solute mass taken up by roots
              !print*,'MAIN sum solute uptake', dummy_root
      End If
           Select Case(uptakeorder)
           Case(1) !exclusion
              sinkE = 0._dp
           Case(2) !passiv uptake, advection only
              sinkE = sink_cube
           Case(3) !active uptake, Michaelis Menten
              sinkE = sink_cube
           Case(4) !passiv uptake, incl. diffusive uptake
              sinkE = fact
              
           Case DEFAULT
              Call stop_program('Sorry, uptakeorder in Input_Partrace_etc. has to be set to a value between 1 and 4. Program terminated.')
           End Select
      Endif
        
      Call setpartracefields(nPt,Vx_old,Vy_old,Vz_old,theta_old,sinkE,fac_deg,nElm,periodic)
      Call runpartrace(timePt,dt)
      Call setConcentration(nElm, SoilSoluteConcentration, SoilSoluteMass)
      If(master) Then
         Print*, 'uptake from partrace', dummy_soil - Sum(SoilSoluteMass), dummy_soil, Sum(SoilSoluteMass)
         Print*, 'uptake by roots', dummy_root
         dummy_diff = (dummy_soil - Sum(SoilSoluteMass)) + dummy_root
         Uptake=dummy_root
         Cum_uptake=Cum_uptake+dummy_root
		   
         If(lSalinity) Call Salinity
      Endif
 
#endif



#ifdef WITH_PARTRACE
     If(savelast) Then
        isavelast=1
     Else
        isavelast=0
     Endif
     Call BroadcastInt(isavelast)
     If(isavelast==0) Then
        savelast=.False.
     Else
        savelast=.True.
     Endif
#endif
     If (savelast) Then
        If(master) Then
           Call OutFEM(t,0)
           If(lvtk) Call OutVTK(kOut)
           Call FlxOut(kOut) !SWMS3D
           If(lOutPartrace) Call PartraceOut(t)
           If (lDou) Then
              Do ipl=1,nplant
              Call OutDou(t,0,ipl)
              End do
              If(lvtk) Then
              Do ipl=1, nplant
              Call OutDouVTK(0,t,ipl)
              End do
              End if
           Endif
           Open (unit=9,File='SimulationTime.out',status='unknown')
        Endif
        Exit TimeLoop
     Endif
     ! end of simulation?:
     If (Abs(t-tMax).Le.0.001_dp*dt) Then
        If(master) Then
           Print *,""
		   Print *, "--------------------------------"
           Print *, "End of your simulation: all went fine! :-)." 
           if (i_noConv.NE.0) Print *, "There were however ",i_noConv," time step(s) without convergence."
           Open (unit=9,File='SimulationTime.out',status='unknown')
           Deallocate(theta_old,hOld,hTemp,hNew,theta,conc)
           Deallocate(Kode,Q)
           Deallocate(conO,con,cap)
           Deallocate(WatIn,SolIn)
        Endif
        Exit TimeLoop
     Endif
  End Do Timeloop

#ifdef WITH_PARTRACE
  If(master) Then
     Deallocate(Vx,Vy,Vz,Vx_old,Vy_old,Vz_old)
  Endif
  Call closePartrace()
#endif

End Program RSWMS_MAIN
