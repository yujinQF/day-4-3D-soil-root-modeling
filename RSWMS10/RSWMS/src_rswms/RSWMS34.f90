!> \file RSWMS34.f90
!! \brief Main RSWMS function.
!!
!!Makefile .   To switch OpenMp on: pgf90 -mp -o executable file.f90
!!                                   export OMP_NUM_THREADS= 2
!!                                   ./executable
!! To use openmp -> switch module OMP_LIB on!!! gfortran compiler can´t read OMP_LIB
!! for gfortran only new compiler (v3) supports openmp
!!******************************************************************************
!!              Mathieu Javaux- Tom Schroeder- Valentin Couvreur-
!!              Natalie Schroeder- Katrin Huber- Felicien meunier -
!! 				Magdalena Landl - Helena Jorda- Nimrod Schwartz -
!!              Axelle Koch 
!!                            Jan Vanderborght
!!                            Fz-Juelich (Germany)
!!
!!
!!
!! R-SWMS_3D :An algorithm for three-dimensional, simultaneous modeling
!! of root growth, transient soil water flow, solute transport, and root
!! water and solute uptake.
!! coupled with RootTyp (Pages et al.)
!!
!! Version 10           August 2020
!!
!! Program Version for Fortan 90
!!
!! based on the code of
!! Francesca Somma & Volker Clausnitzer, University of California, Davis (1998)
!!******************************************************************************

Module RSWMSMod

Contains

  Subroutine RSWMS(t,dt)

    Use typedef
    Use CumData
    Use PlntData
    Use tmctrl
    Use GridData
    Use RootData
    Use TempData
    Use MatData
    Use DoussanMat
    Use SparseMatrix
    Use soldata
    Use ObsData
    Use DomData
    Use ParamData
    Use OMP_LIB
    Use SoluteRootMat
    Use ParticlesInRoot
    Use Doussan, Only: update_Krs, calc_actual_climate, SolveRoot, SetBCroot, update_historics, &
         SetupDou, SetupCou, SSFdis, CalcWnodes
    Use Environmental, Only: Solstr,Temper, ConTox
    Use Output, Only: SubReg, OutDou, OutRoo, OutTrans, OutDouVTK, OutParticleVTK, WriteLog, &
         OutObsProbe, Zprofiles, OutFEM, OutVTK, PartraceOut, FlxOut
    Use MPIutils, Only: stop_program
    Use PlantGrowth, Only: ActTrs, Effncy, Ratio, Leaves
    Use SoluteMod, Only: Solute
    Use Sink, Only: Setsnk
    Use RootGrowth, Only: Root
    Implicit None

    !variable declaration
    Integer(ap) :: ipl,BCtp,icount
    Real(dp) :: BCr
    Real(dp) :: dt,t
    Logical :: ReDo


!> --------------------------- Start of time loop ------------------------------

   solveroot_call=0
   i_noConv=0
!> time step adjustment:
    redo=.True.
    Do While(redo)
       !call outdouvtk(0,t)
       iter=0
       iter_root=0
       If(.Not.Allocated(Q)) Allocate(Q(nPt))
!> Set soil boundary conditions
       Call SetBC(t)
!>signalling affects plant and stomatal conductance
       If (lSign_new)Then
          Call update_Krs(t)
          Call calc_actual_climate(t)
       Endif
!> apply root extraction to FEM-sink term for each node:
       If (.Not.(switchSolve)) Then
          Write (*,'(a)',advance='no') '-> '
!>solve root system and adapt BC for current time	  
          If (lDou) Then
             Call SolveRoot(t,dt,it1,iter_root)
             it1=.False.
!>update sink term	  
          Elseif ((.Not.lno_RWU).And.(.Not.lCalloc)) Then
             If (.Not.lSign_new) Call SetBCroot(t,BCr,BCtp)
             Call Setsnk(t,BCr,BCtp) !Also required when lDou because of solute uptake
          Endif
       End If

!> solve soil water flow equation and calculate pressure heads:
       If(lno_Archi) Call IniMat
       Call Water(t,dt,dtOpt,tOld,ReDo,IAD,IADN,IADD)
     End do! End of redo loop

! Water is solved, solve solute transport now
     If (lSign_new) Then
       Call update_historics(t,dt)
    Endif
!> Calculate hormonal signaling
    If(lSign .Or. lSign_inst) Then
     If (nplant.GE.2) Call stop_program ('Hormonal signaling does not work for multiple roots')
       Call IniSolute(t)
     Do ipl=1,nplant
       Call SoluteRoot(icount,dt,t,ipl)
      End do
    End If
!> solve solute transport equation and calculate concentrations:
    If (lChem) Then
       Call Solute(dt,t,tPulse,dtMaxC,KodCB,IAD,IADN,IADD,icount)
    Endif
    !> calculate water and solute (if included) mass balance
    iCount=1
    Call SubReg(t,iCount)
    !> calculate actual transpiration rate Tact (global variable in MODULE PlntData):
    If (lCalloc.And.(.Not.(lDou))) Call ActTrs
    If (lCalloc) Then
       Do ipl=1,nplant
          If (ipl.Ge.2) Call stop_program('Assimilate allocation currently doesnt work with multiple plants')
          !> translate transpired water into biomass:
          Call Effncy(t,rs,concrs,W)
          delBM=W*Tact(ipl)*dt

          !> partition that goes to new roots is accumulated in dmroot
          !> until next root growth step:
          Call Ratio(t,rs,concrs,RSR)
          dr=RSR/(1.+RSR)*delBM
          dmroot=dmroot+dr

          !> remainder of delBM is partition that goes to shoot:
          dsh=delBM-dr
          mshoot(ipl)=mshoot(ipl)+dsh

          !> calculate current leaf area:
          Call Leaves(t,LAmsh)
          LA=LA+dsh*LAmsh
       End Do
    Endif

!> Write Root related outputs
    If (Abs(t-tOut(kOut)).Le.0.001_dp*dt) Then
       Do ipl=1,nplant
          If(lDou) Then
             Call OutDou(t,kout,ipl)       !> Output for Doussan
             Write(*,'(/a)',advance='yes')'Output files correctly written'
          End If
          If(.Not.lno_root_growth) Call OutRoo(t,kOut,ipl)
          If(continu)  Call OutTrans(kout,ipl)
          If(lvtk .And. (lDou .Or. .Not.lno_root_growth)) Then
             Call OutDouVTK(kout,t,ipl)
             If(lSign .Or. lPartUp)  Call OutParticleVTK(kout,t,ipl)
          End If
       End do
    Endif
	

    !> general output:
    If(.Not.lno_RWU) Then
       Call WriteLog(t)
    Endif

    !> if appropriate, grow root and update the distribution functions 'betaw'
    !> and 'betac' (extraction intensity per volume, [L^-3])
    If (Abs(t-tCallR).Lt.0.001_dp*dt .And. (.Not.lno_root_growth)) Then

       !update soil strength values:
       Call Solstr(hNew)

       !> if temperature input provided, update temperature values:
       If (ltemp) Call Temper(t)

       !> if nutrient deficiency/ion toxicity data are provided, update
       !> the corresponding impedance nodal values:
       If (ltoxi) Then
        If (nplant.Ge.2) Call stop_program('Toxicity does not work for multiple root')
        Do ipl=1,nplant
           Call ConTox(ipl)
        End do
       End if

       !> let root system grow:
       DO ipl=1,nplant
         IF(.not.lUpdate_growth) CALL Root(t,dt,ipl)
       ENDDO
       tCallR=tCallR+dtRoot

       !> update Dousssan weighing factors and matrices
       If (lDou) Then
          !> delete the old vector of list
          Do ipl=1, nplant
             If (nrecold(ipl).Ne.0) Then
                Call SM_delete(plantmatrix(ipl))
             Endif
             If(continu) transroot(nrec(ipl)+1:nrec(ipl)+ngrow(ipl),:,:,1)=transtip(1:ngrow(ipl),:,:,1)
             Call SetupDou(t,dt,ipl)
          End do
       Endif
    Endif

    If (lCou .And. ntimeobs>1) Then
       Call SetupCou(t)
    End If

    !> in case of a static root system with ages smaller than runtime and variable conductivities over time
    If(lUpdate_growth .And. lDou) Then
       !> delete the old vector of list
       Do ipl=1, nplant
       If (nrecold(ipl).Ne.0) Then
             Call SM_delete(plantmatrix(ipl))
       Endif
       End do
       !Deallocate(nBC_irecn)
       Do ipl=1,nplant
          Call SetupDou(t,dt,ipl)
       End do
    End If
!> FEM, root output at specified points in time: (independent of DoussanMat parameters)
    If (Abs(t-tOut(kOut)).Le.0.001_dp*dt) Then
       !< in case macroscopic hydraulic properties need to be saved at each outfem time
       If (lCou.And.(lGap.Or.lAQPc)) Call SSFdis(kOut)
       Call OutFEM(t,kOut)
       If(lvtk) Call OutVTK(kOut)
       Call FlxOut(kOut)
       If(lOutPartrace) Call PartraceOut(t)
       kOut=kOut+1
       If (tOut(kOut).eq.tout(kOut-1)) Then
          Call stop_program('Error in the list of output times (2 times the same)')
       Endif
    Endif

!> ouputs for observation probes
    If (ObsOK) Then
       If (dtprobe.Eq.999) Then
          Call OutObsProbe(t)
       Elseif (Abs(t-tOuProbe(kOuProbe)).Le.0.001_dp*dt) Then
          Call OutObsProbe(t)
          kouProbe=kouProbe+1
       Endif
    Endif

    !> z-profiles
    If (profOK) Then
       If (dtprof.Eq.999) Then
          Call Zprofiles(t)
       Elseif (Abs(t-tOuProf(kOuProf)).Le.0.001_dp*dt) Then
          Call Zprofiles(t)
          kouProf=kouProf+1
       Endif
    Endif

    !> check next BC time for root
    If ((.Not.lno_RWU).And.(.Not.lCalloc)) Then
       If (Abs(t-tBCr(kBCr)).Le.0.001_dp*dt) Then
          kBCr=kBCr+1
       Endif
    Endif

    !> pressure heads for new time level:
    hOld =hNew
    hTemp=hNew

    !> ------------------------ End of time loop ----------------------------------
  End Subroutine rswms

  !**************************************************************************
  !> ### Initializing RSWMS ###
  Subroutine RSWMS_INI(t,dt)
    Use iso_fortran_env

    Use typedef
    Use ConData
    Use CumData
    Use PlntData
    Use tmctrl
    Use GridData
    Use RootData
    Use TempData
    Use MatData
    Use DoussanMat
    Use SparseMatrix
    Use soldata
    Use RhizoData
    Use ObsData
    Use DomData,  Only :zmax
    Use BoundData, Only : tQbcCh
    Use ParamData
    Use SoluteRootMat !, Only: l_degrad    
    Use disToRoot, Only: RStat
    Use Doussan, Only: SetupDou, SSFdis, CalcWnodes
    Use Input
    Use Orthofem, Only: IADMake
    Use Output, Only: OutCouVTK, ObsIni, SubReg
    Use Sink, Only: FedIn, RLDdis, CouIn, BetDis, BetNrm
    Use StrData
    Use GeoData
    !USE OMP_LIB

    Implicit None
    Integer(ap) :: i,ipl,icount
    Integer(ap) :: daytime(8)
    Integer(ap) :: nn
    Integer(ap) :: ttt
    Integer(ap), Allocatable,Dimension (:) :: seed
    Real(dp), Intent (out):: t,dt
    Character form*3,file*8,file2*15
    Integer(8) :: count

    Call RANDOM_Seed(size = nn)
    Allocate(seed(nn))

    !> initiate random number generator:
    Call DATE_AND_Time(values=daytime)
    ttt = Int(daytime(5)*3600+daytime(6)*60+daytime(7))
    Call System_Clock(count)  ! kind=8 gives micro- or nanoseconds
    ttt = Ieor(ttt, Int(count, Kind(ttt)))
    Do i = 1,nn
       seed(i) = lcg(ttt)
    End Do
   ! Print*, 'seed', seed
    Call RANDOM_Seed(put=seed)
    Deallocate(seed)
    Call RANDOM_Number(rand)

    !> open simulation summary file
    Open (UNIT=15,FILE='out/simul_sum.out',STATUS='UNKNOWN')
    Write(15,112)
112 Format('+++++ SIMULATION INPUT SUMMARY +++++ ')
    Close(15)
    ObsOK=.False.

    !> -----------------------------------INPUT---------------------------------------
    !> get input for the specific problem (domain, BC, IC, control parameters):
    Call Applic(dt)
    If (.NOT. lno_Archi .OR. lCou) Call IniRoot
    If (.NOT. lno_Archi .OR. lCou) Call IniPlant
    If ((.NOT. lno_Archi) .AND. (.NOT. lno_root_growth)) Call IniGeo
    If ((.NOT. lno_Archi) .AND. (.NOT. lno_root_growth)) Call IniStrength
    If (lDou .OR. lCou) Call IniDou
    If (lDou .OR. lCou) Call IniSol
    If ((lSomma_growth .OR. lRootBox_growth) .AND. ltoxi) Call IniConc
    If ((.NOT. lno_Archi) .AND. (.NOT. lno_root_growth)) Call IniTemp


    !allocate soil variables (module soldata)
    Allocate(Vx(nPt),Vy(nPt),Vz(nPt))
    Allocate(theta_old(nPt),theta(nPt))
    Allocate(conO(nPt),con(nPt),cap(nPt))
    Allocate(l_elmMacro(1:nElm))
    !> allocate thetaTot's if rhizosphere model.
    If (lRhizo) Allocate(thetaTot(nPt), thetaNonEq(nPt), thetaNonEqOld(nPt), tauTht(nPt), hEqRhizo(nPt))
    l_elmMacro = .False.
    
    Call CalcWnodes()!> calculate Wn (new Javaux)
    If (lChem)  Call ChemIn  !> get solute transport information
    Call SoilIn !> get soil material input and set up K,C-table for interpolation


    !> list for voxels of mixed material
    If(lDou .And. lMacro) Then
       Allocate(MacroList(1:nElm,1:6))
       MacroList = 0
       Allocate(n_neigh(1:nElm))
       n_neigh=0
       Call ListMacro
    End If

    If (.Not.(lno_Archi)) Then
       !> initial root system and growth parameters:
         DO ipl=1,nplant
!          print*,'plant number',ipl
          Call RootIn(t,ipl)
         End do
          !> get Doussan model input information
          If (lDou) Then
           Do ipl=1,nplant
             Call DouIn(ipl)
           End do
          Elseif (lFed) Then
           Do ipl=1,nplant
             Call FedIn
             Call SetupDou(t,dt,ipl)!> Calculates the root system position as compared to the position of the plant collar and of the limits of the periodic domain if so.
             Call RLDdis(t,kout,ipl)!> Calculates the root length & surface density distribution
           End do
             !         CALL BetNrm
          Elseif (lCou) Then
             Write(15,'(/''SSF, Krs and Kcomp generated from RootSys and CondRoot.in'')')
            Do ipl=1,nplant
             Call CouIn(ipl)
             Call DouIn(ipl)
             Call SetupDou(t,dt,ipl)
             Call SSFdis(kout)
             !print*,Krs,size(SSF)
             If (lvtk) Call OutCouVTK(t,ipl)
             If (ldJvL) Call RLDdis(t,kout,ipl)
            End do
          Elseif(.Not. lDou .And. (lSomma_growth.OR.(lRootBox_growth))) Then
             Allocate(transroot(0:maxrec+maxgrw,1:2,1:isubmax,1:nplant))
             transroot=0
             Allocate(transtip(0:maxgrw,1:2,1:isubmax,1:nplant))
             transtip=0
             Allocate(nsub(0:maxgrw+maxrec,1:nplant))
             nsub = 1
          End If
          If (lSign_new) Then
             Call ClimateIn
             Call TardieuIn
             Call ProcessClimate
          Endif
    
    Elseif (lCou) Then
       Call IniMat
       DO ipl=1,nplant
          Call CouIn(ipl)!Initial time determined by PlntIn
       END DO
    Elseif (lFed) Then
       Call FedIn!Initial time determined by PlntIn
    Endif

 !> soil temperature over depth and in time:
    If (ltemp .Or. lCalloc .Or. l_degrad)  Call TempIn
    If (l_degrad) Call ReadDepthDecay
    
    If ((.Not.lno_RWU).And.(.Not.lretry)) Then!> log file created even for RWU models with no root architecture
       !> open log file:
      Do ipl=1,nplant
        Write (file,'(A7)')'out/log'
        Write (file(8:8),'(I1)') ipl
        Open (UNIT=10,FILE=file,STATUS='UNKNOWN')
          If (lDou) Then
             If(lCalloc) Then
                Write (10,90)
             Elseif (lSign .Or. lSign_inst) Then
                Write (10,91)
             Elseif (lPartUp) Then
                Write (10,94)
             Elseif ((lKdrop).and.(nplant.LE.1)) Then
                write(10,922)
             Else
                Write (10,93)
             Endif
          Elseif (lCou) Then
             If (lSign_new) Then
                Write(10,95)
             Else
                Write (10,92)
             Endif
          Endif
90        Format('  Time         Tpot         Tact        grwfac      sAvg         mShoot     mRoot     LA     H_collar     sign_conc    TpLA    TpLA_pot     #stressed_Nodes ')
91        Format('  Time         Tpot         Tact         sAvg       lim. ndes    H_collar   mcol   signal_conc        msign_notrans    csign_notrans    res_time')
92        Format('  Time         Tpot         Tact      lim. ndes     H_collar     H_seqKrsKcomp')
922        Format('  Time         Tpot         Tact        sAvg        lim. ndes    H_collar   Heq_bulk    Heq_sri')
93        Format('  Time         Tpot         Tact        sAvg        lim. ndes    H_collar')
94        Format('  Time         Tpot         Tact         sAvg       lim. ndes    H_collar   mcol  root_mass      soilsolutemass    Total_mass    Uptake     Accumulated_uptake_mass  Error_root_uptake')
95        Format('  Time         Tpot         Tact         sAvg       lim. ndes    H_collar   H_seqKrsgs_t')
          Close (10)
       Enddo
    Else
       t=tQbcCh(1)!initial time =first time of the soil BC.
    Endif

    If (lSign_new) Then
       Write (file2,'(A15)')'out/Tardieu.out'
       Open (UNIT=10,FILE=file2,STATUS='UNKNOWN')
       Write(10,96)
96     Format('TimeKrsKrs_circadampliKrstranspi ABAPHcollarHleafHbundleHcellJxcVcel')
       Close (10)
    Endif

    !> get plant parameter input:
    Call PlntIn(t) !> Reads BCroot.in which is necessary even for RWU models without root architecture, and reads Plant.in which is necessary for assimilate allocation

    If (.Not.lretry) Then
       !> open balance file:
       Open (UNIT=10,FILE='out/balance.out',STATUS='UNKNOWN')
       Write(10,110)
110    Format(' Absolute and relative mass balance error for water and solute transport ')
       Close (10)
       !> open removal file
       Open (UNIT=10,FILE='out/remove.out',STATUS='UNKNOWN')
       Write(10,111)
111    Format(' Total amount of water and solute removed by the root system, ',/,&
            ' by zero and first order reactions, and by drainage at the bottom of the domain. ')
       Close(10)
       !> open Observation node file
       If (ObsOK) Then
          Call ObsIni
       Endif
       If (profOK) Then
          Open (UNIT=121,FILE='out/ProfileTH.out',STATUS='UNKNOWN')
          Open (UNIT=122,FILE='out/ProfilePH.out',STATUS='UNKNOWN')
          Open (UNIT=123,FILE='out/ProfileS.out',STATUS='UNKNOWN')
          Write(form,'(I3)')(nz)
          Write (121,'(/''Averaged water content profile for each time step.'')')
          Write (121,'(/''Time   Z- water content'')')
          Write (121,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
          Write (122,'(/''Averaged water potential profile for each time step.'')')
          Write (122,'(/''Time   Z- water potential'')')
          Write (122,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
          Write (123,'(/''Total Sink profile for each time step.'')')
          Write (123,'(/''Time   Z- sink'')')
          Write (123,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
       Endif
    Endif
    Allocate(IADN(nPt),IADD(nPt))
    Allocate(IAD(maxbnd,nPt))
    Call IADMake(nPt,nElm,maxbnd,IAD,IADN,IADD)
!--------------------------------------------------------------------------------
!> first time step:
    dtOpt=dt
    If(.Not.lno_root_growth) Then
       tCallR=t+dtRoot
    Else
       tCallR = 1e+30_dp
    End If
    dmroot=0.0_dp
    tstart = t

    !> find next output time:
1   kOut=kOut+1
    If (kOut.Le.nOut) Then
       If (tOut(kOut).Le.t) Goto 1
    Endif

    If (lSomma_growth .OR.lRootBox_growth) THEN
2    kaxemg=kaxemg+1
     Do ipl=1,nplant
      If (kaxemg.Le.naxemg(ipl)) Then
       If (tnewax(kaxemg,ipl).Lt.t) Goto 2
      Endif
     End do
    Endif
    !> Z-profiles
    touProf=touProf+t !> to start at initial time < rootin
    If ((ProfOK).And.(dtProf.Ne.999)) Then
3      kouProf=kouProf+1
       If (kouProf.Le.nouProf) Then
          If (touProf(kouProf).Le.t) Goto 3
       Endif
    Endif
    !> probes
    touProbe=touProbe+t !> to start at initial time < rootin
    If ((ObsOK).And.(dtProbe.Ne.999)) Then
4      kouProbe=kouProbe+1
       If (kouProbe.Le.nouProbe) Then
          If (touProbe(kouProbe).Le.t) Goto 4
       Endif
    Endif
    t_begin=t

    !> if Doussan, calculate weighing function and matrices
    If (lDou) Then
       Do ipl=1,nplant
          Call SetupDou(t,dt,ipl)
       !> initialize stress
          stressBC=.False.
          If ((nplant.LE.1).and.(LKdrop)) Then
            Call SSFdis(kout)
          Endif
       End do   
    Elseif (lFed) Then
       Do ipl=1,nplant
       !> get betaw and betac initial distribution:
          Call BetDis(t,ipl)
       !> ...and normalize:
          Call BetNrm
       End do
    Endif
    If ((.Not.lno_RWU).And.(.Not.lCalloc)) Then
       !> find next time change in root BC
5      kBCr=kBCr+1
       If (kBCr.Le.nBCr) Then
          If (tBCr(kBCr).Le.t) Goto 5
       Endif
    Endif

    !> calculate parameters:
    If (lRhizo)  Call RStat(t)
    Call SetMat(0,.false.)!why initiate this? (MJ20)

    Call CalcGeom
    Call CalcWidthSWMS_new
    iCount=0
    Call SubReg(t,iCount)
    !> initialize it1 (first iteration after a root growth)
    it1=.True.
    itMaxRoot=itMax
    If (oldT) Then
       switchsolve=.False.
    Else
       switchsolve=.True.
    Endif

    If(lSign) Then
      Do ipl=1,nplant
       Allocate (l_SignOn(1:nrec(ipl)+1))
       l_SignOn = .False.
      End do
    End If
  Contains
    Function lcg(s)
      Integer(ap) :: lcg
      Integer(ap) :: s
      If (s == 0) Then
         s = 104729
      Else
         s = Mod(s, 4294967)
      End If
      s = Mod(s * 279470, 4294967)
      lcg = Int(Mod(s, Int(Huge(0), ap)), Kind(0))
    End Function lcg

  End Subroutine RSWMS_INI
!******************************************************************************
!> ### adjusts next time step ###
  Subroutine timestep(t,dt)
    Use typedef
    Use tmctrl
    Use ParamData, Only: iter,mxbcch,lPartUp
    Use PlntData, Only: nBCr,tbcr
    Use RootData, Only: lCalloc,lno_RWU
    Use ObsData, Only: ObsOK, ProfOK
    Implicit None

    Real(dp) :: t,dt,dtP,tFemRoo

    tOld=t
    dtOld=dt
    If (kOut.Le.nOut) Then
       tFEMRoo=tOut(kOut)
    Else
       tFEMRoo=tMax
    Endif

    If ((.Not.lno_RWU).And.(.Not.lCalloc)) Then
       If (kBCr.Le.nBCr) Then
          tcBCr=tBCR(kBCr)
       Else
          tcBCr=tMax
       Endif
    Else
       tcBCr=tMax
    Endif
    If ((ProfOK).And.(dtProf.Ne.999).And.(kouProf.Le.nouProf)) Then
       tProf=tOuProf(kouProf)
    Else
       tProf=tMax
    Endif
    If ((ObsOK).And.(dtProbe.Ne.999).And.(kouProbe.Le.nouProbe)) Then
       tProbe=tOuProbe(kouProbe)
    Else
       tProbe=tMax
    Endif
    If(lPartUp) Then
       Call timestep_partrace(dtP)
    Else
       dtP = 1E+30_dp
    End If
    Call TmCont(iter,t,dt,dtOpt,tCallR,tFEMRoo,dtMaxC,tcBCr,tProf,tProbe)

    t=t+dt

  End Subroutine timestep
  !**************************************************************************
  !> ### adjusts partrace time step ###
  Subroutine timestep_partrace(dtP)
    Use typedef
    Use GridData, Only: dxgrid,dygrid,dzgrid,nPt
    Use SolData, Only: Vx,Vy,Vz
    Implicit None

    Integer(ap) :: ii
    Real(dp), Intent(out) :: dtP

    dtP = 1E+30_dp
    Do ii=1,nPt
       dtP = Min(dtP,dxgrid/Abs(Vx(ii)))
       dtP = Min(dtP,dygrid/Abs(Vy(ii)))
       dtP = Min(dtP,dzgrid/Abs(Vz(ii)))
    End Do
    If (dtP .Lt. 1E-20_dp) dtP = 1E+30_dp

  End Subroutine timestep_partrace

  !**************************************************************************

End Module RSWMSMod
