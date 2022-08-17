! ==============================================================================
! Source file ENVIRONMENTAL FUNCTIONS ||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!> \file Environmental.f90
!! \brief module with environmental functions
Module Environmental

Contains
  ! current spatial soil strength distribution:
  SUBROUTINE Solstr(h)
    USE typedef
    USE GridData
    USE GridData2
    USE SolData
    USE WatFun
    USE StrData
    USE RhizoData, ONLY : bulkPara,  lRhizo
    USE RootData, ONLY : l_overburden, ltwo_grids, lHydrotrop, lSoilStrength, lLandl, lBengough1
    USE DomData
    IMPLICIT NONE

    INTEGER(ap) :: i,m,nPt_
    REAL(dp), INTENT(in) ::  h(nPt)
    REAL(dp) :: sat,dummy,por,voidratio,stress,h2(nPt2),Ermax,h_zero
    REAL(dp),   Allocatable, Dimension(:) :: h_use
    INTEGER(ap),   Allocatable, Dimension(:) ::  MatNum_

    IF (ltwo_grids) THEN
         CALL assignement(h, h2)
         nPt_ = nPt2
         ALLOCATE (h_use(nPt_))
         h_use = h2
         ALLOCATE (MatNum_(nPt_))
         MatNum_ = MatNum2
    ELSE 
         nPt_ = nPt
         ALLOCATE (h_use(nPt_))
         h_use = h
         ALLOCATE(MatNum_(nPt_))
         MatNum_ = MatNum
    END IF 

    DO i=1,nPt_
       M=MatNum_(i)
	   
!!!!* assign current soil strength values to each node:
       IF (soiltab) THEN!(Couvreur nov 2011)
          s(i)=((1._dp-(Fth_soiltab(h_use(i),M)-TheTab(nTab,M))/(TheTab(1,M)-TheTab(nTab,M)))**3)*ssMaxTab(M)
       ELSE
          IF (.NOT. lRhizo) THEN
             sat=(FTh(h_use(i),par(:,M))-par(2,M))/(par(3,M)-par(2,M)) 
          ELSE
             sat=(FTh(h_use(i),par(:,M))-bulkPara(1))/(bulkPara(2)-bulkPara(1)) 
          ENDIF
          IF (sat .LE. 0._dp) sat=0.001_dp
          IF (l_overburden) THEN !approach by Gao et al (2016)
             por=1-par(11,M)/2.65 !estimated porosity
             IF (por .EQ. 1._dp) por=0.99_dp
             voidratio=por/(1-por)
             stress=9.81*par(11,M)*ABS(zGrid(i))*0.01 !overburden
             IF (sat .LT. 0.5) sat=0.5_dp
             s(i)=9.81*par(11,M)*((3.560-voidratio)**2/(1+voidratio)*((stress**2.1931)-(h_use(i)*0.1*sat))**0.1846)**2*0.001 !MPa
          ELSE !approach by Whalley et al (2016)
             IF (lSoilStrength) THEN
                dummy=0.35*LOG10(abs(h_use(i))*0.1*sat)+0.93*par(11,M)+1.26 !kPa
             ELSE 
                dummy=1.E-20_dp
             ENDIF
             s(i)=10**dummy*0.001 !MPa
 
          ENDIF
       ENDIF
      ! Change 24.05.22 
       
!!!!* current local impedance factor due to soil strength:
       !local maximum penetrometer resistance, Bengough(2011)
      h_zero = h_use(i)
      IF(.NOT. lHydrotrop) THEN
        h_zero = 0._dp
      ELSEIF(lLandl) THEN
        ssmax = 4._dp-2.33_dp*h_zero*0.0001_dp !MPa
        IF(ssmax.LE.0._dp) ssmax = s(i)
        imps(i) = 1._dp-(s(i)/ssmax)**.5
      ELSE 
        !Elongation rate predicted in eq 1 Bengough 2011  !!!  fitted for maize, should be check for other species (Adrien)
        IF ((s(i).LE. 2.4_dp) .AND. lBengough1) THEN  ! Within the range of validity for eq 1
            Ermax = 31.1_dp -17.8_dp * s(i) + 39.3_dp * h_zero*0.0001_dp + 4._dp*s(i)**2 + 24.7_dp * (h_zero*0.0001)**2 + 3.3_dp*s(i)*(h_zero*0.0001)
            IF (Ermax.LT.0._dp) THEN
                imps(i) = 0._dp
            ELSEIF (Ermax.GT.31.1_dp) THEN ! 31 if EQ 1 Bengough 2011 or 48 if EQ 2 Bengough 2011
                imps(i) = 1._dp
            ELSE 
                imps(i) = Ermax/31.1_dp
            ENDIF
        ELSE ! Outside the range for eq 1 Bengough 2011
            Ermax = 48._dp + 28._dp*h_zero*0.0001_dp-12._dp*s(i)
            IF (Ermax.LT.0._dp) THEN
                imps(i) = 0._dp
            ELSEIF (Ermax.GT.48._dp) THEN ! 31 if EQ 1 Bengough 2011 or 48 if EQ 2 Bengough 2011
                imps(i) = 1._dp
            ELSE 
                imps(i) = Ermax/48._dp
            ENDIF
        ENDIF
      ENDIF
      

    END DO


    RETURN
  END SUBROUTINE Solstr
  !***********************************************************
  ! current spatial temperature distribution in the soil:
  SUBROUTINE Temper(t)
    USE ParamData, ONLY: pi
    USE TempData
    USE GridData
    USE GridData2
    USE RootData, ONLY: nplant,ltwo_grids
    IMPLICIT NONE

    INTEGER(ap) :: i,it,it1,it2,iz1,iz,iz2,ipl, nPt_
    REAL(dp) :: val1,val2,depfac,t,timfac
    REAL(dp), Allocatable, Dimension(:) :: zGrid_

	
!!!!* find the current time interval:
    IF (t.GE.time_S(nt_tempS)) THEN
       it1=nt_tempS
    ELSE
       it=nt_tempS
1      it=it-1
       IF ((t.LT.time_S(it)).AND.(it.GT.1)) GOTO 1
       it1=it
       it2=it+1
       timfac=(t-time_S(it1))/(time_S(it2)-time_S(it1))
    ENDIF
    IF (ltwo_grids) THEN
         nPt_ = nPt2
         ALLOCATE(zGrid_(1:nPt_))
         zGrid_ = zGrid2
    ELSE 
         nPt_ = nPt
         ALLOCATE(zGrid_(1:nPt_)) 
         zGrid_ = zGrid
    ENDIF
    DO  i=1,nPt
!!!!* assign current temperature values to each node --
!!!!* find the appropriate depth interval:
       IF (zGrid(i).LE.zGrid(1)-depth(nz_tempS)) THEN
          iz1=nz_tempS
       ELSE
          iz=nz_tempS
11        iz=iz-1
          IF ((zGrid(i).GT.zGrid(1)-depth(iz)).AND.(iz.GT.1)) GOTO 11
          iz1=iz
          iz2=iz+1
          depfac=(zGrid(i)-(zGrid(1)-depth(iz1)))/((zGrid(1)-depth(iz2))-(zGrid(1)-depth(iz1)))
       ENDIF
!!!!* interpolate along depth- and time-coordinates:
       IF (iz1.EQ.nz_tempS) THEN
          IF (it1.EQ.nt_tempS) THEN
             tem(i)=temtim(nt_tempS,nz_tempS)
          ELSE
             tem(i)=timfac*(temtim(it2,nz_tempS)-temtim(it1,nz_tempS))+temtim(it1,nz_tempS)
          ENDIF
       ELSE
          IF (it1.EQ.nt_tempS) THEN
             tem(i)=depfac*(temtim(nt_tempS,iz2)-temtim(nt_tempS,iz1))+temtim(nt_tempS,iz1)
          ELSE
             val1=depfac*(temtim(it1,iz2)-temtim(it1,iz1))+temtim(it1,iz1)
             val2=depfac*(temtim(it2,iz2)-temtim(it2,iz1))+temtim(it2,iz1)
             tem(i)=timfac*(val2-val1)+val1
          ENDIF
       ENDIF
       !* current local impedance factor due to soil temperature:
       DO ipl=1,nplant
       IF ((tem(i).GE.tempermax(ipl)).OR.(tem(i).LE.tempermin(ipl))) THEN
          impt(i)=0.0_dp
       ELSE
          IF (topt(ipl).LT.tmid(ipl)) THEN
             impt(i)=SIN(pi*((tem(i)-tempermin(ipl))/trange(ipl))**expo(ipl))
          ELSE
             impt(i)=SIN(pi*((tem(i)-tempermax(ipl))/(-trange(ipl)))**expo(ipl))
          ENDIF
       ENDIF
       ENDDO
    END DO
    RETURN
  END SUBROUTINE Temper
  !**********************************************************************
  ! current nodal solute concentration impedance factor values
  SUBROUTINE ConTox(ipl)
    USE GridData
    USE ConData
    USE SolData, ONLY: Conc
    IMPLICIT NONE

    INTEGER(ap) :: i,ipl
    REAL(dp) :: c0,c1,c2,c3

    !* current local impedance factor due to soil water solution concentration:
    c0=cmin(ipl)
    c1=coptmi(ipl)
    c2=coptma(ipl)
    c3=cmax(ipl)
    DO i=1,nPt
       impc(i)=0.0_dp
       IF (Conc(i).GT.c0.AND.Conc(i).LT.c1) impc(i)=(Conc(i)-c0)/(c1-c0)
       IF (Conc(i).GE.c1.AND.Conc(i).LE.c2) impc(i)=1._dp
       IF (Conc(i).GT.c2.AND.Conc(i).LT.c3)impc(i)=(Conc(i)-c3)/(c2-c3)
    END DO
    RETURN
  END SUBROUTINE ConTox
  !***********************************************************************************
  ! current local soil strength value:
  SUBROUTINE StrLoc(corner,sLoc)
    USE Typedef
    USE StrData
    IMPLICIT NONE

    INTEGER(ap), INTENT(in):: corner(8)
    REAL(dp), INTENT(out):: sLoc
	

    sLoc=(s(corner(1))+s(corner(2))+s(corner(3))+s(corner(4))+s(corner(5))+&
         s(corner(6))+s(corner(7))+s(corner(8)))/8._dp
    RETURN
  END SUBROUTINE StrLoc
  !***********************************************************************************
  ! current local temperature value:
  SUBROUTINE TemLoc(corner,tLoc)
    USE TempData
    USE RootData,ONLY: ltemp
    IMPLICIT NONE

    INTEGER(ap):: corner(8),ipl=1
    REAL(dp):: tLoc

    IF (ltemp) THEN
       tLoc=(tem(corner(1))+tem(corner(2))+tem(corner(3))+tem(corner(4))+tem(corner(5))&
            +tem(corner(6))+tem(corner(7))+tem(corner(8)))/8._dp
    ELSE
       tLoc=topt(ipl)
    ENDIF
    RETURN
  END SUBROUTINE TemLoc
  !***********************************************************************************
  !Atmospheric conditions
  SUBROUTINE Atmosphere(t)
    USE TempData
    IMPLICIT NONE

    INTEGER(ap) :: aa
    REAL(dp), INTENT(in) :: t

    !interpolate atmospheric temperature, T_atm for time
    IF (t.GE.time_TA(nt_tempS)) THEN
       Tatm_usr=T_atm(nt_tempS)
    ELSE
       DO aa=1,nt_tempS
          IF (t.GE.time_TA(aa)) THEN
             Tatm_usr=(T_atm(aa+1)-T_atm(aa))/(time_TA(aa+1)-time_TA(aa))*t+T_atm(1)
          ENDIF
       ENDDO
    ENDIF
    Tatm_usr=Tatm_usr+273

    !interpolate atmospheric pressure, P_atm for time
    IF (t.GE.time_PA(nt_presA)) THEN
       Patm_usr=P_atm(nt_presA)
    ELSE
       DO aa=1,nt_presA
          IF (t.GE.time_PA(aa)) THEN
             Patm_usr=(P_atm(aa+1)-P_atm(aa))/(time_PA(aa+1)-time_PA(aa))*t+P_atm(1)
          ENDIF
       ENDDO
    ENDIF

    !interpolate atmospheric pressure deficit, P_diff for time
    IF (t.GE.time_PD(nt_presD)) THEN
       Pdiff_usr=P_diff(nt_presD)
    ELSE
       DO aa=1,nt_presD
          IF (t.GE.time_PD(aa)) THEN
             Pdiff_usr=(P_diff(aa+1)-P_diff(aa))/(time_PD(aa+1)-time_PD(aa))*t+P_diff(1)
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE Atmosphere
!***********************************************************************************
   SUBROUTINE assignement(h, h2)
        ! assigns suitable h values from the coordinates of the coarse grid 
        ! to the coordinates of the fine grid 
     USE Typedef
     USE GridData
     USE GridData2
     IMPLICIT NONE
		
     REAL(dp), INTENT(in) ::  h(nPt)
     REAL(dp), INTENT(out) ::  h2(nPt2)
     REAL(dp),   Allocatable, Dimension(:) :: xarray,yarray, zarray
     INTEGER(ap) :: xlen, ylen, zlen, xyplane, i, xidx,yidx, zidx
     INTEGER(ap) :: idxb(nPt2)
		 
      !make unique arrays for x, y, z 
      CALL unique(xGrid, xarray, xlen)
      CALL unique(yGrid, yarray, ylen)
      CALL unique(zGrid, zarray, zlen)

      xyplane = xlen*ylen 
		
        DO i=1,nPt2
         xidx = MINLOC(ABS(xarray-XGrid2(i)),DIM=1)
         yidx = MINLOC(ABS(yarray-YGrid2(i)),DIM=1)
         zidx = MINLOC(ABS(zarray-ZGrid2(i)),DIM=1)
         idxb(i) = (zidx-1) *xyplane +(yidx-1)*xlen + xidx
        END DO
		
      h2 = h(idxb)
   END SUBROUTINE assignement
!***********************************************************************************
   SUBROUTINE unique(array_in, array_out, k)
    ! make arrays with unique numbers from the x,y,z coordinates of the coarse grid 
     USE Typedef
     IMPLICIT NONE
		
     REAL(dp), Allocatable, Dimension(:), intent(in)  ::  array_in   ! The input
     REAL(dp), Allocatable, Dimension(:), intent(out) ::  array_out  ! The output
     INTEGER(ap), intent(out) :: k                                   ! The number of unique elements
     REAL(dp) :: res(size(array_in))
     INTEGER(ap) :: i, j
 
      k = 1
      res(1) = array_in(1)
      outer: do i=2,size(array_in)
                 do j=1,k
                    if (res(j) == array_in(i)) then
                      ! Found a match so start looking again
                           cycle outer
                    end if
                 end do
			! No match found so add it to the output
                 k = k + 1
                 res(k) = array_in(i)
             end do outer
             array_out = res(1:k)
   END SUBROUTINE unique
! ***********************************************************************************
   ! Calculation of correction factor for solute degradation in soil
   ! will be handed to Partrace and multiplied to first or second order decay factor
   !Dependency on soil moisture, depth and temperature
   !Contains subrouines: 
   !- CalcFieldCapacity: caluclates field capacity for each soil material
   !- CalcMoistFact: calculates moisture factor for each soil node. 
   !- CalcElmTheta: calculates water content on element basis
   !- LocateDFactor: Locates the depth dependent decay factor for each voxel 
   !- InterpSoilTemp: interpolates the soil temperature in space and time and creates element based temp values 

   SUBROUTINE DegradSoil(t)
     USE typedef
     USE GridData, ONLY: nElm
     USE SoluteRootMat, ONLY:fac_temp, fac_dept, fac_wet, fac_deg,decayrate
     Implicit None

     INTEGER(ap) :: iE
     REAL(dp) :: t
     !fac_temp -> temperature factor
     !ElmDfac -> depth factor
     !fm -> moisture factor
     
     Call CalcMoistFactor
     Call TemperElm(t)

     DO iE=1,nElm
        fac_deg(iE)=fac_temp(iE)*fac_wet(iE)*fac_dept(iE)*decayrate
        !print*, fac_deg(iE), iE, 'total factor'
     ENDDO
     
   END SUBROUTINE DegradSoil
!*****************************************************
   SUBROUTINE CalcMoistFactor
     USE typedef
     USE SoluteRootMat, ONLY: thFC, fac_wet,theta_elm
     USE GridData, ONlY: nElm
     Implicit None
     
     REAL(dp):: B
     INTEGER(ap)::iE
     ! IN each time step, call current water content of each node
     ! Calculate Water content of each element
     ! Associate element with material (and therfore field capacity)
     ! calculate moisture factor
     ! Wenn 2 material aufeinander treffen, kann der fc wert nicht gemittelt werden, da nicht lineare zusammehänge exestieren. stattdessen müssen die MVG parameter gemittelt werden?!
     print*, 'Calculate moisture dependent degradation'
     B=0.7   ! B is empirical factor -> FOCUS report 2000, p. 92, PEARL manual, Walker 1974. 
     ! for now only consider one soil material
     CALL CalcElmTheta
     
     DO iE=1,nElm
        
        fac_wet(iE)=(theta_elm(iE)/thFC(iE))**B
        IF (fac_wet(iE).GT.1.0_dp) THEN
           fac_wet(iE)=1.0_dp
        ENDIF
        !print*,theta_elm(iE),thFC(iE), fac_wet(iE)
     END DO
     
     
   END SUBROUTINE CalcMoistFactor
   !********************************
   ! Subroutine uses depthdecay.in 
   ! locate in each soil voxel a depth dependent decay factor
   SUBROUTINE CalcDepthFactor 
     USE typedef
     USE GridData, ONLY: elmnod, nElm, zgrid
     USE SoluteRootMat, ONLY: ndepth, Ddepth, dfactor, fac_dept
     Implicit None

     INTEGER(ap):: iE, nd
     REAL(dp)::  dElm(nElm)
     
     print*, 'Calculate Depth Factor'
        
     ! find for each soil element (center) the soil depth
     
     DO iE=1,nElm 
        dElm(iE)=abs((zgrid(elmnod(1,iE))+zgrid(elmnod(5,iE)))/2)
     ENDDO
     
     
     ! for each soil element find the fitting depth decay factor
     
     DO iE=1,nElm
        DO nd=1, ndepth-1
           IF (dElm(iE).GT.(Ddepth(nd)))THEN
              IF (dElm(iE).LT.(Ddepth(nd+1)))THEN
                 fac_dept(iE)=dfactor(nd)
              ELSEIF (dElm(iE).EQ.Ddepth(nd+1))THEN
                 fac_dept(iE)=dfactor(nd+1)
              ELSEIF(dElm(iE).GE.(Ddepth(ndepth)))THEN
                 fac_dept(iE)=dfactor(ndepth)
              END IF
              
           END IF
        END DO
        !print*, iE,dElm(iE), fac_dept(iE), 'depth factor'
     END DO
     
   END SUBROUTINE CalcDepthFactor
   !****************************************************************************************
   !calculate temperature factor for decay reduction
   !ft=exp(Ea/(Rg*(T-Tr)))
   ! where Ea [J mol-1] is activation energy, Rg is gas constant[J mol-1 K-1], T actual temperature [K], Tr reference temperature  [K]
   
   SUBROUTINE CalcTempFact(t)
     USE typedef
     USE GridData, ONLY: nElm
     USE SoluteRootMat, ONLY: temp_elm,fac_temp
     Implicit None

     REAL(dp) :: Rg, Ea, Tref,t
     INTEGER(ap) :: iE
     
     Rg   = 8.3144598    !gas constant
     Ea   = 54._dp       !activation energy (default in PEARL, Pearl manual p. 51)
     Tref = 20._dp
     
     print*, 'Calculate temperature factor'
     
     Call TemperElm(t)
     
     DO iE=1,nElm
        
        IF (temp_elm(iE).LE.0._dp)THEN
           fac_temp(iE)=0._dp
        ELSE IF (temp_elm(iE).GE.35._dp)THEN
           fac_temp(iE)=exp(-Ea/Rg*(1/(35.0+273.15) - 1/(Tref+273.15)))
        ELSE
           fac_temp(iE)=exp(-Ea/Rg*(1/(temp_elm(iE)+273.15)-1/(Tref+273.15)))
        ENDIF
        !print*,temp_elm(iE), fac_temp(iE), iE, 'TempFActor'
     END DO
     
   END SUBROUTINE CalcTempFact
   !****************************************************************** 
   ! calcualate thFC(iE) for each element
   SUBROUTINE CalcFieldCapacity
     USE typedef
     USE SolData, ONLY: par, nmat, MatNum
     USE SoluteRootMat, ONLY: thFC
     USE GridData, ONLY: elmnod, nElm
     Implicit None

     INTEGER(ap)::mat, iE, matcor(8),cor(8)
     REAL(dp):: thr(nmat), ths(nmat), alp(nmat), nfa(nmat), pF, divi(nElm), high(nElm),thr_elm(nElm),ths_elm(nElm), alp_elm(nElm), nfa_elm(nElm)
       
     pF=2.0_dp     !field capacity is reached at pF=2 (accoring to FOUCS)
     
     DO mat=1,nmat   !Van Genuchten Parameters for each soil material
        thr(mat)=par(2,mat)
        ths(mat)=par(3,mat)
        alp(mat)=par(4,mat)
        nfa(mat)=par(5,mat)
     END DO


     DO iE=1,nElm
        ! find nodes of each element
        cor(1)=elmnod(1,iE)
        cor(2)=elmnod(2,iE)
        cor(3)=elmnod(3,iE)
        cor(4)=elmnod(4,iE)
        cor(5)=elmnod(5,iE)
        cor(6)=elmnod(6,iE)
        cor(7)=elmnod(7,iE)
        cor(8)=elmnod(8,iE)
        
        matcor(1)=MatNum(cor(1))   ! material of each nodes of the cube
        matcor(2)=MatNum(cor(2))
        matcor(3)=MatNum(cor(3))
        matcor(4)=MatNum(cor(4))
        matcor(5)=MatNum(cor(5))
        matcor(6)=MatNum(cor(6))
        matcor(7)=MatNum(cor(7))
        matcor(8)=MatNum(cor(8))
        
        ! find voxel based values for van genuchten parameters
        thr_elm(iE)=(thr(matcor(1))+thr(matcor(2))+thr(matcor(3))+thr(matcor(4))+thr(matcor(5))+thr(matcor(6))+thr(matcor(7))+thr(matcor(8)))/8
        
        ths_elm(iE)=(ths(matcor(1))+ths(matcor(2))+ths(matcor(3))+ths(matcor(4))+ths(matcor(5))+ths(matcor(6))+ths(matcor(7))+ths(matcor(8)))/8
        
        alp_elm(iE)=(alp(matcor(1))+alp(matcor(2))+alp(matcor(3))+alp(matcor(4))+alp(matcor(5))+alp(matcor(6))+alp(matcor(7))+alp(matcor(8)))/8
        
        nfa_elm(iE)=(nfa(matcor(1))+nfa(matcor(2))+nfa(matcor(3))+nfa(matcor(4))+nfa(matcor(5))+nfa(matcor(6))+nfa(matcor(7))+nfa(matcor(8)))/8
     END DO



    
     ! CHECK THIS !!!!!
     DO iE=1,nElm
        
        divi(iE)=(1+(alp_elm(iE)*10**pF)**nfa_elm(iE))
        high(iE)=1-(1/nfa_elm(iE))
        
        thFC(iE)=thr_elm(iE)+((ths_elm(iE)-thr_elm(iE))/(divi(iE)**high(iE)))
     END DO
     
   END SUBROUTINE CalcFieldCapacity
   !**********************
   SUBROUTINE CalcElmTheta
     USE typedef
     USE GridData, ONLY: elmnod, nElm
     USE SolData, ONLY: theta
     USE SoluteRootMat, ONLY: theta_elm
     
     INTEGER(ap)::cor(8), iE
     

     IF (.NOT. ALLOCATED(theta_elm)) THEN
        ALLOCATE(theta_elm(nElm))
        theta_elm = 0._dp
     END IF
     
     DO  iE=1,nElm
        cor(1)=elmnod(1,iE)
        cor(2)=elmnod(2,iE)
        cor(3)=elmnod(3,iE)
        cor(4)=elmnod(4,iE)
        cor(5)=elmnod(5,iE)
        cor(6)=elmnod(6,iE)
        cor(7)=elmnod(7,iE)
        cor(8)=elmnod(8,iE)
        
        theta_elm(iE) = (theta(cor(1))+theta(cor(2))+theta(cor(3))+theta(cor(4))+theta(cor(5))+&
             theta(cor(6))+theta(cor(7))+theta(cor(8)))/8
     ENDDO
     
   END SUBROUTINE CalcElmTheta
   !*********************************
   ! Subroutine uses soiltemp.in
   ! interpolates the temperature values in time and in depth
   ! copied and adapted from Subroutine Temper in environmental.f90, from SWMS code
   ! interpolates temperature values for soil ELEMENTS
   SUBROUTINE TemperElm(t)
     USE typedef
     USE GridData, ONLY: elmnod, nElm, zgrid, nPt
     USE SoluteRootMat, ONLY:  temp_elm
     USE TempData
     Implicit None

     REAL(dp) :: timfac, depfac, temN(nPt)
     REAL(dp) :: t, val1, val2
     INTEGER(ap) :: i,it,it1,it2, iE, cor(2),iz1,iz,iz2
     
       
     print*, 'Interpolate soil Temp in time and space'     
  
     
!!!!* find the current time interval:
     IF (t.GE.time_S(nt_tempS)) THEN
        it1=nt_tempS
     ELSE
        it=nt_tempS      
        IF ((t.LT.time_S(it)).AND.(it.GT.1)) it=it-1
        it1=it
        it2=it+1
        timfac=(t-time_S(it1))/(time_S(it2)-time_S(it1))
     ENDIF

  
     DO  i=1,nPt
!!!!* assign current temperature values to each node --
!!!!* find the appropriate depth interval:
        IF (zGrid(i).LE.zGrid(1)-depth(nz_tempS)) THEN
           iz1=nz_tempS
        ELSE
           iz=nz_tempS
           IF ((zGrid(i).GT.zGrid(1)-depth(iz)).AND.(iz.GT.1)) iz=iz-1
           iz1=iz
           iz2=iz+1
           depfac=(zGrid(i)-(zGrid(1)-depth(iz1)))/((zGrid(1)-depth(iz2))-(zGrid(1)-depth(iz1)))
        ENDIF
        !print*, depfac, 'depth factor'
!!!!* interpolate along depth- and time-coordinates:
        IF (iz1.EQ.nz_tempS) THEN
           IF (it1.EQ.nt_tempS) THEN
              temN(i)=temtim(nt_tempS,nz_tempS)
           ELSE
              temN(i)=timfac*(temtim(it2,nz_tempS)-temtim(it1,nz_tempS))+temtim(it1,nz_tempS)
           ENDIF
        ELSE
           IF (it1.EQ.nt_tempS) THEN
              temN(i)=depfac*(temtim(nt_tempS,iz2)-temtim(nt_tempS,iz1))+temtim(nt_tempS,iz1)
           ELSE
              val1=depfac*(temtim(it1,iz2)-temtim(it1,iz1))+temtim(it1,iz1)
              val2=depfac*(temtim(it2,iz2)-temtim(it2,iz1))+temtim(it2,iz1)
              temN(i)=timfac*(val2-val1)+val1
              !print*, temN(i), 'z-coordinate', zgrid(i)
           ENDIF
        ENDIF        
     ENDDO


! calculate temp for each element, linear interpolation between top and bottom soil voxel planes
 DO  iE=1,nElm
     cor(1)=elmnod(1,iE)
     cor(2)=elmnod(5,iE)
    
     temp_elm(iE) = (temN(cor(1))+temN(cor(2)))/2
     !print*, temp_elm(iE), iE
  ENDDO

END SUBROUTINE TemperElm
End Module Environmental
