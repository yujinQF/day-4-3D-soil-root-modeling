!> \file RootGrowth.f90
!> \brief

!==============================================================================
! Source file ROOT GROWTH ||||||||||||||||||||||||||||||||||||||||||||||||||||||
! contains the following subroutines
! - ROOT
! - Newaxs
! - Brdevl
! - Brnage
! - Uniera
! - Length
! - Ssgcom1
! - Ssgcom2
! - SortMatrix
! - RotMatrix
! - CondDecision
! - Geocom
! - Prfang
! - Nwcomp
! - Angchg
! - Maslen
! - Mkrecd
! - Grow
! - Adjust
! - Spacng
! - Establ
! - Boxlim
! - BoxlimCylinder
! - Compon
! - Remove
! - Update
! - Secondary_Radial_Growth
! - SplitMacro
! - ShiftRoot
! ==============================================================================
!> main subroutine for root system growth
Module RootGrowth

Contains

  Subroutine ROOT(t,dt,ipl)
    Use typedef
    Use ParamData, Only: pi
    Use RootData
    Use GridData, Only: geom,continu,nex,ney,dxgrid,dygrid
    Use tmctrl, Only: kaxemg
    Use SolData, Only: conc,nmat
    Use DoussanMat, Only: transtip, transroot
    Use Environmental, Only: Atmosphere, StrLoc
    Use RootGrowthNeighb, Only: Neighb

    Implicit None
    Real(dp),Intent(in) ::t,dt
    Real(dp) ::  cloc,sloc
    Real(dp) :: newlen,impfac,newmas,sumgrw
    Real(dp) :: dx,dy,dz,space,tempage,drivex,drivey,drivez
    Integer(ap) :: corner_r(8),nrecol
    Integer(ap) :: ngrwnw,iseg,igrow,imin,irecn
    Integer(ap),Intent(in)::ipl
    !> \param t current simulation time
    !> \param dt current time step

    nrecol = nrec(ipl)
    ngrwnw = ngrow(ipl)
    sAvg   = 0.0_dp
    cAvg   = 0.0_dp
    sumgrw = 0.0_dp

    ! Secondary radial growth for all established segments
    If(l_secrad) Call Secondary_Radial_Growth(dt,ipl)

    ! check if it is time to originate new axis(es):
1   If (kaxemg.Le.naxemg(ipl)) Then
       If (t.Gt.tnewax(kaxemg,ipl)) Then
          Call Newaxs(t,sumgrw,ngrwnw,ipl)
          kaxemg = kaxemg+1
          Goto 1
       Endif
    Endif
    ! apply potential growth to all growing branches:
    tiploop:  Do igrow = 1,ngrow(ipl)
       ! develop new branches at all established points that are 'ripe':
       Call Brdevl(t,sumgrw,igrow,ngrwnw,ipl)
       ! find surrounding grid points:
       If(continu) Then
          Call tiptransGrow(xg(igrow,ipl),yg(igrow,ipl),igrow,1,1)
          Call Neighb(xg(igrow,ipl)+transtip(igrow,1,1,1)*nex*dxgrid, &
               yg(igrow,ipl)+transtip(igrow,2,1,1)*ney*dygrid,zg(igrow,ipl),corner_r,imin,ltwo_grids)
       Else
          Call Neighb(xg(igrow,ipl),yg(igrow,ipl),zg(igrow,ipl),corner_r,imin,ltwo_grids)
       End If
       ! calculate representative soil strength value felt by root tip:
       Call StrLoc(corner_r,sLoc)
       sAvg = sAvg+sLoc
       ! calculate branch age:
       Call Brnage(t,igrow,tempage,ipl)
       ! calculate new heading angle and tentative segment length:
       Call Nwcomp(igrow,corner_r,newlen,impfac,dx,dy,dz,drivex,drivey,drivez,tempage,ipl)
       ! calculate tentative segment mass:
       Call Maslen(igrow,sLoc,ordgrw(igrow,ipl),ipl)
       newmas = newlen*MPL(igrow)
       ! 'grow', tentatively:
       Call Grow(igrow,dx,dy,dz,newlen,newmas,t,ordgrw(igrow,ipl),ipl)
       ! add up tentative need for dry mass:
       sumgrw = sumgrw+newmas
       ! make a new, tentative record:
    End Do tiploop

    ! if shoot is considered, increase root mass by accumulated increment 'dmroot' or by potential growth:
    If (lCalloc) Then
       mroot(ipl) = mroot(ipl)+Min(dmroot,sumgrw)
    Else
       mroot(ipl) = mroot(ipl)+sumgrw
    Endif

    ! average soil strength experienced by growing root system:
    sAvg = sAvg/Real(ngrwnw)
    ! calculate growth factor from potential growth and available assimilates:
    If ((lCalloc).And.(sumgrw.Gt.0.0_dp)) Then
       grwfac = dmroot/sumgrw
    Else
       grwfac = 0.0_dp
    Endif

    ! reset 'dmroot':
    dmroot = 0.0_dp
    ALLOCATE(toosml(1:nrec(ipl)+maxgrw))
    toosml = .FALSE.

    ! calculate total branch lengths
    Do igrow = 1,ngrwnw
       If (lCalloc) Then
          ! adjust growth according to available assimilates:
          Call Adjust(igrow,newlen,ipl)
       Else
          newlen = seglen(irecsg(igrow,ipl),ipl)
          brlgth(igrow,ipl) = brlgth(igrow,ipl)+newlen 
          
          If (newlen.GE.1.E-05_dp) Then
             toosml(irecsg(igrow,ipl)) = .False.
          Else
             toosml(irecsg(igrow,ipl))=.True.
          Endif
       Endif

     IF ((.NOT.(stopgr(igrow))).AND.(.NOT.(toosml(irecsg(igrow,ipl)))).AND.(ordgrw(igrow,ipl).LT.norder(ipl))) THEN
          ! calculate branch spacing along new segment:
          Call Spacng(igrow,space,impfac,ipl)
          ! ...and establish new branching points:
          Call Establ(igrow,newlen,space,t,ipl)
       Endif
    End Do

    ! remove all new segments smaller than 1.E-3:
    Call Remove(nrecol,ipl)
    Do igrow = 1,ngrwnw
       If(geom.Eq.3 .And. nmat.Gt.1) Then
          Call Boxlim_Cylinder(igrow,ipl)
       Else
          Call Boxlim(igrow,ipl)
       End If
       ! remove all new segments smaller than 1.E-3:
    End Do

    ! remove branches that have stopped growth from list of growing branches
    ! and update the number of growing tips:
    Call Update(ngrwnw,ipl)
    ! calculate average solute concentration felt by the root system:
    If(ltoxi) Then
       Do iseg = 1,nrec(ipl)
          Call Neighb(xs(iseg,ipl),ys(iseg,ipl),zs(iseg,ipl),corner_r,imin,.FALSE.)
          cLoc = (Conc(corner_r(1))+Conc(corner_r(2))+Conc(corner_r(3))+Conc(corner_r(4))+&
               Conc(corner_r(5))+Conc(corner_r(6))+Conc(corner_r(7))+Conc(corner_r(8)))/8
          cAvg = cAvg+cLoc
       End Do
       cAvg = cAvg/Real(nrec(ipl))
    Else
       cAvg = 0._dp
    End If
    DEALLOCATE(toosml)

    ! continuous domain; determine all the transroot for the new tips
    If(continu .And. .Not.lDou) Then
       Do irecn = 1,nrec(ipl)
          Call roottransGrow(xs(irecn,ipl),ys(irecn,ipl),irecn,1,1)
       End Do
       Do igrow = 1,ngrow(ipl)
          Call tiptransGrow(xg(igrow,ipl),yg(igrow,ipl),igrow,1,1)  
       End Do
       transroot(nrec(ipl)+1:nrec(ipl)+ngrow(ipl),:,1,1)=transtip(1:ngrow(ipl),:,1,1)
    End If

    Return
  End Subroutine ROOT
  !********************************************************************************
  !> originates new axes
  Subroutine Newaxs(t,sumgrw,ngrwnw,ipl)
    Use typedef
    Use RootData
    Use ParamData
    Use tmctrl, Only: kaxemg
    Use DoussanMat, Only: transtip,nsub
    Use GridData, Only: continu,nex,ney,dxgrid,dygrid
    Use Environmental, Only: TemLoc, StrLoc
    Use MPIutils, Only: stop_program
    Use RootGrowthNeighb, Only: Neighb
    Implicit None

    Real(dp),Intent(in) ::t
    Real(dp), Intent(inout) :: sumgrw
    Real(dp) :: newmas,newlen,unimplen,impfac,sLoc  !,MPL
    Real(dp) ::  tloc,alpha,v
    Real(dp) ::  dx,dy,dz,dxini,dyini,dzini,dxstr,dystr,dzstr,drivex,drivey,drivez
    Real(dp), Dimension (3,3) :: condmatrix
    Integer(ap), Intent(inout) :: ngrwnw
    INTEGER(ap),INTENT(in)::ipl
    Integer(ap) :: imin,i,corner_r(8)
    !> \param sumgrw total root dry mass
    !> \param ngrwnw number of root tips (new growth)

    ! generate a total of 'nnewax' new axes:
    Do i = 1,nnewax(kaxemg,ipl)
       ! increment total axes number:
       naxes(ipl) = naxes(ipl)+1
       ! increment total branch number and total number of growing branches:
       nbr(ipl) = nbr(ipl)+1
       If (ngrwnw.Lt.maxgrw) Then
          ngrwnw = ngrwnw+1
       Else
          Call stop_program('Maximum number of growing tips -- PROGRAM TERMINATED.')
       Endif
       ! assign values for axis number, branching order and branch number:
       iaxis(ngrwnw,ipl)  = naxes(ipl)
       ordgrw(ngrwnw,ipl) = 1
       ibrgrw(ngrwnw,ipl) = nbr(ipl)
       ! as of now, the new branch has no estbl. brnch. pts. itself:
       nestbl(ngrwnw,ipl)=0
       ! also, length is zero and no pending brnch. pts. exist:
       brlgth(ngrwnw,ipl)  = 0.0_dp
       ovrtime(ngrwnw,ipl) = -1._dp
       ! no segment behind the tip:
       irecsg(ngrwnw,ipl) = 1
       ! emergence position:
       xg(ngrwnw,ipl) = xs(1,ipl)
       yg(ngrwnw,ipl) = ys(1,ipl)
       zg(ngrwnw,ipl) = zs(1,ipl)
       ! find surrounding grid points:
       If(continu) Then
          Call tiptransGrow(xg(ngrwnw,ipl),yg(ngrwnw,ipl),ngrwnw,1,1)
          nsub(nrec(ipl)+ngrwnw,1) = 1
          Call Neighb(xg(ngrwnw,ipl)+transtip(ngrwnw,1,1,1)*nex*dxgrid, &
               yg(ngrwnw,ipl)+transtip(ngrwnw,2,1,1)*ney*dygrid,zg(ngrwnw,ipl),corner_r,imin,ltwo_grids)
       Else
          Call Neighb(xg(ngrwnw,ipl),yg(ngrwnw,ipl),zg(ngrwnw,ipl),corner_r,imin,ltwo_grids)
       End If
       ! choose approach and get soil strength gradient components:
       If (l_conduc) Then
          Call Ssgcom2(corner_r,condmatrix,ipl)
       Else
          Call Ssgcom1(corner_r,1._dp,strsen(1,ipl),dxstr,dystr,dzstr)
       End If

       ! current local temperature:
       If (ltemp) Call TemLoc(corner_r,tLoc)



       ! calculate the direction of the first sgement of the new axis:
       Call RANDOM_Number(rand)
       alpha = rand*2._dp*pi
       Call Initial(1,alpha,inaxs(i,ipl),dxini,dyini,dzini)
       ! add up all components:
       dx = dxini
       dy = dyini
       dz = dzini

       ! calculate representative soil strength value felt by root tip:
       Call StrLoc(corner_r,sLoc)
       sAvg = sAvg+sLoc
       ! calculate length of the new segment from time left for growth:
       Call Uniera(0._dp,1,ibrgrw(ngrwnw,ipl),v,ipl)
       Call Length(v,t-tnewax(kaxemg,ipl),drivex,drivey,drivez,condmatrix,newlen,unimplen,impfac,corner_r,ipl)

       ! check if maximum length has been reached:
       If (newlen.Ge.brlmax(1,ipl))  newlen = brlmax(1,ipl)
       ! make sure new segment is not too small (prevent division by zero):
       newlen = Max(newlen,1.E-4_dp)
       ! calculate tentative segment mass:
       Call Maslen(ngrwnw,sLoc,ordgrw(ngrwnw,ipl),ipl)
       newmas = newlen*MPL(ngrwnw)
       ! 'grow', tentatively:
       Call Grow(ngrwnw,dx,dy,dz,newlen,newmas,tnewax(kaxemg,ipl),ordgrw(ngrwnw,ipl),ipl)
       ! add up tentative need for dry mass:
       sumgrw = sumgrw+newmas
    End Do

    Return
  End Subroutine Newaxs
  !****************************************************************************
  !> originates new sub-branches from branch 'igrow'
  !> branch needs a certain age to branch and there is a certain inter-branch time
  Subroutine Brdevl(t,sumgrw,igrow,ngrwnw,ipl)
    Use typedef
    Use RootData
    Use ParamData
    Use GridData, Only: continu,nex,ney,dxgrid,dygrid
    Use DoussanMat, Only: transtip,nsub
    Use Environmental, Only: StrLoc
    Use MPIutils, Only: stop_program
    Use RootGrowthNeighb, Only: Neighb
    Implicit None

    Integer(ap) :: irec,iprev,igrow,iest,ngrwnw,kest,ipprev
    INTEGER(ap), INTENT(in):: ipl
    Real(dp) :: newlen,unimplen,impfac,newmas,sumgrw,sLoc   !,MPL
    Real(dp) :: alpha,beta,gamma,delta,deltat,v,gammahdr,thetahdr
    Real(dp) :: t,test,drivex,drivey,drivez
    Real(dp) :: dx,dy,dz,t1,t2,x1,x2,y1,y2,z1,z2,x0,y0,z0,dxpr,dypr,dzpr
    Real(dp), Dimension (3,3) :: condmatrix
    Integer(ap) :: corner_r(8),imin
    !> \param t current time
    !> \param sumgrw total root dry mass
    !> \param igrow current tip or branch number
    !> \param ngrwnw number of root tips (new growth)

    kest = 0
    Do iest = 1,nestbl(igrow,ipl)
       test = timest(igrow,iest)
       If (lRootBox_growth) THEN
          dtbrch = lb(ordgrw(igrow,ipl),ipl)/(vch(1,1,ipl)) !lb is the same for all root orders
       ENDIF
       If ((t-test).Gt.dtbrch(ordgrw(igrow,ipl),ipl)) Then
          ! we have an established point 'ripe' for branching:
          kest = kest+1
          IF (kest.LE.1) THEN
          ! increment total branch number and total number of growing branches:
          nbr(ipl) = nbr(ipl)+1
          If (ngrwnw.Lt.maxgrw) Then
             ngrwnw = ngrwnw+1
          Else
             Call stop_program('Maximum number of growing tips -- PROGRAM TERMINATED.')
          Endif

          ! assign values for axis number, branching order and branch number:
          iaxis(ngrwnw,ipl)  = iaxis(igrow,ipl)
          ordgrw(ngrwnw,ipl) = ordgrw(igrow,ipl)+1
          ibrgrw(ngrwnw,ipl) = nbr(ipl)
          ! as of now, the new branch has no estbl. brnch. pts. itself:
          nestbl(ngrwnw,ipl) = 0
          ! also, length is zero and no pending brnch. pts. exist:
          brlgth(ngrwnw,ipl)  = 0.0_dp
          ovrtime(ngrwnw,ipl) = -1._dp
          ! find the segment where the branching occurs:
          irec=irecsg(igrow,ipl)
          If (timorg(irec,ipl).Le.test) Then
             ! have branching in the most recent segment:
             iprev = irec
             x2    = xg(igrow,ipl)
             y2    = yg(igrow,ipl)
             z2    = zg(igrow,ipl)
             t2    = t
          Else
             ! need to go back through the segment records:
1            If(timorg(irecpr(irec,ipl),ipl).Gt.test) Then
                irec = irecpr(irec,ipl)
                Goto 1
             Endif
             iprev = irecpr(irec,ipl)
             ipprev = irecpr(iprev,ipl)
             x2 = xs(irec,ipl)
             y2 = ys(irec,ipl)
             z2 = zs(irec,ipl)
             t2 = timorg(irec,ipl)
          Endif
          ! this is also the segment behind the tip:
          irecsg(ngrwnw,ipl) = irec
          x1             = xs(iprev,ipl)
          y1             = ys(iprev,ipl)
          z1             = zs(iprev,ipl)
          x0             = xs(ipprev,ipl)
          y0             = ys(ipprev,ipl)
          z0             = zs(ipprev,ipl)
          t1             = timorg(iprev,ipl)
          dx             = x2-x1
          dy             = y2-y1
          dz             = z2-z1
          dxpr           = x1-x0
          dypr           = y1-y0
          dzpr           = z1-z0
          deltat         = t2-t1


          ! --> small error, but origination point of new branch should be at an existing root node for Doussan!
          ! branch originates at x2,y2,z2
          xg(ngrwnw,ipl) = x2
          yg(ngrwnw,ipl) = y2
          zg(ngrwnw,ipl) = z2

          ! find surrounding grid points:
          If(continu) Then
             Call tiptransGrow(xg(ngrwnw,ipl),yg(ngrwnw,ipl),ngrwnw,1,1)
             nsub(nrec(ipl)+ngrwnw,1) = 1
             Call Neighb(xg(ngrwnw,ipl)+transtip(ngrwnw,1,1,1)*nex*dxgrid, &
                  yg(ngrwnw,ipl)+transtip(ngrwnw,2,1,1)*ney*dygrid,zg(ngrwnw,ipl),corner_r,imin,ltwo_grids)
          Else
             Call Neighb(xg(ngrwnw,ipl),yg(ngrwnw,ipl),zg(ngrwnw,ipl),corner_r,imin,ltwo_grids)
          End If

          ! calculate representative soil strength value felt by root tip:
          Call StrLoc(corner_r,sLoc)
          sAvg = sAvg+sLoc

          ! calculate initial heading components:
          Call RANDOM_Number(rand)
          gamma = rand*2._dp*pi
          delta = brnang(ordgrw(igrow,ipl),ipl)
 ! initial heading component in the direction of the lower soil potential 
          IF (lHydropattern) THEN ! Hydropaterning 
                !> calculate a unit vector for the driving force F
                Call Hsgcom1(corner_r,dxpr,dypr,dzpr,thetahdr,gammahdr)
                !gamma = gamma-pi+gammahdr
                gamma =  gammahdr + sin(pi/2._dp-hdrsens*pi/2._dp)* (gamma-gammahdr) ! *abs(sin(thetahdr))
                Call Angchg(dx,dy,dz,alpha,beta,gamma,delta)
          ELSE 
                Call Angchg(dx,dy,dz,alpha,beta,gamma,delta)
          ENDIF

          ! calculate length of the new segment from time left for growth:
          Call Uniera(0.0_dp,ordgrw(ngrwnw,ipl),ibrgrw(ngrwnw,ipl),v,ipl)
          If (lRootBox_growth) THEN
             dtbrch = lb(ordgrw(igrow,ipl),ipl)/v
          ENDIF
          Call Length(v,t-dtbrch(ordgrw(igrow,ipl),ipl)-test,drivex,drivey,drivez,condmatrix,newlen,unimplen,impfac,corner_r,ipl)

          !! make sure new segment is not too small (prevent divison by zero):
          newlen = Max(newlen,1.E-04_dp)

          ! calculate tentative segment mass:
          Call Maslen(ngrwnw,sLoc,ordgrw(ngrwnw,ipl),ipl)
          newmas = newlen*MPL(ngrwnw)

          ! 'grow', tentatively:
          Call Grow(ngrwnw,dx,dy,dz,newlen,newmas,test+dtbrch(ordgrw(igrow,ipl),ipl),ordgrw(ngrwnw,ipl),ipl)
          ! add up tentative need for dry mass:
          sumgrw = sumgrw+newmas
          ENDIF
       Endif
    End Do
    ! remove all developed points from the list of established branching points:
    nestbl(igrow,ipl) = nestbl(igrow,ipl)-kest
    Do iest = 1,nestbl(igrow,ipl)
       timest(igrow,iest) = timest(igrow,iest+kest)
    End Do

    Return
  End Subroutine Brdevl
  !**********************************************************************
  !>  defines branch age
  Subroutine Brnage(t,igrow,tempage,ipl)
    Use typedef
    Use RootData
    Implicit None

    Integer(ap) :: igrow,irec,ipl
    Real(dp), Intent(in) :: t
    Real(dp), Intent(out) :: tempage
    !> \param t current time
    !> \param igrow current branch or tip number
    !> \param tempage branch age (temporal parameter name)

    irec=irecsg(igrow,ipl)
    !> go back through the segment records to first record of branch 'igrow':
1   If (irecpr(irec,ipl).Ne.0) Then
       If (ordseg(irecpr(irec,ipl),ipl).Eq.ordseg(irec,ipl)) Then
          irec = irecpr(irec,ipl)
          Goto 1
       Endif
    Endif
    tempage = t-timorg(irec,ipl)
    If(tempage.Lt.0._dp) Then
       timorg(irec,ipl) = t
       tempage      = 0._dp
    End If

    Return
  End Subroutine Brnage
  !******************************************************************************
  !> unimpeded (= maximum) elongation rate
  Subroutine Uniera(tempage,iorder,ibr,v,ipl)
    Use Typedef
    Use RootData
    Implicit None

    Real(dp), Intent(in) ::  tempage
    Real(dp), Intent(out) :: v
    Integer(ap), Intent(in) :: iorder,ibr,ipl
    Integer(ap):: ivch, iorder_mod,br,countBigLat
    Real(dp) :: brlmaxuse
    !> \param tempage branch age (temporal parameter name)
    !> \param iorder order of current brach
    !> \param v growth velocity
	
    ! calculate unimpeded elongation rate as a function of order, age:
        IF (lSomma_growth) Then
           !! Ax Long and small lateral
           !! Order1 -> Iorder=1 Order2Long -> Iorder=2  Order2Small ->Iorder=3 Order3 ->Iorder=4
           If (iorder.eq.1) Then
              iorder_mod=1
           Elseif (iorder.eq.2) Then
              iorder_mod=3 !By default small lateral (if big lat it will be changed in dowhile)
              br=1
              countBigLat=0
              Do while ((countBigLat .LT. nBigLat(ipl)) .AND. (br .LE. nbr(ipl)))
                 if (ordgrw(br,ipl).EQ.2) Then
                      if (br .EQ. ibr) iorder_mod=2
                 countBigLat=countBigLat+1
                 endif
              br=br+1
              Enddo
           Elseif (iorder.eq.3) Then
              iorder_mod=4
           Endif
           !!
           If (tempage.Ge.agevch(iorder_mod,nvch(iorder_mod,ipl),ipl)) Then
                v = vch(iorder_mod,nvch(iorder_mod,ipl),ipl)
           Else
                ivch = nvch(iorder_mod,ipl)
1               ivch = ivch-1
                If ((tempage.Lt.agevch(iorder_mod,ivch,ipl)).And.(ivch.Gt.1)) Goto 1
                    v = vch(iorder_mod,ivch,ipl)+(tempage-agevch(iorder_mod,ivch,ipl))/(agevch(iorder_mod,ivch+1,ipl)-&
                    agevch(iorder_mod,ivch,ipl))*(vch(iorder_mod,ivch+1,ipl)-vch(iorder_mod,ivch,ipl))
           Endif
        ELSEIF (lRootBox_growth) Then
           !! Ax Long and small lateral
           !! Order1 -> Iorder=1 Order2Long -> Iorder=2  Order2Small ->Iorder=3 Order3 ->Iorder=4
           If (iorder.eq.1) Then
              iorder_mod=1
           Elseif (iorder.eq.2) Then
              iorder_mod=3 !By default small lateral (if big lat it will be changed in dowhile)
              br=1
              countBigLat=0
              Do while ((countBigLat .LT. nBigLat(ipl)) .AND. (br .LE. nbr(ipl)))
                 if (ordgrw(br,ipl).EQ.2) Then
                      if (br .EQ. ibr) iorder_mod=2
                 countBigLat=countBigLat+1
                 endif
              br=br+1
              Enddo
           Elseif (iorder.eq.3) Then
              iorder_mod=4
           Endif
           !!
           IF (brlmax(iorder_mod,ipl).LE.2) Then 
               brlmaxuse = 10 !if brlmax is too small, it leads to errors when continuous soil domain is used 
           Else
               brlmaxuse = brlmax(iorder_mod,ipl)
           Endif        
           v = (vch(iorder_mod,1,ipl)*exp(-(vch(iorder_mod,1,ipl)*tempage/brlmaxuse)))
        ENDIF

    Return
  End Subroutine Uniera
  !******************************************************************************
  !> check whether the root is in a macropore or not
  Subroutine CheckMacro(corner_r,lMacro)
    Use typedef
    Use StrData
    Implicit None

    Logical, Intent(out) :: lMacro
    Integer(ap), Intent(in) ::  corner_r(8)
    Integer(dp) ::  i
    !> \param lMacro
    !> \param corner soil nodes surrounding growing tip

    ! check
    lMacro=.False.
    Do i=1,8
       If (s(corner_r(i)).Le.1.E-1) Then
          lMacro=.True.
          Exit
       Endif
    End Do

    Return
  End Subroutine CheckMacro
  !******************************************************************************
  !> length of new segment
  Subroutine Length(v,tgrow,drivex,drivey,drivez,condmatrix,newlen,unimplen,impfac,corner_r,ipl)
    Use typedef
    Use TempData
    Use ConData
    Use StrData
    Use RootData, Only: ltemp,ltoxi,l_conduc,l_overburden,lSoilStrength
    Implicit None
    Logical :: lMacro
    Real(dp), Intent(in) ::  v,tgrow, drivex, drivey, drivez
    Real(dp), Dimension (3,3) :: condmatrix
    Real(dp) :: factord,drivexun, driveyun, drivezun
    Real(dp) :: kmidx, kmidy, kmidz,seff
    Real(dp) :: cncimp,temimp,strimp
    Real(dp), Intent(out) :: newlen,unimplen, impfac
    Integer(ap), Intent(in) ::  corner_r(8)
    Integer(ap) :: ipl

    !> \param v growth velocity
    !> \param tgrow length of current growth time step
    !> \param newlen length of new segment
    !> \param corner soil nodes surrounding growing tip

    ! calculate length of new segment, taking location into account, approach different for soil strength/conductance approach:

    If (l_conduc) Then
    
       !compute condmatrix
       CALL Ssgcom2(corner_r,condmatrix,ipl)
    
       !> calculate a unit vector for the driving force F
       factord = (Sqrt(drivex**2+drivey**2+drivez**2))

       If ((drivex /= drivex).OR.(drivey /= drivey).OR.(drivez /= drivez)) Then
            drivexun = 1
            driveyun = 1
            drivezun = 1
       else
          drivexun = drivex/factord
          driveyun = drivey/factord
          drivezun = drivez/factord
       Endif


       kmidx = (condmatrix(1,1)*drivexun+condmatrix(2,1)*driveyun+condmatrix(3,1)*drivezun)
       kmidy = (condmatrix(1,2)*drivexun+condmatrix(2,2)*driveyun+condmatrix(3,2)*drivezun)
       kmidz = (condmatrix(1,3)*drivexun+condmatrix(2,3)*driveyun+condmatrix(3,3)*drivezun)

       seff = 1/(Sqrt(kmidx**2+kmidy**2+kmidz**2))
       strimp = 1-(seff/ssmax)**.5 !Landl
    Else
       

       Call CheckMacro(corner_r,lMacro) !> root length shall be max length in macropore
       If (l_overburden) Then
           If (lMacro) Then
               strimp = Max(imps(corner_r(1)),imps(corner_r(2)),imps(corner_r(3)),imps(corner_r(4)),&
               imps(corner_r(5)),imps(corner_r(6)),imps(corner_r(7)),imps(corner_r(8)))
           Endif
           
       Else if (lSoilStrength) Then
           strimp = (imps(corner_r(1))+imps(corner_r(2))+imps(corner_r(3))+imps(corner_r(4))+&
           imps(corner_r(5))+imps(corner_r(6))+imps(corner_r(7))+imps(corner_r(8)))/8._dp
       Else
           strimp = 1._dp
       Endif
    Endif



    If (ltemp) Then
       temimp = (impt(corner_r(1))+impt(corner_r(2))+impt(corner_r(3))+impt(corner_r(4))+&
            impt(corner_r(5))+impt(corner_r(6))+impt(corner_r(7))+impt(corner_r(8)))/8._dp
    Else
       temimp = 1._dp
    Endif
    If (ltoxi) Then
       cncimp = (impc(corner_r(1))+impc(corner_r(2))+impc(corner_r(3))+impc(corner_r(4))+&
            impc(corner_r(5))+impc(corner_r(6))+impc(corner_r(7))+impc(corner_r(8)))/8._dp
    Else
       cncimp = 1._dp
    Endif

    impfac = (cncimp*strimp*temimp)
    newlen = Dble((cncimp*strimp*temimp)*v*tgrow)
    unimplen = Dble(v*tgrow)

    Return
  End Subroutine Length
  !******************************************************************************
  !> soil strength gradient components of new heading vector
  Subroutine Ssgcom1(corner_r,stndrd,senstr,dxstr,dystr,dzstr)
    Use typedef
    Use ParamData, Only: pi
    Use GridData
    Use StrData
    Implicit None

    Real(dp) :: senstr,smax
    Real(dp), Intent(in) :: stndrd
    Real(dp), Intent(out) :: dxstr,dystr,dzstr
    Integer(ap), Intent(in) :: corner_r(8)

    !> \param corner soil nodes surrounding growing tip
    !> \param stndrd length of previous segment
    !> \param senstr sensitivity of growth direction to soil strength
    !> \param dxstr x component of new segment with respect to the soil strength gradient
    !> \param dystr y component of new segment with respect to the soil strength gradient
    !> \param dzstr z component of new segment with respect to the soil strength gradient

    smax = Max(s(corner_r(1)),s(corner_r(2)),s(corner_r(3)),s(corner_r(4)),&
         s(corner_r(5)),s(corner_r(6)),s(corner_r(7)),s(corner_r(8)))


    !> calculate, normalize, and weight soil strength gradient components:
    dxstr = (s(corner_r(1))+s(corner_r(3))+s(corner_r(5))+s(corner_r(7))-&
         s(corner_r(2))-s(corner_r(4))-s(corner_r(6))-s(corner_r(8)))/(4._dp*smax*dxgrid)
    dystr = (s(corner_r(1))+s(corner_r(2))+s(corner_r(5))+s(corner_r(6))-&
         s(corner_r(3))-s(corner_r(4))-s(corner_r(7))-s(corner_r(8)))/(4._dp*smax*dygrid)
    dzstr = -(s(corner_r(1))+s(corner_r(2))+s(corner_r(3))+s(corner_r(4))-&
         s(corner_r(5))-s(corner_r(6))-s(corner_r(7))-s(corner_r(8)))/(4._dp*smax*dzgrid)


    If ((Abs(dxstr).Gt.1.E-10_dp).Or.(Abs(dystr).Gt.1.E-10_dp).Or.(Abs(dzstr).Gt.1.E-10_dp)) Then
       dxstr = dxstr*senstr*stndrd
       dystr = dystr*senstr*stndrd
       dzstr = dzstr*senstr*stndrd
    Endif

    Return
  End Subroutine Ssgcom1
  !******************************************************************************
  !> soil water potential gradient components for a new heading vector
  Subroutine Hsgcom1(corner_r,dxpr,dypr,dzpr,thetahdr,gammahdr)
    Use typedef
    Use ParamData, Only: pi
    Use GridData
    Use RootData
    Use SolData, Only: hnew
    Implicit None

    Real(dp), Intent(out) :: thetahdr,gammahdr
    Integer(ap), Intent(in) :: corner_r(8)
    Real(dp), Intent(in) :: dxpr,dypr,dzpr
    Real(dp) :: dxhdr,dyhdr,dzhdr,normhdr,seglength,obliquity,betahdr,deltahdr,alpha_gen
    Real(dp) :: lamdahdr,sin_lamdahdr,cos_lamdahdr

    !> \param corner soil nodes surrounding growing tip
    !> \param dxpr x component of old segment
    !> \param dypr y component of old segment
    !> \param dzpr z component of old segment
    !> \param dxhdr x component of new segment with respect to the soil water potential gradient
    !> \param dyhdr y component of new segment with respect to the soil water potential gradient
    !> \param dzhdr z component of new segment with respect to the soil water potential gradient

    !> calculate, soil water potential gradient components:
    dxhdr = (hnew(corner_r(1))+hnew(corner_r(3))+hnew(corner_r(5))+hnew(corner_r(7))-&
         hnew(corner_r(2))-hnew(corner_r(4))-hnew(corner_r(6))-hnew(corner_r(8)))/(4._dp)
    dyhdr = (hnew(corner_r(1))+hnew(corner_r(2))+hnew(corner_r(5))+hnew(corner_r(6))-&
         hnew(corner_r(3))-hnew(corner_r(4))-hnew(corner_r(7))-hnew(corner_r(8)))/(4._dp)
    dzhdr = -(hnew(corner_r(1))+hnew(corner_r(2))+hnew(corner_r(3))+hnew(corner_r(4))-&
         hnew(corner_r(5))-hnew(corner_r(6))-hnew(corner_r(7))-hnew(corner_r(8)))/(4._dp)
    
    normhdr = sqrt(dxhdr*dxhdr+dyhdr*dyhdr+dzhdr*dzhdr)
    deltahdr = asin(dzhdr/normhdr) ! the angle between the h gradient and the xy-plane
    seglength = sqrt(dxpr*dxpr+dypr*dypr+dzpr*dzpr)

 ! beta is the angle between the h gradient and the root ecliptic
    thetahdr = acos((dxpr*dxhdr+dypr*dyhdr+dzpr*dzhdr)/(seglength*normhdr))
    if(thetahdr .LT. pi/2._dp) then
       betahdr = pi/2_dp-thetahdr
    else
       betahdr = -(thetahdr-pi/2_dp)
    endif

    obliquity = acos((dzpr)/(seglength))

 ! if the root is growing in the same direction as the differential of water potential gradient then gammahdr should be random across [0-2*pi]
    IF ((thetahdr .LT. 1.E-20_dp) .OR. (thetahdr .GE. 2_dp*pi)) THEN  ! Vertical gradient
       Call RANDOM_Number(rand)
       gammahdr = 2._dp*pi*rand
    ELSE   ! in spherical coordinate if we know dxhdr and theta then we can calculate the second angl over the perpendicular plane
       sin_lamdahdr = (sin(deltahdr)-sin(betahdr)*cos(obliquity))/(cos(betahdr)*sin(obliquity))
       if (sin_lamdahdr .GT. 1._dp) then
           print*, '+',sin_lamdahdr
           sin_lamdahdr = 1._dp
       else if (sin_lamdahdr .LT. -1._dp) then
           print*, '-',sin_lamdahdr
           sin_lamdahdr = -1._dp
       endif
       lamdahdr = asin(sin_lamdahdr)
       alpha_gen = atan2((cos(betahdr)*sin(lamdahdr)*cos(obliquity)-sin(betahdr)*sin(obliquity)), cos(betahdr)*cos(lamdahdr))
       if(alpha_gen .LT. 1.E-20_dp) alpha_gen = alpha_gen + 2._dp*pi
       cos_lamdahdr = cos(alpha_gen)*cos(deltahdr)/cos(betahdr)

       if((alpha_gen .GT. pi/2._dp) .OR. (alpha_gen .LT. 3._dp*pi/2._dp)) then
            gammahdr = -lamdahdr + ((3._dp*pi)/2._dp)
       else
            gammahdr = lamdahdr + (pi/2._dp)
       endif
       IF (dzpr .LT. 0._dp) then
           gammahdr = 2._dp *pi - gammahdr  ! reverse orientation
       ENDIF

    ENDIF
    Return
  End Subroutine Hsgcom1
  
  !******************************************************************************
  !> soil strength gradient components of new heading vector (K-like approach)
  Subroutine Ssgcom2(corner_r,condmatrix,ipl)
    Use typedef
    Use StrData
    Use GridData
    Use RootData
    Use ParamData, Only: pi
    Implicit None

    Real(dp) :: smax,kxx_lin,kyy_lin,kzz_lin
    Real(dp) :: kxx_skew_x,kyy_skew_x,kzz_skew_x
    Real(dp) :: kxx_skew_y,kyy_skew_y,kzz_skew_y
    Real(dp) :: kxx_skew_z,kyy_skew_z,kzz_skew_z
    Real(dp) :: ortho1,ortho2,skew_x1,skew_x2,skew_y1,skew_y2,skew_z1,skew_z2
    Real(dp) :: kxx,kxy,kxz,kyx,kyy,kyz,kzx,kzy,kzz
    Real(dp) :: kx1,ky1,kz1,kx2,ky2,kz2,phi
    Integer(ap) :: i
    Integer(ap), Intent(in) :: corner_r(8),ipl
    Real(dp), Dimension (3,3) :: rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2
    Real(dp), Dimension (3,3) :: condmatrix
    Real(dp), Dimension (3,3) :: skewmatrix_in_x, skewmatrix_in_y, skewmatrix_in_z, linmatrix
    Real(dp), Dimension (3,3) :: skewmatrix_mi
    Real(dp), Dimension (3,3) :: skewmatrix_x, skewmatrix_y, skewmatrix_z,skewmatrix
    Real(dp), Dimension (3,4) :: inmatrix
    !> \param smax = maximum soil strength of one grid cell node
    !> \param kx1,ky1,kz1,kx2,ky2,kz2 = average conductance of each grid cell plane
    !> \param kxx_lin, kyy_lin, kzz_lin = average conductance in direction of the main axes of the Cartesian coordinate system
    !> \param kxx_skew_x,y,z ,kyy_skew_x,y,z, kzz_skew_x,y,z = average conductances for the local coordinate systems rotated by 45° aroudn the x,y or z-axis
    !> \param ortho1,ortho2,skew_x1,skew_x2,skew_y1,skew_y2,skew_z1,skew_z2 = highest & 2nd highest conductance value of each condcutance tensor in the 4 directions (alignes with the C coord system, 45° rotated in x,y,z direction)
    !> \param kxx,kxy,kxz,kyx,kyy,kyz,kzx,kzy,kzz = final values of the conductance tensor
    !> \param phi = rotation angle
    !> \param i = number of soil node surrounding growing tip
    !> \param corner = soil nodes surrounding growing tip
    !> \param skewmatrix_in_x, skewmatrix_in_y, skewmatrix_in_z = conductance matrices int he local coordinate systems
    !> \param rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2 = rotation matrix used to map the local coordinate systems back to the Cartesian coordinate system
    !> \param skewmatrix_mi = intermediate stage of the rotated matrices
    !> \param skewmatrix_x, skewmatrix_y, skewmatrix_z, linmatrix = conductance matrices in all 4 directions
    !> \param inmatrix = matrixdummy used to sort the individual conductance tensors
    !> \param skewmatrix = randomly chosen conductance matrix from skewmatrix_x, skewmatrix_y, skewmatrix_z
    !> \param condmatrix = conductance tensor

    Do i=1,8
       condu(corner_r(i)) = 1/s(corner_r(i))
       If (condu(corner_r(i)).Ge.1.E1_dp) Then
          condu(corner_r(i)) = condMP(ipl)
          s(corner_r(i)) = 1/condMP(ipl)
       Endif
    Enddo

    smax = Max(s(corner_r(1)),s(corner_r(2)),s(corner_r(3)),s(corner_r(4)),&
         s(corner_r(5)),s(corner_r(6)),s(corner_r(7)),s(corner_r(8)))



    !> conductance1: axes of anisotropy aligned with the Cartesian Coordinate systm

    !> calculate averaged & normalized soil conductance per grid plane (harmonic mean of the resistance in each corner point)

    kx1=(condu(corner_r(1))+condu(corner_r(3))+condu(corner_r(5))+condu(corner_r(7)))/(4._dp)
    kx2=(condu(corner_r(2))+condu(corner_r(4))+condu(corner_r(6))+condu(corner_r(8)))/(4._dp)
    ky1=(condu(corner_r(1))+condu(corner_r(2))+condu(corner_r(5))+condu(corner_r(6)))/(4._dp)
    ky2=(condu(corner_r(3))+condu(corner_r(4))+condu(corner_r(7))+condu(corner_r(8)))/(4._dp)
    kz1=(condu(corner_r(1))+condu(corner_r(2))+condu(corner_r(3))+condu(corner_r(4)))/(4._dp)
    kz2=(condu(corner_r(5))+condu(corner_r(6))+condu(corner_r(7))+condu(corner_r(8)))/(4._dp)


    kxx_lin=((2._dp)/(1._dp/kx1+1._dp/kx2))
    kyy_lin=((2._dp)/(1._dp/ky1+1._dp/ky2))
    kzz_lin=((2._dp)/(1._dp/kz1+1._dp/kz2))


    If ((Abs(kxx_lin).Gt.1.E-10_dp).Or.(Abs(kyy_lin).Gt.1.E-10_dp).Or.(Abs(kzz_lin).Gt.1.E-10_dp)) Then
       kxx_lin= Nint((kxx_lin) * 1.E4)/1.E4
       kyy_lin= Nint((kyy_lin) * 1.E4)/1.E4
       kzz_lin= Nint((kzz_lin) * 1.E4)/1.E4
    Endif


    linmatrix = (Reshape((/kxx_lin,0._dp,0._dp,0._dp,kyy_lin,&
         0._dp,0._dp,0._dp,kzz_lin/),Shape(linmatrix)))


    !> conductance2: x is aligned with the Cartesian x-axis, y&z in 45-degree angle

    !> calculate averaged & normalized conductance per cube section (harmonic mean of the conductance in each corner point)

    kx1=(condu(corner_r(1))+condu(corner_r(3))+condu(corner_r(5))+condu(corner_r(7)))/(4._dp)
    kx2=(condu(corner_r(2))+condu(corner_r(4))+condu(corner_r(6))+condu(corner_r(8)))/(4._dp)
    ky1= (condu(corner_r(3))+condu(corner_r(4))+condu(corner_r(1))/2+condu(corner_r(2))/2+condu(corner_r(7))/2+condu(corner_r(8))/2)/(4._dp)
    ky2= (condu(corner_r(5))+condu(corner_r(6))+condu(corner_r(1))/2+condu(corner_r(2))/2+condu(corner_r(7))/2+condu(corner_r(8))/2)/(4._dp)
    kz2= (condu(corner_r(1))+condu(corner_r(2))+condu(corner_r(3))/2+condu(corner_r(4))/2+condu(corner_r(5))/2+condu(corner_r(6))/2)/(4._dp)
    kz2= (condu(corner_r(7))+condu(corner_r(8))+condu(corner_r(3))/2+condu(corner_r(4))/2+condu(corner_r(5))/2+condu(corner_r(6))/2)/(4._dp)

    !> calculate the conductivities along the diagonal principal axes
    kxx_skew_x=((2._dp)/(1._dp/kx1+1._dp/kx2))
    kyy_skew_x=((2._dp)/(1._dp/ky1+1._dp/ky2))
    kzz_skew_x=((2._dp)/(1._dp/kz1+1._dp/kz2))


    !> normalize the 'soil conductivity' with the total conductance (from the axis aligned conductance tensor)
    If ((Abs(kxx_skew_x).Gt.1.E-10_dp).Or.(Abs(kyy_skew_x).Gt.1.E-10_dp).Or.(Abs(kzz_skew_x).Gt.1.E-10_dp)) Then
       kxx_skew_x= Nint((kxx_skew_x)  * 1.E4)/1.E4
       kyy_skew_x= Nint((kyy_skew_x)  * 1.E4)/1.E4
       kzz_skew_x= Nint((kzz_skew_x)  * 1.E4)/1.E4
    Endif


    !>create a conductivity matrix

    skewmatrix_in_x = (Reshape((/kxx_skew_x,0._dp,0._dp,0._dp,kyy_skew_x,&
         0._dp,0._dp,0._dp,kzz_skew_x/),Shape(skewmatrix_in_x)))


    Call RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)

    skewmatrix_mi = Matmul(rotmatrix_x1,skewmatrix_in_x)
    skewmatrix_x = Matmul(skewmatrix_mi,rotmatrix_x2)

    !> conductance3: y is aligned with the Cartesian y-axis, x&z in 45-degree angle


    !> calculate averaged & normalized conductance per cube section (harmonic mean of the conductance in each corner point)

    kx1= (condu(corner_r(1))+condu(corner_r(3))+condu(corner_r(4))/2+condu(corner_r(2))/2+condu(corner_r(5))/2+condu(corner_r(7))/2)/(4._dp)
    kx2= (condu(corner_r(6))+condu(corner_r(8))+condu(corner_r(4))/2+condu(corner_r(2))/2+condu(corner_r(5))/2+condu(corner_r(7))/2)/(4._dp)
    ky1=(condu(corner_r(1))+condu(corner_r(2))+condu(corner_r(5))+condu(corner_r(6)))/(4._dp)
    ky2=(condu(corner_r(3))+condu(corner_r(4))+condu(corner_r(7))+condu(corner_r(8)))/(4._dp)
    kz1= (condu(corner_r(2))+condu(corner_r(4))+condu(corner_r(1))/2+condu(corner_r(3))/2+condu(corner_r(6))/2+condu(corner_r(8))/2)/(4._dp)
    kz2= (condu(corner_r(5))+condu(corner_r(7))+condu(corner_r(1))/2+condu(corner_r(3))/2+condu(corner_r(6))/2+condu(corner_r(8))/2)/(4._dp)

    !> calculate the conductivities along the diagonal principal axes
    kxx_skew_y=((2._dp)/(1._dp/kx1+1._dp/kx2))
    kyy_skew_y=((2._dp)/(1._dp/ky1+1._dp/ky2))
    kzz_skew_y=((2._dp)/(1._dp/kz1+1._dp/kz2))


    !> normalize the 'soil conductivity' with the total strength (from the orthogonal conductivity tensor)
    If ((Abs(kxx_skew_y).Gt.1.E-10_dp).Or.(Abs(kyy_skew_y).Gt.1.E-10_dp).Or.(Abs(kzz_skew_y).Gt.1.E-10_dp)) Then
       kxx_skew_y= Nint((kxx_skew_y)  * 1.E4)/1.E4
       kyy_skew_y= Nint((kyy_skew_y)  * 1.E4)/1.E4
       kzz_skew_y= Nint((kzz_skew_y)  * 1.E4)/1.E4
    Endif


    !>create a conductivity matrix

    skewmatrix_in_y = (Reshape((/kxx_skew_y,0._dp,0._dp,0._dp,kyy_skew_y,&
         0._dp,0._dp,0._dp,kzz_skew_y/),Shape(skewmatrix_in_y)))


    Call RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)

    skewmatrix_mi = Matmul(rotmatrix_y1,skewmatrix_in_y)
    skewmatrix_y = Matmul(skewmatrix_mi,rotmatrix_y2)

    !> conductance3: z is aligned with the Cartesian z-axis, x&y in 45-degree angle


    !> calculate averaged & normalized conductance per cube section (harmonic mean of the conductance in each corner point)

    kx1= (condu(corner_r(1))+condu(corner_r(5))+condu(corner_r(2))/2+condu(corner_r(6))/2+condu(corner_r(3))/2+condu(corner_r(7))/2)/(4._dp)
    kx2= (condu(corner_r(4))+condu(corner_r(8))+condu(corner_r(2))/2+condu(corner_r(6))/2+condu(corner_r(3))/2+condu(corner_r(7))/2)/(4._dp)
    ky1= (condu(corner_r(2))+condu(corner_r(6))+condu(corner_r(1))/2+condu(corner_r(5))/2+condu(corner_r(4))/2+condu(corner_r(8))/2)/(4._dp)
    ky2= (condu(corner_r(3))+condu(corner_r(7))+condu(corner_r(1))/2+condu(corner_r(5))/2+condu(corner_r(4))/2+condu(corner_r(8))/2)/(4._dp)
    kz1=(condu(corner_r(1))+condu(corner_r(2))+condu(corner_r(3))+condu(corner_r(4)))/(4._dp)
    kz2=(condu(corner_r(5))+condu(corner_r(6))+condu(corner_r(7))+condu(corner_r(8)))/(4._dp)

    !> calculate the conductivities along the diagonal principal axes
    kxx_skew_z=((2._dp)/(1._dp/kx1+1._dp/kx2))
    kyy_skew_z=((2._dp)/(1._dp/ky1+1._dp/ky2))
    kzz_skew_z=((2._dp)/(1._dp/kz1+1._dp/kz2))


    !> normalize the 'soil conductivity' with the total strength (from the orthogonal conductivity tensor)
    If ((Abs(kxx_skew_z).Gt.1.E-10_dp).Or.(Abs(kyy_skew_z).Gt.1.E-10_dp).Or.(Abs(kzz_skew_z).Gt.1.E-10_dp)) Then
       kxx_skew_z= Nint((kxx_skew_z)  * 1.E4)/1.E4
       kyy_skew_z= Nint((kyy_skew_z)  * 1.E4)/1.E4
       kzz_skew_z= Nint((kzz_skew_z)  * 1.E4)/1.E4
    Endif


    skewmatrix_in_z = (Reshape((/kxx_skew_z,0._dp,0._dp,0._dp,kyy_skew_z,&
         0._dp,0._dp,0._dp,kzz_skew_z/),Shape(skewmatrix_in_z)))


    Call RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)

    skewmatrix_mi = Matmul(rotmatrix_z1,skewmatrix_in_z)
    skewmatrix_z = Matmul(skewmatrix_mi,rotmatrix_z2)


    !> Find the maximum values of all conductivity tensors
    inmatrix = (Reshape((/kxx_lin,kyy_lin,kzz_lin,kxx_skew_x,kyy_skew_x,kzz_skew_x,&
         kxx_skew_y,kyy_skew_y,kzz_skew_y, kxx_skew_z,&
         kyy_skew_z, kzz_skew_z/),Shape(inmatrix)))

    Call Sortmatrix(inmatrix,ortho1,skew_x1,skew_y1,skew_z1,ortho2,skew_x2,skew_y2,skew_z2)

    !>Decision Rule
    !> 1) take the tensor with the maximum cond in one of the three directions
    If (skew_x1.Gt.Max(ortho1,skew_y1,skew_z1))  Then
       condmatrix = skewmatrix_x
    Elseif (skew_y1.Gt.Max(ortho1,skew_x1,skew_z1)) Then
       condmatrix = skewmatrix_y
    Elseif (skew_z1.Gt.Max(ortho1,skew_x1,skew_y1))Then
       condmatrix = skewmatrix_z
       !> 2) if two or more of these values are the same, take the tensor with the highest second largest value
    Elseif (skew_x2.Gt.Max(ortho2,skew_y2,skew_z2)) Then
       condmatrix = skewmatrix_x
    Elseif (skew_y2.Gt.Max(ortho2,skew_x2,skew_z2)) Then
       condmatrix = skewmatrix_y
    Elseif (skew_z2.Gt.Max(ortho2,skew_x2,skew_y2))Then
       condmatrix = skewmatrix_z
       !> 3) if two of the skew conductivity tensors are the same AND if their maximum values are higher then the maximum value of the orthogonal tensor
    Elseif ((ortho1.Lt.Max(skew_x1,skew_y1,skew_z1)).And.&
         (Max(skew_x1,skew_y1,skew_z1).Ne.Min(skew_x1,skew_y1,skew_z1))) Then
       Call CondDecision1(skew_x1,skew_y1,skew_z1,skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
       condmatrix = skewmatrix
       !> 4) if the three skew conductivity tensors are the same and if their maximum value is hogher than that of the orthogonal tensor
    Elseif (ortho1.Lt.Max(skew_x1,skew_y1,skew_z1)) Then
       Call CondDecision2(skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
       condmatrix = skewmatrix
    Else
       !> in all other cases: take the orthogonal conductivity tensor!
       condmatrix = linmatrix
    Endif

    !> take the k-values out of the conductivity matrix

    kxx = condmatrix(1,1)
    kxy = condmatrix(2,1)
    kxz = condmatrix(3,1)
    kyx = condmatrix(1,2)
    kyy = condmatrix(2,2)
    kyz = condmatrix(3,2)
    kzx = condmatrix(1,3)
    kzy = condmatrix(2,3)
    kzz = condmatrix(3,3)


    Return
  End Subroutine Ssgcom2

  !******************************************************************************
  Subroutine Sortmatrix(inmatrix,ortho1,skew_x1,skew_y1,skew_z1,ortho2,skew_x2,skew_y2,skew_z2)
    Use typedef
    Implicit None

    Real(dp), Dimension (3,4), Intent(in) :: inmatrix
    Real(dp), Dimension (3,4) :: smatrix
    Real(dp), Intent(out) :: ortho1, skew_x1, skew_y1, skew_z1
    Real(dp), Intent(out) :: ortho2, skew_x2, skew_y2, skew_z2
    Real(dp) :: dummy
    Integer(ap) :: hh,ii,jj
    !> \param inmatrix = matrixdummy used to sort the individual conductance tensors
    !> \param smatrix, dummy = dummy matrices
    !> \param hh,ii,jj = indices variables
    !> \param ortho1, skew_x1, skew_y1, skew_z1, ortho2, skew_x2, skew_y2, skew_z2  = highest and 2nd highest conductanec values of each conductance tensor

    smatrix=inmatrix

    Do hh=1,4
       Do ii=1,3
          Do jj=ii,3
             If (smatrix(ii,hh).Le.smatrix(jj,hh)) Then
                dummy = smatrix(jj,hh)
                smatrix(jj,hh) = smatrix(ii,hh)
                smatrix(ii,hh) = dummy
             Endif
          Enddo
       Enddo
    Enddo


    ortho1 =  smatrix(1,1)
    skew_x1 = smatrix(1,2)
    skew_y1 = smatrix(1,3)
    skew_z1 = smatrix(1,4)
    ortho2 =  smatrix(2,1)
    skew_x2 = smatrix(2,2)
    skew_y2 = smatrix(2,3)
    skew_z2 = smatrix(2,4)

  End Subroutine Sortmatrix
  !******************************************************************************
  !> Random decision which of two equivalent conductivities should be chosen
  Subroutine CondDecision1(skew_x1,skew_y1,skew_z1,skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
    Use typedef
    Use RootData, Only: rand
    Implicit None

    Real(dp) :: u,j
    Real(dp), Intent(in):: skew_x1,skew_y1,skew_z1
    Real(dp), Dimension (3,3), Intent(in) :: skewmatrix_x, skewmatrix_y, skewmatrix_z
    Real(dp), Dimension (3,3), Intent(out) :: skewmatrix
    !> \param skew_x1,skew_y1,skew_z1 = highest value of each conductance tensor in the local, by 45° rotated, coordinate system
    !> \param skewmatrix_x, skewmatrix_y, skewmatrix_z = rotated conductance tensors
    !> \param u,j = indices variables
    !> \param skewmatrix = randomly chosen conductance matrix from skewmatrix_x, skewmatrix_y, skewmatrix_z

    Call RANDOM_Number(rand)
    u = rand
    j = Floor(2*u)

    If (skew_x1.Eq.skew_y1) Then
       If (j.Eq.0) Then
          skewmatrix = skewmatrix_x
       Else
          skewmatrix = skewmatrix_y
       Endif
    Elseif (skew_x1.Eq.skew_z1) Then
       If (j.Eq.0) Then
          skewmatrix = skewmatrix_x
       Else
          skewmatrix = skewmatrix_z
       Endif
    Else
       If (j.Eq.0) Then
          skewmatrix = skewmatrix_y
       Else
          skewmatrix = skewmatrix_z
       Endif
    Endif


  End Subroutine CondDecision1
  !******************************************************************************
  !> Random decision which of three equivalent conductivity tensors shall be chosen
  Subroutine CondDecision2(skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
    Use typedef
    Use RootData, Only: rand
    Implicit None

    Real(dp) :: u,j
    Real(dp), Dimension (3,3), Intent(in) :: skewmatrix_x, skewmatrix_y, skewmatrix_z
    Real(dp), Dimension (3,3), Intent(out) :: skewmatrix
    !> \param skewmatrix_x, skewmatrix_y, skewmatrix_z = rotated conductance tensors
    !> \param u,j = indices variables
    !> \param skewmatrix = randomly chosen conductance matrix from skewmatrix_x, skewmatrix_y, skewmatrix_z


    Call RANDOM_Number(rand)
    u = rand
    j = Floor(3*u)

    If (j.Eq.0) Then
       skewmatrix = skewmatrix_x
    Elseif (j.Eq.1) Then
       skewmatrix = skewmatrix_y
    Else
       skewmatrix = skewmatrix_z
    Endif


  End Subroutine CondDecision2
  !******************************************************************************
  !> Rotation Matrices for 45° rotation around the main axes
  Subroutine RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)
    Use typedef
    Use ParamData, Only: pi
    Implicit None

    Real(dp) :: phi
    Real(dp), Dimension (3,3), Intent(out) :: rotmatrix_x1,rotmatrix_y1,rotmatrix_z1
    Real(dp), Dimension (3,3), Intent(out) :: rotmatrix_x2,rotmatrix_y2,rotmatrix_z2

    !> \param phi = rotation angle
    !> \param rotmatrix_x1,rotmatrix_y1,rotmatrix_z1, rotmatrix_x2,rotmatrix_y2,rotmatrix_z2 = rotation matrices

    phi = 45._dp/180._dp*pi

    rotmatrix_x1 = (Reshape((/1._dp,0._dp,0._dp,0._dp,Cos(phi),Sin(phi),0._dp,(-Sin(phi)),Cos(phi)/),Shape(rotmatrix_x1)))
    rotmatrix_x2 = (Reshape((/1._dp,0._dp,0._dp,0._dp,Cos(phi),-Sin(phi),0._dp,(Sin(phi)),Cos(phi)/),Shape(rotmatrix_x2)))
    rotmatrix_y1 = (Reshape((/Cos(phi),0._dp,(-Sin(phi)),0._dp,1._dp,0._dp,Sin(phi),0._dp,Cos(phi)/),Shape(rotmatrix_y1)))
    rotmatrix_y2 = (Reshape((/Cos(phi),0._dp,Sin(phi),0._dp,1._dp,0._dp,(-Sin(phi)),0._dp,Cos(phi)/),Shape(rotmatrix_y2)))
    rotmatrix_z1 = (Reshape((/Cos(phi),Sin(phi),0._dp,(-Sin(phi)),Cos(phi),0._dp,0._dp,0._dp,1._dp/),Shape(rotmatrix_z1)))
    rotmatrix_z2 = (Reshape((/Cos(phi),(-Sin(phi)),0._dp,Sin(phi),Cos(phi),0._dp,0._dp,0._dp,1._dp/),Shape(rotmatrix_z2)))


  End Subroutine RotMatrix
  !******************************************************************************
  !> geotropism components of new heading vector
  Subroutine Geocom(tLoc,iord,iax,stndrd,alpha,dxgeo,dygeo,dzgeo,ipl)
    Use typedef
    Implicit None

    Real(dp), Intent(in)::  alpha
    Real(dp) :: tLoc,betapr,geotrp
    Real(dp), Intent(out) :: dxgeo,dygeo,dzgeo
    Real(dp), Intent(in) :: stndrd
    Integer(ap), Intent(in) :: iord,iax,ipl
    !> \param tLoc average value of the current soil temperature values of the 8 corner nodes
    !> \param iord current branch order
    !> \param iax current axis number
    !> \param stndrd length of previous segment
    !> \param alpha azimuth angle
    !> \param dxstr x component of new segment with respect to the soil strength gradient
    !> \param dystr y component of new segment with respect to the soil strength gradient
    !> \param dzstr z component of new segment with respect to the soil strength gradient

    ! calculate geotropism components
    ! (betapr = preferrential heading angle with xy-plane):
    If (iord.Le.4) Then
       Call Prfang(tLoc,iord,iax,geotrp,betapr,ipl)
       dxgeo = stndrd*geotrp*Cos(betapr)*Cos(alpha)
       dygeo = stndrd*geotrp*Cos(betapr)*Sin(alpha)
       dzgeo = stndrd*geotrp*Sin(betapr)
    Else
       dxgeo = 0.0_dp
       dygeo = 0.0_dp
       dzgeo = 0.0_dp
    Endif

    Return
  End Subroutine Geocom
  !************************************************************!*

  !> initial direction of the new axis
  Subroutine Initial(iord,alpha,inaxs,dxini,dyini,dzini)
    Use typedef
    Use typedef
    Use GeoData
    USE RootData, ONLY: rand
    Implicit None

    Real(dp) ::  alpha, theta
    Real(dp), Intent(in)::  inaxs
    Real(dp), Intent(out) :: dxini,dyini,dzini
    Integer(ap), Intent(in) :: iord
    !> \param tLoc average value of the current soil temperature values of the 8 corner nodes
    !> \param iord current branch order
    !> \param iax current axis number
    !> \param stndrd length of previous segment
    !> \param alpha azimuth angle
    !> \param dxini x component of the first segment of the new axis
    !> \param dyini y component of the first segment of the new axis
    !> \param dzini z component of the first segment of the new axis
    !> \param inaxs angle of the first segment of the new axis (from the horizontal plane)

    CALL RANDOM_NUMBER(rand)
    alpha = rand*2._dp*pi
    CALL RANDOM_NUMBER(rand)
    theta = inaxs !rand*(pi-2*inaxs)+inaxs
	 
    ! calculate components of the initial angle of the new axis:  
     IF (iord.LE.4) THEN
       dxini = SIN(alpha)*COS(theta)
       dyini = COS(alpha)*COS(theta)
       dzini = -1*SIN(theta)
    Else
       dxini = 0.0_dp
       dyini = 0.0_dp
       dzini = 0.0_dp
    Endif


    Return
  End Subroutine Initial
  !************************************************************!*
  !> preferential growth angle with horizontal plane
  Subroutine Prfang(tLoc,iord,iax,geotrp,betapr,ipl)
    Use typedef
    Use GeoData
    Use Rootdata, Only: lSomma_growth
    Implicit None

    Real(dp), Intent(in)::  tLoc
    Real(dp), Intent(out) :: geotrp,betapr
    Integer(ap), Intent(in)::  iord,iax,ipl
    Integer(ap)::  igch
    !> \param tLoc average value of the current soil temperature values of the 8 corner nodes
    !> \param iord current branch order
    !> \param iax current axis number
    !> \param geotrp weighting factor for geotropism angle of growth direction vector
    !> \param betapr preferential vertical growth angle

    ! interpolate along piecewise linear function to get
    ! preferred growth angle as a function of axis#, temperature:
    If (iord.Eq.1) Then
       ! use axis data:
       If (tLoc.Ge.tempax(iax,nangax(iax,ipl),ipl)) Then
          ! temperature greater than greatest T for which betapr is specified --
          ! use values for greatest specified T:
          betapr = angaxs(iax,nangax(iax,ipl),ipl)
       Else
          iGch = nangax(iax,ipl)
1         iGch = iGch-1
          If (iGch.Eq.0) Then
             ! temperature smaller than smallest T for which betapr is specified --
             ! use values for smallest specified T:
             betapr = angaxs(iax,1,ipl)
          Else
             If (tLoc.Lt.tempax(iax,iGch,ipl)) Goto 1
             betapr = angaxs(iax,iGch,ipl)+(tLoc-tempax(iax,iGch,ipl))/&
                  (tempax(iax,iGch+1,ipl)-tempax(iax,iGch,ipl))*&
                  (angaxs(iax,iGch+1,ipl)-angaxs(iax,iGch,ipl))
          Endif
       Endif
          IF (lSomma_growth) Then
             geotrp=geoaxs(iord,ipl)
          ELSE 
             geotrp=sg(iord,ipl)
          ENDIF
    Else
       ! use main lateral data:
       If (tLoc.Ge.templt(nanglt(ipl),ipl)) Then
          ! temperature greater than greatest T for which betapr is specified --
          ! use values for greatest specified T:
          ! addition of random component - for iord=1 already done in Input.f90
          betapr = anglat(nanglt(ipl),ipl)
       Else
          iGch = nanglt(ipl)
2         iGch = iGch-1
          If (iGch.Eq.0) Then
             ! temperature smaller than smallest T for which betapr is specified --
             ! use values for smallest specified T:
             betapr = anglat(1,ipl)
          Else
             If (tLoc.Lt.templt(iGch,ipl)) Goto 2
             betapr = anglat(iGch,ipl)+(tLoc-templt(iGch,ipl))/&
                  (templt(iGch+1,ipl)-templt(iGch,ipl))*(anglat(iGch+1,ipl)-anglat(iGch,ipl))
             betapr = anglat(nanglt(ipl),ipl)
          Endif
       Endif
	   
          IF (lSomma_growth) Then
             geotrp = geolat(ipl)
          ELSE
             geotrp = sg(iord,ipl)
          ENDIF
    Endif

    Return
  End Subroutine Prfang
  !******************************************************************************
  !> length and heading vector components of new segment
  Subroutine Nwcomp(igrow,corner_r,newlen,impfac,dx,dy,dz,drivex,drivey,drivez,tempage,ipl)
    Use typedef
    Use RootData
    Use ParamData
    Use tmctrl, Only: dtroot
    Use Environmental, Only: TemLoc
    Use MPIutils, Only: stop_program
    Implicit None

    Real(dp):: newlen,unimplen,stndrd, impfac
    Real(dp), Intent(out) :: dx,dy,dz,drivex,drivey,drivez
    Real(dp)::  dxstr,dystr,dzstr,dxgeo,dygeo,dzgeo
    Real(dp)::  alpha,beta,tLoc,v
    Real(dp), Dimension (3,3) :: condmatrix
    Real(dp)::  tempage,gamma,stdev,u1,u2,r,theta,delta, normfactor
    Integer(ap), Intent(in):: corner_r(8),igrow
    Integer(ap) ::  itip,iprev,ipl

    !> \param igrow current tip or branch number
    !> \param corner 8 surrounding soil nodes of growing tip
    !> \param newlen length of new segment
    !> \param dx x component of new segment
    !> \param dy y component of new segment
    !> \param dz z component of new segment
    !> \param tempage current branch age


    !> use old segment length as normalizing standard for components:
    stndrd = seglen(irecsg(igrow,ipl),ipl)

    !> calculate the old heading components dx,dy,dz:
    itip = irecsg(igrow,ipl)
    iprev = irecpr(itip,ipl)
    dx = xs(itip,ipl)-xs(iprev,ipl)
    dy = ys(itip,ipl)-ys(iprev,ipl)
    dz = zs(itip,ipl)-zs(iprev,ipl)
    ! dx = xg(igrow)-xs(irecsg(igrow))
    ! dy = yg(igrow)-ys(irecsg(igrow))
    ! dz = zg(igrow)-zs(irecsg(igrow))

    !> change old heading angle by some random amount without changing length:
    Call RANDOM_Number(rand)
    gamma = (2*rand-1)*2._dp*pi

    !> add a random angle with the expected value 0 and a calc stdev to the old heading angle
    !> calc random angle from a random sample from a normal(Gaussian) distribution scaled with the segment length and the range of rdmang
    Call Uniera(tempage,ordgrw(igrow,ipl),ibrgrw(igrow,ipl),v,ipl)
    Call Length(v,dtroot,drivex,drivey,drivez,condmatrix,newlen,unimplen,impfac,corner_r,ipl)
    stdev = (unimplen)**0.5*(rdmang(ordgrw(igrow,ipl),ipl))
    If (stdev.Lt.0.0) Then
       Print*, 'standard deviation must be positive'
    End If

    Call RANDOM_Number(rand)
    u1 = rand
    Call RANDOM_Number(rand)
    u2 = rand
    r = Sqrt( -2.0*Log(u1) )
    theta = 2.0*pi*u2
    delta = stdev*r*Sin(theta)
    Call Angchg(dx,dy,dz,alpha,beta,gamma,delta)

    !> get soil strength gradient components:
    If (l_conduc) Then
       Call Ssgcom2(corner_r,condmatrix,ipl)
    Else
       Call Ssgcom1(corner_r,stndrd,strsen(ordgrw(igrow,ipl),ipl),dxstr,dystr,dzstr)
    End If

    !> current local temperature:
    If(ltemp) Call TemLoc(corner_r,tLoc)

    !< get geotropism components:
    Call Geocom(tLoc,ordgrw(igrow,ipl),iaxis(igrow,ipl),stndrd,alpha,dxgeo,dygeo,dzgeo,ipl)

    !> normalize dx, dy dz before adding them up 
    normfactor = 1/SQRT(dx*dx+dy*dy+dz*dz)
    dx     = normfactor*dx
    dy     = normfactor*dy
    dz     = normfactor*dz
       
    !> add up the different direction components according to the used approach (ssgcom1 or ssgcom2)
    If (l_conduc) Then

       drivex = (dx+dxgeo)
       drivey = (dy+dygeo)
       drivez = (dz+dzgeo)


       dx = (condmatrix(1,1)*drivex+condmatrix(2,1)*drivey+condmatrix(3,1)*drivez)
       dy = (condmatrix(1,2)*drivex+condmatrix(2,2)*drivey+condmatrix(3,2)*drivez)
       dz = (condmatrix(1,3)*drivex+condmatrix(2,3)*drivey+condmatrix(3,3)*drivez)


    Else
       dx = (dxstr+dx+dxgeo)
       dy = (dystr+dy+dygeo)
       dz = (dzstr+dz+dzgeo)


    End If

    !> calculate length of new segment, taking location into account:
    Call Uniera(tempage,ordgrw(igrow,ipl),ibrgrw(igrow,ipl),v,ipl)
    Call Length(v,dtroot,drivex,drivey,drivez,condmatrix,newlen,unimplen,impfac,corner_r,ipl)


    Return
  End Subroutine Nwcomp
  !******************************************************************************
  !> changes the previous azimuth and polar angle by a random component each
  Subroutine Angchg(dx,dy,dz,alpha,beta,gamma,delta)
    Use typedef
    Use paramData
    Implicit None

    Real(dp), Intent(inout) :: dx,dy,dz
    Real(dp), Intent(in) :: gamma,delta
    Real(dp), Intent(out) :: alpha,beta
    Real(dp) :: totlen,horlen,a,r,c,d
    !> \param dx x component of new segment
    !> \param dy y component of new segment
    !> \param dz z component of new segment
    !> \param alpha azimuth of the previous segment
    !> \param beta polar angle of the previous segment
    !> \param gamma perturbation of the azimuth
    !> \param delta perturbation of the polar angle

    totlen = Sqrt(dx*dx+dy*dy+dz*dz)
    horlen = Sqrt(dx*dx+dy*dy)
    r      = Sin(delta)*totlen
    a      = Cos(gamma)*r
    beta   = Asin(Sign(Min(1._dp,Abs(dz/totlen)),dz))
    If (horlen.Gt.1.E-20_dp) Then
       alpha = Acos(Sign(Min(1._dp,Abs(dx/horlen)),dx))
       If (dy.Lt.0.0_dp) alpha = 2._dp*pi-alpha
       c  = Sin(beta)*a
       d  = Sin(gamma)*r
       dx = Cos(delta)*dx-c/horlen*dx-Sin(alpha)*d
       dy = Cos(delta)*dy-c/horlen*dy+Cos(alpha)*d
    Else
       alpha = gamma
       dx    = Cos(alpha)*r
       dy    = Sin(alpha)*r
    Endif
    dz     = Cos(delta)*dz+Cos(beta)*a
    horlen = Sqrt(dx*dx+dy*dy)
    beta   = Asin(Sign(Min(1._dp,Abs(dz/totlen)),dz))
    If (horlen.Gt.1.E-20_dp) Then
       alpha = Acos(Sign(Min(1._dp,Abs(dx/horlen)),dx))
       If (dy.Lt.0.0_dp) alpha = 2._dp*pi-alpha
    Else
       alpha = gamma
    Endif

    Return
  End Subroutine Angchg
  !*****************************************************************
  !> mass per length
  Subroutine Maslen(igrow,sLoc,iorder,ipl)
    Use Typedef
    Use RootData, Only: MPL,MPLch,nMPLch,sMPLch
    Implicit None

    Real(dp), Intent(in) :: sLoc
    Integer(ap), Intent(in) :: igrow,iorder,ipl
    Integer(ap) :: iMPL
    !> \param igrow current growing number of tip or branch
    !> \param sLoc average soil strength of the 8 corner nodes around the growing tip
    !> \param iorder current branch order

    ! calculate mass per length as a function of order, soil strength:
    If (sLoc.Ge.sMPLch(iorder,nMPLch(iorder,ipl),ipl)) Then
       MPL(igrow) = MPLch(iorder,nMPLch(iorder,ipl),ipl)
    Else
       iMPL = nMPLch(iorder,ipl)
1      iMPL = iMPL-1
       If ((sLoc.Lt.sMPLch(iorder,iMPL,ipl)).And.(iMPL.Gt.1)) Goto 1
       MPL(igrow) = MPLch(iorder,iMPL,ipl)+(sLoc-sMPLch(iorder,iMPL,ipl))/&
            (sMPLch(iorder,iMPL+1,ipl)-sMPLch(iorder,iMPL,ipl))*&
            (MPLch(iorder,iMPL+1,ipl)-MPLch(iorder,iMPL,ipl))
    Endif

    Return
  End Subroutine Maslen
  !*****************************************************************
 !> grow last segment record; update tip and segment information
  !> merged old subroutines GROW and MKRECD
  Subroutine Grow(igrow,dx,dy,dz,length,mass,rectime,iorder,ipl)
    Use typedef
    Use RootData
    Use PlntData, Only: SpWgt,rootrad
    Use ParamData, Only: pi
    Use GridData, Only: continu
    Use DoussanMat,Only: nsub
    !Use RootGrowthNeighb, Only: Neighb    
    Use MPIutils, Only: stop_program
    Implicit None

    Real(dp), Intent(in) :: mass
    Real(dp), Intent(inout) :: length    
    Integer(ap) :: igrow,iorder,i,ipl,ordgrw_mod,br,countBigLat
    Real(dp) :: dx,dy,dz,factor,rectime, brldummy, naxdummy(diffnum(ipl)+1)
    !> \param igrow current growing number of tip or branch
    !> \param dx x component of new segment
    !> \param dy y component of new segment
    !> \param dz z component of new segment
    !> \param length length of new segment
    !> \param mass (dry) mass of new segment
    !> \param rectime origination time of new segment

    ! increase number of node records

    If (nrec(ipl).Lt.maxrec) Then
       nrec(ipl) = nrec(ipl)+1
    Else
       Call stop_program('Maximum number of root segments -- PROGRAM TERMINATED.')
    Endif

  If (lRootBox_growth) Then
    naxdummy(diffnum(ipl)+1) = 0
    Do i=1,diffnum(ipl)
       naxdummy(i+1) = naxdummy(i) + numlast(i,ipl)
    EndDo
    
    !check, if length + brlgth exceeds brlmax
    If ((iorder.EQ.1).And.(iaxis(igrow,ipl).GT.(naxtot(ipl)-naxdummy(diffnum(ipl)+1)))) Then 
       Do i=1,diffnum(ipl)
          If (iaxis(igrow,ipl).GT.(naxtot(ipl)-naxdummy(i+1))) Then 
             brldummy = maxlast(diffnum(ipl)-i+1,ipl)
             Exit
          Endif
       Enddo
    Else
          If (ordgrw(igrow,ipl).eq.1) Then
             ordgrw_mod=1
          Elseif (ordgrw(igrow,ipl).eq.2) Then !Difference between long and short laterals
              ordgrw_mod=3 !By default small lateral (if big lat it will be changed in dowhile)
              br=1
              countBigLat=0
              Do while ((countBigLat .LT. nBigLat(ipl)) .AND. (br .LE. nbr(ipl)))
                 if (ordgrw(br,ipl).EQ.2) Then
                      if (br .EQ. igrow) ordgrw_mod=2
                 countBigLat=countBigLat+1
                 endif
              br=br+1
              Enddo
          Else
             ordgrw_mod=4
          Endif
          brldummy = brlmax(ordgrw_mod,ipl)
    Endif    
    
    If (brlgth(igrow,ipl).EQ.brldummy) Then 
       length = 0.0_dp
    Elseif (brlgth(igrow,ipl)+length.GE.brldummy) THEN
       length=brldummy-brlgth(igrow,ipl)
    Endif
  Elseif (lSomma_growth) Then
       If (ordgrw(igrow,ipl).eq.1) Then
           ordgrw_mod=1
       Elseif (ordgrw(igrow,ipl).eq.2) Then !Difference between long and short laterals
              ordgrw_mod=3 !By default small lateral (if big lat it will be changed in dowhile)
              br=1
              countBigLat=0
              Do while ((countBigLat .LT. nBigLat(ipl)) .AND. (br .LE. nbr(ipl)))
                 if (ordgrw(br,ipl).EQ.2) Then
                      if (br .EQ. igrow) ordgrw_mod=2
                 countBigLat=countBigLat+1
                 endif
              br=br+1
              Enddo
       Else
           ordgrw_mod=4
       Endif
       brldummy = brlmax(ordgrw_mod,ipl)   
    
    If (brlgth(igrow,ipl).EQ.brldummy) Then 
       length = 0.0_dp
    Elseif (brlgth(igrow,ipl)+length.GE.brldummy) THEN
       length=brldummy-brlgth(igrow,ipl)
    Endif
  Endif      
    
    ! calculate the actual length of components dx,dy,dz:
    factor = Real(length)/Sqrt(dx*dx+dy*dy+dz*dz)
    dx     = factor*dx
    dy     = factor*dy
    dz     = factor*dz

    !If (dx /= dx) dx=0
    !If (dy /= dy) dy=0
    !If (dz /= dz) dz=0

    ! grow at last root node
    xs(nrec(ipl),ipl) = xg(igrow,ipl)+dx
    ys(nrec(ipl),ipl) = yg(igrow,ipl)+dy
    zs(nrec(ipl),ipl) = zg(igrow,ipl)+dz

    ! make record for new root node
    irecpr(nrec(ipl),ipl) = irecsg(igrow,ipl)
    ordseg(nrec(ipl),ipl) = ordgrw(igrow,ipl)
    ibrseg(nrec(ipl),ipl) = ibrgrw(igrow,ipl)
    seglen(nrec(ipl),ipl) = length
    segmas(nrec(ipl),ipl) = mass
    timorg(nrec(ipl),ipl) = rectime
    If (length.Eq.0._dp) seglen(nrec(ipl),ipl) = 1.E-7_dp
        IF (lSomma_growth) THEN
            segrad(nrec(ipl),ipl) = Sqrt(MPL(igrow)/SpWgt(ipl)/pi)
        ELSE
            segrad(nrec(ipl),ipl) = rootrad(iorder,ipl)
        ENDIF
    !print*,'nrec',nrec,'MPL',MPL(igrow),'rad',segrad(nrec)
    If(segrad(nrec(ipl),ipl).Lt.1.E-7_dp) segrad(nrec(ipl),ipl) = 1.E-7_dp
    segsur(nrec(ipl),ipl)= 2._dp*pi*segrad(nrec(ipl),ipl)*seglen(nrec(ipl),ipl)
    crossSectionSeg(nrec(ipl),ipl) = segrad(nrec(ipl),ipl)*segrad(nrec(ipl),ipl)*pi

    If(continu) Then
       Call roottransGrow(xs(nrec(ipl),ipl),ys(nrec(ipl),ipl),nrec(ipl),1,1)
       nsub(nrec(ipl),1) = 1
    End If

    ! update tip information
    xg(igrow,ipl)     = xs(nrec(ipl),ipl)
    yg(igrow,ipl)     = ys(nrec(ipl),ipl)
    zg(igrow,ipl)     = zs(nrec(ipl),ipl)
    irecsg(igrow,ipl) = nrec(ipl)
    
    If(continu) Then
      Call tiptransGrow(xg(igrow,ipl),yg(igrow,ipl),igrow,1,1)
    End If
    
    Return
  End Subroutine Grow
  !********************************************************************************
  !> in case shoot growth is considered, adjust root growth to not exceed available assimilate
  Subroutine Adjust(igrow,newlen,ipl)
    Use typedef
    Use RootData, Only: stopgr,toosml,segmas,ordgrw,irecsg,seglen,xs,ys,zs,xg,yg,zg,brlgth,brlmax,irecpr,ordseg,grwfac,nbr,nBigLat
    Implicit None

    Integer(ap) :: irec,igrow,ipl, ordgrw_mod,br,countBigLat
    Real(dp), Intent(out) :: newlen
    Real(dp) :: factor
    !> \param igrow number of growing tip or branch
    !> \param newlen length of new segment

    ! grwfac > 1.  indicates more assimilate being sent to root than can be used
    ! for growth under current conditions --
    ! assume that extra assimilate is exudated.
    ! find the segment behind the tip 'igrow':
    irec = irecsg(igrow,ipl)

    ! apply growth factor to tentative segment length:
    newlen = seglen(irec,ipl)*Min(grwfac,1._dp)

    ! check if maximum length has been reached:
     If (ordgrw(igrow,ipl).eq.1) Then
         ordgrw_mod=1
     Elseif (ordgrw(igrow,ipl).eq.2) Then !Difference between long and short laterals
        ordgrw_mod=3 !By default small lateral (if big lat it will be changed in dowhile)
        br=1
        countBigLat=0
        Do while ((countBigLat .LT. nBigLat(ipl)) .AND. (br .LE. nbr(ipl)))
          if (ordgrw(br,ipl).EQ.2) Then
             if (br .EQ. igrow) ordgrw_mod=2
             countBigLat=countBigLat+1
          endif
        br=br+1
        Enddo
     Else
        ordgrw_mod=4
     Endif
    IF (brlgth(igrow,ipl)+newlen.GE.brlmax(ordgrw_mod,ipl)) THEN
        newlen = brlmax(ordgrw_mod,ipl)-brlgth(igrow,ipl)
        stopgr(igrow) = .TRUE.
    ELSE
        stopgr(igrow) = .FALSE.
    ENDIF
    
    ! make sure first branch segments are not too small
    ! after being adjusted (prevent div. by zero):
    If (irecpr(irec,ipl).Eq.0) Then
       newlen = Max(newlen,1.E-07_dp)
    Else
       If (ordseg(irecpr(irec,ipl),ipl).Ne.ordseg(irec,ipl)) newlen = Max(newlen,1.E-07_dp)
    Endif
    If (newlen.Gt.1.E-03_dp) Then
       toosml(irec) = .False.
       ! adjust mass of segment:
       segmas(irec,ipl) = segmas(irec,ipl)*Min(grwfac,1._dp)
       ! calculate length correction factor:
       factor = Real(newlen/seglen(irec,ipl))
       ! calculate exact position of tip:
       xg(igrow,ipl)     = xs(irec,ipl)+factor*(xg(igrow,ipl)-xs(irec,ipl))
       yg(igrow,ipl)     = ys(irec,ipl)+factor*(yg(igrow,ipl)-ys(irec,ipl))
       zg(igrow,ipl)     = zs(irec,ipl)+factor*(zg(igrow,ipl)-zs(irec,ipl))
       brlgth(igrow,ipl) = brlgth(igrow,ipl)+newlen
       ! adjust length of segment:
       seglen(irec,ipl) = newlen
    Else
       ! we have a candidate for removal at the end of time step
       toosml(irec) = .True.
    Endif

    Return
  End Subroutine Adjust
  !******************************************************************************
  !> time difference between successive sub-branch origination points
  Subroutine Spacng(igrow,space,impfac,ipl)
    Use typedef
    Use RootData
    Implicit None

    Integer(ap), Intent(in) :: igrow,ipl
    Real(dp), Intent(in) :: impfac
    Real(dp), Intent(out) :: space
    Real(dp) :: v
    !> \param igrow number of growing tip
    !> \param space age distance between two originating branches

    If (lRootBox_growth) Then
       ! convert length/time units 
       Call Uniera(0.0_dp,ordgrw(igrow,ipl),ibrgrw(igrow,ipl),v,ipl)
       space = brspac(ordgrw(igrow,ipl),ipl)/(v*impfac)
    Else
       space = brspac(ordgrw(igrow,ipl),ipl)
    Endif

    Return
  End Subroutine Spacng
  !******************************************************************************
  !> establish new sub-branch origination points along branch 'igrow'
  Subroutine Establ(igrow,newlen,space,t,ipl)
    Use typedef
    Use RootData
    Use tmctrl, Only: dtroot
    Implicit None

    Real(dp), Intent(in):: t,space
    Real(dp) :: first,deltat,over,smin
    Real(dp), Intent(in):: newlen
    Integer(ap), Intent(in):: igrow
    Integer(ap):: ipl
    !> \param igrow number of growing branch or tip
    !> \param newlen length of new segment
    !> \param space age distance between two originating branches
    !> \param t current simulation time

    over = ovrtime(igrow,ipl)
    If (over.Gt.0._dp) Then
       smin = Min(over,space)
       If (smin.Gt.newlen) Then
          first = -1._dp
          over = smin-dtroot
       Else
          first = smin
       Endif
    Else
       If (space.Gt.dtroot) Then
          first = -1._dp
          over = space-dtroot
       Else
          first = space
       Endif
    Endif
    If ((first.Gt.0._dp).And.(nestbl(igrow,ipl).Lt.maxest)) Then
       deltat                      = first
1      nestbl(igrow,ipl)               = nestbl(igrow,ipl)+1
       timest(igrow,nestbl(igrow,ipl)) = t+deltat
       deltat                      = deltat+space
       If ((deltat.Le.dtroot).And.(nestbl(igrow,ipl).Lt.maxest)) Goto 1
       over = deltat-dtroot
    Endif
    ovrtime(igrow,ipl) = over

    Return
  End Subroutine Establ
  !**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!!!!*
  !> root has to stay within domain limits
  Subroutine Boxlim(igrow,ipl)
    Use typedef
    Use RootData
    Use DomData
    Use GridData, Only: continu
    Use tmctrl, Only: dtroot
    Implicit None

    Integer(ap), Intent(in) :: igrow
    Integer(ap) :: irec,iret,ipl
    Real(dp):: dx,dy,dz,newtime
    Real(dp):: factor,fractn,delmin,del,length,newlen,mass,newmas,surf
    Character wall*4
    !> \param igrow number of growing tip

    ! find the segment behind the tip 'igrow':
    ! then check each of the boudaries (x,y,z)
1   iret = irecsg(igrow,ipl) ! tip node
    irec = irecpr(iret,ipl)  ! prev. node
    If (toosml(iret)) Return

    ! check if tip is outside the domain
    ! and if so, which side wall is intersected first:
    delmin = 1.E+37_dp
    wall   = 'none'

    ! domain walls are the growth boundarys
    If(.Not.continu) Then
       If (xs(iret,ipl).Lt.xmin) Then
          del = (xmin-xs(irec,ipl))/(xs(iret,ipl)-xs(irec,ipl))*seglen(iret,ipl)
          If (del.Lt.delmin) Then
             delmin = del
             wall   = 'xmin'
          Endif
       Endif
       If (xs(iret,ipl).Gt.xmax) Then
          del = (xmax-xs(irec,ipl))/(xs(iret,ipl)-xs(irec,ipl))*seglen(iret,ipl)
          If (del.Lt.delmin) Then
             delmin = del
             wall   = 'xmax'
          Endif
       Endif
       If (ys(iret,ipl).Lt.ymin) Then
          del = (ymin-ys(irec,ipl))/(ys(iret,ipl)-ys(irec,ipl))*seglen(iret,ipl)
          If (del.Lt.delmin) Then
             delmin = del
             wall   = 'ymin'
          Endif
       Endif
       If (ys(iret,ipl).Gt.ymax) Then
          del = (ymax-ys(irec,ipl))/(ys(iret,ipl)-ys(irec,ipl))*seglen(iret,ipl)
          If (del.Lt.delmin) Then
             delmin = del
             wall   = 'ymax'
          Endif
       Endif
    End If
    If (zs(iret,ipl).Lt.zmin) Then
       del = (zmin-zs(irec,ipl))/(zs(iret,ipl)-zs(irec,ipl))*seglen(iret,ipl)
       If (del.Lt.delmin) Then
          delmin = del
          wall   = 'zmin'
       Endif
    Endif
    If (zs(iret,ipl).Gt.zmax) Then
       del = (zmax-zs(irec,ipl))/(zs(iret,ipl)-zs(irec,ipl))*seglen(iret,ipl)
       If (del.Lt.delmin) Then
          delmin = del
          wall   = 'zmax'
       Endif
    Endif

    If (wall.Eq.'none') Return
    ! if we get to here (wall.NE. 'none'), we know we have the tip outside the domain and need
    ! new heading components:
    dx = xs(iret,ipl)-xs(irec,ipl)
    dy = ys(iret,ipl)-ys(irec,ipl)
    dz = zs(iret,ipl)-zs(irec,ipl)
    Call Compon(wall(1:1),irec,dx,dy,dz,ipl)
    If (delmin.Gt.1.E-03_dp) Then
       ! break up segment --
       length = seglen(iret,ipl)
       mass   = segmas(iret,ipl)
       surf = segsur(iret,ipl)
       fractn = delmin/length
       ! part with unchanged heading:
       segmas(iret,ipl) = fractn*mass
       seglen(iret,ipl) = delmin
       segsur(iret,ipl)= fractn*surf

       ! pull back tip to start part with new heading:
       xs(iret,ipl) = xs(irec,ipl)+fractn*(xs(iret,ipl)-xs(irec,ipl))
       ys(iret,ipl) = ys(irec,ipl)+fractn*(ys(iret,ipl)-ys(irec,ipl))
       zs(iret,ipl) = zs(irec,ipl)+fractn*(zs(iret,ipl)-zs(irec,ipl))

       ! to prevent numerical problems, pull back exactly to wall
       If (wall.Eq.'xmin') xs(iret,ipl) = xmin
       If (wall.Eq.'xmax') xs(iret,ipl) = xmax
       If (wall.Eq.'ymin') ys(iret,ipl) = ymin
       If (wall.Eq.'ymax') ys(iret,ipl) = ymax
       If (wall.Eq.'zmin') zs(iret,ipl) = zmin
       If (wall.Eq.'zmax') zs(iret,ipl) = zmax

       ! update tip info
       xg(igrow,ipl) = xs(iret,ipl)
       yg(igrow,ipl) = ys(iret,ipl)
       zg(igrow,ipl) = zs(iret,ipl)

       newlen=length-seglen(iret,ipl)
       newmas=mass-segmas(iret,ipl)
       If (newlen.Gt.1.E-03_dp) Then
          newtime = timorg(iret,ipl)+fractn*dtRoot !convert to sp!
          Call Grow(igrow,dx,dy,dz,newlen,newmas,newtime,ordgrw(igrow,ipl),ipl)
       Endif
    Else
       ! don't split, change heading only --
       ! calculate the actual length of components dx,dy,dz:
       factor = seglen(iret,ipl)/Sqrt(dx*dx+dy*dy+dz*dz)
       dx     = factor*dx
       dy     = factor*dy
       dz     = factor*dz
       ! new position of tip:
       xs(iret,ipl) = xs(irec,ipl)+dx
       ys(iret,ipl) = ys(irec,ipl)+dy
       zs(iret,ipl) = zs(irec,ipl)+dz
       ! update tip info
       xg(igrow,ipl) = xs(iret,ipl)
       yg(igrow,ipl) = ys(iret,ipl)
       zg(igrow,ipl) = zs(iret,ipl)
    End If
    Goto 1

  End Subroutine Boxlim
  !******************************************************************************
  !> for cylindrical soil domains --> growth along cylinder walls
  !> if soil geometry is cylindrical, handle like domain edge
  Subroutine Boxlim_Cylinder(igrow,ipl)
    Use typedef
    Use RootData
    Use DomData
    Use tmctrl, Only: dtroot
    Use GridData, Only:rad_cyl,x_cent,y_cent,dxGrid,dyGrid
    Implicit None

    Integer(ap), Intent(in) :: igrow
    Integer(ap) :: irec,iret,ipl
    Real(dp):: irad,newtime
    Real(dp):: factor,fractn,dx,dy,dz
    Real(dp):: azi,xg_cyl,yg_cyl
    Real(dp):: length, mass, surf,delmin,del,newlen,newmas
    Character wall*4
    !> \param igrow number of growing tip

    ! find the segment behind the tip 'igrow':
    ! loop is needed in case of new growth
    ! pull back tip to smallest possible radius
2   iret = irecsg(igrow,ipl) !> tip node
    irec = irecpr(iret,ipl)  !> prev. node
    If (toosml(iret)) Return
    delmin = 1.E+37_dp
    wall   = 'none'


    ! define orientation of segment and tip node with respect to center, azimuth & quadrant
    irad = Sqrt( (xs(iret,ipl)-x_cent)*(xs(iret,ipl)-x_cent) + (ys(iret,ipl)-y_cent)*(ys(iret,ipl)-y_cent))
    azi  = Atan2(ys(iret,ipl)-y_cent,xs(iret,ipl)-x_cent)

    ! old heading components
    dx     = xs(iret,ipl)-xs(irec,ipl)
    dy     = ys(iret,ipl)-ys(irec,ipl)
    dz     = zs(iret,ipl)-zs(irec,ipl)
    xg_cyl = xs(iret,ipl)
    yg_cyl = ys(iret,ipl)

    ! check if tip node is within cylinder (x and y)
    If (irad .Gt. (rad_cyl)/Sqrt(2._dp)) Then
       xg_cyl = Cos(azi)*(rad_cyl-dxGrid)/Sqrt(2._dp)
       yg_cyl = Sin(azi)*(rad_cyl-dyGrid)/Sqrt(2._dp)
       wall   = 'cyli'
    End If

    ! check in z
    If (zs(iret,ipl) .Lt. zmin) Then
       del = (zmin-zs(irec,ipl))/(zs(iret,ipl)-zs(irec,ipl))*seglen(iret,ipl)
       If (Abs(del) .Lt. delmin) Then
          delmin = del
          wall   = 'zmin'
       End If
    Endif
    If (zs(iret,ipl) .Gt. zmax) Then
       del = (zmax-zs(irec,ipl))/(zs(iret,ipl)-zs(irec,ipl))*seglen(iret,ipl)
       If (Abs(del) .Lt. delmin) Then
          delmin = del
          wall   ='zmax'
       End If
    Endif

    If (wall.Eq.'none') Return

    ! length from segment node to wall
    delmin = Sqrt( (xs(iret,ipl)-xs(irec,ipl))**2 + (ys(iret,ipl)-ys(irec,ipl))**2 + (zs(iret,ipl)-zs(irec,ipl))**2 )

    If (delmin.Gt.1.E-03_dp) Then
       ! break up segment --
       length = seglen(iret,ipl)       
       mass = segmas(iret,ipl)
       surf = segsur(iret,ipl)
       fractn = delmin/seglen(iret,ipl) ! length until wall / total length

       ! part with unchanged heading:
       segmas(iret,ipl) = fractn*mass
       seglen(iret,ipl) = delmin
       segsur(iret,ipl) = fractn*surf

    If (wall.Eq.'cyli') Then
       xs(iret,ipl) = xg_cyl
       ys(iret,ipl) = yg_cyl
       dx       = 0._dp
       dy       = 0._dp
       If(Abs(dz).Lt.1.E-10) Then
          If (Abs(zmax-zs(irec,ipl)).Gt.Abs(zmin-zs(irec,ipl))) Then
             dz = zmax-zs(irec,ipl)
          Else
             dz = zmin-zs(irec,ipl)
          Endif
       End If
    Elseif (wall.Eq.'zmin') Then
       zs(iret,ipl) = zmin+0.0001
       dz = 0._dp
       If(Abs(dx).Lt.1.E-10 .And. Abs(dy).Lt.1.E-10) Then
          dx = x_cent-xs(irec,ipl)
          dy = y_cent-ys(irec,ipl)
       End If
    Elseif (wall.Eq.'zmax') Then
       zs(iret,ipl) = zmax-0.0001
       dz = 0._dp
       If(Abs(dx).Lt.1.E-10 .And. Abs(dy).Lt.1.E-10) Then
          dx = x_cent-xs(irec,ipl)
          dy = y_cent-ys(irec,ipl)
       End If
    End If
       
          ! don´t split, change heading only --
          ! calculate the actual length of components dx,dy,dz:
          !factor = seglen(iret,ipl)/Sqrt(dx*dx+dy*dy+dz*dz)
          !dx     = factor*dx
          !dy     = factor*dy
          !dz     = factor*dz

          ! new position of tip:
          !xs(iret,ipl) = xs(irec,ipl)+dx
          !ys(iret,ipl) = ys(irec,ipl)+dy
          !zs(iret,ipl) = zs(irec,ipl)+dz

          ! update tip info
          xg(igrow,ipl) = xs(iret,ipl)
          yg(igrow,ipl) = ys(iret,ipl)
          zg(igrow,ipl) = zs(iret,ipl)


       newlen = seglen(iret,ipl)-delmin
       newmas = segmas(iret,ipl)*(1-fractn)
       
       If (newlen.Gt.1.E-03_dp) Then
          newtime = timorg(irec,ipl)+(1-fractn)*dtRoot ! convert to sp
          Call Grow(igrow,dx,dy,dz,newlen,newmas,newtime,ordgrw(igrow,ipl),ipl)
       End If
    Else
       ! don´t split, change heading only --
       ! calculate the actual length of components dx,dy,dz:
       factor = seglen(iret,ipl)/Sqrt(dx*dx+dy*dy+dz*dz)
       dx     = factor*dx
       dy     = factor*dy
       dz     = factor*dz

       ! new position of tip:
       xs(iret,ipl) = xs(irec,ipl)+dx
       ys(iret,ipl) = ys(irec,ipl)+dy
       zs(iret,ipl) = zs(irec,ipl)+dz

       ! update tip info
       xg(igrow,ipl) = xs(iret,ipl)
       yg(igrow,ipl) = ys(iret,ipl)
       zg(igrow,ipl) = zs(iret,ipl)

    End If
    Goto 2
  End Subroutine Boxlim_Cylinder
  !******************************************************************************
  !> reorients the part of the segment that was outside the domain back into the domain
  Subroutine Compon(wall,irec,dx,dy,dz,ipl)
    Use typedef
    Use RootData
    Use DomData
    Implicit None

    Integer(ap),Intent(in):: irec,ipl
    Real(dp), Intent(inout):: dx,dy,dz
    Character wall*1
    !> \param wall intersection side
    !> \param irec number of segment
    !> \param dx x component of growing tip
    !> \param dy y component of growing tip
    !> \param dz z component of growing tip

    If (wall.Eq.'x') Then
       dx = 0.0_dp
       If ((Abs(dy).Lt.1.E-10_dp).And.(Abs(dz).Lt.1.E-10_dp)) Then
          If (Abs(ymax-ys(irec,ipl)).Gt.Abs(ymin-ys(irec,ipl))) Then
             dy = ymax-ys(irec,ipl)
          Else
             dy = ymin-ys(irec,ipl)
          Endif
          If (Abs(zmax-zs(irec,ipl)).Gt.Abs(zmin-zs(irec,ipl))) Then
             dz = zmax-zs(irec,ipl)
          Else
             dz = zmin-zs(irec,ipl)
          Endif
       Endif
    Elseif(wall.Eq.'y') Then
       dy = 0.0_dp
       If ((Abs(dx).Lt.1.E-10_dp).And.(Abs(dz).Lt.1.E-10_dp)) Then
          If (Abs(xmax-xs(irec,ipl)).Gt.Abs(xmin-xs(irec,ipl))) Then
             dx = xmax-xs(irec,ipl)
          Else
             dx = xmin-xs(irec,ipl)
          Endif
          If (Abs(zmax-zs(irec,ipl)).Gt.Abs(zmin-zs(irec,ipl))) Then
             dz = zmax-zs(irec,ipl)
          Else
             dz = zmin-zs(irec,ipl)
          Endif
       Endif
    Elseif(wall.Eq.'z') Then
       dz = 0.0_dp
       If ((Abs(dx).Lt.1.E-10_dp).And.(Abs(dy).Lt.1.E-10_dp)) Then
          If (Abs(xmax-xs(irec,ipl)).Gt.Abs(xmin-xs(irec,ipl))) Then
             dx = xmax-xs(irec,ipl)
          Else
             dx = xmin-xs(irec,ipl)
          Endif
          If (Abs(ymax-ys(irec,ipl)).Gt.Abs(ymin-ys(irec,ipl))) Then
             dy = ymax-ys(irec,ipl)
          Else
             dy = ymin-ys(irec,ipl)
          Endif
       Endif
    Endif

    Return
  End Subroutine Compon
  !******************************************************************************
  !> remove all segments from list that are too small
  Subroutine Remove(nrecol,ipl)
    Use typedef
    Use RootData
    Implicit None

    Integer(ap) :: nrecnw,nrecol,ifrom,ito,igrow,ipl
    !> \param nrecol total number of segment at the beginning of the call to ROOT

    ! set nrecnw equal to the current number of segment records
    nrecnw = nrec(ipl)
    ito    = 0
    ifrom  = nrecol
10  ifrom = ifrom+1
    If (toosml(ifrom)) Then
       ! identify and pull back tip that belongs to the record being removed:
       igrow = 0
11     igrow = igrow+1
       If (ibrgrw(igrow,ipl).Ne.ibrseg(ifrom,ipl)) Goto 11
       xg(igrow,ipl)     = xs(ifrom,ipl)
       yg(igrow,ipl)     = ys(ifrom,ipl)
       zg(igrow,ipl)     = zs(ifrom,ipl)
       irecsg(igrow,ipl) = irecpr(ifrom,ipl)
       ! decrease total number of records by '1':
       nrec(ipl) = nrec(ipl)-1
       ! if this the first candidate for removal (to be overwritten),
       ! mark position:
       If (ito.Eq.0) ito = ifrom
    Else If (ito.Gt.0) Then
       ! can move from 'ifrom' because at least one removal candidate
       ! has been previously identified and marked as 'ito':
       toosml(ito) = toosml(ifrom)
       xs(ito,ipl)     = xs(ifrom,ipl)
       ys(ito,ipl)     = ys(ifrom,ipl)
       zs(ito,ipl)     = zs(ifrom,ipl)
       irecpr(ito,ipl) = irecpr(ifrom,ipl)
       ordseg(ito,ipl) = ordseg(ifrom,ipl)
       ibrseg(ito,ipl) = ibrseg(ifrom,ipl)
       seglen(ito,ipl) = seglen(ifrom,ipl)
       segsur(ito,ipl) = segsur(ifrom,ipl)
       segmas(ito,ipl) = segmas(ifrom,ipl)
       timorg(ito,ipl) = timorg(ifrom,ipl)
       ! tip that belongs to moved segment needs the new segment reference, ito:
       igrow = 0
12     igrow = igrow+1
       If (ibrgrw(igrow,ipl).Ne.ibrseg(ito,ipl)) Goto 12
       irecsg(igrow,ipl) = ito
       ito = ito+1
    Endif
    If (ifrom.Lt.nrecnw) Goto 10
    ! now nrec is the actual number of segment records
    Return
  End Subroutine Remove
  !*********************************************************************************
  !> Update root records after removal of segments
  Subroutine Update(ngrwnw,ipl)
    Use typedef
    Use RootData
    Use GridData, Only: continu
    Use DoussanMat, Only: transroot
    Implicit None

    Integer(ap):: ito,ifrom,iest,ngrwnw,ipl
    !> \param ngrwnw number of growing tips

    ! may have more growing branch tips at the end of time step:
    ngrow(ipl)  = ngrwnw
    ito    = 0
    ifrom  = 0
10  ifrom = ifrom+1
  IF (stopgr(ifrom)) THEN
     ngrow = ngrow-1
     ! if this the first candidate for removal (to be overwritten),
     ! mark position:
     IF (ito.EQ.0) ito = ifrom
    ELSEIf (ito.Gt.0) Then
       ! can move from 'ifrom' because at least one removal candidate
       ! has been previously identified and marked as 'ito':
       stopgr(ito)  = stopgr(ifrom)
       xg(ito,ipl)      = xg(ifrom,ipl)
       yg(ito,ipl)      = yg(ifrom,ipl)
       zg(ito,ipl)      = zg(ifrom,ipl)
       iaxis(ito,ipl)   = iaxis(ifrom,ipl)
       irecsg(ito,ipl)  = irecsg(ifrom,ipl)
       ordgrw(ito,ipl)  = ordgrw(ifrom,ipl)
       ibrgrw(ito,ipl)  = ibrgrw(ifrom,ipl)
       brlgth(ito,ipl)  = brlgth(ifrom,ipl)
       ovrtime(ito,ipl) = ovrtime(ifrom,ipl)
       nestbl(ito,ipl)  = nestbl(ifrom,ipl)
       Do iest=1,nestbl(ito,ipl)
          timest(ito,iest) = timest(ifrom,iest)
       End Do
       If(continu) transroot(ito,1:2,1,1) = transroot(ifrom,1:2,1,1)
       ito = ito+1
    Endif
    If (ifrom.Lt.ngrwnw) Goto 10
    ! now ngrow is the actual number of growing tips at beginning of next step
    Return
  End Subroutine Update
  !*********************************************************************************
  !> additional secondary radial growth
  Subroutine Secondary_Radial_Growth(dt,ipl)
    Use typedef
    Use RootData, Only: nrec,segsur,f_rad,ordseg
    Implicit None
    Real(dp),Intent(in) :: dt
    Integer(ap):: irec
    Integer(ap),Intent(in)::ipl

    Do irec = 1,nrec(ipl)
       segsur(irec,ipl) = segsur(irec,ipl) + segsur(irec,ipl)*(f_rad(ordseg(irec,ipl),ipl)-1)*dt
    End Do

  End Subroutine Secondary_Radial_Growth
  !********************************************************************************
  !> ### for continuous domains determine length of the root that is growing back into the
  !> other side of the soil domain ###
  Subroutine roottransGrow(xt,yt,inode,isub,ipl)
    Use typedef
    Use GridData, Only: nex,ney,dxgrid,dygrid
    Use GridData2, Only: nex2,ney2,dxgrid2,dygrid2
    Use DoussanMat, Only: transroot
    USE RootData, ONLY : ltwo_grids
    Use DomData
    Implicit None

    Integer(ap), Intent(in) :: inode,ipl,isub
    Real(dp), Intent(in) :: xt,yt
    Integer(ap) :: nex_,ney_
    Real(dp) :: xd,yd,dxgrid_,dygrid_,xmin_,xmax_,ymin_,ymax_

    If (ltwo_grids) Then 
       nex_ = nex2
       ney_ = ney2
       dxgrid_ = dxgrid2
       dygrid_ = dygrid2
       xmin_ = xmin2
       xmax_ = xmax2
       ymin_ = ymin2
       ymax_ = ymax2
    Else 
       nex_ = nex
       ney_ = ney
       dxgrid_ = dxgrid
       dygrid_ = dygrid
       xmin_ = xmin
       xmax_ = xmax
       ymin_ = ymin
       ymax_ = ymax
    Endif

    xd = xt
    yd = yt
    transroot(inode,1:2,isub,ipl)=0
    Do While (xd.Ge.xmax_) 
       xd=xd-nex_*dxgrid_
       transroot(inode,1,1,ipl)=transroot(inode,1,isub,ipl)-1
    End Do
    Do While (xd.Lt.xmin_)
       xd=xd+nex_*dxgrid_
       transroot(inode,1,1,ipl)=transroot(inode,1,isub,ipl)+1
    End Do
    Do While (yd.Ge.ymax_)
       yd=yd-ney_*dygrid_
       transroot(inode,2,1,ipl)=transroot(inode,2,isub,ipl)-1
    End Do
    Do While (yd.Lt.ymin_)
       yd=yd+ney_*dygrid_
       transroot(inode,2,1,ipl)=transroot(inode,2,isub,ipl)+1
    End Do

  End Subroutine roottransGrow
  !******************************************************************************** 
  !> ### for continuous domains determine length of the root that is growing back into the
  !> other side of the soil domain ###
  Subroutine tiptransGrow(xt,yt,inode,isub,ipl)
    Use typedef
    Use GridData, Only: nex,ney,dxgrid,dygrid
    Use GridData2, Only: nex2,ney2,dxgrid2,dygrid2
    Use DoussanMat, Only: transtip
    USE RootData, ONLY : ltwo_grids
    Use DomData
    Implicit None

    Integer(ap), Intent(in) :: inode,ipl,isub
    Real(dp), Intent(in) :: xt,yt
    Integer(ap) :: nex_,ney_
    Real(dp) :: xd,yd, dxgrid_,dygrid_,xmin_,xmax_,ymin_,ymax_
		
    If (ltwo_grids) Then 
       nex_ = nex2
       ney_ = ney2
       dxgrid_ = dxgrid2
       dygrid_ = dygrid2
       xmin_ = xmin2
       xmax_ = xmax2
       ymin_ = ymin2
       ymax_ = ymax2
    Else 
       nex_ = nex
       ney_ = ney
       dxgrid_ = dxgrid
       dygrid_ = dygrid
       xmin_ = xmin
       xmax_ = xmax
       ymin_ = ymin
       ymax_ = ymax
    Endif

    xd = xt
    yd = yt
    transtip(inode,1:2,isub,ipl)=0
    Do While (xd.Ge.xmax_)
       xd=xd-nex_*dxgrid_
       transtip(inode,1,1,ipl)=transtip(inode,1,isub,ipl)-1
    End Do
    Do While (xd.Lt.xmin_)
       xd=xd+nex_*dxgrid_
       transtip(inode,1,1,ipl)=transtip(inode,1,isub,ipl)+1
    End Do
    Do While (yd.Ge.ymax_)
       yd=yd-ney_*dygrid_
       transtip(inode,2,1,ipl)=transtip(inode,2,isub,ipl)-1
    End Do
    Do While (yd.Lt.ymin_)
       yd=yd+ney_*dygrid_
       transtip(inode,2,1,ipl)=transtip(inode,2,isub,ipl)+1
    End Do

  End Subroutine tiptransGrow
  !********************************************************************************

End Module RootGrowth
