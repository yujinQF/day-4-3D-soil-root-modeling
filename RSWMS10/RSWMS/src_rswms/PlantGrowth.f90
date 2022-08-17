!> \file PlantGrowth.f90
!! \brief current relative stress due to soil strength and solute concentration:

!> Module PlantGrowth
Module PlantGrowth

Contains

  Subroutine Stress(rs,concrs)
    Use ParamData, Only: lChem
    Use PlntData
    Use RootData, Only: sAvg,cAvg
    Use StrData
    Implicit None

    Integer(ap) ::  k
    Real(dp), Intent(out) :: rs,concrs
    ! interpolate along piecewise linear function to get relative stress, rs,
    ! corresponding to actual average soil strength sAvg;
    ! assume rs = 0.0 at zero soil strength;
    ! assume rs = 1.0 at sImp (soil strength at which growth ceases completely):
    If (sAvg.Ge.sImp) Then
       rs=1.0_dp
    Else If (sAvg.Ge.sc(ns)) Then
       rs=rsc(ns)+(sAvg-sc(ns))/(sImp-sc(ns))*(1.-rsc(ns))
    Else
       k=ns
1      k=k-1
       If (k.Eq.0) Then
          ! soil strength smaller than smallest ss for which rs is specified:
          rs=sAvg/sc(1)*rsc(1)
       Else
          If(sc(k).Gt.sAvg) Goto 1
          rs=rsc(k)+(sAvg-sc(k))/(sc(k+1)-sc(k))*(rsc(k+1)-rsc(k))
       Endif
    Endif
    ! now interpolate along piecewise linear function to get relative stress, concrs,
    ! corresponding to actual average solute concentration cAvg;
    ! assume concrs = 0.0 for average solute concentration within the optimal range;
    ! assume concrs = 1.0 for average solute concentration outside the range of
    ! minimum and maximum concentration allowed (at which growth ceases completely):
    If (lChem) Then
       k=0
11     k=k+1
       If (k.Gt.ncnc) Then
          concrs=rscnc(ncnc)
       Else
          If (cAvg.Gt.cncp(k)) Goto 11
          concrs=rscnc(k-1)+(cAvg-cncp(k-1))/(cncp(k)-cncp(k-1))*(rscnc(k)-rscnc(k-1))
       Endif
    Else
       concrs=0.0_dp
    End If
    Return
  End Subroutine Stress
  !***********************************************************************************************
  ! current potential transpiration:
  Subroutine SetTp(t,rs,concrs,ipl)
    Use ParamData, Only: lChem
    Use PlntData
    Use DoussanMat, Only: BCr_usr
    Use TempData
    Use Environmental, Only: Atmosphere
    Implicit None

    Integer(ap), Intent(in) :: ipl
    Real(dp), Intent(in) :: t,rs,concrs
    Integer(ap) :: ifc,k
    Real(dp) :: fcTpLA,fTpLA

    !Call Atmosphere(t)

    ! calculate current TpLA-value:
    ifc=1
301 ifc=ifc+1
    If (ifc.Gt.ntTpLA) Then
       TpLA=TpLAc(ntTpLA)
    Else
       If (t.Gt.tTpLA(ifc)) Goto 301
       TpLA=TpLAc(ifc-1)+(TpLAc(ifc)-TpLAc(ifc-1))*(t-tTpLA(ifc-1))/(tTpLA(ifc)-tTpLA(ifc-1))
    Endif
    ! interpolate along piecewise linear function to get transp. reduction factor,
    ! fTpLA, corresponding to current relative stress due to soil strength, rs;
    ! by definition, fTpLA = 1.0 at zero stress (rs = 0.0)
    If (rs.Ge.sfTpLA(nfTpLA)) Then
       !relative stress greater than greatest rs for which fTpLA is specified:
       fTpLA=fTpLAc(nfTpLA)
    Else
       k=nfTpLA
302    k=k-1
       If (k.Eq.0) Then
          ! relative stress smaller than smallest rs for which fTpLA is specified:
          fTpLA=1.0_dp+rs/sfTpLA(1)*(fTpLAc(1)-1._dp)
       Else
          If (sfTpLA(k).Gt.rs) Goto 302
          fTpLA=fTpLAc(k)+(rs-sfTpLA(k))/(sfTpLA(k+1)-sfTpLA(k))*(fTpLAc(k+1)-fTpLAc(k))
       Endif
    Endif
    ! interpolate along piecewise linear function to get transp. reduction factor,
    ! fcTpLA, corresponding to current relative stress due to solute conc., concrs;
    ! by definition, fcTpLA = 1.0 at zero stress (concrs = 0.0)
    If (lChem) Then
       If (concrs.Ge.scTpLA(ncTpLA)) Then
          ! relative stress greater than greatest concrs for which fTpLA is specified:
          ! fcTpLA=cTpLAc(ncTpLA)
       Else
          k=ncTpLA
303       k=k-1
          If (k.Eq.0) Then
             ! relative stress smaller than smallest concrs for which fTpLA is specified:
             fcTpLA=1.0_dp+concrs/scTpLA(1)*(cTpLAc(1)-1._dp)
          Else
             If (scTpLA(k).Gt.concrs) Goto 303
             fcTpLA=cTpLAc(k)+(concrs-scTpLA(k))/(scTpLA(k+1)-scTpLA(k))*(cTpLAc(k+1)-cTpLAc(k))
          Endif
       Endif
    Else
       fcTpLA=1.0_dp
    Endif
    Tpot(ipl)=TpLA*LA(ipl)
    BCr_usr(ipl)=-Abs(Tpot(ipl)*Min(fTpLA,fcTpLA)) !root BC is negative, while transpiration is positive... 
    Return
  End Subroutine SetTp
  !*********************************************************************************************
  ! integrate uptake over spatial domain:
  Subroutine ActTrs
    Use GridData
    Use PlntData, Only : Tact
    Implicit None

    Integer(ap) :: corner(8),ie,ic
    Real(dp) :: sum,sume

    ! calculate actual overall transpiration rate from soil domain:
    sum=0.0_dp
    Do  iE=1,nElm!,2
       ! assign cuboid corner nodes:
       corner(1)=elmnod(1,iE)
       corner(2)=elmnod(2,iE)
       corner(3)=elmnod(3,iE)
       corner(4)=elmnod(4,iE)
       corner(5)=elmnod(5,iE)
       corner(6)=elmnod(6,iE)
       corner(7)=elmnod(7,iE)
       corner(8)=elmnod(8,iE)
       ! add up average cuboid sink terms, integrate over volume:
       sumE=0.0_dp
       Do ic=1,8
          sumE=sumE+sink(corner(ic))
       End Do
       sum=sum+sumE
    End Do
    Tact=sum*dxGrid*dyGrid*dzGrid/8._dp
    Return
  End Subroutine ActTrs
  !********************************************************************************
  ! current water use efficiency:
  Subroutine Effncy(t,rs,concrs,W)
    Use ParamData, Only: lChem
    Use PlntData
    Implicit None

    Integer(ap) ::  ifc,k
    Real(dp) ::  t,rs,concrs,W,fcw,fw

    ! calculate current W:
    ifc=0
1   ifc=ifc+1
    If (ifc.Gt.ntW) Then
       W=Wc(ntW)
    Else
       If (t.Gt.tW(ifc)) Goto 1
       W=Wc(ifc-1)+(Wc(ifc)-Wc(ifc-1))*(t-tW(ifc-1))/(tW(ifc)-tW(ifc-1))
    Endif
    ! interpolate along piecewise linear function to get water use efficiency reduction
    ! factor, fW, corresponding to current relative stress due to soil strength, rs;
    ! by definition, fW = 1.0 at zero stress (rs = 0.0)
    If (rs.Ge.sfW(nsfW)) Then
       ! relative stress greater than greatest rs for which fW is specified:
       fW=fwc(nsfW)
    Else
       k=nsfW
2      k=k-1
       If (k.Eq.0) Then
          ! relative stress smaller than smallest rs for which fW is specified:
          fW=1.0_dp+rs/sfW(1)*(fWc(1)-1._dp)
       Else
          If (sfW(k).Gt.rs) Goto 2
          fW=fWc(k)+(rs-sfW(k))/(sfW(k+1)-sfW(k)) *(fWc(k+1)-fWc(k))
       Endif
    Endif
    ! interpolate along piecewise linear function to get water use efficiency
    ! reduction factor, fcW, corredponding to current relative stress, concrs;
    ! by definition, fcW = 1.0 at zero stress (concrs = 0.0)
    If (lChem) Then
       If (concrs.Ge.scW(nscW)) Then
          ! relative stress greater than greatest concrs for which fcW is dpecified:
          fcW=cWc(nscW)
       Else
          k=nscW
3         k=k-1
          If (k.Eq.0) Then
             ! relative stress smaller than smallest rs for which fW is dpecified:
             fcW=1.0_dp+concrs/scW(1)*(cWc(1)-1._dp)
          Else
             If (scW(k).Gt.concrs) Goto 3
             fcW=cWc(k)+(concrs-scW(k))/(scW(k+1)-scW(k))*(cWc(k+1)-cWc(k))
          Endif
       Endif
    Else
       fcW=1.0_dp
    Endif
    W=W*Min(fW,fcW)
    Return
  End Subroutine Effncy
!********************************************************************************
! current root/shoot partitioning ratio:
  Subroutine Ratio(t,rs,concrs,RSR)
    Use ParamData, Only: lChem
    Use PlntData
    Implicit None

    Integer(ap) ::  ifc,k
    Real(dp) :: t,rs,concrs,RSR,fcrsr,frsr

    ! calculate current RSR-values for each soil strength:
    ifc=0
1   ifc=ifc+1
    If (ifc.Gt.ntRSR) Then
       RSR=RSRc(ntRSR)
    Else
       If (t.Gt.tRSR(ifc)) Goto 1
       RSR=RSRc(ifc-1)+(RSRc(ifc)-RSRc(ifc-1))*(t-tRSR(ifc-1))/(tRSR(ifc)-tRSR(ifc-1))
    Endif
    ! interpolate along piecewise linear function to get root/shoot ratio reduction
    ! factor, fRSR, corresponding to current relative stress due to soil strength, rs;
    !by definition, fRSR = 1.0 at zero stress (rs = 0.0)
    If (rs.Ge.sfRSR(nsfRSR)) Then
       ! relative stress greater than greatest rs for which fRSR is specified:
       fRSR=fRSRc(nsfRSR)
    Else
       k=nsfRSR
2      k=k-1
       If (k.Eq.0) Then
          ! relative stress smaller than smallest rs for which fRSR is specified:
          fRSR=1._dp+rs/sfRSR(1)*(fRSRc(1)-1._dp)
       Else
          If (sfRSR(k).Gt.rs) Goto 2
          fRSR=fRSRc(k)+(rs-sfRSR(k))/(sfRSR(k+1)-sfRSR(k))*(fRSRc(k+1)-fRSRc(k))
       Endif
    Endif
    !interpolate along piecewise linear function to get root/shoot ratio reduction factor,
    ! fcRSR, corresponding to current relative stress due to solute concentration, concrs;
    ! by definition, fcRSR = 1.0 at zero stress (concrs = 0.0)
    If (lChem) Then
       If (concrs.Ge.scRSR(nscRSR)) Then
          ! relative stress greater than greatest concrs for which fcRSR is specified:
          fcRSR=cRSRc(nscRSR)
       Else
          k=nscRSR
3         k=k-1
          If (k.Eq.0) Then
             ! relative stress smaller than smallest concrs for which fcRSR is specified:
             fcRSR=1.0_dp+concrs/scRSR(1)*(cRSRc(1)-1._dp)
          Else
             If (scRSR(k).Gt.concrs) Goto 3
             fcRSR=cRSRc(k)+(concrs-scRSR(k))/(scRSR(k+1)-scRSR(k))*(cRSRc(k+1)-cRSRc(k))
          Endif
       Endif
    Else
       fcRSR=1.0_dp
    Endif
    RSR=RSR*Max(fRSR,fcRSR)
    Return
  End Subroutine Ratio
!*****************************************************************
! current leaf area increase per new dry shoot mass:
  Subroutine Leaves(t,LAmshv)
    Use PlntData
    Implicit None

    Integer(ap) ::  ifc
    Real(dp) ::  LAmshv,t 

    !calculate current LA/msh:
    ifc=0
1   ifc=ifc+1
    If (ifc.Gt.ntLA) Then
       LAmshv=LAc(ntLA)
    Else
       If (t.Gt.tLA(ifc)) Goto 1
       LAmshv=LAc(ifc-1)+(LAc(ifc)-LAc(ifc-1))*(t-tLA(ifc-1))/(tLA(ifc)-tLA(ifc-1))
    Endif
    Return
  End Subroutine Leaves
!*****************************************************************

End Module PlantGrowth
