!> \file Solute.f90
!> \brief

!> Module Solute
Module SoluteMod
  Implicit None

Contains

  Subroutine SOLUTE(dt,t,tPulse,dtMaxC,KodCB,IAD,IADN,IADD,iCount)
    Use typedef
    Use BoundData
    !USE ParamData, ONLY: maxbnd
    Use GridData, Only: nband,nPt,nElm,elmnod,iL,subN,sink,deter,csink,ci,di,bi,ax,ay,az
    Use SolData, Only: Dispxx,Dispxy,Dispyy,Dispyz,Dispzz,Dispxz,Fc,Gc,conc,Nlevel,matnum,chpar,vx,vy,vz,epsi,theta_old,theta,cono,con,hnew,hold,kode
    Use WatFun
    Use CumData
    Use MatData, Only: As,Bs,A1
    Use Orthofem, Only: Find, OrthoMin, ILU
    Implicit None

    Integer(ap),Intent(inout) :: KodCB(maxbdr),iCount
    Integer(ap):: i,j,k,l,m,iE,iSE,ic,j1,i1,j2,i2,ib,lev,List(4),kk
    Integer(ap):: IAD(nband,nPt),IADN(nPt),IADD(nPt),north
    Real(dp),Intent(in)::dt,t,tpulse
    Real(dp),Intent(inout) :: dtmaxc
    Real(dp) :: VxE(4),VyE(4),VzE(4),S(4,4)
    Real(dp) :: Vzz,Vyy,Vxx
    Real(dp) :: ec1,gce,fmul,rootch,vzee,vyee,vxee,cone,VEl,cayz,s_m
    Real(dp) :: smul2,smul1,AcE,fce,ec6,ec5,ec4,ec3,ec2,CAzz,alf,CAyy
    Real(dp) :: dpom,CAxz,CAxy,CAxx
    Real(dp), Allocatable,Dimension (:) :: Ac
    Real(dp), Allocatable,Dimension (:) :: B1, F ,DS

    !REAL(dp)::Wx(4),Wy(4),Wz(4),W12,W13,W14,W23,W24,W34
    !REAL(dp)::A11,A12,A13,A14,A21,A22,A23,A24,A31,A32,A33,A34,A41,A42,A43,A44

    Allocate(As(nband,nPt))
    Allocate(A1(nband,nPt))
    Allocate(B1(nPt),F(nPt),DS(nPt))
    Allocate(Ac(nPt),Bs(nPt))
    If (.Not.Allocated(Dispxx)) Then
       Allocate(Dispxx(nPt),Dispyy(nPt),Dispzz(nPt))
       Allocate(Dispxy(nPt),Dispxz(nPt),Dispyz(nPt))
    Endif
    If(.Not.Allocated(Fc)) Allocate(Fc(nPt),Gc(nPt))
    If(.Not.Allocated(Qc)) Allocate(Qc(nPt))
    Gc = 0._dp
    Fc = 0._dp
    Ac = 0._dp
    Qc = 0._dp
    F =  0._dp
    DS = 0._dp
    Bs = 0._dp
    AcE=0.0
    !**********************************************************************************
    If (iCount.Eq.0) Call ChInit(dtMaxC,dt)
    ! Initialuization
    alf=1._dp-epsi

    If (t.Gt.tPulse) Then
       cBound=0.0_dp
    Endif
    Bs =0.0_dp
    B1 = Conc
    Qc=0.0_dp
    Do i=1,Npt
       If (epsi.Lt.0.001_dp) Then
          As(iadd(i),i)=0.0_dp
       Else
          Do j=1,nband!2*Nband-1
             As(j,i)=0.0_dp
          End Do
       Endif
    End Do

    Do Lev=1,NLevel
       If (Lev.Eq.NLevel) Then
          Call Veloc
          Call Disper
          Call PeCour(dtMaxC,dt)
          !IF (lUpW) CALL WeFact
       Else
          Call Disper
       Endif

       Do i=1,Npt
          M=MatNum(i)
          If (Lev.Ne.NLevel) Then
             !if(.not.lUpW) then
             DPom=dt/6._dp/(theta_old(i)+ChPar(1,M)*ChPar(5,M))
             Dispxx(i)=Dispxx(i)+Vx(i)*Vx(i)*DPom
             Dispyy(i)=Dispyy(i)+Vy(i)*Vy(i)*DPom
             Dispzz(i)=Dispzz(i)+Vz(i)*Vz(i)*DPom
             Dispxy(i)=Dispxy(i)+Vx(i)*Vy(i)*DPom
             Dispxz(i)=Dispxz(i)+Vx(i)*Vz(i)*DPom
             Dispyz(i)=Dispyz(i)+Vy(i)*Vz(i)*DPom
             !end if
          Else
             Ac(i)=-(theta_old(i)*alf+theta(i)*epsi)-ChPar(1,M)*ChPar(5,M)
             !if(.not.lUpW) then
             DPom=dt/6._dp/(theta(i)+ChPar(1,M)*ChPar(5,M))
             Dispxx(i)=Dispxx(i)-Vx(i)*Vx(i)*DPom
             Dispyy(i)=Dispyy(i)-Vy(i)*Vy(i)*DPom
             Dispzz(i)=Dispzz(i)-Vz(i)*Vz(i)*DPom
             Dispxy(i)=Dispxy(i)-Vx(i)*Vy(i)*DPom
             Dispxz(i)=Dispxz(i)-Vx(i)*Vz(i)*DPom
             Dispyz(i)=Dispyz(i)-Vy(i)*Vz(i)*DPom
             !end if
             Gc(i)=ChPar(8,M)*theta(i)+ChPar(1,M)*ChPar(9,M)
             Fc(i)=ChPar(6,M)*theta(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)-csink(i)+sink(i)
          Endif
       End Do

       Do i=1,Npt
          If (Lev.Eq.NLevel) DS(i)=0.0_dp
       End Do
       ! Loop on elements

       Do iE=1,Nelm
          CAxx=ConAxx(iE)
          CAyy=ConAyy(iE)
          CAzz=ConAzz(iE)
          CAxy=ConAxy(iE)
          CAxz=ConAxz(iE)
          CAyz=ConAyz(iE)
          ! Loop on subelements
          Do iSE=1,5
             i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
             j=elmnod(iL(2,iSE,subN(iE)),iE)
             k=elmnod(iL(3,iSE,subN(iE)),iE)
             l=elmnod(iL(4,iSE,subN(iE)),iE)
             List(1)=i
             List(2)=j
             List(3)=k
             List(4)=l
             VEl=Abs(Deter(iSE,iE))/6._dp!Deter, Ax, Ay, Az from calcgeom (Couvreur oct 2010)
             !            Calculate Velocities
             If (Lev.Eq.NLevel) Then
                Vxx=(Ax(1,iSE,iE)*hNew(i)+Ax(2,iSE,iE)*hNew(j)+Ax(3,iSE,iE)*hNew(k)+Ax(4,iSE,iE)*hNew(l))/Deter(iSE,iE)+CAxz
                Vyy=(Ay(1,iSE,iE)*hNew(i)+Ay(2,iSE,iE)*hNew(j)+Ay(3,iSE,iE)*hNew(k)+Ay(4,iSE,iE)*hNew(l))/Deter(iSE,iE)+CAyz
                Vzz=(Az(1,iSE,iE)*hNew(i)+Az(2,iSE,iE)*hNew(j)+Az(3,iSE,iE)*hNew(k)+Az(4,iSE,iE)*hNew(l))/Deter(iSE,iE)+CAzz
             Else
                Vxx=(Ax(1,iSE,iE)*hOld(i)+Ax(2,iSE,iE)*hOld(j)+Ax(3,iSE,iE)*hOld(k)+Ax(4,iSE,iE)*hOld(l))/Deter(iSE,iE)+CAxz
                Vyy=(Ay(1,iSE,iE)*hOld(i)+Ay(2,iSE,iE)*hOld(j)+Ay(3,iSE,iE)*hOld(k)+Ay(4,iSE,iE)*hOld(l))/Deter(iSE,iE)+CAyz
                Vzz=(Az(1,iSE,iE)*hOld(i)+Az(2,iSE,iE)*hOld(j)+Az(3,iSE,iE)*hOld(k)+Az(4,iSE,iE)*hOld(l))/Deter(iSE,iE)+CAzz
             End If
             If (Lev.Ne.NLevel) Then
                ConE=(ConO(i)+ConO(j)+ConO(k)+ConO(l))/4._dp
                VxE(1)=-ConO(i)*Vxx
                VxE(2)=-ConO(j)*Vxx
                VxE(3)=-ConO(k)*Vxx
                VxE(4)=-ConO(l)*Vxx
                VyE(1)=-ConO(i)*Vyy
                VyE(2)=-ConO(j)*Vyy
                VyE(3)=-ConO(k)*Vyy
                VyE(4)=-ConO(l)*Vyy
                VzE(1)=-ConO(i)*Vzz
                VzE(2)=-ConO(j)*Vzz
                VzE(3)=-ConO(k)*Vzz
                VzE(4)=-ConO(l)*Vzz
             Else
                ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4._dp
                VxE(1)=-Con(i)*Vxx
                VxE(2)=-Con(j)*Vxx
                VxE(3)=-Con(k)*Vxx
                VxE(4)=-Con(l)*Vxx
                VyE(1)=-Con(i)*Vyy
                VyE(2)=-Con(j)*Vyy
                VyE(3)=-Con(k)*Vyy
                VyE(4)=-Con(l)*Vyy
                VzE(1)=-Con(i)*Vzz
                VzE(2)=-Con(j)*Vzz
                VzE(3)=-Con(k)*Vzz
                VzE(4)=-Con(l)*Vzz
             End If
             VxEE=-ConE*Vxx
             VyEE=-ConE*Vyy
             VzEE=-ConE*Vzz

             If (Lev.Eq.1) Then
                RootCh=VEl*dt*(Conc(i)*csink(i)+Conc(j)*csink(j)+Conc(k)*csink(k)+Conc(l)*csink(l))/4._dp
                CumCh0=CumCh0-VEl*dt*(Gc(i)+Gc(j)+Gc(k)+Gc(l))/4._dp
                CumCh1=CumCh1-VEl*dt*((Fc(i)-sink(i))*Conc(i)+(Fc(j)-sink(j))*Conc(j)+&
                     (Fc(k)-sink(k))*Conc(k)+(Fc(l)-sink(l))*Conc(l))/4._dp-RootCh
                CumChR=CumChR+RootCh
             Endif
             FMul=VEl/5._dp
             GcE=(Gc(i)+Gc(j)+Gc(k)+Gc(l))/4._dp
             Ec1=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4._dp
             Ec2=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4._dp
             Ec3=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4._dp
             Ec4=(Dispxy(i)+Dispxy(j)+Dispxy(k)+Dispxy(l))/4._dp
             Ec5=(Dispxz(i)+Dispxz(j)+Dispxz(k)+Dispxz(l))/4._dp
             Ec6=(Dispyz(i)+Dispyz(j)+Dispyz(k)+Dispyz(l))/4._dp
             If (Lev.Eq.NLevel) AcE=(Ac(i)+Ac(j)+Ac(k)+Ac(l))/4._dp
             FcE=(Fc(i)+Fc(j)+Fc(k)+Fc(l))/4._dp
             SMul1=-1._dp/VEl/36._dp
             SMul2=VEl/30._dp
!!$                  if(lUpW) then
!!$                     W12=WeTab(1,iSE)
!!$                     W13=WeTab(2,iSE)
!!$                     W14=WeTab(3,iSE)
!!$                     W23=WeTab(4,iSE)
!!$                     W24=WeTab(5,iSE)
!!$                     W34=WeTab(6,iSE)
!!$                     A11=-2.*W12+2.*W14+2.*W13
!!$                     A12=-2.*W12+W14+W13
!!$                     A13=-W12+W14+2.*W13
!!$                     A14=-W12+2.*W14+W13
!!$                     A21=-W23+2.*W12+W24
!!$                     A22=-2.*W23+2.*W12+2.*W24
!!$                     A23=-2.*W23+W12+W24
!!$                     A24=-W23+W12+2.*W24
!!$                     A31=-W34+W23-2.*W13
!!$                     A32=-W34+2.*W23-W13
!!$                     A33=-2.*W34+2.*W23-2.*W13
!!$                     A34=-2.*W34+W23-W13
!!$                     A41=-2.*W14+W34-W24
!!$                     A42=-W14+W34-2.*W24
!!$                     A43=-W14+2.*W34-W24
!!$                     A44=-2.*W14+2.*W34-2.*W24
!!$                     Wx(1)=VxE(1)*A11+VxE(2)*A12+VxE(3)*A13+VxE(4)*A14
!!$                     Wx(2)=VxE(1)*A21+VxE(2)*A22+VxE(3)*A23+VxE(4)*A24
!!$                     Wx(3)=VxE(1)*A31+VxE(2)*A32+VxE(3)*A33+VxE(4)*A34
!!$                     Wx(4)=VxE(1)*A41+VxE(2)*A42+VxE(3)*A43+VxE(4)*A44
!!$                     Wy(1)=VyE(1)*A11+VyE(2)*A12+VyE(3)*A13+VyE(4)*A14
!!$                     Wy(2)=VyE(1)*A21+VyE(2)*A22+VyE(3)*A23+VyE(4)*A24
!!$                     Wy(3)=VyE(1)*A31+VyE(2)*A32+VyE(3)*A33+VyE(4)*A34
!!$                     Wy(4)=VyE(1)*A41+VyE(2)*A42+VyE(3)*A43+VyE(4)*A44
!!$                     Wz(1)=VzE(1)*A11+VzE(2)*A12+VzE(3)*A13+VzE(4)*A14
!!$                     Wz(2)=VzE(1)*A21+VzE(2)*A22+VzE(3)*A23+VzE(4)*A24
!!$                     Wz(3)=VzE(1)*A31+VzE(2)*A32+VzE(3)*A33+VzE(4)*A34
!!$                     Wz(4)=VzE(1)*A41+VzE(2)*A42+VzE(3)*A43+VzE(4)*A44
!!$           end if
             Do j1=1,4
                i1=List(j1)
                !F(i1)=F(i1)+FMul*(GcE+Gc(i1)/4._dp) + Vel/4._dp * csink_cube(iE)!!f_n
                F(i1)=F(i1)+FMul*(GcE+Gc(i1)/4._dp) !!f_n
                If(Lev.Eq.NLevel) DS(i1)=DS(i1)+FMul*(AcE+Ac(i1)/4._dp)!Q_nm
                S_m=0.0_dp
                Do j2=1,4
                   i2=List(j2)
                   S(j1,j2)=SMul1*(Ec1*bi(j1,iSE,iE)*bi(j2,iSe,iE)+Ec2*ci(j1,iSE,iE)*ci(j2,iSE,iE)+ Ec3*di(j1,iSe,iE)*di(j2,iSe,iE)+&
                        Ec4*(bi(j1,iSE,iE)*ci(j2,iSE,iE)+ci(j1,iSE,iE)*bi(j2,iSE,iE))+Ec5*(bi(j1,iSE,iE)*di(j2,iSE,iE)+di(j1,iSE,iE)*bi(j2,iSE,iE))+&
                        Ec6*(ci(j1,iSE,iE)*di(j2,iSE,iE)+di(j1,iSE,iE)*ci(j2,iSE,iE)))
                   S(j1,j2)=S(j1,j2)-(6*VEl/Deter(iSE,iE))*((bi(j2,iSE,iE)/30._dp)*(VxEE+VxE(j1)/4._dp)+&  !(6*VEl/Deter(iE,iSE))*
                        (ci(j2,iSE,iE)/30._dp)*(VyEE+VyE(j1)/4._dp)+(di(j2,iSE,iE)/30._dp)*(VzEE+VzE(j1)/4._dp))

                   !if(lUpW) S(j1,j2)=S(j1,j2)-(bi(iE,iSE,j2)/240.*Wx(j1)+ci(iE,iSE,j2)/240.*Wy(j1)+di(iE,iSE,j2)/240.*Wz(j1))

                   ic=1
                   If (i1.Eq.i2) ic=2
                   S(j1,j2)=S(j1,j2)+SMul2*ic*(FcE+(Fc(i1)+Fc(i2))/4._dp)
                   !S(j1,j2)=S(j1,j2)+ic*(SMul2*(FcE+(Fc(i1)+Fc(i2))/4._dp) + VEL/24._dp*sink_cube(iE))
                   S_m = S_m -(epsi*S(j1,j2))*Conc(i2)
                   If (Lev.Ne.NLevel) Then
                      Bs(i1)=Bs(i1)-alf*S(j1,j2)*Conc(i2)
                   Else
                      Call Find(i1,i2,kk,nPt,nband,IAD,IADN)
                      iB=kk
                      !iB=iadd_temp(k)
                      As(iB,i1)=As(iB,i1)+epsi*S(j1,j2)
                   Endif
                   If (Lev.Eq.1.And.Kode(i1).Gt.0) Qc(i1)=Qc(i1)-S(j1,j2)*Conc(i2)
                End Do
             End Do
          End Do !sub-element loop
       End Do  !element loop

       Do i=1,Npt
          M=MatNum(i)
          If (Lev.Eq.1.And.Kode(i).Gt.0) Qc(i)=Qc(i)-F(i)
          If (Lev.Ne.NLevel) Then
             Bs(i)=Bs(i)-alf*F(i)
          Else
             As(iadd(i),i)=As(iadd(i),i)+DS(i)/dt
             Bs(i)=Bs(i)+DS(i)/dt*Conc(i)-epsi*F(i)
          Endif
       End Do
    End Do

    ! Boundary condition
    Call C_Bound(KodCB,dt,DS,IADD)
    ! Solve the global matrix equation for transport
    If (epsi.Lt.0.001_dp) Then
       Do i=1,Npt
          Bs(i)=Bs(i)/As(iadd(i),i)
       End Do
    Else
       Call ILU (As,nPt,nband,IAD,IADN,IADD,A1)
       North=4!north=0 (symmetric matrix -> is not the case)
       Call OrthoMin(As,B1,Bs,nPt,nband,nPt,IAD,IADN,IADD,A1,North)
       !CALL SOLVET
    Endif
    Conc=B1
    Bs = B1
    Do i=1,Npt

       !sngl(B(i))
       If (Conc(i).Lt.1.0E-30_dp) Conc(i)=0.0_dp
    End Do
    Call SolInf(dt)
    ! ====================================
    Deallocate(As,A1,Bs)
    Deallocate(B1,F,DS)
    Deallocate(Ac,Fc,Gc,Qc)
    Deallocate(Dispxx,Dispyy,Dispzz)
    Deallocate(Dispxy,Dispxz,Dispyz)
    Return
  End Subroutine SOLUTE
  !*******************************************************************
  Subroutine C_Bound(KodCB,dt,DS,IADD)
    Use typedef
    Use BoundData
    Use GridData, Only: nBCPts
    Use SolData, Only: epsi,Kode,conc
    Use CumData
    Use MatData
    Implicit None

    Integer(ap),Intent(inout) ::KodCB(nPt)
    Integer(ap) :: cKod,i,j
    Integer(ap) :: IADD(nPt)
    Real(dp),Intent(in) :: dt,DS(nPt)
    Real(dp) :: alf,cBnd

    alf=1._dp-epsi
    Do i=1,nPt
       If (Kode(i).Ne.0) Then
          cKod=0
          Do j=1,nBCPts
             If (iBCPt(j).Eq.i) Then

                If (KodCB(j).Gt.0) Then
                   cKod=1
                   cBnd=cBound(KodCB(j))
                Else
                   If (Q(i).Gt.0._dp) Then
                      cKod=3
                      cBnd=cBound(-KodCB(j))
                   Else
                      cKod=2
                   Endif
                Endif
                Goto 12
             Endif
          End Do
12        Continue
          If (cKod.Eq.1) Then
             ! Dirichlet boundary condition
             Qc(i)=Qc(i)+Q(i)*(epsi*cBnd+alf*Conc(i))-DS(i)*(cBnd-Conc(i))/dt
             As(iadd(i),i)=1._dp
             Bs(i)=1._dp*cBnd

!!$             DO 13 j=1,2*NBand-1
!!$                  As(j,i)=0.0_dp
!!$13           CONTINUE
!!$             As(NBand,i)=1._dp
!!$             Bs(i)=cBnd

          Elseif(cKod.Eq.2) Then
             ! Neumann boundary condition
             Qc(i)=Q(i)*Conc(i)

          Elseif (cKod.Eq.3) Then
             ! Cauchy boundary condition
             Bs(i)=Bs(i)-Q(i)*(cBnd-alf*Conc(i))
             As(iadd(i),i)=As(iadd(i),i)-epsi*Q(i)
             Qc(i)=Q(i)*cBnd
          Endif
       Endif
    End Do
    Return
  End Subroutine C_Bound
  !**********************************************************************************
  Subroutine ChInit(dtMaxC,dt)
    Use typedef
    Use GridData
    Use SolData
    Use WatFun
    Implicit None

    Integer(ap) ::  i,M
    Real(dp),Intent(in) :: dt
    Real(dp),Intent(inout) :: dtMaxC

    If (.Not.Allocated(Fc)) Allocate(Fc(nPt),Gc(nPt))

    Do i=1,nPt
       If (NLevel.Eq.2) Then
          M=MatNum(i)
          Gc(i)=ChPar(8,M)*theta(i)+ChPar(1,M)*ChPar(9,M)
          Fc(i)=ChPar(6,M)*theta(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)+(sink(i)-csink(i))
       Endif
    End Do

    Call Veloc
    Call Disper
    Call PeCour(dtMaxC,dt)
    Return
  End Subroutine ChInit
  !*************************************************************************
  Subroutine Veloc
    Use typedef
    Use GridData
    Use SolData
    Use CumData
    Implicit None

    Integer(ap) :: i,j,k,l,m,ise,ie,List(4)
    Real(dp) :: vxx,vyy,vzz,cayz,caxz,caxy,cazz,cayy,caxx, VE, A

    Do i=1,nPt
       Vx(i)=0.0_dp
       Vy(i)=0.0_dp
       Vz(i)=0.0_dp
    End Do
    Do iE=1,nElm
       CAxx=ConAxx(iE)
       CAyy=ConAyy(iE)
       CAzz=ConAzz(iE)
       CAxy=ConAxy(iE)
       CAxz=ConAxz(iE)
       CAyz=ConAyz(iE)
       Do iSE=1,5
          i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          List(1)=i
          List(2)=j
          List(3)=k
          List(4)=l
          VE=Abs(Deter(iSE,iE))/6._dp
          A=1./VE/6
          Vxx=A*(Ax(1,iSE,iE)*hNew(i)+Ax(2,iSE,iE)*hNew(j)+Ax(3,iSE,iE)*hNew(k)+Ax(4,iSE,iE)*hNew(l))!Deter, Ax, Ay and Az are now globals (Couvreur mar 2010)
          Vxx=Vxx+CAxz
          Vyy=A*(Ay(1,iSE,iE)*hNew(i)+Ay(2,iSE,iE)*hNew(j)+Ay(3,iSE,iE)*hNew(k)+Ay(4,iSE,iE)*hNew(l))
          Vyy=Vyy+CAyz
          Vzz=A*(Az(1,iSE,iE)*hNew(i)+Az(2,iSE,iE)*hNew(j)+Az(3,iSE,iE)*hNew(k)+Az(4,iSE,iE)*hNew(l))
          Vzz=Vzz+CAzz
          Do m=1,4
             l=List(m)
             Vx(l)=Vx(l)-Con(l)*Vxx
             Vy(l)=Vy(l)-Con(l)*Vyy
             Vz(l)=Vz(l)-Con(l)*Vzz
          End Do
       End Do
    End Do
    Do i=1,nPt
       Vx(i)=Vx(i)/ListNE(i)
       Vy(i)=Vy(i)/ListNE(i)
       Vz(i)=Vz(i)/ListNE(i)
    End Do
    Return
  End Subroutine Veloc
  !****************************************************************
  Subroutine Disper
    Use typedef
    Use GridData
    Use SolData
    Use WatFun
    Implicit None

    Real(dp) :: Tau,Vabs,ThSati
    Integer(ap) :: i,M

    If (.Not.Allocated(Dispxx)) Then
       Allocate(Dispxx(nPt),Dispyy(nPt),Dispzz(nPt))
       Allocate(Dispxy(nPt),Dispxz(nPt),Dispyz(nPt))
    Endif
    Do i=1,nPt
       M=MatNum(i)
       If (soiltab) Then!(Couvreur nov 2011)
          ThSati=(TheTab(1,M)-TheTab(nTab,M))*Dxy(i)+TheTab(nTab,M)*Exy(i)
       Else
          ThSati=(par(3,M)-par(2,M))*Dxy(i)+par(2,M)*Exy(i)
       Endif
       Tau=theta(i)**(7._dp/3._dp)/ThSati**2
       Vabs=Sqrt(Vx(i)*Vx(i)+Vy(i)*Vy(i)+Vz(i)*Vz(i))
       If (Vabs.Gt.0._dp) Then
          If (soiltab) Then!(Couvreur nov 2011)
             Dispxx(i)=ChPar(3,M)*Vx(i)*Vx(i)/Vabs+ChPar(4,M)*(Vz(i)*Vz(i)+Vy(i)*Vy(i))/Vabs+Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispyy(i)=ChPar(3,M)*Vy(i)*Vy(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vz(i)*Vz(i))/Vabs+Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispzz(i)=ChPar(3,M)*Vz(i)*Vz(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vy(i)*Vy(i))/Vabs+Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
          Else
             Dispxx(i)=ChPar(3,M)*Vx(i)*Vx(i)/Vabs+ChPar(4,M)*(Vz(i)*Vz(i)+Vy(i)*Vy(i))/Vabs+Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
             Dispyy(i)=ChPar(3,M)*Vy(i)*Vy(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vz(i)*Vz(i))/Vabs+Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
             Dispzz(i)=ChPar(3,M)*Vz(i)*Vz(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vy(i)*Vy(i))/Vabs+Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Endif
          Dispxy(i)=(ChPar(3,M)-ChPar(4,M))*Vx(i)*Vy(i)/Vabs
          Dispxz(i)=(ChPar(3,M)-ChPar(4,M))*Vx(i)*Vz(i)/Vabs
          Dispyz(i)=(ChPar(3,M)-ChPar(4,M))*Vy(i)*Vz(i)/Vabs
       Else
          If (soiltab) Then!(Couvreur nov 2011)
             Dispxx(i)=Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispyy(i)=Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispzz(i)=Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
          Else
             Dispxx(i)=Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
             Dispyy(i)=Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
             Dispzz(i)=Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Endif
          Dispxy(i)=0.0_dp
          Dispxz(i)=0.0_dp
          Dispyz(i)=0.0_dp
       Endif
    End Do
  End Subroutine Disper
  !****************************************************************
!!$      SUBROUTINE SolveT
!!$      USE GridData
!!$      USE MatData, ONLY: As, Bs
!!$      IMPLICIT NONE
!!$      REAL(dp) :: P,C,Sum
!!$      INTEGER(ap):: k,N1,i,kk,kc,ii,j,L,jj,M
!!$      N1=Npt-1
!!$      DO 12 k=1,N1
!!$         P=1._dp/As(NBand,k)
!!$         kk=k+1
!!$         kc=NBand
!!$         DO 11 i=kk,Npt
!!$            kc=kc-1
!!$            IF (kc.LE.0) GOTO 12
!!$            C=-P*As(kc,i)
!!$            As(kc,i)=C
!!$            ii=kc+1
!!$            L=kc+NBand-1
!!$            DO 11 j=ii,L
!!$             jj=j+NBand-kc
!!$             As(j,i)=As(j,i)+C*As(jj,k)
!!$11       CONTINUE
!!$12    CONTINUE
!!$      DO 14 i=2,Npt
!!$       jj=NBand+1-i
!!$       ii=1
!!$       IF (jj.LE.0) THEN
!!$          jj=1
!!$          ii=i-NBand+1
!!$       ENDIF
!!$       Sum=0.0_dp
!!$       DO 13 j=jj,NBand-1
!!$          Sum=Sum+As(j,i)*Bs(ii)
!!$          ii=ii+1
!!$13     CONTINUE
!!$       Bs(i)=Bs(i)+Sum
!!$14    CONTINUE
!!$      Bs(Npt)=Bs(Npt)/As(NBand,Npt)
!!$      DO 16 k=1,N1
!!$       i=Npt-k
!!$       jj=i
!!$       m=MIN(2*NBand-1,NBand+k)
!!$       Sum=0.0_dp
!!$       DO 15 j=NBand+1,m
!!$          jj=jj+1
!!$          Sum=Sum+As(j,i)*Bs(jj)
!!$15     CONTINUE
!!$       Bs(i)=(Bs(i)-Sum)/As(NBand,i)
!!$16    CONTINUE
!!$      RETURN
!!$      END SUBROUTINE SolveT
  !****************************************************************
  Subroutine PeCour(dtMaxC,dt)
    Use typedef
    Use GridData
    Use SolData
    Use WatFun
    Implicit None

    Integer(ap)::i,j,k,l,ie,ise
    Real(dp), Intent(in) :: dt
    Real(dp), Intent(out):: dtMaxC
    Real(dp):: pecx,pecy,delx !,xmax,xmin,zmin,zmax,ymin,ymax
    Real(dp):: r1,r2,r3,r4,rmin
    Real(dp):: vzmax,vymax,vxmax,vze,vye,vxe,dze,dye,dxe,delz,dely
    Real(dp):: courx,coury,courz,cour1,cour2,cour3,dt3,dt2,dt1,pecz

    Peclet=0.0_dp
    Courant=0.0_dp
    dtMaxC=1.e+30_dp
    Do iE=1,Nelm
       Do iSE=1,5
          PecX=99999._dp
          PecY=99999._dp
          PecZ=99999._dp
          dt1=1.e+30_dp
          dt2=1.e+30_dp
          dt3=1.e+30_dp

          i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          !            thetai=Fth(hNew(i),par(:,MatNum(i)))
          !            thetaj=Fth(hNew(j),par(:,MatNum(j)))
          !           thetak=Fth(hNew(k),par(:,MatNum(k)))
          !           thetal=Fth(hNew(l),par(:,MatNum(l)))
          !!JAVAUX AMIN1 and AMAX! have been changed to min and max
          !          xmax=max(xGrid(i),xGrid(j),xGrid(k),xGrid(l))
          !          xmin=min(xGrid(i),xGrid(j),xGrid(k),xGrid(l))
          !          ymax=max(yGrid(i),yGrid(j),yGrid(k),yGrid(l))
          !          ymin=min(yGrid(i),yGrid(j),yGrid(k),yGrid(l))
          !          zmax=max(zGrid(i),zGrid(j),zGrid(k),zGrid(l))
          !          zmin=min(zGrid(i),zGrid(j),zGrid(k),zGrid(l))
          !            delX=xmax-xmin
          !            delY=ymax-ymin
          !            delZ=zmax-zmin
          delX=dxGrid !Voxel size is constant and using xmax, xmin, etc causes problems in continuous (Couvreur oct 2010)
          delY=dyGrid
          delZ=dzGrid
          DxE=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4
          DyE=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4
          DzE=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4
          VxE=Abs(Vx(i)+Vx(j)+Vx(k)+Vx(l))/4
          VyE=Abs(Vy(i)+Vy(j)+Vy(k)+Vy(l))/4
          VzE=Abs(Vz(i)+Vz(j)+Vz(k)+Vz(l))/4
          If (DxE.Gt.0._dp) PecX=VxE*delX/DxE
          If (DyE.Gt.0._dp) PecY=VyE*delY/DyE
          If (DzE.Gt.0._dp) PecZ=VzE*delZ/DzE
          !JAVAUX AMIN1 and AMAX! have been changed to min and max
          If (PecX.Ne.99999._dp) Peclet=Max(Peclet,PecX)
          If (PecY.Ne.99999._dp) Peclet=Max(Peclet,PecY)
          If (PecZ.Ne.99999._dp) Peclet=Max(Peclet,PecZ)
          !JAVAUX AMIN1 and AMAX! have been changed to min and max
          Peclet=Min(Peclet,99999._dp)
          VxMax=Max(Abs(Vx(i))/theta(i),Abs(Vx(j))/theta(j),Abs(Vx(k))/theta(k),Abs(Vx(l))/theta(l))
          VyMax=Max(Abs(Vy(i))/theta(i),Abs(Vy(j))/theta(j),Abs(Vy(k))/theta(k),Abs(Vy(l))/theta(l))
          VzMax=Max(Abs(Vz(i))/theta(i),Abs(Vz(j))/theta(j),Abs(Vz(k))/theta(k),Abs(Vz(l))/theta(l))
          R1=1._dp+ChPar(1,MatNum(i))*ChPar(5,MatNum(i))/theta(i)
          R2=1._dp+ChPar(1,MatNum(j))*ChPar(5,MatNum(j))/theta(j)
          R3=1._dp+ChPar(1,MatNum(k))*ChPar(5,MatNum(k))/theta(k)
          R4=1._dp+ChPar(1,MatNum(l))*ChPar(5,MatNum(l))/theta(l)
          RMin=Min(R1,R2,R3,R4)
          CourX=VxMax*dt/(delX*RMin)
          CourY=VyMax*dt/(delY*RMin)
          CourZ=VzMax*dt/(delZ*RMin)
          !JAVAUX AMIN1 and AMAX! have been changed to min and max
          Courant=Max(Courant,CourX,CourY,CourZ)
          Cour1=1.0_dp
          Cour2=1.0_dp
          Cour3=1.0_dp
          !JAVAUX AMIN1 and AMAX! have been changed to min and max
          If(PecX.Ne.99999._dp) Cour1=Min(1._dp,PeCr/Max(0.5_dp,PecX))
          If(PecY.Ne.99999._dp) Cour2=Min(1._dp,PeCr/Max(0.5_dp,PecY))
          If(PecZ.Ne.99999._dp) Cour3=Min(1._dp,PeCr/Max(0.5_dp,PecZ))
          If(VxMax.Gt.0._dp) dt1=Cour1*delX*RMin/VxMax
          If(VyMax.Gt.0._dp) dt2=Cour2*delY*RMin/VyMax
          If(VzMax.Gt.0._dp) dt3=Cour3*delZ*RMin/VzMax
          !JAVAUX AMIN1 and AMAX! have been changed to min and max
          dtMaxC=Min(dtMaxC,dt1,dt2,dt3)
       End Do
    End Do
  End Subroutine PeCour
  !***********************************************************************************
  Subroutine SolInf(dt)
    Use typedef
    Use GridData
    Use CumData
    Use SolData, Only: Kode
    Implicit None

    Integer(ap):: i,j
    Real(dp):: sMean(2),dt

    sMean(1)=0.0_dp
    sMean(2)=0.0_dp
    Do i=1,nPt
       j=iabs(Kode(i))
       If (j.Ne.0) Then
          sMean(j)=sMean(j)-Qc(i) !Qc=nodal value of solute flux, only boundary values!!
       Endif
    End Do
    cCumA=Abs(CumCh0)+Abs(CumCh1)+Abs(CumChR)
    cCumT=CumCh0+CumCh1+CumChR
    Do j=1,2
       ChemS(j)=ChemS(j)+sMean(j)*dt
       cCumT=cCumT+ChemS(j)
       cCumA=cCumA+Abs(ChemS(j))
    End Do
  End Subroutine SolInf

End Module SoluteMod
