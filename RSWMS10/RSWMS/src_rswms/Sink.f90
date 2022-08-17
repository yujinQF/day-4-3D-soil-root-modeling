!> \file Sink.f90
!! \brief Source file SINK.

!> Module Sink
Module Sink
  Implicit None

Contains

  Subroutine BetDis(t,ipl)
    Use typedef
    Use RootData
    Use PlntData
    Use GridData
    Use Doussan, Only: intsec
    Use RootGrowthNeighb, Only: Neighb
    Implicit None

    Integer(ap):: corner(8)
    Integer(ap)::i,ibr,iprv,igrow,iseg,iorder,iurf,iface,ic,irec,imin,ipl,isub,typ
    Real(dp):: betcon(8)
    Real(dp):: x1,x2,xA,xB,y1,y2,yA,yB,z1,z2,zA,zB
    Real(dp):: xInt,yInt,zInt,xCent,yCent,zCent
    Real(dp)::dis,sqrdis,sumB,weight,segage,srface
    Real(dp):: totlen,prvsur,t 
    Real(dp):: blengt
    Logical split

    PrvSur=0.0_dp
    TotSur=0.0_dp
    TotLen=0.0_dp
    split=.False.
    Do iseg=1,nrec(ipl)
       TotLen=TotLen+seglen(iseg,ipl)
       PrvSur=PrvSur+segsur(iseg,ipl)
    End Do
    If (TotLen.Lt.1.E-20_dp) TotLen=1.E-20_dp
    If (PrvSur.Lt.1.E-20_dp) PrvSur=1.E-20_dp
    Do i=1,nPt
       betaw(i)=0.0_dp
       betac(i)=0.0_dp
    End Do
    !> go through each root segment and update the node surface function:
2   Do  ibr=1,nbr(ipl)
       !> find the tip segment of the branch 'ibrnch'
       irec=nrec(ipl)+1
11     irec=irec-1
       If (ibrseg(irec,ipl).Ne.ibr) Goto 11
       If (seglen(irec,ipl).Lt.1.E-20_dp) Then
          iprv=irecpr(irec,ipl)
          If (iprv.Eq.0) Goto 2
          Goto 101
       Else
          Do  igrow=1,ngrow(ipl)
             If (ibrgrw(igrow,ipl).Eq.ibr) Then
                iprv=irec
                If (iprv.Eq.0) Goto 2
                xA=xg(igrow,ipl)
                yA=yg(igrow,ipl)
                zA=zg(igrow,ipl)
                Goto 102
             Endif
          End Do
       Endif
100    iprv=irecpr(irec,ipl)
       If (iprv.Eq.0) Goto 2
101    If (ibrseg(iprv,ipl).Eq.ibrseg(irec,ipl)) Then
          xA=xs(irec,ipl)+xplant(ipl)
          yA=ys(irec,ipl)+yplant(ipl)
          zA=zs(irec,ipl)
       Else
          Goto 2
       Endif
       !> find cuboid around segment:
102    Call Neighb(xA,yA,zA,corner,imin,.FALSE.)
       x1=xgrid(corner(1))
       x2=xgrid(corner(4))
       y1=ygrid(corner(1))
       y2=ygrid(corner(4))
       z1=zgrid(corner(1))
       z2=zgrid(corner(5))
       !> find the other end of the segment:
       xB=(xs(iprv,ipl))+xplant(ipl)
       yB=(ys(iprv,ipl))+yplant(ipl)
       zB=(zs(iprv,ipl))
       !> calculate segment surface:
       blengt=seglen(iprv,ipl)
       srface=segsur(iprv,ipl)
       TotSur=TotSur+segsur(iprv,ipl)
       !> calculate segment weighing factor according to age:
       iorder=0
       iorder=ordseg(iprv,ipl)
       If (maizeroottyp) Then
          If (iorder.Lt.12) Then  !> In RootTyp, maize root types below 12 are principal roots (12 and 13 are lateral)
             typ=1
          Else
             typ=2
          Endif
       Elseif (loliumroottyp) Then
          If (iorder.Lt.3) Then  !> In RootTyp, lolium root types below 3 are principal roots
             typ=1
          Elseif (iorder.Lt.5) Then
             typ=2
          Else
             typ=3
          Endif
       Elseif (wheatroottyp) Then
          If (iorder.Lt.19) Then  !> In RootTyp, wheat root types below 19 are principal roots
             typ=1
          Elseif (iorder.Lt.20) Then
             typ=2
          Else
             typ=3
          Endif
       Else
          If (iorder.Eq.0) Then  !> type 0 is the seed+small segment
             typ=1
          Elseif (iorder.Gt.3) Then !> no more than 3 root types
             typ=3
          Else
             typ=iorder
          Endif
       Endif
       segage=t-timorg(iprv,ipl)
       If (segage.Ge.0.0_dp.And.lUrf) Then
          If (segage.Ge.age(typ,nUrf(typ))) Then
             Weight=Urf(typ,nUrf(typ))
          Else
             iUrf=nUrf(typ)
4            iUrf=iUrf-1
             If ((segage.Lt.age(typ,iUrf)).And.(iUrf.Gt.1)) Goto 4
             Weight=Urf(typ,iUrf)+(segage-age(typ,iUrf))/&
                  (age(typ,iUrf+1)-age(typ,iUrf))*&
                  (Urf(typ,iUrf+1)-Urf(typ,iUrf))
          Endif
       Else
          Weight=1.0_dp
       Endif
       If (intsec(xA,yA,zA,xB,yB,zB,irec,iprv,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) Then
          split=.True.
          ! calculate subsegment length and surface...
41        blengt=Sqrt((xInt-xA)*(xInt-xA)+(yInt-yA)*(yInt-yA)+(zInt-zA)*(zInt-zA))
          srface=blengt*segsur(iprv,ipl)/seglen(iprv,ipl)
          ! calculate subsegment center coordinates...
          xCent=xA+(xInt-xA)/2._dp
          yCent=yA+(yInt-yA)/2._dp
          zCent=zA+(zInt-zA)/2._dp
          ! and calculate the distribution for each node:
          sumB=0.0_dp
          Do ic=1,8
             sqrdis=(xCent-xgrid(corner(ic)))*(xCent-xgrid(corner(ic)))+ &
                  (yCent-ygrid(corner(ic)))*(yCent-ygrid(corner(ic)))+ &
                  (zCent-zgrid(corner(ic)))*(zCent-zgrid(corner(ic)))
             dis=Sqrt(sqrdis)
             If (dis.Lt.1.E-20_dp) dis=1.E-20_dp
             betcon(ic)=1._dp/dis
             sumB=sumB+betcon(ic)
          End Do
          Do ic=1,8
             betaw(corner(ic))=betaw(corner(ic))+betcon(ic)/sumB*Weight*(blengt/TotLen)
             betac(corner(ic))=betac(corner(ic))+betcon(ic)/sumB*Weight*(srface/PrvSur)
          End Do
          xA=xInt
          yA=yInt
          zA=zInt
          If(iFace.Eq.1) zInt=zInt-1.E-5_dp*dzGrid
          If(iFace.Eq.2) zInt=zInt+1.E-5_dp*dzGrid
          If(iFace.Eq.3) yInt=yInt-1.E-5_dp*dyGrid
          If(iFace.Eq.4) yInt=yInt+1.E-5_dp*dyGrid
          If(iFace.Eq.5) xInt=xInt-1.E-5_dp*dxGrid
          If(iFace.Eq.6) xInt=xInt+1.E-5_dp*dxGrid
          Call Neighb(xInt,yInt,zInt,corner,imin,.FALSE.)
          !calculate cuboid´s corners coordinates:
          x1=xgrid(corner(1))
          x2=xgrid(corner(4))
          y1=ygrid(corner(1))
          y2=ygrid(corner(4))
          z1=zgrid(corner(1))
          z2=zgrid(corner(5))
          If (intsec(xA,yA,zA,xB,yB,zB,irec,iprv,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) Then
             Goto 41
          Else
             ! calculate subsegment length and surface...
             split=.False.
             blengt=Sqrt((xB-xA)*(xB-xA)+(yB-yA)*(yB-yA)+(zB-zA)*(zB-zA))
             srface=blengt*segsur(iprv,ipl)/seglen(iprv,ipl)
          Endif
       Endif
       !calculate segment center coordinates:
       xCent=xA+(xB-xA)/2._dp
       yCent=yA+(yB-yA)/2._dp
       zCent=zA+(zB-zA)/2._dp
       !calculate the distribution for each node:
       sumB=0.0
       Do ic=1,8
          sqrdis=(xCent-xgrid(corner(ic)))*(xCent-xgrid(corner(ic)))+&
               (yCent-ygrid(corner(ic)))*(yCent-ygrid(corner(ic)))+&
               (zCent-zgrid(corner(ic)))*(zCent-zgrid(corner(ic)))
          dis=Sqrt(sqrdis)
          If (dis.Lt.1.E-20_dp) dis=1.E-20_dp
          betcon(ic)=1._dp/dis
          sumB=sumB+betcon(ic)
       End Do
       Do ic=1,8
          betaw(corner(ic))=betaw(corner(ic))+betcon(ic)/sumB*Weight*(blengt/TotLen)
          betac(corner(ic))=betac(corner(ic))+betcon(ic)/sumB*Weight*(srface/PrvSur)
       End Do
       irec=iprv
       Goto 100
    End Do
    Return
  End Subroutine BetDis
  !********************************************************************************************
  Subroutine BetNrm
    Use typedef
    Use PlntDATA
    Use RootData, Only: OmegaC,lJarvis,lUrf
    Use DoussanMat, Only: stresfun
    Use GridData
    Implicit None

    Integer(ap):: i,j,k,l,iE,ise,kOut,ipl
    Real(dp):: betwe,betce,VEl,Sbetac,Sbetaw
    Character :: file*17

    kout=0
    ipl=1
    Sbetaw=0.0_dp
    Sbetac=0.0_dp

    Do iE=1,nElm
       Do iSE=1,5

          i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)

          betwE=(betaw(i)+betaw(j)+betaw(k)+betaw(l))/4.
          Sbetaw=Sbetaw+betwE
          betcE=(betac(i)+betac(j)+betac(k)+betac(l))/4.
          Sbetac=Sbetac+betcE
       End Do
    End Do
    VEl=dxGrid*dyGrid*dzGrid/5._dp
    If (Sbetaw.Lt.1.E-20_dp) Return
    Sbetaw=Sbetaw*VEl
    Do  i=1,nPt
       betaw(i)=betaw(i)/Sbetaw
    End Do

    Write (file,'(A10)')'out/Feddes'
    Write (file(11:11),'(I1)') ipl
    Write (file(12:12),'(A1)') '.'
    If (kout.Lt.10) Then
       Write (file(13:13),'(I1)') kout
    Elseif (kout.Lt.100) Then
       Write (file(13:14),'(I2)') kout
    Else
       Write (file(13:15),'(I3)') kout
    Endif
    Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
    Write (8,'(''********** FEDDES ROOT WATER UPTAKE MODEL PARMS **********'')')
    Write (8,*)
    Write (8,'(''Stress function ? (0=no,1=Feddes,2=van Genuchten)'')')
    Write (8,'(I1)') stresfun
    Write (8,*)
    Write (8,'(''Feddes stress function:'')')
    Write (8,'(''h0      h1      h2      h3 [hPa]'')')
    Write (8,'(F5.1,2X,F5.1,2X,F8.1,2X,F9.1)') h0,h1,h2,h3
    Write (8,*)
    Write (8,'(''van Genuchten stress function:'')')
    Write (8,'(''p50[hPa]   h50[hPa]   p1[-]  p2[-]'')')
    Write (8,'(F8.1,2X,F8.1,2X,F5.1,2X,F5.1)') p50,h50,p1,p2
    Write (8,*)
    Write (8,'(''Use of Jarvis (1989) function for compensatory RWU prediction ? (1=yes,2=no)'')')
    If (lJarvis) Then
       Write (8,'(I1)') 1
    Else
       Write (8,'(I1)') 2
    Endif
    Write (8,*)
    Write (8,'(''Omega_c value [-] (1=uncompensated RWU,0=fully compensated RWU)'')')
    Write (8,'(F4.2)') OmegaC
    Write (8,*)
    Write (8,'(''Uptake Reduction Function (URF) [-] with segment age? Only if rRLD not given beside (1=yes,2=no)'')')
    If (lUrf) Then
       Write (8,'(I1)') 1
    Else
       Write (8,'(I1)') 2
    Endif
    Write (8,*)
    Write (8,'(''* Root Length distribution of the root system'')')
    Write (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
    Write (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
    Write (8,'(''    Node#    betaw(RLD)'')')
    Do i=1,nPt
       Write (8,'(I6,3X,F12.10)') i,betaw(nPt+1-i)!nodes order inversed for Hydrus
    Enddo
    Close (8)
    Return
  End Subroutine BetNrm
 !********************************************************************************************
 !> ### root length density distribution ###
  Subroutine RLDdis(t,kout,ipl)
    Use typedef
    Use ParamData, Only: pi
    Use RootData
    Use PlntData
    Use GridData, Only: RLD,RSD,continu,dxGrid,dyGrid,dzGrid,nElm,nex,ney,nez
    Use DoussanMat, Only: transroot,nsub,Intc,cube_i,stresfun
    Implicit None

    Integer(ap):: i,ifol,ibr,igrow,iorder,iurf,irec,ipl,isub,typ,transrootA(2)
    Integer(ap), Intent(in)::kout
    Real(dp):: xA,xB,yA,yB,zA,zB
    Real(dp):: Weight,segage
    Real(dp):: t
    Real(dp):: blengt,srface,TotLen
    Character :: file*17

    Write (*,*) '... Calculating RLD ...'
    TotSur=0.0_dp
    TotLen=0.0_dp
    ipl=1
    Allocate (RLD(1:nElm))
    RLD(1:nElm)=0.0_dp
    Allocate (RSD(1:nElm))
    RSD(1:nElm)=0.0_dp
    !> go through each root segment and update the node surface function:
2   Do ibr=1,nbr(ipl)!> Loop on root branches
       !> find the tip segment of the branch 'ibrnch'
       irec=nrec(ipl)+1
11     irec=irec-1
       If (ibrseg(irec,ipl).Ne.ibr) Goto 11
       Do igrow=1,ngrow(ipl)
          If (ibrgrw(igrow,ipl).Eq.ibr) Then
             ifol=nrec(ipl)+igrow
             If (irec.Eq.0) Goto 2!"if root sys top, next branch!"
             xA=xg(igrow,ipl)+xplant(ipl)
             yA=yg(igrow,ipl)+yplant(ipl)
             zA=zg(igrow,ipl)
             If (continu) Then
                transrootA(1:2)=transroot(nrec(ipl)+igrow,1:2,nsub(nrec(ipl)+igrow,ipl),ipl)
                xA=xA+transrootA(1)*(nex*dxgrid)
                yA=yA+transrootA(2)*(ney*dygrid)
             Endif
             Goto 102
          Endif
       End Do
100    irec=irecpr(ifol,ipl)!Loop on root segments
       If (irec.Eq.0) Goto 2!"if root sys top, next branch!"
       If (ibrseg(irec,ipl).Eq.ibrseg(ifol,ipl)) Then!"on the same brench?"
          xA=Intc(ifol,1,nsub(ifol,ipl),ipl)
          yA=Intc(ifol,2,nsub(ifol,ipl),ipl)
          zA=Intc(ifol,3,nsub(ifol,ipl),ipl)
       Else
          Goto 2!"if branch top, next branch!"By doing that, we commit the error of ommitting the length of the first part of all branches
       Endif
102    If (lUrf) Then
          ! calculate segment weighing factor according to age:
          iorder=0
          iorder=ordseg(irec,ipl)
          If (maizeroottyp) Then
             If (iorder.Lt.12) Then  !In RootTyp, maize root types below 12 are principal roots (12 and 13 are lateral)
                typ=1
             Else
                typ=2
             Endif
          Elseif (loliumroottyp) Then
             If (iorder.Lt.3) Then  !In RootTyp, lolium root types below 3 are principal roots
                typ=1
             Elseif (iorder.Lt.5) Then
                typ=2
             Else
                typ=3
             Endif
          Elseif (wheatroottyp) Then
             If (iorder.Lt.19) Then  !In RootTyp, wheat root types below 19 are principal roots
                typ=1
             Elseif (iorder.Lt.20) Then
                typ=2
             Else
                typ=3
             Endif
          Else
             If (iorder.Eq.0) Then  !type 0 is the seed+small segment
                typ=1
             Elseif (iorder.Gt.3) Then!no more than 3 root types
                typ=3
             Else
                typ=iorder
             Endif
          Endif
          segage=t-timorg(irec,ipl)
          If (segage.Ge.0.0_dp) Then
             If (segage.Ge.age(typ,nUrf(typ))) Then
                Weight=Urf(typ,nUrf(typ))
             Else
                iUrf=nUrf(typ)
4               iUrf=iUrf-1
                If ((segage.Lt.age(typ,iUrf)).And.(iUrf.Gt.1)) Goto 4
                Weight=Urf(typ,iUrf)+(segage-age(typ,iUrf))/&
                     (age(typ,iUrf+1)-age(typ,iUrf))*&
                     (Urf(typ,iUrf+1)-Urf(typ,iUrf))
             Endif
          Else
             Weight=1.0_dp
          Endif
       Else
          Weight=1.0_dp
       Endif
       Do isub=1,nsub(irec,ipl)!Loop on sub-segments
          xB=Intc(irec,1,isub,ipl)
          yB=Intc(irec,2,isub,ipl)
          zB=Intc(irec,3,isub,ipl)
          ! calculate subsegment length and surface...
          If (continu) Then
             blengt=Sqrt((xB-xA+(transrootA(1)-transroot(irec,1,isub,ipl))*nex*dxGrid)**2+(yB-yA+(transrootA(2)-transroot(irec,2,isub,ipl))*ney*dyGrid)**2+(zB-zA)**2)
          Else
             blengt=Sqrt((xB-xA)**2+(yB-yA)**2+(zB-zA)**2)
          Endif
          i=0
42        If (seglen(irec,ipl).Gt.0.) Then
             srface=blengt*segsur(irec+i,ipl)/seglen(irec+i,ipl)
          Else
             i=i-1!Temporary solution for segment whose saved length is 0
             Goto 42
          Endif
          RLD(cube_i(irec,isub,ipl))=RLD(cube_i(irec,isub,ipl))+blengt*Weight
          RSD(cube_i(irec,isub,ipl))=RSD(cube_i(irec,isub,ipl))+srface*Weight
          TotLen=TotLen+blengt*Weight
          TotSur=TotSur+srface*Weight
          xA=xB
          yA=yB
          zA=zB
          transrootA(1:2)=transroot(irec,1:2,isub,ipl)
       End Do!Loop on sub-segments
       ifol=irec
       Goto 100!Loop on root segments
    End Do!Loop on root branches
    RLD=RLD/(dxGrid*dyGrid*dzGrid)!/TotLen
    RSD=RSD/(dxGrid*dyGrid*dzGrid)!/TotSur

    If (lFed) Then
       Write (file,'(A10)')'out/Feddes'
       Write (file(11:11),'(I1)') ipl
       Write (file(12:12),'(A1)') '.'
       If (kout.Lt.10) Then
          Write (file(13:13),'(I1)') kout
       Elseif (kout.Lt.100) Then
          Write (file(13:14),'(I2)') kout
       Else
          Write (file(13:15),'(I3)') kout
       Endif
       Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
       Write (8,'(''********** FEDDES ROOT WATER UPTAKE MODEL PARMS **********'')')
       Write (8,*)
       Write (8,'(''Stress function ? (0=no,1=Feddes,2=van Genuchten)'')')
       Write (8,'(I1)') stresfun
       Write (8,*)
       Write (8,'(''Feddes stress function:'')')
       Write (8,'(''h0      h1      h2      h3 [hPa]'')')
       Write (8,'(F5.1,2X,F5.1,2X,F8.1,2X,F9.1)') h0,h1,h2,h3
       Write (8,*)
       Write (8,'(''van Genuchten stress function:'')')
       Write (8,'(''p50[hPa]   h50[hPa]   p1[-]  p2[-]'')')
       Write (8,'(F8.1,2X,F8.1,2X,F5.1,2X,F5.1)') p50,h50,p1,p2
       Write (8,*)
       Write (8,'(''Use of Jarvis (1989) function for compensatory RWU prediction ? (1=yes,2=no)'')')
       If (lJarvis) Then
          Write (8,'(I1)') 1
       Else
          Write (8,'(I1)') 2
       Endif
       Write (8,*)
       Write (8,'(''Omega_c value [-] (1=uncompensated RWU,0=fully compensated RWU)'')')
       Write (8,'(F4.2)') OmegaC
       Write (8,*)
       Write (8,'(''Uptake Reduction Function (URF) [-] with segment age? Only if rRLD not given beside (1=yes,2=no)'')')
       If (lUrf) Then
          Write (8,'(I1)') 1
       Else
          Write (8,'(I1)') 2
       Endif
       Write (8,*)
       Write (8,'(''* Root Length distribution of the root system'')')
       Write (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
       Write (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
       Write (8,'(''    Element#    RLD'')')
       Do i=1,nElm
          Write (8,'(I6,3X,F12.10)') i,RLD(i)
       Enddo
       Close (8)

    Elseif (ldJvL) Then
       Write (file,'(A8)')'out/dJvL'
       Write (file(9:9),'(I1)') ipl
       Write (file(10:10),'(A1)') '.'
       If (kout.Lt.10) Then
          Write (file(11:11),'(I1)') kout
       Elseif (kout.Lt.100) Then
          Write (file(11:12),'(I2)') kout
       Else
          Write (file(11:13),'(I3)') kout
       Endif
       Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
       Write (8,'(''********** DE JONG VAN LIER APPROACH COUPLED TO COUVREUR ET AL. (2012) RWU MODEL **********'')')
       Write (8,*)
       Write (8,'(''Number of past days considered in the effective sink (if set to zero, only the last sink will be considered)'')')
       Write (8,'(F5.2)') tlim_dJvL
       Write (8,*)
       Write (8,'(''* Distribution of the geometrical parameter for depletion zone characterizing'')')
       Write (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
       Write (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
       Write (8,'(''    Element#    Rho'')')
       Write (*,*) '... Calculating Rho ...'
       Allocate (Rho(1:nElm))
       Do i=1,nElm
          If (RLD(i).Gt.0) Then
             Rho(i)=4/((RSD(i)/RLD(i)/2/pi)**2-0.53**2/pi/RLD(i)+2*(1/(pi*RLD(i))+(RSD(i)/RLD(i)/2/pi)**2)*Log(0.53*2*pi*RLD(i)/RSD(i)/Sqrt(pi*RLD(i))))
          Else
             Rho(i)=999.9999
          Endif
          Write (8,'(I6,3X,F10.6)') i,Rho(i)
       Enddo
       Close (8)
    Endif
    Return
  End Subroutine RLDdis
 !********************************************************************************************
 !> ### Couvreur RWU model input; Couvreur.in ###
  Subroutine CouIn(ipl)!Couvreur feb 2011
    Use typedef
    Use GridData
    Use RootData, Only: Krs,Kcomp,ldJvL,Rho,tlim_dJvL,lno_Archi,lSUF,nrec,ntimeobs,timeobs,Kcomp_mat,Krs_mat
    Use DoussanMat, Only: hx_min,stresfun
    Use Input, Only: MFPC_table
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap):: el,i_dJvL,i_UptakeDomain,i,ipl
    Real(dp):: element

    Open (Unit=10,FILE='in/Couvreur.in',STATUS='OLD',ERR=10)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) stresfun,hx_min,ntimeobs
    Allocate (timeobs(1:ntimeobs))
    Allocate (Krs_mat(1:ntimeobs))
    Allocate (Kcomp_mat(1:ntimeobs))
    Read (10,*)
    Read (10,*)
    Read (10,*) i_UptakeDomain
    If (i_UptakeDomain.Eq.2) lSUF=.True.   ! 2 rootsystem domain
    Read (10,*)
    Read (10,*)
    Read (10,*) i_dJvL   
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    If (i_dJvL.Eq.1.0) ldJvL=.True.
    If (lno_Archi) Then
       Read (10,*)
       Read (10,*)
       Read (10,*) (timeobs(i),i=1,ntimeobs)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) (Krs_mat(i),i=1,ntimeobs)
       Krs=Krs_mat(1)
       Read (10,*)
       Read (10,*)
       Read (10,*) (Kcomp_mat(i),i=1,ntimeobs)
       Kcomp=Kcomp_mat(1)
       !print*,Kcomp,Krs,timeobs(1)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       If (lSUF) Then
          Read (10,*) nrec(ipl)
          Read (10,*)
          Allocate (SSF(1:nrec(ipl)))
          Do el=1,nrec(ipl)
             Read (10,*,ERR=20) element,SSF(el)
          Enddo
       Else
          Read (10,*) nexSSF,neySSF,nezSSF,dxSSF,dySSF,dzSSF
          Read (10,*)
          Allocate (SSF_mat(1:nexSSF*neySSF*nezSSF,1:ntimeobs))
          Do el=1,nexSSF*neySSF*nezSSF
             Read (10,*,ERR=20) element,(SSF_mat(el,i),i=1,ntimeobs)
          Enddo
          Allocate (SSF(1:nexSSF*neySSF*nezSSF))
          SSF(1:nexSSF*neySSF*nezSSF)=SSF_mat(1:nexSSF*neySSF*nezSSF,1)
       Endif
    Endif
    Close (10)
    If (ldJvL) Then
       Call MFPC_Table
       Open (Unit=15,FILE='in/dJvL.in',STATUS='OLD',ERR=30)
       Read (15,*)
       Read (15,*)
       Read (15,*)
       Read (15,*,ERR=40) tlim_dJvL
       If (lno_Archi) Then
          Read (15,*)
          Read (15,*)
          Read (15,*)
          Read (15,*,ERR=40) nexRho,neyRho,nezRho,dxRho,dyRho,dzRho
          Read (15,*)
          Allocate (Rho(1:nElm))
          Do el=1,nexRho*neyRho*nezRho
             Read (15,*,ERR=40) element,Rho(el)
          Enddo
       Endif
       Close (15)
    Endif
    Return
10  Call stop_program('File  < Couvreur.in >  not found -- program terminated.')
20  Call stop_program('Data inconsistency in  < Couvreur.in >  -- program terminated.')
30  Call stop_program('File  < dJvL.in >  not found -- program terminated.')
40  Call stop_program('Data inconsistency in  < dJvL.in >  -- program terminated.')
  End Subroutine CouIn
  !********************************************************************************************
  !> ### Input for Feddes RWU; Feddes.in ###
  Subroutine FedIn
    Use typedef
    Use GridData
    Use RootData, Only: lUrf,lJarvis,OmegaC,lno_Archi
    Use DoussanMat, Only: stresfun
    Use PlntData, Only: h0,h1,h2,h3,p50,h50,p1,p2
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap):: el,element,i_Jarvis,i_Urf

    Open (Unit=10,FILE='in/Feddes.in',STATUS='OLD',ERR=10)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) stresfun
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) h0,h1,h2,h3
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) p50,h50,p1,p2
    Read (10,*)
    Read (10,*)
    Read (10,*) i_Jarvis
    If (i_Jarvis.Eq.1) lJarvis=.True.
    Read (10,*)
    Read (10,*)
    Read (10,*) OmegaC
    Read (10,*)
    Read (10,*)
    Read (10,*) i_Urf
    If (i_Urf.Eq.1) lUrf=.True.
    If (lno_Archi) Then
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) nexRLD,neyRLD,nezRLD,dxRLD,dyRLD,dzRLD
       Read (10,*)
       Allocate (RLD(1:nexRLD*neyRLD*nezRLD))
       Do el=1,nexRLD*neyRLD*nezRLD
          Read (10,*,ERR=20) element,RLD(el)
       Enddo
    Endif
    Close (10)
    Return
10  Call stop_program('File  < Feddes.in >  not found -- program terminated.')
20  Call stop_program('Data inconsistency in  < Feddes.in >  -- program terminated.')
  End Subroutine FedIn
  !**************************************************************************************
  !> ### Calculation of the sink term ###
  Subroutine SetSnk(t,BCr,BCtp)
    Use TypeDef
    Use EnviData
    Use TardieuData
    Use Paramdata,Only: lChem
    Use TempData, Only: tem
    Use PlntData, Only: Tpot,Tact,PHcollar,CMm,VMax,fk,xin,TotSur
    Use GridData, Only: nPt,nElm,betac,sink,csink,zgrid,elmnod,subN,iL,SSF,VElm,HElm,dxGrid,dyGrid,dzGrid,RLD
    Use DoussanMat, Only : sink_cube,betac_cube,csink_cube,stressBC,Lr,Lr_pot,Inv_c1KxKr,Joutr,w_sub,nsub,cube_i,nplant,PHs,loc_Q,w_dis
    Use RootData, Only: lCou,lSUF,lDou,lFed,lJarvis,OmegaC,Krs,Kcomp,lno_RWU,ldJvL,tlim_dJvL,Rho,lPast,tPast,nrec,Hseq,lSign_new,nplant,nrec_m
    Use SolData, Only: hnew,conc,par,theta,MatNum
    Use Watfun, Only: Fh_from_Th,Fh_from_mfp_soiltab,Fmfp_soiltab
    Use Doussan, Only: find_stomatal_aperture, capacitance, Gap, AQPc
    Implicit None

    Integer(ap),Intent(inout) :: BCtp
    Integer(ap) ::i,j,k,l,iE,iSE,isub,irecn,ipl,c_i
    Real(dp), Allocatable, Dimension(:) :: Phi,wPast,Mseq
    !Real(dp) :: alpha
    Real(dp) :: active,active_cube,ConcE,t,HElm2(nElm),J_Loc,Lr_fact(1:nrec_m) !,Lr_fact(1:nrec(nplant))
    Real(dp), Intent(inout) ::BCr
    Real(dp)::gs_loc

    If (lFed) Then
       Tpot(1)=Bcr 
       Do i=1,nElm
          If (RLD(i).Gt.1.E-20_dp) Then
             sink_cube(i)=RLD(i)/Sum(RLD)*Abs(Tpot(1))/VElm(i)*alpha(Fh_from_Th(Sum(theta(elmnod(1:8,i)))/8.0_dp,par(:,MatNum(elmnod(1,i)))),Sum(Conc(elmnod(1:8,i)))/8.0_dp,Sum(tem(elmnod(1:8,i)))/8.0_dp)
          Else
             sink_cube(i)=0.0_dp
          Endif
       End Do
       Tact(1)=DOT_Product(sink_cube(1:nElm),VElm)
       If (lJarvis.And.OmegaC.Lt.1) Then
          sink_cube=sink_cube/Max(OmegaC,Tact(1)/Abs(Tpot(1)))
          Tact(1)=DOT_Product(sink_cube(1:nElm),VElm)
       Endif

    Elseif (lCou.And..Not.lSUF) Then
       Allocate (Phi(1:nElm))
       If (ldJvL) Then
          Allocate (wPast(1:lPast))
          Allocate (Mseq(1:nElm))
          If (lPast.Gt.2) Then
             wPast(1)=((tPast(2)-t)/tlim_dJvL)**2+2*((tPast(2)-t)/tlim_dJvL)+1
             wPast(2:lPast-1)=((tPast(3:lPast)-t)/tlim_dJvL)**2-((tPast(2:lPast-1)-t)/tlim_dJvL)**2+2*((tPast(3:lPast)-t)/tlim_dJvL)-2*((tPast(2:lPast-1)-t)/tlim_dJvL)
             wPast(lPast)=-((tPast(lPast)-t)/tlim_dJvL)**2-2*((tPast(lPast)-t)/tlim_dJvL)
          Elseif (lPast.Eq.2) Then
             wPast(1)=0.0_dp
             wPast(2)=1.0_dp
          Else
             wPast(1)=1.0_dp
          Endif
          Do i=1,nElm
             !> Calculation of the layer's Hseq, here called HElm
             HElm(i)=Fh_from_mfp_soiltab(&
                  Fmfp_soiltab(&
                  Fh_from_Th(&
                  Sum(theta(elmnod(1:8,i)))/8.0_dp,par(:,MatNum(elmnod(1,i)))),&
                  MatNum(elmnod(1,i)),&
                  par(:,MatNum(elmnod(1,i))))-sink_cube(i)/Rho(i),&
                  MatNum(elmnod(1,i)))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_dp
          End Do
       Else
          Do i=1,nElm
             HElm(i)=Fh_from_Th(Sum(theta(elmnod(1:8,i)))/8.0_dp,par(:,MatNum(elmnod(1,i))))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_dp!> water potential of the mean node water content
          End Do
       Endif
       Hseq=DOT_Product(HElm(1:nElm),SSF)
       !print*,Hseq

! use Tardieu and Davies model for stomtal apertuire gs and then run Penmann-Monteith
       If (lSign_new) Then
          Call find_stomatal_aperture(Hseq,Krs,gs_loc)   
          gs_t=gs_loc
          J_loc = 1000000*PM_fun(gs_loc,Sshaded,s_t,Rn_t,rho_t,Cp,ga,VPD_t*1000,lambda_t,gamma_t,conversion_t)! original equation is in mm³/sec/plant
          Tact(1) = J_loc*86.4 !cm/day
          PHcollar=PHcollar_fun(Hseq,Krs,Tact(1))
          Call capacitance
          stressBC=.False.
       Else
          If (BCtp.Eq.2) Then
             PHcollar=-Abs(BCr)/Krs+Hseq
             Tpot(1)=Abs(BCr)![L³/T]
             stressBC=.False.
             Call RootStressCouvreur(PHcollar)  ! Check if root collar < min root collar
          Else
             PHCollar=BCr
          Endif
          If (BCtp.Eq.1.Or.stressBC) Then
             Tact(1)=Krs*(-PHcollar+Hseq)!> Convert the collar pressure head into the actual transpiration
          Else
             Tact(1)=Tpot(1)![L³/T]
          Endif
       Endif
       !print*,Tact(1)
       Phi=Kcomp*(HElm-Hseq)
       sink_cube=Abs(Tact(1))*SSF/dxGrid/dyGrid/dzGrid+Phi*SSF/dxGrid/dyGrid/dzGrid

       If (ldJvL) Then
          !> Calculation of the layer's Hseq, here called HElm
          Do i=1,nElm
             HElm2(i)=Fh_from_mfp_soiltab(&
                  Fmfp_soiltab(&
                  Fh_from_Th(&
                  Sum(theta(elmnod(1:8,i)))/8.0_dp,par(:,MatNum(elmnod(1,i)))),&
                  MatNum(elmnod(1,i)),par(:,MatNum(elmnod(1,i))))-sink_cube(i)/Rho(i),&
                  MatNum(elmnod(1,i)))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_dp
          End Do
          Phi=Kcomp*(HElm-Hseq)
          sink_cube=0.5_dp*sink_cube+0.5_dp*(Abs(Tact(1))*SSF/dxGrid/dyGrid/dzGrid+Phi*SSF/dxGrid/dyGrid/dzGrid)
          HElm=(HElm+HElm2)/2.0_dp
       Endif
    Elseif (lSUF) Then
       Do ipl=1,nplant ! Code not ready for multiple plants yet
          PHs(1:nrec(ipl),ipl)=0._dp
          Do irecn=1,nrec(ipl)
             Do isub=1,nsub(irecn,ipl)
                PHs(irecn,ipl)=PHs(irecn,ipl)+Sum(hnew(loc_Q(irecn,1:8,isub,ipl))*w_dis(irecn,1:8,isub,ipl)/Sum(w_dis(irecn,1:8,isub,ipl)))*w_sub(irecn,isub,ipl)/Sum(w_sub(irecn,1:nsub(irecn,ipl),ipl))
             End Do    !end do-loop over subsegments
             Lr_fact(irecn)=Gap(PHs(irecn,ipl))*AQPc(PHs(irecn,ipl),ipl)
             Lr(irecn,ipl)=Lr_pot(irecn)*Lr_fact(irecn)
          End Do
          Krs=DOT_Product(Inv_c1KxKr(1:nrec(ipl)),Lr_fact(1:nrec(ipl)))
          SSF(1:nrec(ipl))=Inv_c1KxKr(1:nrec(ipl))*Lr_fact(1:nrec(ipl))/Krs
          Kcomp=Krs
          Print *,'SUM(SSF) Krs',Sum(SSF(1:nrec(ipl))),Krs
          SSF(1:nrec(ipl))=SSF(1:nrec(ipl))/Sum(SSF(1:nrec(ipl)))
          Hseq=DOT_Product(PHs(1:nrec(ipl),ipl),SSF)
          If (BCtp.Eq.2) Then
             PHcollar=-Abs(BCr)/Krs+Hseq
             Tpot(ipl)=Abs(BCr)![L³/T]
             stressBC=.False.
             Call RootStressCouvreur(PHcollar)   ! check if collar xylem < collar min
          Else
             PHCollar=BCr
          Endif
          If (BCtp.Eq.1.Or.stressBC) Then
             Tact(1)=Krs*(PHcollar-Hseq) ![L³/T] Convert the collar pressure head into the actual transpiration 
          Else
             Tact(1)=Tpot(ipl) ![L³/T]
          Endif
          Joutr(1:nrec(ipl),ipl)=(Abs(Tact(1))+Kcomp*(PHs(1:nrec(ipl),ipl)-Hseq))*SSF(1:nrec(ipl))
          sink_cube(1:nElm)=0.0
          Do irecn=1,nrec(ipl)
             Do isub=1,nsub(irecn,ipl)
                c_i=cube_i(irecn,isub,ipl)
                sink_cube(c_i)=sink_cube(c_i)+Joutr(irecn,ipl)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
             End Do
          End Do
       End Do ! plant loop
    Elseif (lno_RWU) Then
       Do i=1,nPt
          sink(i)=0.0_dp
       End Do
    Endif
    If (lChem) Then
       Do i=1,nPt
          If (betac(i).Gt.1.E-20_dp) Then
             active=TotSur*betac(i)*(VMax/(CMm+Conc(i))+fk)
          Else
             active=0.0_dp
          End If
          csink(i) = xin*sink(i)+(1-xin)*active
       End Do
    Endif
    Do iE=1,nElm
       If (lDou.And.lChem) Then
          If (betac_cube(iE).Gt.1.E-20_dp) Then
             Do iSE=1,5
                i=elmnod(iL(1,iSE,subN(iE)),iE)
                j=elmnod(iL(2,iSE,subN(iE)),iE)
                k=elmnod(iL(3,iSE,subN(iE)),iE)
                l=elmnod(iL(4,iSE,subN(iE)),iE)
                ConcE = (conc(i) + conc(j) + conc(k) + conc(l) )/4._dp
                active_cube=TotSur*betac_cube(iE)*(VMax/(CMm+ConcE)+fk)
             Enddo
          Else
             active_cube=0.0_dp
          End If
          csink_cube(iE) = xin*sink_cube(iE)+(1-xin)*active_cube
       Endif
    End Do
    Return
  End Subroutine SetSnk
 !**************************************************************************************
 !> ### Calculate osmotic potential [L] from soil solution concentration [M/L3]###
  Function alpha(h,Conc,tem)
    Use typedef
    Use PlntData
    Use DoussanMat, Only: stresfun
    Implicit None

    Real(dp), Intent(in) ::conc,tem,h
    Real(dp):: alpha,alphap,alphah,p,rktemp

    !>     To use the first or the second method simply render active the lines
    !>     where necessary quantities are calculated for the chosen expression and
    !>    comment out the lines related to the other expression; recompile the
    !>     program; beware of units and, if needed, make proper aggiustments.
    !>     Also, adjust for the valency of the chemical species studied if the
    !>     second method is adopted.
    !> **************************
    !>     For very dilute solutions the van´t Hoff expression can be used (T.L.
    !>     Allen and R.M. Keefer, "Chemistry: Experiment and Theory", 1974,
    !>     Harper and Row Publishers, pp 283-285):
    !>                              p=CRT
    !>     where:
    !>     p is the osmotic potential [L], C is the concentration [M/L3], R is a
    !>     constant [(L3*pressure)/(M*temperature)], and T is the temperature
    !>     (K degrees)
    !>
    !>     first transform temperature degrees from Celsius to Kelvin...
    !
    rKtemp=tem+273.15_dp
    !
    !>     1 liter = 1000 cm3; 1 atm = 1033.6 cm H2O;
    !>     R = 0.082 (liter atm K-1 mole-1) = 84811.0144 (cm3 cm K-1 mole-1);
    !
    p=(-Conc)*rKtemp*84811.0144_dp
    !
    !> **********************************
    !>     The following constants are taken from ´Diagnosis and Improvement of
    !>     Saline and Alkali Soils´, USDA Agriculture Handbook no. 60, 1954:
    !>               1 atm = 0.36 millimhos/cm           (Fig 6, p 15);
    !>               1 millimhos/cm = -0.091 meq/liter   (Fig 4, p 12);
    !>     these values are valid for the range of EC that allows plant growth.
    !>    Therefore 1 atm = -0.3276 meq/liter; also:
    !>               1 atm = 1033.6 cm H2O;     and:
    !>               1 meq/liter = 1.0E+06/n M/cm3, being n the valency of the      ion
    !>                                                        in the soil solution.
    !>
    !>      p=Conc*(-0.3276)*1000000/1033.6
    !>
    !> ********************************************************************************
    alpha=0.0_dp
    If (stresfun.Eq.2) Then
       alphah=1._dp/(1._dp+(h/h50)**p1)
       If (p50.Lt.9.E+20_dp) Then
          alphap=1._dp/(1._dp+(p+p50)**p2)
          alpha=alphah*alphap
       Else
          alpha=alphah
       Endif
    Elseif (stresfun.Eq.1) Then
       If ((h.Gt.h3).And.(h.Lt.h2)) Then
          alpha=(h-h3)/(h2-h3)
       Elseif ((h.Ge.h2).And.(h.Le.h1)) Then
          alpha=1.0_dp
       Elseif ((h.Gt.h1).And.(h.Lt.h0)) Then
          alpha=(h-h0)/(h1-h0)
       Else
          alpha=1.E-9_dp
       Endif
    Else
       alpha=1.0_dp
    Endif
    Return
  End Function alpha
 !**************************************************************************
!> ### check if the root collar abs(PH) is larger than abs(hx_min)+tolerance 
!  and adapt the collar BC###
  Subroutine RootstressCouvreur(PHtop)
    Use Typedef
    Use DoussanMat, Only: stressBC,hx_min,stresfun
    Use GridData, Only: RelEps,epslonR,factorRelEps  
    Implicit None

    Real(dp), Intent(inout):: PHtop
    Real(dp):: del_PHr


    If (RelEps) Then 
       del_PHr=-Max(Abs(hx_min/factorRelEps),epslonR)
    Else
       del_PHr=-epslonR
    Endif
    If ((stresfun.Eq.1).And.(PHtop<hx_min+del_PHr)) Then
       !> top node at lower PH than allowed: start of stressed conditions
 !      Write (*,'(/,a,1pe11.3,a,1pe11.3,a)') ' stress in the collar xylem: PHtop=',PHtop,' is lower than criterion ',hx_min
       write(*,'(a)',advance='no') '||' !when stress occurs a vertical bar appears on the screen
       stressBC=.True.
       PHtop=hx_min
    Endif
  End Subroutine RootstressCouvreur
  !********************************************************************************

End Module Sink
