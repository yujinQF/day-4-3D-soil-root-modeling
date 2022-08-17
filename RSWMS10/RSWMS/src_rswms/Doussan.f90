!> \file Doussan.f90


!> Module Doussan
Module Doussan

Contains
  !> calculates where roots intersect which plane of a soil voxel
  Logical Function intsec(xA,yA,zA,xB,yB,zB,ifoln,irecn,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)
    Use Typedef
    Use GridData, Only: dxgrid,dygrid,nex,ney,continu
    Use DoussanMat, Only: transroot,nsub
    Implicit None

    Integer(ap) ::iFace,irecn,ifoln,isub,ipl
    Real(dp) ::xA,yA,zA,xB,xBp,yB,yBp,zB,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,xc
    Real(dp) ::yc,zc,f

    intsec=.False.
    If (continu) Then
       !> if continous domain, reconstruction of the continuous B-Ap segment
       If (isub.Eq.1) Then
          xBp=xB+(transroot(ifoln,1,nsub(ifoln,ipl),ipl)-transroot(irecn,1,1,ipl))*(nex*dxgrid)!ifoln refers to the "Ap" node and irecn to the "B" node
          yBp=yB+(transroot(ifoln,2,nsub(ifoln,ipl),ipl)-transroot(irecn,2,1,ipl))*(ney*dygrid)
       Else
          xBp=xB+(transroot(irecn,1,isub,ipl)-transroot(irecn,1,1,ipl))*(nex*dxgrid)
          yBp=yB+(transroot(irecn,2,isub,ipl)-transroot(irecn,2,1,ipl))*(ney*dygrid)
       Endif
    Else
       xBp=xB
       yBp=yB
    Endif

    If (zB.Lt.z1) Then
       !> may have intersection with x-y-plane of cube at z1
       f=(z1-zA)/(zB-zA)
       xc=xA+f*(xBp-xA)
       yc=yA+f*(yBp-yA)
       If (((xc.Ge.x1).And.(xc.Le.x2).And.yc.Ge.y1).And.(yc.Le.y2)) Then
          xInt=xc
          yInt=yc
          zInt=z1
          iFace=1
          Goto 1
       Endif
    Endif
    If (zB.Gt.z2) Then
       !> may have intersection with x-y-plane of cube at z2
       f=(z2-zA)/(zB-zA)
       xc=xA+f*(xBp-xA)
       yc=yA+f*(yBp-yA)
       If ((xc.Ge.x1).And.(xc.Le.x2).And.(yc.Ge.y1).And.(yc.Le.y2)) Then
          xInt=xc
          yInt=yc
          zInt=z2
          iFace=2
          Goto 1
       Endif
    Endif
    If (yBp.Lt.y1) Then
       !> may have intersection with x-z-plane of cube at y1
       f=(y1-yA)/(yBp-yA)
       xc=xA+f*(xBp-xA)
       zc=zA+f*(zB-zA)
       If ((xc.Ge.x1).And.(xc.Le.x2).And.(zc.Ge.z1).And.(zc.Le.z2)) Then
          xInt=xc
          yInt=y1
          zInt=zc
          iFace=3
          Goto 1
       Endif
    Endif
    If (yBp.Gt.y2) Then
       !> may have intersection with x-z-plane of cube at y2
       f=(y2-yA)/(yBp-yA)
       xc=xA+f*(xBp-xA)
       zc=zA+f*(zB-zA)
       If ((xc.Ge.x1).And.(xc.Le.x2).And.(zc.Ge.z1).And.(zc.Le.z2)) Then
          xInt=xc
          yInt=y2
          zInt=zc
          iFace=4
          Goto 1
       Endif
    Endif
    If (xBp.Lt.x1) Then
       !> may have intersection with y-z-plane of cube at x1
       f=(x1-xA)/(xBp-xA)
       yc=yA+f*(yBp-yA)
       zc=zA+f*(zB-zA)
       If ((yc.Ge.y1).And.(yc.Le.y2).And.(zc.Ge.z1).And.(zc.Le.z2)) Then
          xInt=x1
          yInt=yc
          zInt=zc
          iFace=5
          Goto 1
       Endif
    Endif
    If (xBp.Gt.x2) Then
       !> may have intersection with y-z-plane of cube at x2
       f=(x2-xA)/(xBp-xA)
       yc=yA+f*(yBp-yA)
       zc=zA+f*(zB-zA)
       If ((yc.Ge.y1).And.(yc.Le.y2).And.(zc.Ge.z1).And.(zc.Le.z2)) Then
          xInt=x2
          yInt=yc
          zInt=zc
          iFace=6
          Goto 1
       Endif
    Endif
    Return
1   intsec=.True.
    Return
  End Function intsec
  !===============================================================================
  !> Doussan model implementation through 3 Subroutines
  !>- SetupDou
  !>       - segment
  !>- SetBCroot
  !>       - conductroot
  !>- SolveRoot
  !> and the adaptation of SetSink
 

  !*******************************************************************************
! ### fill up the Doussan matrix based on the root connectivity and root hydraulic
! properties ###  
  Subroutine SetupDou(t,dt,ipl)
    Use TypeDef
    Use Paramdata,Only: ldirect,lChem,lretry,last_out,maxemg 
    Use SparseMatrix
    Use DoussanMat
    Use GridData, Only:xgrid,ygrid,zgrid,dzgrid,dygrid,dxgrid,nex,ney,nElm,continu,nPt,betac,elmnod,subN,iL
    Use RootData, Only: lCalloc,loliumroottyp,maizeroottyp,wheatroottyp,nbr,segsur,ibrseg,seglen,xplant,yplant,ibrgrw,zs,ys,xs,&
         zg,yg,xg,ordseg,timorg,irecpr,nurf,age,urf,lDou,lFed,lCou,lno_Archi,irecsg,nplant
    Use tmctrl, Only: t_begin,tout
    USE SolData, ONLY :nmat
    Use PlntData, Only: TotSur
    Use MPIutils, Only: stop_program
    Use PlantGrowth, Only: Stress, SetTp
    Implicit None

    Integer(ap),INTENT(in) ::ipl
    Integer(ap) ::ibr,irecn,ifoln,igrow,iprvn,iBCn,typ,count,kk,m_tmp,most
    Integer(ap) ::isub,i,j,corner(1:8),iUrf,iseg,iorder,k,l,iE,iSE,mat_seg_tmp(1:8)
    Integer(ap) :: err,ibranch,max_subseg
    Real(dp) ::t,xA,yA,zA,xB,yB,zB,x1,x2,y1,y2,z1,z2,dt,Weight,betce,VEl,Sbetac, segage
    Real(dp) ::PrvSur,w_sub_tmp(1:isubmax)!10 is the max number of subsegemnts
    Real(dp) :: y(1:nLibr),h_mat(1:nLibr),arr(1:nLibr-1),cumsum(1:nLibr-1),PHtemp
    Real(dp) :: theta(1:nLibr),Kh(1:nLibr),C(1:nLibr),temp_phr
    Real(dp) :: rs,concrs
    Logical :: n_apex,run

 
    If (lDou.Or.(lCou.And.(.Not.lno_Archi))) Then  
       !Allocate (nBC_irecn(1:maxemg))
       nBC_irecn=0

       !> initialize the pointer of the first el of the list
       !> initialize matrices
       IF (ipl.EQ.1) Call IniMat
       iBCn=0

       !> shift rootsystem for continous domain - needed if macropores for ConductRoot
       If(continu) Then
          Do irecn=1,nrec(ipl)
             Call roottrans(xs(irecn,ipl),ys(irecn,ipl),irecn,1,1)
          End Do
          Do igrow=1,ngrow(ipl)
             Call tiptrans(xg(igrow,ipl),yg(igrow,ipl),igrow,1,1) 
          End Do
          transroot(nrec(ipl)+1:nrec(ipl)+ngrow(ipl),:,1,1)=transtip(1:ngrow(ipl),:,1,1)
       End If

       !> Current axial and radial conductance matrix
       Call ConductRoot(t,ipl)

       !> prepare calculation of betaw and betac
       PrvSur=0.0_dp
       Do iseg=1,nrec(ipl)+1
          PrvSur=PrvSur+segsur(iseg,ipl)!total segment surface
       Enddo
       If (PrvSur.Lt.1.E-20_dp) PrvSur=1.E-20_dp
       TotSur = PrvSur
       betac=0.0_dp
       !> go through different plants
    Elseif (lFed.And.(.Not.lno_Archi)) Then!For Doussan, that is done in IniMat
       Allocate (transroot(0:nrec(ipl)+ngrow(ipl),1:2,1:isubmax,ipl))
       transroot=0
       Allocate (nsub(0:nrec(ipl)+ngrow(ipl),ipl))
       nsub=0
       Allocate (cube_i(0:nrec(ipl),1:isubmax,ipl))
       cube_i=0
       Allocate (Intc(0:nrec(ipl)+ngrow(ipl),1:3,1:isubmax,ipl))
       Intc=0._dp
       Allocate (l_seg(0:nrec(ipl)))!actually not used
       l_seg=0._dp
       Allocate (beta_weight(0:nrec(ipl),1:isubmax,ipl))!actually not used
       beta_weight=0._dp
       Allocate (w_sub(0:nrec(ipl),1:isubmax,ipl))!actually not used
       w_sub=0._dp
       Allocate (l_sub(0:nrec(ipl),1:isubmax,ipl))!actually not used
       l_sub=0._dp
       Allocate (cent(0:nrec(ipl),1:3,1:isubmax,ipl))!actually not used
       cent=0._dp
       Allocate (sum_dis(0:nrec(ipl),1:isubmax,ipl))
       sum_dis=0._dp
       Allocate (w_dis(0:nrec(ipl),1:8,1:isubmax,ipl))
       w_dis=0._dp
       betaw=0._dp
    Endif

    If (lretry) Then
       t=tOut(last_out)
       Write (*,'(///'' Simulation continues at time = '',F8.3,''.'')') t
    Endif

       !> Current boundary conditions for root
       If (lDou) Then
          If (.Not.lCalloc) Then
             Call SetBCroot(t,curr_BCr(ipl),curr_BCtp(ipl))
          Else
             If (ipl.Ge.2) Call stop_program('Assimilate allocation currently doesnt work with multiple plants')
             Call Stress(rs,concrs)
             Call SetTp(t,rs,concrs,ipl)
             curr_BCr(ipl)=BCr_usr(ipl)     ! -ABS(Tpot) 
             curr_BCtp(ipl)=2
          Endif
          !multiple roots
          err=SM_allocate(plantmatrix(ipl), nrec(ipl)+1, nrec(ipl)+1)
          If(err/=0) Call stop_program('Could not create plantmatrix')
       Endif
       !> go through each root segment and update the node surface function:
       Do ibr=1,nbr(ipl)
          n_apex=.False.
          !> find the tip segment of the branch 'ibr'
          irecn = irecsg(ibr,ipl)

          !> the first one we find is an apex
          If (seglen(irecn,ipl)<1.E-20) Then !> skip this segment too small to be taken
             xA=xs(irecn,ipl)+xplant(ipl)
             yA=ys(irecn,ipl)+yplant(ipl)
             zA=zs(irecn,ipl)
             If (continu) Then
                xA=xA+transroot(irecn,1,1,ipl)*dxgrid*nex
                yA=yA+transroot(irecn,2,1,ipl)*dygrid*ney
             End If

             ifoln=irecn !> following node ID
             irecn=irecpr(irecn,ipl) !> current node ID

          Else !> segment long enough

             n_apex=.True.!> ifoln does not exist if not continu
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==ibr) Then !> apex of branch ibr
                   xA=xg(igrow,ipl)+xplant(ipl)
                   yA=yg(igrow,ipl)+yplant(ipl)
                   zA=zg(igrow,ipl)
                   If (continu) Then
                      xA=xA+transroot(nrec(ipl)+igrow,1,1,ipl)*dxgrid*nex
                      yA=yA+transroot(nrec(ipl)+igrow,2,1,ipl)*dygrid*ney
                      ifoln=nrec(ipl)+igrow                             
                      nsub(nrec(ipl)+igrow,ipl)=1
                   End If
                Endif
             End Do
          Endif

          If (irecn==0) Then !> there exists a branch ibr but not yet any segment!
             run=.False.
          Else
             run=.True.
          Endif

          !> then the rest of the branch up to the seed of the embranchment
          Do While (run)
             !> "upper" node
             iprvn=irecpr(irecn,ipl)
             !> location of the rear or "upper" end
             xB=xs(irecn,ipl)+xplant(ipl)
             yB=ys(irecn,ipl)+yplant(ipl)
             zB=zs(irecn,ipl)
             If (continu) Then
                xB=xB+transroot(irecn,1,1,ipl)*dxgrid*nex
                yB=yB+transroot(irecn,2,1,ipl)*dygrid*ney
             End If
             If (lDou.Or.(lCou.And.(.Not.lno_Archi))) Then
                !> calculate the gravity components z (always positive & maximum at soil surface)
                GH(irecn,ipl)=(zA+zB)/2.
             Endif
             !> calculate segment weighing factor according to age:
             iorder=ordseg(irecn,ipl)
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
            ! Weight for solute
             segage=t-timorg(irecn,ipl)
             If (lChem.And.(segage.Ge.0.0_dp)) Then
                If (segage.Ge.age(typ,nUrf(typ))) Then
                   Weight=Urf(typ,nUrf(typ))
                Else
                   iUrf=nUrf(typ)
                   iUrf=iUrf-1
                   Do While ((segage.Lt.age(typ,iUrf)).And.(iUrf.Gt.1))
                      iUrf=iUrf-1
                   Enddo
                   Weight=Urf(typ,iUrf)+(segage-age(typ,iUrf))/&
                        (age(typ,iUrf+1)-age(typ,iUrf))*&
                        (Urf(typ,iUrf+1)-Urf(typ,iUrf))
                Endif
             Else
                Weight=1.0_dp
             Endif

             !> calculate number of subsegment for node/segment irecn and corresponding weights
             Call segment(xA,yA,zA,xB,yB,zB,ifoln,irecn,ipl,Weight,PrvSur)
			 
             !> associate each segment with a main material and a main voxel for the rhizosphere
             w_sub_tmp=0
             do k=1,nsub(irecn,ipl)
                  w_sub_tmp(k)=w_sub(irecn,k,ipl)
             enddo
             max_subseg=maxloc(w_sub_tmp(1:nsub(irecn,ipl)),1)!biggest subseg of the segment
             vox_seg(irecn,ipl)=cube_i(irecn,max_subseg,ipl)!voxel around the biggest subsegment			 
             if (nMat.eq.1) Then
                  mat_seg(irecn,ipl)=1
             else
                 if (nsub(irecn,ipl).EQ.1) Then !only one voxel or only one material
                   mat_seg_tmp(1:8)=mat_Q(irecn,1:8,1,ipl)!define material around rootsmat_seg(irecn,1:8,ipl)
                Else
                   mat_seg_tmp(1:8)=mat_Q(irecn,1:8,max_subseg,ipl)
                endif
                !> take the most common material
                most=0
                DO k=1,8
                    m_tmp=mat_seg_tmp(k)
                    count=1
                    do kk=k+1,8
                       if(mat_seg_tmp(kk).eq.m_tmp) count=count+1
                    enddo
                    if (most<count) then
                       most=count
                       mat_seg(irecn,ipl)=m_tmp
                    endif
                ENDDO
             Endif

             If (lDou) Then
                !> several changes -> multiple roots
                err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                !> if apex (bottom part of root)
                If (n_apex) Then
                   Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn,ipl)/seglen(irecn,ipl)+Lr(irecn,ipl)*segsur(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Else
                   Err=SM_set(Plantmatrix(ipl), irecn+1,ifoln+1,-Khr(ifoln,ipl)/seglen(ifoln,ipl),.False.) !row, col,value
                   Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn,ipl)/seglen(irecn,ipl)+Khr(ifoln,ipl)/seglen(ifoln,ipl)+Lr(irecn,ipl)*segsur(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Endif

                !> define 1st part of Q (Q=Qi.*PHsoil+Qbc) -> RHS
                Qi(irecn,ipl)=Lr(irecn,ipl)*segsur(irecn,ipl)
                !if(irecn.lt.100) print*,'qi',irecn,qi(irecn,ipl)
             Endif
             !> if reached the seed or the embranchement => change branch
             If (iprvn==0) Then!seed=first segment
                If (lDou) Then
                   GH(0,ipl)=zB
                   iBCn=iBCn+1
                   nBC_irecn(iBCn)=irecn
                   If (curr_BCtp(ipl)==2.And..Not.(ldirect)) Then!flux
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                      !> own additive; position (2,2) changes
                      !> own additive;position (2,1) = 0, changes into zero
                      !> rhs has only 1 entry (first one)
                      Qi(iprvn,ipl)=0._dp
                      Q_bc1(iprvn,ipl)=curr_BCr(ipl) !Q(0,0)
                      !                     Q_bc(irecn,ipl)=curr_BCr
                   Else If ((curr_BCtp(ipl)==2).And.(ldirect)) Then  !flux + DIRECT
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.True.)
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,0._dp,.True.)
                      Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.True.)
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                      !own additive; position (2,2) changes 
                      !own additive;position (2,1) = 0, changes into zero
                      !rhs has only 1 entry (first one)
                      Qi(iprvn,ipl)=0._dp
                      !Q_bc1(iprvn,ipl)=curr_BCr(ipl)!Q(0,0)
                      !                     Q_bc(irecn,ipl)=curr_BCr
                      Do ibranch=1,1  !> all the branches connected to the seed
                         irecn=nBC_irecn(ibranch)
                         temp_phr=-(curr_BCr(ipl)/Khr(irecn,ipl)*seglen(irecn+1,ipl)+200.)

                         Qi(iprvn,ipl)=0._dp !> from seed has to be given, from irecn is already stated in DoussanMat
                         Q_bc1(iprvn,ipl)=temp_Phr*.1+GH(iprvn,ipl) !> use ph from last time step
                         Q_bc2(iprvn,ipl)=temp_Phr*2.+GH(iprvn,ipl)
                         Q_bc1(irecn,ipl)=(temp_Phr*.1+GH(iprvn,ipl))*Khr(irecn,ipl)/seglen(irecn,ipl) !> all the PH for root was given in total head!!
                         Q_bc2(irecn,ipl)=(temp_Phr*2+GH(iprvn,ipl))*Khr(irecn,ipl)/seglen(irecn,ipl)
                      Enddo
                   Else If (curr_BCtp(ipl)==1) Then !PH
                      !iprvn+1,iprvn+1
                      Qi(iprvn,ipl)=0._dp
                      !> true means overwrite (delete first), no addition
                      !iprvn+1,iprvn+1
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.True.) !position (1,1)
                      Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.True.) 
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                      !> first entry in rhs is PH (inlc. gravity)
                      Q_bc1(iprvn,ipl)=curr_BCr(ipl)+GH(iprvn,ipl)!Q(iprv)-BCr*Kh(irecn)/seglen(irecn);
                      !> second entry is PH incl. gravity times these parameters
                      Q_bc1(irecn,ipl)=(curr_BCr(ipl)+GH(iprvn,ipl))*Khr(irecn,ipl)/seglen(irecn,ipl)
                      Q_bc2 = 0._dp
                   Endif
                Endif

                run=.False.

             Elseif (ibrseg(iprvn,ipl).Ne.ibrseg(irecn,ipl)) Then !> start of the branch but not from the seed -> gaat van onder na boven, duz iprvn (nieuwe positie) is nu niet van die ene branch maar van een side branch
                If (lDou) Then
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Endif
                run=.False.
             Endif

             !> definition for the next run of the loop
             ifoln=irecn
             irecn=iprvn
             !> location of the final node
             xA=xB!> A directly redefined with report to B in order not to have to calculate all the translations of A to the inside of the domain
             yA=yB 
             zA=zB
             !> from here, not an apex
             n_apex=.False.
          End Do !> loop on branch nodes
       End Do !> loop on root branches

       If (lDou.Or.(lCou.And.(.Not.lno_Archi))) Then !there is RWU
          If (.Not.(old)) Then
             If ((ave) .Or. (eqDis)) Then !Tom Schroeder soil-based local drop
                Do i=1,nElm !> no_voxels -> number of nodes in  a cuboid (imin)
                   If (no_voxels(i) .Eq. 0) Goto 40
                   Do j=1,no_voxels(i)
                      irecn = voxel_node(i,1,j)
                      isub = voxel_node(i,2,j)
                      corner=loc_Q(irecn,1:8,isub,ipl)
                      !coordinates of voxel
                      x1=xgrid(corner(1))
                      y1=ygrid(corner(1))
                      z1=zgrid(corner(1))
                      x2=xgrid(corner(8))
                      y2=ygrid(corner(8))
                      z2=zgrid(corner(8))
                      If (ave) Then
                         cp_mean(irecn,1,isub,ipl)=(x2+x1)/2
                         cp_mean(irecn,2,isub,ipl)=(y2+y1)/2
                         cp_mean(irecn,3,isub,ipl)=(z2+z1)/2
                      Elseif (eqDis) Then
                         numNodes_voxel(irecn,isub)=no_voxels(i) !> total root nodes of cuboid in no_voxels;  cubiods have same imin
                      Endif
                   Enddo
40                 Continue
                Enddo
             Endif
             Write(*,*)'t=',t
             Write(*,*)'t_begin=',t_begin,t_begin+dt
             Write(*,*)'old=',old
             !if (.not.(old)) then
             If (t .Le. t_begin+dt) Then
                !call cpu_time(t0)
                Write(*,*)'Matric flux potential Library loaded'
                PHtemp = 1e-5
                Call logspace(hx_min,PHtemp,nLibr,y)
                h_mat = -y
                h_mat2=(h_mat(2:Size(h_mat))+h_mat(1:Size(h_mat)-1))/2
                Call setmat_anaLibr(h_mat,theta,Kh,C)
                !> calculate matric flux potential integral K(h)dh -> numerical integration -> to linearize non-linear conductivity
                arr=Abs(h_mat(2:Size(h_mat))-h_mat(1:Size(h_mat)-1))*(Kh(2:Size(Kh))+Kh(1:Size(Kh)-1))/2
                cumsum(Size(arr)) = Sum(arr(1:Size(arr)))
                Phi_mat = cumsum
             Endif
          Endif
          !> total number of embrach. to the seed
          nBCn(ipl)=iBCn
       Endif
  

    If (lChem) Then
       Sbetac=0.0_dp
       VEl=dxGrid*dyGrid*dzGrid/6._dp
       Do iE=1,nElm
          Do iSE=1,5
             i=elmnod(iL(1,iSE,subN(iE)),iE) !> i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
             j=elmnod(iL(2,iSE,subN(iE)),iE)
             k=elmnod(iL(3,iSE,subN(iE)),iE)
             l=elmnod(iL(4,iSE,subN(iE)),iE)
             betcE=(betac(i)+betac(j)+betac(k)+betac(l))/4.
             Sbetac=Sbetac+betcE
          Enddo
       Enddo
       If (Sbetac.Gt.1.E-20_dp) Then
          Sbetac=Sbetac*VEl
          betac=betac/Sbetac
       Else
          Write(*,*)'Sbetac < 10^-20,  Sbetac = ',Sbetac
       Endif
    Endif

  End Subroutine SetupDou
!*******************************************************************************
!>### for a segment irecn between A and B : 
! - calculate the relative length of each subsegment: w_sub(irecn,isub,ipl)
! - calculate the relative surface of each subsegment: beta_weight(irecn,isub,ipl)
! - calculate the inverse dist. btw the center of the subs. and the 8 voxel nodes: w_dis(irecn,ic,isub,ipl) 
! - give the material around each subsegment  mat_Q ###

  Subroutine segment(xA,yA,zA,xB,yB,zB,ifoln,irecn,ipl,Weight,PrvSur)
    Use Typedef
    Use ParamData, Only: lChem
    Use GridData, Only:xgrid,ygrid,zgrid,dzgrid,dygrid,dxgrid,nex,ney,continu,betac,betaw
    Use DoussanMat, Only: nsub,loc_Q,mat_Q,w_sub,sum_dis,w_dis,cent,Intc,no_voxels,voxel_no,&
         voxel_node,indexValue,old,ave,eqdis,l_seg,l_sub,isubmax,transroot,cube_i, beta_weight
    Use RootData, Only: seglen,ibrseg,irecpr,segsur,lDou,lCou,lFed,lno_Archi,lKdrop
    Use MPIutils, Only: stop_program
    Use SolData, Only: MatNum
    Use RootGrowthNeighb, Only: Neighb
    Implicit None

    Integer(ap),Intent(in) :: ipl
    Integer(ap) ::iface,iFaceOld,isub,corner(8),ic,ifoln,irecn,reftransroot(2),imin
    Real(dp),Intent(in) :: xA,yA,zA,xB,yB,zB
    Real(dp):: x1,x2,y1,y2,z1,z2,xAp,yAp,zAp
    Real(dp) ::yint,xint,zint,xCent,yCent,zCent,weight
    Real(dp)::blengt,blengtTot,srface,PrvSur
    Logical :: SPLIT

    !> number of subsegment for node irecn
    !IF(continu) nsub(irecn,ipl)=0
    nsub(irecn,ipl)=0
    If (lDou.Or.((lCou.Or.lFed).And.(.Not.lno_Archi))) Then
       l_seg(irecn)=0
       blengtTot=0.00000001_dp
    Endif
    SPLIT=.True.
    !> initialization
    xAp=xA
    yAp=yA
    zAp=zA
    If(continu) Then
       If(ifoln.Ne.0) Then
          reftransroot=transroot(ifoln,:,nsub(ifoln,ipl),ipl)!> in Intc, transroot, cube_i, etc., the node from rootsys (-> not created by instec) is written in the last columns
       Else
          reftransroot=transroot(ifoln,:,1,ipl)
       Endif
    Endif

    splitloop: Do While(SPLIT)
       nsub(irecn,ipl)=nsub(irecn,ipl)+1
       isub=nsub(irecn,ipl)
       If (isub.Gt.isubmax) Then
          Print *,'Number of subsegments higher than maximum admitted, segment',irecn,'too long with report to the grid resolution'
          Call stop_program('')
       Endif
       !  find cuboid around the apical node of the segment:
       Call Neighb(xAp,yAp,zAp,corner,imin,.FALSE.)!imin : voxel containing the node
       cube_i(irecn,isub,ipl)=Int(imin)
       If (lDou.Or.(lCou.And.(.Not.lno_Archi))) Then
          Loc_q(irecn,1:8,isub,ipl)=corner
          If (((old).and.(lKdrop)).Or.(ave) .Or. (eqdis)) Then
             If (voxel_no(irecn,isub) .Eq. 0) Then !> remove double values from different branchings; if this root node is not processed process it
                voxel_no(irecn,isub)=imin !> voxel ID of a given root node adn subsegemnt
                no_voxels(imin)=no_voxels(imin)+1 !> number of root nodes in a voxel
                voxel_node(imin,1,no_voxels(imin))=irecn
                voxel_node(imin,2,no_voxels(imin))=isub
                If (no_voxels(imin) .Gt. indexValue) Then
                   Print*,'Parameter indexValue in Segment, representing number of nodes in a voxel, has to be set larger'
                   Call stop_program('')
                Endif
             Endif
          Endif
       Endif

       !> calculate cuboid\B4s corners coordinates
       !> Valid in all cases (continu or not) 
       x1=xgrid(corner(5))
       y1=ygrid(corner(5))
       z1=zgrid(corner(5))
       x2=x1+dxgrid
       y2=y1+dygrid
       z2=z1+dzgrid

       !> check if segment is completely included in the cuboid
       If (intsec(xAp,yAp,zAp,xB,yB,zB,ifoln,irecn,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) Then !> true when the root segment intersects the cubo\EFd
          SPLIT=.True.
       Else
          SPLIT=.False. !> just/still 1 loop and then no subsegmentation
          xInt=xB
          yInt=yB
          zInt=zB
       Endif
       !> correction of Ap position 
       If (isub.Gt.1) Then
          Select Case(iFaceOld)
          Case(1,3,5) 
             zAp=zAp+5.E-5_dp*dzGrid
          Case(2,4,6) 
             zAp=zAp-5.E-5_dp*dzGrid
          End Select
       Elseif (continu) Then
          !> The "downward" extremity of a segment is the following node while the other (downward) parts belong to the current node 
          Intc(ifoln,1,nsub(ifoln,ipl),ipl)=xAp
          Intc(ifoln,2,nsub(ifoln,ipl),ipl)=yAp
          Intc(ifoln,3,nsub(ifoln,ipl),ipl)=zAp
       Endif

       If (lDou.Or.((lCou.Or.lFed).And.(.Not.lno_Archi))) Then
          !> calculate (sub-)segment length and surface
          blengt=Sqrt((xInt-xAp)*(xInt-xAp)+(yInt-yAp)*(yInt-yAp)+ (zInt-zAp)*(zInt-zAp))
          blengtTot=blengtTot+blengt
          srface = 0.0_dp
          !TotSur = 0.0_dp
          If (irecn.Gt.0) Then 
             srface=blengt*segsur(irecn,ipl)/seglen(irecn,ipl)
             beta_weight(irecn,isub,ipl)=srface/Prvsur
          Endif


          !> calculate relative length of this (sub)segment
          If (irecn.Eq.0) Then
             w_sub(irecn,isub,ipl)=blengt/seglen(irecn+1,ipl)
          Else
             w_sub(irecn,isub,ipl)=blengt/seglen(irecn,ipl)
          Endif
          If ((.Not.(split)).And.(isub.Eq.1.)) Then
             w_sub(irecn,isub,ipl)=1._dp
          Endif
          If ((w_sub(irecn,isub,ipl).Gt.(1.-1.E-7)).Or.((.Not.(split)).And.(isub.Eq.1.))) Then
             w_sub(irecn,isub,ipl)=1._dp
          Elseif (w_sub(irecn,isub,ipl).Lt.1.E-7) Then
             w_sub(irecn,isub,ipl)=0._dp
          Endif
          l_sub(irecn,isub,ipl)=blengt
          l_seg(irecn)=l_seg(irecn)+blengt

          !> calculate (sub)segment center coordinates...
          xCent=xAp+(xInt-xAp)/2.
          yCent=yAp+(yInt-yAp)/2.
          zCent=zAp+(zInt-zAp)/2.
          cent(irecn,1,isub,ipl)=xCent
          cent(irecn,2,isub,ipl)=yCent
          cent(irecn,3,isub,ipl)=zCent

          !> and calculate the distance to each voxel node:
          sum_dis(irecn,isub,ipl)=0.0_dp
          Do ic=0,1 !Valid in all cases (continu or not)		 
             w_dis(irecn,4*ic+1,isub,ipl)=Sqrt((xCent-xgrid(corner(4*ic+1)))**2+(yCent-ygrid(corner(4*ic+1)))**2+(zCent-zgrid(corner(4*ic+1)))**2)
             w_dis(irecn,4*ic+2,isub,ipl)=Sqrt((xCent-xgrid(corner(4*ic+1))-dxgrid)**2+(yCent-ygrid(corner(4*ic+1)))**2+(zCent-zgrid(corner(4*ic+1)))**2)
             w_dis(irecn,4*ic+3,isub,ipl)=Sqrt((xCent-xgrid(corner(4*ic+1)))**2+(yCent-ygrid(corner(4*ic+1))-dygrid)**2+(zCent-zgrid(corner(4*ic+1)))**2)
             w_dis(irecn,4*ic+4,isub,ipl)=Sqrt((xCent-xgrid(corner(4*ic+1))-dxgrid)**2+(yCent-ygrid(corner(4*ic+1))-dygrid)**2+(zCent-zgrid(corner(4*ic+1)))**2)
          End Do

          Do ic=1,8
             If (w_dis(irecn,ic,isub,ipl).Lt.1.E-20_dp) w_dis(irecn,ic,isub,ipl)=1.E-20_dp
             w_dis(irecn,ic,isub,ipl)=1._dp/w_dis(irecn,ic,isub,ipl)! inverse-distance weight 
             sum_dis(irecn,isub,ipl)=sum_dis(irecn,isub,ipl)+w_dis(irecn,ic,isub,ipl)
             mat_Q(irecn,ic,isub,ipl)=MatNum(corner(ic)) !define material around roots
             If (lFed) betaw(corner(ic))=betaw(corner(ic))+w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*Weight*(blengt/blengtTot)
             If (lChem) betac(corner(ic))=betac(corner(ic))+w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*Weight*(srface/PrvSur)
          End Do
       Endif

       !> save Ap position (inside of the soil domain -> xmax and ymax of the soil boudaries not included)
       If (continu.And.(isub.Gt.1)) Then
          If (xAp.Eq.Minval(xGrid)+nex*dxgrid) Then
             xAp=xAp-nex*dxgrid
             transroot(irecn,1,isub,ipl)=transroot(irecn,1,isub,ipl)-1
          Endif
          If (yAp.Eq.Minval(yGrid)+ney*dygrid) Then
             yAp=yAp-ney*dygrid
             transroot(irecn,2,isub,ipl)=transroot(irecn,2,isub,ipl)-1
          Endif
          Intc(irecn,1,isub-1,ipl)=xAp            
          Intc(irecn,2,isub-1,ipl)=yAp
          Intc(irecn,3,isub-1,ipl)=zAp
       Endif

       If (SPLIT) Then
          !> preparation for next loop
          xAp=xInt
          yAp=yInt
          zAp=zInt
          iFaceOld=iFace
          Select Case (iFace)
          Case(1) 
             zAp=zInt-5.E-5_dp*dzGrid
             If(continu) transroot(irecn,:,isub+1,ipl)=reftransroot
          Case(2) 
             zAp=zInt+5.E-5_dp*dzGrid
             If(continu) transroot(irecn,:,isub+1,ipl)=reftransroot
          Case(3) 
             yAp=yInt-5.E-5_dp*dyGrid
             If (continu) Then
                If (yAp.Lt.Minval(ygrid)) Then
                   yAp=yAp+ney*dygrid
                   transroot(irecn,1,isub+1,ipl)=reftransroot(1)
                   transroot(irecn,2,isub+1,ipl)=reftransroot(2)+1
                   reftransroot=transroot(irecn,:,isub+1,ipl)
                Else
                   transroot(irecn,:,isub+1,ipl)=reftransroot

                Endif
             Endif
          Case(4) 
             yAp=yInt+5.E-5_dp*dyGrid
             If (continu) Then
                If (yAp.Gt.Minval(ygrid)+ney*dygrid) Then
                   yAp=yAp-ney*dygrid
                   transroot(irecn,1,isub+1,ipl)=reftransroot(1)
                   transroot(irecn,2,isub+1,ipl)=reftransroot(2)-1
                   reftransroot=transroot(irecn,:,isub+1,ipl)
                Else
                   transroot(irecn,:,isub+1,ipl)=reftransroot

                Endif
             Endif
          Case(5) 
             xAp=xInt-5.E-5_dp*dxGrid
             If (continu) Then
                If (xAp.Lt.Minval(xgrid)) Then
                   xAp=xAp+nex*dxgrid
                   transroot(irecn,1,isub+1,ipl)=reftransroot(1)+1
                   transroot(irecn,2,isub+1,ipl)=reftransroot(2)
                   reftransroot=transroot(irecn,:,isub+1,ipl)
                Else
                   transroot(irecn,:,isub+1,ipl)=reftransroot

                Endif
             Endif
          Case(6) 
             xAp=xInt+5.E-5_dp*dxGrid
             If (continu) Then
                If (xAp.Gt.Minval(xgrid)+nex*dxgrid) Then
                   xAp=xAp-nex*dxgrid
                   transroot(irecn,1,isub+1,ipl)=reftransroot(1)-1
                   transroot(irecn,2,isub+1,ipl)=reftransroot(2)
                   reftransroot=transroot(irecn,:,isub+1,ipl)
                Else
                   transroot(irecn,:,isub+1,ipl)=reftransroot
                Endif
             Endif
          End Select
       Elseif (continu.And.(isub.Gt.1)) Then
          !> sort transroot in the same order as intc
          transroot(irecn,:,isub+1,ipl)=transroot(irecn,:,1,ipl)
          transroot(irecn,:,1:isub+1,ipl)=transroot(irecn,:,2:isub+2,ipl)

       Endif
    End Do splitloop

    !> Ap was on a face of the cube containig A. If Ap was on a boundary too, it is displaced out of the domain (because of the +/-5.E-5...).
    !> Then -> translation. Ap is now closer to B (on the same side and inside of the domain).

    !> save branch basis position 
    If (irecn.Gt.0) Then
       If (continu.And.(ibrseg(irecpr(irecn,ipl),ipl).Ne.ibrseg(irecn,ipl))) Then
          Intc(irecn,1,isub,ipl)=xB
          Intc(irecn,2,isub,ipl)=yB
          Intc(irecn,3,isub,ipl)=zB
       Endif
    Endif
    If (lDou) Then
       !> correction of inacurrate seglen
       If (nsub(irecn,ipl).Eq.1) Then !nosplit
          w_sub(irecn,1,ipl)=1._dp
       Else
          Do isub=1,nsub(irecn,ipl)
             w_sub(irecn,isub,ipl)=l_sub(irecn,isub,ipl)/l_seg(irecn)
          End Do
       Endif
    Endif
    Return
  End Subroutine segment
  !****************************************************************
  !> ### estimate current BC for root system when using Doussan Sink term ###
  Subroutine SetBCroot(t,BCr,BCtp)
    Use TypeDef
    Use PlntData, Only : tBCr,typeBCr,BCroot,nBCr 
    Implicit None

    Integer(ap),Intent(out) :: BCtp
    Integer(ap) :: ifc
    Real(dp), Intent(in) :: t
    Real(dp), Intent(out) ::BCr

    !> calculate current root UBC-value from input BC(time)-function:
    ifc=0
201 ifc=ifc+1
    If (ifc.Gt.nBCr) Then !> beyond the last user-defined BC value it remains constant
       BCr=BCroot(nBCr)
       BCtp=typeBCr(nBCr)
    Else
       If (t.Ge.tBCr(ifc)) Goto 201 !> try to find the first tBCR wwhich is beyond current t (>=1)
       If (ifc.Eq.1) Then
          BCr=BCroot(ifc)
          BCtp=typeBCr(ifc)
       Elseif ((typeBCr(ifc-1)).Eq.typeBCr(ifc)) Then
          BCr=BCroot(ifc-1)+(BCroot(ifc)-BCroot(ifc-1))*(t-tBCr(ifc-1))/(tBCr(ifc)-tBCr(ifc-1))
          BCtp=typeBCr(ifc-1)
       Else
          BCr=BCroot(ifc-1)
          BCtp=typeBCr(ifc-1)
       Endif
    Endif
    If (BCtp.Eq.2) Then
       BCr=-Abs(BCr) !> always negative flow at the collar
    Endif
    Return
  End Subroutine SetBCroot
  !**************************************************************************
  !> ### Calculates stomatal closure ###
  Subroutine Rootstress(PHtop,BCtop,BCtptop,iterBC,BC_switch,Jcol,ipl)
    Use Typedef
    Use DoussanMat, Only: curr_BCtp,stressBC,hx_min,stresfun,stresval1,stresval2
    Use GridData, Only: RelEps,epslonR,factorRelEps
    Use PlntData, Only: a_r,a_1,a_2
    Use RootData, Only: concol,csign_notrans,PH_crit
    Implicit None

    Integer(ap), Intent(inout):: BCtptop
    Integer(ap), Intent(in):: ipl
    Real(dp), Intent(inout):: BCtop
    Real(dp), Intent(in):: PHtop,Jcol
    Real(dp):: reduc,del_PHr
    Logical, Intent(out) :: BC_switch,iterBC

    !> check if the root collar abs(PH) is larger than abs(hx_min)+tolerance and adapt the collar BC
    If (curr_BCtp(ipl)==2) Then
       If (RelEps) Then 
          del_PHr=-Max(Abs(hx_min/factorRelEps),epslonR)
       Else
          del_PHr=-epslonR
       Endif
       If ((stresfun.Eq.1).And.(PHtop<hx_min+del_PHr)) Then
          !> + stresfun=1: top node at lower PH than allowed: start of stressed conditions
 !         Print *,'stress in the collar xylem: change to PH BC, PHtop=',PHtop,' is lower than criterion ',hx_min,' +',del_PHr
          write(*,'(a)',advance='no') '||' !when stress occurs a vertical bar appears on the screen
          stressBC=.True.
          BCtptop=1
          BCtop=hx_min
          BC_switch=.True.
          iterBC=.True.
       Elseif ((stresfun.Eq.2).And.(PHtop.Le.stresval1)) Then !> + stresfun=2: linear decrease
          reduc=(PHtop-stresval2)/(stresval1-stresval2)
          If (reduc.Gt.1) Then
             reduc=1
          Elseif (reduc.Lt.0) Then
             reduc=0
          Endif
          BCtop=-Abs(Jcol)*reduc
          BCtptop=2
          stressBC=.True.
          iterBC=.True.
          Print *,'stresfun=2; flux reduction of factor',reduc
       Elseif (stresfun.Eq.3) Then !> + stresfun=3: Tuzet function
          reduc=(1+Exp(stresval1*stresval2))/(1+Exp(stresval1*(stresval2-PHtop)))
          BCtop=-Abs(Jcol)*reduc
          BCtptop=2
          Print *,'stresfun=3; flux reduction of factor',reduc
          If (reduc.Lt.1) Then
             BC_switch=.True.
             iterBC=.True.
          Else
             BC_switch=.False.
             iterBC=.False.
          End If
       Elseif (stresfun.Eq.4) Then !> + stresfun=4: additional signaling
          reduc= a_r + (1-a_r)*Exp(-concol*a_1*Exp(a_2*(Abs(PHtop)-Abs(PH_crit))))
          If (reduc.Lt.1) Then
             Print *,'stresfun=4; flux reduction of factor',reduc
             BCtop=-Abs(Jcol)*reduc
             BCtptop=2
             BC_switch=.True.
             iterBC=.True.
          Else
             BC_switch=.False.
             iterBC=.False.
          Endif
       Elseif (stresfun.Eq.5) Then !> + stresfun=5: additional signaling without particle tracker
          reduc= a_r + (1-a_r)*Exp(-csign_notrans*a_1*Exp(a_2*(Abs(PHtop)-Abs(PH_crit))))
          If (reduc.Lt.1) Then
             BCtop=-Abs(Jcol)*reduc
             BCtptop=2
             BC_switch=.True.
             iterBC=.True.
          Else
             BC_switch=.False.
             iterBC=.False.
          Endif
       Else
          BC_switch=.False.
          iterBC=.False.
       Endif

    Elseif ((curr_BCtp(ipl)==1).And.(Abs(Jcol)>1.00001*Abs(BCtop))) Then !factor to avoid oscillation between stress and no stress
       !> check if end of stress (curr_BCtp is PH with PH is hxmin but BCtp_new is still flux)
!       Print *,'end of stress conditions, collar flux=',Abs(Jcol),' is larger than the prescribed flux', Abs(BCtop)
       !  BCtptop=BCnew is already of type 2!
       stressBC=.False.
       BC_switch=.False.
!       Print *,'info',stressBC,BCtptop,BCtop, BC_switch,iterBC

    Else
       BC_switch=.False.
       iterBC=.False.
!       If (stressBC) Print *,'no change: stress=',stressBC
    Endif

  End Subroutine Rootstress
 !*************************************************************
 !> ### allocate radial and axial root conductivities to root segments acc. to their age ###
  Subroutine ConductRoot(t,ipl)
    Use TypeDef
    USE ParamData, ONLY: maxnod
    Use DoussanMat, Only : nLr,Lr,Lrroot,ageLr,nKh,Khr,Khroot,ageKh,Khr_pot,Lr_pot,lcavit,transroot
    Use RootData, Only: nrec,timorg,ordseg,maizeroottyp,loliumroottyp,wheatroottyp,lGap,lAQPc,lKdrop,xs,ys,zs,ltwo_grids
    Use SolData, Only: MatNum, MatNum2,nmat,par,l_elmMacro
    Use GridData, Only: n_neigh,nex,ney,dxgrid,dygrid,continu,nPt
    Use GridData2, Only:nPt2
    Use RootGrowthNeighb, Only: Neighb
    Implicit None

    INTEGER (ap), INTENT(in):: ipl
    Integer(ap) :: irecn,iage,typ,corner(8),imin,iMat,nM,c
    Real(dp), Intent (in) :: t
    Real(dp) :: segage,xC,yC
    Integer(ap),Allocatable,Dimension(:):: MatNum_

    IF (ltwo_grids) THEN
      ALLOCATE (MatNum_(nPt2))
    ELSE
      ALLOCATE(MatNum_(nPt))
    ENDIF

    !> calculate segment lateral and long. conductivity according to age:
    Do irecn=1,nrec(ipl)
       segage=t-timorg(irecn,ipl)
       typ=ordseg(irecn,ipl) 
       If (maizeroottyp) Then
          If (typ.Lt.12) Then  !In RootTyp, maize root types below 12 are principal roots (12 and 13 are lateral) 
             typ=1
          Else
             typ=2
          Endif
       Elseif (loliumroottyp) Then
          If (typ.Lt.3) Then  !In RootTyp, lolium root types below 3 are principal roots
             typ=1
          Elseif (typ.Lt.5) Then
             typ=2
          Else
             typ=3
          Endif
       Elseif (wheatroottyp) Then
          If (typ.Lt.19) Then  !In RootTyp, wheat root types below 19 are principal roots 
             typ=1
          Elseif (typ.Lt.20) Then
             typ=2
          Else
             typ=3
          Endif
       Else
          If (typ.Eq.0) Then  !type 0 is the seed+small segment
             typ=1
          Elseif (typ.Gt.3) Then!no more than 3 root types
             typ=3
          Endif
       Endif
       !> radial conductivity matrix
       iage=1
       Do While (ageLr(typ,iage,ipl).Le.segage)
          iage=iage+1
       Enddo
       If (iage>nLr(typ,ipl)) Then
          Lr(irecn,ipl)=LrRoot(typ,nLr(typ,ipl),ipl)
       Elseif (iage==1) Then
          Lr(irecn,ipl)=LrRoot(typ,iage,ipl)
       Else
          Lr(irecn,ipl)=LrRoot(typ,iage-1,ipl)+(LrRoot(typ,iage,ipl)-LrRoot(typ,iage-1,ipl))*(segage-ageLr(typ,iage-1,ipl))/&
               (ageLr(typ,iage,ipl)-ageLr(typ,iage-1,ipl))
       Endif

       IF (ltwo_grids) THEN
          MatNum_ = MatNum2
       ELSE 
          MatNum_ = MatNum
       END IF 
	   
	   
!> if root segments are within a non-conductive soil material --> set radial conductivity to zero, or scale Lr according to the amount of non-cond. material
       !> macropore
       If(nMat.Gt.1) Then
          Do iMat=1,nMat 
             If (par(6,iMat) .Le. 1.0E-05_dp) Then !non-conductive
                If(continu) Then
                   xC=xs(irecn,ipl)+transroot(irecn,1,1,1)*nex*dxgrid
                   yC=ys(irecn,ipl)+transroot(irecn,2,1,1)*ney*dygrid
                   Call Neighb(xC,yC,zs(irecn,ipl),corner,imin,ltwo_grids)
                Else
                   Call Neighb(xs(irecn,ipl),ys(irecn,ipl),zs(irecn,ipl),corner,imin,ltwo_grids)
                End If
                If (Any(MatNum_(corner) .Eq. iMat))  Lr(irecn,ipl) = 0._dp 
             Elseif(par(11,iMat).Le.1E-5_dp)  Then  !macropore
                If(continu) Then
                   xC=xs(irecn,ipl)+transroot(irecn,1,1,1)*nex*dxgrid
                   yC=ys(irecn,ipl)+transroot(irecn,2,1,1)*ney*dygrid
                   Call Neighb(xC,yC,zs(irecn,ipl),corner,imin,ltwo_grids)
                Else
                   Call Neighb(xs(irecn,ipl),ys(irecn,ipl),zs(irecn,ipl),corner,imin,ltwo_grids)
                End If
                nM=8  !> number of mat.nodes within bulk material
                Do c=1,8
                   If(MatNum_(corner(c)) .Eq. iMat) nM=nM-1
                End Do
                Lr(irecn,ipl) = Lr(irecn,ipl)*nM/8._dp 
                If (l_elmMacro(imin) .And. n_neigh(imin).Lt.1)  Lr(irecn,ipl) = 0._dp
             End If
          End Do
       End If

       !> axial conductance matrix
       iage=1
       Do While (ageKh(typ,iage,ipl).Le.segage)
          iage=iage+1
       Enddo
       If (iage>nKh(typ,ipl)) Then
          Khr(irecn,ipl)=KhRoot(typ,nKh(typ,ipl),ipl)
       Elseif (iage==1) Then
          Khr(irecn,ipl)=Khroot(typ,iage,ipl)
       Else
          Khr(irecn,ipl)=KhRoot(typ,iage-1,ipl)+(KhRoot(typ,iage,ipl)-KhRoot(typ,iage-1,ipl))*(segage-ageKh(typ,iage-1,ipl))/&
               (ageKh(typ,iage,ipl)-ageKh(typ,iage-1,ipl))
       Endif
       If (lcavit) Khr_pot(irecn,ipl)=Khr(irecn,ipl)!only one plant possible for theses processes
       If (lGap.Or.lAQPc.Or.lKdrop) Lr_pot(irecn)=Lr(irecn,ipl)!only one plant possible for theses processes
    End Do

  End Subroutine ConductRoot
 !******************************************************************************
 !> ### solve DOUSSAN flow equations within the xylem ###
  Subroutine SolveRoot(t,dt,it1,iter_root)
    Use TypeDef
    Use SparseMatrix
    Use NumericalRecipes
    Use Watfun
    Use Paramdata,Only: MaxNod, pi, ldirect
    Use GridData, Only: sink,RootSk, sink_cube,betac_cube,dxGrid,dyGrid!,SSF
    Use PlntData, Only :Tpot,Tact
    Use DoussanMat, Only: stressBC,Q_bc1,Q_bc2,Qi,PHs,Qd,PHr,GH,nrecOld,nBCn,bcr_usr,counter2, &
         bctp_usr,curr_bcr,curr_bctp,Khr,NBC_irecn,cent,Lr,ave,old,Joutr,PH_root_sub, &
         cp_mean,oldT,plantmatrix,nplant,lcavit,cavitb,cavitc,Khr_pot,Lr_pot,solveroot_call,&
         KRhizo,PHsri,voxel_node,no_voxels,mat_seg,vox_seg!,w_sub
    Use RootData, Only: nrec,nrec_m,seglen,lCalloc,lGap,lAQPc,nplant,fac_plant,lKdrop,segrad,&
         TypeKdrop,segsur,fac_contact
    Use SolData, Only: par
    Use tmctrl, Only: t_begin, tlevel,dtroot
    Use MPIutils, Only: stop_program
    Use PlantGrowth, Only: Stress, SetTp
    Implicit None
	
    Integer(ap),Intent(in) ::iter_root
    Real(dp) ,Intent(in) :: t,dt
    Logical,Intent(in) :: it1
	
    Integer(ap) ::old_BCtp,iter_loc,itol,itmaxl,ki,imin,isub_tmp,irecn_tmp
    Integer(ap) :: iprvn,irecn,ibranch,BCtp_new,err,ipl,num_elem,nrecmax
    Real(dp)::PHr_tot(1:nrec_m+1,nplant),PH_direct(1:nrec_m+1,nplant), err_loc,Qd2(1:2*(nrec_m+1))
    Real(dp)::maxerr=1.e-10_dp,BCr_new,BCr_old,rootsk2
    Real(dp) :: rs,concrs,phi_b,phi_sri,dist2r,rbulk
    Real(dp), Save::temp_phr
    Logical :: repet,iterBC,BC_switch
    Logical ,Save::firstcall=.True.

!> initialisation
    DO ipl=1,nplant
       If(t.LE.dtRoot) phr_tot(1:nrec(ipl)+1,ipl)=0._dp
    ENDDO
    If (ave)  counter2=0

    solveroot_call=solveroot_call+1     
    nrecmax = MAXVAL(nrec)
    sink=0._dp
    RootSk=0.0_dp
    RootSk2=0.0_dp
    sink_cube=0.0_dp
    betac_cube=0.0_dp

! Loop over different plants
    Do ipl=1,nplant

      If (lGap.Or.lAQPc.Or.(lKdrop.and.((TypeKdrop).eq.1))) Then
             Do irecn=1,nrec(ipl)
                Lr(irecn,ipl)=Lr_pot(irecn)*Gap(PHs(irecn,ipl))*AQPc(PHs(irecn,ipl),ipl)
                !> Gap is 1 in case lGap is .false., same for AQPc. If we wanted to be exact, the total factor should be 
                !> 2*Gap*AQPc/(Gap+AQPc), since the two conductances are in series
                           rbulk=sqrt((dxGrid*dyGrid)/(1*pi))! bulk radius
                If (lKdrop.and.(TypeKdrop.eq.1).and.((Joutr(irecn,ipl)/segsur(irecn,ipl)).GT.0.0000)) Then !only if uptake, not if efflux
                   phi_b=Fmfp_soiltab(PHs(irecn,ipl),mat_seg(irecn,ipl),Par(:,1))
                   phi_sri=phi_b+(Joutr(irecn,ipl)/(segsur(irecn,ipl)*fac_contact))*segrad(irecn,ipl)*log(segrad(irecn,ipl)/rbulk)
                   imin = vox_seg(irecn,ipl)! voxel ID?
                   dist2r=sqrt((dxGrid*dyGrid)/(no_voxels(imin)*pi))! averaged distance btw roots in that voxel
                   Do ki=1,no_voxels(imin)!when more than one root in a voxel, additive approach of Doussan
                       irecn_tmp = voxel_node(imin,1,ki)
                       isub_tmp = voxel_node(imin,2,ki)
                       If ((Joutr(irecn_tmp,ipl).GT.0.000).AND.(irecn_tmp.NE.irecn)) Then
                            phi_sri=phi_sri+(Joutr(irecn_tmp,ipl)/(segsur(irecn_tmp,ipl)*fac_contact))*segrad(irecn_tmp,ipl)*log(dist2r/rbulk)
                       endif
                   Enddo
                   If (Abs(phi_sri-phi_b).GT.0.000) Then
                     If (solveroot_call.EQ.1) Then
                       PHsri(irecn,ipl)=Fh_from_mfp_soiltab(phi_sri,mat_seg(irecn,ipl))
                     Else
                       PHsri(irecn,ipl)=MAX(Fh_from_mfp_soiltab(phi_sri,mat_seg(irecn,ipl)),PHr(irecn+1,ipl)) !!
                     Endif
                       KRhizo(irecn,ipl)=(Joutr(irecn,ipl)/segsur(irecn,ipl))/(PHs(irecn,ipl)-PHsri(irecn,ipl))
                       Lr(irecn,ipl)=Lr(irecn,ipl)*KRhizo(irecn,ipl)/(Lr(irecn,ipl)+KRhizo(irecn,ipl))
                   Else! no significant change
                       PHsri(irecn,ipl)=PHs(irecn,ipl)
                       KRhizo(irecn,ipl)=0
                   Endif
                Else
                   KRhizo(irecn,ipl)=0
                   PHsri(irecn,ipl)=PHs(irecn,ipl)
                Endif
             Enddo
      Else
         Do irecn=1,nrec(ipl)
           PHsri=PHs
         Enddo
      Endif


!> Change of axial/radial conductance due to cavitation/AQP/Kdrop/gaps during the simulation/update once at each time step
       If ((lcavit.Or.lGap.Or.lAQPc.OR.(lKdrop.And.(TypeKdrop.eq.1))) )Then
            call UpdateDou(ipl)
       Endif

!> Set current plant BC
       BCr_old=BCr_usr(ipl)
       If (.Not.lCalloc) Then
          Call SetBCroot(t,BCr_usr(ipl),BCtp_usr(ipl))
          BCr_usr(ipl)=BCr_usr(ipl)*fac_plant(ipl)
       Else
          If (ipl.Ge.2) Call stop_program('Assimilate allocation currently doesnt work with multiple plants')
          Call Stress(rs,concrs)
          Call SetTp(t,rs,concrs,ipl)
          BCr_usr(ipl)=Abs(Tpot(ipl))
          BCtp_usr(ipl)=2
       Endif
       !> calculation of Tpot if BCtp==2
       If (BCtp_usr(ipl)==2) Tpot(ipl)=+Abs(BCr_usr(ipl))
       !> if stress at previous time step, then BCtp new=2 and currBctp=1
       If ((stressBC).And.((BCtp_usr(ipl)==1).Or.(BCr_old/=BCr_usr(ipl)))) Then
          !tocheck!!!with tuzet it will not work
          !there was a stress with BC=2 but in the meantime user BC changed
          stressBC=.False.
       Endif
       BCr_new=BCr_usr(ipl)
       BCtp_new=BCtp_usr(ipl)
       iterBC=.True.
       BC_switch=.False.
       Do While (iterBC)
          !> adapt BC
          If (firstcall.Or.((BCr_new/=curr_BCr(ipl)).And.(.Not.(stressBC))).Or.((BCtp_new/=curr_BCtp(ipl)).And.(.Not.(stressBC))).Or.(BC_switch)) Then
             !> updating BC
             curr_BCr(ipl)=BCr_new
             old_BCtp=curr_BCtp(ipl)
             curr_BCtp(ipl)=BCtp_new
             !> updating the BC part of the matrices Qbc and Cd
             iprvn=0 !> seed=first segment
             If (curr_BCtp(ipl)==2.And.ldirect) Then  !> 1. Flux and direct solver
                If (.Not.tlevel) Then
                   Do ibranch=1,nBCn(ipl)  !> all the branches connected to the seed
                      irecn=nBC_irecn(ibranch) !connected to the seed
                      err=SM_set(plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.True.)!position (1,1)
                      err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,0._dp,.True.)!position (2,1)
                      err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,0._dp,.True.)!position (1,2)
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                   Enddo
                Endif
                Do ibranch=1,nBCn(ipl)  !all the branches connected to the seed
                   irecn=nBC_irecn(ibranch)
                   temp_phr=-(curr_BCr(ipl)/Khr(irecn,ipl)*seglen(irecn+1,ipl)+200.)
                   Qi(iprvn,ipl)=0._dp !> from seed has to be given, from irecn is already stated in DoussanMat
                   Q_bc1(iprvn,ipl)=temp_Phr*.1+GH(iprvn,ipl) !> use ph from last time step
                   Q_bc2(iprvn,ipl)=temp_Phr+GH(iprvn,ipl)
                   Q_bc1(irecn,ipl)=(temp_Phr*.1+GH(iprvn,ipl))*Khr(irecn,ipl)/seglen(irecn,ipl) !> all the PH for root was given in total head!!
                   Q_bc2(irecn,ipl)=(temp_Phr+GH(iprvn,ipl))*Khr(irecn,ipl)/seglen(irecn,ipl)
                Enddo

             Else If (curr_BCtp(ipl)==2.And..Not.ldirect) Then  !> 2. Flux and iter solver
                If (old_BCtp/=curr_BCtp(ipl)) Then
                   Do ibranch=1,nBCn(ipl)  !all the branches connected to the seed
                      irecn=nBC_irecn(ibranch)!connected to the seed
                      If (ibranch==1) Then
                         repet=.True.!delete previous value
                      Else
                         repet=.False.!add to previous value
                      Endif
                      err=SM_set(plantmatrix(ipl),iprvn+1,iprvn+1,Khr(irecn,ipl)/seglen(irecn,ipl),repet)  !position (1,1)
                      err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),repet) !position (1,2)
                      err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),repet) !position (2,1)
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                      Q_bc1(irecn,ipl)=0._dp
                   Enddo
                   !> from seed has to be given, from irecn is already stated in DoussanMat
                   Qi(iprvn,ipl)=0._dp
                Endif
                Q_bc1(iprvn,ipl)=curr_BCr(ipl) !Q(0,0) a flow imposed

             Else If (curr_BCtp(ipl)==1) Then  !PH
                If (old_BCtp/=curr_BCtp(ipl)) Then
                   err=SM_set(plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.True.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Endif
                Qi(iprvn,ipl)=0._dp
                !> all the PH for root was given in total head!!
                Q_bc1(iprvn,ipl)=(curr_BCr(ipl)+GH(iprvn,ipl)) 
                Do ibranch=1,nBCn(ipl)
                   irecn=nBC_irecn(ibranch) !connected to the seed
                   If (old_BCtp/=curr_BCtp(ipl)) Then
                      err=SM_set(plantmatrix(ipl),iprvn+1,iprvn+1,1._dp,.True.)
                      err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,0._dp,.True.)
                      err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,0._dp,.True.)
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                   Endif
                   Q_bc1(irecn,ipl)=(curr_BCr(ipl)+GH(iprvn,ipl))*Khr(irecn,ipl)/seglen(irecn,ipl) !all the PH for root was given in total head!!
                   Q_bc2 = 0._dp
                End Do
             Endif
          Endif

          If (firstcall) firstcall=.False.

          If (oldT) Goto 201
          If (t.Eq.t_begin+dt .And. iter_root .Gt. 1 .Or. t.Gt.t_begin+dt) Then
             Call analytical_approach(ipl)
          Else
             old=.True.
          Endif
          !> for t=tstep. Initially no root flow is calculated, and not all the boundary conditions for the analytical problem are available, so solve
          !> first time step using the averaging way
          !> find current PH soil and create Qd2
201       If (old) Then
             Call average_method(ipl,iter_root)
          Endif
          If (nrec(ipl).Eq.1) Then!only seed
             If (curr_BCtp(ipl)==2) Then
                Q_bc1=Q_bc1 ! this must be checked MJ dec2008
             Endif
          Else
             !> put in format (1:nRec+1) for solving the system
             Qd2(1:nrec(ipl)+1)=Qd(0:nrec(ipl),ipl)
             Qd2(nrec(ipl)+2:) = Qd(nrec(ipl)+1:,ipl)
             PHr_tot(1:nrec(ipl)+1,ipl)=PHr(1:nrec(ipl)+1,ipl)+GH(0:nrec(ipl),ipl) !total water potential
             If (it1) Then
                PHr_tot(1:nrec(ipl)+1,ipl)=PHs(0:nrec(ipl),ipl)   !> initial xylem PH -> could also be zero but then convergence should take longer (this should be the common case)
             Endif
             If (curr_BCtp(ipl)==1) Then
                PHr_tot(1,ipl)=curr_BCr(ipl)+GH(0,ipl) !> bound value should always be incorporated if a PH is imposed
             Endif
             itol=3
             itmaxl=(nrec(ipl)+1)*3 !> maximum iterations for biconjugate gradient-> needed in case of root growth
             !> solve with biconjugated gradient method (Numerical recipies)

             If (ldirect) Then
                num_elem= SM_numberOfElements(plantmatrix(ipl))
                Call SolveRootDirect(ipl,num_elem,Qd2,Ph_direct,temp_phr)
                PHr_tot(1:nrec(ipl)+1,ipl) = ph_direct(1:nrec(ipl)+1,ipl)
             Else
                Call linbcg(nrec(ipl)+1,Qd2(1:nrec(ipl)+1),PHr_tot(1:nrec(ipl)+1,ipl),itol,maxerr,itmaxl,iter_loc,err_loc,ipl)
                If (err_loc.Le.maxerr) Then
                   nrecOld=nrec(ipl)
                Else
                   Print *,'No convergence in the root system'
                   Call stop_program('')
                Endif
             Endif
          Endif
		  !> calculate axial flow and sink for Doussan (MJ RSWMS10)
          call Snk_Doussan(PHr_tot,dt,t,iter_root,ipl)

          !> gravitational component subtracted again to represent data in outRoo.x without gravity
          PHr(1:nrec(ipl)+1,ipl)=PHr_tot(1:nrec(ipl)+1,ipl)-GH(0:nrec(ipl),ipl)
          PHs(0:nrec(ipl),ipl)=PHs(0:nrec(ipl),ipl)-GH(0:nrec(ipl),ipl)      

          If (lcavit) Then
             Khr(1:nrec(ipl),ipl)=Khr_pot(1:nrec(ipl),ipl)*Exp(-(-PHr(2:nrec(ipl)+1,ipl)/cavitb)**cavitc)
          Endif

          If (ave) Then
             PH_root_sub(0:nrec(ipl),:)=PH_root_sub(0:nrec(ipl),:)-cp_mean(0:nrec(ipl),3,:,ipl)
          Elseif (.Not.(ave)) Then
             PH_root_sub(0:nrec(ipl),:)=PH_root_sub(0:nrec(ipl),:)-cent(0:nrec(ipl),3,:,ipl)
          Endif
          If (.Not.lCalloc) Then
             If (((curr_BCtp(ipl).Eq.2).Or.(stressBC)).And.(.Not.(BC_switch))) Then
                !> check stress a the collar in case of BC flux
                !> check flux at the collar in case of PH BC+stress
                !> if swicth at previous iterBC, no check
                Call Rootstress(PHr(1,ipl),BCr_new,BCtp_new,iterBC,BC_switch,Tact(ipl),ipl)
                !> change BCtp_new=currBCtp to 1 if stress occurs
             Else
                iterBC=.False.
                BC_switch=.False.
             Endif
          Endif
       Enddo ! do while iterBC !actual transpiration
    Enddo ! plant loop
    Return
  End Subroutine SolveRoot
 !******************************************************************************
 !> ### calculate Sink term for Doussan model### (MJ2020)
  Subroutine Snk_Doussan(PHr_tot,dt,t,iter_root,ipl)
    Use TypeDef
    !Use SparseMatrix
    !Use NumericalRecipes
    !Use Watfun
	
    Use GridData, Only: betaw,RootSk,dxgrid,dygrid,dzgrid,n_neigh,sink,sink_cube,&
        betac_cube,betac_cube2,MacroList,Wn
    Use PlntData, Only : Tact
    Use DoussanMat, Only: curr_bctp,Khr,PHs,Q_bc1,axialRootFlow,w_dis,nplant,&
        veloRoot,Jintot,Qi,sinkR,nsub,cube_i, w_sub,Joutr,beta_weight,loc_Q,sum_dis
    Use RootData, Only: nrec,nbr,ibrseg,seglen,irecpr,crossSectionSeg,nrec_m
    Use SolData, Only: lMacro,l_elmMacro
    Use tmctrl, Only: t_begin
    Implicit None
	
    Integer(ap) ::ibr,ipl,irecn,ifoln,iprvn,iter_root,c_i,isub,i_dummy
    Integer(ap) :: ic,i_n,corner_ic
    Real(dp),Intent(in) :: PHr_tot(1:nrec_m+1,nplant),dt,t
    Real(dp) :: sink_dummy,betac_dummy
    Real(dp) :: totalAxialFlow=0.0_dp,Voltot,Jcol,checksub,dVolSk
    Logical :: n_apex,run
  
!>Calculate axial root flow: for each branch cumulative flow is calculated
    Branchloop: Do ibr=1,nbr(ipl)
!print*,'ibr',ibr
!print*,'nbr',nbr
        n_apex=.False.
		!> find the tip segment of the branch 'ibr'
        irecn=nrec(ipl)
        Do While (ibrseg(irecn,ipl).Ne.ibr)
           irecn=irecn-1
        End Do
        !> the first one we find is an apex
        If (seglen(irecn,ipl)<1.E-20) Then !> skip this segment too small to be taken
           ifoln=irecn !> following node ID
           irecn=irecpr(irecn,ipl) !> current node ID
        Else
           n_apex=.True.!> ifoln does not exist
        Endif
        If (irecn==0) Then !> there exists a branch ibr but not yet any segment!
           run=.False.
        Else
           run=.True.
        Endif
        !> then the rest of the branch up to the seed or the embranchment
        Do While (run)
           !> "upper" node
           iprvn=irecpr(irecn,ipl)
           axialRootFlow(irecn,ipl) = Khr(irecn,ipl)*(PHr_tot(irecn+1,ipl)-PHr_tot(iprvn+1,ipl))/seglen(irecn,ipl)
           !> if reached the seed or the embranchement => change branch
           If (iprvn==0.Or.(ibrseg(iprvn,ipl).Ne.ibrseg(irecn,ipl))) run=.False.
           If (iprvn .Eq. 0) Then
              If (curr_BCtp(ipl) == 2) Then !> flux
                 axialRootFlow(irecn-1,ipl) = Q_bc1(irecn-1,ipl)
                 If (Q_bc1(irecn-1,ipl) == 0) Then
                    axialRootFlow(irecn,ipl) = Q_bc1(irecn-1,ipl)
                 EndIf
              Endif
              run=.False.
           Endif
           !> get total flow; which is addition of all flows at each large branch connecting to the seedQuid des branches lat\E9rales ? 
           If (iprvn .Eq. 0) Then
               totalAxialFlow = totalAxialFlow + axialRootFlow(irecn,ipl)
           Endif
           !> get velocity for each segement -> needed in soluteRoot
           veloRoot(irecn,ipl) = axialRootFlow(irecn,ipl)/crossSectionSeg(irecn,ipl)
           !> definition for the next run of the loop
           ifoln=irecn
           irecn=iprvn
           !> from here, not an apex
           n_apex=.False.
        End Do !run
    End Do Branchloop


!> calculate SINK term
    totalAxialFlow = 0.0_dp
    betaw=0.0_dp
    Jintot=0.0_dp
    Joutr=0.0_dp
    Jcol=0.0_dp
    Voltot=0.0_dp
    Do irecn=0,nrec(ipl)
        Joutr(irecn,ipl) = Qi(irecn,ipl)*(PHs(irecn,ipl)-PHr_tot(irecn+1,ipl)) !> \param Joutr [L3/T] (soil is still from 0 to nrec while root is from 1 to nrec+1
        sinkR(irecn,ipl)=Joutr(irecn,ipl) !> \param sinkR root flow! [L3/T]
        Jintot=Jintot+Joutr(irecn,ipl)/(dxgrid*dygrid*dzgrid*nsub(irecn,ipl))
        If (t.Eq.t_begin+dt .And. iter_root .Gt. 1 .Or. t.Gt.t_begin+dt) Then
            checksub=0.0_dp
            Do isub=1,nsub(irecn,ipl)
                dVolSk=0.0_dp
                c_i=cube_i(irecn,isub,ipl)
                If(irecn.Eq.0) c_i=cube_i(1,isub,ipl)
                If(lMacro .And. l_elmMacro(c_i)) Then
                !> Uptake from voxels with mixed material is shifted to bulk soil
                    If(n_neigh(c_i).Gt.0) Then
                        sink_dummy = Joutr(irecn,ipl)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Real(n_neigh(c_i))
                        betac_dummy = beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Real(n_neigh(c_i))
                        Do i_n=1,6
                            If(MacroList(c_i,i_n).Gt.0) Then
                                i_dummy = MacroList(c_i,i_n)
                                Sink_cube(i_dummy) = Sink_cube(i_dummy) + sink_dummy
                                betac_cube(i_dummy) = betac_cube(i_dummy) +  betac_dummy
                            End If
                        End Do
                    Else
                        Sink_cube(c_i) = 0._dp
                        betac_cube(c_i) = 0._dp
                    End If
                Else
                    Sink_cube(c_i)=Sink_cube(c_i) + Joutr(irecn,ipl)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                    betac_cube(c_i)=betac_cube(c_i) + beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                Endif
                Do ic=1,8
                    corner_ic=loc_Q(irecn,ic,isub,ipl)
                    If(irecn.Eq.0) corner_ic=loc_Q(1,ic,isub,ipl)
                    dVolSk=Joutr(irecn,ipl)*w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*w_sub(irecn,isub,ipl)!> delta flux extracted
                    sink(corner_ic)=sink(corner_ic)+dVolSk/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)!sink term
                    betaw(corner_ic)=1._dp
                    Voltot=Voltot+w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*w_sub(irecn,isub,ipl)
                End Do
                   checksub=checksub+w_sub(irecn,isub,ipl)
            End Do !>end sub-segments
        Else !first time step
            Do isub=1,nsub(irecn,ipl)
                c_i=cube_i(irecn,isub,ipl)
                If(irecn.Eq.0) c_i=cube_i(1,isub,ipl)
                If(lMacro .And. l_elmMacro(c_i)) Then
                    !Uptake from voxels with mixed material is shifted to bulk soil
                    If(n_neigh(c_i).Gt.0) Then
                        sink_dummy = Joutr(irecn,ipl)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Real(n_neigh(c_i))
                        betac_dummy = beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Real(n_neigh(c_i))
                        Do i_n=1,6
                            If(MacroList(c_i,i_n).Gt.0) Then
                               i_dummy = MacroList(c_i,i_n)
                               Sink_cube(i_dummy) = Sink_cube(i_dummy) + sink_dummy
                               betac_cube(i_dummy) = betac_cube(i_dummy) +  betac_dummy
                            End If
                        End Do
                    Else
                        Sink_cube(c_i) = 0._dp
                        betac_cube(c_i) = 0._dp
                    Endif
                Else
                    Sink_cube(c_i)=Sink_cube(c_i) + Joutr(irecn,ipl)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                    betac_cube(c_i)=betac_cube(c_i) + beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                Endif
                Do ic=1,8
                      corner_ic=loc_Q(irecn,ic,isub,ipl)
                      If(irecn.Eq.0) corner_ic=loc_Q(1,ic,isub,ipl)
                      sink(corner_ic)=sink(corner_ic)+Joutr(irecn,ipl)*w_sub(irecn,isub,ipl)*w_dis(irecn,ic,isub,ipl) &
                           /sum_dis(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                      betaw(corner_ic)=1._dp ! betaw is defined for the Reset subroutine
                End Do
            End Do
        Endif
        Jcol=Jcol+Joutr(irecn,ipl)
    End Do ! nodes loop

    RootSk = Sum(sink_cube*dxgrid*dygrid*dzgrid)
    betac_cube2 = betac_cube/Sum(betac_cube)

    !> actual transpiration
    If (Jcol.Ge.0) Then
        Tact(ipl)=Jcol
    Else
        Tact(ipl)=0.0_dp
    End If
	
 End Subroutine Snk_Doussan
 !**************************************************************************************
 !> ### Estimate the reduction factor due to gap oocurence around roots
 !> function of the local soil pressure head ###
  Function Gap(h)
    Use TypeDef
    Use RootData, Only: lGap,g1,g2
    Implicit None

    Real(dp), Intent(in) ::h
    Real(dp):: Gap

    Gap=1.0_dp
    If (lGap) Then
       If ((h.Lt.g1).And.(h.Gt.g2)) Then
          Gap=(h-g2)/(g1-g2)
       Elseif (h.Le.g2) Then
          Gap=1.E-9_dp
       Endif
    Endif
    Return
  End Function Gap
 !**************************************************************************************
 !> ### Estimate the multiplicative reduction factor due to AQP activity
 !> function of the local soil pressure head ###
  Function AQPc(h,ipl)
    Use TypeDef
    Use RootData, Only: lAQPc,AQPh,AQPv,nAQPc
    Implicit None

    Integer(ap):: ih
    Integer(ap), Intent(in) :: ipl
    Real(dp), Intent(in) ::h
    Real(dp):: AQPc

    AQPc=1.0_dp
    If (lAQPc) Then
       ih=1
       Do While (AQPh(ih,ipl).Ge.h)
          ih=ih+1
       Enddo
       If (ih.Gt.nAQPc) Then
          AQPc=AQPv(nAQPc,ipl)
       Elseif (ih.Eq.1) Then
          AQPc=AQPv(1,ipl)
       Else
          AQPc=AQPv(ih-1,ipl)+(AQPv(ih,ipl)-AQPv(ih-1,ipl))*(h-AQPh(ih-1,ipl))/(AQPh(ih,ipl)-AQPh(ih-1,ipl))
       Endif
    Endif
    Return
  End Function AQPc

 !******************************************************************************
 !### calcualte SSF, Krs for Couvreur et al. (2012) model###
  Subroutine SetupCou(t)
    Use TypeDef
    Use GridData
    Use RootData
    Implicit None

    Integer(ap) :: pos,i,iElm
    Real(dp) :: t,min_time,diff_time

    pos=1
    min_time=Abs(t-timeobs(1))
    Do i=2,ntimeobs
       diff_time=Abs(t-timeobs(i))
       If (diff_time .Lt. min_time) Then
          pos=i
          min_time=diff_time
       End If
    End Do


    If (timeobs(pos) .Gt. t) Then
       Krs=(Krs_mat(pos)-Krs_mat(pos-1))/(timeobs(pos)-timeobs(pos-1))*(t-timeobs(pos-1))+Krs_mat(pos-1)
       Kcomp=(Kcomp_mat(pos)-Kcomp_mat(pos-1))/(timeobs(pos)-timeobs(pos-1))*(t-timeobs(pos-1))+Kcomp_mat(pos-1)
       Do iElm=1,nElm
          SSF(iElm)=(SSF_mat(iElm,pos)-SSF_mat(iElm,pos-1))/(timeobs(pos)-timeobs(pos-1))*(t-timeobs(pos-1))+SSF_mat(iElm,pos-1)
       End Do
    Elseif (timeobs(pos) .Eq. t) Then
       Krs=Krs_mat(pos)
       Kcomp=Kcomp_mat(pos)
       Do iElm=1,nElm
          SSF(iElm)=SSF_mat(iElm,pos)
       End Do
    Elseif ((timeobs(pos) .Lt. t) .And. (pos .Lt. Size(timeobs))) Then
       Krs=(Krs_mat(pos+1)-Krs_mat(pos))/(timeobs(pos+1)-timeobs(pos))*(t-timeobs(pos))+Krs_mat(pos)
       Kcomp=(Kcomp_mat(pos+1)-Kcomp_mat(pos))/(timeobs(pos+1)-timeobs(pos))*(t-timeobs(pos))+Kcomp_mat(pos)
       Do iElm=1,nElm
          SSF(iElm)=(SSF_mat(iElm,pos+1)-SSF_mat(iElm,pos))/(timeobs(pos+1)-timeobs(pos))*(t-timeobs(pos))+SSF_mat(iElm,pos)
       End Do
    Else
       Krs=Krs_mat(pos)
       Kcomp=Kcomp_mat(pos)
       Do iElm=1,nElm
          SSF(iElm)=SSF_mat(iElm,pos)
       End Do
    End If
  End Subroutine SetupCou
 !-----------------------------------------------------------------
 !> ### calculates Couvreur et al. (2012) RWU parameters
 !> solves DOUSSAN flow equations within the xylem ###
  Subroutine SSFdis(kout)
    Use TypeDef
    Use SparseMatrix
    Use NumericalRecipes
    Use GridData, Only: dxgrid,dygrid,dzgrid,zGrid,nElm,nPt,nex,ney,nez,SSF,Wn,nx,ny,nz
    Use DoussanMat, Only: GH,w_dis,nsub,loc_Q,Khr,w_sub,sum_dis,Lr,Lr_pot,hx_min,plantmatrix,cube_i,Inv_c1KxKr,PHs,stresfun
    Use RootData, Only: nrec,seglen,segsur,nbr,ibrseg,irecpr,Krs,Kcomp,ldJvL,lSinkCube,lSUF,nplant,vol_root,lDou
    Use TardieuData, Only: Krs0,Kcomp0
    Use SolData, Only: hNew
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap), Intent(in)::kout
    Integer(ap), Allocatable,Dimension(:) :: ind
    Integer(ap) :: iter_loc,itol,itmaxl,corner_ic,ibr,c_i,i,j,k,nk
    Integer(ap) :: iprvn,irecn,ifoln,ic,isub,err,ipl
    Real(dp) :: PHr_tot1(1:nrec(nplant)+1),PHr_tot2(1:nrec(nplant)+1),err_loc,Q1(1:nrec(nplant)+1),Q2(1:nrec(nplant)+1),Joutr1(1:nrec(nplant)),Joutr2(1:nrec(nplant))
    Real(dp) :: PHs1(0:nrec(nplant)),PHs2(0:nrec(nplant)),maxerr=1.e-10_dp,Jcol1,Jcol2,SSFt,SSFcum,sumPhiHsi,sumPhi,sumPhi2,sumHsi,sumHsi2,Beta(1:nx*ny*nz)
    Real(dp) :: Hsi,r,Phi,PhiSUF(1:nrec(nplant)),Hseq,Lr_fact(1:nrec(nplant))
    Real(dp), Allocatable,Dimension(:,:) ::SSFunstockNew,SSFunstock
    Real(dp), Allocatable,Dimension(:) ::sink_cube1,sink_cube2,sink1,sink2,Hs_cube2
    Logical :: n_apex,run
    Character :: file*17
    

    Write (*,*) '... Calculating SSF and Krs ...'
    If (.Not.lSUF.Or.kout.Eq.0) Then
       If (lDou) lSUF=.TRUE.
       Write (*,*) '... Solving Doussan system ...'
       Allocate (sink_cube1(1:nElm))
       Allocate (sink_cube2(1:nElm))
       Allocate (Hs_cube2(1:nElm))
       sink_cube1=0._dp
       sink_cube2=0._dp
       Hs_cube2=0._dp
       Allocate (sink1(1:nPt))
       Allocate (sink2(1:nPt))
       sink1=0._dp
       sink2=0._dp
       Do ipl=1,nplant
          !> PHs homogeneous when characterizing SSF and Krs
          PHs1(0:nrec(ipl))=-300.0_dp
          !> PHs heterogeneous when characterizing Kcomp
          PHs2(0:nrec(ipl))=-150.0_dp+2*GH(0:nrec(ipl),ipl)
       End Do
       Do ipl=1,nplant      
          err=SM_allocate(plantmatrix(ipl), nrec(ipl)+1, nrec(ipl)+1)
          If(err/=0) Call stop_program('Could not create plantmatrix')
          Do ibr=1,nbr(ipl) !all plants have same number of roots!
             run=.True.
             n_apex=.False.
             !* find the tip segment of the branch 'ibr'
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.ibr)
                irecn=irecn-1
             End Do
             If (seglen(irecn,ipl)<1.E-20) Then ! skip this segment too small to be taken
                ifoln=irecn ! following node ID
                irecn=irecpr(irecn,ipl) !current node ID
             Else
                n_apex=.True.!ifoln does not exist if not continu
             Endif

             !* then the rest of the branch up to the seed of the embranchment
             Do While (run)
                !* "upper" node
                iprvn=irecpr(irecn,ipl)
                Q1(irecn+1)=Lr(irecn,ipl)*segsur(irecn,ipl)*PHs1(irecn)
                Q2(irecn+1)=Lr(irecn,ipl)*segsur(irecn,ipl)*PHs2(irecn)
                err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                !if apex (bottom part of root)
                If (n_apex) Then
                   Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn,ipl)/seglen(irecn,ipl)+Lr(irecn,ipl)*segsur(irecn,ipl),.False.)!Here Lr corresponds to Lr_pot if Solveroot was not called yet
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Else
                   Err=SM_set(Plantmatrix(ipl), irecn+1,ifoln+1,-Khr(ifoln,ipl)/seglen(ifoln,ipl),.False.) !row, col,value
                   Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn,ipl)/seglen(irecn,ipl)+Khr(ifoln,ipl)/seglen(ifoln,ipl)+Lr(irecn,ipl)*segsur(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Endif
                If (iprvn==0) Then!seed=first segment
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                   Q1(iprvn+1)=-100.0_dp
                   Q2(iprvn+1)=-100.0_dp
                   run=.False.
                Elseif (ibrseg(iprvn,ipl).Ne.ibrseg(irecn,ipl)) Then!start of the branch but not from the seed
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn,ipl)/seglen(irecn,ipl),.False.)!iprvn+1,iprvn+1
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)!iprvn+1,irecn+1
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                   run=.False.
                Endif
                ! definition for the next run of the loop
                ifoln=irecn
                irecn=iprvn
                n_apex=.False.
             End Do ! loop on branch nodes
          End Do ! loop on branches

          PHr_tot1(1:nrec(ipl)+1)=PHs1(0:nrec(ipl))   ! initial xylem PH -> could also be zero but then convergence should take longer (this should be the common case)
          PHr_tot2(1:nrec(ipl)+1)=PHs2(0:nrec(ipl))
          itol=3
          itmaxl=(nrec(ipl)+1)*3 !maximum iterations for biconjugate gradient-> needed in case of root growth
          ! solve with biconjugated gradient method (Numerical recipies)
          Call linbcg(nrec(ipl)+1,Q1(1:nrec(ipl)+1),PHr_tot1,itol,maxerr,itmaxl,iter_loc,err_loc,ipl)           
          If (err_loc.Gt.maxerr) Then
             Print *,'No convergence in the root system, err_loc > maxerr',err_loc,maxerr
             Call stop_program('')
          Endif
          Call linbcg(nrec(ipl)+1,Q2(1:nrec(ipl)+1),PHr_tot2,itol,maxerr,itmaxl,iter_loc,err_loc,ipl)           
          If (err_loc.Gt.maxerr) Then
             Print *,'No convergence in the root system, err_loc > maxerr',err_loc,maxerr
             Call stop_program('')
          Endif

          !====================== end solve root system ====================================
          !***********************************************************************
          ! calculate SINK term
          !**********************************************************************
          ! initialisation
          Jcol1=0.0_dp
          Jcol2=0.0_dp
          Do irecn=1,nrec(ipl)
             Joutr1(irecn) = Lr(irecn,ipl)*segsur(irecn,ipl)*(PHs1(irecn)-PHr_tot1(irecn+1)) ![L3/T] (soil is still from 0 to nrec while root is from 1 to nrec+1)
             Joutr2(irecn) = Lr(irecn,ipl)*segsur(irecn,ipl)*(PHs2(irecn)-PHr_tot2(irecn+1))
             If (.Not.lSUF) Then
                Do isub=1,nsub(irecn,ipl)
                   c_i=cube_i(irecn,isub,ipl)
                   Sink_cube1(c_i)=Sink_cube1(c_i)+Joutr1(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                   Sink_cube2(c_i)=Sink_cube2(c_i)+Joutr2(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                   If (.Not.lSinkCube) Then
                      Do ic=1,8
                         corner_ic=loc_Q(irecn,ic,isub,ipl)
                         sink1(corner_ic)=sink1(corner_ic)+Joutr1(irecn)*w_sub(irecn,isub,ipl)*w_dis(irecn,ic,isub,ipl) &
                              /sum_dis(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                         sink2(corner_ic)=sink2(corner_ic)+Joutr2(irecn)*w_sub(irecn,isub,ipl)*w_dis(irecn,ic,isub,ipl) &
                              /sum_dis(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                      End Do
                   Endif
                End Do
             Endif
             Jcol1=Jcol1+Joutr1(irecn)
             Jcol2=Jcol2+Joutr2(irecn)
          End Do ! nodes loop

          Krs=Jcol1/(-300.0_dp-PHr_tot1(1))
          If (.Not.lSUF) Then
             If(.Not.Allocated(SSF)) Allocate (SSF(1:nElm))
             SSF(1:nElm)=Sink_cube1(1:nElm)*dxgrid*dygrid*dzgrid/Jcol1
             If (.Not.lSinkCube) Then
                Beta(1:nx*ny*nz)=sink1(1:nx*ny*nz)/Jcol1
             Endif
          Else!> if lSUF (and kout=0)
             Allocate (SSF(1:nrec(ipl)))
             Allocate (Inv_c1KxKr(1:nrec(ipl)))
             SSF(1:nrec(ipl))=Joutr1(1:nrec(ipl))/Jcol1
             Inv_c1KxKr(1:nrec(ipl))=SSF(1:nrec(ipl))*Krs
          Endif
       End Do !plant loop
    Else!> if lSUF and kout>0
       Do ipl=1,nplant !Actually currently not made to work with several plants
          Write (*,*) '... Using Inv_c1 ...'
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
          SSFcum=Sum(SSF(1:nrec(ipl)))
       End Do
    Endif
    Krs0=Krs   
    Print *,'... Calculating Kcomp ...'
    Do ipl=1,nplant 
       If (.Not.lSUF.Or.kout.Eq.0) Then
          If (.Not.lSUF) Then
             Allocate (SSFunstock(1:nElm,1:2))
             SSFunstock(1:nElm,2)=SSF!> Will progressively unstock SSFi
             Do i=1,nElm
                SSFunstock(i,1)=i!> Contains element numbers
                Hs_cube2(i)=-150.0_dp+(zGrid(i)+zGrid(i+nx*ny)) !> Works for both XY periodic and non-periodic domains
             End Do
             nk=nElm!> \param nk=number of SSFi stocked in SSFunstock
             SSFt=0.0001_dp
             !> \param SSFt Threshold value above which SSFi is considered in the estimation of Kcomp
             SSFcum=0.0_dp!> \param SSFcum Cumulation of "accepted" SSFi
             j=0!> j=number of SSFi already cumulated in SSFcum
             Do While (SSFcum.Lt.0.9_dp)
                SSFt=0.5_dp*SSFt+0.5_dp*Sum(SSFunstock(1:nk,2))/nk
                Allocate (SSFunstockNew(1:nk,1:2))
                SSFunstockNew=0.0_dp
                k=0!> k=number of SSFi stocked in SSFunstockNew
                Do i=1,nk
                   If (SSFunstock(i,2).Ge.SSFt) Then
                      SSFcum=SSFcum+SSFunstock(i,2)
                      j=j+1
                      Phi=Sink_cube2(Int(SSFunstock(i,1)))*dxgrid*dygrid*dzgrid/SSFunstock(i,2)-Jcol2
                      Hsi=Hs_cube2(Int(SSFunstock(i,1)))+200.0_dp
                      sumPhiHsi=sumPhiHsi+Phi*Hsi
                      sumHsi=sumHsi+Hsi
                      sumHsi2=sumHsi2+Hsi**2
                      sumPhi=sumPhi+Phi
                      sumPhi2=sumPhi2+Phi**2
                   Elseif (SSF(i).Gt.0.0_dp) Then
                      SSFunstockNew(k+1,1)=SSFunstock(i,1)
                      SSFunstockNew(k+1,2)=SSFunstock(i,2)
                      k=k+1
                   Endif
                End Do
                nk=k
                Deallocate (SSFunstock)
                Allocate (SSFunstock(1:nk,1:2))
                SSFunstock=SSFunstockNew(1:nk,1:2)
                Deallocate (SSFunstockNew)
             End Do
             Kcomp=(j*sumPhiHsi-sumHsi*sumPhi)/(j*sumHsi2-sumHsi**2)
             r=(sumPhiHsi-sumHsi*sumPhi/j)/Sqrt((sumHsi2-sumHsi**2/j)*(sumPhi2-sumPhi**2/j))
          Else!if lSUF
             Allocate (ind(1:nrec(ipl)))
             Hseq=DOT_Product(SSF(1:nrec(ipl)),PHs2(1:nrec(ipl)))
             k=0
             Do i=1,nrec(ipl)
                If (SSF(i).Gt.0.5/nrec(ipl)) Then
                   PhiSUF(i)=Joutr2(i)/SSF(i)-Jcol2
                   k=k+1
                   ind(k)=i
                Else
                   PhiSUF(i)=0.0
                Endif
             End Do
             Kcomp=DOT_Product(PHs2(ind(1:k))-Hseq,PhiSUF(ind(1:k)))/DOT_Product(PHs2(ind(1:k))-Hseq,PHs2(ind(1:k))-Hseq)
             r=(DOT_Product(PhiSUF(ind(1:k)),PHs2(ind(1:k)))-Sum(PHs2(ind(1:k)))*Sum(PhiSUF(ind(1:k)))/k)/Sqrt((DOT_Product(PhiSUF(ind(1:k)),PhiSUF(ind(1:k)))-Sum(PhiSUF(ind(1:k)))**2/k)*(DOT_Product(PHs2(ind(1:k)),PHs2(ind(1:k)))-Sum(PHs2(ind(1:k)))**2/k))
             SSFcum=Sum(SSF(ind(1:k)))
          Endif
          Write (*,*) 'Krs [cm\B3/d]',Krs
          Write (*,*) 'Kcomp [cm\B3/d]',Kcomp
          Write (*,*) 'Correlation coefficient [-]',r
          Write (*,*) 'Percentage of cumulated SSF considered in the calculation of Kcomp [%]',SSFcum*100
       Else!if lSUF and kout>0
          PHs(1:nrec(ipl),ipl)=0._dp
          Do irecn=1,nrec(ipl)
             Do isub=1,nsub(irecn,ipl)
                PHs(irecn,ipl)=PHs(irecn,ipl)+Sum(hnew(loc_Q(irecn,1:8,isub,ipl))*w_dis(irecn,1:8,isub,ipl)/Sum(w_dis(irecn,1:8,isub,ipl)))*w_sub(irecn,isub,ipl)/Sum(w_sub(irecn,1:nsub(irecn,ipl),ipl))
             End Do    !end do-loop over subsegments
             Lr_fact(irecn)=Gap(PHs(irecn,ipl))*AQPc(PHs(irecn,ipl),ipl)
             If (.NOT. lDou) Then
               Lr(irecn,ipl)=Lr_pot(irecn)*Lr_fact(irecn)
             Endif
          End Do
          Krs=DOT_Product(Inv_c1KxKr(1:nrec(ipl)),Lr_fact(1:nrec(ipl)))
          SSF(1:nrec(ipl))=Inv_c1KxKr(1:nrec(ipl))*Lr_fact(1:nrec(ipl))/Krs
          Kcomp=Krs
       Endif
       Kcomp0=Kcomp
       !> writes Couvreur.out
       Write (file,'(A12)')'out/Couvreur'
       Write (file(13:13),'(I1)') ipl
       Write (file(14:14),'(A1)') '.'
       If (kOut.Lt.10) Then
          Write (file(15:15),'(I1)') kOut
       Elseif (kOut.Lt.100) Then
          Write (file(15:16),'(I2)') kOut
       Else
          Write (file(15:17),'(I3)') kOut
       Endif
       Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
       Write (8,'(''********** COUVREUR ET AL. (2012) ROOT WATER UPTAKE MODEL PARAMETERS **********'')')
       Write (8,*)
       Write (8,'(''Stress function ? (1=yes,2=no)/ Threshold value of the stress function [hPa]'')')
       If (stresfun.Eq.1) Then
          Write (8,'(I1,3X,F10.2)') 1,hx_min
       Else
          Write (8,'(I1,3X,I1)') 2,0!> non isohydric plants (Tuzet, etc. cannot currently be considered when using Couvreur RWU model)
       Endif
       Write (8,*)
       Write (8,'(''Dimension of the root water uptake model (1=soil domain,2=root system domain)'')')
       If (lSUF) Then
          Write (8,'(I1)') 2
       Else
          Write (8,'(I1)') 1
       Endif
       Write (8,*)
       Write (8,'(''Use of modified de Jong van Lier function for Hseq prediction ? (1=yes,2=no)'')')
       If (ldJvL) Then
          Write (8,'(I1)') 1
       Else
          Write (8,'(I1)') 2
       Endif
       Write (8,*)
       Write (8,'(''Root volume [cm³]'')')
       Write (8,'(6pE20.12)') vol_root
       Write (8,*)
       Write (8,'(''Equivalent conductance of the whole root system [cm\B3/hPa/day]'')')
       Write (8,'(F10.8)') Krs
       Write (8,*)
       Write (8,'(''Compensatory RWU conductance of the root system [cm\B3/hPa/day]'')')
       Write (8,'(F10.8)') Kcomp
       Write (8,*)
       If (lSUF) Then
          Write (8,'(''* Standard Uptake Fraction distribution'')')
          Write (8,'(''Number of root segments [-], root volume [cm³]'')')
          Write (8,'(I6,3X,F12.10)') nrec(ipl),vol_root
          Write (8,'(''  SegID#    SUF'')')
          Do i=1,nrec(ipl)
             Write (8,'(I6,3X,F12.10)') i,SSF(i)!> Note that both SUF and SSF have the same name (SSF) in the script since they cannot be used in the same time
          Enddo
          Close (8)
       Elseif (lSinkCube) Then
          Write (8,'(''* Standard Sink Fraction distribution'')')
          Write (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
          Write (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
          Write (8,'(''  Element#    SSF'')')
          Do i=1,nElm
             Write (8,'(I6,3X,F12.10)') i,SSF(i)
          Enddo
          Close (8)
       Else
          Write (8,'(''* Standard Sink Fraction distribution'')')
          Write (8,'(''Number of nodes in each dimension [-] and elements size [cm]'')')
          Write (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nx,ny,nz,dxGrid,dyGrid,dzGrid
          Write (8,'(''  Node#    BetaSSF'')')
          Do i=1,nx*ny*nz
             Write (8,'(I6,3X,F14.10)') i,Beta(nx*ny*nz+1-i)!> Nodes # inverted for Hydrus
          Enddo
          Close (8)
          Call stop_program('SSF generated for Hydrus!')
       Endif
    Enddo ! plant loop
    Return
  End Subroutine SSFdis
  !***************************************************************************
  Subroutine analytical_approach(ipl)
    !      USE OMP_LIB
    Use TypeDef
    Use Paramdata,Only: MaxNod, pi,lSalinity
    Use GridData, Only: xgrid,ygrid,zgrid,dxgrid,dygrid,dzgrid
    Use DoussanMat, Only :old,oldT,ave,eqDis,nsub,PHs,loc_q,cp_mean,voxel_no,no_voxels, &
         voxel_node,numnodes_voxel,h_mat2,hx_min,phr_sub,n,indexvalue,counter2,phi_mat, &
         PH_micro2,ph_root_sub,Qd,w_dis,w_sub,cent,Intc,PHr,Lr,sum_dis,switchcriterion,Qi, &
         tcheck,tcheck2,count_nodes,k_ave,B,Q_bc1,Q_bc2,cube_i,PHs_osmotic,PHo
    Use RootData, Only: nrec,seglen,segsur,sigma
    Use SolData, Only: hnew, matnum
    Use RootGrowthNeighb, Only: FlowNeighb
    Implicit None

    Integer(ap) :: checker1=0,blaat,irecnTemp,isubTemp,ki
    Integer(ap) ::j,counter,l,w,Flow_corner(1:8,1:6)
    Integer(ap) ::corner(1:8),var1,kl,imin,ipl
    Integer(ap) :: irecn,i,isub,pos,ja=0,irecn_tmp,isub_tmp
    Real(dp) :: beta
    Real(dp) :: y(1:n),h_mat(1:n),arr(1:n-1),Phi(1:n-1)
    Real(dp) :: K_k,C_k,theta_k,temp
    Real(dp) :: cumsum(1:n-1),Phi_rend
    Real(dp) :: z1,z2,frac1,frac2,xCent,yCent,Zcent,xInt,yInt,zInt
    Real(dp) :: testx(2),testy(2),testz(2),xgrid_max,xgrid_min,ygrid_max,ygrid_min,mXY(2),point(1:4,1:2)
    Real(dp) :: zgrid_max,zgrid_min,mZX(2),mZY(2)
    Real(dp) :: PH_micro(1:n-1)
    Real(dp) :: theta(1:n),K(1:n),C(1:n),rho,r0,r(1:n-1),rend
    Real(dp) ::sumW,psub,q_root,q_out
    Real(dp) :: pos_cent,testvar
    Real(dp),Dimension(1:4,1:4) :: dis,Flowpoint_rel,PHpoint_rel
    Real(dp),Dimension(1:4) :: mFlows,mPHs,PH_inter,Flow_inter
    Real(dp),Dimension(1:8) ::h_matK,temp2
    Real(dp),Dimension(1:8) ::pond_i,q_tot_c
    Real(dp) :: meanPHS,meanFlows
    Real(dp) :: z_angle,x_angle,y_angle,x1,x2,y1,y2,q_root_tmp(1:indexValue),pond_subTemp
    Logical :: zalign,yalign,xalign,extra_checker=.False.,countja

    !> parameters for counting nodes that are under limiting conditions (ks < lr*r0)
    countja=.False.
    count_nodes=0
    Do irecn=0,nrec(ipl) !> loop over segments
       old=oldT
       x_angle=0
       y_angle=0
       z_angle=0
       !> initialize soil PH matrix
       PHs(irecn,ipl)=0._dp
       Do isub=1,nsub(irecn,ipl) !> loop over subsegments
          dis=0
          Flowpoint_rel=0
          PHpoint_rel=0
          mFlows=0
          mPHs=0
          PH_inter=0
          Flow_inter=0
          !> ---------------------------calculate the flux in each soil node with Darcy\B4s law--------
          If (.Not.eqDis) Then
             Do i=1,8 !> for each corner of a cube
                corner(i)=loc_Q(irecn,i,isub,ipl)
             End Do
             !> determine nodes surrounding a soil node to determine the flow
             Call FlowNeighb(corner,Flow_corner)
             Call flux_node(Flow_corner,corner,q_tot_c)
          Elseif (eqDis) Then !method C
             q_tot_c(1:8)=0
          Endif
          !> ---------------------------end flow calculation-----------------------------------
          !> initialization and r0 definement
          checker1=0 !> reset for each cp
          Do i=1,8
             pond_i(i)=w_dis(irecn,i,isub,ipl)
             corner(i)=loc_Q(irecn,i,isub,ipl)
          End Do
          zalign=.False.
          yalign=.False.
          xalign=.False.
          !> not always true: adjust accordingly
          If (irecn .Eq. 0) Then
             r0=5e-2 !set to minimal radius for root collar ID only-> not important
          Else
             If (segsur(irecn,ipl) .Eq. 0 .Or. seglen(irecn,ipl) .Eq. 0) Then
                r0 = 5e-2 !> minimal root radius (in cm)
             Else
                r0 = 1/(2*pi)*segsur(irecn,ipl)/seglen(irecn,ipl)
             End If
          Endif
          !-----------------------------------------------------------------------
          !
          !> REDUCTION from 3D -> 2D plane, currently averaging procedure is taken
          !> in axial direction (i.e. direction perpendicular to radial flow)  which is sufficient (checked)
          !----------------------------------------------------------------------
          !> for the one_root procedure; the center point of a segment is the middle of a soil voxel
          If (ave) Then
             xCent = cp_mean(irecn,1,isub,ipl)
             yCent = cp_mean(irecn,2,isub,ipl)
             zCent = cp_mean(irecn,3,isub,ipl)
             imin = voxel_no(irecn,isub) !voxel number
             If (counter2(imin) .Ne. 0) Then !> check of conductivity drop in this voxel is already calculated
                irecnTemp = voxel_node(imin,1,counter2(imin))
                isubTemp = voxel_node(imin,2,counter2(imin))
                pond_subTemp = w_sub(irecnTemp,isubTemp,ipl)
                extra_checker=.True.
                Goto 210 !> to calculate analytical solution only once for the large root node
             Elseif (counter2(imin) .Eq. 0 .And. no_voxels(imin) .Ne. 0) Then
                Do blaat=1,no_voxels(imin)
                   If (voxel_node(imin,1,blaat) .Eq. irecn .And. voxel_node(imin,2,blaat) &
                        .Eq. isub .And. w_sub(irecn,isub,ipl).Ne.0) counter2(imin)=blaat
                Enddo
             Endif
          Elseif (.Not.(ave)) Then !> if not one root method each root node is a center point
             xCent = cent(irecn,1,isub,ipl)
             yCent = cent(irecn,2,isub,ipl)
             zCent = cent(irecn,3,isub,ipl)
          Endif
          !> interception point x,y,zInt with soil cube
          xInt = Intc(irecn,1,isub,ipl)
          yInt = Intc(irecn,2,isub,ipl)
          zInt = Intc(irecn,3,isub,ipl)
          !> ---------------------------derive x,y,z-allignment and segment angle------------------
          !> Get aligment of a root in the current voxel. If an intersection point of a root with the cube is 
          !> identical to the center point of the root, pick a direction dependent on the surface it crosses
          !> otherwise calculate the segment angle. first from a z-directional point of view, if the branching 
          !> angle is smaller then a given range then look if it is
          !> aligned in x-direction (using a projection of the z-coordinate on the 2D generated plain)
          !> if this angle doesnt comply to a criterion then the allignement is performed in y-direction
          !> weigh factor Wb is calculated but not needed in remainder calculations
          !> dimensions of the voxel, maxima and minima in each direction
          zgrid_max = Maxval(zgrid(corner(1:8)))
          zgrid_min = Minval(zgrid(corner(1:8)))
          xgrid_max = Maxval(xgrid(corner(1:8)))
          xgrid_min = Minval(xgrid(corner(1:8)))
          ygrid_max = Maxval(ygrid(corner(1:8)))
          ygrid_min = Minval(ygrid(corner(1:8)))
          !> get alignment direction
          If (xCent-xInt .Eq. 0 .And. yCent-yInt .Eq. 0 .And. zCent-zInt .Eq. 0) Then
             If (zCent .Eq. zgrid_max .Or. zCent .Eq. zgrid_min) Then
                zalign=.True.
             Elseif (yCent .Eq. ygrid_max .Or. yCent .Eq. ygrid_min) Then
                yalign=.True.
             Elseif (xCent .Eq. xgrid_max .Or. xCent .Eq. xgrid_min) Then
                xalign=.True.
             Endif
          Else
             !> geometry; beta is radial angle betwen z difference and the length of the segment-> convert to degree
             beta = Asin(Abs(zInt-zCent)/Sqrt((xInt-xCent)*(xInt-xCent)+ &
                  (yInt-yCent)*(yInt-yCent)+(zInt-zCent)*(zInt-zCent)))
             z_angle=Abs(beta/(2*pi))*360
             If (z_angle .Le. 90.1 .And. z_angle .Ge. 45) Then !> certain error deviatin given to criterion
                zalign=.True.
             Else
                zInt=zCent !> projection in 2d plane to get angle in x-or y direction
                beta = Acos(Abs(xInt-xCent)/Sqrt((xInt-xCent)*(xInt-xCent)+ &
                     (yInt-yCent)*(yInt-yCent)+(zInt-zCent)*(zInt-zCent)))
                x_angle=Abs(beta/(2*pi))*360
                x_angle=90-x_angle !> becoz of cosine, angle interpreted from 0 to 45 and angle has to be adjusted for weight factor
                If (x_angle .Le. 90 .And. x_angle .Ge. 45) Then
                   xalign=.True.
                Else
                   yalign=.True.
                Endif
             Endif
          Endif
          !> -----------------------------------end allignment-----------------------------------
          !> for each allocated root direction in a voxel a checker1 is set which determines where the cp is located, 
          !> in the voxel or at the edges. If it is located in the  ! voxel then the analytical solution is applied, 
          !> otherwise the old method is sufficient
          !> set checker1
          If (zalign) Then
             If (Abs(xCent-xgrid_max) .Gt. r0 .And. Abs(xCent-xgrid_min).Gt.r0 &
                  .And. Abs(yCent-ygrid_max) .Gt. r0 .And. Abs(yCent-ygrid_min).Gt.r0 ) Then !in plane
                checker1 = 1
             Else !> if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
                checker1=0
             End If
          Elseif (xalign) Then
             If (Abs(yCent-ygrid_max) .Gt. r0 .And. Abs(yCent-ygrid_min).Gt. r0 &
                  .And. Abs(zCent-zgrid_max) .Gt. r0 .And. Abs(zCent-zgrid_min).Gt.r0 ) Then !in plane
                checker1 = 1
             Else !> if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
                checker1=0
             End If
          Elseif (yalign) Then
             If (Abs(xCent-xgrid_max) .Gt. r0 .And. Abs(xCent-xgrid_min).Gt.r0 &
                  .And. Abs(zCent-zgrid_max) .Gt. r0 .And. Abs(zCent-zgrid_min).Gt.r0 ) Then !in plane
                checker1 = 1
             Else !> if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
                checker1=0
             End If
          Endif
          !> Average for the specific direction chosen from a 3D to a 2D mesh (for application of the 2D analytical solution). 
          !> In direction plane is orientated
          !> To be clearer; using the allignment of the root node the PH and Flow boundary conditions
          !> are transferred to the 4 corner nodes spanning the 2D domain in which the
          !> analytical approach is applied
          !>
          !>              C1----------------C2
          !>|  |
          !>|  |
          !>|  |
          !>|o  |
          !>|  |
          !>|  |
          !>|  |
          !>C3----------------C4
          !>
          !> z-align, x-align, y-align are discussed; basicly the same routines;
          If (checker1.Eq.1 .And. (zalign)) Then
             Do i=1,4
                z1 = zgrid(corner(i))
                z2 = zgrid(corner(i+4))
                frac1 = 1/Abs(zCent-z1)
                frac2 = 1/Abs(zCent-z2)
                If (Abs(zCent-z1) .Eq. 0) Then
                   PH_inter(i) = hnew(corner(i))
                   !> ph bound. cond.
                   Flow_inter(i) = q_tot_c(i)
                   !> flow bound. cond.
                Elseif (Abs(zCent-z2) .Eq. 0) Then
                   PH_inter(i) = hnew(corner(i+4))
                   Flow_inter(i) = q_tot_c(i+4)
                Else
                   PH_inter(i) = (frac1*hnew(corner(i)) + frac2*hnew(corner(i+4)))/(frac1+frac2)
                   Flow_inter(i) = (frac1*q_tot_c(i) + frac2*q_tot_c(i+4))/(frac1+frac2)
                   !> corner -> inter (2D plane)
                Endif
             Enddo
             testx(1) = Abs(xgrid_max-xCent)
             testx(2) = Abs(xgrid_min-xCent)
             mXY(1) = Minval(testx) !> minimum distance in x-dimension
             testy(1) = Abs(ygrid_max-yCent)
             testy(2) = Abs(ygrid_min-yCent)
             mXY(2) = Minval(testy) !> minimum distance in y-direction
             !> get minimum radius from root node to a 2D edge; to be used as outer radius
             !> for the analytical approach
             rend = Minval(mXY)     !> outer radius of the soil cylinder
             If (eqDis) Then !> get an equal radius dependent on the number of root nodes in this voxel
                testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
                rend = 1/2d0*Sqrt((dzgrid)/Floor(testvar))
             Endif
             If (rend .Lt. r0) Then
                ja=1
                Goto 203
             Endif
             !> the dashed arrow line is set as the minimal distance in this example; i.e. rend (outer radius)
             !>             C1------0---------C2
             !>||  |
             !>||  |
             !>||  |
             !>0<----->o---------0
             !>||  |
             !>||  |
             !>||  |
             !>C3------0---------C4
             !> new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
             !> so the outer radius has intersection points with the 2D plane as drawn above (the zeros\B4s (0) are the
             !> intersection points with the 2D plane
             !> these intersection points are calculated here; for method B,C (see MS)
             If (eqDis) Then !> middle of voxel as point of departure for bnd.cond.
                point(1,1) = (xgrid_max+xgrid_min)/2+rend
                point(1,2) = (ygrid_max+ygrid_min)/2
                point(2,1) = (xgrid_max+xgrid_min)/2
                point(2,2) = (ygrid_max+ygrid_min)/2+rend
                point(3,1) = (xgrid_max+xgrid_min)/2-rend
                point(3,2) = (ygrid_max+ygrid_min)/2
                point(4,1) = (xgrid_max+xgrid_min)/2
                point(4,2) = (ygrid_max+ygrid_min)/2-rend
             Else
                point(1,1) = xCent+rend
                point(1,2) = yCent
                point(2,1) = xCent
                point(2,2) = yCent+rend
                point(3,1) = xCent-rend
                point(3,2) = yCent
                point(4,1) = xCent
                point(4,2) = yCent-rend
             Endif
             !> calculate the pressure head, based on a distance average, for each new coordinate surrounding the 
             !> CP(center ponint, i.e. the root node)---- 2D!!!!!!!! z-alignment-> x and y-direction is evaluated
             !> so outer boundary conditions are mapped on the 4 points that represent the outer circle
             !> The new PH in each intersection point (zero) with 2D domain is mPHs, here only part of that solution is given
             !> due to distance saciling!! so "_rel" is intermediate value for PH or Flow (later these values are put to 
             !> the correct valuesin mPHs,mFlow
             Do j=1,Size(point,1)
                Do i=1,4
                   dis(j,i) = Sqrt((point(j,1)-xgrid(corner(i)))**2+(point(j,2)-ygrid(corner(i)))**2)
                   PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i) !> inter values being the values in the 4 corner nodes of the 2D domain
                   Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                   !> inter -> point (distance)
                End Do
             End Do
          Elseif (checker1.Eq.1 .And. (xalign)) Then
             kl=1
             Do i=1,7,2
                x1 = xgrid(corner(i))
                x2 = xgrid(corner(i+1))
                frac1 = 1/Abs(xCent-x1)
                frac2 = 1/Abs(xCent-x2)
                If (Abs(xCent-x1) .Eq. 0) Then
                   PH_inter(kl) = hnew(corner(i))
                   Flow_inter(kl) = q_tot_c(i)
                Else If (Abs(xCent-x2) .Eq. 0) Then
                   PH_inter(kl) = hnew(corner(i+1))
                   Flow_inter(kl) = q_tot_c(i+1)
                Else
                   PH_inter(kl) = (frac1*hnew(corner(i)) + frac2*hnew(corner(i+1)))/(frac1+frac2)
                   Flow_inter(kl) = (frac1*q_tot_c(i) + frac2*q_tot_c(i+1))/(frac1+frac2)
                   !> corner -> inter (2D plane)
                Endif
                kl=kl+1
             End Do
             testz(1) = Abs(zgrid_max-zCent)
             testz(2) = Abs(zgrid_min-zCent)
             mZY(1) = Minval(testz) !> minimum distance in z-dimension
             testy(1) = Abs(ygrid_max-yCent)
             testy(2) = Abs(ygrid_min-yCent)
             mZY(2) = Minval(testy) !> minimum distance in y-direction
             rend = Minval(mZY)     !> outer radius of the soil cylinder
             If (eqDis) Then
                testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
                rend = 1/2d0*Sqrt((dxgrid)/Floor(testvar))
             Endif
             If (rend .Lt. r0) Then
                ja=1
                Goto 203
             Endif
             !> new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
             If (eqDis) Then !> middle of voxel as point of departure for bnd.cond.
                point(1,1) = (zgrid_max+zgrid_min)/2+rend
                point(1,2) = (ygrid_max+ygrid_min)/2
                point(2,1) = (zgrid_max+zgrid_min)/2
                point(2,2) = (ygrid_max+ygrid_min)/2+rend
                point(3,1) = (zgrid_max+zgrid_min)/2-rend
                point(3,2) = (ygrid_max+ygrid_min)/2
                point(4,1) = (zgrid_max+zgrid_min)/2
                point(4,2) = (ygrid_max+ygrid_min)/2-rend
             Else
                point(1,1) = zCent+rend
                point(1,2) = yCent
                point(2,1) = zCent
                point(2,2) = yCent+rend
                point(3,1) = zCent-rend
                point(3,2) = yCent
                point(4,1) = zCent
                point(4,2) = yCent-rend
             Endif
             Do j=1,Size(point,1)
                Do i=1,4
                   If (i.Eq.1) kl=1 !> y and z-coord of these nodes needed
                   If (i.Eq.2) kl=3
                   If (i.Eq.3) kl=5
                   If (i.Eq.4) kl=7
                   dis(j,i) = Sqrt((point(j,1)-zgrid(corner(kl)))**2+(point(j,2)-ygrid(corner(kl)))**2)
                   PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i)
                   Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                   !> inter -> point (distance)
                End Do
             End Do
          Elseif (checker1.Eq.1 .And. (yalign)) Then
             kl=1
             Do i=1,4
                var1=i
                If (i.Gt.2) var1=var1+2
                y1 = xgrid(corner(var1))
                y2 = xgrid(corner(var1+2))
                frac1 = 1/Abs(yCent-y1)
                frac2 = 1/Abs(yCent-y2)
                If (Abs(yCent-y1) .Eq. 0) Then
                   PH_inter(kl) = hnew(corner(var1))
                   Flow_inter(kl) = q_tot_c(var1)
                Elseif (Abs(yCent-y2) .Eq. 0) Then
                   PH_inter(kl) = hnew(corner(var1+2))
                   Flow_inter(kl) = q_tot_c(var1+2)
                Else
                   PH_inter(kl) = (frac1*hnew(corner(var1)) + frac2*hnew(corner(var1+2)))/(frac1+frac2)
                   Flow_inter(kl) = (frac1*q_tot_c(var1) + &
                        frac2*q_tot_c(var1+2))/(frac1+frac2)
                   !> corner -> inter (2D plane)
                Endif
                kl=kl+1
             Enddo
             testz(1) = Abs(zgrid_max-zCent)
             testz(2) = Abs(zgrid_min-zCent)
             mZX(1) = Minval(testz) !> minimum distance in z-dimension
             testx(1) = Abs(xgrid_max-xCent)
             testx(2) = Abs(xgrid_min-xCent)
             mZX(2) = Minval(testx) !> minimum distance in y-direction
             rend = Minval(mZX)
             If (eqDis) Then
                testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
                rend = 1/2d0*Sqrt((dygrid)/Floor(testvar))
             Endif
             If (rend .Lt. r0) Then
                ja=1
                Goto 203
             Endif
             !> new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
             If (eqDis) Then !> middle of voxel as point of departure for bnd.cond.
                point(1,1) = (zgrid_max+zgrid_min)/2+rend
                point(1,2) = (xgrid_max+xgrid_min)/2
                point(2,1) = (zgrid_max+zgrid_min)/2
                point(2,2) = (xgrid_max+xgrid_min)/2+rend
                point(3,1) = (zgrid_max+zgrid_min)/2-rend
                point(3,2) = (xgrid_max+xgrid_min)/2
                point(4,1) = (zgrid_max+zgrid_min)/2
                point(4,2) = (xgrid_max+xgrid_min)/2-rend
             Else
                point(1,1) = zCent+rend
                point(1,2) = xCent
                point(2,1) = zCent
                point(2,2) = xCent+rend
                point(3,1) = zCent-rend
                point(3,2) = xCent
                point(4,1) = zCent
                point(4,2) = xCent-rend
             Endif
             Do j=1,Size(point,1)
                Do i=1,4
                   If (i.Eq.1) kl=1
                   If (i.Eq.2) kl=2
                   If (i.Eq.3) kl=5
                   If (i.Eq.4) kl=6
                   dis(j,i) = Sqrt((point(j,1)-zgrid(corner(kl)))**2+(point(j,2)-xgrid(corner(kl)))**2)
                   PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i)
                   Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                   !> inter -> point (distance)
                End Do
             End Do
          Endif
          !> IF CP is located in the plane (most of the cases)
210       Continue   !> one root approach continues here if applicable
          If (checker1 .Eq. 1 .Or. (extra_checker)) Then
             If (extra_checker) Then
                extra_checker=.False.
                Goto 200
             Endif
             !-----------------------------------------------------------------------------------------------------------------
             !> Determine the soil cylinder
             !> -----------------------------------------------------------------------------------------------------------------
             !> after outer boundary conditions at the new circle surrounding the cp are known for each allignment, 
             !> calculate the average PH and flow for these boundary conditions
             !> these average PH and Flow values (meanPHs and meanFlows) at the outer radius are the boundary conditions
             !> for the analytical approach
             Do j=1,4
                mPHs(j) = Sum(PHpoint_rel(j,:))/Sum(1/dis(j,:))
                mFlows(j) = Sum(Flowpoint_rel(j,:))/Sum(1/dis(j,:))
             End Do
             meanPHs = (Sum(mPHs(1:4))/4)
             meanFlows = (Sum(mFlows(1:4))/4)
             !> fractioned value should be used (subsegmentation)
             meanPHs = meanPHs*w_sub(irecn,isub,ipl)
             meanFlows = meanFlows*w_sub(irecn,isub,ipl)
             !-----------------------------------------------------------------------------------------------------------------
             !> partitioning the radii used for the anlytical approach (in log scale); closer
             !> radii near the root, larger radii at outer radius
             !> -----------------------------------------------------------------------------------------------------------------
             !> spatial grid
             Call logspace(r0,rend,n-1,y)
             r=y(1:n-1)
             rho = rend/r0
             !-----------------------------------------------------------------------
             !
             !> Calculate BOUNDARY conditions for analytical solution: Phi_rend(based on meanPHs),q_end,q_root
             !>
             !> -----------------------------------------------------------------------------------------------------------------
             !>  Calculate Phi_rend -> obtained from the PH boundary condition (the PH at outer radius of microscopic model)
             !> -----------------------------------------------------------------------------------------------------------------
             Call logspace(hx_min,meanPHs,n,y)
             h_mat = -y
             Call SetMat_ana(h_mat,theta,K,C)
             !> calculate matric flux potential integral K(h)dh -> numerical integration
             !> cumsum = cumulative sum
             arr=Abs(h_mat(2:Size(h_mat))-h_mat(1:Size(h_mat)-1))*(K(2:Size(K))+K(1:Size(K)-1))/2
             Do i=1,Size(arr)
                cumsum(i) = Sum(arr(1:i))
             End Do
             !> last value from integral is at r=rend
             Phi_rend=cumsum(Size(cumsum))
             !--------------q_root------------------------
             !> get current voxel number; for a loop over all the nodes in current voxel, get their irecn and isub values; calculate q_root for all root nodes
             If (ave) Then
                imin = voxel_no(irecn,isub)
                Do ki=1,no_voxels(imin)
                   irecn_tmp = voxel_node(imin,1,ki)
                   isub_tmp = voxel_node(imin,2,ki)
                   PHr_sub(irecn_tmp,isub_tmp) = w_sub(irecn_tmp,isub_tmp,ipl)*PHr(irecn_tmp+1,ipl)!you should not wietugh with w_sub
                   q_root_tmp(ki)= Lr(irecn_tmp,ipl) * (PH_root_sub(irecn_tmp,isub_tmp) &
                        - PHr_sub(irecn_tmp,isub_tmp))
                Enddo
                !===============================================================================
                !> chose the way how q_root is defined:
                !> 1) sum of all root flows, consider as one sucking root
                !> PHr_sub is not used, then only to check wheter or not limiting conditions have been reached -> see PHi
                PHr_sub(irecn,isub) = w_sub(irecn,isub,ipl)*PHr(irecn+1,ipl)
                !1)
                q_root = Sum(q_root_tmp(1:ki)) !> total uptake -> one big root
                !!-------------------------------------------------------------------------------
             Elseif (.Not.(ave)) Then
                !> PHxylem = PHr; PHr_sub = subsegment xylemPH
                PHr_sub(irecn,isub) = w_sub(irecn,isub,ipl)*PHr(irecn+1,ipl)
                !> ALCULATE BOUNDARY CONDITION q_root
                !> should be fractioned
                q_root= Lr(irecn,ipl) * (PH_root_sub(irecn,isub) - PHr_sub(irecn,isub))
             Endif
             !> CALCULATE BOUNDARY CONDITION q_out
             !> steady state -> condition for doussan matrix
             If (meanFlows .Le. q_root*(r0/rend)) Then
                q_out = meanFlows
             Else
                q_out = q_root*(r0/rend)
             Endif
             !de willigen, q_out=zero
             !if (eqDis) q_out=0
             !**************************************************************************
             !> ------------------------------ calculate PHI -----------------------------
             !> **************************************************************************
             !-----------------------------------------------------
             !> actual PHi calculation; hx_min equals lower boundary integral from flux density -> Phi_r_root = 0
             !> -----------------------------------------------------
             If (PHr_sub(irecn,isub) .Le. hx_min*w_sub(irecn,isub,ipl)) Then 
                !> if PHsoil-root interface = PHxylem limiting -> flow is zero (ultimate case limiting root PH
                Phi = ((Phi_rend + q_out*rend/2*Log(1/rho) ) / &
                     ( rho**2-1 + 2*rho**2*Log(1/rho) ) ) * (r**2/r0**2-1 + 2*rho**2*Log(r0/r)) &
                     + rend*q_out*Log(r/r0)
             Else !> always true, also when hx_min is reached for PHxylem. If finally PHsoil-root inter = PHxylem (=hx_min), see above
                Phi = Phi_rend + (q_root*r0 - q_out*rend)* &
                     (r**2/r0**2/(2*(1-rho**2)) + rho**2/(1-rho**2)*(Log(rend/r)-1/2.)) &
                     + q_out*rend*Log(r/rend)
             End If
             counter=1
             !> determine PH from Phi in microscopic domain
             Do l=1,Size(Phi)
                Do w=counter,Size(Phi_mat)
                   If (Phi_mat(w) .Lt. Phi(l)) Then
                      counter = w;
                   Elseif (Phi_mat(w) .Gt. Phi(l)) Then
                      Exit
                   End If
                End Do
                PH_micro(l)=h_mat2(counter)
             End Do
             PH_micro2(irecn,1:Size(Phi),isub)=PH_micro
             !> water potential criterion; check B*K_ave (soil) with Lr*r0 (see MS)
             If (PH_micro(n-1)-PH_micro(1) .Eq. 0) Then
                k_ave(irecn,isub) = 0
                ja=1
                Goto 203
             Else
                k_ave(irecn,isub) = (Phi(n-1)-Phi(1))/(PH_micro(n-1)-PH_micro(1))
             Endif
             B(irecn,isub)=2*(1-rho**2)/(-2*rho**2*(Log(rho)-1/2.)-1)
             If (k_ave(irecn,isub)*B(irecn,isub) .Lt. Lr(irecn,ipl)*r0) Then
                countja=.True.
             Endif
             !> set switchcriterion; if B*Ks > Lr*r0
             !> set switchcriterion; if B*Ks < Lr*r0
             !> set switchcriterion; if B*Ks +_ 5% < > Lr*r0
             If ((ave) .And. switchcriterion.Eq.1) PH_micro(1)=meanPHs+ &
                  r0/k_ave(irecn,isub)*q_out*rho*Log(1/rho)+r0/ &
                  (k_ave(irecn,isub)*B(irecn,isub))*q_out*rho
             If (k_ave(irecn,isub)*B(irecn,isub) .Lt. Lr(irecn,ipl)*r0+Lr(irecn,ipl)*r0*5/100 .And. k_ave(irecn,isub)*B(irecn,isub) .Gt. &
                  Lr(irecn,ipl)*r0-Lr(irecn,ipl)*r0*5/100) Then
                switchCriterion=2
                PH_micro(1)=meanPHs/2+PHr(irecn,ipl)/2+1/(Lr(irecn,ipl)*2)*q_out*rho*Log(1/rho)+ &
                     1/(Lr(irecn,ipl)*2)*q_out*rho
             Elseif (k_ave(irecn,isub)*B(irecn,isub) .Lt. Lr(irecn,ipl)*r0-Lr(irecn,ipl)*r0*5/100) Then
                switchCriterion=3
                PH_micro(1)=PHr(irecn,ipl)+B(irecn,isub)/Lr(irecn,ipl)*q_out*rho*Log(1/rho)+1/Lr(irecn,ipl)*q_out*rho
             Endif
             If (switchCriterion .Ne. 1 .And. (.Not.(tcheck)) ) Then
                tcheck2=.True.
                tcheck=.True.
             Endif
             If (PH_micro(1).Gt.-0.8) Then
                ja=1
                Goto 203
             Endif
             !>if one big root approach is taken:
             !>
             !> calculate the radius r at which the real centerpoint is located. Do not take the z-height
             !> into account. Find at which position in the radius array this value matches (sort of)
             !> take that index value in PH_micro and assign that PH as interface PH of that center point of a root segment.
200          If (ave) Then
                pos_cent=Sqrt((cent(irecn,1,isub,ipl)-cp_mean(irecn,1,isub,ipl))*(cent(irecn,1,isub,ipl)-cp_mean(irecn,1,isub,ipl))&
                     + (cent(irecn,2,isub,ipl)-cp_mean(irecn,2,isub,ipl))*(cent(irecn,2,isub,ipl)-cp_mean(irecn,2,isub,ipl))) 
                !> calculate radius from original cp to center of a voxel in 2D plane
                Do i=1,Size(r)-1
                   If (pos_cent .Eq. 0_dp) Then
                      pos=1
                      Exit
                   Elseif (pos_cent .Le. r(i+1) .And. pos_cent .Ge. r(i)) Then
                      pos=i+1
                      Exit
                   Endif
                Enddo
                If (checker1.Eq.1) Then
                   PH_root_sub(irecn,isub) = PH_micro2(irecn,pos,isub)
                Else
                   !> denote PH_int to all root nodes within this voxel independent of distance
                   PH_root_sub(irecn,isub) = PH_micro2(irecnTemp,pos,isubTemp)/pond_subTemp*w_sub(irecn,isub,ipl) 
                   !> take difference in w_sub between cps into account
                   If (segsur(irecnTemp,ipl) .Eq. 0 .Or. seglen(irecnTemp,ipl) .Eq. 0) Then
                      r0 = 5e-2 !> minimal root radius (in cm)
                   Else
                      r0 = 1/(2*pi)*segsur(irecnTemp,ipl)/seglen(irecnTemp,ipl)
                   End If
                   If (k_ave(irecnTemp,isubTemp)*B(irecnTemp,isubTemp) .Lt. Lr(irecnTemp,ipl)*r0) Then
                      countja=.True.  !> to check how many root nodes are under "limiting conditions"
                   Endif
                Endif
             Endif
             If (.Not.(ave)) Then
                PH_root_sub(irecn,isub) = PH_micro(1)
             Endif
             !-----------------------------------------------------------------------
             !> Give every soil-root interface node its correct PH value
             !> -----------------------------------------------------------------------
             PH_root_sub(irecn,isub) = PH_root_sub(irecn,isub) + zCent*w_sub(irecn,isub,ipl)!calc PH at interface with gravity
          Else!checker 
             !> some segments are divided in subsegments in which a cp is located at the edge of the plane (and/or cube),use average approach of original R-SWMS
             sumW=sum_dis(irecn,isub,ipl)
             psub=w_sub(irecn,isub,ipl)
             PH_root_sub(irecn,isub) = ((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
                  (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
                  zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
                  (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
                  zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
                  (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
          End If
          !---------------------------------------------------------------------------
          !> ph at interface
          PHs(irecn,ipl)= PHs(irecn,ipl) + PH_root_sub(irecn,isub)

          !> add salinity potential 
          If (lSalinity) Then
             PHs(irecn,ipl)=PHs(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
             PHo(irecn,ipl)=PHo(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
          Endif

203       If (ja.Eq.1) Then !> average method used in some cases
             If (irecn .Eq. 0) Then
                r0=5e-2
             Else
                r0 = 1/(2*pi)*segsur(irecn,ipl)/seglen(irecn,ipl)
             Endif
             temp=Lr(irecn,ipl)*r0
             Do j=1,8
                h_matK(j)=hnew(corner(j))
                Call SetMat_singlevalue(h_matK(j),MatNum(corner(j)),theta_k,K_k,C_k,corner(j))
                temp2(j)=K_k
             Enddo
             If (Sum(temp2)/8 .Lt. temp) Then
                countja=.True.
             Endif
             ja=0
             sumW=sum_dis(irecn,isub,ipl)
             psub=w_sub(irecn,isub,ipl)
             PH_root_sub(irecn,isub)=((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
                  (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
                  zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
                  (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
                  zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
                  (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
             PHs(irecn,ipl)=PHs(irecn,ipl)+PH_root_sub(irecn,isub)

             If(lSalinity) Then
                PHs(irecn,ipl)=PHs(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
                PHo(irecn,ipl)=PHo(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
             Endif

          Endif
       End Do !> end do-loop over subsegments
       !> calculate matrix Q (doussan matrix)
       If (countja) count_nodes=count_nodes+1
       countja=.False.
    End Do

   Qd(0:nrec(ipl),ipl) = Qi(0:nrec(ipl),ipl)*PHs(0:nrec(ipl),ipl)+Q_bc1(0:nrec(ipl),ipl)
   Qd(nrec(ipl)+1:2*nrec(ipl)+1,ipl) = Qi(0:nrec(ipl),ipl)*PHs(0:nrec(ipl),ipl)+Q_bc2(0:nrec(ipl),ipl)
  End Subroutine analytical_approach
  !=====================================================================================
  Subroutine flux_node(Flow_corner,corner,q_tot_c)
    Use typedef
    Use GridData, Only: xgrid,ygrid,zgrid,dxgrid,dygrid,dzgrid
    Use SolData, Only: MatNum, hnew
    Implicit None

    Integer(ap) :: i,j,Flow_corner(1:8,1:6)
    Integer(ap) ::corner(1:8)
    Real(dp) :: q_tot_c_t
    Real(dp) :: theta_single,K_single,C_single,Kn,Kcn,Kc,dz
    Real(dp),Intent(out) :: q_tot_c(1:8)

    q_tot_c(1:8) = 0
    !> loop over all soil nodes from cubic (sub cubic -> subsegment)
    Do i=1,8
       !> calculate soil hydraulic conductivity in the corner node and the surrounding node (in next loop)
       Call SetMat_singlevalue(hnew(corner(i)),MatNum(corner(i)),theta_single,K_single,C_single, corner(i))
       Kc = K_single
       !> loop over all surrounding nodes
       Do j=1,6
          If (Flow_corner(i,j).Eq.0) Then
             q_tot_c_t = 0
             Goto 40
          Endif
          Call SetMat_singlevalue(hnew(Flow_corner(i,j)),MatNum(corner(j)),theta_single,K_single,C_single,Flow_corner(i,j))
          Kn = K_single
          Kcn = (Kc+Kn)/2._dp !> arithmetic mean soil hydraulic conductivity
          !> distance between the two nodes evaluated
          If (xgrid(corner(i))-xgrid(Flow_corner(i,j)) .Ne. 0) Then
             dz=dxgrid
          Else If (ygrid(corner(i))-ygrid(Flow_corner(i,j)) .Ne. 0) Then
             dz=dygrid
          Else If (zgrid(corner(i))-zgrid(Flow_corner(i,j)) .Ne. 0) Then
             dz=dzgrid
             !> -> q_tot_c_t = 0 (no flow from z-direction taken into account
          End If
          !> flow from corner node to surrounding node
          q_tot_c_t = - Kcn * (hnew(corner(i)) - hnew(Flow_corner(i,j)))/dz
          !> take the sum over all the calculates flows from and to the corner node
          !> q_tot_c has size 8 (corner nodes of cube)
40        Continue
          q_tot_c(i) = q_tot_c(i) + q_tot_c_t
       End Do
    End Do
    Return
  End Subroutine flux_node
 !=====================================================================================
  Subroutine average_method(ipl,iter_root)
  !### when lOld is true (typical case), estimate Soil PH
  ! including or not osmotic potential###
    Use typedef
    Use Watfun
    Use Paramdata,Only: MaxNod,pi,lSalinity
    Use doussanmat, Only: PHs,nsub,w_sub,sum_dis,PH_root_sub,Qi,Qd,Q_bc1,Q_bc2,w_dis,&
         loc_Q,Lr,count_nodes,cube_i,PHs_osmotic,PHo,Joutr,PHsri,voxel_node,no_voxels,&
         mat_seg,vox_seg,w_sub,PHbulk
    Use GridData, Only: zgrid, dxGrid, dyGrid
    Use RootData, Only: nrec,segsur,seglen,sigma,lKdrop,segrad,TypeKdrop
    Use SolData, Only: MatNum, hnew,par

    Implicit None

    Integer(ap) :: isub,j,irecn_tmp,isub_tmp,ki,imin,iter_root
    Integer(ap) :: corner(1:8)
    Integer(ap) :: irecn,ipl,M
    Real(dp) :: sumW,psub
    Real(dp),Dimension(1:8) :: pond_i
    Real(dp),Dimension(1:8) :: h_mat,temp2
    Real(dp) :: K,C,theta,r0,temp,phi_sri,phi_b,dist2r
    Logical :: countja

    countja=.False.
    count_nodes=0
    Do irecn=1,nrec(ipl)
       If (irecn .Eq. 0) Then
          r0=5e-2
       Else
          r0 = 1/(2*pi)*segsur(irecn,ipl)/seglen(irecn,ipl)
       Endif
       temp=Lr(irecn,ipl)*r0
       !> initialize soil PH matrix
       PHs(irecn,ipl)=0._dp
       PHo(irecn,ipl)=0._dp
       Do isub=1,nsub(irecn,ipl)
          pond_i(:)=w_dis(irecn,:,isub,ipl)
          corner(:)=loc_Q(irecn,:,isub,ipl)
          Do j=1,8
             h_mat(j)=hnew(corner(j))
             M=MatNum(corner(j))
             Call SetMat_singlevalue(h_mat(j),M,theta,K,C,corner(j))!associate k
             temp2(j)=K
          Enddo
          If (Sum(temp2)/8 .Lt. temp) Then!>Ksoil<Lr
             countja=.True.
          Endif
          sumW=sum_dis(irecn,isub,ipl) !> total inverse distance to node 
          psub=w_sub(irecn,isub,ipl)   !> subsegment length/total segment length ratio
          PH_root_sub(irecn,isub)=((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
               (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
               zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
               (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
               zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
               (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
          PHs(irecn,ipl)=PHs(irecn,ipl)+PH_root_sub(irecn,isub)   
          if ((lKdrop).and.((TypeKdrop).eq.2).and.(iter_root.GT.1)) Then !update of PHsri at each iteration
              if (Joutr(irecn,ipl)/segsur(irecn,ipl).GT.0.001) Then !from previous it
                  phi_b=Fmfp_soiltab(PHs(irecn,ipl),mat_seg(irecn,ipl),Par(:,1))
                  phi_sri=phi_b+(Joutr(irecn,ipl)/segsur(irecn,ipl)*segrad(irecn,ipl)*log(segrad(irecn,ipl)/dxGrid))
                  imin = vox_seg(irecn,ipl)! voxel correspondant
                  dist2r=sqrt((dxGrid*dyGrid)/(no_voxels(imin)*pi))! averaged distance btw roots in that voxel
                  Do ki=1,no_voxels(imin)!when more than one root in a voxel, additive approach of Doussan
                      irecn_tmp = voxel_node(imin,1,ki)
                      isub_tmp = voxel_node(imin,2,ki)
                      If ((Joutr(irecn_tmp,ipl).GT.0.001).AND.(irecn_tmp.NE.irecn)) Then
                           phi_sri=phi_sri+(Joutr(irecn_tmp,ipl)/segsur(irecn_tmp,ipl))*w_sub(irecn_tmp,isub_tmp,ipl)*segrad(irecn_tmp,ipl)*log(dist2r/dxGrid)
                      Endif
                  Enddo
                  If (Abs(phi_sri-phi_b).GT.0.0001) Then
                     PHbulk(irecn,ipl)=PHs(irecn,ipl)
                     PHsri(irecn,ipl)=Fh_from_mfp_soiltab(phi_sri,mat_seg(irecn,ipl))
! print*,irecn,no_voxels(imin),imin, 'PHsri',PHsri(irecn,ipl),'PHk', PHbulk(irecn,ipl),dist2r,phi_sri
!print *, iter, iter_root, 
                     PHs(irecn,ipl)=PHsri(irecn,ipl)
                  Endif
              Endif
          Else!PHs remains the same=> PHS=PHbulk=PHSri
             PHbulk(irecn,ipl)=PHs(irecn,ipl)
             !PHsri(irecn,ipl)=PHs(irecn,ipl)
          Endif
          If(lSalinity) Then
             PHs(irecn,ipl)=PHs(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
             PHo(irecn,ipl)=PHo(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
          Endif
       End Do
       If (countja) count_nodes=count_nodes + 1
       countja=.False.
    End Do
!> calculate matrix Qd
    Qd(0:nrec(ipl),ipl) = Qi(0:nrec(ipl),ipl)*PHs(0:nrec(ipl),ipl)+Q_bc1(0:nrec(ipl),ipl)
    Qd(nrec(ipl)+1:2*nrec(ipl)+1,ipl) = Qi(0:nrec(ipl),ipl)*PHs(0:nrec(ipl),ipl)+Q_bc2(0:nrec(ipl),ipl)

  End Subroutine average_method
 
!***********************************************************************************
!### create a log distributed array ###
  Subroutine logspace(d1,d2,n,y)
    Use typedef
    Implicit None

    Integer(ap),Intent(in) :: n
    Integer(ap) :: i
    Real(dp) :: qq(1:n),Y_(1:n)
    Real(dp) :: d1, d2, d31, d32
    Real(dp), Intent(out) :: y(1:n)

    !>logspace algorithm
    d31 = Log10(Abs(d1))
    d32 = Log10(Abs(d2))

    Do i=1,n-1
       qq(i) =d31+(i-1)*(d32-d31)/(n-1)
    End Do

    !>put in one array
    y_(1:n-1) = qq(1:n-1)
    y_(n) = d32

    !> convert logscale to normal values
    y = (10)**y_(1:n)
    Return
  End Subroutine logspace
 !***********************************************************************
  Subroutine SetMat_ana(hTemp,theta,K,C)
    Use typedef
    Use GridData
    Use SolData, Only: par
    Use WatFun
    Use Doussanmat, Only : n
    Use RhizoData, Only : lRhizo
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap):: i,M
    Real(dp), Intent(in) :: hTemp(1:n)
    Real(dp), Intent(out) :: theta(1:n),C(1:n),K(1:n)
    Real(dp):: Ti(1:n)

    If (lRhizo) Then
       Call stop_program('Someone callSetMat_ana while RhizoModel')
    Endif
    Do i=1,n
       M=nMat
       !> nodal conductivity values:
       K(i)=FKP(hTemp(i),par(:,M),1)*Bxy(i)
       C(i)=FCP(hTemp(i),par(:,M))*Dxy(i)/Axy(i)
       Ti(i)=Fth(hTemp(i),par(:,M))
       theta(i)=par(2,M)*Exy(i)+(Ti(i)-par(2,M))*Dxy(i) 
    End Do

    Return
  End Subroutine SetMat_ana
 !***********************************************************************
  Subroutine SetMat_anaLibr(hTemp,theta,K,C)
    Use typedef
    Use GridData
    Use SolData, Only: par
    Use WatFun
    Use Doussanmat, Only: nLibr
    Use RhizoData, Only : lRhizo
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: i,M
    Real(dp), Intent(in) :: hTemp(1:nLibr)
    Real(dp), Intent(out) :: theta(1:nLibr),C(1:nLibr),K(1:nLibr)

    If (lRhizo) Then
       Call stop_program('Someone call SetMat_anaLibr while RhizoModel')
    Endif
    M=nMat
    Do i=1,nLibr
       !> nodal conductivity values:
       K(i)=FKP(hTemp(i),par(:,M),1)
       C(i)=FCP(hTemp(i),par(:,M))
       theta(i)=Fth(hTemp(i),par(:,M))
    End Do

    Return
  End Subroutine SetMat_anaLibr
  !***********************************************************************
  Subroutine SetMat_singlevalue(hTemp,M,theta,K,C,soilNodeAddress)
    Use typedef
    Use GridData
    Use SolData, Only: par,soiltab
    Use WatFun
    Implicit None

    Integer(ap), Intent(in) :: M
    Integer(ap) :: soilNodeAddress
    Real(dp), Intent(in) :: hTemp
    Real(dp), Intent(out) :: theta,C,K

    !> nodal conductivity values:
    If (soiltab) Then
       K=FKP_soiltab(hTemp,M)
       C=FCP_soiltab(hTemp,M)
       theta=Fth_soiltab(hTemp,M)
    Else
       K=FKP(hTemp,par(:,M), soilNodeAddress)
       C=FCP(hTemp,par(:,M))
       theta=Fth(hTemp,par(:,M)) 
    End If
    Return
  End Subroutine SetMat_singlevalue
  !**************************************************************************************
  ! weigth for each node in fucntion of their neighbours
  Subroutine CalcWnodes()
    Use typedef
    Use GridData
    Implicit None

    Integer(ap) :: iE

    !> calculate actual overall transpiration rate from soil domain:
    Wn=0.0_dp
    Do iE=1,nElm
       !> assign to Wn the number of times a cuboid corner node is called
       Wn(elmnod(1,iE))=Wn(elmnod(1,iE))+1
       Wn(elmnod(2,iE))=Wn(elmnod(2,iE))+1
       Wn(elmnod(3,iE))=Wn(elmnod(3,iE))+1
       Wn(elmnod(4,iE))=Wn(elmnod(4,iE))+1
       Wn(elmnod(5,iE))=Wn(elmnod(5,iE))+1
       Wn(elmnod(6,iE))=Wn(elmnod(6,iE))+1
       Wn(elmnod(7,iE))=Wn(elmnod(7,iE))+1
       Wn(elmnod(8,iE))=Wn(elmnod(8,iE))+1
    End Do
    Wn=Wn/8._dp

    Return
  End Subroutine CalcWnodes
  !********************************************************************************
  !> ### for continuous domains determine length of the root that is growing back into the
  !> other side of the soil domain ###
  Subroutine roottrans(xt,yt,inode,isub,ipl)
    Use typedef
    Use GridData, Only: nex,ney,dxgrid,dygrid
    Use DoussanMat, Only: transroot
    Use DomData
    Implicit None

    Integer(ap), Intent(in) :: inode,ipl,isub
    Real(dp), Intent(in) :: xt,yt
    Real(dp) :: xd,yd

    xd = xt
    yd = yt
    transroot(inode,1:2,isub,ipl)=0
    Do While (xd.Ge.xmax) 
       xd=xd-nex*dxgrid
       transroot(inode,1,1,ipl)=transroot(inode,1,isub,ipl)-1
    End Do
    Do While (xd.Lt.xmin)
       xd=xd+nex*dxgrid
       transroot(inode,1,1,ipl)=transroot(inode,1,isub,ipl)+1
    End Do
    Do While (yd.Ge.ymax)
       yd=yd-ney*dygrid
       transroot(inode,2,1,ipl)=transroot(inode,2,isub,ipl)-1
    End Do
    Do While (yd.Lt.ymin)
       yd=yd+ney*dygrid
       transroot(inode,2,1,ipl)=transroot(inode,2,isub,ipl)+1
    End Do

  End Subroutine roottrans
  !******************************************************************************** 
  !> ### for continuous domains determine length of the root that is growing back into the
  !> other side of the soil domain ###
  Subroutine tiptrans(xt,yt,inode,isub,ipl)
    Use typedef
    Use GridData, Only: nex,ney,dxgrid,dygrid   
    Use DoussanMat, Only:transtip 
    Use DomData
    Implicit None

    Integer(ap), Intent(in) :: inode,ipl,isub
    Real(dp), Intent(in) :: xt,yt
    Real(dp) :: xd,yd
    		
    xd = xt
    yd = yt
    transtip(inode,1:2,isub,ipl)=0
    Do While (xd.Ge.xmax)
       xd=xd-nex*dxgrid
       transtip(inode,1,1,ipl)=transtip(inode,1,isub,ipl)-1
    End Do
    Do While (xd.Lt.xmin)
       xd=xd+nex*dxgrid
       transtip(inode,1,1,ipl)=transtip(inode,1,isub,ipl)+1
    End Do
    Do While (yd.Ge.ymax)
       yd=yd-ney*dygrid
       transtip(inode,2,1,ipl)=transtip(inode,2,isub,ipl)-1
    End Do
    Do While (yd.Lt.ymin)
       yd=yd+ney*dygrid
       transtip(inode,2,1,ipl)=transtip(inode,2,isub,ipl)+1
    End Do

  End Subroutine tiptrans
  !********************************************************************************
  Subroutine calc_actual_climate(t)
    Use typedef
    Use EnviData
    Implicit None

    Real(dp), Intent(in)::t
    Real(dp) :: temp_PPFD,temp_VPD,temp_temperature
    Real(dp) :: temp_s,temp_rho,temp_gamma,temp_conversion,temp_lambda

    Call interp_linear (nclimaticdata,time_climate,PPFD,t,temp_PPFD)
    Call interp_linear (nclimaticdata,time_climate,VPD,t,temp_VPD)
    Call interp_linear (nclimaticdata,time_climate,Temperature,t,temp_temperature)

    Call interp_linear (nclimaticdata,time_climate,s_f,t,temp_s)
    Call interp_linear (nclimaticdata,time_climate,rho_f,t,temp_rho)
    Call interp_linear (nclimaticdata,time_climate,gamma_f,t,temp_gamma)
    Call interp_linear (nclimaticdata,time_climate,conversion_f,t,temp_conversion)
    Call interp_linear (nclimaticdata,time_climate,lambda_f,t,temp_lambda)


    VPD_t=temp_VPD
    PPFD_t=temp_PPFD
    Rn_t=PPFD_t/2.
    temperature_t=temp_temperature

    s_t=temp_s
    rho_t=temp_rho
    gamma_t=temp_gamma
    conversion_t=temp_conversion
    lambda_t=temp_lambda
    !print*,temp_t,VPD_t,PPFD_t,temperature_t

  End Subroutine calc_actual_climate
  !********************************************************************************
  Subroutine interp_linear (data_num, t_data, p_data, &
       t_interp, p_interp )
    Use typedef
    Implicit None

    Integer(ap) :: data_num
    Integer(ap) :: left
    Integer(ap) :: right
    Real(dp) :: p_data(data_num)
    Real(dp) :: p_interp
    Real(dp) :: t
    Real(dp) :: t_data(data_num)
    Real(dp) :: t_interp



    t = t_interp
    !
    !  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
    !  nearest to, TVAL.
    !
    Call r8vec_bracket ( data_num, t_data, t, left, right )

    p_interp = &
         ( ( t_data(right) - t                ) * p_data(left)   &
         + (                 t - t_data(left) ) * p_data(right) ) &
         / ( t_data(right)     - t_data(left) )


    Return
  End Subroutine interp_linear
  !********************************************************************************
  Subroutine r8vec_bracket(n,x,xval,left,right)
    Use typedef

    Implicit None

    Integer(ap), Intent(in) :: n
    Integer(ap), Intent(out) :: left
    Integer(ap), Intent(out) :: right
    Integer(ap) :: i
    Real(dp), Intent(in) :: x(n)
    Real(dp), Intent(in) :: xval

    Do i = 2, n - 1

       If ( xval < x(i) ) Then
          left = i - 1
          right = i
          Return
       End If

    End Do

    left = n - 1
    right = n

    Return
  End Subroutine r8vec_bracket
!*************************************************************************************
!### Update Doussan matrix at each tome solveroot is called in case gaps,
!  AQP, Kdrop are activated ### 
  Subroutine UpdateDou(ipl)
    Use TypeDef
    Use SparseMatrix
    Use DoussanMat, Only: plantmatrix,curr_bctp,Khr,Lr,Qi
    Use RootData, Only: nrec,nbr,ibrseg,irecpr,seglen,segsur
    Use MPIutils, Only: stop_program

    Implicit None
    Integer(ap) :: ibr,iprvn,irecn,ifoln,err,ipl
    Logical :: n_apex,run

  
          err=SM_allocate(plantmatrix(ipl), nrec(ipl)+1, nrec(ipl)+1)
          If(err/=0) Call stop_program('Could not create plantmatrix')
          Do ibr=1,nbr(ipl) !> all plants have same number of root segments!
             n_apex=.False.
             !> find the tip segment of the branch 'ibr'
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.ibr)
                irecn=irecn-1
             End Do
             If (seglen(irecn,ipl)<1.E-20) Then ! skip this segment too small to be taken
                ifoln=irecn ! following node ID
                irecn=irecpr(irecn,ipl) !current node ID
             Else
                n_apex=.True.!ifoln does not exist if not continu
             Endif
             If (irecn==0) Then !there exists a branch ibr but not yet any segment!
                run=.False.
             Else
                run=.True.
             Endif
             !> then the rest of the branch up to the seed or the embranchment
             Do While (run)
                !> "upper" node
                iprvn=irecpr(irecn,ipl)
                !> several changes -> multiple roots
                err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                !> if apex (bottom part of root)
                If (n_apex) Then
                   Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn,ipl)/seglen(irecn,ipl)+Lr(irecn,ipl)*segsur(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Else
                   Err=SM_set(Plantmatrix(ipl), irecn+1,ifoln+1,-Khr(ifoln,ipl)/seglen(ifoln,ipl),.False.) !row, col,value
                   Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn,ipl)/seglen(irecn,ipl)+Khr(ifoln,ipl)/seglen(ifoln,ipl)+Lr(irecn,ipl)*segsur(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                Endif
                !> define 1st part of Q (Q=Qi.*PHsoil+Qbc) -> RHS
                Qi(irecn,ipl)=Lr(irecn,ipl)*segsur(irecn,ipl)
                If (iprvn==0) Then!> seed=first segment
                   If (curr_BCtp(ipl)==2) Then !> flux BC
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                   Else If (curr_BCtp(ipl)==1) Then !> PH BC
                      Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.True.)!position (1,1)
                      ! true means overwrite (delete first), no addition
                      Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.True.)!irecn+1,iprvn+1
                      If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                   Endif
                   run=.False.
                Elseif (ibrseg(iprvn,ipl).Ne.ibrseg(irecn,ipl)) Then !start of the branch but not from the seed -> gaat van onder na boven, duz iprvn (nieuwe positie) is nu niet van die ene branch maar van een side branch
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                   Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn,ipl)/seglen(irecn,ipl),.False.)
                   If(err/=0) Call stop_program('Could not insert element into plantmatrix')
                   run=.False.
                Endif
                !> definition for the next run of the loop
                ifoln=irecn
                irecn=iprvn
                !> from here, not an apex
                n_apex=.False.
             End Do ! loop on branch nodes
          End Do ! loop on root branches
		  
  End Subroutine UpdateDou
  
 !*************************************************************************************  
    Subroutine update_historics(t,dt)
    Use typedef
    Use EnviData
    Use TardieuData
    Use PlntData, Only: Tact,PHcollar
    Use tmctrl, Only : dtmin,tout
    Use RootData, Only : Krs
    Implicit None

    Integer(ap) :: current_day
    Real(dp), Intent(in)::t,dt
    Real(dp) :: temp_day
    Character :: file*15

    If (.Not.Allocated(Global_transpi_and_Psi)) Then
       Allocate(Global_transpi_and_Psi(Nint((Maxval(tout)-t_initial)/dtmin),4))
       Global_transpi_and_Psi=0
    Endif
    !print*,NINT((maxval(tout)-t_initial)/dtmin),dtmin,t_initial,maxval(tout)
    Global_transpi_and_Psi(count_glob,1)=t
    Global_transpi_and_Psi(count_glob,2)=Tact(1)
    Global_transpi_and_Psi(count_glob,3)=PHcollar
    Global_transpi_and_Psi(count_glob,4)=dt
    count_glob=count_glob+1
    !print*,'test',maxval(Global_transpi_and_Psi(:,2))

    Call interp_linear (nclimaticdata,time_climate,days,t,temp_day)
    current_day=Nint(temp_day)  

    Write (file,'(A15)')'out/Tardieu.out'
    Open (UNIT=10,FILE=file,STATUS='OLD',POSITION='APPEND')
    Write (10,'(2(1pE11.4,1X),12(1pE11.4,1X))')&
         t,Krs,Krs_circad,ampli(current_day),Krs_transpi,ABA_collar,PHcollar,Hleaf,Hbundle,Hcell,Jxc_t,Vcel
    Close (10)

  End Subroutine update_historics

  !*************************************************************************************
  Subroutine update_Krs(t)
    Use typedef
    Use TardieuData
    Use EnviData
    Use ParamData, Only : pi
    Use RootData, Only: Krs,Hseq,Kcomp
    Use PlntData, Only: Tact,PHcollar
    Implicit None

    Integer(ap) :: i,current_day
    Real(dp), Intent(in)::t
    Real(dp) :: sum_dt,sum_Tact,hourphoto,temp_day,Jminus40=0.  
    Real(dp) :: Gxl_transpi,Gxl_circad,Gc_transpi,Gc_circad,Gstem_transpi,Gstem_circad
    Real(dp) :: factor_ABA

    ! calc ABA in the collar
    ABA_collar=ABA_collar_fun(Hseq,Tact(1),ABAconstit,ABAa,ABAb)

    ! check the current day of simulation
    Call interp_linear (nclimaticdata,time_climate,days,t,temp_day)
    current_day=Nint(temp_day)  

    ! calc photosynthesis hour
    hourphoto = (t-t_initial) - datad(current_day,2)
    !print*,current_day,hourphoto


    If (current_day.Gt.2) Then
       ampli(current_day) = Psixylmax(current_day-2) - Psixylmin(current_day-2)
    Endif

    Krs_circad = TaucircadGr*ampli(current_day)*Cos(hourphoto*24*pi/12)
    Gxl_circad= TaucircadGxl*ampli(current_day)*Cos(hourphoto*24*pi/12)
    Gc_circad= TaucircadGc*ampli(current_day)*Cos(hourphoto*24*pi/12)
    Gstem_circad= TaucircadGstem*ampli(current_day)*Cos(hourphoto*24*pi/12)      

    If (.Not.(Allocated(Global_transpi_and_Psi))) Then
       Krs_transpi = Grmin
       Gxl_transpi = Gxlmin
       Gc_transpi = Gcmin
       Gstem_transpi = Gstemmin
    Else
       If (t.Lt.(Global_transpi_and_Psi(1,1)+1.5/24.)) Then
          Jminus40=Global_transpi_and_Psi(count_glob-1,2)
          Krs_transpi = Min(Max(Grmin,TautranspiGr*Jminus40),Grmax)
          Gxl_transpi = Min(Max(Gxlmin,TautranspiGxl*Jminus40),Gxlmax)
          Gc_transpi = Min(Max(Gcmin,TautranspiGc*Jminus40),Gcmax)
          Gstem_transpi = Min(Max(Gstemmin,TautranspiGstem*Jminus40),Gstemmax)
       Else
          sum_dt=0.
          sum_Tact=0.
          Do i=1,count_glob-1
             If ((Global_transpi_and_Psi(i,1).Gt.(t-60./60./24.)).And.(Global_transpi_and_Psi(i,1).Lt.(t-20./60./24.))) Then
                sum_dt=sum_dt+Global_transpi_and_Psi(i,4)
                sum_Tact=sum_Tact+Global_transpi_and_Psi(i,2)*Global_transpi_and_Psi(i,4)
             Endif
          Enddo
          Jminus40=sum_Tact/sum_dt
          Krs_transpi = Min(Max(Grmin,TautranspiGr*Jminus40),Grmax)
          Gxl_transpi = Min(Max(Gxlmin,TautranspiGxl*Jminus40),Gxlmax)
          Gc_transpi = Min(Max(Gcmin,TautranspiGc*Jminus40),Gcmax)
          Gstem_transpi = Min(Max(Gstemmin,TautranspiGstem*Jminus40),Gstemmax)
       Endif
    Endif
    factor_ABA=(1+tauABAGr*Log(ABA_collar/20))
    Krs=Min(Max(Grmin,(Krs0 + Krs_transpi + Krs_circad)*factor_ABA),Grmax)
    Kcomp=Min(Max(Grmin,(Kcomp0 + Krs_transpi + Krs_circad )*factor_ABA),Grmax)
    Gxl=Min(Max(Gxlmin,(Gxl0 + Gxl_transpi + Gxl_circad)*factor_ABA),Gxlmax)
    Gc=Min(Max(Gcmin,(Gc0 + Gc_transpi + Gc_circad)*factor_ABA),Gcmax)
    Gstem=Min(Max(Gstemmin,(Gstem0 + Gstem_transpi + Gstem_circad)*factor_ABA),Gstemmax)

    !print*,Gxl,Gc,Gstem 
    If (Abs(PHcollar).Gt.0) Then
       Psixylmax(current_day) = Max(Psixylmax(current_day), Abs(PHcollar));
       Psixylmin(current_day) = Min(Psixylmin(current_day), Abs(PHcollar));
    Endif


  End Subroutine update_Krs

 !*********************************************************************************
  !> ### estimate stomatal conductance gs as a function of the Tardieu-Davies model
  !>  Tardieu et al., 2015: https://doi.org/10.1093/jxb/erv039 ###

  Subroutine find_stomatal_aperture(Hseq,Krs,gs_loc)    
    Use typedef
    Use TardieuData
    Use EnviData
    Implicit None

    Integer(ap) :: j 
    Real(dp), Intent(in) :: Hseq
    Real(dp),Intent(inout) :: gs_loc
    Real(dp) :: g1,g2,difcurrent,difinit,ABA_loc,Hcollar_loc,Hbundle_loc
    Real(dp) :: J_loc
    Real(dp) :: Tact_loc,Krs

    g1 = 0.01
    difcurrent = 1
    Do j=1,500     
       gs_loc = g1/conversion_t
       J_loc = 1000000*PM_fun(gs_loc,Sshaded,s_t,Rn_t,rho_t,Cp,ga,VPD_t*1000,lambda_t,gamma_t,conversion_t)
       Tact_loc = J_loc*86.4
       Hcollar_loc=PHcollar_fun(Hseq,Krs,Tact_loc)
       Hbundle_loc=Hcollar_loc-Tact_loc/Gstem-Tact_loc/Gxl
       ABA_loc = ABA_collar_fun(Hseq,Tact_loc,ABAconstit,ABAa,ABAb)
       g2=2*gs_fun(ABA_loc,Hbundle_loc,gsmin,gsalpha,gsbeta,gsdelta)
       difinit = difcurrent
       difcurrent=(g1-g2)**2
       !print*,j,difcurrent,difinit
       ! calcul de dif
       If (difcurrent.Lt.difinit)Then
          g1 = g1 + 0.001
       Else
          Exit
       Endif
       gs_loc = g2
    Enddo

  End Subroutine find_stomatal_aperture

  !***************************************************************************************

  Subroutine capacitance
    Use typedef
    Use EnviData
    Use TardieuData
    Use PlntData, Only: Tact,PHcollar
    Use RootData, Only: Hseq,Krs

    Implicit None

    Real(dp) :: dt_current
    Real(dp) :: Hcellcurrent,RWC

    If (.Not.Allocated(Global_transpi_and_Psi)) Then
       dt_current=0.
    Else
       dt_current=Global_transpi_and_Psi(count_glob-1,4)
    Endif

    If (count_glob.Eq.1) Then
       Hleaf=PHcollar-Tact(1)/Gstem
       Hbundle=PHcollar-Tact(1)/Gstem-Tact(1)/Gxl
       Hcell=Hbundle
       Vcel=Vres+ (Vsat - Vres)* (1/(1+(alphacap*(-Hcell*10/10000))**ncap))**(1-1/ncap)
       !print*,-Hcellcurrent+Hbundle,Jxc_t,dt_current,Vcel,'here'          
    Else
       Hcellcurrent=Hcell
       Hleaf=PHcollar-Tact(1)/Gstem
       Hbundle=PHcollar-Tact(1)/Gstem-Tact(1)/Gxl

       Jxc_t = (Gc * (-Hcellcurrent + Hbundle)) / (1 + Gc/Krs + Gc/Gstem + Gc/Gxl) 
       PHcollar=PHcollar_fun(Hseq,Krs,Tact(1)+Jxc_t)
       Hleaf = PHcollar-(Tact(1)+Jxc_t)/Gstem
       Hbundle = PHcollar-(Tact(1)+Jxc_t)/Gstem - (Tact(1)+Jxc_t)/Gxl
       Vcel = Max(0.0_dp,Vcel + Jxc_t*dt_current)
       !print*,-Hcellcurrent+Hbundle,Jxc_t,dt_current,Vcel
       If (Vcel.Ge.Vsat) Then
          Vcel=Vsat
          Hcellcurrent=0.
       Else
          RWC = (Vcel-Vres)/(Vsat-Vres)
          If (RWC.Lt.0.00001)Then
             Hcellcurrent=Hbundle
          Else
             Hcellcurrent = (-0.1*((1-RWC**(ncap/(ncap-1)))**(1/ncap)/(alphacap*RWC**(1/(ncap-1)))))*10000
          Endif
       Endif
       Hcell = Hcellcurrent

    Endif
    !print*,Hcell,Jxc_t*dt,Vcel

  End Subroutine capacitance


End Module Doussan
