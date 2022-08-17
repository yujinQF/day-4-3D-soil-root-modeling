!> \file RootTyp.f90
!! \brief ROOT GROWTH with RootTyp

!> Module RootTyp
Module RootTyp

Contains
  
  Subroutine RunRootTip(t,ipl)
    ! run RootTip for a given time (in days) from 0
    ! return naxes:= primary roots in RooTyp
    ! call C-functions, which return arrays of variables
    ! most of these are then saved under RootData module
    Use Typedef
    Use RootData
    Use GridData, Only: continu
    Use Output, Only: OutRoo
    Implicit None

    Integer(ap) :: n,simtime, n_nodes, n_meris,n_br_crea,n_br_del, n_axes,t
    Integer(ap) :: n_nodecorr,ipl !runroottip currently doesn´t work with more than one plant. It should be updated... (Couvreur dec 2009)
    Integer(ap), Allocatable :: node(:)
    Integer(ap), Allocatable :: prev_node(:)
    Integer(ap), Allocatable :: ordn(:)
    Integer(ap), Allocatable :: orda(:)
    Integer(ap), Allocatable :: agerec(:)
    Integer(ap), Allocatable :: axeF(:)
    Integer(ap), Allocatable :: prev_nodea(:)
    Integer(ap), Allocatable :: axenf(:)
    Real(dp) :: origin(3)=(/ 0.0,0.0,0.0 /)
    Real(dp), Allocatable :: xrec(:)
    Real(dp), Allocatable :: yrec(:)
    Real(dp), Allocatable :: zrec(:)
    Real(dp), Allocatable :: xa(:)
    Real(dp), Allocatable :: ya(:)
    Real(dp), Allocatable :: za(:)
    Real(dp), Allocatable :: diamrec(:)
    Real(dp), Allocatable :: diama(:)
    Real(dp), Allocatable :: agea(:)
    Real(dp) :: treal

    !initialize root function
    Call init_RootType1(origin)
    Print *,'RootTyp is running ...'
    Do simtime=1,t,1
       Call iterate_roottype1(simtime)
    End Do
    Print *,'RootTyp converged'
    Call number_of_nodes(n_nodes, n_br_crea,n_br_del) !C-function
    n_meris=n_br_crea-n_br_del
    !n_meris= all meristems
    !n_nodes= number maximum of nodes: some of them have been deleted...
    !print *,n_br_crea,n_br_del,n_meris,n_nodes
    !print *,'------------------------------'
    Allocate (node(n_nodes))
    Allocate (prev_node(n_nodes))
    Allocate (ordn(n_nodes))
    Allocate (xrec(n_nodes))
    Allocate (yrec(n_nodes))
    Allocate (zrec(n_nodes))
    Allocate (diamrec(n_nodes))
    Allocate (agerec(n_nodes)) !in days
    Allocate (axenf(n_nodes))
    Allocate (axeF(n_meris))
    Allocate (orda(n_meris))
    Allocate (xa(n_meris))
    Allocate (ya(n_meris))
    Allocate (za(n_meris))
    Allocate (diama(n_meris))
    Allocate (agea(n_meris))
    Allocate (prev_nodea(n_meris))

    ! extract nodes C-function < interface
    Call extract_nodes(node,prev_node,ordn,xrec,yrec,zrec,diamrec,agerec,axeF,orda,xa,&
         & ya,za,diama,agea,prev_nodea,n_axes,axenf,n_nodecorr) 
    nrec(ipl)=n_nodecorr

    !meaning of node is not clear?
    xs(1:nrec(ipl),ipl)=xrec(1:nrec(ipl))
    ys(1:nrec(ipl),ipl)=yrec(1:nrec(ipl))
    zs(1:nrec(ipl),ipl)=-zrec(1:nrec(ipl))!roottip in mm
    !Z must decrease downward!!
    Do n=2,nrec(ipl)
       If (zs(n,ipl)>0.0_dp) zs(n,ipl)=-zs(n,ipl) !in case there is roots >soil surface
       If ((n.Ne.1).And.(axenf(n).Ne.(axenf(n-1)))) Then
          If (axenf(n).Ne.(axenf(n-1)+1)) Then
             Print *,'axe',axenf(n-1)+1,'is missing'
          Endif
       Endif
    Enddo

    ! create R-SWMS variables for root
    irecpr(1:nrec(ipl),ipl)=prev_node(1:nrec(ipl))
    ordseg(1:nrec(ipl),ipl)=ordn(1:nrec(ipl))+1 !in Roottip from 0 to 7, here from 1 to 8
    ibrseg(1:nrec(ipl),ipl)=axenf(1:nrec(ipl))
    timorg(1:nrec(ipl),ipl)=agerec(1:nrec(ipl)) !day of creation+time real=age
    segdiam(1:nrec(ipl))=diamrec(1:nrec(ipl))

    !growing apices= "meristems with axes" in RootTyp
    ngrow=n_axes
    nbr=n_axes !number of apex=number of branhces
    xg(1:n_axes,ipl)=xa(1:n_axes)
    yg(1:n_axes,ipl)=ya(1:n_axes)
    zg(1:n_axes,ipl)=-za(1:n_axes)
    irecsg(1:n_axes,ipl)=prev_nodea(1:n_axes)
    ordgrw(1:n_axes,ipl)=orda(1:n_axes)+1
    ibrgrw(1:n_axes,ipl)=axeF(1:n_axes)
    Allocate(connex(1:nrec(ipl)))

    !check root nodes outside of the soil
    If (zg(1,ipl)>0.0_dp) zg(1,ipl)=-zg(1,ipl)
    Do n=2,ngrow(ipl)
       If (zg(n,ipl)>0.0_dp) zg(n,ipl)=-zg(n,ipl) !in case there is roots >soil surface
       If (axenf(n).Ne.(axenf(n-1))) Then
          If (axenf(n).Ne.(axenf(n-1)+1)) Then
             Print *,'axe',axenf(n-1)+1,'is missing'
          Endif
       Endif
    Enddo

    !check root type (Couvreur nov 2010)
    If (Maxval(ordseg).Eq.13) Then!Recognizes old roots created by RootTyp
       maizeroottyp=.True.
    Elseif (Maxval(ordseg).Eq.5) Then
       loliumroottyp=.True.
    Elseif (Maxval(ordseg).Eq.20) Then
       wheatroottyp=.True.
    Endif

    !estimate length and surface of segments and get naxes
    Call Estimate_Seglen(ipl)
    !   CALL OutRoo(real(0),naxes,0,0,0,0,0)
    !simplifiy system
    Do n=1,2
       Call SimplifyRoot(ipl)
       Call Estimate_Seglen2(ipl)
    Enddo
    !   CALL OutRoo(real(999),naxes,0,0,0,0,0)
    !   CALL AdaptOutRoot
    !   CALL Estimate_Seglen2
    Call close_c_interface() !C-function
    Call finish_RootType1() !C-function
    treal=t
    Print *,'number of nodes after',n,' iterations at time t= ', treal,' is ',nrec(ipl)
    Call OutRoo(treal,0,ipl)!kOut is 0
    If (.Not.continu) Call CheckSize(1)!No need to check size if continuous domain. ipl is 1 because RootTyp can currently only simulate 1 root system (Couvreur feb 2010)
  End Subroutine RunRootTip
  !********************************************************************************
  Subroutine Estimate_Seglen(ipl)
    Use typedef
    Use ParamData
    Use RootData
    Implicit None

    Real(dp) :: xend,yend,zend,xtop,ytop,ztop
    Integer(ap) :: inode,n,inode2,num
    INTEGER(ap),INTENT(in)::ipl
    
    If(lrrt) connex(:)=.False.
    inode=1
    naxes(ipl)=0
    Do n=1,nbr(ipl) !for all branches
       brlgth(n,ipl)=0.0_dp
       inode=irecsg(n,ipl) !inode =prec
       xtop=xs(inode,ipl)!end note=apex
       ytop=ys(inode,ipl)
       ztop=zs(inode,ipl)
       xend=xg(n,ipl)!Length associated with a node is the length of the segment under the node (Couvreur feb 2010)
       yend=yg(n,ipl)
       zend=zg(n,ipl)
       seglen(inode,ipl)=Sqrt((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
       If (seglen(inode,ipl)==0) Then!apical node "on" the meristem
          xs(inode,ipl)=xs(irecpr(inode,ipl),ipl)+(xtop-xs(irecpr(inode,ipl),ipl))/2.0_dp
          ys(inode,ipl)=ys(irecpr(inode,ipl),ipl)+(ytop-ys(irecpr(inode,ipl),ipl))/2.0_dp
          zs(inode,ipl)=zs(irecpr(inode,ipl),ipl)+(ztop-zs(irecpr(inode,ipl),ipl))/2.0_dp
          xtop=xs(inode,ipl)
          ytop=ys(inode,ipl)
          ztop=zs(inode,ipl)
          seglen(inode,ipl)=Sqrt((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
       Endif
       brlgth(n,ipl)=brlgth(n,ipl)+seglen(inode,ipl)
       num=1
       If (lrrt) Then
          segsur(inode,ipl)=seglen(inode,ipl)*pi*segdiam(inode)!segment surface
       Else 
          segsur(inode,ipl)=seglen(inode,ipl)*pi*2._dp*segrad(inode,ipl) 
       End If
       Do While (ibrseg(irecpr(inode,ipl),ipl)==n.And.irecpr(inode,ipl).Gt.0)
          inode2=irecpr(inode,ipl)
          xend=xtop
          yend=ytop
          zend=ztop
          xtop=xs(inode2,ipl)
          ytop=ys(inode2,ipl)
          ztop=zs(inode2,ipl)
          seglen(inode2,ipl)=Sqrt((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
          If (seglen(inode2,ipl)==0) Then!superposed nodes
             xs(inode2,ipl)=xs(irecpr(inode2,ipl),ipl)+(xtop-xs(irecpr(inode2,ipl),ipl))/2.0_dp
             ys(inode2,ipl)=ys(irecpr(inode2,ipl),ipl)+(ytop-ys(irecpr(inode2,ipl),ipl))/2.0_dp
             zs(inode2,ipl)=zs(irecpr(inode2,ipl),ipl)+(ztop-zs(irecpr(inode2,ipl),ipl))/2.0_dp
             xtop=xs(inode2,ipl)
             ytop=ys(inode2,ipl)
             ztop=zs(inode2,ipl)
             seglen(inode2,ipl)=Sqrt((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
          Endif
          If (lrrt) Then
             segsur(inode2,ipl)=seglen(inode2,ipl)*pi*segdiam(inode2)!segment surface
          Else
             segsur(inode2,ipl)=seglen(inode2,ipl)*pi*2._dp*segrad(inode2,ipl)
          End If
          brlgth(n,ipl)=brlgth(n,ipl)+seglen(inode2,ipl)
          num=num+1;
          inode=inode2
       Enddo
       br_rec(n)=irecpr(inode,ipl)!records at which node (on the father axis) this axis is connected
       num_seg(n,ipl)=num
       If (maizeroottyp) Then
          If (ordgrw(n,ipl).Le.11) naxes(ipl)=naxes(ipl)+1 !number of axes ("principal roots")
       Elseif (loliumroottyp) Then
          If (ordgrw(n,ipl).Lt.3) naxes(ipl)=naxes(ipl)+1
       Elseif (wheatroottyp) Then
          If (ordgrw(n,ipl).Le.18) naxes(ipl)=naxes(ipl)+1
       Elseif (ordgrw(n,ipl).Lt.3) Then!Other roots (Couvreur nov 2010)
          naxes=naxes+1 !number of axes ("principal roots")
       Endif
       If (lrrt .And. (br_rec(n).Ne.0)) connex(br_rec(n))=.True. !logical which defines whether there is a connection to that node
    Enddo
    Return
  End Subroutine Estimate_Seglen
  !****************************************************************************
  Subroutine SimplifyRoot(ipl)!Redesigned for the new RootTyp with senescence (Couvreur may 2010)
    Use typedef
    Use RootData, Only : irecpr,nrec,ibrseg,irecsg,connex,xs,ys,zs,timorg,segdiam,ordseg,seglen,nplant
    Use GridData, Only : dxGrid,dyGrid,dzGrid
    Implicit None

    Real(dp) :: xrec(1:nrec(nplant)),yrec(1:nrec(nplant)),zrec(1:nrec(nplant)),agerec(1:nrec(nplant))
    Real(dp) :: diamrec(1:nrec(nplant)),ltot,resolution
    Integer(ap) :: ordn(1:nrec(nplant)),axenf(1:nrec(nplant))
    Integer(ap) :: inew,iold,iprec_old,oldi(nrec(nplant))
    Integer(ap) :: old2new(nrec(nplant)),iprec(1:nrec(nplant))
    INTEGER(ap),INTENT(in) :: ipl

    !initialisation
    xrec=xs(1:nrec(ipl),ipl)
    yrec=ys(1:nrec(ipl),ipl)
    zrec=zs(1:nrec(ipl),ipl)
    axenf=ibrseg(1:nrec(ipl),ipl)
    agerec=timorg(1:nrec(ipl),ipl)
    diamrec=segdiam(1:nrec(ipl))
    ordn=ordseg(1:nrec(ipl),ipl)
    iprec=irecpr(1:nrec(ipl),ipl)
    resolution=Min(2.0_dp,dxGrid,dyGrid,dzGrid) !Adapt root simplification to grid resolution (threshold = 2 cm because root properties evolve with distance of that order) (Couvreur mar 2010)
    inew=1
    iold=1!start with node 2==proximal node (Couvreur may 2010)
    old2new(1)=1
    !check each root
    Do While (iold<nrec(ipl))
       inew=inew+1
       iold=iold+1
       iprec_old=iprec(iold)!node before
       irecpr(inew,ipl)=old2new(iprec(iold))
       ltot=seglen(iprec_old,ipl) !length
       Do While ((axenf(iprec_old)==axenf(iold)).And.(connex(iold).Eqv..False.)&
            &.And.(seglen(iold,ipl)<0.8*resolution).And.(seglen(iold,ipl)+ltot.Le.resolution)&
            &.And.((iold).Ne.irecsg(ibrseg(iold,ipl),ipl))) !this is not the apical node
          !same branches + no connection to iold
          !(one more condition : the new segment is shorter than the grid resolution) (Couvreur feb 2010)
          old2new(iold)=999999
          ltot=ltot+seglen(iold,ipl)
          iold=iold+1
          iprec_old=iprec(iold)!node before
       Enddo
       ibrseg(inew,ipl)=axenf(iold) !br#
       xs(inew,ipl)=xrec(iold)
       ys(inew,ipl)=yrec(iold)
       zs(inew,ipl)=zrec(iold)
       ordseg(inew,ipl)=ordn(iold)
       timorg(inew,ipl)=agerec(iold) !orig time
       segdiam(inew)=diamrec(iold) !diam
       oldi(inew)=iold !keep in mind the previous numerotation
       old2new(iold)=inew
       If (irecsg(ibrseg(inew,ipl),ipl)==oldi(inew)) Then
          ! if the previous node of the apex of the branch where node inew is is inew, then itmust be also updated
          irecsg(ibrseg(inew,ipl),ipl)=inew
       Endif
    Enddo
    nrec(ipl)=inew
  End Subroutine SimplifyRoot
  !****************************************************************************
  !adapt root description to RSWMS
  Subroutine AdaptOutRoot(ipl)
    Use TypeDef
    Use RootData, Only : irecpr,nrec,ibrseg,irecsg,xs,ys,zs,timorg,segdiam,num_seg,nbr,nplant
    Implicit None

    Real(dp):: xrec(1:nrec(nplant)),yrec(1:nrec(nplant)),zrec(1:nrec(nplant)),agerec(1:nrec(nplant))
    Real(dp) :: diamrec(1:nrec(nplant))
    Integer(ap):: prev(1:nrec(nplant)),axenf(1:nrec(nplant)),prec(1:nrec(nplant))
    Integer(ap) :: ibr2,ibr1,nn,i,n, old2new(nrec(nplant))
    INTEGER(ap),INTENT(in) :: ipl

    !initialisation
    xrec=xs(1:nrec(ipl),ipl)
    yrec=ys(1:nrec(ipl),ipl)
    zrec=zs(1:nrec(ipl),ipl)
    axenf=ibrseg(1:nrec(ipl),ipl)
    agerec=timorg(1:nrec(ipl),ipl)
    prev=irecpr(1:nrec(ipl),ipl)
    diamrec=segdiam(1:nrec(ipl))
    ibr1=1
    Do n=1,nbr(ipl) !for all branches
       !get the number of nodes for that branch
       nn=num_seg(n,ipl)
       ibr1=ibr1 !ibr1=ibr2+1 at the run  of teh next loop!!
       ibr2=ibr1+nn-1
       !adapt previous node to the meristem
       irecsg(n,ipl)=ibr2
       Do i=1,nn
          xs(ibr2,ipl)=xrec(ibr1)
          ys(ibr2,ipl)=yrec(ibr1)
          zs(ibr2,ipl)=zrec(ibr1)
          ibrseg(ibr2,ipl)=axenf(ibr1)
          timorg(ibr2,ipl)=agerec(ibr1) !orig time
          segdiam(ibr2)=diamrec(ibr1) !diam
          prec(ibr2)=prev(ibr1)
          old2new(ibr1)=ibr2
          ibr1=ibr1+1
          ibr2=ibr2-1
       Enddo
    Enddo

    !correct prec matrix
    Do i=1,nrec(ipl)
       If (prec(i).Ne.0_dp) Then
          irecpr(i,ipl)=old2new(prec(i))
       Else
          irecpr(i,ipl)=0
       Endif
    Enddo
  End Subroutine AdaptOutRoot
  !*********************************************************************
  Subroutine Estimate_Seglen2(ipl)
    Use typedef
    Use ParamData
    Use RootData, Only : seglen,segsur,brlgth,segdiam,irecpr,nbr, &
         ibrseg,zs,ys,xs,irecsg,num_seg,connex,zg,yg,xg,naxes
    Implicit None

    Real(dp) :: xend,yend,zend,xtop,ytop,ztop
    Integer(ap) :: inode,n,brn,num,connected2
    INTEGER(ap),INTENT(in)::ipl

    connex(:)=.False.
    inode=1
    naxes=0
    Do n=1,nbr(ipl) !for all branches
       inode=irecsg(n,ipl) !inode =apex<-meristem
       xend=xg(n,ipl)!end note=apex
       yend=yg(n,ipl)
       zend=zg(n,ipl)
       brn=n !branch number
       brlgth(brn,ipl)=0.0_dp
       num=0
       Do While (brn==n.And.inode.Gt.0)
          xtop=xs(inode,ipl)
          ytop=ys(inode,ipl)
          ztop=zs(inode,ipl)
          seglen(inode,ipl)=Sqrt((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
          !length of node is the length of the segment located on the top of it!Wrong, at the basis (Couvreur feb 2010)
          !length of node 1=0
          If (seglen(inode,ipl)==0) Then!superposed nodes
             xs(inode,ipl)=xs(irecpr(inode,ipl),ipl)+(xtop-xs(irecpr(inode,ipl),ipl))/2.0_dp
             ys(inode,ipl)=ys(irecpr(inode,ipl),ipl)+(ytop-ys(irecpr(inode,ipl),ipl))/2.0_dp
             zs(inode,ipl)=zs(irecpr(inode,ipl),ipl)+(ztop-zs(irecpr(inode,ipl),ipl))/2.0_dp
             xtop=xs(inode,ipl)
             ytop=ys(inode,ipl)
             ztop=zs(inode,ipl)
             seglen(inode,ipl)=Sqrt((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
          Endif
          segsur(inode,ipl)=seglen(inode,ipl)*pi*segdiam(inode)!segment surface
          brlgth(brn,ipl)=brlgth(brn,ipl)+seglen(inode,ipl)
          num=num+1
          xend=xtop
          yend=ytop
          zend=ztop
          inode=irecpr(inode,ipl)
          connected2=inode
          brn=ibrseg(inode,ipl)
       Enddo
       num_seg(n,ipl)=num
       If (connected2.Ne.0) connex(connected2)=.True.
    Enddo
  End Subroutine Estimate_Seglen2
  !****************************************************************************
  Subroutine CheckSize(ipl)!Checksize has to know wich plant is being checked in order to place it correctly with report to the grid (Couvreur dec 2009)
    ! check max/min position of roots as compared to the soil grid
    Use RootData
    Use GridData
    Use MPIutils, Only: stop_program
    Implicit None

    Real(dp)::maxX,maxY,maxZ,maxXs,maxYs,maxZs,minX,minY,minZ,minXs,minYs,minZs
    Integer(ap)::ipl

    maxX=Minval(xGrid)+nex*dxgrid !Adapted to continuous and non continuous soil domain (Couvreur dec 2009)
    minX=Minval(xGrid)
    maxY=Minval(YGrid)+ney*dygrid
    minY=Minval(YGrid)
    maxZ=Maxval(ZGrid)
    minZ=Minval(ZGrid)
    maxXs=Maxval(xs(1:nrec(ipl),ipl))+xplant(ipl)
    minXs=Minval(xs(1:nrec(ipl),ipl))+xplant(ipl)
    maxYs=Maxval(Ys(1:nrec(ipl),ipl))+yplant(ipl)
    minYs=Minval(Ys(1:nrec(ipl),ipl))+yplant(ipl)
    maxZs=Maxval(Zs(1:nrec(ipl),ipl))
    minZs=Minval(Zs(1:nrec(ipl),ipl))
    If (maxXs>maxX) Then
       Print *,'X root too large'
       Goto 20
    Endif
    If (maxYs>maxY) Then
       Print *,'Y root too large'
       Goto 20
    Endif
    If (maxZs>maxZ) Then !
       Print *,'Upper root node (=',maxZs,') is higher than soil max. z (=',maxZ,')'
       Goto 20
    Endif
    If (minXs<minX) Then
       Print *,'X root too small'
       Goto 20
    Endif
    If (minYs<minY) Then
       Print *,'Y root too small'
       Goto 20
    Endif
    If (minZs<minZ) Then
       Print *,'Lower root node (=',minZs,') is deeper than soil min. z (=',minZ,')'
       Goto 20
    Endif
    Return
20  Print *,'root max:',maxXs,maxYs,maxZs,'root min:',minXs,minYs,minZs,'soil:',maxX,maxY
    Call stop_program('Please re-run R-SWMS/RootTyp')
  End Subroutine CheckSize
  !****************************************************************************

End Module RootTyp
