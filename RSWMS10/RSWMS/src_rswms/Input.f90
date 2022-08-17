!> \file Input.f90
!! \brief Module for Input.

!> Module Input
Module Input

Contains

!> ### reads main inputs ###
  Subroutine Applic(dt)
    !> \param nPt = total number of nodal points
    !> \param nBCpts = total number of nodal points with specified BC
    !> \param nElm = total number of elements
    !> \param ne* = number of elements (half-cuboids) in '*'-direction
    !> \param nel = nex*ney (number of elements per horizontal element-layer)
    Use iso_c_binding
    Use Typedef
    Use ParamData, Only: ldirect,lvtk,lOutPartrace,lChem,lretry,last_out,&
         lSalinity,maxmat,maxIrrig,maxbdr,mxBcCh,lPartUp, lClimate
    Use RootData, Only: lno_root_growth,lno_RWU,ltemp,lCalloc,lrrs,lrrt,&
         lRootTyp_growth,lSomma_growth,lRootBox_growth,lUpdate_growth,lCou,lDou,lno_Archi,&
         lFed,lSinkCube,ltwo_grids,lSign_new,nplant
    Use GridData
    Use GridData2
    Use tmctrl
    Use BoundData
    Use SolData, Only: hold, htemp, hnew, conc, Kode, Kcell, matNum, matNum2,nmat,KodCB,lTab,lCelia,lRoot_explicit
    Use CumData
    Use DomData
    Use ObsData
    Use doussanmat, Only : old,oldT,ana_aan
    Use RhizoData, Only: lRhizo
    Use SoluteRootMat, Only: uptakeorder,decayrate,l_degrad
    Use EnviData
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: nh,nh_tot,nQ,nI,nFrDr,l,j,ise,ie,ip,i,idum,idum2,k,err,ii,nl,nl2,aa,bb,dummy
    Integer(ap) :: i2,assCode,hfun,nQ1,decayorder
    Integer(ap),dimension(8) :: values
    Integer(ap) :: xIrrig(maxIrrig),yIrrig(maxIrrig),zIrrig(maxIrrig)
    Integer(ap), Allocatable, Dimension(:) :: node_temp
    Real(dp) :: Qbcrec(mxBcCh)
    Real(dp) :: A11,A22,A33,A12,A13,A23,C11,C22,C33,ini,ixmin,ixymin,imin,dt,xqmin,xqmax,xhmin,xhmax,xminh,xmaxh,zminh,zmaxh
    Logical :: lOrt=.False.,firstOK=.False.,ltop=.True.,lFrDr=.False.
    Character text*5
    Character infile*14
    Character date*8,time*10,zone*5

    Allocate(node_temp(nPt))
    nMat=0
    xmax=-1.E+30_dp
    xmin=+1.E+30_dp
    ymax=-1.E+30_dp
    ymin=+1.E+30_dp
    zmax=-1.E+30_dp
    zmin=+1.E+30_dp

    IF (ltwo_grids) THEN
    xmax2=-1.E+30_dp
    xmin2=+1.E+30_dp
    ymax2=-1.E+30_dp
    ymin2=+1.E+30_dp
    zmax2=-1.E+30_dp
    zmin2=+1.E+30_dp
    END IF 

    profOK=.False.
    Open(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')

!> #### reads control.in ####
    Open (Unit=10,FILE='in/control.in',STATUS='OLD',ERR=10)
    Write(15,'(//''++ General Information ++'')',advance='no')
    Write(15,'(/''-------------------------'')',advance='no')
    CALL DATE_AND_TIME(DATE, TIME, ZONE, VALUES)
    Write(15,'(/''Run with RSWMS10 on '',a,''/'',a,''/'',a,'' at '', a,''hr'',a,''.'')',advance='no')date(7:8),date(5:6),date(1:4),time(1:2),time(3:4) 
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) LnUnit,TmUnit,MsUnit,CnUnit
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) itMax,itMaxRoot
    Read (10,*)
    Read (10,*)
    Read (10,*) RelEps,factorRelEps
    Read (10,*)
    Read (10,*)
    Read (10,*) epslonPH,epslonWC,epslonR,epslonS
    If (factorRelEps.Eq.0) Then
       Call stop_program('Control.in: factor for relative criterium may not be zero')
    Endif
    If (RelEps)   Write(15,'(/''* Relative tolerance is used for WC, PH and Sink with value of'',f5.3)',advance='no')1/factorRelEps

    Read (10,*)
    Read (10,*)
    Read (10,*) dt,dtMin,dtMax,FacInc,FacDec,dtRoot
    Read (10,*)
    Read (10,*)
    Read (10,*) lretry,last_out,nplant
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) nOut
    If (nOut.Gt.mxOut) Call stop_program('Input from  < control.in >  exceeds parameter "mxOut". Program terminated.')
    Read (10,*)
    Read (10,*)
    Read (10,*) (tOut(i),i=1,nOut)
    tmax=tOut(nOut)
    Read (10,*)
    Read (10,*)
    Read (10,*) lvtk,lOutpartrace,profOK
    Read (10,*)
    Read (10,*)
    Read (10,*) dtprof
    If(profOK) Then
       If (dtprof<999) Then
          nouProf=0
          ini=0
          Do  While (ini.Lt.tmax)
             nouProf=nouProf+1
             tOuProf(nouProf)=ini+dtprof
             ini=tOuProf(nouProf)
             If (nouProf.Gt.mxProf) Then
                Print *,'too small time step in z-profiles: only ',mxProf,' profiles will be kept'
                Goto 77
             Endif
          End Do
       End If
    End If
77  Read(10,*)
    Read(10,*)
    Read(10,*)
    Read(10,*)
! Root water uptake modeling
    Read(10,*) lno_RWU,lFed,lDou,lCou,lSinkCube !
    Read(10,*)
    Read(10,*)
    Read(10,*)
    Read(10,*) lno_Archi,lrrs,lrrt
    Read(10,*)
    Read(10,*)
    Read(10,*)
! Root growth
    Read(10,*) lno_root_growth,lRootTyp_growth,lSomma_growth,lRootBox_growth,lUpdate_growth
    Read(10,*)
    Read(10,*)
    Read(10,*)
    Read(10,*) lCalloc,lChem,ltemp,continu,lSalinity,lPartUp, lRhizo, lclimate, ltwo_grids
    Read(10,*)
    Read(10,*)
    Read(10,*)
    Read(10,*) ldirect,old,ana_aan,ltab,lCelia,lRoot_explicit
    Close (10)

    oldT=old

!> check for compatibility of processes
    If((.Not.lno_RWU) .And. lno_archi .And. .Not.lCou)  Call stop_program('RWU models should always be associated with a root architecture')
    If(.Not.lCou .And. .Not.lSinkCube)  Call stop_program('Sink nodes currently only programmed to generate SSF for Hydrus -> uses Couvreur & archi from RootSys')
    If(.Not.lno_RWU .And. .Not.lFed .And. .Not.lDou .And. .Not.lCou) Call stop_program('Please set at least one RWU option in control.in to .TRUE.')
    If(lRootTyp_growth) Call stop_program('RootTyp root growth during the scenario run is not yet implemented')
    If(lCalloc .And. .Not. ((lRootBox_growth.OR.lSomma_growth).And.lDou)) Call stop_program('Assimilate allocation can currently only be associated with Somma/RootBox root growth and Doussan root water uptake')
    If(continu) Write(15,'(/''* XY domain periodicity for soil water flow, root system architecture and root water flow'')',advance='no')
    If (continu.And.(.Not.old)) Call stop_program('Domain continuity currently only works with soil Averaging method')
    If(ldirect .And. .Not. lDou) Call stop_program('Direct Doussan can only be used with Doussan RWU model')
    If(ldirect .And. .Not. lno_root_growth) Call stop_program('Direct Doussan can only be used without root growth')
    If(ana_aan .And. .Not. old) Call stop_program('Cannot use microscopic model and memory reduction simultaneous')
    If(.Not. lno_root_growth.and. .Not. (lRootBox_growth.OR.lSomma_growth.OR.lRootTyp_growth.OR.lUpdate_growth)) Call stop_program('You should define a growth model')!lRootTyp_growth	lSomma_growth	lRootBox_growth	lUpdate_growth

!> #### reads grid.in ####
    Open (Unit=10,FILE='in/grid.in',STATUS='OLD',ERR=20)
    Write(15,'(//''++ Soil Geometry ++'')',advance='no')
    Write(15,'(/''--------------------'')',advance='no')
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) nPt,nElm
    Allocate(hOld(nPt))
    Allocate(hTemp(nPt))
    Allocate(hNew(nPt))
    Allocate(conc(nPt))
    Allocate(Kode(nPt))
    Allocate(Kcell(nPt))
    Allocate(MatNum(nPt))
    If (nPt.Gt.maxnod) Call stop_program('Input from  < grid.in >  exceeds parameter "maxnod". Program terminated.')
    If (nElm.Gt.maxelm) Call stop_program('Input from  < grid.in >  exceeds parameter "maxelm". Program terminated.')
    lOrt=.True.
    Read (10,*)
    Read (10,*)
    Read (10,*) nx,ny,nz,nex,ney,nez,dxGrid,dyGrid,dzGrid
    nel=nex*ney
    nl=nx*ny
    Call IniGrid !> initialize xgrid, ygrid,zgrid and scaling factors
    ! elements:
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Do i=1,nElm
       Read(10,*) iDum,(elmnod(k,i),k=1,8),subN(i),A11,A22,A33,A12,A13,A23,C11,C22,C33
       ConAxx(i)=C11*A11*A11+C22*A12*A12+C33*A13*A13
       ConAyy(i)=C11*A12*A12+C22*A22*A22+C33*A23*A23
       ConAzz(i)=C11*A13*A13+C22*A23*A23+C33*A33*A33
       ConAxy(i)=C11*A11*A12+C22*A12*A22+C33*A13*A23
       ConAxz(i)=C11*A11*A13+C22*A12*A23+C33*A13*A33
       ConAyz(i)=C11*A12*A13+C22*A22*A23+C33*A23*A33
    End Do
    Write(15,'(/''* Number of elements in x, y and z directions:'',3i5)',advance='no')nex,ney,nez
    Write(15,'(/''* Number of nodes in x, y and z directions:'',3i5)',advance='no')nx,ny,nz
    Write(15,'(/''* Length of elements in x, y and z directions:'',3f7.2)',advance='no')dxGrid,dyGrid,dzGrid
    Close (10)

    IF (ltwo_grids) THEN
!> #### reads grid2.in ####
    Open (Unit=10,FILE='in/grid2.in',STATUS='OLD',ERR=20)
    Write(15,'(//''++ Soil Geometry of second grid ++'')',advance='no')
    Write(15,'(/''--------------------'')',advance='no')
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) nPt2,nElm2
    Allocate(MatNum2(nPt2))
    If (nPt2.Gt.maxnod) Call stop_program('Input from  < grid2.in >  exceeds parameter "maxnod". Program terminated.')
    If (nElm2.Gt.maxelm) Call stop_program('Input from  < grid2.in >  exceeds parameter "maxelm". Program terminated.')
    lOrt=.True.
    Read (10,*)
    Read (10,*)
    Read (10,*) nx2,ny2,nz2,nex2,ney2,nez2,dxGrid2,dyGrid2,dzGrid2
    nel2=nex2*ney2
    nl2=nx2*ny2
    IF ((nex*dxGrid.NE.nex2*dxGrid2).OR.(nex*dxGrid.NE.nex2*dxGrid2).OR.(nex*dxGrid.NE.nex2*dxGrid2)) THEN
        Call stop_program('total size of grid 1 and grid 2 must be in agreement') 
    ENDIF 
    Call IniGrid2 !> initialize xgrid, ygrid,zgrid and scaling factors
    ! elements:
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Do i=1,nElm2
       Read(10,*) iDum2,(elmnod2(k,i),k=1,8) 
    End Do
    Write(15,'(/''* grid two, Number of elements in x, y and z directions:'',3i5)',advance='no')nex2,ney2,nez2
    Write(15,'(/''* grid two, Number of nodes in x, y and z directions:'',3i5)',advance='no')nx2,ny2,nz2
    Write(15,'(/''* grid two, Length of elements in x, y and z directions:'',3f7.2)',advance='no')dxGrid2,dyGrid2,dzGrid2
    END IF 
    Close (10)
	
!> #### read nodes.in ####
    If (.Not.lretry) Then
       Open (Unit=10,FILE='in/nodes.in',STATUS='OLD',ERR=30)
3      Read (10,'(A5)') text
       If (text.Ne.'Node#') Goto 3
       Do i=1,nPt
          Read (Unit=10,iostat=err) iDum,MatNum(i),xGrid(i),yGrid(i), zGrid(i),hOld(i),Conc(i),Axy(i),Bxy(i),Dxy(i),Exy(i)
          If (err.Ne.0) Then ! in case there are no scaling factors
             Read (10,*) iDum,MatNum(i),xGrid(i),yGrid(i),zGrid(i), hOld(i),Conc(i)
             Axy(i)=1.
             Bxy(i)=1.
             Dxy(i)=1.
             Exy(i)=1.
          Endif
          Kode(i)=0
          hNew(i)=hOld(i)
          hTemp(i)=hOld(i)
          nMat=Max(nMat,MatNum(i))
          xmin=Min(xmin,xGrid(i))
          ymin=Min(ymin,yGrid(i))
          zmin=Min(zmin,zGrid(i))
          xmax=Max(xmax,xGrid(i))
          ymax=Max(ymax,yGrid(i))
          zmax=Max(zmax,zGrid(i))
       End Do
    Else
       If (last_out.Gt.0.And.last_out.Lt.10) Then
          Write(infile,'(A11,I1)') 'out/outfem.',last_out
       Elseif (last_out.Lt.100) Then
          Write(infile,'(A11,I2)') 'out/outfem.',last_out
       Elseif (last_out.Lt.1000) Then
          Write(infile,'(A11,I3)') 'out/outfem.',last_out
       Else
          Call stop_program('Last output # should be between 1 and 999, when retrying a simulation from an existing output')
       Endif
       Open (Unit=10,FILE=infile,STATUS='OLD',ERR=30)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Do i=1,nPt ! in case there are no scaling factors
          Read (10,*) iDum,MatNum(i),xGrid(i),yGrid(i),zGrid(i), hOld(i),Conc(i)
          Axy(i)=1.
          Bxy(i)=1.
          Dxy(i)=1.
          Exy(i)=1.
          Kode(i)=0
          hNew(i)=hOld(i)
          hTemp(i)=hOld(i)
          nMat=Max(nMat,MatNum(i))
          xmin=Min(xmin,xGrid(i))
          ymin=Min(ymin,yGrid(i))
          zmin=Min(zmin,zGrid(i))
          xmax=Max(xmax,xGrid(i))
          ymax=Max(ymax,yGrid(i))
          zmax=Max(zmax,zGrid(i))
       End Do
    Endif
    Do i=1,nx
       xCol(i)=xmin+(i-1)*dxGrid
    Enddo
    Do i=1,ny
       yCol(i)=ymin+(i-1)*dyGrid
    Enddo
    Close (10)
    If (nmat.Gt.maxmat) Call stop_program('Input from nodes file exceeds parameter "maxmat". Program terminated.')
    Do iP=1,nPt
       ListNE(iP)=0
    End Do
    nBand=1
    !> nBand calculated like in SWMS3D
    Do  iE=1,nElm
       Do  iSE=1,5
          i=elmnod(iL(1,iSE,subN(iE)),iE)
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          ListNE(i)=ListNE(i)+1
          ListNE(j)=ListNE(j)+1
          ListNE(k)=ListNE(k)+1
          ListNE(l)=ListNE(l)+1
          If(Abs(i-j).Gt.nBand) nBand=Abs(i-j)
          If(Abs(i-k).Gt.nBand) nBand=Abs(i-k)
          If(Abs(i-l).Gt.nBand) nBand=Abs(i-l)
          If(Abs(j-k).Gt.nBand) nBand=Abs(j-k)
          If(Abs(j-l).Gt.nBand) nBand=Abs(j-l)
          If(Abs(k-l).Gt.nBand) nBand=Abs(k-l)
       End Do
    End Do
    nBand=nBand+1
	
!> #### reads nodes2.in ####
    If (ltwo_grids) Then
       Open (Unit=10,FILE='in/nodes2.in',STATUS='OLD',ERR=30)
399      Read (10,'(A5)') text
       If (text.Ne.'Node#') Goto 399
       Do i=1,nPt2
          Read (Unit=10,iostat=err) iDum2,MatNum2(i),xGrid2(i),yGrid2(i), zGrid2(i)
          If (err.Ne.0) Then ! in case there are no scaling factors
             Read (10,*) iDum2,MatNum2(i),xGrid2(i),yGrid2(i),zGrid2(i)
          Endif
          xmin2=Min(xmin2,xGrid2(i))
          ymin2=Min(ymin2,yGrid2(i))
          zmin2=Min(zmin2,zGrid2(i))
          xmax2=Max(xmax2,xGrid2(i))
          ymax2=Max(ymax2,yGrid2(i))
          zmax2=Max(zmax2,zGrid2(i))
       End Do
       If (nmat.Gt.maxmat) Call stop_program('Input from nodes2 file exceeds parameter "maxmat". Program terminated.')
    ENDIF
    Close (10)
    
	
!> #### reads mesh.in for additional geometrical information ####
    Open(Unit=10,FILE='in/mesh.in',STATUS='OLD',ERR=60)
    Do i=1,9
       Read (10,*)
    End Do
    Read (10,*) geom
    If (geom.Eq.3) Then
       Read (10,*)
       Read (10,*) rad_cyl
       ! grid center
       x_cent = (nex*dxGrid/2+xmin)
       y_cent = (ney*dyGrid/2+ymin)
    End If
    Close(10)

    If(lclimate)Then
       Call ClimateIN
    Endif

!> #### reads bc.in; soil boundary conditions #### ! Si 2 materiaux (puis voir water.f90)
  IF ((geom.Eq.4) .And.(nMat.Eq.2)) THEN! rhizotron case
    Open (Unit=10,FILE='in/bcRhizotron.in',STATUS='OLD',ERR=40)
    Write(15,'(//''++ Soil Boundary Conditions ++'')',advance='no')
    Write(15,'(/''------------------------------'')',advance='no')
    Read(10,*)
    Read(10,*)
    Read(10,*)
       Write(15,'(/''* Rhizotron specific BC'')',advance='no')
    nBCPts=0
    ii=1
    nh_tot=0
    Do l=1,9 
      Read(10,*) dummy,xminh,xmaxh,zminh,zmaxh
       If ((xminh.Lt.xmin).Or.(xmaxh.Gt.xmax).Or.(zminh.Lt.zmin).Or.(zmaxh.Gt.zmax)) Call stop_program('xmin, xmax, zmin or zmax for PH bc in <bcRhizotron.in> not within soil domain. Program terminated.')
       xminh=Floor(xminh/dxgrid)*dxgrid
       xmaxh=Floor((xmaxh+dxgrid)/dxgrid)*dxgrid !incl. right nodes...
       zminh=Floor(zminh/dzgrid)*dzgrid
       zmaxh=Floor((zmaxh+dzgrid)/dzgrid)*dzgrid
       nh=NINT(((xmaxh-xminh)/dxgrid)*((zmaxh-zminh)/dzgrid))
       If (nh.Gt.maxbdr) Call stop_program('Input from  < bc.in >  exceeds parameter "maxbdr". Program terminated.')
       aa=NINT(Abs(xmaxh-xminh)/dxgrid)
       bb=NINT(Abs(zmaxh-zminh)/dzgrid)
        
       Do j=1,bb
         Do k=1,aa 
            iBCPt(ii)=NINT(((zmax-zmaxh)/dzgrid)*nx*ny+(xminh-xmin)/dxgrid+k+(j-1)*nx*ny)
            ii=ii+1
         End Do
       End Do
         Do ii=nh_tot+1,nh_tot+nh
            Kode(iBCPt(ii))=+3
            If (l.Eq.1) Kcell(iBCPt(ii))=1
            If (l.Eq.2) Kcell(iBCPt(ii))=2
            If (l.Eq.3) Kcell(iBCPt(ii))=3
            If (l.Eq.4) Kcell(iBCPt(ii))=4
            If (l.Eq.5) Kcell(iBCPt(ii))=5
            If (l.Eq.6) Kcell(iBCPt(ii))=6
            If (l.Eq.7) Kcell(iBCPt(ii))=7
            If (l.Eq.8) Kcell(iBCPt(ii))=8
            If (l.Eq.9) Kcell(iBCPt(ii))=9
         End Do
         nh_tot=nh_tot+nh
         nBCPts=nBCPts+nh
    End Do
    Read (10,*)
    Read (10,*)
    Read (10,*) nhbcCh
    If (nhbcCh.Gt.mxBcCh) Call stop_program('Input from  < bcRhizotron.in >  exceeds parameter "mxBcCh". Program terminated.')
    Read (10,*)
    If (nhbcCh.Eq.0) Then
       Read (10,*)
    Else
       Read (10,*) (thbcCh(i),hbc1(i),hbc2(i),hbc3(i),hbc4(i),hbc5(i),hbc6(i),hbc7(i),hbc8(i),hbc9(i),i=1,nhbcCh)
    End If
 ELSE
    Open (Unit=10,FILE='in/bc.in',STATUS='OLD',ERR=40)
    Write(15,'(//''++ Soil Boundary Conditions ++'')',advance='no')
    Write(15,'(/''------------------------------'')',advance='no')
! flux BC: (t)-BC [L^3/T]  (Kode(i)=-1)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) qfun,ltop
    Read (10,*)

    Select Case(qfun)
    Case(0) !> + case0: no flux
       Read (10,*)
       Write(15,'(/''* No flux'')',advance='no')
       nQ=0
       noBCflux=.true.
    Case(1) !> + case1: flux BC on total top surface
       Read(10,*)
       Write(15,'(/''* Infiltration top BC'')',advance='no')
       nQ=nx*ny
       If (ltop) Then
          Do i=1,nQ
             iBCPt(i)=i
          Enddo
       Else
          ii=nPt-nQ+1
          Do i=1,nQ
             iBCPt(i)=ii
             ii=ii+1
          Enddo
       End If
       Do  i=1,nQ
          Kode(iBCPt(i))=-1
       End Do
    Case(2) !> + case2: top domain is partially irrigated on a strip between xqmin and xqmax
       Write(15,'(/''* Irrigation on a strip'')',advance='no')
       Read(10,*) xqmin, xqmax
       If ((xqmin.Lt.xmin).Or.(xqmax.Gt.xmax)) Call stop_program('xmin or xmax for Q bc in <bc.in> not within soil domain. Program terminated.')
       xqmin=Floor(xqmin/dxgrid)*dxgrid
       xqmax=Floor((xqmax+dxgrid)/dxgrid)*dxgrid !incl. right nodes...
       nQ=NINT((xqmax-xqmin)/dxgrid*ny)
       aa=NINT(Abs(xqmax-xqmin)/dxgrid)
       i=1
       Do j=1,ny
          Do k=1,aa
             iBCPt(i)=NINT((xqmin-xmin)/dxgrid+k+(j-1)*nx)
             i=i+1
          End Do
       End Do
       Do i=1,nQ
          Kode(iBCPt(i))=-1
       End Do
    Case(3) !> + case3: top domain is partially irrigated on two strips between xqmin1 and xqmax1 and from xqmin2 and xqmax2
       Write(15,'(/''* Irrigation on two strips'')',advance='no')
       homogene = 2 !non-homogeneous distribution of BC -- time and space dependent changes of flux
       Read(10,*) xqmin1, xqmax1, xqmin2, xqmax2
       If ((xqmin1.Lt.xmin).Or.(xqmin2.Lt.xmin).Or.(xqmax1.Gt.xmax).Or.(xqmax2.Gt.xmax)) Call stop_program('boundaries for partial root zone irrigation for Q bc in <bc.in> not within soil domain. Program terminated.')
       If((xqmin2.Lt.xqmin1).Or.(xqmax1.Gt.xqmax2)) Call stop_program('xmin1 & xmax1 should be smaller than xmin2 & xmax2, Program terminated from <bc.in> ')
       If(xqmin2.Le.xqmax1) Call stop_program('no double allocation of the same nodes possible for partial root zone irrigation; Program terminated from <bc.in>')
       xqmin1=Floor(xqmin1/dxgrid)*dxgrid
       xqmax1=Floor(xqmax1/dxgrid)*dxgrid
       xqmin2=Floor((xqmin2+dxgrid)/dxgrid)*dxgrid
       xqmax2=Floor((xqmax2+dxgrid)/dxgrid)*dxgrid
       nQ1=NINT((xqmax1-xqmin1)/dxgrid*ny)
       aa=NINT(Abs(xqmax1-xqmin1)/dxgrid)
       i=1
       Do j=1,ny
          Do k=1,aa
             iBCPt(i)=NINT((xqmin1-xmin)/dxgrid+k+(j-1)*nx)
             i=i+1
          End Do
       End Do
       !=nQ1+1start allocating nodes from the end of the first part
       nQ=nQ1+NINT((xqmax2-xqmin2)/dxgrid*ny)
       aa=NINT(Abs(xqmax2-xqmin2)/dxgrid)
       Do j=1,ny
          Do k=1,aa
             iBCPt(i)=NINT((xqmin2-xmin)/dxgrid+k+(j-1)*nx)
             i=i+1
          End Do
       End Do
       Do i=1,nQ
          Kode(iBCPt(i))=-1
       End Do
    Case(4) !> + case4: irrigation only on material #1 (e.g. cylindrical domain where outer filling material = #2)
       Write(15,'(/''* Irrigation on material 1 only'')',advance='no')
       Read(10,*)
       dummy = nx*ny
       nQ = 0
       Do i=1,dummy
          If(MatNum(i).Eq.1) Then
             nQ = nQ+1
             iBCPt(nQ)=i
          End If
       Enddo
       Do  i=1,nQ
          Kode(iBCPt(i))=-1
       End Do
    Case DEFAULT
       Call stop_program('Sorry, qfun in <bc.in> can only take a value between 0 and 4. Program terminated.')
    End Select
    Read (10,*)
    Read (10,*) nQbcCh
    If (nQbcCh.Gt.mxBcCh) Call stop_program('Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.')

    If((lclimate).And.(.Not.lSign_new))Then
       Write(15,'(/''* Climate BC at soil surface'')',advance='no')
       Do i=1,nclimaticdata
          tQbcCh(i)=time_climate(i)
          Qbcrec(i)=Precip(i) !(Abs(Precip(i))-Abs(E_pot(i)))
          nQbcCh=nclimaticdata
       End Do
       nQ=0
       nQ=nx*ny
       Do i=1,nQ
             iBCPt(i)=i
       Enddo
       Do  i=1,nQ
          Kode(iBCPt(i))=-1
       End Do
       Qbc(:,1)=Qbcrec
       Read (10,*)
       Read (10,*)
    Else
       Read (10,*)
       If(nQbcCh.Eq.0) Then
          READ(10,*)
       Else
          If (homogene.Eq.1) Then
             Do i=1,nQbcCh
                Read (10,*) tQbcCh(i),Qbcrec(i)
             End Do
             Qbc(:,1)=Qbcrec
          Else
             If(qfun.Eq.3) Then
                Do j=1,nQbcCh
                   Read (10,*) (Qbcrec(i),i=1,3) !only 2 compartments possible
                   tQbcCh(j)=Qbcrec(1)    !> \param tQbcCh time of a top flux BC
                   Qbc(j,1:nQ1)=Qbcrec(2) !> \param Qbc top flux [L T-1]
                   Qbc(j,nQ1+1:nQ)=Qbcrec(3)
                Enddo
             Else
                Do j=1,nQbcCh
                   Read (10,*) (Qbcrec(i),i=1,nQ+1) !on peut enlever une dimension à qbctot!
                   tQbcCh(j)=Qbcrec(1)
                   Qbc(j,1:nQ)=Qbcrec(2:nQ+1)
                Enddo
             Endif
          Endif
       Endif
    Endif


    If (nQ.Eq.0) Then
       Write(15,'(/''* No flux B.C.'')',advance='no')
    Else
       Write(15,'(/''* Nodes with flux B.C.:'',i5)',advance='no') nQ
       Write(15,'(/''* Number of time imposed flux B.C.:'',i5)',advance='no')nQbcCh
    Endif
    !> Irrigators [L/T]  (Kode(i)=-3):
    !> flux will be multiplied by surface of 1 voxel
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) nI
    If (nI.Gt.maxIrrig) Call stop_program('Input from  < bc.in >  exceeds parameter "maxIrrig". Program terminated.')
    Read (10,*)
    If(nI.Eq.0)Then
       Read(10,*)
    Else
       Write(15,'(/''* Irrigators at soil surface'')',advance='no')
       homogene=2
       Do i=1,nI
          Read (10,*) xIrrig(i),yIrrig(i),zIrrig(i)
       End Do
    End If
    Read (10,*)
    Read (10,*) nIbcCh
    If (nIbcCh.Gt.mxBcCh) Call stop_program('Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.')
    Read (10,*)
    If(nIBcCh.Eq.0)Then
       Read(10,*)
    Else
       Do i=1,nI
          Read (10,*) (tIbcCh(i,j),Ibc(i,j),j=1,nIbcCh)
       End Do
    End If
    Do i=1,nI
       If ((xIrrig(i).Gt.xmax).Or.(xIrrig(i).Lt.xmin)) Call stop_program('Irrigator out of soil domain (X). Program terminated.')
       If ((yIrrig(i).Gt.ymax).Or.(yIrrig(i).Lt.ymin)) Call stop_program('Irrigator out of soil domain (Y). Program terminated.')
       If ((zIrrig(i).Gt.zmax).Or.(zIrrig(i).Lt.zmin)) Call stop_program('Irrigator out of soil domain (Z). Program terminated.')
       ixmin=1+((xIrrig(i)-xmin)/dxgrid)
       ixymin=ixmin+((yIrrig(i)-ymin)/dygrid)*nx
       imin=ixymin+(Abs(zIrrig(i)-zmax)/dzgrid)*nl
       If (ixmin.Ne.Floor(ixmin)) Call stop_program('Irrigator position not equal to a node position (X). Program terminated.')
       If (ixymin.Ne.Floor(ixymin)) Call stop_program('Irrigator position not equal to a node position (Y). Program terminated.')
       If (imin.Ne.Floor(imin)) Call stop_program('Irrigator position not equal to a node position (Z). Program terminated.')
       iBCPt(nQ+i)=Int(imin)
    Enddo
    Do i=1,nI
       Kode(iBCPt(nQ+i))=-3
    Enddo
!> Pressure head BC; h(t)-BC [L]  (Kode(i)=+1):
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) hfun,ltop
    Read (10,*)
    If ((ltop .And. (qfun.Gt.0)) .And. (hfun.Gt.0)) Call stop_program('Top boundary condition can either be flux or pressure head. Please redefine hfun and/or qfun in <bc.in>. Program terminated.')
    Select Case(hfun)
    Case(0) !> + Case0: no PH BC
       Write(15,'(/''* No PH top BC'')',advance='no')
       Read(10,*)
       nh = 0
    Case (1) !> + Case1: total surface BC
       Write(15,'(/''* PH top BC'')',advance='no')
       Read (10,*)
       If (ltop) Then
          nh=nx*ny
          Do i=1,nh
             iBCPt(nQ+nI+i)=i
          Enddo
          Do i=1,nh
             Kode(iBCPt(nQ+nI+i))=+1
          End Do
       Else
          If(qfun.Eq.4) Then !irrigation on cylinder --> also cylindrical bottom bc
             dummy = nx*ny
             nh = 0
             ii = nPt-dummy
             Do i=1,dummy
                ii=ii+1
                If(MatNum(ii).Eq.1) Then
                   nh = nh+1
                   iBCPt(nQ+nI+nh) = ii
                   Kode(iBCPt(nQ+nI+nh)) = +2
                End If
             Enddo
          Else
             nh=nx*ny
             ii=nPt-nh+1
             Do i=1,nh
                iBCPt(nQ+nI+i)=ii
                ii=ii+1
             End Do
             Do i=1,nh
                Kode(iBCPt(nQ+nI+i))=+2
             End Do
          End If
       Endif
    Case(2) !> + Case2: PH BC on partial surface area between xmin and xmax
       Write(15,'(/''* PH top BC on partial surface area'')',advance='no')
       Read (10,*) xhmin,xhmax
       If ((xhmin.Lt.xmin).Or.(xhmax.Gt.xmax)) Call stop_program('xmin or xmax for PH bc in <bc.in> not within soil domain. Program terminated.')
       xhmin=Floor(xhmin/dxgrid)*dxgrid
       xhmax=Floor(xhmax+dxgrid/dxgrid)*dxgrid !incl. right nodes...
       nh=NINT((xhmax-xhmin)/dxgrid*ny)
       If (nh.Gt.maxbdr) Call stop_program('Input from  < bc.in >  exceeds parameter "maxbdr". Program terminated.')
       aa=NINT(Abs(xhmax-xhmin)/dxgrid)
       If (ltop) Then
          i=1
          Do j=1,ny
             Do k=1,aa
                iBCPt(nQ+nI+i)=NINT((xhmin-xmin)/dxgrid+k+(j-1)*nx)
                i=i+1
             End Do
          End Do
          Do i=1,nh
             Kode(iBCPt(nQ+nI+i))=+1
          End Do
       Else
          i=1
          Do j=1,ny
             Do k=1,aa
                iBCPt(nQ+nI+i)=NINT((nPt-nx*ny)+((xhmin-xmin)/dxgrid+k+(j-1)*nx))
                i=i+1
             End Do
          End Do
          Do i=1,nh
             Kode(iBCPt(nQ+nI+i))=+2
          End Do
       End If
    Case DEFAULT
       Call stop_program('Sorry, hfun in <bc.in> can only take a value of 0, 1 or 2. Program terminated.')
    End Select

    Read (10,*)
    Read (10,*) nhbcCh
    If (nhbcCh.Gt.mxBcCh) Call stop_program('Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.')
    Read (10,*)
    If (nhbcCh.Eq.0) Then
       Read (10,*)
    Else
       Read (10,*) (thbcCh(i),hbc(i),i=1,nhbcCh)
    End If
    If (nh.Eq.0) Then
       Write(15,'(/''* No head B.C.'')',advance='no')
    Else
       Write(15,'(/''* Nodes with head B.C.:'',i5)',advance='no')nh
       Write(15,'(/''* Number of time imposed head B.C.:'',i5)',advance='no')nhbcCh
    Endif
    !> free drainage-BC  (Kode(i)=-2):
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) lFrDr

    If (.Not.ltop .And. lFrdr) Call stop_program('Bottom BC is already defined as pressure head -- no free drainage possible')
    If(lFrDr) Then
       nFrDr = nx*ny
       ii = nx*ny*nz - nFrDr +1
       Do i=1,nFrDr
          iBCPt(nQ+nI+nh+i)=ii
          ii=ii+1
       End Do
       Do i=1,nFrDr
          Kode(iBCPt(nQ+nI+nh+i))=-2
       End Do
    Else
       nfrdr = 0
    End If

    If (.Not. lFrDr) Then
       Write(15,'(/''* No free drainage nodes.'')',advance='no')
    Else
       Write(15,'(/''* Bottom B.C. is free drainage'')',advance='no')
    Endif
!!$  DO i=1,nFrdr
!!$     Kode(iBCPt(nQ+nI+nh+i))=-2
!!$  END DO
    nBCPts=nQ+nI+nh+nFrdr
    Read (10,*)
    Read (10,*)
    Read (10,*)
    !> solute transport BC:
    !> the type of BC to be invoked for solute transport (KodCB>0 first type; KodCB<0 third type)
    Read (10,*) assCode
    KodCB=0
    Do i=1,nBCPts/2
       KodCB(i)=assCode
    Enddo
    Read (10,*)
    Read (10,*)
    Read (10,*) nCBnd1, nCBnd2
    If (nCBnd1.Gt.mxBcCh) Call stop_program('Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.')
    If (nCBnd2.Gt.mxBcCh) Call stop_program('Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.')
    Read (10,*)
    If (nCBnd1.Gt.0) Then
       Read (10,*) (tCBnd1(i),CBnd1(i),i=1,nCBnd1)
    Else
       Read (10,*)
    End If
    Read (10,*)
    If (nCBnd2.Gt.0) Then
       Read (10,*) (tCBnd2(i),CBnd2(i),i=1,nCBnd2)
    Else
       Read (10,*)
    End If
    Read (10,*)
    Read (10,*)
    Read (10,*) tPulse
    If ((nCBnd1==0).And.(nCBnd2==0)) Then
       Write(15,'(/''* No solute transport BCs'')',advance='no')
    Elseif (nCBnd1.Ne.0) Then
       Write(15,'(/''* 1st time dependent solute transport B.C.: '',i5)',advance='no')nCBnd1
    Elseif (nCBnd2.Ne.0) Then
       Write(15,'(/''* 2nd time dependent solute transport B.C.: '',i5)',advance='no')nCBnd2
    Endif
    Close (10)
 END IF

    ! reads decayrate and uptake order from partrace input
    If(lPartUp.Or.lSalinity) Then
       Open (Unit=10,FILE='Input_PARTRACE.pTraceInpV12',STATUS='OLD',ERR=70)
       Do i=1,82
          Read (10,*)
       End Do
       Read (10,*) decayorder
       Read (10,*) decayrate
       Read (10,*)
       Read (10,*) uptakeorder
       Close (10)
    End If
    If (decayorder.EQ.3) l_degrad = .TRUE.
    
!> #### if input file is available, reads probes.in ####
    Open (Unit=10,FILE='in/Probes.in',STATUS='OLD',ERR=111)
    ObsOK=.True.
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) npr,dtProbe
    Write(15,'(/''* Number of sampling probes:'',i5)')npr
    Call IniProbes
    Read (10,*)
    Read (10,*)
    Do i=1,npr
       Read (10,*) Pr(i),Pt(i),CrP(i)
    Enddo
    Read (10,*)
    Read (10,*)
    Read (10,*) (VarP(i),i=1,npr)
    Read (10,*)
    Read (10,*)
    Read (10,*) (distrP(i),i=1,npr)
    firstOK=.True.
    Do i=1,npr
       If (Pt(i)==4) Then !user define node
          If (firstOK) Then
             Read (10,*)
             Read (10,*)
             firstOK=.False.
          Endif
          NodebyPr(i)=Crp(i)
          Read (10,*) (node_temp(i2),i2=1,CrP(i)+1)
          Do i2=1,CrP(i)
             NodePr(i,i2)=node_temp(i2+1)
          Enddo
       Endif
    Enddo
    !time
    If (dtprobe.Ne.999) Then
       nouProbe=0
       ini=0
       Do While (ini.Lt.tmax)
          nouProbe=nouProbe+1
          tOuProbe(nouProbe)=ini+dtprobe
          ini=tOuProbe(nouProbe)
          If (nouProbe.Gt.mxProf) Then
             Print *,'too small time step in probe.in: only ',mxProf,' profiles will be kept'
             Goto 78
          Endif
       Enddo
    Endif
78  Close (10)

    !run first calculations
    Call NodebyProbe
    Deallocate(node_temp)
    Return
10  Call stop_program('File  < control.in >  not found -- program terminated.')
20  Call stop_program('File  < grid.in >  not found -- program terminated.')
30  Call stop_program('Input file for nodal values not found -- program terminated.')
40  Call stop_program('File  < bc.in >  not found -- program terminated.')
60  Call stop_program('File  < mesh.in >  not found -- program terminated.')
70  Call stop_program(' File  <Input_PARTRACE.pTraceInpV12>  not found --')
111 Write (15,'(/''* File  <Probes.in>  not found: no defined observation locations'')')
  End Subroutine Applic
  !******************************************************************************
  !> ### reads Doussan RWU related inputs ###
  Subroutine DouIn(ipl)
    Use Typedef
    Use TardieuData, Only : gsmin,gsalpha,gsbeta,gsdelta
    Use DoussanMat, Only : ageLr,LrRoot,nLr,ageKh,Khroot,nKh,hx_min,ave,eqDis,stresval1,&
         stresval2,lcavit,cavitb,cavitc,stresfun
    Use PlntData, Only : a_r,a_1,a_2
    Use RootData, Only : lGap,lAQPc,nAQPc,AQPh,AQPv,lSign,lSign_new,lSign_inst,&
        delta_h,nplant,sign_in,PH_crit,size_buff,vol_buff,sign_in,lKdrop,TypeKdrop,&
        PH_crit,g1,g2,TR1,seglen,segsur,vol_root,fac_contact
    Use ParamData, Only : pi,lPartUp,lSalinity,n_MFP
    Use SoluteRootMat, Only : nPerm,agePr,PrRoot,nPass,agePs,PsRoot,sorp,l_linSorb,l_freundSorb,&
         theta_R,rho_R,numPa
    Use MPIutils, Only: stop_program
    !Use WatFun, Only : hTab_MFP,MFPTab!,hcheck,mfpcheck,hnewcheck
    Implicit None

    Integer(ap) :: i,i_ave,i_eqDis
    INTEGER(ap),INTENT(in) :: ipl
    CHARACTER :: infile*24

! Open simulation summary file
    Open(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
    Write(15,'(//''++ Hydraulic Data for roots ++'')',advance='no')
    Write(15,'(/''-----------------------------'')',advance='no')

 !> #### reads CondRoot.in ####
    WRITE(infile, "(A14,I1)")'in/CondRoot.in',ipl
    If (nplant.GT.1) Then  
        Open (Unit=10,FILE=infile,STATUS='OLD',ERR=3001)
    Else 
        Open (Unit=10,FILE='in/CondRoot.in',STATUS='OLD',ERR=3001)
    End If
    Write(15,'(/''* The Doussan algorithm will be used for the root water uptake (lDou=TRUE)'')',advance='no')
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) (nLr(i,ipl),i=1,3)
    Read (10,*)
    !> \param LrRoot radial root hydraulic conductivity [L/T/L]; per root order (must be defined for 3 orders)
  READ (10,*) (ageLr(1,i,ipl),LrRoot(1,i,ipl),i=1,nLr(1,ipl))  !main axes
  READ (10,*) (ageLr(2,i,ipl),LrRoot(2,i,ipl),i=1,nLr(2,ipl))  !secondary axes
  READ (10,*) (ageLr(3,i,ipl),LrRoot(3,i,ipl),i=1,nLr(3,ipl))  !tertiary axes
  READ (10,*)
  READ (10,*)
  READ (10,*) (nKh(i,ipl),i=1, 3)
  READ (10,*)
  !> \param KhRoot axial root hydraulic conductance [L4/T/L]; per root order (must be defined for 3 orders)
  READ (10,*) (ageKh(1,i,ipl),KhRoot(1,i,ipl),i=1,nKh(1,ipl))
  READ (10,*) (ageKh(2,i,ipl),KhRoot(2,i,ipl),i=1,nKh(2,ipl))
  READ (10,*) (ageKh(3,i,ipl),KhRoot(3,i,ipl),i=1,nKh(3,ipl))
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) (nPerm(i,ipl),i=1, 3)
  READ (10,*)
!> \param PrRoot radial root solute permeability [L/T]; per root order (must be defined for 3 orders)
  IF(lPartUp.OR.lSalinity)THEN
     READ (10,*) (agePr(1,i,ipl),PrRoot(1,i,ipl),i=1,nPerm(1,ipl))
     READ (10,*) (agePr(2,i,ipl),PrRoot(2,i,ipl),i=1,nPerm(2,ipl))
     READ (10,*) (agePr(3,i,ipl),PrRoot(3,i,ipl),i=1,nPerm(3,ipl))
     READ (10,*)
     READ (10,*)
     READ (10,*) (nPass(i,ipl),i=1, 3)
     READ (10,*)
!> \param PrRoot radial root solute permeability [L/T]; per root order (must be defined for 3 orders)
     READ (10,*) (agePs(1,i,ipl),PsRoot(1,i,ipl),i=1,nPass(1,ipl))
     READ (10,*) (agePs(2,i,ipl),PsRoot(2,i,ipl),i=1,nPass(2,ipl))
     READ (10,*) (agePs(3,i,ipl),PsRoot(3,i,ipl),i=1,nPass(3,ipl))
       Read (10,*)
       Read (10,*)
       Read (10,*) l_linSorb, l_freundSorb
       Read (10,*)
       If(l_linSorb) Then
          Read (10,*) sorp(1)
       Elseif (l_freundSorb) Then
          Read (10,*) sorp(1), sorp(2)
       Else
          Read(10,*)
       End If
          Read (10,*)
          Read (10,*)
          Read (10,*)theta_R,rho_R
          Read (10,*)
          Read (10,*)
          Read (10,*)numPa
    Else
       Do i=1,21
          Read(10,*)
       End Do
    End If
    !averaging method
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) i_ave, i_eqDis
    If (i_ave .Eq. 0) Then
       ave=.False.
    Elseif (i_ave .Eq. 1) Then
       ave=.True.
       Write(15,'(/''* Averaging method: treated as one root (ave= TRUE)'')',advance='no')
    Else
       Call stop_program('Set AVERAGE properly in Condroot.in')
    Endif
    If (i_eqDis .Eq. 0) Then
       eqDis=.False.
       If (i_ave .Eq. 0) Then
          Write(15,'(/''* No averaging method will be used below voxel scale (ave= FALSE and eqDis= FALSE)'')',advance='no')
       Endif
    Elseif (i_eqDis .Eq. 1) Then
       eqDis=.True.
       Write(15,'(/''* Averaging method: use equal distance function (eqDis= TRUE)'')',advance='no')
    Else
       Call stop_program('Set eqDis parameter properly in Condroot.in')
    Endif
    !> 1. stress functions
    Read (10,*,ERR=20)
    Read (10,*,ERR=20)
    Read (10,*,ERR=20) stresfun
    Read (10,*,ERR=20)
    Select Case(stresfun)
!> + case0: No plant water stress responses considered
    Case(0) 
       Read (10,*,ERR=20)
       Write(15,'(/''* No plant water stress is considered '')')
!> + case1: Hydraulic regulation of stomata at a limiting stress value (isohydric)
    Case(1) 
       Read (10,*,ERR=20) hx_min !> \param hx_min limitin PH_collar
       hx_min=-Abs(hx_min)
       Write(15,'(/''* Hydraulic regulation of stomata with a limiting stress value = '',1pE11.4,1pE11.4)',advance='no') &
            hx_min
!> + case2: Hydraulic regulation of stomata; linear stress function
    Case(2) 
       Read (10,*,ERR=20) stresval1,stresval2
       Write(15,'(/''* A linear stress function will be used with values = '',1pE11.4,1pE11.4)',advance='no') &
            stresval1,stresval2
!> + case3: Hydraulic regulation of stomata; Tuzet et al. (2003) function
    Case(3) 
       Read (10,*,ERR=20) stresval1,stresval2
       Write(15,'(/''* The Tuzet stress function will be used with values  = '',1pE11.4,1pE11.4)',advance='no') &
            stresval1,stresval2
 !> + case4: Hydraulic and hormonal regulation of stomata; Tardieu et al. (1993) function; add. hormone transport
    Case(4)
       lSign=.True.
       Read (10,*,ERR=20) a_r,a_1,a_2,PH_crit,delta_h,sign_in,size_buff
       Write(15,'(/''* Tardieu&Davies stomatal conductance will be calculated with signal transport; '',2(1pE11.4,1X))')  &
            a_r,a_1,a_2
       If(size_buff.Gt.0) Then
          vol_buff = vol_root(ipl)*size_buff
       Else
          vol_buff = segsur(1,ipl)*segsur(1,ipl)/seglen(1,ipl)/4/pi
       End If
       TR1 = Abs(PH_crit)-Abs(delta_h)
!> + case5: Hydraulic and hormonal regulation of stomata; Tardieu et al. (1993) function; NO hormone transport
    Case(5) 
       lSign_inst=.True.
       Read (10,*,ERR=20) a_r,a_1,a_2,PH_crit,delta_h,sign_in,size_buff
       Write(15,'(/''* Tardieu&Davies stomatal conductance will be calculated with instantaneous signaling; '',2(1pE11.4,1X))')  &
            a_r,a_1,a_2
       If(size_buff.Gt.0) Then
          vol_buff = vol_root(ipl)*size_buff
       Else
          vol_buff = segsur(1,ipl)*segsur(1,ipl)/seglen(1,ipl)/4/pi
       End If
       TR1 = Abs(PH_crit)-Abs(delta_h)
!> + case6:! Tardieu and Davies Updated
    Case(6) 
       lSign_new=.True.
       Read (10,*,ERR=20) gsmin,gsalpha,gsbeta,gsdelta
       Write(15,'(/''* New Tardieu and Davies stomatal conductance will be ; '',4(1pE11.4,1X))')  &
            gsmin,gsalpha,gsbeta,gsdelta
    Case DEFAULT
       Call stop_program('Sorry, stresfun in <bc.in> can only take values from 0-6. Program terminated.')
    End Select
	
!> \Additional impacts on radial root conductivity
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=60) lKdrop !> 1. Soil conductivity drop around roots
    If (lKdrop) Then
       Write (*,*) 
       Write(15,'(/''* Soil conductivity drop around roots will be considered (lKdrop=TRUE)'')',advance='no')
       Read (10,*)
       Read (10,*,ERR=25) TypeKdrop !defines how K drop is modeled
       Read (10,*)
       Read (10,*,ERR=26) fac_contact ! fraction of root surface in contact with the rhizosphere
       Call MFPC_table
! !this is just a check for matric flux potential claulations
       ! Open(UNIT=16,FILE='out/MFPTable.out',STATUS='UNKNOWN')
       ! Write(16,'(//''++ Matric Flux potential Table ++'')',advance='no')
       ! Write(16,'(/''-------------------------'')',advance='no')
       ! Write(16,'(/''   PH    PHI '')')
       ! Do i=1,n_MFP
          ! Write (16,'(5(1X,1pE10.3))') hTab_MFP(i),MFPTab(i,1),hnewcheck(i),mfpcheck(i,1),hcheck(i)
       ! End Do 
       ! Close (16)
       ! Write(*,'('' lKdrop table calculated'')')
    Else
        Write(15,'(/''* Soil conductivity drop around roots will not be considered (lKdrop=false)'')',advance='no')
        Read (10,*)
        Read (10,*)
        Read (10,*)
        Read (10,*)
    Endif
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=30) lcavit !> 2. cavitational effect on root axial conductance
    If (lcavit) Then
       Write (*,*) 
       Write(15,'(/''* Cavitation effect on root axial conductance will be considered (lcavit=TRUE)'')',advance='no')
       Read (10,*)
       Read (10,*,ERR=30) cavitb,cavitc
    Else
       Write(15,'(/''* No cavitation effect on root axial conductance will be considered (lcavit=false)'')',advance='no')
       Read (10,*)
       Read (10,*)
    Endif
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=40) lGap !> 3. Air gap / Rhizosphere gap effect on radial root conductivity
    If (lGap) Then
       Write(15,'(/''* Air gap / rhizosphere hydrophobicity effect on radial conductivity will be considered (lgap=TRUE)'')',advance='no')
       Read (10,*)
       Read (10,*,ERR=40) g1,g2
    Else
       Write(15,'(/''* No Air gap / rhizosphere hydrophobicity effect on radial conductivity will be considered (lgap=false)'')',advance='no')
       Read (10,*)
       Read (10,*)
    Endif
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=50) lAQPc, nAQPc !> 4. Aquaporin effect on radial root conductivity
    If (lAQPc) Then
       Write(15,'(/''* Aquaporin status effect on radial conductivity will be considered (lAQPc=TRUE)'')',advance='no')
       Read (10,*)
       Read (10,*,ERR=50) (AQPh(i,ipl),AQPv(i,ipl),i=1,nAQPc)
    Else
       Write(15,'(/''* No aquaporin status effect on radial conductivity will be considered (lAQPc=false)'')',advance='no')
       Read (10,*)
       Do i = 1,nAQPc
          Read (10,*)
       Enddo
    Endif
    If ((ave) .And. (eqDis)) Call stop_program('set one root and equidistant approach correctly; may not be both switched on')
    Close (15)
    Return
    Close (10)
3001 Call stop_program(' File  <CondRoot.in>  not found --')
20  Call stop_program(' Water stress information not found in <CondRoot.in>. -- program terminated.')
25  Call stop_program(' Kdrop not found in <CondRoot.in>. -- program terminated.')
26  Call stop_program(' fraction not found in <CondRoot.in>. -- program terminated.')
30  Call stop_program(' Cavitation information not found in <CondRoot.in>. -- program terminated.')
40  Call stop_program(' Gap / hydrophobicity information not found in <CondRoot.in>. -- program terminated.')
50  Call stop_program(' AQP behavior information not found in <CondRoot.in>. -- program terminated.')
60  Call stop_program(' Kdrop information not found in <CondRoot.in>. -- program terminated.')
  End Subroutine DouIn
  !*******************************************************************************
  ! > Read climate.in
  Subroutine ClimateIn
    Use Typedef
    Use Envidata
    Use RootData, Only : lSign_new
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: i
    Real(dp) :: dummy

    Open(Unit=10,File='in/Climate.in',STATUS='OLD',ERR=10)
    Read (10,*)
    Read (10,*) nclimaticdata
    If (lSign_new) Then !Stressfun in condroot.in=6
       Read (10,*)
       Read (10,*)
       Allocate(time_climate(nclimaticdata))
       Allocate(T_pot(nclimaticdata))
       Allocate(Precip(nclimaticdata))
       Allocate(PPFD(nclimaticdata))
       Allocate(VPD(nclimaticdata))
       Allocate(Temperature(nclimaticdata))
       Read (10,*,ERR=20) (dummy,time_climate(i),T_pot(i),Precip(i),PPFD(i),Temperature(i),VPD(i),i=1,nclimaticdata)
    Else
       Allocate(time_climate(nclimaticdata))
       Allocate(T_pot(nclimaticdata))
       Allocate(Precip(nclimaticdata))
       Read (10,*)
       Read (10,*)
       Read(10,*,ERR=20) (dummy, time_climate(i),T_pot(i),Precip(i),i=1,nclimaticdata)
       !time_climate=t_initial+time_climate/24.
       ! print*,'CLIMATE DATA', dummy, time_climate, Precip, E_pot, T_pot
    End If
    Close (10)
    Print*,'... Reading Climate.in ...'
    Return
10  Call stop_program('Climate.in not found, program terminated')
20  Call stop_program('Problem in reading Climate.in')
  End Subroutine ClimateIn
  !*******************************************************************************
  ! > Treat Climate data
  Subroutine ProcessClimate
    Use Typedef
    Use Envidata
    Use TardieuData, Only : ampliinit
    Implicit None

    Integer(ap) :: count_day=1,i,ndays,j,j_init
    Real(dp) :: sum_daily_PPFD

    Allocate(s_f(nclimaticdata))
    Allocate(rho_f(nclimaticdata))
    Allocate(gamma_f(nclimaticdata))
    Allocate(lambda_f(nclimaticdata))
    Allocate(conversion_f(nclimaticdata))
    Allocate(PPFDcum(nclimaticdata))
    Allocate(Rn(nclimaticdata))

    Allocate(hours(nclimaticdata))
    Allocate(days(nclimaticdata))

    Print*,'... Processing Climate ...'
    ! calc coefficient temperature-dependent
    s_f = 0.1849 * (Temperature** 2) + 1.005 * Temperature + 50.87  ;
    rho_f = 0.00002 * (Temperature**2) - 0.0047 * Temperature +  1.2921  ;
    gamma_f = 0.0643 * Temperature + 64.9  ;
    lambda_f = (-0.0024 * Temperature + 2.5008) * 1000000 ;
    conversion_f  = -0.1424 * Temperature + 43.917;

    Rn = PPFD/2.
    ! calc coefficien temperature-independent
    Cp  =  1012;

    ! calc hours/days
    Do i=1,nclimaticdata
       If (Floor((time_climate(i)-(t_initial))).Lt.count_day) Then
          days(i)=count_day
          hours(i)=(time_climate(i)-(t_initial+count_day-1.))*24.
       Else
          count_day=count_day+1
          days(i)=count_day
          hours(i)=(time_climate(i)-(t_initial)-(days(i)-1.))*24.
       Endif
       !print*,i,count_day,days(i),hours(i),time_climate(i)
    Enddo

    ndays=Maxval(Nint(days))

    Allocate(datad(ndays,4))
    Allocate(ampli(ndays))
    Allocate(psixylmax(ndays))
    Allocate(psixylmin(ndays))

    ampli=ampliinit
    psixylmax=Abs(-500.-ampliinit)
    psixylmin=Abs(-500.)

    Do i=1,ndays
       j=1
       Do While (days(j).Lt.i)
          j=j+1
       Enddo
       j_init=j
       Do While ((PPFD(j).Lt.10).And.(hours(j).Lt.8))
          j=j+1
       Enddo
       If (j.Gt.j_init) Then
          datad(i,1)=Maxval(time_climate(j_init:j-1))-t_initial
          datad(i,2)=datad(i,1)
       Else
          datad(i,2)=datad(i-1,2)+1
          datad(i,1)=datad(i-1,1)
       Endif
       !print*,datad(i,1),datad(i,2)
    Enddo

    Do i=1,ndays
       j=1
       Do While (days(j).Lt.i)
          j=j+1
       Enddo
       j_init=j
       sum_daily_PPFD=0.
       Do While ((j.Le.Size(days)).And.(days(Min(j,Size(days))).Le.i))
          If (j.Gt.j_init) Then
             sum_daily_PPFD=sum_daily_PPFD+PPFD(j)*(time_climate(j)-time_climate(j-1))*24.*60*60./1000000.
             PPFDcum(j)=PPFDcum(j-1)+PPFD(j)*(time_climate(j)-time_climate(j-1))*24.*60*60./1000000.
          Else
             PPFDcum(j)=0.
          Endif
          !print*,PPFDcum(j)
          j=j+1
       Enddo
       If (j.Gt.j_init) Then
          datad(i,3)=Maxval(PPFD(j_init:j-1))
          datad(i,4)=sum_daily_PPFD
       Else
          datad(i,3)=0
          datad(i,4)=0
       Endif
       !print*,datad(i,3),datad(i,4)
    Enddo
  End Subroutine ProcessClimate

  !*******************************************************************************
  ! > Read Tardieu.in
  Subroutine TardieuIn
    Use Typedef
    Use TardieuData
    Use MPIutils, Only: stop_program
    Implicit None
    Open(Unit=10,File='in/Tardieu.in',STATUS='OLD',ERR=10)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) ampliinit,TaucircadGr,TaucircadGxl,TaucircadGstem,TaucircadGc
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) TautranspiGr,TautranspiGxl,TautranspiGstem,TautranspiGc
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) TauABAGr,TauABAGxl,TauABAGstem,TauABAGc
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) Vres,Vsat,ncap,alphacap
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) Grmin,Grmax,Gxlmin,Gxlmax,Gstemmin,Gstemmax,Gcmin,Gcmax,Gxl0,Gstem0,Gc0
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) ABAconstit,ABAa,ABAb
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) S,shading,ga

    Sshaded=S*shading
    Close (10)
    Print*,'... Reading Tardieu.in ...'
    Return
10  Call stop_program('Tardieu.in not found, program terminated')
20  Call stop_program('Problem in reading Tardieu.in')
  End Subroutine TardieuIn

  !*******************************************************************************
  !> ### Reads solute related inputs ###
  Subroutine ChemIn
    Use Typedef
    Use SolData
    Use PlntData, Only: cMm,Vmax,xin,fk
    Use RootData, Only: nUrf,age,Urf
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: i,k

    Open(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
    Write(15,'(//''++ Solute Transport Information ++'')',advance='no')
    Write(15,'(/''----------------------------------'')')
    Open(Unit=10,File='in/chem.in',STATUS='OLD',ERR=10)
    Read (10,*)
    Read (10,*)
    NLevel=1
    Read (10,*) epsi,PeCr
    If (epsi.Lt.0.999_dp) NLevel=2
    Read (10,*)
    Read (10,*)
    Do k=1,nMat
       Read (10,*,ERR=20) (ChPar(i,k),i=1,9)
    End Do
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) CMm,VMax,xin,fk
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=20) (nUrf(i),i=1,3)
    Read (10,*,ERR=20)
    Do  i=1,3
       Read (10,*,ERR=20) (age(i,k),Urf(i,k),k=1,nUrf(i))
    End Do
    Close (10)
    Write (15,'('' -- simulation will include solute transport.'')',advance='no')
    Close (15)
    Return
10  Call stop_program('File  < chem.in >  not found  -- program terminated')
20  Call stop_program('Insufficient data in  < chem.in >  -- program terminated.')
  End Subroutine ChemIn
 !****************************************************************************
 !> ### reads soil related inputs ###
  Subroutine SoilIn
    Use Typedef
    Use SolData, Only: matNum,par,nMat,lMacro,soiltab,ssMaxTab
    Use WatFun
    Use RhizoData, Only: rhizoModel, lRhizo, bulkPara
    Use MPIutils, Only: stop_program
    Use GridData, Only: nPt,Axy,denomSe
    Implicit None

    Integer(ap) ::i,k,m,ii
    Real(dp) ::alh,h1,hn,htab1,htabn,ths,thr

    Open(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
    Write(15,'(//''++ Soil Information ++'')',advance='no')
    Write(15,'(/''----------------------'')',advance='no')
!> Check which  Rhizosphere model is considered.
    If (lRhizo) Then
       Open(Unit=11,FILE='in/Rhizo.in',STATUS='OLD',ERR=60)
       Write(15,'(//''* Rhizosphere mucilage properties are considered'')',advance='no')
       Read(11,*)
       Read(11,*)
       Read(11,*)
       Read(11,*,ERR=60) rhizoModel
       Close(11)
    Endif
!>  reads soil.in 
    Open(Unit=10,FILE='in/soil.in',STATUS='OLD',ERR=10)
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*,ERR=30)nMat,h1,hN,nTab
    If (h1*hN.Eq..0_dp) Then
       Write (*,'('' Input error -- please use non-zero values for hTab1 and hTab2 in  < soil.in >.'')')
       Call stop_program('Program terminated.')
    Endif
    Write(15,'(/''* Tabulated values between: '',1pE11.4,'' and '',1pE11.4)',advance='no') h1,hN
    Write(15,'(/''* Number of soil layers:'',1i7)',advance='no') nMat
    hTab1=-Min(Abs(h1),Abs(hN))
    hTabN=-Max(Abs(h1),Abs(hN))
    Read (10,*)
    Read (10,*)
! read soil hydraulic parameters
    Do k=1,nMat
       Read (10,*,ERR=20)(par(i,k),i=1,11)
    EndDo
!> Call for rhizosphere properties
    If (lRhizo)  Call RhizoIn
    Close (10)
	
!create a vector with nodal denominator value for relative water content
    do ii=1,nPt
        m=Matnum(ii)
        If (.NOT. lRhizo) Then
            thr=par(2,m)
            ths=par(3,m)
        else 
            thr=bulkPara(1)
            ths=bulkPara(2)
        endif
        denomSe(ii)=(ths-thr)/Axy(ii)
    Enddo

!> define if the soil types are contrasting (e.g. macropore)
    If (nMat.Gt.1) Then
       Do k=1,nMat
          If(par(11,k).Lt.1E-5_dp) lMacro=.True.
       End Do
    End If
	

!> generate table for interpolation:	
    Call IniTab
    If (Abs(hTabN-hTab1).Gt.0._dp) Then
       alh1=Log10(-hTab1)
       dlh=(Log10(-hTabN)-alh1)/Real(nTab-1)
       Do i=1,nTab
          alh=alh1+Real(i-1)*dlh
          hTab(i)=-(10._dp)**alh
       End Do
       Do k=1,nMat
          Write(15,'(/''   PH    KH    C    TH '')')
          Do i=1,nTab
             ConTab(i,k)=FKP(hTab(i),par(:,k),i)
             CapTab(i,k)=FCP(hTab(i),par(:,k))
             TheTab(i,k)=FTh(hTab(i),par(:,k))
             Write (15,'(4(1X,1pE10.3))') hTab(i),TheTab(i,k),CapTab(i,k),ConTab(i,k)
          End Do
       End Do
    Else
       Write (*,'('' Input error --  hTab1  and  hTab2  in  < >  cannot be equal.'')')
       Call stop_program('Program terminated.')
    Endif
    Return
!> #### if available, reads soiltab.in ####
10  Open(Unit=10,FILE='in/soiltab.in',STATUS='OLD',ERR=40)
    soiltab=.True.
    Read (10,*,ERR=50)
    Read (10,*,ERR=50)
    Read (10,*,ERR=50)
    Read (10,*,ERR=50) nMat, nTab
    Read (10,*,ERR=50)
    Read (10,*,ERR=50)
    Call IniTab

    !> tabulated soil hydraulic properties
    !> col1=PH, col2=WC, col3=CH, col4=KH
    Do i=1,nTab
       Read (10,*,ERR=50) hTab(i),(TheTab(i,k),CapTab(i,k),ConTab(i,k),k=1,nMat)
    Enddo
    Read (10,*,ERR=50)
    Read (10,*,ERR=50) (ssMaxTab(k),k=1,nMat)
    alh1=Log10(Abs(hTab(1)))
    dlh=(Log10(Abs(hTab(nTab)))-alh1)/Real(nTab-1)
    Write(15,'(/''* Tabulated soil input file'')',advance='no')
    Write(15,'(/''* Number of soil layers:'',1i7)',advance='no') nMat
    Write(15,'(/''* Number of tabulations:'',1i7)',advance='no') nTab
    Close (10)
    Return
20  Call stop_program('Insufficient data in  < soil.in >  -- program terminated.')
30  Call stop_program('Soil input file has been changed (25Feb2011): now, add at the 4th line a 4th number giving the number of tabulations (by default, put 100) !')
40  Call stop_program('File  < soil.in > or < soiltab.in > not found -- program terminated.')
50  Call stop_program('Error in < soiltab.in > -- program terminated.')
60  Call stop_program('ERROR in <rhizo.in -- program terminated')
  End Subroutine SoilIn
  !***************************************************************************
  Subroutine RhizoIn
    !> A subroutine that read the rhizosphere parameters and calculate
    !RhizoStaticPara
    Use Typedef
    Use RhizoData, Only: rhizoModel, bulkPara, RhizoPara, RhizoSpatialPara, StaticRhizoPara
    Use RhizoStat, Only: RetCurve, ERR
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: i,j, itNum=0 , fitLen=15000, maxIter=10000, inx_hcr
    Real(dp) :: thtS, thtR, lambd, hcr, omega, beta, cTot, rhob, rhow
    Real(dp) ::  tht1, tht2, f1, f2, rt, dx, thtMid, sumERR, preERR
    Real(dp) :: lambd_stat, epsi
    Real(dp), Allocatable, Dimension(:) :: tempTht, tempH, err1
    Logical :: converge = .False.

    Allocate(tempTht(fitLen), tempH(fitLen))
    tempH = -1.
    epsi  = 0.0001
    Open(Unit=11,FILE='in/Rhizo.in',STATUS='OLD',ERR=20)
    Read(11,*)
    Read(11,*)
    Read(11,*)
    Read(11,*,ERR=20) rhizoModel
    Read(11,*)
    Read(11,*)
    Read(11,*,ERR=20) (bulkPara(i), i=1,6)
    Read(11,*)
    Read(11,*)
    Read(11,*,ERR=20) (RhizoPara(i), i=1,9)
    Read(11,*)
    Read(11,*)
    Read(11,*,ERR=20) (RhizoSpatialPara(i), i=1,2)
    Close(11)
    thtR = bulkPara(1); thtS = bulkPara(2); lambd = bulkPara(3); hcr = bulkPara(4)
    omega = RhizoPara(1); beta = RhizoPara(2); cTot = RhizoPara(3); rhob = RhizoPara(8); rhow = RhizoPara(9)
    !> A procedure to calculate Rhizosphere static parameter. Fitst we use the
    !Secant method to generate tempTht(tempH) relation and after we will fit that
    !data into the brooks and cory model to find lambd_stat, hcr_stat. The static
    !is also used for the initial condition.
    tht1 = thtR
    tht2 = thtS
    Do i=1,fitLen
       tempH(i) =  - i
       f1 = ERR(tht1, tempH(i))
       f2 = ERR(tht2, tempH(i))
       If (f1 * f2 .Ge. 0.) Then
          tempTht(i) = thtS
       Else
          If (f1 < 0.) Then
             rt = tht1
             dx = tht2 - tht1
          Else
             rt = tht2
             dx = tht1 - tht2
          Endif
          Do While(.Not. converge )
             dx = dx*0.5
             thtMid = rt + dx
             f2 = ERR(thtMid, tempH(i))
             If (f2 .Le. 0.) rt = thtMid
             If (Abs(dx) .Lt. epsi) Then
                tempTht(i) = thtMid
                converge = .True.
             Endif
             itNum = itNum + 1
             If (itNum .Ge. maxIter) Then
                Call stop_program('Max Iteration in the secant method reached. see RhizoIn subroutine')
             Endif
          Enddo
       Endif
       !    write(*,*) tempH(i), tempTht(i), RetCurve(tempH(i),hcr,lambd), itNum, abs(dx)
       converge = .False.
       itNum = 0
       tht1 = thtR
       tht2 = thtS

    Enddo
    !< Finding hcr of the static Rhizo
    Do i=1,fitLen
       If (tempTht(i) .Ne. tempTht(i+1)) Then
          inx_hcr = i+1
          StaticRhizoPara(2) = (tempH(i+1) + tempH(i))/2.
          Goto 13
       Endif
       If (i .Eq. fitLen)  Call stop_program('Did not found hcr Static - program terminated ')

    Enddo
13  Continue
    !< A simple iterative procedure to find lambda rhizo. Note that lambda rhizo is
    !always lower then lambda bulk and higher than zero.
    f1=0.
    lambd_stat = lambd
    converge = .False.
    itNum = 0
    Allocate(err1(fitLen-inx_hcr))
    Do While (.Not. converge)
       Do j=1,(fitLen-inx_hcr)
          f1 = RetCurve(tempH(j+inx_hcr),hcr,lambd_stat)
          err1(j) = (f1-tempTht(j+inx_hcr))**2.0
       Enddo
       preERR = sumERR
       sumERR = Sum(err1)
       If (itNum .Eq. 0)  preERR = sumERR + 1.
       If (sumERR .Lt. epsi .Or. preERR .Lt. sumERR) Then
          StaticRhizoPara(1) = lambd_stat
          converge = .True.
       Endif
       lambd_stat = lambd_stat - 0.001
       itNum = itNum + 1
    Enddo
    If (sumERR .Ge. 1.0) Then
       Call stop_program('can not converge in RhizoIn.')
    Endif
    Return
20  Call stop_program('ERROR in <Rhizo.in> -- Program terminated')
  End Subroutine RhizoIn
 !***************************************************************************
 ! ### numerical calculation of matric flux potential PHI (MFP)###
  Subroutine MFPC_table
    Use Typedef
    Use ParamData, Only : n_MFP
    Use SolData, Only : nmat, par
    Use WatFun, Only : hmax_MFP,hTab_MFP,MFPTab,FKP,hmin_MFP,hcheck,mfpcheck,hnewcheck,Fmfp_soiltab,Fh_from_mfp_soiltab,dexp_MFP
    Use RhizoData, Only: lRhizo
    Use MPIutils, Only: stop_program
    Use Doussan, Only: logspace
    Implicit None

    Integer(ap):: k,i
    real(dp):: y(n_MFP)

    If(lRhizo) Then
       Call stop_program('MFPC_table can not be use in rhizosphere modelling')
    Endif
	
!initialisation	
    Call logspace(hmax_MFP,hmin_MFP,n_MFP,y)
    hTab_MFP=-y
    dexp_MFP=(log10(abs(hmin_MFP))-log10(abs(hmax_MFP)))/(n_MFP-1)!increment of htab_MFP

    Do  k=1,nMat
       MFPTab(n_MFP,k)=0.0_dp
    End Do
!vector building	
    Do  i=1,n_MFP-1
       Do  k=1,nMat
          MFPTab(n_MFP-i,k)=MFPTab(n_MFP+1-i,k)+(FKP(hTab_MFP(n_MFP-i),par(:,k),i)+FKP(hTab_MFP(n_MFP+1-i),par(:,k),i))*(hTab_MFP(n_MFP-i)-hTab_MFP(n_MFP+1-i))/2.0_dp
       End Do
    End Do
!	%check
    Do  i=1,n_MFP
      Do  k=1,nMat
        hnewcheck(i)=hTab_MFP(i)-0.5
        mfpcheck(i,k)=Fmfp_soiltab(hnewcheck(i),k,par(:,k))
        hcheck(i)=Fh_from_mfp_soiltab(mfpcheck(i,k),k)
      End Do
    End Do
    ! hTab_MFP(1501)=-15010.0_dp
    ! Do  k=1,nMat
       ! MFPTab(1501,k)=0.0_dp
    ! End Do
    ! Do  i=1,1500
       ! hTab_MFP(1501-i)=-15010.0_dp+10.0_dp*i
       ! Do  k=1,nMat
          ! MFPTab(1501-i,k)=MFPTab(1502-i,k)+(FKP(hTab_MFP(1501-i),par(:,k),i)+FKP(hTab_MFP(1502-i),par(:,k),i))*10.0_dp/2.0_dp
       ! End Do
    !End Do
    Return
  End Subroutine MFPC_Table
  !***************************************************************************
  !> ### reads root related inputs ###
  Subroutine RootIn(t,ipl)
    Use Typedef
    Use tmctrl
    Use ParamData
    Use EnviData
    Use TempData
    Use RootData
    Use PlntData
    Use GridData
    Use GeoData
    Use ConData
    Use StrData
    Use MPIutils, Only: stop_program
    Use RootTyp, Only: RunRootTip, CheckSize
    Implicit None

    INTEGER(ap),INTENT(in):: ipl
    Integer(ap) :: ifg(maxemg)=0,dumI(maxemg)=0
    Integer(ap) :: iestbl,i,j,irec,igrow,ifive,timeRT,linlim,trash
    Real(dp) :: t
    Real(dp) :: maxZ,dumR(maxemg),dumR1(maxemg,mxpnts)=0,dumR2(maxemg,mxpnts)=0
    Character :: infile*12

!> #### if RootTyp is used, reads param.txt ####
    If (lrrt) Then
       Open (UNIT=10,FILE='in/param.txt',STATUS='OLD',ERR=50)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) timeRT
       !> no multiple roots with RootTyp
       nplant=1
       xplant(1)=0
       yplant(1)=0
       Call RunRootTip(timeRT,ipl)
       t=timeRT !implicit change of type INT->real
       Write(15,'(/''RootTyp was used for generating the root system'')')
       Write(15,'(/''Initial Time for simulation is '',1pE11.4)',advance='no') t
!> #### reads RootSys ####
    Elseif (lrrs) Then!Read RootSys and root.in --> root growth with Somma or RootBox
       WRITE (infile, "(A10,I1)") 'in/RootSys',ipl!only valid for less than 10 plants!
       If (nplant.GT.1) Then
            Open (UNIT=10,FILE=infile,STATUS='OLD',ERR=60)
       Else
            Open (UNIT=10,FILE='in/RootSys',STATUS='OLD',ERR=60)
       End If
       Write(15,'(/''* Root system imported from file '',1A12)',advance='no') infile
       Read (10,*,ERR=40)
       !> 1. general parameters
       Read (10,*,ERR=40) t
       Write(15,'(/''* Initial Time for simulation is '',1pE11.4)',advance='no') t
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) 
       Write(15,'(/''* Number of different root systems : '',I5)',advance='no') nplant
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) trash,xplant(ipl),yplant(ipl)
       Write(15,'(/''* Seed locations for plant '',1I5, '' is '', 2F9.3)',advance='no') ipl,xplant(ipl),yplant(ipl)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) mroot(ipl),mshoot(ipl),LA(ipl)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) sAvg,cAvg
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) naxes(ipl)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) nbr(ipl)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) nrec(ipl)
       If ((nbr(ipl).Lt.1).Or.(nrec(ipl).Lt.1))  Goto 20
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       !> 2. segment information
       Do i=1,nrec(ipl)
          Read (10,*,ERR=40) irec,xs(irec,ipl),ys(irec,ipl),zs(irec,ipl),irecpr(irec,ipl),ordseg(irec,ipl),&
               ibrseg(irec,ipl),seglen(irec,ipl),segsur(irec,ipl),segmas(irec,ipl)
          Read (10,*,ERR=40) timorg(irec,ipl)
          !Cross Section of each segment - later needed for flow velocities & particle tracker
          segrad(irec,ipl) = segsur(irec,ipl)/2._dp/pi/seglen(irec,ipl)
          crossSectionSeg(irec,ipl) = pi*segrad(irec,ipl)**2
          segvol(irec,ipl) =  segrad(irec,ipl)**2*pi*seglen(irec,ipl) !root volume per root segment
          tot_len(ipl) = tot_len(ipl) + seglen(irec,ipl) !total root length
          vol_root(ipl) = vol_root(ipl) +  segvol(irec,ipl) !total root volume
       End Do
       If (Maxval(ordseg).Eq.13) Then!Recognizes old roots created by RootTyp (Couvreur nov 2010)
          maizeroottyp=.True.
       Elseif (Maxval(ordseg).Eq.5) Then
          loliumroottyp=.True.
       Elseif (Maxval(ordseg).Eq.20) Then
          wheatroottyp=.True.
       Endif
       ibrseg(0,ipl)=0
       segsur(0,ipl) =0.
       If ((zs(1,ipl).Gt.zgrid(1)).Or.(zs(1,ipl).Lt.zgrid(npt))) Then
          Write (*,'(//'' Inconsistency -- Position of first root segment not within the spatial domain'')')
          Write (*,'('' as defined in  < nodes.in >.''/)')
          Call stop_program('Program terminated.')
       Endif
       If (.Not.(continu)) Then
       Call CheckSize(ipl)
       Endif
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40) ngrow(ipl)
       If (ngrow(ipl).Lt.1) Goto 20
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       Read (10,*,ERR=40)
       !> 3. tip information
8      Read (10,*,ERR=40) igrow,xg(igrow,ipl),yg(igrow,ipl),zg(igrow,ipl),irecsg(igrow,ipl),ordgrw(igrow,ipl),ibrgrw(igrow,ipl),&
            brlgth(igrow,ipl), iaxis(igrow,ipl)
       Read (10,*,ERR=40) ovrtime(igrow,ipl),nestbl(igrow,ipl)
       If (nestbl(igrow,ipl).Gt.0) Then
          ifive=0
81        linlim=Min(ifive+5,nestbl(igrow,ipl))
          Read (10,*,ERR=40) (timest(igrow,iestbl),iestbl=ifive+1,linlim)
          ifive=ifive+5
          If (linlim.Lt.nestbl(igrow,ipl)) Goto 81
       Endif
       maxZ=Maxval(ZGrid)
       If (igrow.Lt.ngrow(ipl)) Goto 8
       zg(1:ngrow(ipl),ipl)=zg(1:ngrow(ipl),ipl)-maxZ !??? maxZ?
       Close (10)
    End If

!> #### if Somma growth, reads root.in  ####
    If (lSomma_growth) Then
       Open(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
       Write(15,'(//''++ Root System ++'')',advance='no')
       Write(15,'(/''-----------------'')',advance='no')
       IF (nplant.GT.1) Then
            WRITE (infile,"(A10,I1)") 'in/root.in',ipl
       Else
            WRITE (infile,"(A10,I1)") 'in/root.in'
       End If

       Open (UNIT=10,FILE=infile,STATUS='OLD',ERR=10)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) naxemg(ipl)
       If (naxemg(ipl).Gt.maxemg) Call stop_program('Input from  < root.in >  exceeds parameter "maxemg". Program terminated.')
       Read (10,*,ERR=20)
       Do i=1,naxemg(ipl)
          Read (10,*,ERR=20) tnewax(i,ipl),nnewax(i,ipl), dumR(i)
       End Do
       naxtot(ipl)=naxes(ipl)
       ifg(1:naxtot(ipl))=1 !initial axes aways belong to axis group #1
       Do i=1,naxemg(ipl)
          ifg(naxtot(ipl)+1:naxtot(ipl)+nnewax(i,ipl))=i
          naxtot(ipl)=naxtot(ipl)+nnewax(i,ipl)
       End Do
       Do i=1,naxtot(ipl)
          inaxs(i,ipl)=dumR(ifg(i))/180._dp*pi
       End Do
       dumR = 0.
       If (naxtot(ipl).Gt.maxemg) Call stop_program('Input from  < root.in >  exceeds parameter "maxemg". Program terminated.')
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (dumR(i),i=1,naxemg(ipl))
       Do i=1,naxtot(ipl)
          geoaxs(i,ipl)=dumR(ifg(i))
       End Do
       dumR = 0.
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (dumI(i),i=1,naxemg(ipl))
       Do i=1,naxtot(ipl)
          nangax(i,ipl)=dumI(ifg(i))
       End Do
       Read (10,*,ERR=20)
       Do  i=1,naxemg(ipl)
          Read (10,*,ERR=20) (dumR1(i,j),dumR2(i,j), j=1, dumI(i))
       End Do
       Do i=1,naxtot(ipl)
          tempax(i,1:nangax(i,ipl),ipl) = dumR1(ifg(i),1:nangax(i,ipl))
          angaxs(i,1:nangax(i,ipl),ipl) = dumR2(ifg(i),1:nangax(i,ipl))/180._dp*pi
       End Do
       dumI = 0
       dumR1 = 0.
       dumR2 = 0.
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) geolat(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) nanglt(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (templt(i,ipl),anglat(i,ipl),i=1,nanglt(ipl))
       Do j=1,nanglt(ipl)
          anglat(j,ipl)=anglat(j,ipl)/180._dp*pi
       End Do
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) norder(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) nBigLat(ipl)
       If (norder(ipl).Gt.maxord) Call stop_program('Input from  < root.in >  exceeds parameter "maxord". Program terminated.')
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) 
       Read (10,*,ERR=20) (nVch(i,ipl),i=1,norder(ipl)+1)
       Read (10,*,ERR=20)
       Do i=1,(norder(ipl)+1)
           Read (10,*,ERR=20) (ageVch(i,j,ipl),Vch(i,j,ipl),j=1,nVch(i,ipl))
       End Do
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) SpWgt(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (nMPLch(i,ipl),i=1,norder(ipl))
       Read (10,*,ERR=20)
       Do  i=1,norder(ipl)
          Read (10,*,ERR=20) (sMPLch(i,j,ipl),MPLch(i,j,ipl),j=1,nMPLch(i,ipl))
       End Do
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) l_conduc,lLandl
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) condMP(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) l_overburden,lBengough1,lSoilStrength
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (strsen(i,ipl),i=1,norder(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (rdmang(i,ipl),i=1,norder(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) tempermin(ipl),topt(ipl),tempermax(ipl)
       trange(ipl)=tempermax(ipl)-tempermin(ipl)
       tmid(ipl)=(tempermin(ipl)+tempermax(ipl))/2._dp
       If (topt(ipl).Lt.tmid(ipl)) Then
          expo(ipl)=Log(.5_dp)/Log((topt(ipl)-tempermin(ipl))/trange(ipl))
       Else
          expo(ipl)=Log(.5_dp)/Log((topt(ipl)-tempermax(ipl))/(-trange(ipl)))
       Endif
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) ltoxi
       Read (10,*,ERR=20)
       If (ltoxi) Then
          Read (10,*,ERR=20) cmin(ipl),coptmi(ipl),coptma(ipl),cmax(ipl)
       Else
          Read (10,*,ERR=20)
       Endif
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (brlmax(i,ipl),i=1,norder(ipl)+1)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (brspac(i,ipl),i=1,norder(ipl)-1)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (brnang(i,ipl),i=1,norder(ipl)-1)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (dtbrch(i,ipl),i=1,norder(ipl)-1)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) l_secrad
       Read (10,*,ERR=20)
       If(l_secrad) Then
          Read (10,*,ERR=20) (f_rad(i,ipl),i=1,norder(ipl))
       Else
          Read (10,*,ERR=20)
       EndIf
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) lHydrotrop, lHydropattern, lXerobranch
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) hdrsens, Xerosens   
       Close (10)
       Do i=1,norder(ipl)-1
          brnang(i,ipl)=brnang(i,ipl)/180._dp*pi
          rdmang(i,ipl)=rdmang(i,ipl)/180._dp*pi
       End Do
       rdmang(norder(ipl),ipl)=rdmang(norder(ipl),ipl)/180._dp*pi
       ELSEIF (lRootBox_growth) Then
       Open(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
       Write(15,'(//''++ Root System ++'')',advance='no')
       Write(15,'(/''-----------------'')',advance='no')

       WRITE (infile,"(A10,I1)") 'in/rootBox.in',ipl
       If (nplant.GT.1) Then
            Open (UNIT=10,FILE=infile,STATUS='OLD',ERR=10)
       Else
            Open (UNIT=10,FILE='in/rootBox.in',STATUS='OLD',ERR=10)
       End If
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)   
       Read (10,*,ERR=20) norder(ipl)
       Read (10,*,ERR=20)   
       Read (10,*,ERR=20) nBigLat(ipl)
       If (norder(ipl).Gt.maxord) Call stop_program('Input from  < rootBox.in >  exceeds parameter "maxord". Program terminated.')
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (brlmax(i,ipl), i=1,norder(ipl)+1)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) diffnum(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (maxlast(i,ipl),i=1,diffnum(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (numlast(i,ipl),i=1,diffnum(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (Vch(i,1,ipl), i=1,norder(ipl)+1)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (sg(i,ipl), i=1,norder(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (rdmang(i,ipl), i=1,norder(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (rootrad(i,ipl), i=1,norder(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (lb(i,ipl), i=1,norder(ipl)-1)   
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (brspac(i,ipl), i=1,norder(ipl)-1) 
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (brnang(i,ipl), i=1,norder(ipl)-1)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) naxemg(ipl)
       If (naxemg(ipl).Gt.maxemg) Call stop_program('Input from  < root.in >  exceeds parameter "maxemg". Program terminated.')
       Read (10,*,ERR=20)
       Do i=1,naxemg(ipl)
          Read (10,*,ERR=20) tnewax(i,ipl),nnewax(i,ipl), dumR(i)
       End Do
       naxtot(ipl)=naxes(ipl)
       ifg(1:naxtot(ipl))=1 !initial axes aways belong to axis group #1
       Do i=1,naxemg(ipl)
          ifg(naxtot(ipl)+1:naxtot(ipl)+nnewax(i,ipl))=i
          naxtot(ipl)=naxtot(ipl)+nnewax(i,ipl)
       End Do
       Do i=1,naxtot(ipl)
          inaxs(i,ipl)=dumR(ifg(i))/180._dp*pi
       End Do
       dumR = 0.
       If (naxtot(ipl).Gt.maxemg) Call stop_program('Input from  < root.in >  exceeds parameter "maxemg". Program terminated.')
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (dumI(i),i=1,naxemg(ipl))
       Do i=1,naxtot(ipl)
          nangax(i,ipl)=dumI(ifg(i))
       End Do
       Read (10,*,ERR=20)
       Do  i=1,naxemg(ipl)
          Read (10,*,ERR=20) (dumR1(i,j),dumR2(i,j), j=1, dumI(i))
       End Do
       Do i=1,naxtot(ipl)
          tempax(i,1:nangax(i,ipl),ipl) = dumR1(ifg(i),1:nangax(i,ipl))
          angaxs(i,1:nangax(i,ipl),ipl) = dumR2(ifg(i),1:nangax(i,ipl))/180._dp*pi
       End Do
       dumI = 0
       dumR1 = 0.
       dumR2 = 0.
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) nanglt(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (templt(i,ipl),anglat(i,ipl),i=1,nanglt(ipl))
       Do j=1,nanglt(ipl)
          anglat(j,ipl)=anglat(j,ipl)/180._dp*pi
       End Do
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (nMPLch(i,ipl),i=1,norder(ipl))
       Read (10,*,ERR=20)
       Do  i=1,norder(ipl)
          Read (10,*,ERR=20) (sMPLch(i,j,ipl),MPLch(i,j,ipl),j=1,nMPLch(i,ipl))
       End Do
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) l_conduc
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) condMP(ipl)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) l_overburden
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) (strsen(i,ipl),i=1,norder(ipl))
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) tempermin(ipl),topt(ipl),tempermax(ipl)
       trange(ipl)=tempermax(ipl)-tempermin(ipl)
       tmid(ipl)=(tempermin(ipl)+tempermax(ipl))/2._dp
       If (topt(ipl).Lt.tmid(ipl)) Then
          expo(ipl)=Log(.5_dp)/Log((topt(ipl)-tempermin(ipl))/trange(ipl))
       Else
          expo(ipl)=Log(.5_dp)/Log((topt(ipl)-tempermax(ipl))/(-trange(ipl)))
       Endif
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) ltoxi
       Read (10,*,ERR=20)
       If (ltoxi) Then
          Read (10,*,ERR=20) cmin(ipl),coptmi(ipl),coptma(ipl),cmax(ipl)
       Else
          Read (10,*,ERR=20)
       Endif
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20)
       Read (10,*,ERR=20) l_secrad
       Read (10,*,ERR=20)
       If(l_secrad) Read (10,*,ERR=20) (f_rad(i,ipl),i=1,norder(ipl))
       Close (10)
       Do i=1,norder(ipl)-1
          brnang(i,ipl)=brnang(i,ipl)/180._dp*pi
          rdmang(i,ipl)=rdmang(i,ipl)/180._dp*pi
       End Do
       rdmang(norder(ipl),ipl)=rdmang(norder(ipl),ipl)/180._dp*pi
    Endif

    If (.Not.lretry) Then
       Write (15,'(///'' Simulation starts at time = '',F8.3,''.'')') t
       Write (15,'(/''The root system consists now of '',I5,'' branch(es)[including '',I3,'' axis(es)]'')')nbr,naxes
       Write (15,'('' or of a total of '',I5,'' segment(s) and '',I5,'' growing branch tip(s).'')')nrec,ngrow
       Write (15,'(/'' Total root mass is '',F9.3,'', total shoot mass '',F9.3,''.'')')mroot,mshoot
       Write (15,'('' Leaf area is '',F9.3,''.'')')LA
    Endif
    Return

10  Call stop_program('Data inconsistency in  < root.in or rootBox.in > not found -- program terminated.')
20  Call stop_program('Data inconsistency in  < root.in or rootBox.in>  -- program terminated.')
    !    30 call stop_program('Root input file does not exist -- program terminated.')
40  Call stop_program('Data error in root system file  -- program terminated.')
50  Call stop_program('File param.txt not found --program terminated')
60  Call stop_program('File RootSys not found --program terminated')
  End Subroutine RootIn
  !**************************************************************************************************
  !> ### plant / shoot relevant inputs ###
  Subroutine PlntIn(t)
    Use typedef
    Use ParamData, Only: lChem,Pi,lretry,last_out, lclimate
    Use PlntData
    Use RootData,Only : lDou,lCou,lFed,lSomma_growth,lRootBox_growth,lno_RWU,lCalloc,lno_Archi,ltoxi,timeobs,nplant,fac_plant
    Use DoussanMat, Only: stresfun
    Use tmctrl
    Use EnviData
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: i,j,funBC,aa,n,nt,typed(mxbcch)
    Integer(ap), Allocatable, Dimension(:) :: typeBCrt
    Real(dp) :: t,BCd(mxbcch)
    Real(dp), Allocatable, Dimension(:) :: tBCrt1,dummy

    ALLOCATE(fac_plant(1:nplant))
    !read BCRoot.in
    !IF  (((.NOT.lno_RWU).AND.(.NOT.lSign_new))) THEN!BCroot.in now read as long as the RWU process and no Carbon allocation are simulated, Tpot.in erased
    ! file with boundary conditions given for root in therms of PH or of flux
    If  (.Not.lno_RWU) Then

!> #### reads BCroot.in ####
       Open (Unit=10,FILE='in/BCroot.in',STATUS='OLD',ERR=3001)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) funBC,(fac_plant(i),i=1,nplant)
       Read (10,*)
       Read (10,*) nBCr

       If (nBCr.Gt.mxBcCh) Then
          Call stop_program(' the number of boundary conditions for root exceeded the maximum. Simulation stopped.')
       Endif
       Read (10,*)

       If(lclimate)Then
          nBCr=nclimaticdata
          If (funBC.Le.2) Then !> free format or constant value
              Do i=1,nclimaticdata
                typeBCr(i) = 2 !Constant value
                tBCr(i)=time_climate(i)
                BCroot(i) = T_pot(i)
                !   print*, tBCr(i), BCroot(i), 'time and value'
              Enddo
          Else
              Allocate (tBCrt1(1:nBCr+1))
              Allocate (dummy(1:nBCr+1))
              Allocate (typeBCrt(1:nBCr+1))     
              
              Do i=1,nclimaticdata
                typeBCrt(i) = 2 !Constant value
                tBCrt1(i)=time_climate(i)
                dummy(i) = T_pot(i)
                !   print*, tBCr(i), BCroot(i), 'time and value'
              Enddo
             
             !> if transpiration is not constant, interpolate daily values
             nt = 0
             Do i=2,nBCr
                aa = FLOOR((tBCrt1(i)-tBCrt1(i-1))+0.5)
                Do j = 1,aa
                   nt = nt+1
                   If (funBC .LT. 3) Then
                   BCd(nt) = dummy(i-1) + j*(dummy(i)-dummy(i-1))/(tBCrt1(i)-tBCrt1(i-1))
                   Else
                   BCd(nt) = dummy(i-1) !if sinusoidal, no interpolation 
                   Endif
                   Print*, 'Bcd, nt', BCd(nt), nt
                   typed(nt) = typeBCrt(i-1)
                End Do
             End Do
             n = 1
             Do i=1,nt
                Do j=1,25 !> sub-daily interpolation steps (25/day)
                   tBCr(n) = tBCrt1(1)+(i-1)+Real(j)/25
                   typeBCr(n) = typed(i)
                   If (funBC .Eq. 3) Then
                      BCroot(n) = BCd(i)*(Sin(2*pi*(tBCr(n)-0.25))+1)
                      n = n+1
                   Else If (funBC .Eq. 4) Then
                      BCroot(n) = BCd(i)*pi*Sin((tBCr(n)-0.25)*2*pi)
                      If (BCroot(n).Lt.0) BCroot(n)=0.0
                      n = n+1
                   Else
                      Call stop_program(' ... ERROR in input file BCroot.in ... (funBC id#)')
                   End If
                End Do
             End Do
             nBCr = n-1
          Endif   
       Elseif (funBC.Le.4) Then 
          If (funBC.Le.2) Then !> free format or constant value
             Read (10,*) (tBCr(i),typeBCr(i),BCroot(i),i=1,nBCr)
             If (lFed) Then
                Do i=1,nBCr
                   If (typeBCr(i).Ne.2) Call stop_program('With Feddes RWU model, the upper BC at the root collar should be of flux type')
                Enddo
             Endif
          Else
             Allocate (tBCrt1(1:nBCr+1))
             Allocate (dummy(1:nBCr+1))
             Allocate (typeBCrt(1:nBCr+1))
             Read (10,*) (tBCrt1(i),typeBCrt(i),dummy(i),i=1,nBCr)
             If (lFed) Then
                Do i=1,nBCr
                   If (typeBCrt(i).Ne.2) Call stop_program('With Feddes RWU model, the upper BC at the root collar should be of flux type')
                Enddo
             Endif
             !> if runtime is longer than last tBCrin, the last BCr will be used
             If (tBCrt1(nBCr).Lt.tMax) Then
                nBCr = nBCr+1
                tBCrt1(nBCr) = tMax
                dummy(nBCr) = dummy(nBCr-1)
             End If
             !> if transpiration is not constant, interpolate daily values
             nt = 0
             Do i=2,nBCr
                aa = FLOOR((tBCrt1(i)-tBCrt1(i-1))+0.5)
                Do j = 1,aa
                   nt = nt+1
                   If (funBC .LT. 3) Then
                        BCd(nt) = dummy(i-1) + j*(dummy(i)-dummy(i-1))/(tBCrt1(i)-tBCrt1(i-1))
                   Else
                        BCd(nt) = dummy(i-1) !if sinusoidal, no interpolation 
                   Endif                   
                   typed(nt) = typeBCrt(i-1)
                End Do
             End Do

             n = 1
             Do i=1,nt
                Do j=1,25 !> sub-daily interpolation steps (25/day)
                   tBCr(n) = tBCrt1(1)+(i-1)+Real(j)/25
                   typeBCr(n) = typed(i)
                   If (funBC .Eq. 3) Then
                      BCroot(n) = BCd(i)*(Sin(2*pi*(tBCr(n)-0.25))+1)
                      n = n+1
                   Else If (funBC .Eq. 4) Then
                      BCroot(n) = BCd(i)*pi*Sin((tBCr(n)-0.25)*2*pi)
                      If (BCroot(n).Lt.0) BCroot(n)=0.0
                      n = n+1
                   Else
                      Call stop_program(' ... ERROR in input file BCroot.in ... (funBC id#)')
                   End If
                End Do
             End Do
             nBCr = n-1
          Endif
       Endif
!	   print *, 'nt', nt,'tBCrt1(1)',tBCrt1(1),'n',n
!       Write (*,'(/'' BC-data for root found.'')'))
       Close (10)

       If (lno_Archi) Then
          If (lretry) Then
             t=tOut(last_out)
             Write (*,'(///'' Simulation continues at time = '',F8.3,''.'')') t
          Elseif (lcou)Then
             t=timeObs(1)
          Else
             t=(tBCr(1)) !First time step is the first time of BCroot.in when RootSys or RootTyp architectures are not used
          Endif
       Endif
       If (lretry) Then
           t=tOut(last_out)
       Endif          
       If ((tBCR(1).Lt.t).And.(.Not.lretry)) Then
          Print *,'tBCR',tBCRt1(1),tBCR(1),tBCR(2),tBCR(3),tBCR(4),tBCR(5),'t',t
          Call stop_program(' initial time for root BC lower than initial simulation time')
       Endif
    Endif
    !> #### if shoot growth is considered, reads plant.in ####
    If (lCalloc) Then
       Open (Unit=10,FILE='in/plant.in',STATUS='OLD',ERR=11)
       !    --- all plant input is expressed as piecewise linear functions ---
       ! no-stress potential transpiration rate per leaf area as a function of time:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) ntTpLA
       If (ntTpLA.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (tTpLA(i),TpLAc(i),i=1,ntTpLA)
       !> no-stress water use efficiency (dry mass gained per water mass transpired)
       !> as f(time):
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) ntW
       If (ntW.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (tW(i),Wc(i),i=1,ntW)
       !> no-stress root/shoot ratio as a f(time):
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) ntRSR
       If (ntRSR.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (tRSR(i),RSRc(i),i=1,ntRSR)
       !> transpiration reduction factor as a function of relative stress due soil strength:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) nfTpLA
       If (nfTpLA.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (sfTpLA(i),fTpLAc(i),i=1,nfTpLA)
       !> water use efficiency factor as a function of relative stress due soil strength:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) nsfW
       If (nsfW.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (sfW(i),fWc(i),i=1,nsfW)
       !> root/shoot-ratio factor as a function of relative stress due soil strength:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) nsfRSR
       If (nsfRSR.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (sfRSR(i),fRSRc(i),i=1,nsfRSR)
       !> relative stress as a function of soil strength:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) ns
       If (ns.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (sc(i),rsc(i),i=1,ns)
       !> transpiration reduction factor as a function of relative stress due solute concentration:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) ncTpLA
       If (ncTpLA.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (scTpLA(i),cTpLAc(i),i=1,ncTpLA)
       !> water use efficiency factor as a function of relative stress due solute concentration:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) nscW
       If (nscW.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (scW(i),cWc(i),i=1,nscW)
       !> root/shoot-ratio factor as a function of relative stress due solute concentration:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) nscRSR
       If (nscRSR.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (scRSR(i),cRSRc(i),i=1,nscRSR)
       !> relative stress as a function of solute concentration:
       Read (10,*)
       Read (10,*)
       Read (10,*)
       Read (10,*) ncnc
       If (ncnc.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (cncp(i),rscnc(i),i=1,ncnc)
       Read (10,*)
       Read (10,*)
       Read (10,*)
       !> leaf area per dry shoot mass as a function of time:
       Read (10,*) ntLA
       If (ntLA.Gt.mxBcCh) Goto 99
       Read (10,*)
       Read (10,*) (tLA(i),LAc(i),i=1,ntLA)
       Close (10)
    End If
    !> Options declaration; which processes are modeled
    If ((lSomma_growth).OR.(lRootBox_growth)) Then
       If (.Not.lno_RWU) Then
          If (lChem) Then
             Write (15,'(/'' Interaction between root water - solute uptake / root growth'')')
             Write (15,'(''     and soil water content / soil strength will be simulated.'')')
          Else
             Write (15,'(/'' Interaction between root water uptake / root growth'')')
             Write (15,'(''     and soil water content / soil strength will be simulated.'')')
          End If
          If (lCalloc) Then
             Write (15,'(/'' Shoot growth and assimilate distribution will be considered. '')')
          Else
             Write (15,'(/'' Shoot growth and assimilate distribution will not be considered. '')')
          Endif
       Else
          Write (15,'('' -- will simulate root growth as affected by soil strength only.'')')
          Write (15,'(/'' Water and solute uptake, shoot growth and assimilate      '')')
          Write (15,'(''     distribution will not be considered.'')')
       Endif
       If (ltoxi) Then
          If (lChem) Then
             Write (15,'(//'' Nutrient deficiency and/or ion toxicity effects will be included'')')
             Write (15,'(''     in the simulation. '')')
          Else
             Write (15,'(//'' Nutrient deficiency and/or ion toxicity effects will be included'')')
             Write (15,'(''     in the simulation using initial concentration.'')')
          End If
       Else
          Write (15,'(//'' Nutrient deficiency and/or ion toxicity effects will not be '')')
          Write (15,'(''     included in the simulation. '')')
       Endif
    Endif
    If (lFed) Then
       If (stresfun.Eq.2) Then
          If (p50.Lt.9.E+20_dp) Then
             Write (*,'(/'' The van Genuchten (1978) expression will be used as water'')')
             Write      (*,'(''     extraction function. '')')
          Else
             Write (*,'(/'' The van Genuchten (1978) expression will be used for the water'')')
             Write (*,'(''     extraction function, but osmotic potential effects will be'')')
             Write (*,'(''     neglected. '')')
          Endif
       Elseif (stresfun.Eq.1) Then
          Write (*,'(//'' The Feddes et al. (1978) expression will be used for the water'')')
          Write (*,'(''     extraction function. '')')
       Else
          Write (15,'(//'' The water extraction function will be directly proportional to the relative'')')
          Write (15,'(''     root length density. '')')
       Endif
    Elseif (lDou) Then
       Write (15,'(//'' The Doussan et al. (1998) expression will be used as water extraction function'')')
    Elseif (lCou) Then
       Write (15,'(//'' The Couvreur et al. (2012) expression will be used as water extraction function'')')
    Endif

    Return
99  Call stop_program('Input from  < plant.in >  exceeds parameter "mxBcCh". Program terminated.')
3001 Call stop_program(' File  <BCroot.in>  not found --')
11  Call stop_program(' File  < plant.in >  not found --')
  End Subroutine PlntIn
  !***************************************************************************************
  !> ### if temperature effect is considered, reads temp.in ###
  Subroutine TempIn
    Use Typedef
    Use TempData
    !Use RootData, Only: ltemp,lSign
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: ii, jj

    Open (Unit=10,FILE='in/temp.in',STATUS='OLD',ERR=100)
    !> soil temperature as a function of time and depth
    !> (piecewise linear function):
    Read (10,*)
    Read (10,*)
    Read (10,*)
    Read (10,*) nz_tempS,nt_tempS
    If (nz_tempS.Gt.mxdpth) Call stop_program('Input from  < temp.in >  exceeds parameter"mxdpth". Program terminated.')
    Read (10,*)
    Read (10,*)
    Read (10,*) (depth(ii),ii=1,nz_tempS)
    depth=abs(depth)
    If (Abs(depth(1)).Gt.1.E-30_dp) Then
       Write (*,'(//'' Inconsistency - Second depth value in  < temp.in >  corresponds to soil surface'')')
       Write (*,'('' and must be zero.''/)')
       Call stop_program('Program terminated.')
    Endif
    Read (10,*)
    Read (10,*)
    Read (10,*) 
    Do ii=1, nt_tempS
       Read (10,*) time_S(ii),(temtim(ii,jj),jj=1,nz_tempS)
    End Do
   ! Read (10,*)
   ! Read (10,*)
   ! Read (10,*)
   ! Read (10,*) nt_tempA,nt_presA,nt_presD
   ! Read (10,*)
   ! Read (10,*)
   ! Read (10,*) (time_TA(ii),T_atm(ii),ii=1,nt_TempA)
   ! Read (10,*)
   ! Read (10,*)
   ! Read (10,*) (time_PA(ii),P_atm(ii),ii=1,nt_presA)
   ! Read (10,*)
   ! Read (10,*)
   ! Read (10,*) (time_PD(ii),P_diff(ii),ii=1,nt_presD)
    Close(10)
    Return
100 Call stop_program(' File  < temp.in >  not found -- needed if temperature influence on root growth or signaling should be modeled!')
  End Subroutine TempIn
  !********************************************************************
    !> ### reads depthdecay.in ###
Subroutine ReadDepthDecay
  Use typedef
  Use SoluteRootMat, Only: ndepth,Ddepth,dfactor
  Use MPIutils, Only: stop_program
  Implicit None
  
  Integer(ap) :: ii


  Open (Unit=10, File='in/depthdecay.in',Status='Old',Err=10)
  Read (10,*)
  Read (10,*)
  Read (10,*)
  Read (10,*) ndepth  
  Read (10,*)
  Read (10,*)
  Do ii=1,ndepth
     Read (10,*) Ddepth(ii), dfactor(ii)
     if (Ddepth(ii).lt.0._dp) Ddepth(ii)=abs(Ddepth(ii))
  End Do

  Return

10 Call stop_program('File depthdecay.in not found -- needed for solute decay in soil')
End Subroutine ReadDepthDecay
  !********************************************************************
    !> ### defines the number of nodes located in one observation plane / probe window ###
  Subroutine NodebyProbe
    Use Typedef
    Use ObsData
    Use GridData
    Use DomData
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: ip,delta,in,imin,i,iz,ix,iy

    !> observation planes
    Do ip=1,npr
       If (Pt(ip)==3) Then !horiz. plane perp. to Z
          nodebyPr(ip)=nx*ny
          delta=Int((Abs(CrP(ip)-zmax))/dzgrid)
          imin=nodebyPr(ip)*delta+1
          in=imin
          i=0
          Do i=1,nodebyPr(ip)
             NodePr(ip,i)=in
             in=in+1
          Enddo
       Elseif (Pt(ip)==2) Then !vert. plane perp. to Y
          nodebyPr(ip)=nx*nz
          delta=Int((Abs(CrP(ip)-ymin))/dygrid)
          imin=delta*nx+1
          i=0
          in=imin
          Do iz=1,nz
             Do ix=1,nx
                i=i+1
                NodePr(ip,i)=in
                in=in+1
             Enddo
             in=in+(nx)*(ny-1)
          Enddo
       Elseif (Pt(ip)==1) Then ! vert plane perp to X
          nodebyPr(ip)=ny*nz
          delta=Int((Abs(CrP(ip)-xmin))/dxgrid)
          imin=delta+1
          i=0
          in=imin
          Do iz=1,nz
             Do iy=1,ny
                i=i+1
                NodePr(ip,i)=in
                in=in+nx
             Enddo
          Enddo
       Elseif (Pt(ip).Ne.4) Then
          Call stop_program('Wrong plane direction in Probes.in')
       Endif
       Write(*,*)'plane',ip,'has ',nodebyPr(ip),'nodes.'
       Write(*,*) Pt(ip),ip,nodebyPr(ip),delta,imin,in
       If (nodebyPr(ip)>1000) Then
          Call stop_program( 'too many elements in one observation plane (>1000)')
       Endif
    Enddo
    Return
  End Subroutine NodebyProbe
  !*********************************************************************
  !> ### Defines list with elements that are of mixed material ###
  Subroutine ListMacro
    Use typedef
    Use GridData
    Use SolData, Only: MatNum,l_elmMacro,par,nMat
    Use RootGrowthNeighb, Only: ElemNeighb
    Implicit None

    Integer(ap) :: iElm,i
    Integer(ap) :: NeighElm(1:6),cif,k,imat

    cif=0  !number of mixed material voxels
    ! identify macropore material
    Do k=1,nMat
       If(par(11,k).Lt.1E-5_dp) imat=k
    End Do
    ! identify mixed material elements
    Do iElm=1,nElm
       If(Any(MatNum(elmnod(:,iElm)).Eq.imat)) Then
          l_elmMacro(iElm) = .True.
          cif=cif+1
          Call ElemNeighb(iElm,NeighElm)
          Do i = 1,6
             If(NeighElm(i).Gt.0)Then
                If(Any(MatNum(elmnod(:,NeighElm(i))).Eq.imat)) Then
                   MacroList(iElm,i)=0
                Else
                   MacroList(iElm,i)=NeighElm(i)
                   n_neigh(iElm) = n_neigh(iElm)+1
                End If
             End If
          End Do
       End If
    End Do
    !print*, 'number of mixed material voxels', cif
  End Subroutine ListMacro

End Module Input
