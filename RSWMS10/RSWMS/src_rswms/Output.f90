!> \file Output.f90
!! \brief Module for output.

!> Ouput module
Module Output

Contains
  !> generates output for observation nodes
  Subroutine ObsIni    
    Use Typedef
    Use ObsData
    Use GridData, Only : xgrid,ygrid,zgrid
    Implicit None

    Integer(ap) :: ip,minl,i
    Character filename*13,form*3,form2*2,form1*1

    Do ip=1,nPr
       minl=Min(nodebyPr(ip),500)
       !> define file names
       Write (filename,'(A13)')'out/Probe.  '
       If (ip.Lt.10) Then
          Write (filename(11:11),'(I1)') ip
       Else
          Write (filename(11:12),'(I2)') ip
       Endif
       Open (UNIT=10,FILE=filename,STATUS='UNKNOWN')
       !if more than 500 nodes
       If (minl.Ne.nodebyPr(ip)) Then
          If (ip.Lt.10) Then
             Write (filename(12:12),'(A1)')'b'
          Else
             Write (filename(13:13),'(A1)')'b'
          Endif
          Open(UNIT=11,FILE=filename,STATUS='UNKNOWN')
       Endif
       ! write titles
       Write (10,'(/''Observation nodes for each time step.'')')
       If (Pt(ip)==1) Write(10,'(/''Cross-section perpendicular to X axis at X = '',1pE11.4)')xgrid(NodePr(ip,1))
       If (Pt(ip)==2) Write(10,'(/''Cross-section perpendicular to Y axis at Y = '',1pE11.4)')ygrid(NodePr(ip,1))
       If (Pt(ip)==3) Write(10,'(/''Cross-section  perpendicular to Z axis at Z = '',1pE11.4)')zgrid(NodePr(ip,1))
       If (Pt(ip)==4) Write(10,'(/''Probe location defined by the user'')')
       If (DistrP(ip).Eq.1)  Then
          Write(10,'(/''This probe contains'',I4,'' nodes'')') nodebyPr(ip)
          If (varP(ip)==2) Write (10,'(/''      Time   Averaged water content'')')
          If (varP(ip)==1) Write (10,'(/''      Time   Averaged water potential'')')
          If (varP(ip)==3) Write (10,'(/''      Time   Aver. water potential   Aver. water content'')')
       Else
          If (minl.Ne.nodebyPr(ip)) Then
             Write(10,'(/''This probe contains'',I4,'' nodes but only 500 are given in this file, the rest is in the file b.'')')nodebyPr(ip)
             Write(form,'(I3)')500
          Else
             Write(10,'(/''This plane contains'',I4,'' nodes'')') nodebyPr(ip)
             If (nodebyPr(ip).Gt.99) Then
                Write(form,'(I3)')nodebyPr(ip)!only valid for larger than 100 numbers!!!!
             Elseif (nodebyPr(ip).Gt.9) Then
                Write(form2,'(I2)')nodebyPr(ip)
             Else
                Write(form1,'(I1)')nodebyPr(ip)
             Endif
          Endif
          Write(10,'(/''Nodes ID'')')
          If (nodebyPr(ip).Gt.99) Then
             Write(10,'('//form//'(1X,I5))') (NodePr(ip,i),i=1,minl)
          Elseif (nodebyPr(ip).Gt.9) Then
             Write(10,'('//form2//'(1X,I5))') (NodePr(ip,i),i=1,minl)
          Else
             Write(10,'('//form1//'(1X,I5))') (NodePr(ip,i),i=1,minl)
          Endif
          If (varP(ip)==2) Write(10,'(/''      Time   Node water content'')')
          If (varP(ip)==1) Write(10,'(/''      Time   Node water potential'')')
          If (varP(ip)==3) Write(10,'(/''      Time   Node water potential (line 1) and node water content (line 2)'')')
          !if more than 500 nodes
          If (minl.Ne.nodebyPr(ip)) Then
             Write (11,'(/''Observation nodes for each time step.'')')
             If (Pt(ip)==1) Write(11,'(/''plane perpendicular to X axis at X = '',1pE11.4)')Crp(ip)
             If (Pt(ip)==2) Write(11,'(/''plane perpendicular to Y axis at Y = '',1pE11.4)')Crp(ip)
             If (Pt(ip)==3) Write(11,'(/''plane perpendicular to Z axis at Z = '',1pE11.4)')Crp(ip)
             If (Pt(ip)==4) Write(11,'(/''probe  location defined by the user'')')
             Write(11,'(/''This probe contains'',I4,'' nodes. Here are given the'',I4,'' nodes larger than 500'')') &
                  nodebyPr(ip),(NodebyPr(ip)-500)
             Write (11,'(/''Nodes ID'')')
             Write(form,'(I3)')(nodebyPr(ip)-500)
             Write(11,'('//form//'(1X,I5))') (NodePr(ip,i),i=501,nodebyPr(ip))
             If (varP(ip)==2) Write (11,'(/''      Time   Node water content'')')
             If (varP(ip)==1) Write (11,'(/''      Time   Node water potential'')')
             If (varP(ip)==3) Write (11,'(/''      Time   Node water potential (line 1) and node water content (line 2)'')')
          Endif
       Endif
    Enddo
  End Subroutine ObsIni
 !****************************************************************************************
 !> generates log file; output of log depends on modelled processes
  Subroutine WriteLog(t)
    Use Typedef
    Use PlntData
    Use RootData
    Use EnviData
    Use Doussanmat,Only : count_nodes,PHr,nplant,PHs,PHsri
    Use GridData,Only: HElm,SSF
    Use ParamData, Only: lPartUp
    Use SolData
    Implicit None

    Integer(ap) ::ipl
    Real(dp):: t
    Character :: na*10,file*8
    !> \param t current time
	
DO ipl=1,nplant 
    Write (file,'(A7)')'out/log'
    na='    n/a   '
       Write (file(8:8),'(I1)') ipl
       Open (UNIT=10,FILE=file,STATUS='OLD',POSITION='APPEND')
       If (lFed) Then
          Write(10,'(3(1pE10.3,1X),A10,1X,2(1pE10.3,1X),A10,1X,1pE10.3)') &
               t,Tpot(ipl),Tact(ipl),na,sAvg,cAvg,na,mroot(ipl)
       Elseif (lCalloc) Then 
          Write (10,'(12(1pE11.4,1X),I6)')&
               t,Tpot(ipl),Tact(ipl),grwfac,sAvg,mshoot(ipl),mroot(ipl),LA(ipl),PHr(1,ipl),&
               concol,TpLA,TpLA_pot,count_nodes
       Elseif (lDou .And. .Not.lCalloc) Then
          If (lSign .Or. lSign_inst) Then
             Write (10,'(1pE13.6,1X,3(1pE11.4,1X),I6,1X,9(1pE11.4,1X))')&
                  t,Tpot(ipl),Tact(ipl),sAvg,count_nodes,PHr(1,ipl),mcol,concol,&
                  msign_notrans,csign_notrans,res_t
          Elseif(lPartUp) Then
             Write (10,'(4(1pE17.10,1X),I6,1X,13(1pE15.8,1X))')&
                  t,Tpot(ipl),Tact(ipl),sAvg,count_nodes,PHr(1,ipl),mcol,root_mass,&
                  Sum(SoilSoluteMass),(mcol+root_mass+Sum(SoilSoluteMass)),&
                  Uptake,Cum_uptake,(mcol+root_mass+Cum_uptake)
          Else
             If ((nplant.LE.1).and.(lKdrop)) Then
               Write (10,'(1pE13.6,1X,3(1pE11.4,1X),I6,6X,3(1pE11.4,1X))')&
                  t,Tpot(ipl),Tact(ipl),sAvg,count_nodes,PHr(1,ipl),DOT_Product(PHs(1:nrec_m,1),SSF),DOT_Product(PHsri(1:nrec_m,1),SSF)
             Else
               Write (10,'(1pE13.6,1X,3(1pE11.4,1X),I6,6X,(1pE11.4,1X))')&
                  t,Tpot(ipl),Tact(ipl),sAvg,count_nodes,PHr(1,ipl)
             Endif
          Endif
       Elseif (lCou) Then
          If (lSign_new) Then
             Write (10,'(1pE13.6,1X,3(1pE11.4,1X),I6,1X,5(1pE11.4,1X))')&
                  t,Tpot(ipl),Tact(ipl),sAvg,count_nodes,PHcollar,Hseq,Krs,gs_t,conversion_t
          Else
             Write (10,'(3(1pE11.4,1X),I2,1X,1pE11.4,1X,1pE11.4,1X,1pE11.4,1X,1pE11.4)')&
                  t,Tpot(ipl),Tact(ipl),0,PHcollar,DOT_Product(HElm,SSF),Krs,Kcomp
          Endif
       Endif
 End Do
    Return
  End Subroutine WriteLog

  !****************************************************************************************
  !> generates soil grid output, format is equal to nodes.in
  Subroutine OutFEM(t,kout)
    Use GridData
    Use SolData
    Use WatFun
    Use RhizoData
    Implicit None

    Integer(ap), Intent(in)::kout
    Integer(ap):: corner(8),ie,m,ic,i
    Real(dp),Intent(in)::t
    Real(dp):: sum,sumc,sume,sumce,tht_(nPt),solt,tottht,totsol,tht
    Character file*14
    !> \param t current time
    !> \param kout current number for the FEM output

    Write (file,'(A11)')'out/outfem.'
    If (kout.Lt.10) Then
       Write (file(12:12),'(I1)') kout
    Elseif (kout.Lt.100) Then
       Write (file(12:13),'(I2)') kout
    Else
       Write (file(12:14),'(I3)') kout
    Endif
    Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
    ! calculate total water and solute volume within domain:
    sum=0.0_dp
    sumC=0.0_dp
    Do 1 iE=1,nElm!,2
       ! assign cuboid corner nodes:
       corner(1)=elmnod(1,iE)
       corner(2)=elmnod(2,iE)
       corner(3)=elmnod(3,iE)
       corner(4)=elmnod(4,iE)
       corner(5)=elmnod(5,iE)
       corner(6)=elmnod(6,iE)
       corner(7)=elmnod(7,iE)
       corner(8)=elmnod(8,iE)
       ! average cuboid water and solute content:
       sumE=0.0_dp
       sumCE=0.0_dp
       Do 11 ic=1,8
          i=corner(ic)
          M=MatNum(i)
          If (soiltab) Then
             tht=Fth_soiltab(hNew(i),M)
             solt=(Fth_soiltab(hNew(i),M)+ChPar(1,M)*ChPar(5,M))*Conc(i)
          Else
             tht=Fth(hNew(i),par(:,M))
             solt=(Fth(hNew(i),par(:,M))+ChPar(1,M)*ChPar(5,M))*Conc(i)
          Endif
          sumE=sumE+tht
          sumCe=sumCE+solt
11     Enddo
       sum=sum+sumE
       sumC=sumC+sumCE
1   Enddo
    tottht= sum*dxGrid*dyGrid*dzGrid/8._dp
    totsol=sumC*dxGrid*dyGrid*dzGrid/8._dp
    Write(8,'(/''Total Water Volume at Time  '',F12.4,1X,A5,'' is'',1pE15.8,     ''.'')')&
         t,TmUnit,tottht
    Write(8,'(/''Total Solute Volume at Time '',F12.4,1X,A5,'' is'',1pE15.8,     ''.'')')&
         t,TmUnit,totsol
    Write (8,'(/''Length Unit is '',A5)')LnUnit
    If (.Not. lRhizo) Then
       Write (8,'(/4X,''Node#'',2X,''Mater.#'',9X,''x'',12X,''y'',12X,''z'',12X,''h'',10X,''conc.'',8X,''theta'',8X,''wsink'',8X,''csink'')')
       If (soiltab) Then
          Do  i=1,nPt
             M=MatNum(i)
             tht_(i)=Fth_soiltab(hNew(i),M)
             Write (8,'(I7,5X,I2,6X,3(2X,1pE11.4),5(2X,3pE11.4))')i,M,xGrid(i),yGrid(i),zGrid(i),hNew(i),Conc(i),tht_(i),sink(i),csink(i)
          Enddo
       Else
          Do  i=1,nPt
             M=MatNum(i)
             tht_(i)=Fth(hNew(i),par(:,M))
             ! if (hNew(i).lt.(-1*10**(99))) Then
                ! hNew(i)=-1*10**(99)
             ! endif
             Write (8,'(I7,5X,I2,6X,3(2X,1pE11.4),2(2X,3pE11.4),3(2X,1pE11.4))')i,M,xGrid(i),yGrid(i),zGrid(i),hNew(i),Conc(i),tht_(i),sink(i),csink(i)
          Enddo
       Endif
    Else
       Write (8,'(/4X,''Node#'',2X,''Mater.#'',8X,''x'',9X,''y'',9X,''z'',11X,''h'',9X,''conc.'',7X,''theta'',7X,''wsink'',7X,''csink'',7X,''thetaTot'',7X,''cond'')')
       Do i=1,nPt
          Write (8,'(I7,5X,I2,6X,3(2X,1pE11.4),7(2X,1pE11.4))')i,M,xGrid(i),yGrid(i),zGrid(i),hNew(i),Conc(i),theta(i),sink(i),csink(i),&
               thetaTot(i), con(i)
       Enddo
    Endif
    Close (8)
    Return
  End Subroutine OutFEM
  !****************************************************************************************
  !> in case of a continous domain, saves transroot for the basic root nodes -> draw the root in matlab
  Subroutine OutTrans(kout,ipl)
    Use typedef
    Use DoussanMat, Only: transroot,transtip,nsub
    Use RootData, Only: nrec,ngrow
    Implicit None

    Integer(ap),Intent(in)::kout
    Integer(ap) :: i,ipl
    Character file*17    
    !> \param t current time
    !> \param kout current number for the FEM output

    Write (file,'(A14)')'out/transroot.'
    Write (file(13:13),'(I1)') ipl
    If (kout.Lt.10) Then
       Write (file(15:15),'(I1)') kout
    Elseif (kout.Lt.100) Then
       Write (file(15:16),'(I2)') kout
    Else
       Write (file(15:17),'(I3)') kout
    Endif
    Open (UNIT=444,FILE='out/Transroot.out',STATUS='UNKNOWN')
    Write (444,'(/''Basic root nodes translations for each plant.'')')
       Write (444,'(/'' '')')
       Write (444,'(/''Plant '',i1)') ipl
       Write (444,'(/'' X  Y (times)'')')
       Do i=1,nrec(ipl)
          Write (444,'(i3, i3)') transroot(i,1,nsub(i,ipl),ipl), transroot(i,2,nsub(i,ipl),ipl)
       End Do
       Write (444,'(/'' '')')
       Write (444,'(/'' '')')
       Do i=1,ngrow(ipl)
          Write (444,'(i3, i3)') transtip(i,1,nsub(i,ipl),ipl), transtip(i,2,nsub(i,ipl),ipl)
       End Do
    Close(444)
    Return
  End Subroutine OutTrans
  !*************************************************************************
  !> if appropriate, writes depths profiles 
  Subroutine Zprofiles(t)
    Use GridData
    Use SolData
    Use WatFun
    Implicit None

    Integer(ap):: iz,M,prof,ixy,i
    Real(dp),Intent(in) :: t
    Real(dp):: thz(1000),phz(1000),sz(1000)
    Character form*5
    !> \param t current time

    Write(form,'(I3)')(nx*ny)
    ! open files
    Open (UNIT=121,FILE='out/ProfileTH.out',STATUS='OLD',POSITION='APPEND')
    Open (UNIT=122,FILE='out/ProfilePH.out',STATUS='OLD',POSITION='APPEND')
    Open (UNIT=123,FILE='out/ProfileS.out',STATUS='OLD',POSITION='APPEND')
    ! calculate average/sum
    prof=0
    Do iz=1,nz
       thz(iz)=0.
       phz(iz)=0.
       sz(iz)=0.
       Do ixy=1,nx*ny
          i=ixy+prof
          M=MatNum(i)
          If (soiltab) Then
             thz(iz)=thz(iz)+Fth_soiltab(hNew(i),M)
          Else
             thz(iz)=thz(iz)+Fth(hNew(i),par(:,M))
          Endif
          phz(iz)=phz(iz)+hNew(i)
          sz(iz)=sz(iz)+sink(i)
       Enddo
       thz(iz)=thz(iz)/(nx*ny)
       phz(iz)=phz(iz)/(nx*ny)
       prof=prof+nx*ny
    Enddo
    Write (121,'(1X,F15.6)') t,(thz(i),i=1,nz)
    Write (122,'(1X,F15.6)') t,(phz(i),i=1,nz)
    Write (123,'(1X,F15.6)') t,(sz(i),i=1,nz)
    Close(121)
    Close(122)
    Close(123)
  End Subroutine Zprofiles
  !*************************************************************************
  !> output for observation probes
  Subroutine OutObsProbe(t)
    Use GridData
    Use SolData
    Use WatFun
    Use ObsData
    Implicit None

    Integer(ap)::ip,n,M,minl,i
    Real(dp),Intent(in) :: t
    Real(dp):: tht(1000),pht(1000),AverPH,AverTH
    Character file*13,form*5
    !> \param t current time

    Do ip=1,nPr
       ! calculations, theta and PH
       tht=0
       pht=0
       Do n=1,nodebyPr(ip)
          M=MatNum(NodePr(ip,n))
          If (soiltab) Then
             tht(n)=Fth_soiltab(hNew(NodePr(ip,n)),M)
          Else
             tht(n)=Fth(hNew(NodePr(ip,n)),par(:,M))
          Endif
          pht(n)=hNew(NodePr(ip,n))
       Enddo

       ! define file names
       Write (file,'(A13)')'out/Probe.   '
       If (ip.Lt.10) Then
          Write (file(11:11),'(I1)') ip
       Else
          Write (file(11:12),'(I2)') ip
       Endif
       Open (UNIT=10,FILE=file,STATUS='OLD',POSITION='APPEND')
       If (DistrP(ip).Eq.1) Then
          !> average of each plane
          AverTH=Sum(tht)/nodebyPr(ip)
          AverPH=Sum(pht)/nodebyPr(ip)
          If (VarP(ip)==1) Write (10,'(2(1X,F12.4))') t,AverPH
          If (VarP(ip)==2) Write (10,'(2(1X,F12.4))') t,AverTH
          If (VarP(ip)==3) Write (10,'(3(1X,F12.4))') t,AverPH,AverTH
          Close(10)
       Else
          ! complete distribution asked
          minl=Min(nodebyPr(ip),500)
          Write(form,'(I3)')(minl+1)
          If (VarP(ip)==1) Write (10,'('//form//'(1X,F12.4))') t,(pht(i),i=1,minl)
          If (VarP(ip)==2) Write (10,'('//form//'(1X,F12.4))') t,(tht(i),i=1,minl)
          If (VarP(ip)==3) Then
             Write (10,'('//form//'(1X,F12.4))') t,(pht(i),i=1,minl)
             Write (10,'('//form//'(1X,F12.4))') t,(tht(i),i=1,minl)
          Endif
          Close(10)
          !IF more than 500 nodes
          If (minl.Ne.nodebyPr(ip)) Then
             If (ip.Lt.10) Then
                Write (file(12:12),'(A1)')'b'
             Else
                Write (file(13:13),'(A1)')'b'
             Endif
             Open(UNIT=11,FILE=file,STATUS='OLD',POSITION='APPEND')
             Write(form,'(I3)')(nodebyPr(ip)-500)+1
             If (VarP(ip)==1) Write (11,'('//form//'(1X,F12.4))') t,(pht(i),i=501,nodebyPr(ip))
             If (VarP(ip)==2) Write (11,'('//form//'(1X,F12.4))') t,(tht(i),i=501,nodebyPr(ip))
             If (VarP(ip)==3) Then
                Write (11,'('//form//'(1X,F12.4))') t,(pht(i),i=501,nodebyPr(ip))
                Write (11,'('//form//'(1X,F12.4))') t,(tht(i),i=501,nodebyPr(ip))
             Endif
             Close(11)
          Endif
       Endif
    Enddo
    Return
  End Subroutine OutObsProbe
  !*************************************************************************
  !> in case of Somma or RootBox root growth, writes root output; format is the same as RootSys 
  Subroutine OutRoo(t,kOut,ipl)
    Use typedef
    Use ParamData, Only : pi
    Use RootData
    Use PlntData, Only : LA
    Implicit None

    Integer(ap) :: irec,igrow!,ifive,iestbl,linlim
    Integer(ap), Intent(in) :: kOut,ipl
    Real(dp), Intent(in) :: t
    Character outfile*16
	
    Write (outfile,"(A11,I1,A1)")'out/RootSys',ipl,'.'
    If (kout.Lt.10) Then
       Write (outfile(14:14),'(I1)') kout
    Elseif (kout.Lt.100) Then
       Write (outfile(14:15),'(I2)') kout
    Else
       Write (outfile(14:16),'(I3)') kout
    Endif
!before version 10, tehse files had only the time as name
!    Write (outfile,"(A4,I1)")'out/',ipl
!   Write (outfile(5:12),'(F7.3)')t
    Open (UNIT=8,FILE=outfile,STATUS='UNKNOWN')
    Write (8,'(''Time:'')')
    Write (8,*) t
    Write (8,*)
    Write (8,'(''Number of seeds'')')
    Write (8,*) nplant
    Write (8,*)
    Write (8,'(''ID, X and Y coordinates of the seeds (one per line)'')')
    Write (8,'(I5,3(1X,1pE9.2))') ipl,xplant(ipl),yplant(ipl)
    Write (8,*)
    Write (8,'(''Root DM, shoot DM, leaf area:'')')
    Write (8,*) mroot(ipl),mshoot(ipl),LA(ipl)
    Write (8,*)
    Write (8,'(''Average soil strength and solute concentration experienced by root system:'')')
    Write (8,*) sAvg,cAvg
    Write (8,*)
    Write (8,'(''Total # of axes:'')')
    Write (8,*) naxes(ipl)
    Write (8,*)
    Write (8,'(''Total # of branches, including axis(es):'')')
    Write (8,*) nbr(ipl)
    Write (8,*)
    Write (8,'(''Total # of segment records:'')')
    Write (8,*) nrec(ipl)
    Write (8,*)
    Write (8,'(''segID#'',6X,''x'',11X,''y'',11X,''z'',6X,''prev or '','' br#  length   surface  mass'')')
    Write (8,'(''origination time'')')
    ! write list of all segment records:
    Do irec=1,nrec(ipl)
       Write (8,'(I6,3(1X,1pE12.5),1X,I6,1X,I2,1X,I5,3(1X,1pE11.4))')&
            irec,xs(irec,ipl),ys(irec,ipl),zs(irec,ipl),irecpr(irec,ipl),ordseg(irec,ipl),ibrseg(irec,ipl),&
            seglen(irec,ipl),segsur(irec,ipl),segmas(irec,ipl)
       Write (8,'(1pE11.4)') timorg(irec,ipl)
    End Do
    Write (8,*)
    Write (8,'(''Total # of growing branch tips:'')')
    Write (8,*) ngrow(ipl)
    Write (8,*)
    Write (8,'(''tipID#'',4X,''xg'',10X,''yg'',10X,''zg'',6X,''sg.bhd.tp. '',''ord  br#  tot.br.lgth. axs#'')')
    Write (8,'(''overlength'',2X,''# of estblished points'')')

    Write (8,'(''time of establishing (-->)'')')
    ! write list of all growing tips:
    Do igrow=1,ngrow(ipl)
       Write(8,'(I6,3(1X,1pE11.4),1X,I6,6X,I2,1X,I5,1X,1pE11.4,3X,I3)')&
            igrow,xg(igrow,ipl),yg(igrow,ipl),zg(igrow,ipl),irecsg(igrow,ipl),ordgrw(igrow,ipl),&
            ibrgrw(igrow,ipl),brlgth(igrow,ipl),iaxis(igrow,ipl)
       Write (8,'(1pE11.4,1X,I5)')ovrtime(igrow,ipl),nestbl(igrow,ipl)
     !  If (nestbl(igrow,ipl).Gt.0) Then
     !     ifive=0
!51   !     linlim=Min(ifive+5,nestbl(igrow,ipl))
     !     Write (8,'(5(1X,1pE13.6))')(timest(igrow,iestbl),iestbl=ifive+1,linlim)
     !     ifive=ifive+5
     !     If (linlim.Lt.nestbl(igrow,ipl)) Goto 51
     !  Endif
    End Do
    Close(8)
    Return
  End Subroutine OutRoo
  !*********************************************************************
  !> writes a list of all root segments and their corresponding radial and axial
  !> fluxes, conductivities, etc.: OutRoot.XX 
  Subroutine OutDou(t,kout,ipl)
    Use typedef
    Use RootData, Only: lKdrop,TypeKdrop,irecpr,xs,ys,zs,segrad,nrec,ibrseg,&
        yplant,xplant
    Use DoussanMat, Only: Lr,Khr,PHs,PHr,PHsri,sinkR,axialRootFlow,Qi,Qd,Q_bc1,&
       veloRoot,PHbulk,PHsri
    Use ParamData, Only: pi
    Use SoluteRootMat, Only : segsorb,l_linSorb,l_freundSorb
    Implicit None
    
    Integer(ap) :: irec,kout,ipl
    Real(dp), Intent(in) :: t
    Character file*15

    Write (file,'(A12)')'out/outRoot.'
    !> \param t current time
    !> \param kout current number for the root output

       Write (file(11:11),'(I1)') ipl
       If (kout.Lt.10) Then
          Write (file(13:13),'(I1)') kout
       Elseif (kout.Lt.100) Then
          Write (file(13:14),'(I2)') kout
       Else
          Write (file(13:15),'(I3)') kout
       Endif
       Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
       Write (8,'(''Time:'')')
       Write (8,*) t
       Write (8,*)
       Write (8,'(''Total # of segment records:'')')
       Write (8,*) nrec(ipl)
       Write (8,*)
       Write (8,90)
90     Format(' segID       x          y         z        br#   prev     Lr       Kx         PHbulk       PHinter       PHxylem   radialRootFlow axialRootFlow      Qi           Qd           Q_bc         radius        veloRoot       segconc ')
       If(l_linSorb.Or.l_freundSorb)Then
          Do irec=1,nrec(ipl)
             Write (8,'(I6,3(1X,1pE10.3),I5,1X,I6,1X,E10.3E3,1X,E10.3E3,12(1X,E14.5E3))') &
                  irec,xs(irec,ipl)+xplant(ipl),ys(irec,ipl)+yplant(ipl),zs(irec,ipl),ibrseg(irec,ipl),irecpr(irec,ipl),&
                  Lr(irec,ipl),Khr(irec,ipl),PHs(irec,ipl),PHsri(irec,ipl),PHr(irec+1,ipl),sinkR(irec,ipl),&
                  axialRootFlow(irec,ipl),Qi(irec,ipl),Qd(irec,ipl),Q_bc1(irec,ipl),segrad(irec,ipl),veloRoot(irec,ipl),&
                  segsorb(irec)
          End Do
          Close(8)
       Else
          Do irec=1,nrec(ipl)
             if(lKdrop.and.(TypeKdrop.eq.2))Then ! PHsri is estimated at each time step
                Write (8,'(I6,3(1X,1pE10.3),I5,1X,I6,3X,E10.3E3,1X,E10.3E3,12(1X,E14.5E3))') &
                  irec,xs(irec,ipl)+xplant(ipl),ys(irec,ipl)+yplant(ipl),zs(irec,ipl),ibrseg(irec,ipl),irecpr(irec,ipl),&
                  Lr(irec,ipl),Khr(irec,ipl),PHbulk(irec,ipl),PHsri(irec,ipl),PHr(irec+1,ipl),sinkR(irec,ipl),&
                  axialRootFlow(irec,ipl),Qi(irec,ipl),Qd(irec,ipl),Q_bc1(irec,ipl),segrad(irec,ipl),veloRoot(irec,ipl),&
                  segsorb(irec)
             Else
                Write (8,'(I6,3(1X,1pE10.3),I5,1X,I6,3X,E10.3E3,1X,E10.3E3,12(1X,E14.5E3))') &
                  irec,xs(irec,ipl)+xplant(ipl),ys(irec,ipl)+yplant(ipl),zs(irec,ipl),ibrseg(irec,ipl),irecpr(irec,ipl),&
                  Lr(irec,ipl),Khr(irec,ipl),PHs(irec,ipl),PHsri(irec,ipl),PHr(irec+1,ipl),sinkR(irec,ipl),&
                  axialRootFlow(irec,ipl),Qi(irec,ipl),Qd(irec,ipl),Q_bc1(irec,ipl),segrad(irec,ipl),veloRoot(irec,ipl),&
                  segsorb(irec)
             EndIf
          End Do
          Close(8)
       End If
    Return
  End Subroutine OutDou
  !***************************************************************************
  !> writes balance.out and remove.out 
  Subroutine SubReg(t,iCount)
    Use ParamData, Only: lChem,lretry,last_out
    Use GridData
    Use SolData
    Use WatFun
    Use CumData
    Use Doussanmat, Only : SinkR
    Use Rootdata, Only : lDou,lCou,lFed,lno_RWU,lno_Archi,ldJvL,tlim_dJvL,lPast,tPast,sinkredDPast
    Use tmctrl, Only : tOut
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap) :: i,j,k,l,Mi,Mj,Mk,Ml,iSE,iE,iCount,icut
    Integer(ap) :: xE,yE,zE
    Real(dp),Intent(in) :: t
    Real(dp) :: cbalr=0,cc,cbalt
    Real(dp) ::cnewe,cel,wel,DeltC,time
    Real(dp) ::WatVol,ConVol,ww,wBalT,DeltW,WNewE,VE
    Real(dp), Allocatable,Dimension (:) :: tPast_old
    Real(dp), Allocatable,Dimension (:,:) :: sinkredD,sinkredDPast_old
    Real(dp), Allocatable,Dimension (:) :: SSF_old,RLD_old
    !> \param t current time
    !> \param iCount at the first timestep iCount = 0 

    If (.Not.Allocated(watin)) Allocate(WatIn(nElm),SolIn(nElm))
    If (ldJvL) Allocate (sinkredD(1:nElm,1))
    If (icount.Eq.0) Then
       VElm(1:nElm)=dxGrid*dyGrid*dzGrid
       If (lCou.And.lno_Archi.And.(nexSSF.Ne.nex.Or.neySSF.Ne.ney.Or.nezSSF.Ne.nez)) Then
          If (iCount.Eq.0) Then
             Allocate (SSF_old(1:nexSSF*neySSF*nezSSF))
             SSF_old=SSF
             Deallocate (SSF)
             Allocate (SSF(1:nElm))
             SSF=0._dp
             Do zE=1,nezSSF
                Do yE=1,neySSF
                   Do xE=1,nexSSF
                      SSF(Int(Floor(((Real(xE)-1)*nex)/nexSSF)+1+Floor(((Real(yE)-1)*ney)/neySSF)*nex+Floor(((Real(zE)-1)*nez)/nezSSF)*nex*ney))=SSF(Int(Floor(((Real(xE)-1)*nex)/nexSSF)+1+Floor(((Real(yE)-1)*ney)/neySSF)*nex+Floor(((Real(zE)-1)*nez)/nezSSF)*nex*ney))+SSF_old(xE+(yE-1)*nexSSF+(zE-1)*nexSSF*neySSF)
                   End Do
                End Do
             End Do
          Endif
       Elseif (lFed.And.lno_Archi.And.(nexRLD.Ne.nex.Or.neyRLD.Ne.ney.Or.nezRLD.Ne.nez)) Then
          If (iCount.Eq.0) Then
             Allocate (RLD_old(1:nexRLD*neyRLD*nezRLD))
             RLD_old=RLD
             Deallocate (RLD)
             Allocate (RLD(1:nElm))
             RLD=0._dp
             Do zE=1,nezRLD
                Do yE=1,neyRLD
                   Do xE=1,nexRLD
                      RLD(Int(Floor(((Real(xE)-1)*nex)/nexRLD)+1+Floor(((Real(yE)-1)*ney)/neyRLD)*nex+Floor(((Real(zE)-1)*nez)/nezRLD)*nex*ney))=(RLD(Int(Floor(((Real(xE)-1)*nex)/nexRLD)+1+Floor(((Real(yE)-1)*ney)/neyRLD)*nex+Floor(((Real(zE)-1)*nez)/nezRLD)*nex*ney))*dxGrid*dyGrid*dzGrid+RLD_old(xE+(yE-1)*nexRLD+(zE-1)*nexRLD*neyRLD)*dxRLD*dyRLD*dzRLD)/dxGrid*dyGrid*dzGrid
                   End Do
                End Do
             End Do
          Endif
       Endif
    Endif

    ! write balance.out
    Open (UNIT=10,FILE='out/balance.out',STATUS='OLD',POSITION='APPEND')

    ! initializing some variables
    WatVol=0.0_dp
    DeltW=0.0
    ConVol=0.0
    DeltC=0.0

    Do iE=1,nElm ! loop over elements
       wEl=0.0
       cEl=0.0
       Do iSE=1,5  
          i=elmnod(iL(1,iSE,subN(iE)),iE) ! i,j,k,l are corner nodes of the tetrahedal element, ise counts 5 which is equal to five tedrahedrals which is equal to a cubic.
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          Mi=MatNum(i)
          Mj=MatNum(j)
          Mk=MatNum(k)
          Ml=MatNum(l)
          VE=Abs(Deter(iSE,iE))/6.
          WNewE=VE*(theta(i)+theta(j)+theta(k)+theta(l))/4 !V -> water mass balance (elements)
          WatVol=WatVol+WNewE !sum of V
          wEl=wEl+WNewE
          If (lChem) Then
             CNewE=VE*((theta(i)+ChPar(1,Mi)*ChPar(5,Mi))*Conc(i)+(theta(j)+&
                  ChPar(1,Mj)*ChPar(5,Mj))*Conc(j)+(theta(k)+ChPar(1,Mk)*ChPar(5,Mk))*Conc(k)+&
                  (theta(l)+ChPar(1,Ml)*ChPar(5,Ml))*Conc(l))/4  ! solute mass balance (elements)
             ConVol=ConVol+CNewE !sum of solute mass
             cEl=cEl+CNewE
          Endif
          If (iSE.Eq.5) Then
             If (iCount.Eq.0) Then
                WatIn(iE)=wEl
                If (lChem) SolIn(iE)=cEl
             Else
                DeltW=DeltW+Abs(WatIn(iE)-wEl)
                If (lChem) DeltC=DeltC+Abs(SolIn(iE)-cEl)
             Endif
          Endif
       End Do
       If (ldJvL) sinkredD(iE,1)=sink_cube(iE)
    End Do
    If (ldJvL) Then
       If (iCount.Eq.0) Then
          lPast=3!Number of times saved in tPast, sinkRedD, etc.
          Allocate (tPast(1:lPast))
          tPast=(/0._dp,t-0.0002_dp,t-0.0001_dp/)!Times saved in sinkRedD
          Allocate (sinkredDPast(Size(sinkredD,1),1:lPast))!Dimensions are (number of soil elements ; number of times saved)
          sinkredDPast(:,1)=sinkredD(:,1)!Contains the past values of sink_cube for the latest times, 1st column for the oldest time step.
          sinkredDPast(:,2)=sinkredD(:,1)
          sinkredDPast(:,3)=sinkredD(:,1)
       Else
          Allocate (tPast_old(1:lPast))
          tPast_old=tPast
          Allocate (sinkredDPast_old(Size(sinkredD,1),1:lPast))
          sinkredDPast_old=sinkredDPast!Keep previous information in mind
          icut=0
14        icut=icut+1!tPast and sinkredD are going to be pruned, beginning by the oldest time steps. "icut-1" is the number of pruned time steps.
          If ((tPast(icut)-t).Lt.(-tlim_dJvL)) Then
             lPast=lPast-1
             If (lPast.Gt.0) Then
                Goto 14
             Else
                lPast=lPast+1
             Endif
          Else
             lPast=lPast+1
          Endif
          Deallocate (tPast,sinkredDPast)
          Allocate (tPast(1:lPast))
          Allocate (sinkredDPast(Size(sinkredD,1),1:lPast))
          If (lPast.Gt.1) Then
             tPast=(/tPast_old(icut:lPast+icut-2),t/)!Here we already use the "new" value of lPast.
             sinkredDPast(:,1:lPast-1)=sinkredDPast_old(:,icut:lPast+icut-2)
             sinkredDPast(:,lPast)=sinkredD(:,1)
          Else
             tPast=t
             sinkredDPast(:,1)=sinkredD(:,1)
          Endif
       Endif
    Endif

    ! Mass balance calculation
    If (iCount.Eq.0) Then
       If (.Not.lretry) Then
          wVolI=WatVol
          If (lChem) cVolI=ConVol
          If (lChem) Then 
             Write(10,130)
          Elseif (lno_RWU) Then
             Write(10,135)
          Elseif (lDou) Then
             Write(10,136)
          Else
             Write(10,137)
          Endif
          If (lChem) Then
             Write(10,140) t,WatVol,ConVol
          Else
             Write(10,150) t,WatVol
          Endif
       Else
          Close (10)
          Open (Unit=20,FILE='out/balance.out',STATUS='OLD',ERR=10)
          Read (20,*)
          Read (20,*)
          Read (20,*)
          If (.Not.lChem) Then
             Read (20,*) time,wVolI
          Else
             Read (20,*) time,wVolI,cVolI
          Endif
3         Read (20,*) time,WatVol,wBalT,wBalR,ww,wCumT,WatVol,wCumA,deltW
          If (time.Ne.tOut(last_out)) Goto 3
          Close (20)
          ! write remove.out
          Open (Unit=20,FILE='out/remove.out',STATUS='OLD',ERR=20)
          Read (20,*)
          Read (20,*)
          Read (20,*)
          Read (20,*)
4         Read (20,*) time,CumCh0,CumCh1,CumChR,ChemS(1),ChemS(2),CumRt,CumQ(1),CumQ(2)
          If (time.Ne.tOut(last_out)) Goto 4
          Close (20)
       Endif
    Else
       wBalT=WatVol-wVolI+wCumT
       ww=Max(DeltW,wCumA)
       If (ww.Ge.1.e-25_dp) wBalR=Abs(wBalT)/ww*100  
       If (lChem) Then
          cBalT=ConVol-cVolI+cCumT
          cc=Max(DeltC,cCumA)
          If (cc.Ge.1.e-25_dp) cBalR=Abs(cBalT)/cc*100
          Write(*,*)'cBalR=',cBalR
          Write(10,160) t,WatVol,wBalT,wBalR,ConVol,cBalT,cBalR,Peclet,Courant,wcumT,wcumA,deltW,deltc 
       Elseif (lno_RWU) Then
          Write(10,171) t,WatVol,wBalT,wBalR,peclet,courant,wcumT,WatVol-wvolI,wcumA,deltW
       Elseif (lDou) Then
          Write(10,170) t,WatVol,wBalT,wBalR,ww,wcumT,WatVol-wvolI,wcumA,deltW,RootSk,Sum(SinkR)
       Else
          Write(10,172) t,WatVol,wBalT,wBalR,ww,wcumT,WatVol-wvolI,wcumA,deltW,RootSk
       Endif
    Endif
    Close (10)
    WatVolOld=WatVol
130 Format(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]  CncVol [M]   ',&
         'CncBalT [M]  CncBalR [%]  Peclet       Courant      ',&
         'wCumT [V]    wCumA[V]   DeltaW[V]   DeltaC[V]')
135 Format(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%] ' ,&
         'Peclet    Courant',&
         'wCumT [V]  WatVol-wVolI[V]   wCumA[V]    Deltaw[V] ')
136 Format(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]  ww [V]  ',&
         'wCumT [V]  WatVol-wVolI[V]   wCumA[V]    Deltaw[V]   TotalRootFlow [V/T]',&
         '  TotRadialFlow [V/T]')
137 Format(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]  ww [V]  ',&
         'wCumT [V]  WatVol-wVolI[V]   wCumA[V]    Deltaw[V]   TotalRootFlow [V/T]')
140 Format(f12.4,1X,1pE12.5,27x,1pE12.5)
150 Format(f12.4,1X,1pE12.5)
160 Format(f12.4,1X,14(1pE12.5,1X))
170 Format(f12.4,1X,10(1pE12.5,1X))
171 Format(f12.4,1X,9(1pE12.5,1X))
172 Format(f12.4,1X,9(1pE12.5,1X))
    Open (UNIT=10,FILE='out/remove.out',STATUS='OLD',POSITION='APPEND')
    If (iCount.Eq.0) Then
       If (.Not.lretry) Then
          Write(10,230)
       Endif
    Else
       Write(10,240)  t,CumCh0,CumCh1,CumChR,ChemS(1),ChemS(2),CumRt ,CumQ(1) ,CumQ(2)
    Endif
230 Format(/,'    Time [T]   CumCh0 [M]   CumCh1 [M]   CumChR [M] ChemS(1) [M]',&
         ' ChemS(2) [M]    CumRt [V]  CumQ(1) [V]  CumQ(2) [V]')
240 Format(f12.4,1X,8(1pE12.5,1X))
    Close(10)
    Return
10  Call stop_program('< out/balance.out > not found  -- program terminated.')
20  Call stop_program('< out/remove.out > not found  -- program terminated.')
  End Subroutine SubReg
  !************************************************************************
  !> writes veloci.XX, snkElm.XX, and betElm.XX 
  Subroutine FlxOut(kOut)
    Use typedef
    Use ParamData, Only: lOutPartrace,pi
    Use Soldata
    Use GridData
    Use DoussanMat, Only :betac_cube2
    Use PlntData, Only:TotSur
    Use RootData, Only:lDou,lCou,ldJvL,lFed
    Use Watfun, Only: Fh_from_Th
    Use SoluteMod, Only: Veloc
    Implicit None

    Integer(ap), Intent(in)::kout
    Integer(ap):: i
    Character file*14,  file2*14,  file3*14
    !> \param kOut current number for the FEM output  

    Write (file,'(A11)')'out/veloci.'
    Write (file2,'(A11)')'out/snkElm.'
    Write (file3,'(A11)')'out/betElm.'

    If (kOut.Lt.10) Then
       Write (file(12:12),'(I1)') kOut
       Write (file2(12:12),'(I1)') kOut
       Write (file3(12:12),'(I1)') kOut
    Elseif (kout.Lt.100) Then
       Write (file(12:13),'(I2)') kOut
       Write (file2(12:13),'(I2)') kOut
       Write (file3(12:13),'(I2)') kOut
    Else
       Write (file(12:14),'(I3)') kOut
       Write (file2(12:14),'(I3)') kOut
       Write (file3(12:14),'(I3)') kOut
    Endif

    Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
    Write (8,'(/3X,''Node#'',10X,''x'',11X,''y'',11X,''z'',11X,''Vx'',10X,''Vy'',10X,''Vz'',10X,''Width'')')
    Call Veloc

    ! write nodal values:
    Do i=1,nPt
       Write (8,'(I5,8X,3(1X,1pE11.4),4(1X,1pE11.4))')i,xGrid(i),yGrid(i),zGrid(i),Vx(i),Vy(i),Vz(i),width(i)
    End Do
    Close(8)

    If (lOutPartrace) Then
       Open (UNIT=50,FILE=file3,STATUS='UNKNOWN')
       Write (50,'(/''Element Nr.'',9X,''BetaElm'')' )
       Do i=1,nElm
          Write (50,*) i,  betac_cube2 (i)*TotSur
       End Do
       Close(50)
    End If
    If (lOutPartrace.Or.lDou) Then
       Open (UNIT=90,FILE=file2,STATUS='UNKNOWN')
       Write (90,'(/''Element Nr.'',9X,''SinkElm'')' )
       Do i=1,nElm
          Write (90,*) i,Sink_cube(i)
       End Do
       Close(90)
    Elseif ((lCou.And..Not.ldJvL).Or.lFed) Then
       Open (UNIT=90,FILE=file2,STATUS='UNKNOWN')
       Write (90,'(/''Element Nr.'',9X,''SinkElm'',9X,''HElm'')' )
       Do i=1,nElm
          Write (90,*) i,Sink_cube(i),HElm(i)
       End Do
       Close(90)
    Elseif (lCou.And.ldJvL) Then
       Open (UNIT=90,FILE=file2,STATUS='UNKNOWN')
       Write (90,'(/''Element Nr.'',9X,''SinkElm'',9X,''Hint'',9X,''Hbulk'')' )
       Do i=1,nElm
          Write (90,*) i,Sink_cube(i),HElm(i),Fh_from_Th(Sum(theta(elmnod(1:8,i)))/8.0_dp,par(:,MatNum(elmnod(1,i))))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_dp
       End Do
       Close(90)
    End If
    Return
  End Subroutine FlxOut
  !****************************************************************************
  !> writes .vtk output for the soil grid; veloci.XX 
  Subroutine OutVTK(kOut)
    Use typedef
    Use Soldata
    Use GridData
    Use ParamData, Only: lchem,lSalinity,lPartUp
    Use RootData, Only: lSomma_growth,lRootBox_growth,lCou
    Use StrData, Only: s
    Use SoluteRootMat, Only: SoilSoluteUptake,SoilUptake_adv,SoilUptake_diff
    Use DoussanMat, Only: PHs_osmotic
    Use RhizoData
    Use SoluteMod, Only: Veloc
    Implicit None

    Integer(ap), Intent(in)::kout
    Integer(ap):: i, k, j, ix, iy, iz
    Character file*22
    !> \param kOut current number for the FEM output

    Write (file,'(A15)')'out/vtk/veloci.'

    If (kOut.Lt.10) Then
       Write (file(16:16),'(I1)') kOut
       Write (file(17:20),'(A4)') '.vtk'
    Elseif (kOut.Lt.100) Then
       Write (file(16:17),'(I2)') kOut
       Write (file(18:21),'(A4)') '.vtk'
    Else
       Write (file(16:18),'(I3)') kOut
       Write (file(19:22),'(A4)') '.vtk'
    Endif

    Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
    Write (8,'(A26)')'# vtk DataFile Version 3.0' 
    Write (8,'(A12)')'model R-SWMS'
    Write (8,'(A5)')'ASCII'
    Write (8,'(A24)')'DATASET RECTILINEAR_GRID'

    If (continu) Then

       Write (8,'(A11,1X,3I6)')'DIMENSIONS ', nx+1, ny+1, nz

       Write (8,'(A14,1X,I6,1X,A5)')'X_COORDINATES ',nx+1, 'float' 
       Do i=1,nx
          Write (8,'(1X,1pE11.4)',advance="no") xgrid(i)
       End Do
       Write (8,'(1X,1pE11.4)',advance="no") xgrid(Size(xgrid))+dxgrid

       Write (8,'(/A14,1X,I6,1X,A5)')'Y_COORDINATES ',ny+1, 'float'
       Do i=1,ny*nx,nx
          Write (8,'(1pE11.4)',advance="no") ygrid(i)
       End Do
       Write (8,'(1X,1pE11.4)',advance="no") ygrid(Size(ygrid))+dygrid

       Write (8,'(/A14,1X,I6,1X,A5)')'Z_COORDINATES ',nz, 'float' 
       Do i=1,Size(zgrid),nx*ny
          Write (8,'(1pE11.4)',advance="no") zgrid(i)
       End Do

       Write (8,'(/A10,1X,I7)')'POINT_DATA', (nx+1)*(ny+1)*nz
       Write (8,'(A24)')'SCALARS velocity float 3'
       Write (8,'(A20)')'LOOKUP_TABLE default'

       Call Veloc

       k=0
       j=0
       Do  iz=1,nz
          Do  iy=1,ny+1
             Do  ix=1,nx+1     
                k = ix + (iy-1)*nx + (iz-1)*nx*ny
                j = j +1
                If (ix.Eq.nx+1) Then
                   k = k - nx
                End If
                If  (iy.Eq.ny+1) Then
                   k = k - (nx*ny)
                End If
                Write (8,'(3(1X,1pE11.4))')Vx(k),Vy(k),Vz(k)
             End Do
          End Do
       End Do

       Write (8,'(A27)')'SCALARS pressurehead float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       k=0
       j=0
       Do  iz=1,nz
          Do  iy=1,ny+1
             Do  ix=1,nx+1     
                k = ix + (iy-1)*nx + (iz-1)*nx*ny
                j = j +1
                If (ix.Eq.nx+1) Then
                   k = k - nx
                End If
                If  (iy.Eq.ny+1) Then
                   k = k - (nx*ny)
                End If
                Write (8,'(1X,1pE11.4)',advance="no")hnew(k)
             End Do
          End Do
       End Do
       Write (8,*)''
       Write (8,'(A17)')'SCALARS wc float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       k=0
       j=0
       Do  iz=1,nz
          Do  iy=1,ny+1
             Do  ix=1,nx+1     
                k = ix + (iy-1)*nx + (iz-1)*nx*ny
                j = j +1
                If (ix.Eq.nx+1) Then
                   k = k - nx
                End If
                If  (iy.Eq.ny+1) Then
                   k = k - (nx*ny)
                End If
                Write (8,'(1X,1pE11.4)',advance="no")theta(k)
             End Do
          End Do
       End Do
       Write (8,*)''
       Write (8,'(A17)')'SCALARS mat integer '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       k=0
       j=0
       Do  iz=1,nz
          Do  iy=1,ny+1
             Do  ix=1,nx+1     
                k = ix + (iy-1)*nx + (iz-1)*nx*ny
                j = j +1
                If (ix.Eq.nx+1) Then
                   k = k - nx
                End If
                If  (iy.Eq.ny+1) Then
                   k = k - (nx*ny)
                End If
                Write (8,'(1X,I3)',advance="no") MatNum(k)
             End Do
          End Do
       End Do
       Write (8,*)''
       Write (8,'(A17)')'SCALARS Kode integer '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       k=0
       j=0
       Do  iz=1,nz
          Do  iy=1,ny+1
             Do  ix=1,nx+1     
                k = ix + (iy-1)*nx + (iz-1)*nx*ny
                j = j +1
                If (ix.Eq.nx+1) Then
                   k = k - nx
                End If
                If  (iy.Eq.ny+1) Then
                   k = k - (nx*ny)
                End If
                Write (8,'(1X,I3)',advance="no") Kode(k)
             End Do
          End Do
       End Do

       If(lchem) Then
          Write (8,*)''
          Write (8,'(A19)')'SCALARS conc float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          k=0
          j=0
          Do  iz=1,nz
             Do  iy=1,ny+1
                Do  ix=1,nx+1     
                   k = ix + (iy-1)*nx + (iz-1)*nx*ny
                   j = j +1
                   If (ix.Eq.nx+1) Then
                      k = k - nx
                   End If
                   If  (iy.Eq.ny+1) Then
                      k = k - (nx*ny)
                   End If
                   Write (8,'(1X,1pE11.4)',advance="no")conc(k)
                End Do
             End Do
          End Do
       End If
       Write (8,*)''
       If (lRhizo) Then
          Write (8,'(A20)')'SCALARS thtTot float'
          Write (8,'(A20)')'LOOKUP_TABLE default'
          k=0
          j=0
          Do iz=1,nz
             Do iy=1,ny+1
                Do ix=1,nx+1
                   k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                   j = j + 1
                   If (ix .Eq. nx+1) Then
                      k = k -nx
                   Endif
                   If (iy .Eq. ny + 1) Then
                      k = k - (nx*ny)
                   Endif
                   Write (8, '(1X,1pE11.4)', advance="no") thetaTot(k)
                Enddo
             Enddo
          Enddo
       Endif
       If (lRhizo) Then
          Write (8,'(A19)')'SCALARS Rnorm float'
          Write (8,'(A20)')'LOOKUP_TABLE default'
          k=0
          j=0
          Do iz=1,nz
             Do iy=1,ny+1
                Do ix=1,nx+1
                   k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                   j = j + 1
                   If (ix .Eq. nx+1) Then
                      k = k -nx
                   Endif
                   If (iy .Eq. ny + 1) Then
                      k = k - (nx*ny)
                   Endif
                   Write (8, '(1X,1pE11.4)', advance="no") Rnorm(k)
                Enddo
             Enddo
          Enddo
       Endif
       If (lRhizo) Then
          Write (8,'(A18)')'SCALARS cTot float'
          Write (8,'(A20)')'LOOKUP_TABLE default'
          k=0
          j=0
          Do iz=1,nz
             Do iy=1,ny+1
                Do ix=1,nx+1
                   k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                   j = j + 1
                   If (ix .Eq. nx+1) Then
                      k = k -nx
                   Endif
                   If (iy .Eq. ny + 1) Then
                      k = k - (nx*ny)
                   Endif
                   Write (8, '(1X,1pE11.4)', advance="no") cTot_r(k)
                Enddo
             Enddo
          Enddo
       Endif

       If (lRhizo) Then
          Write (8,'(A18)')'SCALARS cond float'
          Write (8,'(A20)')'LOOKUP_TABLE default'
          k=0
          j=0
          Do iz=1,nz
             Do iy=1,ny+1
                Do ix=1,nx+1
                   k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                   j = j + 1
                   If (ix .Eq. nx+1) Then
                      k = k -nx
                   Endif
                   If (iy .Eq. ny + 1) Then
                      k = k - (nx*ny)
                   Endif
                   Write (8, '(1X,1pE11.4)', advance="no") con(k)
                Enddo
             Enddo
          Enddo
       Endif


       If((lSomma_growth.OR.(lRootBox_growth))) Then
          Write (8,'(A27)')'SCALARS SoilStrength float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          k=0
          j=0
          Do  iz=1,nz
             Do  iy=1,ny+1
                Do  ix=1,nx+1     
                   k = ix + (iy-1)*nx + (iz-1)*nx*ny
                   j = j +1
                   If (ix.Eq.nx+1) Then
                      k = k - nx
                   End If
                   If  (iy.Eq.ny+1) Then
                      k = k - (nx*ny)
                   End If
                   Write (8,'(1X,1pE11.4)',advance="no") s(k)
                End Do
             End Do
          End Do
       End If

       If(lchem) Then
          Write (8,*)''
          Write (8,'(A20)')'SCALARS csink float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          k=0
          j=0
          Do  iz=1,nz
             Do  iy=1,ny+1
                Do  ix=1,nx+1     
                   k = ix + (iy-1)*nx + (iz-1)*nx*ny
                   j = j +1
                   If (ix.Eq.nx+1) Then
                      k = k - nx
                   End If
                   If  (iy.Eq.ny+1) Then
                      k = k - (nx*ny)
                   End If
                   Write (8,'(1X,1pE11.4)',advance="no")csink(k)
                End Do
             End Do
          End Do
       End If

    Else

       Write (8,'(A11,1X,3I6)')'DIMENSIONS ', nx, ny, nz

       Write (8,'(A14,1X,I6,1X,A5)')'X_COORDINATES ',nx, 'float' 
       Do i=1,nx
          Write (8,'(1X,1pE11.4)',advance="no") xgrid(i)
       End Do

       Write (8,'(/A14,1X,I6,1X,A5)')'Y_COORDINATES ',ny, 'float'
       Do i=1,ny*nx,nx
          Write (8,'(1pE11.4)',advance="no") ygrid(i)
       End Do

       Write (8,'(/A14,1X,I6,1X,A5)')'Z_COORDINATES ',nz, 'float' 
       Do i=1,Size(zgrid),nx*ny
          Write (8,'(1pE11.4)',advance="no") zgrid(i)
       End Do

       Write (8,'(/A15,1X,I7)')'POINT_DATA', nPt
       Write (8,'(A24)')'SCALARS velocity float 3'
       Write (8,'(A20)')'LOOKUP_TABLE default'

       Call Veloc
       Do i=1,nPt
          Write (8,'(3(1X,1pE11.4))')Vx(i),Vy(i),Vz(i)
       End Do

       Write (8,'(A27)')'SCALARS pressurehead float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE11.4)',advance="no") hnew(1:nPt)


       Write (8,'(A17)')'SCALARS wc float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE11.4)',advance="no") theta(1:nPt)

      If(lChem) then
        Write (8,'(A19)')'SCALARS conc  float'
        Write (8,'(A20)')'LOOKUP_TABLE default'
        Write (8,'(1X,1pE11.4)',advance="no") conc(1:nPt)

        Write (8,'(A20)')'SCALARS csink float '
        Write (8,'(A20)')'LOOKUP_TABLE default'
        Write (8,'(1X,1pE11.4)',advance="no") csink(1:nPt)
     End If
       
       Write (8,'(A19)')'SCALARS mat integer'
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,I3)',advance="no") MatNum(1:nPt)
       
       WRITE (8,'(A19)')'SCALARS width float'
       WRITE (8,'(A20)')'LOOKUP_TABLE default'
       WRITE (8,'(1X,1pE11.4)',advance="no") width(1:nPt) 

       If (lRhizo) Then
          Write (8, '(A20)')'SCALARS thtTot float'
          Write (8, '(A20)')'LOOKUP_TABLE default'
          Write (8, '(1X,1pE11.4)', advance="no") thetaTot(1:nPt)

          Write (8, '(A19)')'SCALARS Rnorm float'
          Write (8, '(A20)')'LOOKUP_TABLE default'
          Write (8, '(1X,1pE11.4)', advance="no") Rnorm(1:nPt)

          Write (8, '(A18)')'SCALARS cTot float'
          Write (8, '(A20)')'LOOKUP_TABLE default'
          Write (8, '(1X,1pE11.4)', advance="no") cTot_r(1:nPt)

          Write (8, '(A18)')'SCALARS cond float'
          Write (8, '(A20)')'LOOKUP_TABLE default'
          Write (8, '(1X,1pE11.4)', advance='no') con(1:nPt)
       Endif

       Write (8,'(A20)')'SCALARS Kode integer '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,I3)',advance="no") Kode(1:nPt)

       If((lSomma_growth.OR.(lRootBox_growth))) Then
          Write (8,'(A27)')'SCALARS SoilStrength float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") s(1:nPt)
       End If

    Endif

    Write (8,*)''
    Write (8,'(/A9,1X,I7)')'CELL_DATA', nElm
    Write (8,'(A24)')'SCALARS sinkElm float'
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE11.4)',advance="no") sink_cube(1:nElm)

    If (lCou) Then
       Write (8,'(A24)')'SCALARS SSF float'
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE11.4)',advance="no") SSF(1:nElm)
    Endif

    Write (8,'(A24)')'SCALARS csinkElm float'
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE11.4)',advance="no") csink_cube(1:nElm)

    Write (8,'(A24)')'SCALARS cmassElm float'
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE11.4)',advance="no") cmass_cube(1:nElm)  


    If(lPartUp) Then
       Write (8,'(A37)')'SCALARS SoilSoluteConcentration float'
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE13.4)',advance="no") SoilSoluteConcentration(1:nElm)
       Write (8,*)''

       Write (8,'(A30)')'SCALARS SoilSoluteUptake float'
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE13.4)',advance="no") SoilSoluteUptake(1:nElm)
       Write (8,*)''

       Write (8,'(A30)')'SCALARS SoilUptake_adv float'
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE13.4)',advance="no") SoilUptake_adv(1:nElm)
       Write (8,*)''

       Write (8,'(A30)')'SCALARS SoilUptake_diff float'
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE13.4)',advance="no") SoilUptake_diff(1:nElm)
       Write (8,*)''

    End If

    If (lsalinity) Then
       Write (8,*)''
       Write (8,'(A27)')'SCALARS osmoticHead float'
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE11.4)',advance="no") PHs_osmotic(1:nElm)
       Write (8,*)''
    End If

    Close (8)

    Return
  End Subroutine OutVTK
  !************************************************************************
  !> generates .vtk output for the root outRoot.XX.XXX 
  Subroutine OutDouVTK(kOut,t,ipl)
    Use typedef
    Use ParamData, Only: pi,lPartUp
    Use RootData
    Use DoussanMat, Only: Lr,Khr,PHs,PHr,sinkR,axialRootFlow,Qi,Qd,Q_bc1,nsub,transroot,&
         transtip,nplant
    Use GridData, Only: dxgrid,dygrid,nx,ny,continu
    Use SoluteRootMat, Only: segsorb,l_linSorb,l_freundSorb,seg_upt
    Implicit None

    Integer(ap):: j, sizen, kk
    Integer(ap) :: irec,kout,ipl,trans_i(2), trans_j(2)
    Integer(ap):: igrow,irecn
    Real(dp) ::  xr(nrec(nplant)+ngrow(nplant)), yr(nrec(nplant)+ngrow(nplant)),t
    Character file*29 
    !> \param t current simulation time
    !> \param kOut current number for the FEM output

       Write (file,'(A18)')'out/vtk/outRootXX.'


       If (ipl.Lt.10) Then
          Write (file(16:16),'(I1)') 0
          Write (file(17:17),'(I1)') ipl
       Else
          Write (file(16:17),'(I2)') ipl
       End If
       If (kOut.Lt.10) Then
          Write (file(19:19),'(I1)') kOut
          Write (file(20:24),'(A4)') '.vtk'
       Elseif (kOut.Lt.100) Then
          Write (file(19:20),'(I2)') kOut
          Write (file(21:25),'(A4)') '.vtk'
       Else
          Write (file(19:21),'(I3)') kOut
          Write (file(22:26),'(A4)') '.vtk'
       Endif

       Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
       Write (8,'(A26)')'# vtk DataFile Version 3.0' 
       Write (8,'(A12)')'model R-SWMS'
       Write (8,'(A5)')'ASCII'
       Write (8,'(A16)')'DATASET POLYDATA'
       Write (8,'(A7,1X,I6,1X,A5)')'POINTS ', nrec(ipl)+ngrow(ipl), 'float'


       If (continu) Then

          ! write out points root segments
          Do  irec=1,nrec(ipl)
             xr(irec) = xs(irec,ipl)+xplant(ipl)+1.*transroot(irec,1,nsub(irec,ipl),ipl)*dxgrid*nx
             yr(irec) = ys(irec,ipl)+yplant(ipl)+1.*transroot(irec,2,nsub(irec,ipl),ipl)*dygrid*ny
             Write (8,'(1X,3(1X,E12.5))',advance="no") xr(irec),yr(irec),zs(irec,ipl)
             Write (8,*)''
          End Do

          ! write out points root tips
          Do  igrow=1,ngrow(ipl)
             xr(nrec(ipl)+igrow) = xg(igrow,ipl)+xplant(ipl)+1.*transtip(igrow,1,nsub(1,ipl),ipl)*dxgrid*nx
             yr(nrec(ipl)+igrow) = yg(igrow,ipl)+yplant(ipl)+1.*transtip(igrow,2,nsub(1,ipl),ipl)*dygrid*ny
             Write (8,'(1X,3(1X,E12.5))',advance="no") xr(nrec(ipl)+igrow),yr(nrec(ipl)+igrow),zg(igrow,ipl)
             Write (8,*)''
          End Do

          kk=0
          Do j=1,nrec(ipl) !segments
             trans_i = transroot(irecpr(j,ipl),:,nsub(irecpr(j,ipl),ipl),ipl) 
             trans_j = transroot(j,:,nsub(j,ipl),ipl)
             If (All(trans_i.Eq.trans_j)) Then
                If(irecpr(j,ipl).Ne.0) Then
                   kk = kk +1
                End If

             Endif
          End Do

          Do j=1,nbr(ipl)! tips
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   trans_i = transtip(igrow,:,nsub(nrec(ipl)+igrow,ipl),ipl)
                   trans_j = transroot(irecn,:,nsub(irecn,ipl),ipl)
                   If (All(trans_i.Eq.trans_j)) Then
                      kk = kk +1
                   End If
                Endif
             End Do
          End Do


          Write (8,'(A6,1X,I8,1X,I9)')'LINES ',kk, kk*3

          Do j=1,nrec(ipl) !tips
             trans_i = transroot(irecpr(j,ipl),:,nsub(irecpr(j,ipl),ipl),ipl)
             trans_j = transroot(j,:,nsub(j,ipl),ipl)
             If (All(trans_i.Eq.trans_j)) Then
                If(irecpr(j,ipl).Ne.0) Then
                   Write (8,'(3(1X,I8))',advance="no") 2, irecpr(j,ipl)-1, j-1
                   Write (8,*)''
                End If

             Endif
          Enddo


          Do j=1,nbr(ipl) !segments
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   trans_i =  transtip(igrow,:,nsub(nrec(ipl)+igrow,ipl),ipl)
                   trans_j = transroot(irecn,:,nsub(irecn,ipl),ipl)
                   If (All(trans_i.Eq.trans_j)) Then
                      Write (8,'(3(1X,I8))',advance="no") 2, irecn-1, nrec(ipl)+igrow-1
                      Write (8,*)''
                   End If
                Endif
             End Do
          End Do



       Else

          ! write points segments
          Do  irec=1,nrec(ipl)
             Write (8,'(1X,3(1X,1pE10.3))',advance="no") xs(irec,ipl)+xplant(ipl),ys(irec,ipl)+yplant(ipl),zs(irec,ipl)
             Write (8,*)''
          End Do

          ! write points tips
          Do  igrow=1,ngrow(ipl)
             Write (8,'(1X,3(1X,E12.5))',advance="no") xg(igrow,ipl)+xplant(ipl),yg(igrow,ipl)+yplant(ipl),zg(igrow,ipl)
             Write (8,*)''
          End Do

          sizen = 0
          Do j=1,nrec(ipl) !segments
             If(irecpr(j,ipl).Ne.0) Then
                sizen = sizen +1
             Endif
          Enddo

          Do j=1,nbr(ipl) !tips
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   sizen = sizen +1
                End If
             End Do
          End Do



          Write (8,'(A6,1X,I8,1X,I9)')'LINES ', sizen, sizen*3

          Do j=1,nrec(ipl) !segments
             If(irecpr(j,ipl).Ne.0) Then
                Write (8,'(3(1X,I8))',advance="no") 2, irecpr(j,ipl)-1, j-1
                Write (8,*)''
             End If
          Enddo

          Do j=1,nbr(ipl) !tips
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   Write (8,'(3(1X,I8))',advance="no") 2, irecn-1, nrec(ipl)+igrow-1
                   Write (8,*)''
                Endif
             End Do
          End Do

       Endif !end if continu

       Write (8,'(A11,1X,I6)')'POINT_DATA ', nrec(ipl)+ngrow(ipl)
       If(.Not.lno_RWU) Then
          Write (8,'(A17)')'SCALARS xs float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") xs(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") xs(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A17)')'SCALARS ys float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") ys(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") ys(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A17)')'SCALARS zs float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") zs(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") zs(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A17)')'SCALARS Lr float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") Lr(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Lr(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A18)')'SCALARS Khr float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") Khr(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Khr(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A21)')'SCALARS PHinter float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") PHs(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") PHs(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A21)')'SCALARS PHxylem float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") PHr(2:nrec(ipl)+1,ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Phr(irecsg(1:ngrow(ipl),ipl)+1,ipl)
          Write (8,*)''

          Write (8,'(A28)')'SCALARS radialRootFlow float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") sinkR(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") sinkR(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A28)')'SCALARS axialRootFlow float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") axialRootFlow(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") axialRootFlow(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A17)')'SCALARS Qi float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") Qi(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Qi(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A17)')'SCALARS Qd float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") Qd(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Qd(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A17)')'SCALARS Qb float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") Q_bc1(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Q_bc1(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

        If(lPartUp) Then
          Write (8,'(A25)')'SCALARS segconc float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") segconc(1:nrec(ipl))
          Write (8,'(1X,1pE11.4)',advance="no") segconc(irecsg(1:ngrow(ipl),ipl))
          Write (8,*)''

          Write (8,'(A25)')'SCALARS seg_upt float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") seg_upt(1:nrec(ipl))
          DO igrow=1,ngrow(ipl)
                          Write (8,'(1X,1pE11.4)',advance="no") 0E+00
          END DO
          Write (8,*)''

          If(l_linSorb .Or. l_freundSorb) Then
             Write (8,'(A21)')'SCALARS segsorb float '
             Write (8,'(A20)')'LOOKUP_TABLE default'
             Write (8,'(1X,1pE11.4)',advance="no") segsorb(1:nrec(ipl))
             Do igrow=1,ngrow(ipl)
                Write (8,'(1X,1pE11.4)',advance="no") 0E+00
             End Do
             Write (8,*)'' 
          End If
        End if
       End If

       Write (8,'(A19)')'SCALARS age float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE11.4)',advance="no") t-timorg(1:nrec(ipl),ipl)
       Do irec=1,ngrow(ipl)
          Write (8,'(1X,1pE11.4)',advance="no") 0E+00
       End Do
       Write (8,*)''

       Write (8,'(A28)')'SCALARS segmentRadius float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Do  irec=1,nrec(ipl)
          Write (8,'(1X,1pE14.8)',advance="no") segrad(irec,ipl) 
       End Do
       Do  irec=1,ngrow(ipl)
          Write (8,'(1X,1pE14.8)',advance="no") 0E+00 
       End Do
       Write (8,*)''

       Close(8)
       Return

  End Subroutine OutDouVTK
  !************************************************************************
  Subroutine OutCouVTK(t,ipl)
    Use typedef
    Use ParamData, Only: pi
    Use RootData
    Use DoussanMat, Only: Lr,Khr,nsub,transroot,transtip
    Use GridData, Only: dxgrid,dygrid,nx,ny,continu
    Implicit None

    Integer(ap) :: j, sizen,  kk
    Integer(ap) :: irec, ipl, trans_i(2), trans_j(2)
    Integer(ap) :: igrow,irecn
    Real(dp) ::  xr(nrec(nplant)+ngrow(nplant)), yr(nrec(nplant)+ngrow(nplant)),t
    Character file*24 
    !> \param t current simulation time
    !> \param kOut current number for the FEM output

       Write (file,'(A18)')'out/vtk/outRootXX.'


       If (ipl.Lt.10) Then
          Write (file(16:16),'(I1)') 0
          Write (file(17:17),'(I1)') ipl
       Else
          Write (file(16:17),'(I2)') ipl
       End If
       Write (file(19:19),'(I1)') 0
       Write (file(20:24),'(A4)') '.vtk'


       Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
       Write (8,'(A26)')'# vtk DataFile Version 3.0' 
       Write (8,'(A12)')'model R-SWMS'
       Write (8,'(A5)')'ASCII'
       Write (8,'(A16)')'DATASET POLYDATA'
       Write (8,'(A7,1X,I6,1X,A5)')'POINTS ', nrec(ipl)+ngrow(ipl), 'float'


       If (continu) Then

          ! write out points root segments
          Do  irec=1,nrec(ipl)
             xr(irec) = xs(irec,ipl)+xplant(ipl)+1.*transroot(irec,1,nsub(irec,ipl),ipl)*dxgrid*nx
             yr(irec) = ys(irec,ipl)+yplant(ipl)+1.*transroot(irec,2,nsub(irec,ipl),ipl)*dygrid*ny
             Write (8,'(1X,3(1X,E12.5))',advance="no") xr(irec),yr(irec),zs(irec,ipl)
             Write (8,*)''
          End Do

          ! write out points root tips
          Do  igrow=1,ngrow(ipl)
             xr(nrec(ipl)+igrow) = xg(igrow,ipl)+xplant(ipl)+1.*transtip(igrow,1,nsub(1,ipl),ipl)*dxgrid*nx
             yr(nrec(ipl)+igrow) = yg(igrow,ipl)+yplant(ipl)+1.*transtip(igrow,2,nsub(1,ipl),ipl)*dygrid*ny
             Write (8,'(1X,3(1X,E12.5))',advance="no") xr(nrec(ipl)+igrow),yr(nrec(ipl)+igrow),zg(igrow,ipl)
             Write (8,*)''
          End Do

          kk=0
          Do j=1,nrec(ipl) !segments
             trans_i = transroot(irecpr(j,ipl),:,nsub(irecpr(j,ipl),ipl),ipl)
             trans_j = transroot(j,:,nsub(j,ipl),ipl)
             If (All(trans_i.Eq.trans_j)) Then
                If(irecpr(j,ipl).Ne.0) Then
                   kk = kk +1
                End If

             Endif
          End Do

          Do j=1,nbr(ipl)! tips
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   trans_i = transtip(igrow,:,nsub(nrec(ipl)+igrow,ipl),ipl)
                   trans_j = transroot(irecn,:,nsub(irecn,ipl),ipl)
                   If (All(trans_i.Eq.trans_j)) Then
                      kk = kk +1
                   End If
                Endif
             End Do
          End Do


          Write (8,'(A6,1X,I8,1X,I9)')'LINES ',kk, kk*3

          Do j=1,nrec(ipl) !tips
             trans_i = transroot(irecpr(j,ipl),:,nsub(irecpr(j,ipl),ipl),ipl)
             trans_j = transroot(j,:,nsub(j,ipl),ipl)
             If (All(trans_i.Eq.trans_j)) Then
                If(irecpr(j,ipl).Ne.0) Then
                   Write (8,'(3(1X,I8))',advance="no") 2, irecpr(j,ipl)-1, j-1
                   Write (8,*)''
                End If

             Endif
          Enddo


          Do j=1,nbr(ipl) !segments
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   trans_i = transroot(nrec(ipl)+igrow,:,nsub(nrec(ipl)+igrow,ipl),ipl)
                   trans_j = transroot(irecn,:,nsub(irecn,ipl),ipl)
                   If (All(trans_i.Eq.trans_j)) Then
                      Write (8,'(3(1X,I8))',advance="no") 2, irecn-1, nrec(ipl)+igrow-1
                      Write (8,*)''
                   End If
                Endif
             End Do
          End Do



       Else

          ! write points segments
          Do  irec=1,nrec(ipl)
             Write (8,'(1X,3(1X,1pE10.3))',advance="no") xs(irec,ipl)+xplant(ipl),ys(irec,ipl)+yplant(ipl),zs(irec,ipl)
             Write (8,*)''
          End Do

          ! write points tips
          Do  igrow=1,ngrow(ipl)
             Write (8,'(1X,3(1X,E12.5))',advance="no") xg(igrow,ipl)+xplant(ipl),yg(igrow,ipl)+yplant(ipl),zg(igrow,ipl)
             Write (8,*)''
          End Do

          sizen = 0
          Do j=1,nrec(ipl) !segments
             If(irecpr(j,ipl).Ne.0) Then
                sizen = sizen +1
             Endif
          Enddo

          Do j=1,nbr(ipl) !tips
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   sizen = sizen +1
                End If
             End Do
          End Do



          Write (8,'(A6,1X,I8,1X,I9)')'LINES ', sizen, sizen*3

          Do j=1,nrec(ipl) !segments
             If(irecpr(j,ipl).Ne.0) Then
                Write (8,'(3(1X,I8))',advance="no") 2, irecpr(j,ipl)-1, j-1
                Write (8,*)''
             End If
          Enddo

          Do j=1,nbr(ipl) !tips
             irecn=nrec(ipl)
             Do While (ibrseg(irecn,ipl).Ne.j)
                irecn=irecn-1
             End Do
             Do igrow=1,ngrow(ipl)
                If (ibrgrw(igrow,ipl)==j) Then
                   Write (8,'(3(1X,I8))',advance="no") 2, irecn-1, nrec(ipl)+igrow-1
                   Write (8,*)''
                Endif
             End Do
          End Do

       Endif !end if continu

       Write (8,'(A11,1X,I6)')'POINT_DATA ', nrec(ipl)+ngrow(ipl)
       If(.Not.lno_RWU) Then
          Write (8,'(A17)')'SCALARS Lr float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") Lr(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Lr(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''

          Write (8,'(A18)')'SCALARS Khr float '
          Write (8,'(A20)')'LOOKUP_TABLE default'
          Write (8,'(1X,1pE11.4)',advance="no") Khr(1:nrec(ipl),ipl)
          Write (8,'(1X,1pE11.4)',advance="no") Khr(irecsg(1:ngrow(ipl),ipl),ipl)
          Write (8,*)''
       End If

       Write (8,'(A19)')'SCALARS age float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Write (8,'(1X,1pE11.4)',advance="no") t-timorg(1:nrec(ipl),ipl)
       Do irec=1,ngrow(ipl)
          Write (8,'(1X,1pE11.4)',advance="no") 0E+00
       End Do
       Write (8,*)''

       Write (8,'(A28)')'SCALARS segmentRadius float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       Do  irec=1,nrec(ipl)
          Write (8,'(1X,1pE14.8)',advance="no") segrad(irec,ipl) 
       End Do
       Do  irec=1,ngrow(ipl)
          Write (8,'(1X,1pE14.8)',advance="no") 0E+00 
       End Do
       Write (8,*)''

       Close(8)

       Return

  End Subroutine OutCouVTK
  !************************************************************************
  !> generates output in case partrace was running
  Subroutine PartraceOut(t)
    Use typedef
    Use Soldata
    Use GridData
    Use SoluteMod, Only: Veloc
    Implicit None

    Integer(ap):: i,j,k,l,ix,iy,iz
    Real(dp),Intent(in)::t
    Logical,Save::lfirst=.True.
    Character file*31,file2*19,file3*36
    !> \param t current simulation time

    If (lfirst) Then
       Open (UNIT=8,FILE='out/Partrace/out.velocity.00001',STATUS='UNKNOWN')
       Open (UNIT=10,FILE='out/Partrace/out.water_content.00001',STATUS='UNKNOWN')

       Write (10,'(a14)') '#water content'
       Write (10,'(a15)') '#CPU: 1 from: 1'

       Write (8,'(a9)') '#velocity'
       Write (8,'(a15)') '#CPU: 1 from: 1'
       If (continu) Then
          Write (10,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', (nx+1)*(ny+1)*nz, '  local: ',(nx+1)*(ny+1)*nz
          Write (10,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx+1,'  from: ',nx+1
          Write (10,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny+1,'  from: ',ny+1
          Write (10,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz

          Write (8,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', (nx+1)*(ny+1)*nz, '  local: ',(nx+1)*(ny+1)*nz
          Write (8,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx+1,'  from: ',nx+1
          Write (8,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny+1,'  from: ',ny+1
          Write (8,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz
       Else
          Write (10,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', nx*ny*nz, '  local: ',nx*ny*nz
          Write (10,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx,'  from: ',nx
          Write (10,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny,'  from: ',ny
          Write (10,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz

          Write (8,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', nx*ny*nz, '  local: ',nx*ny*nz
          Write (8,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx,'  from: ',nx
          Write (8,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny,'  from: ',ny
          Write (8,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz
       Endif
       Close(8)
       Close(10)
       lfirst = .False.
    Endif

    Write (file,'(A31)')'out/Partrace/out.velocity.00001'
    Write (file2,'(A19)')'out/Partrace/geoPar'
    Write (file3,'(A36)')'out/Partrace/out.water_content.00001'

    Open (UNIT=8,FILE=file,STATUS='OLD',POSITION='APPEND')
    Open (UNIT=9,FILE=file2,STATUS='UNKNOWN')
    Open (UNIT=10,FILE=file3,STATUS='OLD',POSITION='APPEND')      

    Write (10,'(A19,F0.3)') '#simulation time:  ',t-100
    Write (8,'(A19,F0.3)') '#simulation time:  ',t-100

    Call Veloc
    !> write nodal soil values:
    If (continu) Then
       k=0
       j=0
       Do  iz=nz,1,-1
          Do  iy=1,ny+1
             Do  ix=1,nx+1     
                k = ix + (iy-1)*nx + (iz-1)*nx*ny
                j = j +1
                If (ix.Eq.nx+1) Then
                   k = k - nx
                End If
                If  (iy.Eq.ny+1) Then
                   k = k - (nx*ny)
                End If
                Write (8,'(I7,8X,4(1X,1pE11.4))')j,Vx(k)/theta(k),Vy(k)/theta(k),Vz(k)/theta(k),theta(k)
                Write (9,'(I7,8X,3(1X,1pE11.4))')j,xGrid(k),yGrid(k),zGrid(k)
                Write (10,'(I7,8X,1(1X,1pE11.4))')j,theta(k)
             End Do
          End Do
       End Do
    Else
       l=0
       Do j=1,nz
          Do i=1,nx*ny
             k = nx*ny*nz - j*nx*ny +i
             l = l+1
             Write (8,'(I7,8X,4(1X,1pE11.4))')l,Vx(k)/theta(k),Vy(k)/theta(k),Vz(k)/theta(k),theta(k)
             Write (9,'(I7,8X,3(1X,1pE11.4))')l,xGrid(k),yGrid(k),zGrid(k)
             Write (10,'(I7,8X,1(1X,1pE11.4))')l,theta(k)
          End Do
       End Do
    Endif


    Close(8)
    Close(9)
    Close(10)
    Return
  End Subroutine PartraceOut
  !****************************************************************************
  !> UNUSED AT THE MOMENT 
  Subroutine OutFEM_VTK_Sub(kOut)
    Use typedef
    Use GridData
    Use SolData
    Use WatFun
    Implicit None

    Integer(ap), Intent(in)::kout
    Integer(ap)::ie,i,ise,global
    Character file*20
    !> \param kout current number for the FEM output

    Write (file,'(A20)')'out/outfem_sub00.vtk'
    If (kout.Lt.10) Then
       Write (file(16:16),'(I1)') kout
    Else
       Write (file(15:16),'(I2)') kout
    Endif


    Open (UNIT=8,FILE=file,STATUS='UNKNOWN')
    Write (8,'(A26)')'# vtk DataFile Version 3.0' 
    Write (8,'(A12)')'model R-SWMS'
    Write (8,'(A5)')'ASCII'
    Write (8,'(A25)')'DATASET UNSTRUCTURED_GRID'

    Write (8,'(/A6,1X,I7, 1X, A5)')'POINTS', nPt, 'float'
    Do i=1,nPt
       Write (8,'(3(1X,1pE12.5))',advance="no") xgrid(i), ygrid(i) ,zGrid(i)
       Write (8,*)''
    End Do


    Write (8,'(/A5,1X,I7,1X,I8)')'CELLS', nElm*5, 5*nElm*5
    Do iE=1,nElm
       Do iSE=1,5
          Write (8,'((1X,I6))',advance="no") 4
          Do global=1,4
             Write (8,'((1X,I6))',advance="no") elmnod(iL(global,iSE,subN(iE)), iE)-1
          End Do
          Write (8,*)''
       End Do
    Enddo

    Write (8,'(/A10,1X,I7)')'CELL_TYPES', nElm*5
    Do iE=1,nElm
       Do iSe=1,5
          Write (8,'(8(1X,I2))',advance="no") 10
       End Do
    Enddo

    Write (8,*)''
    Write (8,'(A11,1X,I7)')'POINT_DATA ', nPt
    Write (8,'(A30)')'SCALARS pressurehead float '
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE15.6)',advance="no") hnew(1:nPt)

    Write (8,*)''
    Write (8,'(A30)')'SCALARS wc float '
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE15.6)',advance="no") theta(1:nPt)

    Write (8,*)''
    Write (8,'(A30)')'SCALARS conc float '
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE15.6)',advance="no") conc(1:nPt)

    Write (8,*)''
    Write (8,'(A30)')'SCALARS sink float '
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE15.6)',advance="no") sink(1:nPt)

    Write (8,*)''
    Write (8,'(A30)')'SCALARS csink float '
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Write (8,'(1X,1pE15.6)',advance="no") csink(1:nPt)

    Write (8,*)''
    Write (8,'(A11,1X,I7)')'CELL_DATA ', nElm*5
    Write (8,'(A30)')'SCALARS subN float '
    Write (8,'(A20)')'LOOKUP_TABLE default'
    Do iE=1,nElm
       Do iSe=1,5
          Write (8,'(1X,I7)',advance="no") subN(iE)
       End Do
    Enddo


    Close (8)
    Return
  End Subroutine OutFEM_VTK_Sub
  !****************************************************************************
  !> generates .vtk outputs for particles within roots, in case of hormone transport; ParRoot.XX 
  Subroutine OutParticleVTK(kout,t,ipl)
    Use typedef
    Use RootData, Only: seglen,nrec,ngrow,nbr,xg,yg,zg,xs,ys,zs,xplant,yplant,irecpr,ibrseg,ibrgrw,nplant
    Use DoussanMat, Only:  nsub,transroot
    Use GridData, Only: dxgrid,dygrid,nx,ny,continu
    Use SoluteRootMat, Only: Particle,firstP,totalParticleNum

    Implicit None

    Integer(ap) :: irec,kout,ipl
    Integer(ap):: igrow,irecn,j
    Real(dp), Intent(in) :: t
    Real(dp) :: xr(nrec(nplant)+ngrow(nplant)), yr(nrec(nplant)+ngrow(nplant)), zr(nrec(nplant)+ngrow(nplant))
    Real(dp) :: seg_tipx(nrec(nplant)),seg_tipy(nrec(nplant)),seg_tipz(nrec(nplant))
    Character file*29 
    Type(Particle) ,Pointer  :: P
    !> \param t current simulation time
    !> \param rOuR current number for the root output


       Write (file,'(A18)')'out/vtk/ParRootXX.'


       If (ipl.Lt.10) Then
          Write (file(16:16),'(I1)') 0
          Write (file(17:17),'(I1)') ipl
       Else
          Write (file(16:17),'(I2)') ipl
       End If
       If (kout.Lt.10) Then
          Write (file(19:19),'(I1)') kout
          Write (file(20:24),'(A4)') '.vtk'
       Elseif (kout.Lt.100) Then
          Write (file(19:20),'(I2)') kout
          Write (file(21:25),'(A4)') '.vtk'
       Else
          Write (file(19:21),'(I3)') kout
          Write (file(22:26),'(A4)') '.vtk'
       Endif

       Open (UNIT=8,FILE=file,STATUS='UNKNOWN')

       If (continu) Then

          ! write out points root segments
          Do  irec=1,nrec(ipl)
             xr(irec) = xs(irec,ipl)+xplant(ipl)+1.*transroot(irec,1,nsub(irec,ipl),ipl)*dxgrid*nx
             yr(irec) = ys(irec,ipl)+yplant(ipl)+1.*transroot(irec,2,nsub(irec,ipl),ipl)*dygrid*ny
             zr(irec) = zs(irec,ipl)
             !WRITE (8,'(1X,3(1X,E12.5))',advance="no") xr(irec),yr(irec),zs(irec)
             !WRITE (8,*)''
          End Do

          ! write out points foot tips
          Do  igrow=1,ngrow(ipl)
             xr(nrec(ipl)+igrow) = xg(igrow,ipl)+xplant(ipl)+1.*transroot(nrec(ipl)+igrow,1,nsub(1,ipl),ipl)*dxgrid*nx
             yr(nrec(ipl)+igrow) = yg(igrow,ipl)+yplant(ipl)+1.*transroot(nrec(ipl)+igrow,2,nsub(1,ipl),ipl)*dygrid*ny
             zr(nrec(ipl)+igrow) = zg(igrow,ipl)
             !WRITE (8,'(1X,3(1X,E12.5))',advance="no") xr(nrec+igrow),yr(nrec+igrow),zg(igrow)
             !WRITE (8,*)''
          End Do



       Else !if not continue

          ! points segments
          Do  irec=1,nrec(ipl)
             xr(irec) = xs(irec,ipl)+xplant(ipl)
             yr(irec) = ys(irec,ipl)+yplant(ipl)
             zr(irec) = zs(irec,ipl)
          End Do

          ! points tips
          Do  igrow=1,ngrow(ipl)
             xr(nrec(ipl)+igrow) = xg(igrow,ipl)+xplant(ipl)
             yr(nrec(ipl)+igrow) = yg(igrow,ipl)+yplant(ipl)              
             zr(nrec(ipl)+igrow) = zg(igrow,ipl)        
          End Do

       End If !continu 

       Do j=1,nrec(ipl) !segments
          If(irecpr(j,ipl).Ne.0) Then
             seg_tipx(irecpr(j,ipl)) = (xr(j)-xr(irecpr(j,ipl)))/seglen(irecpr(j,ipl),ipl)
             seg_tipy(irecpr(j,ipl)) = (yr(j)-yr(irecpr(j,ipl)))/seglen(irecpr(j,ipl),ipl)
             seg_tipz(irecpr(j,ipl)) = (zr(j)-zr(irecpr(j,ipl)))/seglen(irecpr(j,ipl),ipl)
          End If
       Enddo

       Do j=1,nbr(ipl) !tips
          irecn=nrec(ipl)
          Do While (ibrseg(irecn,ipl).Ne.j)
             irecn=irecn-1
          End Do
          Do igrow=1,ngrow(ipl)
             If (ibrgrw(igrow,ipl)==j) Then

                seg_tipx(irecn) = (xr(nrec(ipl)+igrow)-xr(irecn))/seglen(irecn,ipl)
                seg_tipy(irecn) = (yr(nrec(ipl)+igrow)-yr(irecn))/seglen(irecn,ipl)
                seg_tipz(irecn) = (zr(nrec(ipl)+igrow)-zr(irecn))/seglen(irecn,ipl)

             Endif
          End Do
       End Do

       Write (8,'(A26)')'# vtk DataFile Version 3.0' 
       Write (8,'(A12)')'model R-SWMS'
       Write (8,'(A5)')'ASCII'
       Write (8,'(A25)')'DATASET UNSTRUCTURED_GRID'
       Write (8,'(A7,1X,I8,1X,A5)')'POINTS ',totalParticleNum  , 'float'   
       p => firstP
       Do While(Associated(P))
          Write (8,'(3(1X,1pE11.4))',advance="no") xr(p%segnum)+seg_tipx(p%segnum)*(seglen(p%segnum,ipl)-p%position), &
               yr(p%segnum)+seg_tipy(p%segnum)*(seglen(p%segnum,ipl)-p%position), &
               zr(p%segnum)+seg_tipz(p%segnum)*(seglen(p%segnum,ipl)-p%position)
          Write (8,*)''
          P => P%NEXT
       End Do

       Write (8,'(A8,1X,I8,1X,I8)')'CELLS ',totalParticleNum ,totalParticleNum*2
       Do j=1,totalParticleNum 
          Write (8,'(2(1X,I7))',advance="no") 1,j-1
          Write (8,*)''
       End Do

       Write (8,'(A16,1X,I8)')'CELL_TYPES ',totalParticleNum 
       Do j=1,totalParticleNum 
          Write (8,'(1X,I7)',advance="no") 1
       End Do
       Write (8,*)''

       Write (8,'(A11,1X,I8)')'POINT_DATA ', totalParticleNum
       Write (8,'(A17)')'SCALARS ID float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       p=> firstP
       Do While(Associated(P))
          Write (8,'((1X,I8))',advance="no") p%ID
          P => P%NEXT
       End Do
       Write (8,*)''

       Write (8,'(A19)')'SCALARS Mass float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       p=> firstP
       Do While(Associated(P))
          Write (8,'((1X,1pE11.4))',advance="no") p%mass
          P => P%NEXT
       End Do
       Write (8,*)''

       Write (8,'(A21)')'SCALARS SegNum float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       p=> firstP
       Do While(Associated(P))
          Write (8,'((1X,I8))',advance="no") p%segnum
          P => P%NEXT
       End Do
       Write (8,*)''

       Write (8,'(A21)')'SCALARS PartAge float '
       Write (8,'(A20)')'LOOKUP_TABLE default'
       p=> firstP
       Do While(Associated(P))
          Write (8,'((1X,1pE11.4))',advance="no") t-p%partOrig
          P => P%NEXT
       End Do
       Write (8,*)''

       Close(8)

       Return
  End Subroutine OutParticleVTK
  !****************************************************************************
End Module Output
