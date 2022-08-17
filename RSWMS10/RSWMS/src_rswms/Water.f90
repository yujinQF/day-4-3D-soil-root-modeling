!>File Water.f90 RSWMS9 - 22/04/2020

! Main subroutine WATER calculates water flow in the soil
  Subroutine WATER(t,dt,dtOpt,tOld,ReDo,IAD,IADN,IADD)
    Use Typedef
    Use ParamData, Only: maxbnd,maxplant,iter,iter_tot,i_noConv
    Use tmctrl, Only: dtMin
    Use GridData, Only : epslonPH,epslonR,epslonS,RelEps,factorRelEps,&
         nPt,itMax,sink_cube,sink_cubeOld,HElm,HElmOld,epslonWC,denomSe
    Use MatData
    Use DoussanMat, Only : PHr,PHrOld,sinkRold,SinkR,PHrtemp, &
         savelast,stressBC,JoutrOld,Joutr
    Use RootData, Only :nrec,lCou,lDou,lFed,ldJvL,nplant,lno_RWU,lno_root_growth!,nrec_m
    Use WatFun
    Use SolData, Only: lCelia, hOld, hTemp, hNew, theta, theta_old,Cap,lRoot_explicit
    Use RhizoData
    Use Doussan, Only: SolveRoot, SetBCroot
    Use MPIutils, Only: stop_program
    Use Sink, Only: SetSnk
    Implicit None

    Integer(ap) :: IAD(maxbnd,nPt),IADN(nPt),IADD(nPt),iter_root,iter_soil!,irecn
    Integer(ap) :: ipl,BCtp,i,Passage,Maxpassage
    Real(dp),Intent(in):: told
    Real(dp),Intent(inout)::dt,t,dtopt
    Real(dp):: tol_R,tol_S,tol_PH,tol_WC,ThNew(1:nPt)
    Real(dp):: deltaPHr,deltaPHint,deltaS,BCr,deltaTH!,deltaJoutr
    Real(dp):: deltaPHs, maxPhrtot, maxPhr(1:maxplant)
!    Real(dp):: Div,OLD(1:nrec_m),NEW(1:nrec_m),DiffOLDNEW(1:nrec_m)
    Real(dp), Allocatable,Dimension (:) ::sol
    Logical :: ItCrit,stressBC_old,Explic=.false.,l_updtD=.False.
    Logical,Intent(out) :: ReDo

    Allocate(sol(nPt),B(nPt))
    !> \param t current time
    !> \param dt time step size
    !> \param dtOpt optimal time step size calculated in the time step control TmCont
    !> \param tOld time before last time step
    !> \param ReDo logical to rerun the model from a given time step
    !> \param IAD adjacency matrix (nodal connections)
    !> \param IADN number of adjacent nodes in IAD (self-inclusive)
    !> \param IADD position of diagonal in adjacency matrix

!----------------------- initialisation ------------------------------------

   iter_root=0
    ReDo=.False.
    theta_old=theta
    If (lDou) Then
       Do ipl=1,nplant
          JoutrOld(0:nrec(ipl),ipl)=Joutr(0:nrec(ipl),ipl)
          PHrOld(1:nrec(ipl)+1,ipl)=PHr(1:nrec(ipl)+1,ipl)
          sinkRold(0:nrec(ipl),ipl)=sinkR(0:nrec(ipl),ipl)
       Enddo
    Endif
    iter_soil = iter
    tol_PH=epslonPH
    tol_R=epslonR
    tol_S=epslonS
    tol_WC=epslonWC

240 Continue
    iter = 0

Maxpassage=5
Passage=0
! Explicit root solving
241  If (lRoot_explicit) Then
Passage=Passage+1
       iter_root=iter_root+1
       Write(*,'(a)',advance='no') 'R'
       sink_cubeOld=sink_cube
       Call SolveRoot(t,dt,.False.,iter_root)
     EndIf

! ----------------------- Start of FEM iteration loop -------------------------
250 Continue

    If (lDou) stressBC_old=stressBC
! NonEquilibrium in the rhizosphere
    If (lRhizo .And. RhizoModel .Eq. 3) Then ! Only Relevant for the Dynamic Scenario 
       Do i=1,nPt
          If (Rnorm(i) .Ne. 1.) Then
             Call thtNonEq(hnew(i), dt, i)
          Else
             thetaNonEq(i) = 0.0
             thetaTot(i) = theta(i)
             hEqRhizo(i) = hNew(i)
             tauTht(i) = 1.0
          Endif
       Enddo
    Endif

! Calculate theta at k, Conductivity at k and capacity	
    Call SetMat(1,Explic)
    If (lCelia) Then
       Call ResetCelia(dt,IAD,IADN,IADD)
     Else
        Call Reset(dt,IAD,IADN,IADD)
    Endif
    Call Dirich(IADD)
    Call SolveIlut(sol)
    hTemp=hNew
    hNew=sol
    ThNew=theta+Cap*(hNew-hTemp)/denomSe! Updated theta at k+1 for error quantification
	
    iter=iter+1
    iter_tot=iter_tot+1
 
! Test for convergence:
    ItCrit=.False.! Stop iterating when ItCrit remains false
! When relative error is on, we take the largest tolerance on the maximum PH
    If (RelEps) Then
       If (epslonPH .Lt. Abs(Maxval(hNew(1:nPt)/factorRelEps))) tol_PH= Abs(Maxval(hNew(1:nPt)/factorRelEps))
    Endif

! Error in PH and WC
    deltaTH=Maxval(Abs(ThNew(1:nPt)-theta(1:nPt))) ! theta is k and ThNew is k+1
    deltaPHs=Maxval(Abs(hNew(1:nPt)-hTemp(1:nPt))) ! htemp is k and hNew is k+1

! Check solution accuracy
    If ((deltaPHs.GT.tol_PH).OR.(deltaTH.GT.tol_WC)) Then ! Strictest condition wins
       If (deltaPHs.GT.tol_PH) Then
          Write (*,'(a)',advance='no') 'h'
       Else
          Write (*,'(a)',advance='no') 'w'
       EndIf
       ItCrit=.True.! Error on the soil solution too large
    Endif
 !!Ax 2022 
   If (lRoot_explicit) Then
      !    Do ipl=1,nplant
      !       PHrTemp(1:nrec(ipl)+1,ipl)=PHr(1:nrec(ipl)+1,ipl)!PHrtemp is at k
      !    Enddo
      !    If (RelEps) Then 
      !       Do ipl=1,nplant
      !          maxPhr(ipl)=Maxval(PHr(1:nrec(ipl)+1,1:nplant)/factorRelEps)
      !       Enddo
      !       maxPhrtot=Maxval(maxPhr(1:nplant))
      !       tol_R=epslonR
      !       If (epslonR .Lt. Abs(maxPhrtot)) tol_R= Abs(maxPhrtot) ! we take the largest tolerance
      !    Endif
! Check for all nodes if PHr is lower than relative or absolute criterion
      !    maxPhr=0
      !    Do ipl=1,nplant
      !       maxPhr(ipl)=Maxval(Abs(PHr(1:nrec(ipl)+1,ipl)-PHrTemp(1:nrec(ipl)+1,ipl)))
      !    Enddo
      !    deltaPHr=Maxval(maxPhr(1:nplant))
      !    deltaS=Maxval(Abs(sink_cubeOld-sink_cube))
      !    l_updtD=.False.
      !    If (deltaPHr.Gt.tol_R) Then
      !       Write(*,'(a)',advance='no') 'r'
      !       ItCrit=.True.! no convergence for root-> rerun soil and root.
      !       l_updtD=.True.
      !    Elseif (deltaS.Gt.tol_S) Then
      !       Write(*,'(a)',advance='no') 's'
      !       ItCrit=.True.! no convergence for root-> rerun soil and root.
      !       l_updtD=.True.
      !    Endif
      !If (iter.Lt.itMax) Then
      !    If (ItCrit) GOTO 241
      !Endif
      If (Passage.LT.Maxpassage) GOTO 241
   Endif


! In case there is convergence for soil (itCrit==.false.), we check RWU is still OK
    If (.Not.(itCrit).and.(.not.(lRoot_explicit))) Then 
       sink_cubeOld=sink_cube
       If (lDou) Then
          Do ipl=1,nplant
             PHrTemp(1:nrec(ipl)+1,ipl)=PHr(1:nrec(ipl)+1,ipl)!PHrtemp is at k
          Enddo
! Solving Doussan equation -> PHr at k+1
             iter_root=iter_root+1
             Call SolveRoot(t,dt,.False.,iter_root)
! Relative or absolute error
          If (RelEps) Then 
             Do ipl=1,nplant
                maxPhr(ipl)=Maxval(PHr(1:nrec(ipl)+1,1:nplant)/factorRelEps)
             Enddo
             maxPhrtot=Maxval(maxPhr(1:nplant))
             tol_R=epslonR
             If (epslonR .Lt. Abs(maxPhrtot)) tol_R= Abs(maxPhrtot) ! we take the largest tolerance
          Endif
! Check for all nodes if PHr is lower than relative or absolute criterion
          maxPhr=0
          Do ipl=1,nplant
             maxPhr(ipl)=Maxval(Abs(PHr(1:nrec(ipl)+1,ipl)-PHrTemp(1:nrec(ipl)+1,ipl)))
          Enddo
          deltaPHr=Maxval(maxPhr(1:nplant))
          deltaS=Maxval(Abs(sink_cubeOld-sink_cube))
          l_updtD=.False.
          If (deltaPHr.Gt.tol_R) Then
             Write(*,'(a)',advance='no') 'r'
             ItCrit=.True.! no convergence for root-> rerun soil and root. !Ax 2021
             l_updtD=.True.
          Elseif (deltaS.Gt.tol_S) Then
             Write(*,'(a)',advance='no') 's'
             ItCrit=.True.! no convergence for root-> rerun soil and root. !Ax 2021
             l_updtD=.True.
          Endif
       Elseif (lFed.Or.lCou) Then 
          If (ldJvL) Then
             HElmOld=HElm
          Endif
          Call SetBCroot(t,BCr,BCtp)
          Call SetSnk(t,BCr,BCtp)! estimate sink_cube
          If (ldJvL) Then
             deltaPHint=Maxval(Abs(HElmOld-HElm))
             If (deltaPHint.Gt.tol_R) Then
                Write(*,'(a)',advance='no') '*'
                ItCrit=.True.! no convergence for soil-root interface water potential-> rerun soil and root.
             Endif
          Endif
          deltaS=Maxval(Abs(sink_cubeOld-sink_cube))
          If (deltaS.Gt.tol_S) Then
             Write(*,'(a)',advance='no') '+'
             ItCrit=.True.! no convergence for root-> rerun soil and root.
          Endif
       Endif
    Endif

    If(Explic) ItCrit=.false. !no check of convergence anymore -> end of this routine
    If ((savelast) .And. (.Not.(iTCrit)) ) Return

    If (ItCrit) Then
       If (iter.Lt.itMax) Then
          Goto 250
! No convergence and minimum dt ->explicit solution
       Else If (dt.Le.dtMin) Then
          Write(*,'(a)',advance='no') 'No convergence-> explicit solution'
          Explic=.true.
          i_noConv=i_noConv+1
          hNew =hOld
          hTemp=hOld
          Goto 250
! No convergence -> dt reduction
       Else 
          Write (*,'(a)') '!'
          hNew =hOld
          hTemp=hOld
          If (lDou) Then
             stressBC=stressBC_old
             Do ipl=1,nplant
                 PHr(1:nrec(ipl)+1,:)=PHrOld(1:nrec(ipl)+1,:)
             End do
          Endif
          dt=Max(dt/3._dp,dtMin)
          dtOpt=dt
          t=tOld+dt
          Write (*,'(/,a,1pe12.5,a,1pe11.3,a,i6)') ' t = ',t,' dt = ',dt,' iter = ',iter_tot
          Write (*,'(a)',advance='no') '->'
! When there is root growth or solute transport, dt impacts other processes.
          If ((lno_RWU).OR.(lno_root_growth)) Then
              Goto 240              
          else
              ReDo=.True.
              Goto 200
           Endif
       Endif
    Endif
!iteration is OK just make a last update of soil WC& PH with the right sink term
    Call WatInf(dt)
! Update theta 
    Call SetMat(1,Explic)
! --------------------- End of FEM iteration loop -----------------------------

200 Deallocate(sol,B)
    Return

  End Subroutine WATER

  
  
  !*******************************************************************************
  !> solve the linear equation system Ax=b - preconditioner and iterative solver can be chosen here    
  Subroutine SOLVEILUT(sol)
    Use Typedef
    Use MatData
    Use GridData, Only : nPt

    Integer(ap) :: lfil,maxits=1000,ierr
    Integer(ap) :: ipar(16),nwk
    Real(dp) :: droptol
    Real(dp):: fpar(16)
    Real(dp), Intent(out) ::sol(nPt)

    Integer(ap) , Allocatable, Dimension(:) :: jwork,ju
    Real(dp), Allocatable,Dimension (:) ::work,xran,ww
    Real(dp), Allocatable,Dimension (:,:) ::vv
    External cgnr, fom, ilut,ilu0,ulid
    External cg,bcg,dbcg,bcgstab,tfqmr,gmres,fgmres,dqgmres

    Allocate (jwork(2*nPt))
    Allocate (xran(nPt))
    Allocate (work(nPt+1))
    Allocate (vv(nPt,20))
    Allocate (ww(nPt*40))
    Allocate (ju(nPt))
    Allocate (alu(numnz))
    Allocate (jlu(numnz))

    !> \param sol solution vector of the linear equation system -> pressure head

    !set the parameters for the iterative solvers
    ipar(2) = 1      !1 is left preconditioning,2 right,3 both, 0 = no precon
    ipar(3) = 1      !1 -> convergence test based on residualnorms; 2->based on change in approximate solution
    ipar(4) = nPt*40      !# elements in array w
    ipar(5) = 16      !size of krylov subspace
    ipar(6) = maxits      !max number of matrix-vector operations
    fpar(1) = 1.0E-8   !relative tolerance
    fpar(2) = 1.0E-12  !absolute tolerance
    lfil = 1
    droptol =1.0E-2
    nwk = 2*numNZ

    !preconditioner
    !IF (MOD(time_step,10).EQ.0) THEN
    Call ilut (nPt,A_dparse,JCOL,IROW,lfil,droptol,alu,jlu,ju,nwk,work,jwork,ierr)
    !END IF
    time_step = time_step +1
    !solve
    xran=0.0 !initial guess
    Call runrc(nPt,B,sol,ipar,fpar,ww,xran,A_dparse,JCOL,IROW,alu,jlu,ju,bcgstab)

    Deallocate (jwork)
    Deallocate (xran)
    Deallocate (work)
    Deallocate (vv)
    Deallocate (ww)
    Deallocate (ju)
    Deallocate (alu)
    Deallocate (jlu)

  End Subroutine SOLVEILUT
  !**************************************************************************
  !>     the actual tester. It starts the iterative linear system solvers
  !>     with a initial guess suppied by the user.
  !>
  !>     The structure {au, jau, ju} is assumed to have the output from
  !>     the ILU* routines in ilut.f.   
  Subroutine runrc(n,rhs,sol,ipar,fpar,wk,guess,a,ja,ia,&
       au,jau,ju,solver)
    Use Typedef
    Implicit None

    Integer(ap) :: n,ipar(16),ia(n+1),ja(*),ju(*),jau(*)
    Real(dp) :: fpar(16),rhs(n),sol(n),guess(n),wk(*),a(*),au(*)
    External solver

    !c-----------------------------------------------------------------------
    !c     local variables
    Integer(ap) i, its
    Real(dp) res, dnrm2
    !c     external dtime
    External dnrm2
    Save its,res

    !ipar(2) can be 0, 1, 2, please don´t use 3
    !      IF (ipar(2).GT.2) THEN
    !         PRINT *, 'I can not do both left and right preconditioning.'
    !         RETURN
    !      ENDIF

    !normal execution

    its = 0
    res = 0._dp
    Do i = 1, n
       sol(i) = guess(i)
    Enddo
    ipar(1) = 0

10  Call solver(n,rhs,sol,ipar,fpar,wk)

    !output the residuals

    If (ipar(7).Ne.its) Then
       its = ipar(7)
    Endif
    res = fpar(5)
    !c
    If (ipar(1).Eq.1) Then
       Call amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
       Goto 10
    Else If (ipar(1).Eq.2) Then
       Call atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
       Goto 10
    Else If (ipar(1).Eq.3 .Or. ipar(1).Eq.5) Then
       Call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
       Goto 10
    Else If (ipar(1).Eq.4 .Or. ipar(1).Eq.6) Then
       Call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
       Goto 10
    Else If (ipar(1).Le.0) Then
       If (ipar(1).Eq.0) Then
          !            WRITE (*,'(''*'',$)')
          !Iterative solver has satisfied convergence test.
       Else If (ipar(1).Eq.-1) Then
          Print *, 'Iterative solver has iterated too many times.'
       Else If (ipar(1).Eq.-2) Then
          Print *, 'Iterative solver was not given enough work space.'
          Print *, 'The work space should at least have ', ipar(4),&
               ' elements.'
       Else If (ipar(1).Eq.-3) Then
          Print *, 'Iterative sovler is facing a break-down.'
       Else
          Print *, 'Iterative solver terminated. code =', ipar(1)
       Endif
    Endif
    Return
  End Subroutine runrc
  !c-----end-of-runrc
  Function distdot(n,x,ix,y,iy)
    Use Typedef

    Integer(ap) :: n, ix, iy
    Real(dp) :: distdot, x(*), y(*), ddot

    External ddot
    distdot = ddot(n,x,ix,y,iy)
    Return
  End Function distdot
!*********************************************************************
!> sets the boundary conditions in the matrix and in the right-hand side vector 
!> of the linear equation system
  Subroutine SetBC(t)
    Use typedef
    Use BoundData
    Use GridData
    Use SolData, Only: hnew,Kode,Kcell,nmat
    Use CumData
    Implicit None

    Real(dp) :: t,head,Qtot,wtot,m,b
    Real(dp) :: head1,head2,head3,head4,head5,head6,head7,head8,head9 ! useful for rhizotron
    Real(dp) , Dimension (:) :: volflw(1+(homogene-1)*mxBcCh),irrigflw(nBcPts)
    Integer(ap) :: i,jj,n


    Q = 0._dp
    ! calculate current value for Q-BC´s, if any:
    ! first check kode is -1 (check for just one node). If kode is -1, then follow loop and calculate flow
    If (noBCflux) Goto 23 !new tp avoid error message
    If (Kode(iBCPt(1)).Eq.-1) Then
       i=0
22     i=i+1   
       If (t.Ge.tQbcCh(nQbcCh)) Then
          If (homogene.Eq.1) Then
             volflw=Qbc(nQbcCh,1)
          Else
             volflw(1:mxBcCh)=Qbc(nQbcCh,:)
       Endif
       Else
          If (t.Gt.tQbcCh(i)) Goto 22
             If (homogene.Eq.1) Then
                volflw=Qbc(i-1,1)+(t-tQbcCh(i-1))/(tQbcCh(i)-tQbcCh(i-1))*(Qbc(i,1)-Qbc(i-1,1))
             Else
                volflw(1:mxBcCh)=Qbc(i-1,:)+(t-tQbcCh(i-1))/(tQbcCh(i)-tQbcCh(i-1))*(Qbc(i,:)-Qbc(i-1,:))
             Endif
          Endif
    Endif


! Calculate current value for I-BC´s, if any:
23    If (nIbcCh.Eq.0) Goto 30  !mark 25 was removed here because it gave a warning
    irrigflw(:)=0
! Loop over irrigators
    n=0
    Do i=1,nBcPts
       If (Kode(iBcPt(i)).Eq.-3) Then
          n=n+1      !irrigator number
          If(t.Ge.tIBcCh(n,nIBcCh)) Then
             irrigflw(i)=Ibc(n,nIBcCh)
          Else
             jj=nIbcCh
27           jj=jj-1
             If((t.Lt.tIBcCh(n,jj)).And.(jj.Gt.1)) Goto 27
             m=(Ibc(n,jj+1)-Ibc(n,jj))/(tIBcCh(n,jj+1)-tIBcCh(n,jj))
             b=Ibc(n,jj)-tIBcCh(n,jj)*m
             irrigflw(i)=m*t+b
          Endif
       Endif
    Enddo

! Calculate current value for h-BC´s, if any:
30  If (nhbcCh.Eq.0) Goto 40
    i=1
33  i=i+1

    If (i.Gt.nhbcCh) Then
       If ((geom.Eq.4) .And. (nmat.Eq.2)) Then ! rhizotron
       head1=hbc1(nhbcCh)
       head2=hbc2(nhbcCh)
       head3=hbc3(nhbcCh)
       head4=hbc4(nhbcCh)
       head5=hbc5(nhbcCh)
       head6=hbc6(nhbcCh)
       head7=hbc7(nhbcCh)
       head8=hbc8(nhbcCh)
       head9=hbc9(nhbcCh)
       Else 
       head=hbc(nhbcCh)
       Endif
    Else
       If (t.Gt.thbcCh(i)) Goto 33
          If ((geom.Eq.4) .And. (nmat.Eq.2))  Then ! rhizotron
          head1=hbc1(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc1(i)-hbc1(i-1))
          head2=hbc2(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc2(i)-hbc2(i-1))
          head3=hbc3(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc3(i)-hbc3(i-1))
          head4=hbc4(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc4(i)-hbc4(i-1))
          head5=hbc5(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc5(i)-hbc5(i-1))
          head6=hbc6(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc6(i)-hbc6(i-1))
          head7=hbc7(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc7(i)-hbc7(i-1))
          head8=hbc8(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc8(i)-hbc8(i-1))
          head9=hbc9(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc9(i)-hbc9(i-1))
          Else
          head=hbc(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc(i)-hbc(i-1))
          Endif
    Endif

! Calculate QBC*surface
40  Qtot=0
    wtot=0
    Do i=1,nBCPts
       If ((Kode(iBCPt(i)).Eq.+1).Or.(Kode(iBCPt(i)).Eq.+2)) Then 
           hnew(iBCPt(i))=head
       Elseif (Kode(iBCPt(i)).Eq.+3) Then ! rhizotron
           If (Kcell(iBCPt(i)).Eq.1) Then
              hnew(iBCPt(i))=head1
          Elseif (Kcell(iBCPt(i)).Eq.2) Then
              hnew(iBCPt(i))=head2
          Elseif (Kcell(iBCPt(i)).Eq.3) Then
              hnew(iBCPt(i))=head3
          Elseif (Kcell(iBCPt(i)).Eq.4) Then
              hnew(iBCPt(i))=head4
          Elseif (Kcell(iBCPt(i)).Eq.5) Then
              hnew(iBCPt(i))=head5
          Elseif (Kcell(iBCPt(i)).Eq.6) Then
              hnew(iBCPt(i))=head6
          Elseif (Kcell(iBCPt(i)).Eq.7) Then
              hnew(iBCPt(i))=head7
          Elseif (Kcell(iBCPt(i)).Eq.8) Then
              hnew(iBCPt(i))=head8
          Elseif (Kcell(iBCPt(i)).Eq.9) Then
              hnew(iBCPt(i))=head9
          Endif
       Elseif (Kode(iBCPt(i)).Eq.-1) Then
          If (homogene.Eq.1) Then
             Q(iBCPt(i))=volflw(1)*Width(iBCPt(i))
          Else
             Q(iBCPt(i))=volflw(i)*width(iBCPt(i)) 
          Endif
          Qtot=Qtot+Q(iBCPt(i))
          wtot=wtot+Width(iBCPt(i))
       Elseif (Kode(iBCPt(i)).Eq.-3) Then
          Q(iBCPt(i))=irrigflw(i)*dxgrid*dygrid
          Qtot=Qtot+Q(iBCPt(i))
          wtot=wtot+Width(i)
       Endif
    End Do
! Calculate current value for solute transport boundary conditions, if any:
    If (nCBnd1.Eq.0) Goto 50
    i=0
55  i=i+1
    If (i.Ge.nCBnd1) Then
       cBound(1)=CBnd1(nCBnd1)
    Else
       If (t.Gt.tCBnd1(i)) Goto 55
       cBound(1)=CBnd1(i-1)+(t-tCBnd1(i-1))/(tCBnd1(i)-tCBnd1(i-1))*(CBnd1(i)-CBnd1(i-1))
    Endif
50  If (nCBnd2.Eq.0) Goto 60
    i=0
66  i=i+1
    If (i.Ge.nCBnd2) Then
       cBound(2)=CBnd2(nCBnd2)
    Else
       If (t.Gt.tCBnd2(i)) Goto 66
       cBound(2)=CBnd2(i-1)+(t-tCBnd2(i-1))/(tCBnd2(i)-tCBnd2(i-1))*(CBnd2(i)-CBnd2(i-1))
    Endif
60  Return
  End Subroutine SetBC
  !*********************************************************************************
  !> calculates the width at the boundary faces similar to the SWMS_3D model    
  Subroutine CalcWidthSWMS_new
    Use Typedef
    Use GridData
    Use DomData
    Use BoundData, Only: xqmax1,xqmin2,qfun
    Use SolData, Only: MatNum
    Implicit None

    Integer(ap):: iB,ix,iy,iCord,iDiv,iRest,i,bord1,bord2
    Real (dp):: WidthX,WidthY

    iB=nPt-(nx*ny)
    i=0
    Do iy=1,ny
       Do ix=1,nx
          i=i+1
          iB=iB+1
          If(qfun.Eq.4) Then  !for cylindrical/arbitrary setup              
             !X direction
             If(ix.Ne.1 .And. ix.Ne.nx) Then
                If(MatNum(i-1).Eq.1 .And. MatNum(i+1).Eq.1) Then
                   WidthX =  xCol(ix+1)-xCol(ix-1)
                Else If(MatNum(i-1).Ne.1) Then !left border
                   WidthX = xCol(ix+1)-xCol(ix)
                Else If(MatNum(i+1).Ne.1) Then !right border
                   WidthX = xCol(ix)-xCol(ix-1)
                End If
             Else If(ix.Eq.1) Then
                WidthX = xCol(ix+1)-xCol(1)
             Else If(ix.Eq.nx) Then
                WidthX = xCol(nx)-xCol(nx-1)
             End If
             !Y direction
             If(iy.Ne.1 .And. iy.Ne.ny) Then
                If(MatNum(i-ny).Eq.1 .And. MatNum(i+ny).Eq.1) Then
                   WidthY =  yCol(iy+1)-yCol(iy-1)
                Else If(MatNum(i-ny).Ne.1) Then !lower border
                   WidthY = yCol(iy+1)-xCol(iy)
                Else If(MatNum(i+ny).Ne.1) Then !upper border
                   WidthY = yCol(iy)-yCol(iy-1)
                End If
             Else If(iy.Eq.1) Then
                WidthY = yCol(iy+1)-yCol(1)
             Else If(iy.Eq.ny) Then
                WidthY = yCol(ny)-yCol(ny-1)
             End If

          Else               
             If(ix.Ne.1.And.ix.Ne.nx) Then
                If(qfun.Eq.3) Then  !for split setup introducing a second and third border at xqmax1,xqmin2
                   bord1=NINT((xqmax1-xmin)/dxgrid)
                   bord2=NINT((xqmin2-xmin)/dxgrid+1)
                   If(ix.Ne.bord1.And.ix.Ne.bord2) Then
                      WidthX=(xCol(ix+1)-xCol(ix-1))
                   Else If(ix.Eq.bord1) Then
                      WidthX=xcol(bord1)-xcol(bord1-1)
                   Else If(ix.Eq.bord2) Then
                      WidthX=xcol(bord2+1)-xcol(bord2)
                   End If
                Else
                   WidthX=(xCol(ix+1)-xCol(ix-1))
                End If
             Else If(ix.Eq.1) Then
                WidthX=(xCol(ix+1)-xCol(1))
                If(continu) WidthX=(xCol(2)-xCol(nx)+nx*dxGrid)
             Else If(ix.Eq.nx) Then
                WidthX=(xCol(nx)-xCol(nx-1))
                If(continu) WidthX=(xCol(1)-xCol(nx-1)+nx*dxGrid)
             End If
             If(iy.Ne.1.And.iy.Ne.ny) Then
                WidthY=(yCol(iy+1)-yCol(iy-1))
             Else If(iy.Eq.1) Then
                WidthY=(yCol(iy+1)-yCol(1))
                If(continu) WidthY=(yCol(2)-yCol(ny)+ny*dyGrid)
             Else If(iy.Eq.ny) Then
                WidthY=(yCol(ny)-yCol(ny-1))
                If(continu) WidthY=(yCol(1)-yCol(ny-1)+ny*dyGrid)
             End If
          End If
          iCord=ix+iy
          iRest=Mod(iCord,2)
          iDiv=3
          If(iRest.Eq.1) iDiv=6 
          Width(i)=WidthX*WidthY/iDiv !top surface
          iCord=ix+iy+nz
          iRest=Mod(iCord,2)
          iDiv=3
          If(iRest.Eq.0) iDiv=6
          Width(iB)= WidthX*WidthY/iDiv !bottom surface
       End Do
    End Do
    Return
  End Subroutine CalcWidthSWMS_new
!****************************************************************
!> calculates the soil conductivity - the used model is  defined in the input files
!K is only important for setmatRhizo. Does not change anything for setmat  
  Subroutine SetMat(K,Explic)
    Use Typedef
    Use GridData, Only : nPt,Axy,Bxy,Dxy,Exy
    Use SolData
    Use WatFun
    Use RhizoData
    Implicit None

    Real(dp):: ci,cpi,Ti!,t0,t1
    Real(dp):: sl,s2,him,hi2,hi1
    Integer(ap):: it,i,K,m
    Logical Explic

    If (lRhizo) Then
       Call SetMatRhizo(K)
    Else
       !$OMP PARALLEL DO private(i,M,hi1,hi2,hiM,iT,Sl,S2,ci,Cpi,Ti)
       Do i=1,nPt
! Calculate nodal conductivity values:
          If(K.Eq.1) conO(i)=con(i)
          M=MatNum(i)
          hi1=hTemp(i)/Axy(i) !amin1(hSat(M),hTemp(i)/Axz(i))
          hi2=hNew(i)/Axy(i)!amin1(hSat(M), hNew(i)/Axz(i))
          if(Explic) hi2=hi1
          hiM=0.1_dp*hi1+0.9_dp*hi2
          If (soiltab) Then
             ci=FKP_soiltab(hiM,M)
          Elseif (lTab.And.(hiM.Ge.hTab(nTab)).And.(hiM.Le.hTab(1))) Then
             iT=Int((Log10(-hiM)-alh1)/dlh)+1
             Sl=(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))
             ci=ConTab(iT,M)+Sl*(hiM-hTab(iT))
          Else
             ci=FKP(hiM,par(:,M),i)
          Endif
          con(i)=ci*Bxy(i)
          If(K.Eq.0) conO(i)=con(i)

! Calculate nodal capacity values at k+1, except if explicit:
          hi1=hOld(i)/Axy(i)
          hi2=hNew(i)/Axy(i)
          If(Explic) hi2=hi1
          If (soiltab) Then
             Cpi=FCP_soiltab(hi2,M)
             Ti=Fth_soiltab(hi2,M)
          Elseif (lTab.And.(hi2.Ge.hTab(nTab)).And.(hi2.Le.hTab(1))) Then
             iT=Int((Log10(-hi2)-alh1)/dlh)+1
             Sl =(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))
             S2 =(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))
             Cpi=CapTab(iT,M)+Sl*(hi2-hTab(iT))
             Ti=TheTab(iT,M)+S2*(hi2-hTab(iT))
          Else
             Cpi=FCP(hi2,par(:,M))
             Ti=Fth(hi2,par(:,M))
          Endif
          Cap(i)=Cpi*Dxy(i)/Axy(i)
          If (soiltab) Then
             theta(i)=TheTab(nTab,M)*Exy(i)+(Ti-TheTab(nTab,M))*Dxy(i)
          Else
             theta(i)=par(2,M)*Exy(i)+(Ti-par(2,M))*Dxy(i)
          Endif
       Enddo
       !$OMP END PARALLEL DO
    Endif

    Return
  End Subroutine SetMat
!************************************************************************************
  !> SetMat for the Rhizosphere
  !> Calculate soil conductivity and capacity
  Subroutine SetMatRhizo(K)
    Use Typedef
    Use GridData, Only : nPt
    Use SolDAta
    Use WatFun
    Use RhizoData
    Implicit None

    Real(dp)    :: him,hi2,hi1
    Integer(ap) :: i,K,M

    M = MatNum(1)
    Do i=1,nPt
       !! In this section the conductivity is calculated
       If(K .Eq. 1) conO(i)=con(i)
       hi1 = hTemp(i)
       hi2 = hNew(i)
       hiM = 0.1_dp*hi1 + 0.9_dp*hi2
       theta(i) = Fth(hiM,par(:,M))
       ! Calculate the initial thetaTot water contet
       If (K.Eq. 0) Call IniRhizoTheta(hiM,i)
       ! If rhizoStatic, calculate the ThetaTot
       If (K .Eq. 1 .And. RhizoModel .Eq. 2) Call IniRhizoTheta(hiM,i)
       con(i) = FKP(hiM,par(:,M),i)
       If (K .Eq. 0) conO(i)=con(i)

       !! In this section the capacity is calculated
       !! For the capacity hiM is calculated differently. 
       hi1=hOld(i)
       hiM = (hi2+hi1)/2
       ! Capacity
       Cap(i) = FCP(hiM,par(:,M))
       ! The water content is recalculated based on the hiM
       theta(i) = Fth(hiM,par(:,M))
       If (K .Eq. 0) Call IniRhizoTheta(hiM,i)
       If (K .Eq. 1 .And. RhizoModel .Ne. 1) Call IniRhizoTheta(hiM,i)
    Enddo
    Return
  End Subroutine SetMatRhizo


  !************************************************************************************
  !> calculate the geometry, shape_functions etc needed for the finite element method   
  !> takes care of splitting the cube into 5 tetrahedrals 
  Subroutine CalcGeom
    Use Typedef
    Use GridData
    Use CumData
    Use DomData
    Use MatData
    Implicit None

    Integer(ap) :: iE,iSE,i,j,k,l,ii,jj
    Integer(ap) :: ex,ey
    Real(dp) :: xGridi,xGridj,xGridk,xGridl,yGridi,yGridj,yGridk,yGridl
    Real(dp) :: Cxx,Czz,Cyy,Cxy,Cxz,Cyz,VE,ai(1:4)

    ! Loop on elements:
    ex=0
    ey=1
    Do iE=1,nElm  ! loop over all cubes
       Cxx=ConAxx(iE)
       Cyy=ConAyy(iE)
       Czz=ConAzz(iE)
       Cxy=ConAxy(iE)
       Cxz=ConAxz(iE)
       Cyz=ConAyz(iE)
       ex=ex+1!Count the position of the element in the grid (Couvreur dec 2009)
       If (ex.Gt.nex) Then
          ex=1
          ey=ey+1
          If (ey.Gt.ney) ey=1
       Endif
       ! Loop on subelements: 5 tetrahedrals per soil cube
       Do iSE=1,5     
          i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          !if master node then coefficient calculated as normal; secondly signs are needed for orientation
          If (((ex.Eq.nex).Or.(ey.Eq.ney)).And.continu) Then!"Bridge sub-element?" (Couvreur dec 2009)
             xGridi=xGrid(i)
             xGridj=xGrid(j)
             xGridk=xGrid(k)
             xGridl=xGrid(l)
             yGridi=yGrid(i)
             yGridj=yGrid(j)
             yGridk=yGrid(k)
             yGridl=yGrid(l)
             If (ex.Eq.nex) Then
                Select Case(subN(iE))
                Case(1) 
                   Select Case(iSE)
                   Case (1,3) 
                      xGridj = xGridj + nex*dxgrid
                   Case (2,4)
                      xGridj = xGridj + nex*dxgrid
                      xGridk = xGridk + nex*dxgrid
                      xGridl = xGridl + nex*dxgrid                           
                   Case (5)
                      xGridj = xGridj + nex*dxgrid
                      xGridl = xGridl + nex*dxgrid
                   End Select
                Case (2)
                   Select Case(iSE)
                   Case (1,3) 
                      xGridj = xGridj + nex*dxgrid
                   Case (2) 
                      xGridj = xGridj + nex*dxgrid
                      xGridk = xGridk + nex*dxgrid
                      xGridl = xGridl + nex*dxgrid                           
                   Case (4)
                      xGridk = xGridk + nex*dxgrid
                      xGridj = xGridj + nex*dxgrid  
                   Case (5) 
                      xGridj = xGridj + nex*dxgrid
                      xGridi = xGridi + nex*dxgrid
                      xGridl = xGridl + nex*dxgrid  
                   End Select
                End Select
             End If

             If (ey.Eq.ney) Then
                Select Case (subN(iE))
                Case (1) 
                   Select Case (iSE)
                   Case (1) 
                      yGridj = yGridj + ney*dygrid
                      yGridk = yGridk + ney*dygrid
                      yGridl = yGridl + ney*dygrid
                   Case (2,3) 
                      yGridk = yGridk + ney*dygrid
                   Case (4) 
                      yGridi = yGridi + ney*dygrid
                      yGridk = yGridk + ney*dygrid
                      yGridl = yGridl + ney*dygrid
                   Case (5) 
                      yGridk = yGridk + ney*dygrid
                      yGridl = yGridl + ney*dygrid
                   End Select
                Case(2) 
                   Select Case(iSE)
                   Case(1,2) 
                      yGridk = yGridk + ney*dygrid  
                   Case (3,5)
                      yGridk = yGridk + ney*dygrid
                      yGridj = yGridj + ney*dygrid
                      yGridl = yGridl + ney*dygrid
                   Case (4)
                      yGridk = yGridk + ney*dygrid
                      yGridl = yGridl + ney*dygrid
                   End Select
                End Select
             End If

             bi(1,iSE,iE)=-(yGridk-yGridj)*(zGrid(l)-zGrid(j))+(yGridl-yGridj)*(zGrid(k)-zGrid(j))!bi, ci and di are now globals 
             bi(2,iSE,iE)=+(yGridl-yGridk)*(zGrid(i)-zGrid(k))-(yGridi-yGridk)*(zGrid(l)-zGrid(k))
             bi(3,iSE,iE)=-(yGridi-yGridl)*(zGrid(j)-zGrid(l))+(yGridj-yGridl)*(zGrid(i)-zGrid(l))
             bi(4,iSE,iE)=+(yGridj-yGridi)*(zGrid(k)-zGrid(i))-(yGridk-yGridi)*(zGrid(j)-zGrid(i))
             ci(1,iSE,iE)=+(xGridk-xGridj)*(zGrid(l)-zGrid(j))-(xGridl-xGridj)*(zGrid(k)-zGrid(j))
             ci(2,iSE,iE)=-(xGridl-xGridk)*(zGrid(i)-zGrid(k))+(xGridi-xGridk)*(zGrid(l)-zGrid(k))
             ci(3,iSE,iE)=+(xGridi-xGridl)*(zGrid(j)-zGrid(l))-(xGridj-xGridl)*(zGrid(i)-zGrid(l))
             ci(4,iSE,iE)=-(xGridj-xGridi)*(zGrid(k)-zGrid(i))+(xGridk-xGridi)*(zGrid(j)-zGrid(i))
             di(1,iSE,iE)=-(xGridk-xGridj)*(yGridl-yGridj)+(xGridl-xGridj)*(yGridk-yGridj)
             di(2,iSE,iE)=+(xGridl-xGridk)*(yGridi-yGridk)-(xGridi-xGridk)*(yGridl-yGridk)
             di(3,iSE,iE)=-(xGridi-xGridl)*(yGridj-yGridl)+(xGridj-xGridl)*(yGridi-yGridl)
             di(4,iSE,iE)=+(xGridj-xGridi)*(yGridk-yGridi)-(xGridk-xGridi)*(yGridj-yGridi)
             ! coefficient a-> shape_function = a + bx +cy +dz
             ai(1)=xGridj*yGridk*zGrid(l) + xGridk*yGridl*zGrid(j) + xGridl*yGridj*zGrid(k) - &
                  xGridl*yGridk*zGrid(j) - xGridj*yGridl*zGrid(k) - xGridk*yGridj*zGrid(l)
             ai(2)=xGridi*yGridl*zGrid(k) + xGridk*yGridi*zGrid(l) + xGridl*yGridk*zGrid(i) - &
                  xGridi*yGridk*zGrid(l) - xGridk*yGridl*zGrid(i) - xGridl*yGridi*zGrid(k)
             ai(3)=xGridi*yGridj*zGrid(l) + xGridj*yGridl*zGrid(i) + xGridl*yGridi*zGrid(j) - &
                  xGridi*yGridl*zGrid(j) - xGridj*yGridi*zGrid(l) - xGridl*yGridj*zGrid(i)
             ai(4)=- xGridi*yGridj*zGrid(k) - xGridj*yGridk*zGrid(i) - xGridk*yGridi*zGrid(j) + &
                  xGridi*yGridk*zGrid(j) + xGridj*yGridi*zGrid(k) + xGridk*yGridj*zGrid(i)

             Deter(iSE,iE)=(xGridl-xGridi)*bi(4,iSE,iE)+(yGridl-yGridi)*ci(4,iSE,iE)+(zGrid(l)-zGrid(i))*di(4,iSE,iE)!Deter, Ax, Ay, Az and B1fact are now globals (Couvreur mar 2010)
             Ax(:,iSE,iE)=Cxx*bi(:,iSE,iE)+Cxy*ci(:,iSE,iE)+Cxz*di(:,iSE,iE)
             Ay(:,iSE,iE)=Cxy*bi(:,iSE,iE)+Cyy*ci(:,iSE,iE)+Cyz*di(:,iSE,iE)
             Az(:,iSE,iE)=Cxz*bi(:,iSE,iE)+Cyz*ci(:,iSE,iE)+Czz*di(:,iSE,iE)
          Else
             bi(1,iSE,iE)=-(yGrid(k)-yGrid(j))*(zGrid(l)-zGrid(j))+(yGrid(l)-yGrid(j))*(zGrid(k)-zGrid(j))
             bi(2,iSE,iE)=+(yGrid(l)-yGrid(k))*(zGrid(i)-zGrid(k))-(yGrid(i)-yGrid(k))*(zGrid(l)-zGrid(k))
             bi(3,iSE,iE)=-(yGrid(i)-yGrid(l))*(zGrid(j)-zGrid(l))+(yGrid(j)-yGrid(l))*(zGrid(i)-zGrid(l))
             bi(4,iSE,iE)=+(yGrid(j)-yGrid(i))*(zGrid(k)-zGrid(i))-(yGrid(k)-yGrid(i))*(zGrid(j)-zGrid(i))
             ci(1,iSE,iE)=+(xGrid(k)-xGrid(j))*(zGrid(l)-zGrid(j))-(xGrid(l)-xGrid(j))*(zGrid(k)-zGrid(j))
             ci(2,iSE,iE)=-(xGrid(l)-xGrid(k))*(zGrid(i)-zGrid(k))+(xGrid(i)-xGrid(k))*(zGrid(l)-zGrid(k))
             ci(3,iSE,iE)=+(xGrid(i)-xGrid(l))*(zGrid(j)-zGrid(l))-(xGrid(j)-xGrid(l))*(zGrid(i)-zGrid(l))
             ci(4,iSE,iE)=-(xGrid(j)-xGrid(i))*(zGrid(k)-zGrid(i))+(xGrid(k)-xGrid(i))*(zGrid(j)-zGrid(i))
             di(1,iSE,iE)=-(xGrid(k)-xGrid(j))*(yGrid(l)-yGrid(j))+(xGrid(l)-xGrid(j))*(yGrid(k)-yGrid(j))
             di(2,iSE,iE)=+(xGrid(l)-xGrid(k))*(yGrid(i)-yGrid(k))-(xGrid(i)-xGrid(k))*(yGrid(l)-yGrid(k))
             di(3,iSE,iE)=-(xGrid(i)-xGrid(l))*(yGrid(j)-yGrid(l))+(xGrid(j)-xGrid(l))*(yGrid(i)-yGrid(l))
             di(4,iSE,iE)=+(xGrid(j)-xGrid(i))*(yGrid(k)-yGrid(i))-(xGrid(k)-xGrid(i))*(yGrid(j)-yGrid(i))
             ! coefficient a-> shape_function = a + bx +cy +dz
             ai(1)=xGrid(j)*yGrid(k)*zGrid(l) + xGrid(k)*yGrid(l)*zGrid(j) + xGrid(l)*yGrid(j)*zGrid(k) - &
                  xGrid(l)*yGrid(k)*zGrid(j) - xGrid(j)*yGrid(l)*zGrid(k) - xGrid(k)*yGrid(j)*zGrid(l)
             ai(2)=xGrid(i)*yGrid(l)*zGrid(k) + xGrid(k)*yGrid(i)*zGrid(l) + xGrid(l)*yGrid(k)*zGrid(i) - &
                  xGrid(i)*yGrid(k)*zGrid(l) - xGrid(k)*yGrid(l)*zGrid(i) - xGrid(l)*yGrid(i)*zGrid(k)
             ai(3)=xGrid(i)*yGrid(j)*zGrid(l) + xGrid(j)*yGrid(l)*zGrid(i) + xGrid(l)*yGrid(i)*zGrid(j) - &
                  xGrid(i)*yGrid(l)*zGrid(j) - xGrid(j)*yGrid(i)*zGrid(l) - xGrid(l)*yGrid(j)*zGrid(i)
             ai(4)=- xGrid(i)*yGrid(j)*zGrid(k) - xGrid(j)*yGrid(k)*zGrid(i) - xGrid(k)*yGrid(i)*zGrid(j) + &
                  xGrid(i)*yGrid(k)*zGrid(j) + xGrid(j)*yGrid(i)*zGrid(k) + xGrid(k)*yGrid(j)*zGrid(i)

             Deter(iSE,iE)=(xGrid(l)-xGrid(i))*bi(4,iSE,iE)+(yGrid(l)-yGrid(i))*ci(4,iSE,iE)+(zGrid(l)-zGrid(i))*di(4,iSE,iE)
             Ax(:,iSE,iE)=Cxx*bi(:,iSE,iE)+Cxy*ci(:,iSE,iE)+Cxz*di(:,iSE,iE)
             Ay(:,iSE,iE)=Cxy*bi(:,iSE,iE)+Cyy*ci(:,iSE,iE)+Cyz*di(:,iSE,iE)
             Az(:,iSE,iE)=Cxz*bi(:,iSE,iE)+Cyz*ci(:,iSE,iE)+Czz*di(:,iSE,iE)
          Endif
          VE=Abs(Deter(iSE,iE))/6.
          If (Deter(iSE,iE).Eq.0 .Or. VE.Eq.0) Then
             Write(*,*)'Ve=0',VE,'deter=',Deter(iSE,iE),'element nr.',iE
             Print*,'coo i',xgrid(i),ygrid(i),zgrid(i)
             Print*,'coo j',xgrid(j),ygrid(j),zgrid(j)
             Print*,'coo k',xgrid(k),ygrid(k),zgrid(k)
             Print*,'coo l',xgrid(l),ygrid(l),zgrid(l)
             Read(*,*)
          Endif

          Do ii=1,4
             B1fact(ii,iSE,iE)=Cxz*bi(ii,iSE,iE)+Cyz*ci(ii,iSE,iE)+Czz*di(ii,iSE,iE)
             Do jj=1,4
                E(jj,ii,iSE,iE) = Cxx*bi(ii,iSE,iE)*bi(jj,iSE,iE)+Cyy*ci(ii,iSE,iE)*ci(jj,iSE,iE)+Czz*di(ii,iSE,iE)*di(jj,iSE,iE) &
                     +Cxy*(bi(ii,iSE,iE)*ci(jj,iSE,iE)+bi(jj,iSE,iE)*ci(ii,iSE,iE))+Cxz*(bi(ii,iSE,iE)*di(jj,iSE,iE)+bi(jj,iSE,iE)*di(ii,iSE,iE)) &
                     +Cyz*(ci(ii,iSE,iE)*di(jj,iSE,iE)+ci(jj,iSE,iE)*di(ii,iSE,iE))
             End Do
          End Do
       End Do
    End Do
  End Subroutine CalcGeom
!****************************************************************************
  !> Assembles the matrix and the right-hand side of the linear equation system
  !> based on the standard galerkin finite element method
  !> use of sparse format 
  Subroutine Reset(dt,IAD,IADN,IADD)
    Use Typedef
    Use ParamData, Only: maxbnd
    Use GridData, Only: nPt,nElm,rootSkOld,subN,elmnod,iL,deter,sink_cube,B1fact, iadd_temp,e,width
    Use CumData
    Use SolData
    Use DomData
    Use MatData
    Use tmctrl, Only: tlevel_soil
    Use RhizoData, Only: Rnorm, tauTht, hEqRhizo, RhizoModel, lRhizo
    Use Orthofem, Only: Find
    Implicit None

    Integer(ap) :: ind,ind1
    Integer(ap) :: iSE,i,j,k,l,iE,n,ii,jj,kk,iG,jG,kkk
    Integer(ap) :: IAD(maxbnd,nPt),IADN(nPt),IADD(nPt)
    Real(dp) :: cone,cape  
    Real(dp) :: fmul,bmul
    Real(dp) :: amul,QN
    Real(dp) :: dt
    Real(dp) :: SinkE
    Real(dp) :: VE
    Real(dp), Allocatable,Dimension (:) :: B1,F,DS

    Allocate (B1(nPt),F(nPt),DS(nPt))

! Initialisation
    a_dparse(1:irow(npt+1)-1) = 0.0d0
    B(1:npt)=0.0
    B1=0.0
    DS=0.0
    F=0.0
    RootSkold=0.0

!save matrix index for sparse format in the first step -> tlevel_soil = false
    !> \param IAD adjacency matrix (nodal connections)
    !> \param IADN number of adjacent nodes in IAD (self-inclusive)
    !> \param IADD position of diagonal in adjacency matrix
    k= 0
    If(.Not.tlevel_soil) Then
! Loop on cubic elements:
       Do  iE=1,nElm
! Loop on subelements: 5 tetrahedrals per soil cube
          Do  iSE=1,5
! Loop on nodes: 4 nodes per tetraheders  
             Do  ii=1,4
                iG=elmnod(iL(ii,iSE,subN(iE)),iE)
                Do  jj=1,4
                   jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                   Call FIND (iG,jG,kk,nPt,maxbnd,IAD,IADN)
                   k = k +1
                   iadd_temp(k) = kk
                Enddo
             Enddo
          Enddo
       Enddo
       tlevel_soil = .True.
    Endif
    
    kkk = 0 !index for matrix index: iadd_temp
! Loop on cubic elements:
    Do iE=1,nElm
! Loop on subelements: 5 tetrahedrals per soil cube
       Do iSE=1,5 
          i=elmnod(iL(1,iSE,subN(iE)),iE)! i,j,k,l are corner nodes of the tetrahedal element
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          CapE=(Cap(i)+Cap(j)+Cap(k)+Cap(l))/4.
          ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4.
          VE=Abs(Deter(iSE,iE))/6.!Deter and B1fact are now globals 
          AMul=ConE/VE/36.![/L²/T]
          BMul=ConE/6.![L/T]
          FMul=VE/20.![L³]
! DS-term (sink)
          SinkE=FMul*5*sink_cube(iE) ![L³/T]
          DS(i)=DS(i)+SinkE ![L³/T]
          DS(j)=DS(j)+SinkE
          DS(k)=DS(k)+SinkE
          DS(l)=DS(l)+SinkE
!FC-term
          F(i)=F(i)+FMul*(4*CapE+Cap(i)) ![L²]
          F(j)=F(j)+FMul*(4*CapE+Cap(j))
          F(k)=F(k)+FMul*(4*CapE+Cap(k)) 
          F(l)=F(l)+FMul*(4*CapE+Cap(l)) 

          Do ii=1,4
             iG=elmnod(iL(ii,iSE,subN(iE)),iE)
             B1(iG)=B1(iG)+BMul*B1fact(ii,iSE,iE)![L³/T](B1fact [L²])
             Do jj=1,4
                jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                kkk = kkk+1
!A_dparse (CSR format)
                A_dparse(IROW(iG)-1+iadd_temp(kkk))=A_dparse(IROW(iG)-1+iadd_temp(kkk))+AMul*E(jj,ii,iSE,iE)              
             End Do
          End Do
       End Do
    End Do

! Determine Boundary fluxes, loop over all nodes
    Do N=1,nPt
! PH boundary nodes
       If (Kode(N).Gt.0) Then
          QN=B1(N)+DS(N) ! previously: +F(N)*(theta(N)-theta_old(N))/dt ![L³/T]
          ind=IROW(N)
          ind1=IROW(N+1)-1
          Do j=ind,ind1
             QN=QN+A_dparse(j)*hNew(JCOL(j))
             !if a row is part of a bound.cond. then get the values in A_dparse and multiply by the corresponding column values in hNew
          Enddo
          Q(N)=QN
! Free drainage       
       Elseif (Kode(N).Eq.-2) Then
          Q(N)=-con(N)*Width(N)
       Endif

! Construction of effective matrix:
       If(.Not. lRhizo .Or. RhizoModel .Eq. 1 .Or. RhizoModel .Eq. 2) Then
! with equilibrium
          A_dparse(IROW(N)-1+IADD(N)) = A_dparse(IROW(N)-1+IADD(N)) +  F(N)/dt! use row index - 1 to get the number of previous entries, add up the entry number of the diagonal value to find the location of the diagonal value in A_dparse
! with non equilibrium 
       Else
          A_dparse(IROW(N)-1+IADD(N)) = A_dparse(IROW(N)-1+IADD(N)) + Rnorm(N)*F(N)/dt + (1. - Rnorm(N))*FMul/tauTht(N) 
       Endif
    End Do

! Complete construction of RHS vector(no i-loop neeeded):
    If(.Not. lRhizo .Or. RhizoModel .Eq. 1 .Or. RhizoModel .Eq. 2) Then
       B=Q-B1+hOld*F/dt-DS
    Else
       B = Q - B1 + Rnorm*hOld*F/dt-DS + (1. - Rnorm)*hEqRhizo*FMul/tauTht
    Endif

    Deallocate (B1,F,DS)
  End Subroutine Reset
!****************************************************************************
  !> assambles the matrix and the right-hand side of the linear equation system
  !> based on the standard galerkin finite element method
  !> used the mass lumping algorith from Celia et al 1990 (similar to SWMS_3D)
  Subroutine ResetCelia(dt,IAD,IADN,IADD)
    Use Typedef
    Use ParamData, Only: maxbnd,NoGrav
    Use GridData, Only: nPt,nElm,rootSkOld,subN,elmnod,iL,deter,sink_cube,b1fact, iadd_temp,e,width
    Use CumData
    Use SolData
    Use DomData
    Use MatData
    Use WatFun
    Use tmctrl, Only: tlevel_soil
    Use Orthofem, Only: Find
    Implicit None

    Integer(ap) :: ind,ind1
    Integer(ap) :: iSE,i,j,k,l,iE,n,ii,jj,kk,iG,jG,kkk
    Integer(ap) :: IAD(maxbnd,nPt),IADN(nPt),IADD(nPt)
    Real(dp) :: cone,cape
    Real(dp) :: fmul,bmul
    Real(dp) :: amul,QN,VE
    Real(dp) :: dt
    Real(dp) :: SinkE
    Real(dp), Allocatable,Dimension (:) ::B1, F ,DS 

    Allocate (B1(nPt),F(nPt),DS(nPt))

! Initialisation
    a_dparse(1:irow(npt+1)-1) = 0.0d0
    B(1:npt)=0.0
    B1=0.0
    DS=0.0
    F=0.0
    RootSkold=0.0

!save matrix index for sparse format in the first step -> tlevel_soil = false    !> \param dt time step size
    !> \param IAD adjacency matrix (nodal connections)
    !> \param IADN number of adjacent nodes in IAD (self-inclusive)
    !> \param IADD position of diagonal in adjacency matrix
    k= 0
    If(.Not.tlevel_soil) Then
! Loop on cubic elements (initial loop):
       Do  iE=1,nElm
! Loop on subelements: 5 tetrahedrals per soil cube
          Do  iSE=1,5
! Loop on nodes: 4 nodes per tetraheders
             Do  ii=1,4
                iG=elmnod(iL(ii,iSE,subN(iE)),iE)
                Do  jj=1,4
                   jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                   Call FIND (iG,jG,kk,nPt,maxbnd,IAD,IADN)
                   k = k +1
                   iadd_temp(k) = kk
                Enddo
             Enddo
          Enddo
       Enddo
       tlevel_soil = .True.
    Endif

    kkk=0
! Loop on cubic elements:
    Do iE=1,nElm
! Loop on subelements: 5 tetrahedrals per soil cube
       Do iSE=1,5
          i=elmnod(iL(1,iSE,subN(iE)),iE) !i,j,k,l are corner nodes of the tetraedal element
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          CapE=(Cap(i)+Cap(j)+Cap(k)+Cap(l))/4
          ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4
          VE=Abs(Deter(iSE,iE))/6.
          AMul=ConE/VE/36.
          If (NoGrav) Then
             BMul=0.
          Else
             BMul=ConE/6.
          Endif
          FMul=VE/20.
! D-term 
          SinkE=FMul*5*sink_cube(iE)![L³/T]
          DS(i)=DS(i)+SinkE![L³/T]
          DS(j)=DS(j)+SinkE
          DS(k)=DS(k)+SinkE
          DS(l)=DS(l)+SinkE
 
!flow leaving soil is average sink times their volume over all the elements
          RootSkold=RootSkold+VE*SinkE
          Do ii=1,4
             iG=elmnod(iL(ii,iSE,subN(iE)),iE)
! F-term
             F(iG)=F(iG)+FMul*5
             B1(iG)=B1(iG)+BMul*B1fact(ii,iSE,iE)
             Do jj=1,4
                jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                kkk = kkk+1
! A_dparse (CSR format)
                A_dparse(IROW(iG)-1+iadd_temp(kkk))=A_dparse(IROW(iG)-1+iadd_temp(kkk))+AMul*E(jj,ii,iSE,iE)
             End Do
          End Do
       End Do
    End Do

! Determine Boundary fluxes :
    Do N=1,nPt
! PH boundary node
       If (Kode(N).Gt.0) Then !hbc
          QN=B1(N)+DS(N)+F(N)*(theta(N)-theta_old(N))/dt
          ind=IROW(N)
          ind1=IROW(N+1)-1
          Do j=ind,ind1
             QN=QN+A_dparse(j)*hNew(JCOL(j))
             !if a row is part of a bound.cond. then get the values in A_dparse and multiply by the corresponding column values in hNew
          Enddo
          Q(N)=QN
! Free drainage
       Elseif (Kode(N).Eq.-2) Then
          Q(N)=-con(N)*Width(N)
       Endif
! Form effective matrix:
    A_dparse(IROW(N)-1+IADD(N)) = A_dparse(IROW(N)-1+IADD(N)) + F(N)*Cap(N)/dt! use row index - 1 to get the number of previous entries, add up the entry number of the diagonal value to find the location of the diagonal value in A_dparse!
    End do

! Complete construction of RHS vector
    B=F*Cap*hNew/dt-F*(theta-theta_old)/dt+ Q-B1-DS

    Deallocate (B1,F,DS)
  End Subroutine ResetCelia
  !***********************************************************************
  !> adding Dirichlet boundary condition to the sparse matrix 
  Subroutine Dirich(IADD)
    Use Typedef
    Use GridData, Only:nPt
    Use MatData, Only:A_dparse,B,IROW
    Use SolData, Only: Kode, hnew
    Implicit None

    Integer(ap) :: N,IADD(nPt),ind
    !> \param IADD position of diagonal in adjacency matrix
    !usage of compressed sparse row format
    ! IA(1:nPt+1) represent the rows + 1. The values in IA correspond to the positions in JA
    ! (column indices) and AA (nonzero values). The indicators ind and ind1 can be defined
    ! which show for one row the indices for the arrays JA and AA.
    !
    ! example (from Youcef Saad: numerical methods for large eigenvalue problems 1992)
    !
    ! AA = 1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
    ! JA = 1 4 1 2 4 1 3 4 5 3 4 5
    ! IA = 1 3 6 10 12 13
    !
    ! A = |   1. 0. 0.  2.  0.   |
    !     |   3. 4. 0.  5.  0.   |
    !     |   6. 0. 7.  8.  9.   |
    !     |   0. 0. 10. 11. 0.   |
    !     |   0. 0. 0.  0.  12.  |
    Do N=1,nPt
       If (Kode(N).Ge.1) Then
       !ind is row index; to find where diagonal position in A_dparse is
       !located take row index -1 (gives number of previous entries in
       !A_dparse and add up the number which corresponds to the diagonal value
       !, i.e., IADD(N)
       ind=IROW(N)
       A_dparse(ind-1+IADD(N))=10.d30
       B(N)=10.d30*hNew(N)
       End if
    End Do
  End Subroutine Dirich
!*****************************************************************************
!> time step control: calculates the next time step size for  
!> water, solute and root growth
  Subroutine TmCont(iter,t,dt,dtOpt,tCallR,tFEMRoo,dtMaxC,tcBCr,tProf,tProbe)
    Use typedef
    Use tmctrl, Only: dtMax,dtMin,facdec,facinc,tmax
    Implicit None

    Integer(ap) :: iter
    Real(dp), Intent(out) :: dt
    Real(dp), Intent(in) :: t,dtMaxC,tCallR,tFEMRoo,tcBCr,tProf,tProbe
    Real(dp), Intent(inout) :: dtOpt
    Real(dp)::tfix
    !> \param iter
    !> \param t current simulation time
    !> \param dt current time step size
    !> \param dtOpt time step when there has been a decrease in Water.
    !> \param tCallR
    !> \param tFEMRoo output time for the next FEM/SOIL and ROOT output
    !> \param dtMaxC
    !> \param tcBCr
    !> \param tProf
    !> \param tProbe

    dtMax=Min(dtMax,dtMaxC)
    tFix=Min(tcallr,tFEMRoo,tMax,tcBCr,tProf,tProbe)
!    Write (*,'(6(1X,1pe9.3),1X,I4)')  tFix-t, dtMax, dtMin, dtOpt, tcallr,tFEMRoo,iter
    If (iter.Le.4.And.(tFix-t).Ge.FacInc*dtOpt)  dtOpt=Min(dtMax,FacInc*dtOpt)!3
    If (iter.Ge.7)  dtOpt=Max(dtMin,FacDec*dtOpt) !7
    dt=Min(dtOpt,tFix-t)
    dt=Min((tFix-t)/Anint((tFix-t)/dt),dtMax)
    If(tFix-t.Ne.dt.And.dt.Gt.(tFix-t)/2._dp) dt=(tFix-t)/2._dp
 ! Write (*,'(2(1X,1pe9.3))')  dtOpt, dt

  End Subroutine TmCont
  !****************************************************************
  Subroutine WatInf(dt)
    Use Typedef
    Use GridData, Only: dxGrid,dyGrid,dzGrid,nElm,nPt,RootSk,sink_cube
    Use CumData
    Use RootData, Only: lFed, lCou
    Use SolData, Only: Kode
    Implicit None

    Integer(ap)::i,j
    Real(dp):: vMean(3),dt

    vMean=0.0_dp

    Do i=1,nPt
       j=Abs(Kode(i))
       wCumA=wCumA+Abs(Q(i))*dt
       If (j.Ne.0) Then
          vMean(j)=vMean(j)-Q(i)
       Endif
    End Do
    If (lCou.Or.lFed) Then
       RootSk=0.0_dp
       Do i=1,nElm
          RootSk=RootSk+sink_cube(i)*dxGrid*dyGrid*dzGrid
       Enddo
    Endif
    wCumA=wCumA+(RootSk*dt)!abs(RootSk*dt)
    CumRt=CumRt+RootSk*dt
    VolSink = RootSk*dt
    wCumT=CumRt
    CumQ(1)=CumQ(1)+(vMean(1)+vMean(3))*dt
    CumQ(2)=CumQ(2)+vMean(2)*dt
    wCumT=wCumT+CumQ(1)+CumQ(2)
  End Subroutine WatInf
  !****************************************************************
  !< SUBROUTINE TO CALCULATE THE INITIAL WATER CONTENT IN THE RHIZOSPEHRE
  Subroutine IniRhizoTheta(h0,nodeID)
    Use Typedef
    Use RhizoData, Only: bulkPara, StaticRhizoPara, thetaTot, Rnorm
    Use solData,   Only: theta
    Implicit None

    Integer(ap)  :: nodeID
    Real(dp) :: h0, thetaEqIni=0.
    Real(dp) :: lambda_s, hcr_s, thtR, thtS

    thtR = bulkPara(1); thtS = bulkPara(2)
    lambda_s = StaticRhizoPara(1)
    hcr_s    = StaticRhizoPara(2)

    If (Rnorm(nodeID) .Lt. 1.) Then
       If (h0 .Lt. hcr_s) Then
          thetaEqIni = thtR +(thtS-thtR)*(hcr_s/h0)**lambda_s
       Else
          thetaEqIni = thtS
       Endif
       thetaTot(nodeID) = Rnorm(nodeID)*theta(nodeID) + (1. - Rnorm(nodeID))*thetaEqIni
    Else
       thetaTot(nodeID) = theta(nodeID)
    Endif
  End Subroutine IniRhizoTheta
  !***************************************************************
  !< SUBROUTINE TO CALCULATE heq (equilibrium pressure head in the rhizopshere
  Subroutine hEq(thtIn, hIn, i)
    ! This subroutine calculate hEq base on theta and location (i.e. cTot)
    Use Typedef
    Use RhizoData, Only: bulkPara, RhizoPara, hEqRhizo, cTot_r
    Implicit None

    Integer(ap)  :: i
    Real(dp) :: thtR, thtS, hcr, lambda, cw, rhow, rhob, beta, omega
    Real(dp) :: thtIn, hIn, h1

    thtR = bulkPara(1); thtS = bulkPara(2); lambda = bulkPara(3); hcr = bulkPara(4)
    omega = RhizoPara(1); beta = RhizoPara(2); rhob = RhizoPara(8); rhow = RhizoPara(9)

    If (thtIn .Lt. thtS) Then
       h1 = hcr/((thtIn-thtR)/(thtS-thtR))**(1./lambda)
    Else
       h1 = hIn
    Endif
    cw = cTot_r(i)*rhob/(rhow*thtIn)
    hEqRhizo(i) = h1 - omega*cw**beta
  End Subroutine hEq
  !**************************************************************
  !< SUBROUTINE TO CALCULATE thetaNonEquilibrium
  Subroutine thtNonEq(hNew, dt, i)
    Use Typedef
    Use SolData, Only: theta
    Use RhizoData, Only: bulkPara, RhizoPara, thetaNonEq, thetaNonEqOld, thetaTot, hEqRhizo, thetaTot, Rnorm, tauTht
    Implicit None

    Integer(ap)  :: i, maxIter, iter
    Real(dp) :: hNew, dt
    Real(dp) :: tolerance, diffTht
    Real(dp) :: gamm, tau0, thtS

    ! Initiate parameters
    tau0 = RhizoPara(6)
    gamm = RhizoPara(7)
    thtS = bulkPara(2)
    maxIter = 500
    tolerance = 1e-4
    iter = 0
    diffTht = 1.0
    thetaNonEq(i) = theta(i)*0.5
    ! The main loop that calculate nonequilibrium parameters: 
    ! At each iteration the the equilibrium presure head at the rhizosphere 
    ! is calculate. Following, tau, the mucilage water content and 
    ! the total water content  are calculated. 
    ! The loop is ended when difference in the local water content between two successive
    ! iteration is lower than the tolerance or the iteration exceed the maximum
    ! allowablw iteration. 
    Do While (Abs(diffTht) .Gt. tolerance .And. iter .Lt. maxIter)
       Call hEq(thetaNonEq(i), hnew, i) ! Calculating the new equilibrium pressure head
       diffTht = thetaTot(i)
       tauTht(i) = thetaNonEq(i)**(-gamm)*tau0
       ! d(theta)\dt = 1/tau * h-heq ! 
       thetaNonEq(i) = thetaNonEqOld(i) + (1./tauTht(i))*(hnew - hEqRhizo(i))*dt
       If (thetaNonEq(i) .Ge. thtS) Then
          thetaNonEq(i) = thtS
          thetaTot(i) = Rnorm(i)*theta(i) + (1. - Rnorm(i))*thetaNonEq(i)
          Goto 251
       Endif
       thetaTot(i) = Rnorm(i)*theta(i) + (1. - Rnorm(i))*thetaNonEq(i)
       diffTht = diffTht - thetaTot(i) ! difference between thetaTot within two iterations.
       iter = iter + 1
       If (iter .Ge. maxIter .And. Abs(diffTht) .Gt. tolerance) Then
          Write(*,*) '!!!-!!!'
          Write(*,*) 'rhizoERR = ' , diffTht
          Goto 251
       Endif
    Enddo
251 Continue
    !    WRITE(*,*) iter, diffTht, thetaTot(i), thetaNonEq(i), theta(i)
  End Subroutine thtNonEq
!End Module Water
