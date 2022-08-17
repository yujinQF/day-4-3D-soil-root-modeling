 MODULE mumps_mod

CONTAINS

  SUBROUTINE MUMPS(NN,NNZ,irow,jcol,aij,rhs,sol,ipl)

    USE typedef
    USE tmctrl, ONLY:  tlevel
    USE DoussanMat, ONLY: curr_BCtp, nBCn, nBC_irecn, Khr, curr_BCr
    USE RootData, ONLY: seglen 
    IMPLICIT NONE
    include 'dmumps_struc.h'

    INTEGER(ap), INTENT(IN)::NN, NNZ, ipl
    INTEGER(ap), INTENT(INOUT)::irow(NNZ),jcol(NNZ)
    INTEGER(ap) ::iprvn, irecn, ibranch
    REAL(dp), INTENT(INOUT)::aij(NNZ),rhs(2*NN),sol(NN)
    REAL(dp)::sol1(NN),sol2(NN),Q1(NN),Q2(NN), alpha
    TYPE (DMUMPS_STRUC), SAVE:: mumps_par !take care: mumps_par is save!

    !if first iter step
    IF (.NOT.tlevel)THEN

       mumps_par%JOB = -1
       mumps_par%SYM = 0
       mumps_par%PAR = 1

       !job = -1 -> initializes an instance of the package
       CALL DMUMPS(mumps_par)

       !supress output
       mumps_par%ICNTL(1) = -5
       mumps_par%ICNTL(2) = -5
       mumps_par%ICNTL(3) = -5
       mumps_par%ICNTL(4) = -5

       !matrix entries
       mumps_par%N = NN
       mumps_par%NZ = NNZ
       ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
       ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
       ALLOCATE( mumps_par%A( mumps_par%NZ ) )


       !ini matrix in coordinate format
       mumps_par%IRN = irow
       mumps_par%JCN = jcol
       mumps_par%A = aij

       !use AMD minimal degree ordering
       mumps_par%ICNTL(7) = 0

       !job = 4 -> analysis and factorisation of the matrix
       mumps_par%JOB = 4
       CALL DMUMPS(mumps_par)

       !after first step set tlevel = true
       tlevel=.TRUE.

    ENDIF

    !if flux bc, two rhs´es are needed for the solution
    IF (curr_BCtp(ipl)==2) THEN  !Flux
       mumps_par%NRHS = 2
       mumps_par%LRHS = NN
       ALLOCATE( mumps_par%RHS ( mumps_par%N * mumps_par%NRHS ) )

       !initialize rhs 
       !take care to create two right hand sides 
       mumps_par%RHS = rhs


       !solve again with same matrix
       !analysis and factorisation is stored
       mumps_par%JOB = 3
       CALL DMUMPS(mumps_par)

       sol1 = mumps_par%RHS(1:NN)  !solution first rhs
       sol2 = mumps_par%RHS(NN+1:) !solution second rhs


       iprvn = 0
       DO ibranch=1,nBCn(ipl)  !all the branches connected to the seed
          irecn=nBC_irecn(ibranch)!connected to the seed

          Q1(irecn) = Khr(irecn,ipl)/seglen(irecn,ipl)*(sol1(iprvn+1)-sol1(irecn+1))
          Q2(irecn) = Khr(irecn,ipl)/seglen(irecn,ipl)*(sol2(iprvn+1)-sol2(irecn+1))
       ENDDO

       alpha =  (curr_BCr(ipl) -  Q2(1) ) / ( Q1(1) -  Q2(1))
       sol = sol1*alpha + (1-alpha)*sol2

       !if ph bc -> solve only with one rhs
    ELSEIF (curr_BCtp(ipl)==1) THEN  !PH

       mumps_par%NRHS = 1
       mumps_par%LRHS = NN
       ALLOCATE( mumps_par%RHS ( mumps_par%N * mumps_par%NRHS ) )

       !initialize rhs 
       !take care to use only first part of rhs
       mumps_par%RHS = rhs(1:NN)

       !solve again with same matrix
       !analysis and factorisation is stored
       mumps_par%JOB = 3
       CALL DMUMPS(mumps_par)

       !solution stored in rhs of mumps struct
       sol = mumps_par%RHS  
    ENDIF

    DEALLOCATE( mumps_par%RHS )

  END SUBROUTINE MUMPS

END MODULE mumps_mod
!******************************************************************************
SUBROUTINE SolveRootDirect(ipl,num_elem,rhs,solution)
  !solve DOUSSAN flow equations within the xylem
  !with direct scheme and save inverse matrix
  USE TypeDef
  USE SparseMatrix
  USE RootData, ONLY :nrec,nplant
  USE DoussanMat, ONLY: jao,iao,ao,iro,jco,aij,plantmatrix !sa,ija,
  USE tmctrl, ONLY:  tlevel
  USE mumps_mod
  IMPLICIT NONE

  INTEGER(ap), INTENT(in):: ipl
  INTEGER(ap)::num_elem
  REAL(dp):: rhs(1:2*nrec(nplant)+2)
  REAL(dp), INTENT(out):: solution(nrec(nplant)+1)

  !CALL SM_write(plantmatrix(ipl))

  IF (.NOT.tlevel)THEN !first iteration step 

     ALLOCATE(iro(num_elem))
     ALLOCATE(jco(num_elem))
     ALLOCATE(aij(num_elem))
     !convert to COO sparse matrix format
     CALL SM_convert(plantmatrix(ipl),iro,jco,aij,num_elem,ipl) 

     ALLOCATE(jao(num_elem))
     ALLOCATE(ao(num_elem))
     ALLOCATE(iao(nrec(ipl)+2))
     !convert from COO to CRS sparse matrix format -> sparskit routine!
     CALL coocsr(nrec(ipl)+1,num_elem,aij,iro,jco,ao,jao,iao) 

     !CALL pspltm(nrec+1,nrec+1,0,jao,iao,'title',1,15.,'cm',0,0,12) !print in fort.12 -> ps file

  END IF

  CALL  MUMPS(nrec(ipl)+1,num_elem,iro,jco,aij,rhs,solution,ipl)


END SUBROUTINE SolveRootDirect
!******************************************************************************
SUBROUTINE SM_convert(mat,irow,jcol,aij,num_elem,ipl)
  USE typedef
  USE SparseMatrix
  USE RootData, ONLY : nrec,nplant
  IMPLICIT NONE

  INTEGER(ap), INTENT(IN) ::ipl, num_elem
  INTEGER(ap), INTENT(INOUT) ::irow(num_elem), jcol(num_elem)
  INTEGER(ap):: jcol_diag(nrec(nplant)+1),irow_diag(nrec(nplant)+1)
  INTEGER(ap)::jcol_off(num_elem-nrec(nplant)-1),irow_off(num_elem-nrec(nplant)-1)
  INTEGER(ap) :: i, k, kk
  REAL(dp), INTENT(INOUT) ::aij(num_elem)
  REAL(dp) ::aij_diag(nrec(nplant)+1),aij_off(num_elem-nrec(nplant)-1)
  TYPE(SparseMatrixType), INTENT(IN) :: mat
  TYPE(SparseMatrixNode), POINTER :: p

  !WRITE(*,*) 'Convert sparse matrix with ',mat%nlines,' lines and ',mat%ncols,' columns'
  !PRINT *, 'Non Zero Elements: ', num_elem
  k = 1
  kk = 1
  DO i=1,mat%nlines
     !Write(*,*) 'line ',i,' has ',mat%lines(i)%numnodes,' elements'      
     p=>mat%lines(i)%first
     DO WHILE(ASSOCIATED(p))
        !Write(*,*) '   column=',p%colno,' data=',p%data
        IF (i.EQ.p%colno)THEN !diagonal elements
           irow_diag(kk) = i
           jcol_diag(kk) = p%colno
           aij_diag(kk) = p%data
           kk = kk + 1      
        ELSE ! off-diagonal elements
           irow_off(k) = i
           jcol_off(k) = p%colno
           aij_off(k) = p%data
           k = k + 1      
        END IF
        p=>p%next
     ENDDO
  ENDDO

  ! put both arrays together
  ! first diagonal elements - then off diagonal elements
  irow(1:nrec(ipl)+1) = irow_diag
  jcol(1:nrec(ipl)+1) = jcol_diag
  aij(1:nrec(ipl)+1) = aij_diag

  irow(nrec(ipl)+2:) = irow_off
  jcol(nrec(ipl)+2:) = jcol_off
  aij(nrec(ipl)+2:) = aij_off

END SUBROUTINE SM_CONVERT
!******************************************************************************
