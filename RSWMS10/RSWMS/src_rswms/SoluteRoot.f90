! ==============================================================================
! Source file SOLUTE ROOT ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!>  \file SoluteRoot.f90
!!  \brief module for particle tracking inside the roots
!> Solute root model implementation through 9 Subroutines
!>      - SoluteRoot (main)
!>- MoveParticles
!>      - MoveSingleParticle
!>- InsertParticle
!>      - RemoveParticle
!>- GenerateHormones
!>      - UptakePassiv
!>      - OutParticleVTK (Output.f90)
!>      - UpdateRootSys
!>      - ReduceParticleMass

Module ParticlesInRoot

Contains

  !***************************************************************************
  !> main subroutine of solute transport inside root
  !> calculates for one time step the input and the moving
  !> of the particles
  Subroutine SoluteRoot(iCount,dt,t,ipl)
    Use TypeDef
    Use ParamData, Only:  maxParticle, pi, lPartUp
    Use RootData
    Use SoluteRootMat
    Use PlntData, Only: Tact
    Implicit None

    Integer(ap),Intent(in) :: icount
    Integer(ap) :: irec, numVic
    Integer(ap):: ipl
    Real(dp),Intent(in) :: dt,t
    Type(Particle) ,Pointer  :: P
    !> \param icount icount=0 at first time step
    !> \param dt time step size
    !> \param t current simulation time

    !* insert new particles
    m_in = 0._dp
    numVic = 0
    res_t = 0._dp
    !If ((icount.Eq.0) .Or. (totalParticleNum .Eq. 0)) Nullify(firstP)
    If(lSign_inst.Or.lSign) Then
       Call GenerateHormones(1_ap,dt,t)
    Elseif(lPartUp) Then
       Call UptakePassiv(1_ap,dt,t)
    End If

    !* move particles
    If (Associated(firstP)) Call MoveParticles(dt,icount,numVic,t)

    !* calculate retardation factor for sorption
    If(l_linSorb .Or. l_freundSorb) Call CalculateRootRetard(ipl)

    !* calculate mass per root segment
    p=> firstP
    segSoluteMass = 0.
    Do While(Associated(P))
       segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass
       P => P%NEXT
    End Do

    root_mass=sum(segSoluteMass(1:nrec(ipl)))

    !* calculate concentration in root segment
    Do irec=1,nrec(ipl)
       If (seglen(irec,ipl) .Ge. 1.E-20) Then ! skip this segment too small to be taken
          segconc(irec) = segSoluteMass(irec)/segvol(irec,ipl)  !segconc in [M/L³tot]
       End If
    End Do


    !* first root segment acts as buffer
    If(lPartUp) mcol = mcol + m_in !mass gain of buffer
    !print*,'mcol',mcol

    !* calculate cumulated conc. at root collar (for the last time step)
    If(lSign) Then
       concol = mcol/vol_buff ! calculate  conc. at root collar (for the last time step)
       mcol = Max(mcol - concol*Tact(1)*dt,0._dp) !* mass loss of buffer
    End If

  End Subroutine SoluteRoot

  !***************************************************************************
  !> move a single particle by a given velocity
  !> calculates the new position, new segment number and
  !> removes particle form list if it leaves the root (collar)
  Subroutine MoveSingleParticle(P,ipl,dt,numVic,t)
    Use TypeDef
    Use DoussanMat, Only: veloRoot
    Use RootData, Only: irecpr,seglen
    Use SoluteRootMat, Only: Particle, totalParticleNum,irecfollow,numfollow,retard
    Use ParamData, Only: lPartUp
    Implicit None

    Integer(ap) ::ipl,numVic
    Real(dp) ,Intent(in) ::dt,t
    Real(dp):: moveVec,tm,lm1,lmt
    Logical :: lforw
    Type(Particle) ,Pointer ,Intent(inout) :: P!,pFirst
    Type(Particle) ,Pointer :: Pprev
    !> \param pFirst pointer to first particle in particlelist
    !> \param P pointer to the current particle
    !> \param ipl current root system (not used)
    !> \param dt time step size
    !> \param numVic number of particle (victims) which will be removed in the next time step
    !> \param t current simulation time


    tm  = 0. ! time a particle already was moved
    lm1 = 0. ! moved length in segment
    lmt = 0. ! moved length within one time step

    Do While((tm.Lt.dt ))
       moveVec    =  veloRoot(P%segNum,ipl) *retard(P%segNum) * (dt-tm)
       P%position = P%positionOld + moveVec

       If (moveVec .Ge. 0) Then !* forward movement
          lforw = .True.
          lm1 = seglen(P%segNum,ipl) - P%positionOld
          lmt = lmt + lm1
          tm  = tm + lm1/(veloRoot(P%segNum,ipl)*retard(P%segNum))

          If (P%position .Gt. seglen(P%segNum,ipl)) Then ! Change segment
             P%segNum = irecpr(P%segNum,ipl)

             If (P%segNum.Eq.0)  Then
                Call RemoveParticle(P,Pprev,lforw,numVic,t)
                totalParticleNum =  totalParticleNum -1
                P => Pprev
                Exit
             End If

             P%segLen = seglen(P%segNum,ipl)
             P%positionOld = 0.

          End If


       Else !* backward movement
          lforw = .False.
          lm1   = P%positionOld
          lmt   = lmt + lm1
          tm    = tm + lm1/Abs(veloRoot(P%segNum,ipl)*retard(P%segNum))

          If(P%position .Lt. 0) Then !Change segment

             If (numfollow(P%segNum).Gt.0) Then
                P%segNum = irecfollow(P%segNum,1) !CHANGE HERE FOR BIFURCATION!
                P%positionOld = seglen(P%segNum,ipl)
             Else !root tip reached
                If(lPartUp) Then
                   !* Particles get stuck in end of segment
                   p%position    = 0.0
                   P%positionOld = 0.0
                   P%segNum      = P%segNum
                   p%segLen      = segLen(P%segNum,ipl)
                   tm = dt
                Else
                   !* hormonal signalling
                   Call RemoveParticle(P,Pprev,lforw,numVic,t)
                   totalParticleNum =  totalParticleNum -1
                   P => Pprev
                   Exit
                End If
             End If
          End If

       End If

    End Do

  End Subroutine MoveSingleParticle

  !***************************************************************************
  !> go over list of particles and move each particle
  !> with a given velocity
  Subroutine MoveParticles(dt,icount,numVic,t)
    Use TypeDef
    Use ParamData, Only:  maxParticle, pi
    Use DoussanMat, Only: nplant
    Use RootData, Only: segSoluteMass,lno_root_growth,m_in,res_t
    Use SoluteRootMat
    Implicit None

    Integer(ap) :: ipl,numVic
    Integer(ap), Intent(in) :: icount
    Real(dp) ,Intent(in) :: dt,t
    !Type(Particle) ,Pointer ,Intent(inout) :: pFirst

    !> \param pFirst pointer to first particle in particlelist
    !> \param icount icount=0 at first time step
    !> \param dt time step size
    !> \param numVic number of particle (victims) which will be removed in the next time step
    !> \param t current simulation time

    m_in=0._dp
    numVic = 0
    res_t = 0._dp

    If(icount.Eq.0)  Then
       Call MakeFollowList(ipl)
    Elseif(.Not.lno_root_growth) Then
       Call MakeFollowList(ipl)
    End If


    Do ipl=1,nplant
       segSoluteMass = 0.
       pParticle => firstP !point to beginning of the list
       Do While (Associated(pParticle)) !end of list reached

          Call MoveSingleParticle(pParticle,ipl,dt,numVic,t)!move one particle
          If (Associated(pParticle)) Then
             pParticle%positionOld = pParticle%position
             pParticle => pParticle%next !go to next link in the list
          Else
             pParticle => firstP !pFirst
          End If
       End Do
    End Do

  End Subroutine MoveParticles

  !***************************************************************************
  !> this routine creates a new particle (new entry in linked list)
  !> which is added at the beginning (left hand side) of the list
  Subroutine InsertParticle(irec,pos,mass,t,ipl)
    Use TypeDef
    Use RootData, Only: seglen
    Use SoluteRootMat, Only: Particle,firstP
    Implicit None

    Integer(ap),Intent(in):: irec
    Integer(ap),Save :: particleID=0
    Integer(ap):: ipl
    Real(dp),Intent(in):: t
    Real(dp),Intent(in):: mass
    Real(dp),Intent(in):: pos
    !Type(Particle) ,Pointer,Intent(inout) :: pFirst
    Type(Particle) ,Pointer :: pParticle
    !> \param pFirst pointer to first particle in particlelist
    !> \param irec current root segment
    !> \param pos location of particle within root segment
    !> \param mass mass of current particle
    !> \param t current simulation time

    particleID =  particleID +1
    Allocate(pParticle)      ! new particle

    Nullify(pParticle%prev)  ! prev of new particle points to Null
    pParticle%next => firstP !pFirst ! next of new particle is pFirst

    ! prev of pFirst points to new particle
    !If(Associated(pFirst)) pFirst%prev => pParticle
    If(Associated(firstP)) FirstP%prev => pParticle
    FirstP => pParticle      ! pFirst is now new particle

    pParticle%ID = particleID   ! ID
    pParticle%segNum = irec     ! segNum
    pParticle%positionOld = pos ! old position (same as new)
    pParticle%position = pos    ! new position
    pParticle%mass = mass      ! particle mass
    pParticle%segLen = seglen(irec,ipl)
    pParticle%partOrig = t      !origination time of particle


  End Subroutine InsertParticle

  !***************************************************************************
  !> this routine removes a particle (entry in linked list) and
  !> links the pointers of the prev and next particle (entry in list)
  !> adds up particle masses that are removed at the collar
  !> calculates average residence time of particles
  Subroutine RemoveParticle(pVictim,pVictim_prev,lforw,numVic,t)
    Use TypeDef
    Use SoluteRootMat, Only: Particle,firstP
    Use RootData, Only: m_in,res_t
    Implicit None

    Integer(ap) :: numVic
    Real(dp),Intent(in) :: t
    Logical :: lforw
    !Type(Particle) ,Pointer ,Intent(inout) :: pFirst
    Type(Particle) ,Pointer ,Intent(inout) :: pVictim
    Type(Particle) ,Pointer ,Intent(out) :: pVictim_prev
    !> \param pFirst pointer to first particle in particlelist
    !> \param pVictim particle to be removed from the system
    !> \param pVictim_prev previous particle from pVictim
    !> \param lforw true if upward transport/flow
    !> \param t current simulation time
    !> \param numVic number of particles to be removed from the system

    If (lforw) Then
       m_in=m_in+pVictim%mass
       res_t=(res_t*numVic+t-pVictim%partOrig)/(numVic+1)
       numVic=numVic+1
    End If

    If (Associated(pVictim%prev)) Then
       ! pVictim has a prev entry
       pVictim%prev%next => pVictim%next
       pVictim_prev => pVictim%prev
    Else
       !if pVictim = pFirst -> there is no prev: set pFirst to next entry
       !pVictim_prev => pFirst%prev
       !pFirst => pVictim%next
       pVictim_prev => firstP%prev
       Firstp => pVictim%next
    End If

    If (Associated(pVictim%next)) Then
       ! pVictim has a next entry
       pVictim%next%prev => pVictim%prev
    Else
       ! if last entry of the list -> set next pointer from prev to Null
       ! but only if it is not the only (very last) particle
       If (Associated(pVictim%prev)) Nullify(pVictim%prev%next)
    End If

    !* remove pVictim
    Deallocate(pVictim)

  End Subroutine RemoveParticle

  !***************************************************************************
  !> this routine calculates the hormone mass which is generated by the root segments
  !> and accordingly inserts particles in the roots
  Subroutine GenerateHormones(ipl,dt,t)
    Use TypeDef
    Use RootData
    Use SoluteRootMat
    Use PlntData, Only: Tact
    Use DoussanMat, Only: Phr,w_sub,nsub
    Implicit None

    Integer(ap):: irecn,ipl,ip,isub
    Real(dp):: xp, posi
    Real(dp), Intent(in):: dt,t
    !Type(Particle) ,Pointer ,Intent(out) :: pFirst
    !> \param pFirst pointer to first particle in particlelist
    !> \param ipl current root system (not used)
    !> \param dt time step size
    !> \param t current simulation time

    Do irecn=1,nrec(ipl)
       mhorm=0._dp
       msign_notrans=0._dp
       csign_notrans=0._dp
       
       If(Abs(Phr(irecn+1,ipl)).Gt. TR1)  mhorm(irecn,ipl) = &
            sign_in*(Abs(Phr(irecn+1,ipl))-TR1)*dt*segmas(irecn,ipl)
    End Do
    
  
    msign_notrans=Sum(mhorm)
    If (msign_notrans .Gt. 0._dp) Then
       mcol_i = mcol_i + msign_notrans
       csign_notrans = mcol_i/vol_buff
       mcol_i = mcol_i - csign_notrans*Tact(1)*dt
    End If

    If(lSign_inst) Return

   
    Do irecn=1,nrec(ipl)
       If(mhorm(irecn,ipl).Gt.1E-22) Then
          Do isub=1,nsub(irecn,ipl)
             If (isub.Eq.1) Then
                xp = 0.
             Else
                xp = xp + seglen(irecn,ipl)*w_sub(irecn,isub-1,ipl)
             End If
             Do ip=1,numPa
                ! position of particle in root segement
                posi = (ip-1)*(seglen(irecn,ipl)*w_sub(irecn,isub,ipl)/numPa)+xp
                Call InsertParticle(irecn,posi,mhorm(irecn,ipl)/numPa,t,ipl)
                totalParticleNum =  totalParticleNum +1
             End Do
          End Do
       End If
    End Do

  End Subroutine GenerateHormones

  !***************************************************************************
  !> calculates the mass which is taken up passively (advection+diffusion) by a root segment
  !> and accordingly inserts particles inside the root
  Subroutine UptakePassiv(ipl,dt,t)
    Use TypeDef
    Use DoussanMat, Only: nsub,w_sub,cube_i,Joutr,isubmax
    Use RootData
    Use GridData, Only: dzgrid,dygrid,dxgrid, nElm, elmnod
    Use SoluteRootMat
    Use SolData, Only: theta,SoilSoluteMass
    Implicit None

    Integer(ap):: ibr,irecn,ifoln,iprvn,ip,isub,ipl,cor(8),iE
    Real(dp):: xp, posi, VoxVol, RootUptake_adv,RootUptake_diff
    Real(dp):: delta_ru,dummy
    Real(dp),ALLOCATABLE,DIMENSION(:,:) :: RootUptake
    Real(dp),ALLOCATABLE,DIMENSION(:) :: RS_uptake
    Real(dp), Intent(in):: dt,t
    Logical:: n_apex,run, corr_root
    Type(Particle) ,Pointer :: P
    !Type(Particle) ,Pointer ,Intent(inout) :: pFirst
    !> \param pFirst pointer to first particle in particlelist
    !> \param ipl current root system (not used)
    !> \param dt time step size
    !> \param t current simulation time

    ! uptake by roots - can change in case of root growth
    IF (.NOT. ALLOCATED(RootUptake) .OR. .NOT. lno_root_growth) THEN
       ALLOCATE(RootUptake(nrec(ipl),isubmax))
    END IF
    RootUptake = 0._dp
    IF (.NOT. ALLOCATED(RS_uptake) .OR. .NOT. lno_root_growth) THEN
       ALLOCATE(RS_uptake(nelm))
    END IF
    RS_uptake = 0._dp
    
    !SoilSoluteUptake = 0._dp
    seg_upt=0._dp
    RootUptake_adv = 0._dp
    RootUptake_diff = 0._dp
    SoilUptake_adv = 0._dp
    SoilUptake_diff = 0._dp
    delta_ru = 0._dp
    
    corr_root = .False.

    VoxVol=dzgrid*dygrid*dxgrid

    !* in case of root growth update matrices
    If(.Not.lno_root_growth) Call UpdateRootSys(ipl)

    !* Soil Voxel water content for solute concentration
    Do iE=1,nElm
       cor(1)=elmnod(1,iE)
       cor(2)=elmnod(2,iE)
       cor(3)=elmnod(3,iE)
       cor(4)=elmnod(4,iE)
       cor(5)=elmnod(5,iE)
       cor(6)=elmnod(6,iE)
       cor(7)=elmnod(7,iE)
       cor(8)=elmnod(8,iE)

       theta_elm(iE) = (theta(cor(1))+theta(cor(2))+theta(cor(3))+theta(cor(4))+theta(cor(5))+&
            theta(cor(6))+theta(cor(7))+theta(cor(8)))/8
       ccube(iE) = SoilSoluteMass(iE)/(VoxVol*theta_elm(iE))
    End Do

    !* Get mass in root segments
    P=> firstP     !pointer from SoluteRootMat
    segSoluteMass=0._dp
    Do While(Associated(P))
       segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass
       P => P%NEXT
    End Do
    !print*,'segsolutemass in', segsolutemass(1:nrec)

    Do ibr=1,nbr(ipl)
       n_apex=.False.
       !* find the tip segment of the branch 'ibr'
       irecn=nrec(ipl)
       Do While (ibrseg(irecn,ipl).Ne.ibr)
          irecn=irecn-1
       End Do
       !* the first one is an apex
       If (seglen(irecn,ipl)<1.E-20) Then ! skip this segment too small to be taken
          ifoln=irecn ! following node ID
          irecn=irecpr(irecn,ipl) !current node ID
       Else
          n_apex=.True.!ifoln does not exist
       Endif
       If (irecn==0) Then !there exists a branch ibr but not yet any segment!
          run=.False.
       Else
          run=.True.
       Endif
       !* then the rest of the branch up to the seed or the embranchment
       Do While (run)
          !* "upper" node
          iprvn=irecpr(irecn,ipl)
          ! if reached the seed or the embranchement => change branch
          If (iprvn.Eq.0 .Or. (ibrseg(iprvn,ipl).Ne.ibrseg(irecn,ipl))) run=.False.

          xp = 0.
          Do isub=1,nsub(irecn,ipl)
             iE=cube_i(irecn,isub,ipl) !element number

             If (isub.Eq.1) Then
                xp = 0.
             Else
                xp = xp + seglen(irecn,ipl)*w_sub(irecn,isub-1,ipl)
             End If

			 !* advective uptake
             If(Joutr(irecn,ipl).Ge.0) Then !* uptake
                RootUptake_adv = frac_pass(irecn)*Joutr(irecn,ipl)*ccube(iE)*w_sub(irecn,isub,ipl)*dt
             Elseif (Joutr(irecn,ipl).Lt.0) Then  !* exudation
                RootUptake_adv = frac_pass(irecn)*Joutr(irecn,ipl)*(segconc(irecn)-segsorb(irecn))/theta_R*&
                     w_sub(irecn,isub,ipl)*dt
             End If

			 !* diffusive uptake
             RootUptake_diff = (ccube(iE)-(segconc(irecn)-segsorb(irecn))/theta_R)*segsur(irecn,ipl)* &
                  w_sub(irecn,isub,ipl)*Perm(irecn)*dt
             !if (RootUptake_diff<0._dp) RootUptake_diff = 0._dp !switch off exudation 
             !ccube is in [M/L³water]
             !segconc and segsorb is in [M/L³tot] -> divide by theta_R to receive cw=[M/L³wat]

             !* potential uptake for subsegment
             RootUptake(irecn,isub) = RootUptake_adv + RootUptake_diff

             SoilUptake_adv(iE) = SoilUptake_adv(iE) - RootUptake_adv
             SoilUptake_diff(iE) = SoilUptake_diff(iE) - RootUptake_diff
             RootUptake_adv = 0._dp
             RootUptake_diff = 0._dp
             
             !* Total root surface that takes up solutes per soil element
             if(RootUptake(irecn,isub).GT.0._dp) RS_uptake(iE) = RS_uptake(iE)+segsur(irecn,ipl)*w_sub(irecn,isub,ipl)
             
          End Do !sub segments


          ! definition for the next run of the loop
          ifoln=irecn
          irecn=iprvn
          !* from here, not an apex
          n_apex=.False.
       End Do !nodes
    End Do ! branches



    !* potential soil solute uptake is calculated
    ! positive for exudation, negative for uptake   
    DO iE = 1,nelm
       SoilSoluteUptake(iE) = SoilUptake_adv(iE)+SoilUptake_diff(iE)
    END DO

    !print*,'sum root uptake pot', sum(RootUptake_pot(1:nrec,1:nsub(nrec,1)))
    !print*,'sum soil uptake pot', sum(soilsoluteuptake(1:nelm)) 
    !print*, 'soil solute mass sum', sum(SoilSoluteMass(1:nelm))
    !print*, 'sum root solute mass', sum(segSoluteMass(1:nrec))


    !* Check that root exudation is not limited by root solute mass
    DO irecn=1,nrec(ipl)
       DO isub=1,nsub(irecn,ipl)
          DO iE=1,nelm
             IF(iE.EQ.cube_i(irecn,isub,ipl)) THEN
                IF ((RootUptake(irecn,isub).LT.0._dp)) THEN !exudation
                   delta_ru = segSoluteMass(irecn)*w_sub(irecn,isub,ipl)+RootUptake(irecn,isub)
                   IF (delta_ru .LT. 0._dp) THEN !exudation > available root mass
                      !* update soilsoluteuptake and soil solute mass for this element
                      !print*,irecn,isub,'SSU in',SoilSoluteUptake(iE),'delta_ru',delta_ru,'soilsm in',soilsolutemass(ie)
                      SoilSoluteUptake(iE) = SoilSoluteUptake(iE) + delta_ru
                      !SoilSoluteMass(iE) = SoilSoluteMass(iE) + delta_ru
                      !print*, 'soilsm out',soilsolutemass(ie)
                      !* limit root uptake per subsegment
                      RootUptake(irecn,isub) = -segSoluteMass(irecn)*w_sub(irecn,isub,ipl)
                      !print*,'ru out',rootuptake(irecn,isub)
                   END IF
                END IF
             END IF
          END DO
       END DO
    END DO
   
    !* Check that there is no negative soil solute mass...
    ! soil solute uptake <0 when solutes are removed (taken up by roots)
    dummy=0._dp
    DO iE=1,nelm
       IF(SoilSoluteMass(ie).GT.0._dp)THEN
          IF (SoilSoluteMass(iE)+SoilSoluteUptake(iE).LT.0._dp) THEN               
             !* limit soilsoluteuptake for roots that take up solutes (RootUptake>0)
             !* Scale the available soil mass per root subsegment with the relative root surface, based on the total uptake surface
             !print*,'soilsm in', Soilsolutemass(ie),'soilupt in',soilsoluteuptake(ie)
             SoilSoluteUptake(iE) = -(1._dp-(1._dp/1000._dp))*SoilSoluteMass(iE)
             DO irecn=1,nrec(ipl)
                DO isub=1,nsub(irecn,ipl)
                   IF(cube_i(irecn,isub,ipl).EQ.iE) THEN
                      IF(RootUptake(irecn,isub).GT.0._dp) THEN
                         !print*,'rootupt in',rootuptake(irecn,isub)
                         RootUptake(irecn,isub) = (1._dp-(1._dp/1000._dp))*SoilSoluteMass(iE)*segsur(irecn,ipl)*w_sub(irecn,isub,ipl)/RS_uptake(iE)
                         dummy = dummy+RootUptake(irecn,isub)
                         !print*,'rootupt out',rootuptake(irecn,isub),'soilsm',SoilSoluteMass(ie),'scaling',segsur(irecn)*w_sub(irecn,isub,ipl)/RS_uptake(iE),'dummy',dummy
                      END IF
                   END IF
                END DO
             END DO
             !print*,'dummy',dummy,'SoilSoluteUptake + dummy',SoilSoluteUptake(iE)+dummy
             dummy=0._dp
          END IF
       END IF
    END DO
    
    
    !*update segupt
    seg_upt = 0._dp
    DO irecn=1,nrec(ipl)
       DO isub=1,nsub(irecn,ipl)
          seg_upt(irecn) = seg_upt(irecn)+RootUptake(irecn,isub)
       END DO
    END DO
        
    !* Calculate fact to submit to ParTrace
    !* with initial soil mass
    Do iE=1, nElm
       fact(iE)=1._dp
       If(Abs(SoilSoluteUptake(iE)).Ge.1e-22 .AND. SoilSoluteMass(iE).GT.0._dp) Then
             fact(iE) = (1+SoilSoluteUptake(iE)/SoilSoluteMass(iE)) 
          Else
             fact(iE) =1._dp
       End If
    End Do
    !print*,'fact',fact(1:nelm)
    print*,'sum soilsupt out',sum(soilsoluteuptake(1:nelm)) 
       
	!* Calculate root solute mass
    !* insert particles in root
    Do irecn=1,nrec(ipl)
       !print*,'IN irecn',irecn,'seg solute mass', segsolutemass(irecn),'segupt',seg_upt(irecn)
       If (seg_upt(irecn).Gt.1e-22) Then
          p=>firstp
          !print*, 'insertParticle'
          Do ip=1,numPa
             ! position of particle in root segement
             posi = (ip-1)*(seglen(irecn,ipl)*w_sub(irecn,isub,ipl)/numPa)+xp
             Call InsertParticle(irecn,posi,seg_upt(irecn)/numPa,t,ipl)
             totalParticleNum =  totalParticleNum +1
          End Do

       Else If  (seg_upt(irecn).Le.-1e-22) Then   ! mass leaves the root segment
          !Print*, 'ReduceMass'
          If (segSoluteMass(irecn).Le.Abs(seg_upt(irecn)))  corr_root=.True.
          Call ReduceParticleMass(seg_upt(irecn), irecn, corr_root)
       End If
    End Do

    !* update root solute mass
    P=> firstP     !* pointer from SoluteRootMat
    segSoluteMass=0._dp
    Do While(Associated(P))
       segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass
       P => P%NEXT
    End Do
    !print*,'sum root solute mass out', sum(segsolutemass(1:nrec))
     
     IF(.NOT. lno_root_growth) DEALLOCATE(RootUptake)

  End Subroutine UptakePassiv
   !***************************************************************************
  !> this routine stores a list of the connection of root segments
  Subroutine MakeFollowList(ipl)
    Use TypeDef
    Use RootData, Only: nrec, nbr, ibrseg, irecpr, ibrseg, seglen
    Use SoluteRootMat, Only: irecfollow,numfollow
    Implicit None

    Integer(ap) ::ibr,irecn,ifoln,iprvn,a,ipl
    Logical :: n_apex,run

    irecFollow = 0
    numFollow = 0

    Do ibr=1,nbr(ipl)
       n_apex=.False.
       !* find the tip segment of the branch 'ibr'
       irecn=nrec(ipl)
       Do While (ibrseg(irecn,ipl).Ne.ibr)
          irecn=irecn-1
       End Do
       !* the first one we find is an apex
       If (seglen(irecn,ipl)<1.E-20) Then ! skip this segment too small to be taken
          ifoln=irecn ! following node ID
          irecn=irecpr(irecn,ipl) !current node ID
       Else
          n_apex=.True.!ifoln does not exist
       Endif
       If (irecn==0) Then !there exists a branch ibr but not yet any segment!
          run=.False.
       Else
          run=.True.
       Endif
       !* then the rest of the branch up to the seed or the embranchment
       Do While (run)
          !* "upper" node
          iprvn=irecpr(irecn,ipl)
          ! if reached the seed or the embranchement => change branch
          If (iprvn==0.Or.(ibrseg(iprvn,ipl).Ne.ibrseg(irecn,ipl)))Then !run=.FALSE.
             Exit
          Endif
          !number of following root segments (branches)
          a = numFollow(irecpr(irecn,ipl))+1
          !list which root segment follows
          irecfollow(irecpr(irecn,ipl),a) = irecn
          numFollow(irecpr(irecn,ipl)) = a
          ! definition for the next run of the loop
          ifoln=irecn
          irecn=iprvn
          !* from here, not an apex
          n_apex=.False.
       End Do !nodes
    End Do ! branches

  End Subroutine MakeFollowList

  !***************************************************************************
  !> calculates a retardation factor for the solute particles in case of sorption
  Subroutine CalculateRootRetard(ipl)
    Use typedef
    Use RootData, Only: nrec,segconc
    Use SoluteRootMat, Only: retard, sorp, l_linSorb, l_freundSorb, segsorb , rho_R, theta_R, segsolv
    Implicit None

    Integer(ap) :: ii,irecn,ipl
    Real(dp) :: epsilon = 1e-5, cw, r, c, rdata, k, n

    !> calculates retardation within roots due to sorption
    !> assumes fully saturated roots; kd is not scaled with water content
    !> references in Partrace / particleclassSorp.cpp and dissertation of
    !> O. Neuendorf

    If(l_linSorb) Then
       Do irecn=1,nrec(ipl)
          k = sorp(1)
          retard(irecn) = 1/(theta_R + rho_R*k)
          segsorb(irecn)=(1-retard(irecn))*segconc(irecn)
       End Do

    Elseif(l_freundSorb) Then

       Do irecn=1,nrec(ipl)
          c = segconc(irecn)! in [M/L³tot]
          rdata = retard(irecn)
          k = sorp(1)/rho_R*theta_R
          n = sorp(2)


          !* change to water volume based concentrations

          c=c/theta_R  !in [M/L³wat]
          !print*, c, rdata, retard(irecn)

          If (rdata > 0._dp) Then
             cw = rdata * c
          Else
             cw = 0.5 * c
          End If

          Do ii=1,50 !max. 50 iterations
             r = cw-(cw+k*cw**n-c)/(1._dp+n*k*cw**(n-1));

             If(r<0) Then
                cw = cw/2
             Else
                If(Abs((cw-r)/cw)<epsilon) Then
                   cw = r
                   Goto 10
                Else
                   cw = r
                End If
             End If
          End Do ! iterations

10        If(ii==50) Then
             rdata = 1._dp

          Elseif  (c==0._dp) Then
             rdata = 1._dp

          Else
             rdata = cw/c !rdata in [-] (M/L³wat)
          End If

          retard(irecn) = rdata*theta_R !retard in [-] [M/L³tot]

          segsorb(irecn) = (c-cw)*theta_R !segsorb in [M/L³tot]
          segsolv(irecn)= cw !segsolv in [M/L³tot]

       End Do ! loop over nrec

    End If ! sorption type

  End Subroutine CalculateRootRetard

  !******************************
  Subroutine UpdateRootSys(ipl)
    Use typedef
    Use RootData, Only: seglen, segsur, nrec, segvol, vol_root,nplant
    Use ParamData, Only: Pi
    Implicit None

    Integer(ap) :: irec,ipl
    Real(dp) :: crossSectionSeg(nrec(nplant),nplant) , segrad(nrec(nplant),nplant)
    !> When root growth, update segment information:
    !> - segment volume
    !> - buffer volume

    Do irec=1,nrec(ipl)
       segrad(irec,ipl) = segsur(irec,ipl)/2._dp/pi/seglen(irec,ipl)
       crossSectionSeg(irec,ipl) = segrad(irec,ipl)**2*pi
       segvol(irec,ipl) =  segrad(irec,ipl)**2*pi*seglen(irec,ipl) !root volume per root segment
       vol_root(ipl) = vol_root(ipl) +  segvol(irec,ipl) !total root volume
    End Do


  End Subroutine UpdateRootSys

  !***************************************************************************
  Subroutine ReduceParticleMass(m_red, irecn, corr_root)
    Use typedef
    Use RootData, Only:segSoluteMass, nrec,nplant
    Use SoluteRootmat, Only:Particle, firstP
    Use MPIutils, Only: stop_program
    Implicit None

    Integer(ap):: irecn
    Real(dp)::red_fact(nrec(nplant)) , m_red
    Logical:: corr_root
    Type(Particle) ,Pointer :: P

    !* Calculate Reduction Factor red_fact
    red_fact=1._dp

    P => firstP
    segSoluteMass=0._dp
    Do While(Associated(P))
       segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass   !segsolutemass length of maxseg
       P => P%NEXT

    End Do

    If (corr_root) Then
       red_fact(irecn)=0._dp
    Else If (segSoluteMass(irecn).Ge.1e-22) Then
       red_fact(irecn)=1+(m_red/segSoluteMass(irecn))

    Elseif(segSoluteMass(irecn).Lt.1e-22) Then
       red_fact(irecn)=1.0_dp  ! produces mass!

    Endif

    P=>firstP
    segSoluteMass=0._dp
    Do While(Associated(P))
       !print*, P%Mass, 'particle mass', P%segNum, 'ParticleSegNum', red_fact(P%segNum), 'red_fact'
       P%Mass=P%Mass*red_fact(P%segNum)
       segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass
       P => P%Next
    End Do
    !print*,'reduce2',' segsolutemass',segsolutemass(irecn),'irecn',irecn
    If(red_fact(irecn).Eq.1) Call stop_program('red_fact(irecn).EQ.1')

  End Subroutine ReduceParticleMass

  !********************************************************************************
End Module ParticlesInRoot
