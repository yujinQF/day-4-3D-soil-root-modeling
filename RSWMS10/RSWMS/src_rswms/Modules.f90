!> \file Modules.f90
!> general information, definition of max. matrix sizes for pre-allocation
    

MODULE ParamData
  USE typedef
  IMPLICIT NONE
  INTEGER(ap),PARAMETER :: maxnod=662661,mxpnts=8,maxplant=5 
  INTEGER(ap),PARAMETER :: maxbdr=13000,mxBcCh=13000,maxIrrig=5,mxBcCh1=50
  INTEGER(ap),PARAMETER :: mxtime=1000,mxdpth=100,maxelm= 640000
  INTEGER(ap),PARAMETER :: maxgrw=100000,maxrec=100000,maxest=5000
  INTEGER(ap),PARAMETER :: maxemg=50,maxord=3,maxmat=10,maxbnd=19
  INTEGER(ap),PARAMETER :: maxobs=30,n_MFP=10000!200000
  INTEGER(ap),PARAMETER :: maxParticle=10000 !solute particles inside root
  INTEGER(ap),PARAMETER:: maxAMG=5000000,maxamg2=5000000 !array size for sparse A matrix and sparse row index matrix
  INTEGER(ap) :: iter_root,iter_tot=0,iter=5,last_out,i_noConv
  REAL(dp), PARAMETER :: pi=3.14159265358979323
  LOGICAL :: lvtk=.FALSE., lOutPartrace=.FALSE., ldirect=.FALSE.,lChem=.FALSE.,lretry=.False.
  LOGICAL :: lOrt, lSalinity=.true.,lPartUp=.FALSE., lclimate=.FALSE., lTard=.FALSE.,NoGrav=.False.
END MODULE ParamData
!************************************************************************
!> data needed for root
MODULE RootData
  USE typedef
  USE paramData, ONLY: maxord,maxplant,maxrec,maxgrw,maxord,maxest,mxpnts,maxemg
  IMPLICIT NONE
  INTEGER(ap) :: nrec_m,ngrow_m,nplant,TypeKdrop
  INTEGER(ap) :: nAQPc,ntimeobs,lPast
  INTEGER(ap), ALLOCATABLE, DIMENSION(:) :: nUrf
  INTEGER(ap), ALLOCATABLE, DIMENSION(:) :: nbr,nrec,ngrow,naxemg,naxes,norder,diffnum,naxtot,nBigLat
  INTEGER(ap), ALLOCATABLE, DIMENSION(:) :: br_rec
  INTEGER(ap), ALLOCATABLE, DIMENSION(:,:):: ordseg,ibrseg,irecpr
  INTEGER(ap), ALLOCATABLE, DIMENSION(:,:):: iaxis,irecsg,nestbl,ordgrw,ibrgrw,num_seg
  INTEGER(ap), ALLOCATABLE, DIMENSION(:,:):: nvch,nMPLch
  INTEGER(ap), ALLOCATABLE, DIMENSION(:,:):: nnewax
  REAL(dp) :: rand !random number created in RSWMS34.f90
  REAL(dp) :: Krs,Kcomp
  REAL(dp) :: sign_in,PH_crit,sAvg,cAvg,vol_buff,size_buff,m_in
  REAL(dp) :: g1,g2,fac_contact,hdrsens=0._dp,Xerosens=-1500._dp
  REAL(dp) :: MPL(maxgrw)=0._dp
  REAL(dp) :: tlim_dJvL,OmegaC
  REAL(dp) :: delbm,LAmsh
  REAL(dp) :: rsr,dr,dsh,grwfac,w,dmroot,rs,concrs,Hseq
  REAL(dp) :: sigma=0.  
  REAL(dp) :: concol,mcol=0,msign_notrans,csign_notrans,res_t,delta_h,mcol_i=0._dp,TR1,root_mass=0
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: xplant,yplant,condMP,vol_root,mroot,mshoot,tot_len
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: segdiam,segSoluteMass,segconc
  REAL(dp), ALLOCATABLE, DIMENSION(:,:):: age,Urf,tnewax,inaxs,timest
  REAL(dp), ALLOCATABLE, DIMENSION(:,:):: brlmax,f_rad,strsen,rdmang,brspac,brnang,dtbrch,lb,maxlast,numlast
  REAL(dp), ALLOCATABLE, DIMENSION(:,:):: xs,ys,zs,timorg,seglen,segsur,segmas,segrad,crossSectionSeg,segvol,mhorm
  REAL(dp), ALLOCATABLE, DIMENSION(:,:):: xg,yg,zg,ovrtime,brlgth,AQPh,AQPv
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)::agevch,vch,sMPLch,MPLch
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: Rho,tPast,fac_plant,timeobs,Krs_mat,Kcomp_mat
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sinkredDPast,distseg2soil
  LOGICAL :: maizeroottyp=.FALSE.,loliumroottyp=.FALSE.,wheatroottyp=.FALSE.,lrrt=.FALSE.,lrrs=.FALSE.,lSomma=.FALSE.,lRootBox=.FALSE.
  LOGICAL :: lDou=.FALSE.,lCou=.FALSE.,lSUF=.FALSE.,lFed=.FALSE.,ldJvL=.FALSE.,lJarvis=.FALSE.,lUrf=.FALSE.,lKdrop=.FALSE.
  LOGICAL :: lGap=.FALSE.,lAQPc=.FALSE.,l_secrad=.FALSE.,l_conduc=.FALSE.,l_overburden=.FALSE.,lno_RWU=.FALSE.,stopgr(maxgrw)=.FALSE.
  LOGICAL :: lno_Archi=.FALSE.,lno_root_growth=.FALSE.,lRootTyp_growth=.FALSE.,lSomma_growth=.FALSE.,lRootBox_growth=.FALSE.,lUpdate_growth=.FALSE.
  LOGICAL :: lCalloc=.FALSE.,ltemp=.FALSE.,lSign=.FALSE.,lSign_new=.FALSE.,lSign_inst=.FALSE.,ltoxi=.FALSE.,lSinkCube=.FALSE.,ltwo_grids=.FALSE.
  LOGICAL :: it1,lHydrotrop=.FALSE.,lXerobranch=.FALSE.,lHydropattern=.FALSE.,lSoilStrength=.FALSE.,lLandl=.TRUE.,lBengough1=.FALSE.
  LOGICAL(dp), ALLOCATABLE, DIMENSION(:) :: toosml, connex,l_SignOn
  Save mcol

CONTAINS

    SUBROUTINE IniRoot

    IMPLICIT NONE

    ALLOCATE(nUrf(maxord))
    ALLOCATE(nbr(nplant))
    ALLOCATE(nrec(nplant))
    ALLOCATE(ngrow(nplant))
    ALLOCATE(naxemg(nplant))
    ALLOCATE(naxes(nplant))
    ALLOCATE(norder(nplant))
    ALLOCATE(diffnum(nplant))
    ALLOCATE(naxtot(nplant))
    ALLOCATE(nBigLat(nplant))
    ALLOCATE(br_rec(maxgrw))
    ALLOCATE(ordseg(maxrec,nplant))
    ALLOCATE(ibrseg(0:maxrec,nplant))
    ALLOCATE(irecpr(maxrec,nplant))
    ALLOCATE(iaxis(maxgrw,nplant))
    ALLOCATE(num_seg(maxgrw,nplant))
    ALLOCATE(irecsg(maxgrw,nplant))
    ALLOCATE(nestbl(maxgrw,nplant))
    ALLOCATE(ordgrw(maxgrw,nplant))
    ALLOCATE(ibrgrw(maxgrw,nplant))
    ALLOCATE(nvch(maxord+1,nplant))
    ALLOCATE(nMPLch(maxord,nplant))
    ALLOCATE(nnewax(maxemg,nplant))
    ALLOCATE(xplant(nplant))
    ALLOCATE(yplant(nplant))
    ALLOCATE(condMP(nplant))
    ALLOCATE(vol_root(nplant))
    ALLOCATE(mroot(nplant))
    ALLOCATE(mshoot(nplant))
    ALLOCATE(tot_len(nplant))
    ALLOCATE(segdiam(maxrec))
    ALLOCATE(segSoluteMass(maxrec))
    ALLOCATE(segconc(maxrec))
    ALLOCATE(age(maxord,mxpnts))
    ALLOCATE(Urf(maxord,mxpnts))
    ALLOCATE(timest(maxgrw,maxest))
    ALLOCATE(tnewax(maxemg,nplant))
    ALLOCATE(inaxs(maxemg,nplant))
    ALLOCATE(f_rad(maxord,nplant))
    ALLOCATE(strsen(maxord,nplant))
    ALLOCATE(rdmang(maxord,nplant))
    ALLOCATE(brlmax(maxord+1,nplant))
    ALLOCATE(maxlast(maxord,nplant))
    ALLOCATE(numlast(maxord,nplant))
    ALLOCATE(xs(maxrec,nplant))
    ALLOCATE(ys(maxrec,nplant))
    ALLOCATE(zs(maxrec,nplant))
    ALLOCATE(timorg(maxrec,nplant))
    ALLOCATE(seglen(maxrec,nplant))
    ALLOCATE(segsur(0:maxrec,nplant))!from 0?
    ALLOCATE(segmas(maxrec,nplant))
    ALLOCATE(segrad(maxrec,nplant))
    ALLOCATE(crossSectionSeg(maxrec,nplant))
    ALLOCATE(segvol(maxrec,nplant))
    ALLOCATE(mhorm(maxrec,nplant))
    ALLOCATE(xg(maxgrw,nplant))
    ALLOCATE(yg(maxgrw,nplant))
    ALLOCATE(zg(maxgrw,nplant))
    ALLOCATE(ovrtime(maxgrw,nplant))
    ALLOCATE(brlgth(maxgrw,nplant))
    ALLOCATE(agevch(maxord+1,mxpnts,nplant))
    ALLOCATE(vch(maxord,mxpnts,nplant))
    ALLOCATE(sMPLch(maxord,mxpnts,nplant))
    ALLOCATE(MPLch(maxord,mxpnts,nplant))
    ALLOCATE(brspac(maxord-1,nplant))
    ALLOCATE(brnang(maxord-1,nplant))
    ALLOCATE(dtbrch(maxord-1,nplant))
    ALLOCATE(lb(maxord-1,nplant))
    Allocate (AQPh(20,nplant))
    Allocate (AQPv(20,nplant))
    END SUBROUTINE IniRoot
END MODULE RootData
!************************************************************************
!> data needed for plant
!> root characteristics and Doussan input variables
MODULE PlntData
  USE typedef
  USE ParamData, ONLY: mxBcCh, mxBcCh1, maxplant,maxord
  IMPLICIT NONE
  INTEGER(ap) :: nTpot,ntTpLA,nfTpLA,ncTpLA,ntLA,ntRSR,nscRSR,t_ini,t_dev,t_mid,t_end
  INTEGER(ap) :: ntW,nsfW,nscW,ncnc,ns,nBCr,nsfRSR
  REAL(dp) :: hlim,sf
  REAL(dp) :: h50,p50,p1,p2,CMm,VMax,fk,xin,h0,h1,h2,h3 
  REAL(dp) :: TpLA, TpLA_pot
  REAL(dp) :: a_r,a_1,a_2
  REAL(dp) :: TotSur=0._dp,PHcollar,b_1,b_2 

 INTEGER(ap), ALLOCATABLE,DIMENSION(:):: typeBCr
 REAL(dp), ALLOCATABLE,DIMENSION(:):: tTpot,Tpotc,tTpLA,TpLAc,sfTpLA,fTpLAc,scTpLA,tLA,LAc,cTpLAc,tRSR,RSRc,tBCr,BCroot
 REAL(dp), ALLOCATABLE,DIMENSION(:):: sfRSR,fRSRc,cncp,rscnc,tW,Wc,cRSRc,sfW,fWc,scW,cWc,scRSR,sc,rsc
 REAL(dp), ALLOCATABLE,DIMENSION(:):: Tpot,Tact,LA,SpWgt
 REAL(dp), ALLOCATABLE,DIMENSION(:,:):: a,rootrad

 CONTAINS

  SUBROUTINE IniPlant
    USE RootData, ONLY:nplant
    IMPLICIT NONE

    ALLOCATE(typeBCr(mxBcCh))
    ALLOCATE(tTpot(mxBcCh))
    ALLOCATE(Tpotc(mxBcCh))
    ALLOCATE(tTpLA(mxBcCh))
    ALLOCATE(TpLAc(mxBcCh))
    ALLOCATE(sfTpLA(mxBcCh))
    ALLOCATE(fTpLAc(mxBcCh))
    ALLOCATE(scTpLA(mxBcCh))
    ALLOCATE(tLA(mxBcCh))
    ALLOCATE(LAc(mxBcCh))
    ALLOCATE(cTpLAc(mxBcCh))
    ALLOCATE(tRSR(mxBcCh))
    ALLOCATE(RSRc(mxBcCh))
    ALLOCATE(tBCr(mxBcCh))
    ALLOCATE(BCroot(mxBcCh))
    ALLOCATE(sfRSR(mxBcCh1))
    ALLOCATE(fRSRc(mxBcCh1))
    ALLOCATE(cncp(mxBcCh1))
    ALLOCATE(rscnc(mxBcCh1))
    ALLOCATE(tW(mxBcCh1))
    ALLOCATE(Wc(mxBcCh1))
    ALLOCATE(cRSRc(mxBcCh1))
    ALLOCATE(sfW(mxBcCh1))
    ALLOCATE(fWc(mxBcCh1))
    ALLOCATE(scW(mxBcCh1))
    ALLOCATE(cWc(mxBcCh1))
    ALLOCATE(scRSR(mxBcCh1))
    ALLOCATE(sc(mxBcCh1))
    ALLOCATE(rsc(mxBcCh1))

    ALLOCATE(Tpot(nplant))
    ALLOCATE(Tact(nplant))
    ALLOCATE(LA(nplant))
    ALLOCATE(SpWgt(nplant))
    ALLOCATE(a(maxord,nplant))
    ALLOCATE(rootrad(maxord,nplant))
  END SUBROUTINE IniPlant
END MODULE PlntData
!************************************************************************
!> Data needed for Tardieu and Davies model (environnement)
MODULE EnviData
  USE typedef
  IMPLICIT NONE
  INTEGER(ap) :: nclimaticdata
  REAL(dp) :: Cp
  REAL(dp) :: t_initial=0.
  REAL(dp) :: PPFD_t,VPD_t,Temperature_t,Rn_t
  REAL(dp) :: s_t, rho_t, gamma_t, conversion_t,lambda_t,gs_t,Jxc_t
  REAL(dp) :: Hleaf,Hbundle,Hcell,Vcel
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: time_climate,PPFD,VPD,Temperature,PPFDcum,Rn, Precip, T_pot
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: s_f, rho_f, gamma_f, conversion_f,lambda_f
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: hours,days,ampli,psixylmax,psixylmin
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: datad

CONTAINS

  REAL(dp) FUNCTION gs_fun(ABA,Psi,gsmin,gsalpha,gsbeta,gsdelta)
    REAL(dp),INTENT(IN):: ABA,Psi,gsmin,gsalpha,gsbeta,gsdelta
    gs_fun = gsmin + gsalpha*EXP(gsbeta*ABA*EXP(gsdelta*Psi))
    !print*,gsbeta,ABA,(gsdelta*Psi)
    !print*,gsmin,gsalpha*EXP(gsbeta*ABA*EXP(gsdelta*Psi))
    RETURN
  END FUNCTION gs_fun

! ###Penmann-Monteith Equation as in Tardieu et al., 2015: https://doi.org/10.1093/jxb/erv039
  REAL(dp) FUNCTION PM_fun(gs,Sshaded,s,Rn,rho,Cp,ga,VPD,lambda,gamma,conversion)

    REAL(dp),INTENT(IN):: gs,Sshaded,s,Rn,rho,Cp,ga,VPD,lambda,gamma,conversion
    PM_fun = Sshaded*(s*Rn + rho *Cp*ga*VPD)/(lambda*(s+gamma*(1+ga*conversion/gs)))
    RETURN
  END FUNCTION PM_fun

  REAL(dp) FUNCTION PHcollar_fun(PHeq,Krs,Tact)
    REAL(dp),INTENT(IN)::PHeq
    REAL(dp),INTENT(IN):: Tact,Krs
    PHcollar_fun=PHeq-Tact/Krs
    !print*,PHcollar_fun
  END FUNCTION PHcollar_fun

  REAL(dp) FUNCTION ABA_collar_fun(Hseq,Tact,ABAconstit,ABAa,ABAb)
    REAL(dp),INTENT(IN)::Hseq,ABAconstit,ABAa,ABAb
    REAL(dp),INTENT(IN)::Tact
    ABA_collar_fun=ABAconstit+ ABS(ABAa*Hseq)/(ABAb+ABS(Tact)) 
    !print*,ABA_collar_fun
  END FUNCTION ABA_collar_fun

END MODULE EnviData
!************************************************************************
!> Data needed for Tardieu and Davies model (parameters)
MODULE TardieuData
  USE typedef
  IMPLICIT NONE
  INTEGER(ap) :: count_glob=1
  REAL(dp) :: gsmin,gsalpha,gsbeta,gsdelta
  REAL(dp) :: ampliinit,TaucircadGr,TaucircadGxl,TaucircadGstem,TaucircadGc
  REAL(dp) :: TautranspiGr,TautranspiGxl,TautranspiGstem,TautranspiGc
  REAL(dp) :: TauABAGr,TauABAGxl,TauABAGstem,TauABAGc
  REAL(dp) :: Vres,Vsat,ncap,alphacap
  REAL(dp) :: Grmin,Grmax,Gxlmin,Gxlmax,Gstemmin,Gstemmax,Gcmin,Gcmax,Gxl0,Gstem0,Gc0
  REAL(dp) :: Gc,Gxl,Gstem
  REAL(dp) :: ABAconstit,ABAa,ABAb
  REAL(dp) :: S,shading,ga,Sshaded
  Real(dp) :: ABA_collar=0.,Krs0,Kcomp0
  REAL(dp) :: Krs_transpi,Krs_circad
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: Global_transpi_and_Psi
END MODULE TardieuData

!************************************************************************
!> Data needed for the soil grid 
MODULE GridData
  USE ParamData, ONLY: maxElm
  USE typedef
  IMPLICIT NONE
  INTEGER(ap) ::nPt,nElm,nel,nBand,nBCPts,nex,ney,nez,itMax,itMaxRoot,nx,ny,nz,nexSSF,neySSF
  INTEGER(ap) :: nezSSF,nexRLD,neyRLD,nezRLD,nexRho,neyRho,nezRho,geom
  INTEGER(ap), ALLOCATABLE,DIMENSION (:,:) :: elmnod
  INTEGER(ap), ALLOCATABLE,DIMENSION (:) :: subN,n_neigh
  INTEGER(ap), ALLOCATABLE, DIMENSION (:) :: IADN,IADD
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:) :: IAD,MacroList  
  REAL(dp) :: dxGrid,dyGrid,dzGrid,epslonPH,epslonWC,epslonR,factorRelEps,epslonS,dxSSF
  REAL(dp) :: dySSF,dzSSF,dxRLD,dyRLD,dzRLD,dxRho,dyRho,dzRho
  REAL(dp) ::xCol(1000),yCol(1000),x_cent=0,y_cent=0
  REAL(dp) :: RootSk,checkSink,RootSkOld,rad_cyl=0
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: xgrid,ygrid,zgrid,Axy,Bxy,Dxy,Exy,betaw,betac,Width
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: Vn,VElm,HElm,HElmOld,RLD,RSD,denomSe
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: SSF
  REAL(dp), ALLOCATABLE,DIMENSION(:) :: Wn,sink,csink,sinkOld,sink_cube,sink_cubeOld,csink_cube,betac_cube,betac_cube2,cmass_cube
  REAL(dp), ALLOCATABLE,DIMENSION(:,:) :: Deter
  REAL(dp), ALLOCATABLE,DIMENSION(:,:) :: SSF_mat
  REAL(dp), ALLOCATABLE,DIMENSION(:,:,:) :: B1fact,Ax,Ay,Az
  REAL(dp), ALLOCATABLE,DIMENSION(:,:,:) :: bi,ci,di
  REAL(dp), ALLOCATABLE,DIMENSION(:,:,:,:) :: E
  LOGICAL :: RelEps,continu
  CHARACTER LnUnit*5,TmUnit*5,MsUnit*5,CnUnit*5
  !subelement order of first cube->  subN(iE) = 1
  !subelement order of second  cube->  subN(iE) = 2
  !iL(corner_tetraheda,subelements,subN(iE))
  INTEGER(dp), PARAMETER :: iL(1:4,1:5,1:2) = RESHAPE([1, 4, 7, 3, &
       1, 6, 4, 2, &
       5, 6, 7, 1, &
       7, 6, 8, 4, &
       1, 6, 7, 4, &
       5, 2, 3, 1, &
       5, 6, 8, 2, &
       5, 8, 7, 3, &
       5, 2, 8, 3, &
       2, 8, 3, 4 ],SHAPE(iL))
  INTEGER(ap),SAVE::iadd_temp(maxElm*80)

CONTAINS
  SUBROUTINE IniGrid
    IMPLICIT NONE
    ALLOCATE (xgrid(1:nPt))
    ALLOCATE (ygrid(1:nPt))
    ALLOCATE (zgrid(1:nPt))
    ALLOCATE (Vn(1:nPt))
    ALLOCATE (VElm(1:nElm))
    ALLOCATE (HElm(1:nElm))
    ALLOCATE (HElmOld(1:nElm))
    ALLOCATE (subN(1:nPt))
    ALLOCATE (Axy(1:nPt))
    ALLOCATE (denomSe(1:nPt))
    ALLOCATE (Bxy(1:nPt))
    ALLOCATE (Dxy(1:nPt))
    ALLOCATE (Exy(1:nPt))
    ALLOCATE (betaw(1:nPt))
    ALLOCATE (betac(1:nPt))
    ALLOCATE (sink(1:nPt))
    sink=0._dp
    ALLOCATE (sinkOld(1:nPt))
    sinkOld=0._dp
    ALLOCATE (sink_cube(1:nElm))
    sink_cube=0._dp
    ALLOCATE (sink_cubeOld(1:nElm))
    sink_cubeOld=0._dp
    ALLOCATE (csink_cube(1:nElm))
    csink_cube=0._dp
    ALLOCATE (cmass_cube(1:nElm))
    cmass_cube=0._dp
    ALLOCATE (betac_cube(1:nElm))
    betac_cube=0._dp
    ALLOCATE (betac_cube2(1:nElm))
    betac_cube2=0._dp
    ALLOCATE (csink(1:nPt))
    ALLOCATE (Wn(1:nPt))
    ALLOCATE (Width(1:nPt))
    ALLOCATE (elmnod(1:8,1:nElm))
    ALLOCATE (Deter(1:5,1:nElm))
    ALLOCATE (B1fact(1:4,1:5,1:nElm))
    ALLOCATE (Ax(1:4,1:5,1:nElm))
    ALLOCATE (Ay(1:4,1:5,1:nElm))
    ALLOCATE (Az(1:4,1:5,1:nElm))
    ALLOCATE (bi(1:4,1:5,1:nElm))
    ALLOCATE (ci(1:4,1:5,1:nElm))
    ALLOCATE (di(1:4,1:5,1:nElm))
    ALLOCATE (E(1:4,1:4,1:5,1:nElm))

  END SUBROUTINE IniGrid
END MODULE GridData
!*************************************************************************
!> Data needed for the soil domain grid 
MODULE GridData2
  USE ParamData, ONLY: maxElm
  USE typedef
  IMPLICIT NONE
  INTEGER(ap) ::nPt2,nElm2,nel2,nex2,ney2,nez2,nx2,ny2,nz2
  INTEGER(ap), ALLOCATABLE,DIMENSION (:,:) :: elmnod2
  INTEGER(ap), ALLOCATABLE,DIMENSION (:) :: n_neigh2
  REAL(dp) :: dxGrid2,dyGrid2,dzGrid2
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: xgrid2,ygrid2,zgrid2
  !subelement order of first cube->  subN(iE) = 1
  !subelement order of second  cube->  subN(iE) = 2
  !iL(corner_tetraheda,subelements,subN(iE))

CONTAINS
    SUBROUTINE IniGrid2
    IMPLICIT NONE
    ALLOCATE (xgrid2(1:nPt2))
    ALLOCATE (ygrid2(1:nPt2))
    ALLOCATE (zgrid2(1:nPt2))
    ALLOCATE (elmnod2(1:8,1:nElm2))
    END SUBROUTINE IniGrid2
	
END MODULE GridData2
!*************************************************************************
!> Data needed for the Doussan root water uptake module
MODULE DoussanMat
  USE typedef
  USE paramData, ONLY: mxBcCh1,maxplant,maxemg
  USE GridData
  USE RootData, ONLY: ngrow,nrec,maxrec,maxgrw,lSomma_growth,lRootBox_growth,nplant,lKdrop
  USE SparseMatrix !multiple roots
  IMPLICIT NONE

  ! Doussan matrices
  INTEGER(ap) ::nmax,isubmax=100,solveroot_call=0
  INTEGER(ap) :: nLibr=10000,indexValue=10000,switchcriterion=1,n=50 
  INTEGER(ap) :: count_nodes, stresfun
  INTEGER(ap) :: nBC_irecn(maxemg)
  INTEGER(ap), ALLOCATABLE, DIMENSION (:) :: nBCn,nrecOld
  INTEGER(ap), ALLOCATABLE, DIMENSION (:) :: counter2,BCtp_usr,curr_BCtp
  INTEGER(ap), ALLOCATABLE, DIMENSION (:) :: nBC_iprvn,no_voxels,segPerCube !ija,sa,
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:) :: nKh,nLr,mat_seg,vox_seg
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:) :: numNodes_voxel,voxel_no,nsub
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:,:) :: voxel_node
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:,:,:) :: loc_Q,mat_Q,transroot,transtip
  INTEGER(ap), ALLOCATABLE, DIMENSION (:) ::iro,jco,jao,iao
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:,:) :: cube_i
  REAL(dp) :: Jintot,sinktot
  REAL(dp) :: hx_min,stresval1,stresval2,cavitb,cavitc
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: PHr_sub,PH_root_sub,axialRootFlow,veloRoot
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: k_ave,B,Qd,Qi,Q_bc1,Q_bc2,Q_bc,Lr,Khr,GH,Khr_pot
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: qroot,l_seg,Lr_pot,Inv_c1KxKr,Inv_ciiKrKr
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: Phi_mat,h_mat2
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: delta2,delta2old,sinkRtemp
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: tempSinkR,curr_BCr,BCr_usr
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: PHs,PHr,PHrOld,PHrTemp,PHsTemp,PHo,PHsri,PHbulk
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: sinkR,SinkROld,Joutr,KRhizo,JoutrOld
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: Phs_osmotic !osmotic head
  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: PH_micro2
  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: KhRoot,LrRoot,ageKh,ageLr
  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: w_dis,cent,cp_mean,Intc
  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: w_sub,l_sub,sum_dis, beta_weight
  REAL(dp), ALLOCATABLE, DIMENSION (:)::aij,ao

  LOGICAL :: loop1=.TRUE.,stressBC=.FALSE.,ave,old,oldT,eqDis,ItCrit_root=.TRUE.,once=.TRUE.,lcavit=.false.
  LOGICAL :: moment=.FALSE.,savelast=.FALSE.,switchSolve,tcheck=.FALSE.
  LOGICAL :: tcheck2=.FALSE.,ana_aan
  TYPE(SparseMatrixType), ALLOCATABLE, DIMENSION (:) :: plantmatrix !multiple roots
				   
CONTAINS

  SUBROUTINE IniDou
    USE RootData, ONLY :nplant
    IMPLICIT NONE
    ALLOCATE (nBCn(nplant))
    ALLOCATE (nrecOld(nplant))
    ALLOCATE (nKh(1:3,nplant))
    ALLOCATE (nLr(1:3,nplant))
    ALLOCATE (KhRoot(1:3,mxBcCh1,nplant))
    ALLOCATE (LrRoot(1:3,mxBcCh1,nplant))
    ALLOCATE (ageKh(1:3,mxBcCh1,nplant))
    ALLOCATE (ageLr(1:3,mxBcCh1,nplant))
  END SUBROUTINE IniDou

  SUBROUTINE IniMat
    USE RootData, ONLY :nrec,nrec_m,ngrow,ngrow_m,nplant,lno_Archi
    IMPLICIT NONE
    INTEGER(ap) :: ipl2
       
nrec_m=0 !new
ngrow_m=0 !new
IF (.NOT. lno_Archi) THEN 
 DO ipl2=1,nplant !new
  IF (nrec(ipl2) .GT. nrec_m) THEN !new
  nrec_m=nrec(ipl2) !new
  ENDIF !new
  IF (ngrow(ipl2) .GT. ngrow_m) THEN !new
  ngrow_m=ngrow(ipl2) !new
  ENDIF !new
 ENDDO !new
ENDIF
    
    IF (loop1) THEN
       loop1=.FALSE.
    ELSE
       DEALLOCATE(voxel_no)
       DEALLOCATE(mat_seg)
       DEALLOCATE(vox_seg)
       DEALLOCATE(no_voxels)
       DEALLOCATE(numNodes_voxel)
       DEALLOCATE(voxel_node)
       DEALLOCATE(PH_micro2)
       DEALLOCATE(Phi_mat)
       DEALLOCATE(h_mat2)
       DEALLOCATE(delta2)
       DEALLOCATE(delta2old)
       DEALLOCATE(w_dis)
       DEALLOCATE(cube_i)
       DEALLOCATE(cent)
       DEALLOCATE(cp_mean)
       DEALLOCATE(Intc)
       DEALLOCATE(tempSinkR)
       DEALLOCATE(nsub)
       DEALLOCATE(w_sub)
       DEALLOCATE(beta_weight)
       DEALLOCATE(l_sub)
       DEALLOCATE(l_seg)
       DEALLOCATE(GH)
       DEALLOCATE(Joutr)
       DEALLOCATE(JoutrOld)
       DEALLOCATE(KRhizo)
       DEALLOCATE(veloRoot)
       DEALLOCATE(axialRootFlow)
       DEALLOCATE(sinkR)
       DEALLOCATE(sinkRtemp)
       DEALLOCATE(Lr)
       DEALLOCATE(Lr_pot)
       DEALLOCATE(Khr)
       DEALLOCATE(Khr_pot)
       DEALLOCATE(PHs)
       DEALLOCATE(PHsri)
       DEALLOCATE(PHbulk)
       DEALLOCATE(PHo)
       DEALLOCATE(PHs_osmotic)
       DEALLOCATE(PHr)
       DEALLOCATE(PHrOld)
       DEALLOCATE(SinkROld)
       DEALLOCATE(PHrTemp)
       DEALLOCATE(PHr_sub)
       DEALLOCATE(PHsTemp)
       DEALLOCATE(PH_root_sub)
       DEALLOCATE(loc_Q)
       DEALLOCATE(mat_Q)
       DEALLOCATE(Qi)
       DEALLOCATE(Q_bc)
       DEALLOCATE(Q_bc1)
       DEALLOCATE(Q_bc2)
       DEALLOCATE(Qd)
       DEALLOCATE(sum_dis)
       DEALLOCATE(transroot)
       DEALLOCATE(transtip)
       DEALLOCATE(BCtp_usr)
       DEALLOCATE(curr_BCtp)
       DEALLOCATE(curr_BCr)
       DEALLOCATE(BCR_usr)
       DEALLOCATE(k_ave)
       DEALLOCATE(B)
       DEALLOCATE(plantmatrix)
       DEALLOCATE(segPerCube)
       IF (ave) THEN
          DEALLOCATE(counter2)
       ENDIF
    ENDIF
    IF (ana_aan) THEN
       IF (ave) THEN
          ALLOCATE(counter2(1))
          counter2=0
       ENDIF
       ALLOCATE(PH_micro2(1,1,1)) !1:500 -> length r is lower than 500, so no. of compartments lower
       PH_micro2=0
       ALLOCATE (Phi_mat(1))
       Phi_mat = 0._dp
       ALLOCATE (h_mat2(1))
       h_mat2 = 0._dp
    ELSE
       IF (ave) THEN
          ALLOCATE(counter2(1:Nelm))

          counter2=0
       ENDIF
       ALLOCATE(PH_micro2(0:nrec_m,1:500,1:isubmax)) !1:500 -> length r is lower than 500, so no. of compartments lower
       PH_micro2=0
       ALLOCATE(numNodes_voxel(0:nrec_m,1:isubmax))
       numNodes_voxel=0
       ALLOCATE (Phi_mat(1:nLibr-1))
       Phi_mat = 0._dp
       ALLOCATE (h_mat2(1:nLibr-1))
       h_mat2 = 0._dp
    ENDIF
    IF (((old.and.lKdrop)).OR.(ave) .OR. (eqdis)) THEN
          ALLOCATE (voxel_no(0:nrec_m,1:isubmax))
          voxel_no=0
          ALLOCATE (no_voxels(1:nElm))
          no_voxels=0
          ALLOCATE (voxel_node(1:nElm,1:2,indexValue))
          voxel_node=0
          ALLOCATE(numNodes_voxel(0:nrec_m,1:isubmax))
          numNodes_voxel=0
       ELSE
          ALLOCATE(numNodes_voxel(1,1))
          numNodes_voxel=0
          ALLOCATE (voxel_no(1,1))
          voxel_no=0
          ALLOCATE (no_voxels(1))
          no_voxels=0
          ALLOCATE (voxel_node(1,1,1))
          voxel_node=0
       ENDIF

    ALLOCATE(delta2(1:nrec_m+1))
    delta2=0
    ALLOCATE(delta2old(1:nrec_m+1))
    delta2old=0
    ALLOCATE(tempSinkR(1:nrec_m+1))
    tempSinkR=0
    ALLOCATE (l_seg(0:nrec_m))
    l_seg=0._dp
    ALLOCATE (w_sub(0:nrec_m,1:isubmax,1:nplant))
    w_sub=0._dp
    ALLOCATE (beta_weight(0:nrec_m,1:isubmax,1:nplant))
    beta_weight=0._dp
    ALLOCATE (l_sub(0:nrec_m,1:isubmax,1:nplant))
    l_sub=0._dp
    ALLOCATE (PHr_sub(0:nrec_m,1:isubmax))
    PHr_sub=0._dp
    ALLOCATE (PH_root_sub(0:nrec_m,1:isubmax))
    PH_root_sub=0._dp
    ALLOCATE (GH(0:nrec_m,1:nplant))
    GH=0._dp
    ALLOCATE (sinkR(0:nrec_m,1:nplant))
    sinkR=0._dp
    ALLOCATE (BCtp_usr(1:nplant))
    BCtp_usr=0
    ALLOCATE (curr_BCtp(1:nplant))
    curr_BCtp=0
    ALLOCATE (curr_BCr(1:nplant))
    curr_BCr=0
    ALLOCATE (BCr_usr(1:nplant))
    BCr_usr=0
    ALLOCATE (sinkRtemp(0:nrec_m))
    sinkRtemp=0._dp
    ALLOCATE (Joutr(0:nrec_m,1:nplant))
    Joutr=0._dp
    ALLOCATE (JoutrOld(0:nrec_m,1:nplant))
    JoutrOld=0._dp
    ALLOCATE (veloRoot(0:nrec_m,1:nplant))
    veloRoot=0._dp
    ALLOCATE (axialRootFlow(0:nrec_m,1:nplant))
    axialRootFlow=0._dp
    ALLOCATE (Lr(0:nrec_m,1:nplant))
    Lr=0._dp
    ALLOCATE (KRhizo(0:nrec_m,1:nplant))
    KRhizo=0._dp
    ALLOCATE (Lr_pot(0:nrec_m))
    Lr_pot=0._dp
    ALLOCATE (Khr(0:nrec_m,1:nplant))
    Khr=0._dp
    ALLOCATE (Khr_pot(0:nrec_m,1:nplant))
    Khr_pot=0._dp
    ALLOCATE (PHs(0:nrec_m,1:nplant))
    PHs=0._dp
    ALLOCATE (PHsri(0:nrec_m,1:nplant))
    PHsri=0._dp
    ALLOCATE (PHbulk(0:nrec_m,1:nplant))
    PHbulk=0._dp
    ALLOCATE (PHo(0:nrec_m,1:nplant))
    PHo=0._dp
    ALLOCATE (PHs_osmotic(nElm))
    PHs_osmotic=0._dp
    ALLOCATE (PHr(1:nrec_m+1,1:nplant))
    PHr=0._dp
    ALLOCATE (PHrOld(1:nrec_m+1,1:nplant))
    PHrOld=0._dp
    ALLOCATE (SinkROld(0:nrec_m,1:nplant))
    SinkROld=0._dp
    ALLOCATE (PHrTemp(1:nrec_m+1,1:nplant))
    PHrTemp=0._dp
    ALLOCATE (PHsTemp(0:nrec_m,1:nplant))
    PHsTemp=0._dp
    IF((lSomma_growth).OR.(lRootBox_growth)) THEN
       ALLOCATE (transroot(0:maxrec+maxgrw,1:2,1:isubmax,1:nplant))
       transroot=0
       ALLOCATE (transtip(0:maxgrw,1:2,1:isubmax,1:nplant))
       transtip=0
       ALLOCATE (nsub(0:maxrec+maxgrw,1:nplant))
       nsub=1
    ELSEIF(.NOT. (lSomma_growth.OR.(lRootBox_growth))) THEN
       ALLOCATE (nsub(0:nrec_m+ngrow_m,1:nplant))
       nsub=1
       ALLOCATE (transroot(0:nrec_m+ngrow_m,1:2,1:isubmax,1:nplant))
       transroot=0
       ALLOCATE (transtip(0:ngrow_m,1:2,1:isubmax,1:nplant))
       transtip=0
    END IF
    ALLOCATE (loc_Q(0:nrec_m,1:8,1:isubmax,1:nplant))
    loc_Q=0
    ALLOCATE (mat_Q(0:nrec_m,1:8,1:isubmax,1:nplant))
    mat_Q=0
    ALLOCATE (mat_seg(0:nrec_m,1:nplant))
    mat_seg=0
    ALLOCATE (vox_seg(0:nrec_m,1:nplant))
    vox_seg=0
    ALLOCATE (w_dis(0:nrec_m,1:8,1:isubmax,1:nplant))
    w_dis=0._dp
    ALLOCATE (cube_i(0:nrec_m,1:isubmax,1:nplant))
    cube_i=0
    ALLOCATE (cent(0:nrec_m,1:3,1:isubmax,1:nplant))
    cent=0._dp
    ALLOCATE (cp_mean(0:nrec_m,1:3,1:isubmax,1:nplant))
    cp_mean=0._dp
    ALLOCATE (Intc(0:nrec_m+ngrow_m,1:3,1:isubmax,1:nplant))
    Intc=0._dp
    ALLOCATE (Qi(0:nrec_m,1:nplant))
    Qi=0._dp
    ALLOCATE (Q_bc(0:nrec_m,1:nplant))
    Q_bc=0._dp
    ALLOCATE (Q_bc1(0:nrec_m,1:nplant))
    Q_bc1=0._dp
    ALLOCATE (Q_bc2(0:nrec_m,1:nplant))
    Q_bc2=0._dp
    ALLOCATE (Qd(0:2*nrec_m+1,1:nplant))
    Qd=0._dp
    ALLOCATE (sum_dis(0:nrec_m,1:isubmax,1:nplant))
    sum_dis=0._dp
    ALLOCATE (k_ave(0:nrec_m,1:isubmax))
    k_ave=0._dp
    ALLOCATE (B(0:nrec_m,1:isubmax))
    B=0._dp
    ALLOCATE (plantmatrix(1:nplant))
    ALLOCATE (segPerCube(1:nElm))
    segPerCube = 1
  END SUBROUTINE IniMat
END MODULE DoussanMat

!*******************************************************************
!> Data needed for the solute transport inside the roots
MODULE SoluteRootMat
  USE typedef
  USE ParamData, only:  maxParticle,maxrec,mxBcCh1,mxdpth,maxplant
  USE RootData, ONLY: nrec
  IMPLICIT NONE 

  Type Particle
     INTEGER(ap) :: ID
     INTEGER(ap) :: segNum
     REAL(dp) :: positionOld
     REAL(dp) :: position
     REAL(dp) :: mass
     REAL(dp) :: segLen(maxplant)
     REAL(dp) :: partOrig
     TYPE(Particle), POINTER :: prev
     TYPE(Particle), POINTER :: next
  END type Particle

  INTEGER(ap)  :: totalParticleNum,uptakeorder
  INTEGER(ap)  :: segNumParticle(maxParticle),numPa,ndepth
  INTEGER(ap), ALLOCATABLE, DIMENSION (:):: numfollow
  INTEGER(ap), ALLOCATABLE, DIMENSION(:,:)::irecfollow
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:):: nPerm,nPass
  REAL(dp) :: particelMass(maxParticle),theta_R,rho_R, segsorb(maxrec),timestepfactor
  REAL(dp) :: sorp(2)=0, seg_upt(maxrec),Ddepth(mxdpth),dfactor(mxdpth),decayrate
  REAL(dp), allocatable, DIMENSION (:) :: Perm,retard,frac_pass, segsolv,theta_elm,Tottransfer
  REAL(dp), allocatable, DIMENSION (:) :: ccube,SoilUptake_adv,SoilUptake_diff,SoilSoluteUptake,fact
  REAL(dp), allocatable, DIMENSION (:) :: thFC,fac_wet,fac_temp,fac_dept,fac_deg,temp_elm
  REAL(dp), allocatable, DIMENSION (:,:,:) ::agePr,agePs,PrRoot,PsRoot

  LOGICAL :: loop2=.TRUE.,l_linSorb=.false.,l_freundSorb=.false.,l_degrad=.false.
  TYPE(Particle), POINTER :: firstP, pParticle

CONTAINS
  SUBROUTINE IniSol
    USE RootData, ONLY: nplant
    IMPLICIT NONE
    ALLOCATE(nPerm(1:3,nplant))
    ALLOCATE(nPass(1:3,nplant))
    ALLOCATE(agePr(1:3,mxBcCh1,nplant))
    ALLOCATE(agePs(1:3,mxBcCh1,nplant))
    ALLOCATE(PrRoot(1:3,mxBcCh1,nplant))
    ALLOCATE(PsRoot(1:3,mxBcCh1,nplant))
  END SUBROUTINE IniSol

  SUBROUTINE IniSolute(t)
    USE typedef
    USE RootData, ONLY: nrec_m,timorg,ordseg,nplant
    USE GridData, ONLY: nElm
    IMPLICIT NONE

    INTEGER(ap) :: iage,irecn,typ,ipl
    REAL(dp), INTENT(in) :: t
    REAL(dp) :: segage

    IF (.NOT. ALLOCATED(ccube)) THEN
       ALLOCATE(ccube(nElm))
       ccube = 0._dp
    END IF

    IF (.NOT. ALLOCATED(fact)) THEN
       ALLOCATE(fact(nElm))
       fact = 1._dp
    END IF

    IF (.NOT. ALLOCATED(Tottransfer)) THEN
       ALLOCATE(Tottransfer(1))
       Tottransfer=0._dp
    END IF

    IF (.NOT. ALLOCATED(SoilSoluteUptake)) THEN
       ALLOCATE(SoilSoluteUptake(nElm))
       SoilSoluteUptake = 0._dp
    END IF
    IF (.NOT. ALLOCATED(SoilUptake_adv)) THEN
       ALLOCATE(SoilUptake_adv(nElm))
       SoilUptake_adv = 0._dp
    END IF
    IF (.NOT. ALLOCATED(SoilUptake_diff)) THEN
       ALLOCATE(SoilUptake_diff(nElm))
       SoilUptake_diff = 0._dp
    END IF
    IF (.NOT. ALLOCATED(theta_elm)) THEN
       ALLOCATE(theta_elm(nElm))
       theta_elm = 0._dp
    END IF
    IF(.NOT. ALLOCATED(retard)) THEN
       ALLOCATE(retard(maxrec))
       retard = 1._dp
    END IF
    IF(.NOT. ALLOCATED(segsolv)) THEN !!nrec
       ALLOCATE(segsolv(nrec_m))
       segsolv = 0._dp
    END IF
        If(.Not. Allocated(fac_deg)) Then
       Allocate(fac_deg(nElm))
       fac_deg = 1._dp
    End If
    If (l_degrad) Then
       If(.Not. Allocated(thFC)) Then
          Allocate(thFC(nElm))
          thFC = 0._dp
       End If
       If(.Not. Allocated(temp_elm)) Then
          Allocate(temp_elm(nElm))
          temp_elm = 0._dp
       End If
       If(.Not. Allocated(fac_wet)) Then
          Allocate(fac_wet(nElm))
          fac_wet = 1._dp
       End If
       If(.Not. Allocated(fac_temp)) Then
          Allocate(fac_temp(nElm))
          fac_temp = 1._dp
       End If
       If(.Not. Allocated(fac_dept)) Then
          Allocate(fac_dept(nElm))
          fac_dept = 1._dp
       End If
       If(.Not. Allocated(numfollow)) Then
          Allocate(numfollow(nrec_m))
       End If
       If(.Not. Allocated(irecfollow)) Then
          Allocate(irecfollow(nrec_m,4))
       End If
    End If


    ! root permeabilities 
    IF(ALLOCATED(Perm)) DEALLOCATE(Perm)
    ALLOCATE(Perm(nrec_m))
    Perm = 0._dp
    IF(ALLOCATED(frac_pass)) DEALLOCATE(frac_pass)
    ALLOCATE(frac_pass(nrec_m))
    frac_pass = 0._dp

Do ipl=1,nplant
    !> root permeability matrices
    DO irecn=1,nrec(ipl)
       segage = t-timorg(irecn,ipl)
       typ = ordseg(irecn,ipl)
       iage=1
       DO WHILE (agePr(typ,iage,ipl).le.segage)
          iage=iage+1
       ENDDO
       IF (iage>nPerm(typ,ipl)) THEN
          Perm(irecn)=PrRoot(typ,nPerm(typ,ipl),ipl)
       ELSEIF (iage==1) THEN
          Perm(irecn)=PrRoot(typ,iage,ipl)
       ELSE
          Perm(irecn)=PrRoot(typ,iage-1,ipl)+(PrRoot(typ,iage,ipl)-PrRoot(typ,iage-1,ipl))*(segage-agePr(typ,iage-1,ipl))/&
               (agePr(typ,iage,ipl)-agePr(typ,iage-1,ipl))
       ENDIF
    END DO

    DO irecn=1,nrec(ipl)
       segage = t-timorg(irecn,ipl)
       typ = ordseg(irecn,ipl)
       iage=1
       DO WHILE (agePs(typ,iage,ipl).le.segage)
          iage=iage+1
       ENDDO
       IF (iage>nPass(typ,ipl)) THEN
          frac_pass(irecn)=PsRoot(typ,nPass(typ,ipl),ipl)
       ELSEIF (iage==1) THEN
          frac_pass(irecn)=PsRoot(typ,iage,ipl)
       ELSE
          frac_pass(irecn)=PsRoot(typ,iage-1,ipl)+(PsRoot(typ,iage,ipl)-PsRoot(typ,iage-1,ipl))*(segage-agePs(typ,iage-1,ipl))/&
               (agePs(typ,iage,ipl)-agePs(typ,iage-1,ipl))
       ENDIF
    END DO
   End Do

  END SUBROUTINE IniSolute

END MODULE SoluteRootMat
!*******************************************************************
MODULE NumericalRecipes
  USE Typedef

CONTAINS

  SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err,ipl)
    ! solve inverse matrix with preconditionned biconjugate gradient method      
    ! code taken from Numerical recipes in fortran 77, p. 79
    USE SparseMatrix
    USE DoussanMat, ONLY: plantmatrix
    IMPLICIT NONE

    INTEGER(ap), INTENT(in) :: itmax,itol,n,ipl
    INTEGER(ap), INTENT(out) :: iter
    INTEGER(dp) :: j
    REAL(dp), INTENT(inout):: tol,b(:),x(:)
    REAL(dp), INTENT(out) ::err
    REAL(dp) :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm,EPS
    REAL(dp) :: p(n),pp(n),r(n),rr(n),z(n),zz(n)
    PARAMETER (EPS=1.E-14_dp)

    iter=0     
    CALL SM_multiply(plantmatrix(ipl), x, r)
    Do j=1,n
       r(j)=b(j)-r(j)
       rr(j)=r(j)
    End Do
    znrm=1._dp
    bnrm=0
    IF (itol.EQ.1) THEN
       bnrm=snrm(n,b,itol)
    ELSE IF (itol.EQ.2) THEN
       CALL SM_divide_by_diagonal(plantmatrix(ipl), b,z)
       bnrm=snrm(n,z,itol)
    ELSE IF (itol.EQ.3.OR.itol.EQ.4) THEN
       CALL SM_divide_by_diagonal(plantmatrix(ipl), b,z)
       bnrm=snrm(n,z,itol)
       CALL SM_divide_by_diagonal(plantmatrix(ipl), r,z)
       znrm=snrm(n,z,itol)
    ELSE
       PRINT *,'illegal itol in linbcg'
    ENDIF
    CALL SM_divide_by_diagonal(plantmatrix(ipl), r,z)     

    bkden=0
100 If (iter.Le.itmax) Then
       iter=iter+1
       zm1nrm=znrm
       CALL SM_divide_by_diagonal(plantmatrix(ipl), rr,zz)
       bknum=0._dp
       Do j=1,n
          bknum=bknum+z(j)*rr(j)
       End Do
       IF(iter.EQ.1) THEN
          Do j=1,n
             p(j)=z(j)
             pp(j)=zz(j)
          End Do
       ELSE
          bk=bknum/bkden
          Do j=1,n
             p(j)=bk*p(j)+z(j)
             pp(j)=bk*pp(j)+zz(j)
          End Do
       ENDIF
       bkden=bknum
       ! call atimes(n,p,z,0)
       CALL SM_multiply(plantmatrix(ipl), p, z)
       akden=0._dp
       Do j=1,n
          akden=akden+z(j)*pp(j)
       End Do
       ak=bknum/akden
       !         call atimes(n,pp,zz,1)

       CALL SM_multiply_transpose(plantmatrix(ipl), pp, zz)
       Do j=1,n
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
       End Do
       Call SM_divide_by_diagonal(plantmatrix(ipl), r,z)
       IF(itol.EQ.1.OR.itol.EQ.2)THEN
          znrm=1._dp
          err=snrm(n,r,itol)/bnrm
       ELSE IF(itol.EQ.3.OR.itol.EQ.4)THEN
          znrm=snrm(n,z,itol)
          IF(ABS(zm1nrm-znrm).GT.EPS*znrm) THEN
             dxnrm=ABS(ak)*snrm(n,p,itol)
             err=znrm/ABS(zm1nrm-znrm)*dxnrm
          ELSE
             err=znrm/bnrm
             GOTO 100
          ENDIF
          xnrm=snrm(n,x,itol)
          IF(err.LE.0.5_dp*xnrm) THEN
             err=err/xnrm
          ELSE
             err=znrm/bnrm
             GOTO 100
          ENDIF
       ENDIF
       IF(err.GT.tol) THEN
          GOTO 100
       ENDIF
    ENDIF
    RETURN      
  END SUBROUTINE linbcg
  !  (C) Copr. 1986-92 Numerical Recipes Software '%12'%*Sim+).
  !***************************************************************************
  FUNCTION snrm(n,sx,itol)
    !> compute one or two norms of a vector sx(1:n), as signaled by itol. used by linbcg
    !> from Numerical Recipes in F., p.81      
    INTEGER(ap):: n,itol,i,isamax
    REAL(dp):: sx(:),snrm

    IF (itol.LE.3)THEN
       snrm=0._dp
       Do i=1,n
          snrm=snrm+sx(i)**2
       End Do
       snrm=SQRT(snrm)
    ELSE
       isamax=1
       Do i=1,n
          IF(ABS(sx(i)).GT.ABS(sx(isamax))) isamax=i
       End Do
       snrm=ABS(sx(isamax))
    ENDIF
    RETURN
  END FUNCTION snrm

END MODULE NumericalRecipes
!========================================================================
!> Data needed for time control
MODULE tmctrl
  USE typedef
  IMPLICIT NONE
  INTEGER(ap), PARAMETER :: mxOut=1000,mxProf=1000
  INTEGER(ap) :: nOut,nouProf,nouProbe
  INTEGER(ap) :: kout=0,kouprof=0,kouprobe=0,kbcr
  INTEGER(ap) ::kaxemg=0
  REAL(dp) :: tOut(mxOut),touProf(mxProf),touProbe(mxProf)
  REAL(dp) :: dtroot=1000.,dtMin,dtMax,FacInc,FacDec,tmax,t_begin,dtProf,dtProbe
  REAL(dp) :: dtMaxC=1.E+30_dp,tcallr,told,dtold,tpulse,dtopt
  REAL(dp) :: tcBCr,tProf,tProbe
  LOGICAL :: tlevel,tlevel_soil
END MODULE tmctrl
!*******************************************************************
!> Concentration data
MODULE ConData
  USE typedef
  USE ParamData, ONLY: maxnod,maxplant
  IMPLICIT NONE
  REAL(dp):: impc(maxnod)
  REAL(dp),ALLOCATABLE,DIMENSION(:)::coptma,cmax,cmin,coptmi
  
CONTAINS
  SUBROUTINE IniConc
  USE RootData, ONLY:nplant
  IMPLICIT NONE

  ALLOCATE(coptma(nplant))
  ALLOCATE(cmax(nplant))
  ALLOCATE(cmin(nplant))
  ALLOCATE(coptmi(nplant))
  END SUBROUTINE IniConc
END MODULE ConData
!********************************************************************
!> data for boundary conditions
MODULE BoundData
  USE typedef
  USE ParamData, ONLY: mxBcCh, maxbdr, maxIrrig
  IMPLICIT NONE
  INTEGER(ap) :: iBCPt(maxbdr+maxIrrig),nQbcCh,nIbcCh,nCBnd1,nhbcCh,nCBnd2,qfun,homogene=1
  REAL(dp) :: cBound(2),CBnd2(mxBcCh)
  REAL(dp) :: tQbcCh(mxBcCh),thbcCh(mxBcCh),hbc(mxBcCh),hbc1(mxBcCh),hbc2(mxBcCh),hbc3(mxBcCh),hbc4(mxBcCh),hbc5(mxBcCh),hbc6(mxBcCh),hbc7(mxBcCh),hbc8(mxBcCh),hbc9(mxBcCh)
  REAL(dp) :: tCBnd1(mxBcCh),CBnd1(mxBcCh),tCBnd2(mxBcCh)
  REAL(dp) :: xqmin1,xqmin2,xqmax1,xqmax2
  REAL(dp), DIMENSION (:,:) :: Qbc(mxBcCh,mxBcCh),tIbcCh(maxIrrig,mxBcCh),Ibc(maxIrrig,mxBcCh)
  LOGICAL :: noBCflux=.false.
END MODULE BoundData
!********************************************************************
!> data for temperature
MODULE TempData
  USE ParamData, ONLY: mxtime, mxdpth, maxnod,maxplant
  USE typedef
  IMPLICIT NONE
  INTEGER(ap):: nz_tempS,nt_tempS,nt_tempA,nt_presA,nt_presD
  REAL(dp):: time_S(mxtime),depth(mxdpth),time_TA(mxtime),time_PA(mxtime),time_PD(mxtime),tstart
  REAL(dp):: tem(maxnod),impt(maxnod),T_atm(mxtime),P_atm(mxtime),P_diff(mxtime)
  REAL(dp):: Tatm_usr,Patm_usr,Pdiff_usr
  REAL(dp),ALLOCATABLE,DIMENSION(:)::tempermin,topt,tempermax,trange,tmid,expo
  REAL(dp), DIMENSION (:,:) :: temtim(mxtime,mxdpth)
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: Tdepth,Ttime

CONTAINS
  SUBROUTINE IniTemp
  USE RootData, ONLY:nplant
  IMPLICIT NONE
  ALLOCATE(tempermin(nplant))
  ALLOCATE(topt(nplant))
  ALLOCATE(tempermax(nplant))
  ALLOCATE(trange(nplant))
  ALLOCATE(tmid(nplant))
  ALLOCATE(expo(nplant))
  END SUBROUTINE IniTemp 
END MODULE TempData
!********************************************************************
MODULE CumData
  USE Typedef
  USE Paramdata, ONLY: maxelm, maxnod
  IMPLICIT NONE
  INTEGER(ap) :: ListNE(maxnod)
  REAL(dp) :: cumCh0=0.,cumChr=0.,cumCh1=0.,CumRt=0.
  REAL(dp) :: wCumA=0.,wCumT=0.,cCumA,cCumT,wVolI,cVolI
  REAL(dp) :: VolSink=0._dp, VolQ=0._dp, WatVolOld=0._dp,wBalR
  REAL(dp) :: ConAxx(maxelm),ConAyy(maxelm),ConAzz(maxelm)
  REAL(dp):: ConAxy(maxelm),ConAxz(maxelm),ConAyz(maxelm)
  REAL(dp), DIMENSION(2):: CumQ=(/0.,0./),ChemS=(/0.,0./)
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: Qc, Q
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: WatIn,SolIn
END MODULE CumData
!********************************************************************
!> Data for soil domain
MODULE DomData
  USE TypeDef
  IMPLICIT NONE
  REAL(dp):: xmin,xmax,ymin,ymax,zmin,zmax
  REAL(dp):: xmin2,xmax2,ymin2,ymax2,zmin2,zmax2
END MODULE DomData
!********************************************************************

MODULE GeoData
  USE TYPEDEF
  USE ParamData
  IMPLICIT NONE
  INTEGER(ap),ALLOCATABLE,DIMENSION(:)::nanglt
  INTEGER(ap),ALLOCATABLE,DIMENSION(:,:)::nangax
  REAL(dp),ALLOCATABLE,DIMENSION(:)::geolat
  REAL(dp),ALLOCATABLE,DIMENSION(:,:)::anglat,templt
  REAL(dp),ALLOCATABLE,DIMENSION(:,:)::geoaxs
  REAL(dp),ALLOCATABLE,DIMENSION(:,:)::sg
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:)::tempax,angaxs

CONTAINS
  SUBROUTINE IniGeo
  USE RootData, ONLY: nplant
  IMPLICIT NONE
  ALLOCATE(nanglt(nplant))
  ALLOCATE(nangax(maxemg,nplant))
  ALLOCATE(geolat(nplant))
  ALLOCATE(anglat(mxpnts,nplant))
  ALLOCATE(templt(mxpnts,nplant))
  ALLOCATE(geoaxs(maxemg,nplant))
  ALLOCATE(sg(maxord,nplant))
  ALLOCATE(tempax(maxemg,mxpnts,nplant))
  ALLOCATE(angaxs(maxemg,mxpnts,nplant))
  END SUBROUTINE IniGeo
END MODULE GeoData
!********************************************************************
MODULE MatData
  USE TypeDef
  USE GridData, ONLY :nband,nPt 
  IMPLICIT NONE
  INTEGER(ap) :: NumNZ, time_step=0
  INTEGER(ap), ALLOCATABLE, DIMENSION (:) ::jlu,IROW,JCOL
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: alu
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: A
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: A_dparse
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: As
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: A1
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: B
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: Bs
END MODULE MatData
!********************************************************************
!> output variables for observation probes
MODULE ObsData
  USE TypeDef
  USE Paramdata, ONLY: maxobs
  ! npr: number of probes
  ! Pr: Probe ID
  ! Pt: plane direction (perp to X (1), to Y (2) or to Z (3)) or probe given by the user (4))
  ! CrP: crossing point or number of user nodes (if Pt=4)
  ! VarP: variable (WC=1, PH=2, both=3)
  ! distrP: 1=all, 2=average
  IMPLICIT NONE
  INTEGER(ap) :: npr,nprof
  INTEGER(ap) :: Pl(maxobs),varProf(maxobs)
  INTEGER(ap), ALLOCATABLE, DIMENSION(:):: Pt, nodebyPr,Pr,VarP,distrP,CrP
  INTEGER(ap), ALLOCATABLE, DIMENSION(:,:):: NodePr
  LOGICAL :: ObsOK,profOK

CONTAINS
  SUBROUTINE IniProbes
  IMPLICIT NONE
  ALLOCATE(Pt(npr))
  ALLOCATE(nodebyPr(npr))
  ALLOCATE(Pr(npr))
  ALLOCATE(VarP(npr))
  ALLOCATE(distrP(npr))
  ALLOCATE(CrP(npr))
  ALLOCATE(NodePr(npr,1000)) !1000=max number of node for a given plane
  END SUBROUTINE IniProbes
END MODULE ObsData
!********************************************************************
!> Soil data
MODULE SolData
  USE Typedef
  USE ParamData, ONLY: maxnod, maxIrrig, maxbdr, maxmat
  IMPLICIT NONE
  INTEGER(ap) :: NLevel,nMat
  INTEGER(ap) ::KodCB(maxbdr+maxIrrig)
  INTEGER(ap), ALLOCATABLE, DIMENSION(:)::Kode,Kcell,MatNum,MatNum2
  REAL(dp) :: ChPar(10,maxmat),par(11,maxmat),epsi,Peclet,Courant,ssMaxTab(maxmat)
  REAL(dp) :: PeCr
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: Vx,Vy,Vz,Vx_old,Vy_old,Vz_old
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: conO,con,Cap
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: Fc, Gc
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: Dispxx,Dispyy,Dispzz
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: Dispxy,Dispxz,Dispyz
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: hOld,hTemp,hNew,Conc
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: theta, theta_old
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: SoilSoluteConcentration
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: SoilSoluteMass
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: SoilSoluteMass_new
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: Uptake
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: Cum_uptake 
  LOGICAL :: soiltab=.false.,lTab=.false.,lMacro=.FALSE., lCelia=.true.,lRoot_explicit=.false.
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: l_elmMacro
END MODULE SolData
!********************************************************************
!> Rhizosphere Data
MODULE RhizoData
  USE Typedef
  IMPLICIT NONE
  INTEGER(ap) :: rhizoModel
  REAL(dp) :: bulkPara(6), StaticRhizoPara(2)
  REAL(dp) :: RhizoPara(9)
  REAL(dp) :: RhizoSpatialPara(2)
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: thetaTot, thetaNonEq, thetaNonEqOld
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: cTot_r, Rnorm, tauTht, hEqRhizo 
  LOGICAL :: lRhizo
END MODULE RhizoData
!********************************************************************
!> Data for root growth components influenced by soil strength
MODULE StrData
  USE TypeDef
  USE ParamData, ONLY: maxnod
  IMPLICIT NONE
  REAL (dp)::ssmax,simp,refgrd
  REAL (dp), ALLOCATABLE, DIMENSION(:)::s,condu,imps

CONTAINS
  SUBROUTINE IniStrength
  USE GridData, ONLY:nPt
  USE GridData2, ONLY:nPt2
  USE RootData, ONLY:ltwo_grids
  IMPLICIT NONE
  IF (ltwo_grids) THEN
  nPt=nPt2
  ENDIF
  ALLOCATE(s(nPt))
  ALLOCATE(condu(nPt))
  ALLOCATE(imps(nPt))
  END SUBROUTINE IniStrength
END MODULE StrData
!********************************************************************
MODULE WatFun
  USE Typedef
  Use ParamData, Only : n_MFP
  USE SolData, ONLY :nmat
  USE RhizoData
  IMPLICIT NONE

  INTEGER(ap) :: nTab
  REAL(dp),ALLOCATABLE, DIMENSION (:) :: hcheck,hTab,hTab_MFP,hnewcheck
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: TheTab,ConTab,CapTab,MFPTab,mfpcheck
  REAL(dp) ::alh1,dlh,dexp_MFP, hmax_MFP=-0.000001, hmin_MFP=-500000000
CONTAINS

! initialize matrix sizes for tabulated soil characteristics
  SUBROUTINE IniTab
    IMPLICIT NONE
    ALLOCATE (hTab(1:nTab))
    ALLOCATE (TheTab(1:nTab,1:nmat))
    ALLOCATE (ConTab(1:nTab,1:nmat))
    ALLOCATE (CapTab(1:nTab,1:nmat))
    ALLOCATE (MFPTab(n_MFP,1:nmat))
    ALLOCATE (mfpcheck(n_MFP,1:nmat))
    ALLOCATE (hTab_MFP(1:n_MFP))
    ALLOCATE (hcheck(1:n_MFP))
    ALLOCATE (hnewcheck(1:n_MFP))
  END SUBROUTINE IniTab

! current local theta (volumetric soil water content) values:
  REAL(dp) FUNCTION Fth(h,Parloc)
    REAL(dp) :: thr,ths,THT, hcr
    REAL(dp),INTENT(in) :: h,Parloc(:)
    ! if rhizosphere then thr and ths are from rhizobulk
    IF (.NOT. lRhizo) THEN
       thr=Parloc(2)
       ths=Parloc(3)
       IF (h.LT.0.0_dp) THEN
          THT=FThNrm(h,Parloc)
          Fth=thr+THT*(ths-thr)
       ELSE
          Fth=ths
       ENDIF
       ! In case Rhizosphere is considered   
    ELSE              
       hcr = bulkPara(4)
       thr=bulkPara(1)
       ths=bulkPara(2)
       IF (h .LT. hcr) THEN
          THT = FthNrm(h,Parloc)
          Fth = thr + THT*(ths-thr)
       ELSE
          Fth = ths
       ENDIF
    ENDIF
    RETURN
  END FUNCTION Fth

  ! current local THETA (0...1 normalized theta) values:
  REAL(dp) FUNCTION FThNrm(h,Parloc)
    Use MPIutils, Only: stop_program
    INTEGER(ap) :: MatMod
    REAL(dp) :: a,n,m,a2,n2, m2,w1,w2
    REAL(dp) :: lambda=0, hcr=0
    REAL(dp),INTENT(in) :: h,Parloc(:)

    FthNrm=1.0_dp
    a=Parloc(4)
    n=Parloc(5)
    m=1._dp-1._dp/n
    MatMod=INT(Parloc(1))
    IF (lRhizo) THEN ! Rhizosphere is considered
       lambda = bulkPara(3)
       hcr    = bulkPara(4)
       MatMod = 3          
    ENDIF
    SELECT CASE(MatMod)
    CASE (1) ! case1: Mualem Van Genuchten
       IF (h.LT.0.0_dp) THEN
          FthNrm=(1._dp+(a*ABS(h))**n)**(-m)
       ENDIF
    CASE (2) ! case2: Dual porosity model
       w2=parloc(8)
       a2=parloc(9)
       n2=parloc(10)
       w1=1._dp-w2
       m2=1._dp-1._dp/n2
       IF (h.LT.0.0_dp) THEN
          FthNrm=w1*(1._dp+(a*ABS(h))**n)**(-m) + w2*(1._dp+(a2*ABS(h))**n2)**(-m2)
       ENDIF
    CASE (3) ! Case 3, Brooks and Corey as defined in Rhizo.in
       IF (h .LT. hcr) THEN
          FthNrm = (hcr/h)**lambda
       ENDIF
    CASE DEFAULT
       call stop_program('Sorry, only three soil models possible at the moment. Please check <soil.in>.')
    END SELECT
    RETURN
  END FUNCTION FThNrm

!***********************************************************
!###find soil water content from h values
! (from tabulated data)###

  REAL(dp) FUNCTION Fth_soiltab(h,M)
    REAL(dp),INTENT(in) :: h
    INTEGER(ap),INTENT(in) :: M
    INTEGER(ap) :: iT
    IF (h.LT.0.0_dp) THEN
       iT=1
       DO WHILE ((h.LT.hTab(iT+4)).AND.(iT.LE.nTab-8))
          iT=iT+4
       END DO
       DO WHILE ((h.LT.hTab(iT+1)).AND.(iT.LE.nTab-2))
          iT=iT+1
       END DO
       Fth_soiltab=TheTab(iT,M)+(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))*(h-hTab(iT))
    ELSE
       Fth_soiltab=TheTab(1,M)
    ENDIF
    RETURN
  END FUNCTION Fth_soiltab

!***********************************************************
!###find water potential from a theta values###
  REAL(dp) FUNCTION Fh_from_Th(Th,Parloc)
    Use MPIutils, Only: stop_program
    REAL(dp) :: thr,ths,a,n,m
    REAL(dp),INTENT(in) :: Th,Parloc(:)
    INTEGER(ap) :: MatMod

    Fh_from_Th=0.0_dp
    thr=Parloc(2)
    ths=Parloc(3)
    a=Parloc(4)
    n=Parloc(5)
    m=1._dp-1._dp/n
    MatMod=INT(Parloc(1))

    SELECT CASE(MatMod)
    CASE (1) !> + case1: Mualem Van Genuchten
       IF (Th.LT.ths) THEN
          Fh_from_Th=-(((((Th-thr)/(ths-thr))**(-1.0_dp/m))-1.0_dp)**(1.0_dp/n))/a
       ENDIF
    CASE (2)  !> + case2: Dual porosity model
       call stop_program('Can t currently use Fh_from_Th with dual porosity model')
    CASE (3)  !> + case3: Brooks and Corey
       call stop_program('Can t currently use Fh_from_Th with rhizosphere')
    END SELECT
    RETURN
  END FUNCTION Fh_from_Th
!*********************************************************************
!> ### Finds current local matric flux potential value phi from h ###
  REAL(dp) FUNCTION Fmfp_soiltab(h,M,Parloc)
    REAL(dp),INTENT(in) :: h,Parloc(:)
    INTEGER(ap),INTENT(in) :: M
    INTEGER(ap) :: iT
    REAL (dp):: iTR
	
    IF (h.LT.0.0_dp) THEN
      iTR=(log10(abs(h))-log10(abs(hmax_MFP)))/dexp_MFP
      iT=INT(iTR)
      Fmfp_soiltab=MFPTab(iT,M)+(MFPTab(iT+1,M)-MFPTab(iT,M))/(hTab_MFP(iT+1)-hTab_MFP(iT))*(h-hTab_MFP(iT))
    ELSE
       Fmfp_soiltab=MFPTab(1,M)+h*Parloc(6)
    ENDIF
    ! IF (h.LT.0.0_dp) THEN
       ! iT=1
       ! ! DO WHILE ((h.LT.hTab_MFP(iT+1000)).AND.(iT.LE.n_MFP-2000))
          ! ! iT=iT+1000
       ! ! END DO
       ! DO WHILE ((h.LT.hTab_MFP(iT+100)).AND.(iT.LE.n_MFP-200))
          ! iT=iT+100
       ! END DO
       ! DO WHILE ((h.LT.hTab_MFP(iT+10)).AND.(iT.LE.n_MFP-20))
          ! iT=iT+10
       ! END DO
       ! DO WHILE ((h.LT.hTab_MFP(iT+1)).AND.(iT.LE.n_MFP-2))
          ! iT=iT+1
       ! END DO
       ! Fmfp_soiltab=MFPTab(iT,M)+(MFPTab(iT+1,M)-MFPTab(iT,M))/(hTab_MFP(iT+1)-hTab_MFP(iT))*(h-hTab_MFP(iT))
    ! ELSE
       ! Fmfp_soiltab=MFPTab(1,M)+h*Parloc(6)
    ! ENDIF
    RETURN
  END FUNCTION Fmfp_soiltab
 !*********************************************************************
 !###find water potential from a phi value
  REAL(dp) FUNCTION Fh_from_mfp_soiltab(mfp,M)
    REAL(dp),INTENT(in) :: mfp
    INTEGER(ap) :: iT
    INTEGER(ap),INTENT(in) :: M
	
    ! If (mfp.GT.0.0_dp) THEN
	   ! mfp_temp=MFPTab(:,M)
       ! iT=minloc(abs(mfp_temp-mfp))
       ! if (MFPTab(iT,M).EQ.mfp) Then
          ! Fh_from_mfp_soiltab=hTab_MFP(iT)
       ! elseif
          ! iT=iT-1
          ! Fh_from_mfp_soiltab=hTab_MFP(iT)+(hTab_MFP(iT+1)-hTab_MFP(iT))/(MFPTab(iT+1,M)-MFPTab(iT,M))*(mfp-MFPTab(iT,M))
       ! endif
    ! Else
       ! Fh_from_mfp_soiltab=hmin_MFP
    ! Endif
	
    IF (mfp.GT.0.0_dp) THEN
       iT=1
       ! DO WHILE ((mfp.LT.MFPTab(iT+1000,M)).AND.(iT.LE.n_MFP-2000))
          ! iT=iT+1000
       ! END DO
       ! DO WHILE ((mfp.LT.MFPTab(iT+100,M)).AND.(iT.LE.n_MFP-200))
          ! iT=iT+100
       ! END DO
       DO WHILE ((mfp.LT.MFPTab(iT+10,M)).AND.(iT.LE.n_MFP-20))
          iT=iT+10
       END DO
       DO WHILE ((mfp.LT.MFPTab(iT+1,M)).AND.(iT.LE.n_MFP-2))
          iT=iT+1
       END DO
       Fh_from_mfp_soiltab=hTab_MFP(iT)+(hTab_MFP(iT+1)-hTab_MFP(iT))/(MFPTab(iT+1,M)-MFPTab(iT,M))*(mfp-MFPTab(iT,M))
    ELSE
       Fh_from_mfp_soiltab=hmin_MFP
    ENDIF
    RETURN
  END FUNCTION Fh_from_mfp_soiltab
  
 !********************************************************************
 !>### Finds current local soil conductivity for a given h
 ! (analytical solution) ###
  REAL(dp) FUNCTION FKP(h,Parloc,soilNode) 
    INTEGER(ap) :: MatMod, soilNode
    REAL(dp) ::  Ks,THT,lambda,a,a2,n,n2,m,m2,w1,w2,Sv1,Sv2,rNumer,rDenom
    REAL(dp),INTENT(in) :: h,Parloc(:)
    REAL(dp) :: L, ni,d, thtR,thtS, rhob, rhow,cw,  hcr
    REAL(dp) :: Kb,S, thtM, mu, thtBulk
    LOGICAL :: Ralloc =.TRUE.

    FKP=0
    a=Parloc(4)
    n=Parloc(5)
    Ks=Parloc(6)
    lambda=Parloc(7)
    m=1._dp-1._dp/n
    MatMod=INT(Parloc(1))
    !> FKP in case of rhizosphere model
    IF (lRhizo) MatMod = 3
    SELECT CASE(MatMod)
    CASE (1) ! Mualem Van Genuchten
       IF (h.LT.0.0_dp) THEN
          THT=(1._dp+(a*ABS(h))**n)**(-m)
          FKP=Ks*(THT**lambda)*(1._dp-(1._dp-THT**(1._dp/m))**m)**2
       ELSE
          FKP=Ks
       ENDIF
    CASE (2) !Dual porosity model
       w2=parloc(8)
       a2=parloc(9)
       n2=parloc(10)
       w1=1.0-w2
       m2=1._dp-1._dp/n2
       IF (h.LT.0.0_dp) THEN
          Sv1=(1._dp+(a*ABS(h))**n)**(-m)
          Sv2=(1._dp+(a2*ABS(h))**n2)**(-m2)
          THT=w1*Sv1+w2*Sv2
          rNumer=w1*a*(1._dp-(1._dp-Sv1**(1._dp/m))**m) +&
               w2*a2*(1._dp-(1._dp-Sv2**(1._dp/m2))**m2) 
          rDenom=w1*a+w2*a2
          FKP=Ks*(THT**lambda)*(rNumer/rDenom)**2
       ELSE
          FKP=Ks
       ENDIF
    CASE (3) ! Conductivity for Rhizosphere model
       thtR =   bulkPara(1)
       thtS =   bulkPara(2)
       lambda = bulkPara(3)
       hcr =    bulkPara(4)
       Ks =     bulkPara(5)
       L  =     bulkPara(6)
       ni =     rhizoPara(4)
       d  =     rhizoPara(5)
       rhob =   rhizoPara(8)
       rhow =   rhizoPara(9)
       ! Total Saturation
       S = (thetaTot(soilNode)-thtR)/(thtS-thtR)
       ! Bulk conductivity
       IF (S .LT. 1.) THEN
          Kb = Ks*S**(2./lambda + 2. + L)
          ! mucilage water content thtM
          IF (.NOT. ALLOCATED(Rnorm)) THEN
             Ralloc = .false.
             ALLOCATE (Rnorm(soilNode))
             Rnorm = 1
          ENDIF
          IF (Rnorm(soilNode) .LT. 1.) THEN
             IF (h .LT. hcr) THEN
                thtBulk = (thtS-thtR)*(hcr/h)**lambda + thtR
             ELSE
                thtBulk = thtS
             ENDIF
             thtM = ((thetaTot(soilNode)-thtBulk*Rnorm(soilNode)))/&
                  (1. - Rnorm(soilNode))
             cw = (cTot_r(soilNode)*rhob)/(rhow*thtM)
             mu = 1. + ni*cw**d
             FKP = 1./mu * Kb
          ELSE
             FKP = Kb
          ENDIF
       ELSE
          FKP = Ks
       ENDIF
    END SELECT
    IF (.NOT. Ralloc) THEN
       DEALLOCATE(Rnorm)
       Ralloc = .TRUE.
    ENDIF
    RETURN
  END FUNCTION FKP

 !********************************************************************
 !> ### Finds current local soil water conductivity value from h 
 ! (from tabulated values)### 
  REAL(dp) FUNCTION FKP_soiltab(h,M)
    REAL(dp),INTENT(in) :: h
    INTEGER(ap),INTENT(in) :: M
    INTEGER(ap) :: iT
    IF (h.LT.0.0_dp) THEN
       iT=1
       DO WHILE ((h.LT.hTab(iT+4)).AND.(iT.LE.nTab-8))
          iT=iT+4
       END DO
       DO WHILE ((h.LT.hTab(iT+1)).AND.(iT.LE.nTab-2))
          iT=iT+1
       END DO
       FKP_soiltab=(ConTab(iT,M)+(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))*(h-hTab(iT)))
    ELSE
       FKP_soiltab=ConTab(1,M)
    ENDIF
    RETURN
  END FUNCTION FKP_soiltab

 !********************************************************************
 !> ### Finds current local soil water capacity values from h 
 ! (analytical solution)###
  REAL(dp) FUNCTION FCP(h,Parloc)
    REAL(dp), INTENT(in) :: h,Parloc(:)
    REAL(dp) :: thr,ths,a,n,m,C2a,C2b,m2,w1,W2,a2,n2
    REAL(dp) :: hcr, lambda 
    Integer(ap) :: MatMod

    FCP=0
    MatMod=INT(Parloc(1))
    thr=parloc(2)
    ths=parloc(3)
    a=parloc(4)
    n=parloc(5)
    m=1._dp-1._dp/n
    ! Rhizosphere
    IF (lRhizo) MatMod = 3
    SELECT CASE(MatMod)
    CASE (1) ! Mualem Van Genuchten
       IF (h.LT.0.0_dp) THEN
          FCP=(ths-thr) *(a*n*m*((a*ABS(h))**(n-1._dp)))/((1._dp+(a*ABS(h))**n)**(m+1._dp))
       ENDIF
    CASE (2) !Dual porosity model
       w2=parloc(8)
       a2=parloc(9)
       n2=parloc(10)
       w1=1._dp-w2
       m2=1._dp-1._dp/n2
       IF (h.LT.0.0_dp) THEN
          C2a=(ths-thr) *(a*n*m*((a*ABS(h))**(n-1._dp)))/((1._dp+(a*ABS(h))**n)**(m+1._dp))
          C2b=(ths-thr) *(a2*n2*m2*((a2*ABS(h))**(n2-1._dp)))/((1._dp+(a2*ABS(h))**n2)**(m2+1._dp))
          FCP= w1*C2a+w2*C2b
       ENDIF
    CASE (3) ! Rhizosphere model
       thr =    bulkPara(1)
       ths =    bulkPara(2)
       lambda = bulkPara(3)
       hcr =    bulkPara(4)
       IF (h .LT. hcr) THEN
          FCP = -(ths-thr)*lambda*(hcr/h)**lambda/h
       ENDIF

    END SELECT
    RETURN
  END FUNCTION FCP
  
 !********************************************************************
 !> ### Finds current local soil water capacity values from h 
 ! (from tabulated data)###  
  REAL(dp) FUNCTION FCP_soiltab(h,M)
    REAL(dp),INTENT(in) :: h
    INTEGER(ap),INTENT(in) :: M
    INTEGER(ap) :: iT
    IF (h.LT.0.0_dp) THEN
       iT=1
       DO WHILE ((h.LT.hTab(iT+4)).AND.(iT.LE.nTab-8))
          iT=iT+4
       END DO
       DO WHILE ((h.LT.hTab(iT+1)).AND.(iT.LE.nTab-2))
          iT=iT+1
       END DO
       FCP_soiltab=(CapTab(iT,M)+(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))*(h-hTab(iT)))
    ELSE
       FCP_soiltab=0.0_dp
    ENDIF
    RETURN
  END FUNCTION FCP_soiltab
END MODULE WatFun

!********************************************************************
MODULE disToRoot
  USE typedef
  USE RootData, ONLY: xs, ys, zs, nrec, timorg,nplant
  USE GridData, ONLY: xgrid, ygrid, zgrid, nPt, dxGrid, dyGrid, dzGrid
  USE RhizoData, ONLY: RhizoPara, RhizoSpatialPara, cTot_r, Rnorm
  USE DoussanMat, ONLY: loc_Q
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) ::  rhizoMatrix

CONTAINS
  !> RStat is a Subroutine that calculate Rnorm and Ctot_r for a static root
  !system
  SUBROUTINE RStat(t)
    IMPLICIT NONE
    INTEGER(ap)  ::  i, k, jj, count1,ipl
    INTEGER(ap), allocatable, dimension(:) :: rhizoNodes
    REAL(dp) ::t
    REAL(dp) :: alpha, dist, cTot,  RnormEXP
    ALLOCATE(cTot_r(nPt))
    ALLOCATE(Rnorm(nPt))
    ALLOCATE(rhizoNodes(nPt))
    !> Initalized parameters
    rhizoNodes = 0  
    Rnorm(1:nPt) = 1.
    !> Finding and storing soil nodes near the root
    DO i=0,size(loc_Q(:,1,1,1))-1
       DO jj=1,size(loc_Q(1,1,:,1))
          DO k=1,8
             IF(loc_Q(i,k,jj,1) .GT. 0) THEN
                rhizoNodes(loc_Q(i,k,jj,1)) = loc_Q(i,k,jj,1) ! Rhizo node true
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !> RhizoMatrix contain information of the soil node (:,1), min distance from the root (:,2), and the age of the closest root (:,3)
    count1=0
    DO i=1,size(rhizoNodes)
       IF(rhizoNodes(i) .GT. 0) count1 = count1+1
    ENDDO
    ALLOCATE(rhizoMatrix(count1,3))
    ! Initialized rhizoMatrix !
    rhizoMatrix(:,1) = 0
    rhizoMatrix(:,2) = 10000.
    rhizoMatrix(:,3) = 10000.
    k=1
    DO ipl=1,nplant
    DO i=1,npt
       IF(rhizoNodes(i) .GT. 0) THEN
          rhizoMatrix(k,1) = rhizoNodes(i)
          DO jj=1,nrec(ipl)
             dist = SQRT((xs(jj,ipl)-xgrid(i))**2 +(ys(jj,ipl)-ygrid(i))**2 + (zs(jj,ipl)-zgrid(i))**2)
             IF (dist .LT. rhizoMatrix(k,2)) THEN
                rhizoMatrix(k,2) = dist
                rhizoMatrix(k,3) = timorg(jj,ipl)
             ENDIF
          ENDDO
          k = k +1
       ENDIF
    ENDDO
    ENDDO

    cTot  = RhizoPara(3)
    alpha = RhizoSpatialPara(1)
    RnormEXP = RhizoSpatialPara(2)
    cTot_r = 0.
    !> CTot_r is a vector that contain the concentration of mucilage at each soil
    !node. It is calculated based on the Monod equation and has the Form:
    !Ctot_r = ctot * (timeOrigine/(maxTime/alpha + timeorigine)). 
    !alpha is an empirical parameter. The bigger alpha the concentration
    !approach cTot
    !faster. 

    DO i=1,count1
       cTot_r(INT(rhizoMatrix(i,1))) =cTot*rhizoMatrix(i,3)/(MAXVAL(rhizoMatrix(:,3))/alpha + rhizoMatrix(i,3))
       Rnorm(INT(rhizoMatrix(i,1))) = EXP(-RnormEXP*cTot_r(INT(rhizoMatrix(i,1))))
    ENDDO

    call GrowRhizo(t)
    OPEN (UNIT=15,FILE='out/Rnorm.out',STATUS='UNKNOWN')
    WRITE(15,*) 'soilNode   minDist    age cTot_r, Rnorm'
    DO i=1,size(rhizoMatrix(:,1))
       WRITE(15,*) rhizoMatrix(i,1), rhizoMatrix(i,2), rhizoMatrix(i,3),cTot_r(INT(rhizoMatrix(i,1))),Rnorm(INT(rhizoMatrix(i,1)))
    END DO
    CLOSE(15)

  END SUBROUTINE RStat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !** A subroutine to calculate the concentration and Rnorm for growing root. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  SUBROUTINE GrowRhizo(t)
    USE Typedef
    IMPLICIT NONE

    INTEGER(ap)  :: i
    REAL(dp) :: t, relConc, timeOrigin
    REAL(dp) :: alpha, cTot, RnormEXP

    cTot = RhizoPara(3)
    alpha = RhizoSpatialPara(1)
    RnormEXP = RhizoSpatialPara(2)

    DO i=1,size((rhizoMatrix(:,3)))
       timeOrigin = rhizoMatrix(i,3)
       relConc = relC(t, timeOrigin,alpha)
       cTot_r(INT(rhizoMatrix(i,1))) = relConc * cTot !*rhizoMatrix(i,3)/(MAXVAL(rhizoMatrix(:,3))/alpha + rhizoMatrix(i,3))
       Rnorm(INT(rhizoMatrix(i,1))) = EXP(-RnormEXP*cTot_r(INT(rhizoMatrix(i,1))))
    ENDDO

  END SUBROUTINE GrowRhizo

  FUNCTION relC(t,tOrig,alpha)
    ! A function that calculate the relative concentration of mucilage. 
    ! The function get the current simulation time and the time of origin
    ! and return the relative concentration fot that time. 

    REAL(dp) relC ! Relative concentration of mucilage
    REAL(dp) tDec, t, tOrig, alpha ! Time from which concentration start to decay
    tDec = tOrig + 5.0

    IF (t .LE. tOrig) THEN
       relC = 0.0
    ELSEIF (t .GT. tOrig .AND. t .LT. tDec) THEN
       relC = erf(t)
    ELSE
       relC = erf(t) -erf((t-tDec)/alpha)
    ENDIF

    RETURN
  END FUNCTION relC

END MODULE disToRoot
!********************************************************************
!> A Module to calculate the rhizopshrer static parameters. This is also used
!for the rhizosphere initial condition. 
MODULE RhizoStat
  USE typedef
  USE RhizoData, ONLY: bulkPara, RhizoPara

CONTAINS
  !< Secant Methodi
  !*************************************************
  FUNCTION ERR(tht1, h1)
    REAL(dp) :: omega, rhob, rhow, cTot, beta, hcr, lambda
    REAL(dp) :: tht1, heq, h1, thtModel, ERR
    omega = RhizoPara(1); beta = RhizoPara(2); cTot = RhizoPara(3)
    rhob = RhizoPara(8); rhow = RhizoPara(9)
    lambda = bulkPara(3); hcr = bulkPara(4)
    heq = h1 + omega * (rhob*cTot/(rhow*tht1))**beta
    thtModel = RetCurve(heq,hcr,lambda)
    ERR = thtModel - tht1
    RETURN
  END FUNCTION ERR

  !> RetCurve return theta based on brooks and corey    
  FUNCTION RetCurve(h, hcr, lambda)
    REAL(dp) :: h, S, RetCurve
    REAL(dp) :: hcr, lambda, thtR, thtS
    thtR = bulkPara(1); thtS = bulkPara(2)
    IF (h .LT. hcr) THEN
       S = (hcr/h)**lambda
    ELSE
       S = 1.
    ENDIF
    RetCurve = S*(thtS-thtR) + thtR
    RETURN
  END FUNCTION RetCurve
  !*************************************************

END MODULE RhizoStat
