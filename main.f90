PROGRAM WATER_1
USE IFPORT
USE mt19937
USE PARAM_WATER
USE FORCE_1 ,ONLY : INIT_FORCE , FORCE_LJ_EL, ENERGY_LJ_EL &
                       ,ENERGY_BOND, FORCE_BOND &
                       ,CG_ITER, CG_MAXITER, CG_MINITER, CG_CONV &
                       ,SETPB,GETPROP,MPOUT,CONFOUT &
                       ,FREEING_UP,FORCE_ELEC_THz,FORCE_ELEC_RAMAN
IMPLICIT NONE

TYPE(property) :: prop
TYPE(property) :: proptmp
REAL(KIND=double) :: ntemp

!Quantity
REAL(KIND=double) :: eel, ebond, elj, epol
REAL(KIND=double) :: keng, ktemp, envt, enpt
REAL(KIND=double) :: volume,tmass
REAL(KIND=double) :: pint

!phase
INTEGER :: natom
INTEGER :: nmol
REAL(KIND=double),ALLOCATABLE :: cor(:)
REAL(KIND=double),ALLOCATABLE :: com(:)
REAL(KIND=double),ALLOCATABLE :: ptt(:)
REAL(KIND=double),ALLOCATABLE :: for_s(:)
REAL(KIND=double),ALLOCATABLE :: for_f(:)
REAL(KIND=double),ALLOCATABLE :: rot(:,:,:)

!For d/dt
REAL(KIND=double),ALLOCATABLE :: cor_tmp(:)
REAL(KIND=double),ALLOCATABLE :: ptt_tmp(:)
REAL(KIND=double),ALLOCATABLE :: for_stmp(:)
REAL(KIND=double),ALLOCATABLE :: for_ftmp(:)

!For NEMD
REAL(KIND=double) :: pm
REAL(KIND=double),ALLOCATABLE :: cor_eq0(:)
REAL(KIND=double),ALLOCATABLE :: ptt_eq0(:)
REAL(KIND=double),ALLOCATABLE :: for_eq0(:)

!Chaege, Dipole and Polarization
REAL(KIND=double),ALLOCATABLE :: cp(:)
REAL(KIND=double),ALLOCATABLE :: ci(:)
REAL(KIND=double),ALLOCATABLE :: dipole(:)
REAL(KIND=double),ALLOCATABLE :: dipole_ind(:)
REAL(KIND=double),ALLOCATABLE :: polv(:,:)
REAL(KIND=double),ALLOCATABLE :: polv_perm(:,:)
REAL(KIND=double),ALLOCATABLE :: mutmp(:,:)

REAL(KIND=double) :: dipt(3), dipit(3), polt(9)
REAL(KIND=double) :: diptp(3), dipitp(3), poltp(9)
REAL(KIND=double) :: polpt(9),polptp(9)

REAL(KIND=double) :: qHperm,qOperm
REAL(KIND=double) :: qHind,qOind
REAL(KIND=double) :: qHpav,qHiav
REAL(KIND=double) :: qOpav,qOiav
REAL(KIND=double) :: dppav, dppinst
REAL(KIND=double) :: dipav, dipinst
REAL(KIND=double) :: popav, popinst
REAL(KIND=double) :: polav, polinst
REAL(KIND=double) :: anpopav, anpopinst
REAL(KIND=double) :: anpolav, anpolinst


!For newton eq (6th sympletic)
REAL(KIND=double) :: om(7)

!Bath and Pressure
REAL(KIND=double),ALLOCATABLE :: snvt(:)
REAL(KIND=double),ALLOCATABLE :: pnvt(:)
REAL(KIND=double),ALLOCATABLE :: ggt(:)
REAL(KIND=double)  :: gkbt, pnpt, snpt, gnpt, gnpt_f, gnpt_s
REAL(KIND=double)  :: gkt,gkr,gk1,vp,qt,qt1,qtp,odnf,wpt,kbt
INTEGER :: nchain
REAL(KIND=double),ALLOCATABLE :: snvtm(:,:)
REAL(KIND=double),ALLOCATABLE :: pnvtm(:,:)
REAL(KIND=double),ALLOCATABLE :: psin1(:,:)
REAL(KIND=double),ALLOCATABLE :: psin2(:,:)
REAL(KIND=double),ALLOCATABLE :: mass(:)
REAL(KIND=double) :: Qsin1, Qsin2, gsin

!Condition
INTEGER :: maxstep
INTEGER :: nfst !Multitime scale
REAL(KIND=double) :: dt, edt, dtf
REAL(KIND=double) :: pbox,pbox0
REAL(KIND=double) :: rtemp
REAL(KIND=double) :: pext
REAL(KIND=double) :: cut1 , cut2 
REAL(KIND=double) :: cut(2)
REAL(KIND=double) :: cgconv, seng, sengtmp
REAL(KIND=double) :: ewk,epsrf
INTEGER           :: newh
INTEGER           :: ninit , nout , ntrj ,npol, nttime, ntinit, ncheck
CHARACTER(LEN=5)   :: cjob
CHARACTER(LEN=10)  :: cname
CHARACTER(LEN=256) :: mfile, nmfile, polfile, trjfile, hessfile,dipfile
CHARACTER(LEN=256) :: outfile, trjfjle,checkfile,chafile
LOGICAL :: fmin,feng
LOGICAL :: ftrjmin, fpol, fdpdt,fbtrj
LOGICAL :: EqMD

!+NEMD
CHARACTER(LEN=32)  :: jobtype
CHARACTER(LEN=256) :: neqfile
REAL(KIND=double)  :: Elec
INTEGER :: mstep   
INTEGER :: nneq
CHARACTER(LEN=32)  :: NEMDtype
LOGICAL :: Noneq, NoneqTHz, NoneqRAMAN
LOGICAL :: NoneqTHzchk, NoneqRAMANchk

!Init energy
REAL(KIND=double) :: H0

!Configuration
REAL(KIND=double) :: roh,rhh,th
REAL(KIND=double) :: rohav,rhhav,thav
INTEGER :: nconf

INTEGER :: ierr, iter
INTEGER :: ndip,ncha
INTEGER :: n, m , i , j , k , l , ii, kk, jj
INTEGER :: ix, iy ,iz ,jx, jy ,jz ,kx ,ky ,kz ,lx ,ly ,lz
INTEGER :: iaix, iaiz

REAL(KIND=double) :: tmp
REAL(KIND=double) :: phi,c,s,r11,r13,r33,r31

REAL(KIND=double),PARAMETER  :: kcalm2j   = 1.660540d-21*4.184d0
REAL(KIND=double),PARAMETER  :: kcalmolang2bar = kcalm2j*1d25
REAL(KIND=double), ALLOCATABLE :: gvec(:)
REAL(KIND=double), ALLOCATABLE :: hvec(:)
REAL(KIND=double), ALLOCATABLE :: xivec(:)
REAL(KIND=double), ALLOCATABLE :: pcom(:)
REAL(KIND=double), ALLOCATABLE :: xicom(:)
REAL(KIND=double), ALLOCATABLE :: df(:)
REAL(KIND=double), ALLOCATABLE :: xtmp(:)
INTEGER :: nseed
LOGICAL :: fread
LOGICAL :: exist

REAL(KIND=double) :: ct,ct1,ct2
INTEGER :: time_d,time_h,time_m,time_s,time_ms
INTEGER :: v1(8)

REAL(KIND=double) :: eta
REAL(KIND=double) :: lgam
REAL(KIND=double) :: nu
INTEGER :: sc1

NAMELIST/ctr/ maxstep , dt, rtemp, pbox, pext , nout ,ntrj, npol  &
&              , cjob, ninit, cut1, cut2 &
&              , cname, ewk, newh, nmol                      &
&              , nchain, vp, wpt, ntinit , cgconv, epsrf &
&              , fmin, ftrjmin, fpol,feng,seng , fdpdt, fbtrj ,ncha     &
&              , mfile,nmfile,trjfile,polfile,outfile,dipfile,chafile   &
&              , mstep,Elec,neqfile,nneq,jobtype,pbox0,nfst


CALL DATE_AND_TIME(values=v1)

CALL SYSTEM_CLOCK( sc1 )
CALL init_genrand_real(mod(sc1,10000000))
eta = genrand_real()
iset = 0

INQUIRE(FILE='ctrfile.dat',EXIST=exist)

IF(exist /= .TRUE.) THEN
  WRITE(*,'(A)') "Error : ctrfile not found."
  READ(*,*)
  STOP
ENDIF

OPEN(10,FILE='ctrfile.dat',FORM="FORMATTED")
READ(10,ctr)
CLOSE(10)


nseed = INT(SECNDS(0.0))
!bar -> kcal/mol/ang^3
pext = pext/kcalmolang2bar


kbt = rtemp * 1.98624d-3
lgam = 0.005d0 * vp
nu = 0.005d0 * vp

Qsin1 = kbt*10d0**2/vp**2
Qsin2 = kbt*10d0**2/vp**2
gsin = 0.1d0

natom = nmol * 3
vp = vp * 0.00284968d0
ndip = 0

dtf = dt / dble(nfst) 

rohav = 0d0
rhhav = 0d0
thav  = 0d0
nconf = 0

tmass = (mass_h*dble(nmol*2) + mass_o*dble(nmol))*1.6605402d-24/masscoef

IF(cjob == 'NVT'.or.cjob == 'BANAB')THEN
  gkt = dble(3 * natom - 3) * rtemp * 1.98624d-3
ELSE IF(cjob == 'MASSB'.or.cjob == 'MASSN')THEN
  gkt = rtemp * 1.98624d-3
ELSE IF( cjob == 'NPT')THEN
  gkt = dble(3 * natom - 3 - 1) * rtemp * 1.98624d-3
  odnf = 1d0 + 3d0/dble(3*natom - 3)
ENDIF

gk1 = rtemp * 1.98624d-3
qt  = gkt/vp**2
qt1 = gk1/vp**2
qtp = wpt

ALLOCATE(cor(3*natom))
ALLOCATE(com(3*nmol))
ALLOCATE(ptt(3*natom))
ALLOCATE(for_s(3*natom))
ALLOCATE(for_f(3*natom))

ALLOCATE(cor_tmp(3*natom))
ALLOCATE(ptt_tmp(3*natom))
ALLOCATE(for_stmp(3*natom))
ALLOCATE(for_ftmp(3*natom))

ALLOCATE(cor_eq0(3*natom))
ALLOCATE(ptt_eq0(3*natom))
ALLOCATE(for_eq0(3*natom))

ALLOCATE(cp(3*nmol))
ALLOCATE(ci(3*nmol))
ALLOCATE(dipole(3*nmol))
ALLOCATE(dipole_ind(3*nmol))
ALLOCATE(polv(11*nmol,3))
ALLOCATE(polv_perm(11*nmol,3))
ALLOCATE(mutmp(11*nmol,8))

ALLOCATE(pnvt(nchain))
ALLOCATE(snvt(nchain))
ALLOCATE(ggt(nchain))
ALLOCATE(pnvtm(nchain,3*natom))
ALLOCATE(snvtm(nchain,3*natom))
ALLOCATE(psin1(nchain,3*natom))
ALLOCATE(psin2(nchain,3*natom))
ALLOCATE(mass(3*natom))


DO i = 1, natom*3
   IF(mod(i,9)>=1.and.mod(i,9)<=6) THEN
      mass(i) = mass_h
   ELSE
      mass(i) = mass_o
   END IF
END DO


cor = 0d0
ptt = 0d0
ntemp = 0
cut(1) = cut1
cut(2) = cut2
EqMD = .TRUE.

OPEN(9,FILE=outfile,FORM="FORMATTED")
OPEN(10,FILE="energy.dat")

IF(ntrj > 0)THEN
  IF(fbtrj)THEN
  OPEN(15,FILE=trjfile,FORM="UNFORMATTED")
  ELSE
  OPEN(15,FILE=trjfile,FORM="FORMATTED")
  ENDIF
ENDIF

WRITE(9,'(a)') " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"
WRITE(9,'(2a)') " -- JOBNAME : ",cname
WRITE(9,'(a)') " -- MOLECULES : WATER (flexible model)"
WRITE(9,'(a)') " -- POTENTIAL : POLARIZATION MODEL"
WRITE(9,'(a,a)') " -- Job Type : ",jobtype
IF(cjob=='NVT'.or.cjob=='BANAB')WRITE(9,'(a)') " -- NVT algorithm  : NHC"
IF(cjob=='NVT')WRITE(9,'(a)') " -- NHC integrator  : side-NHC"
IF(cjob=='BANAB')WRITE(9,'(a)') " -- NHC integrator  : middle-NHC"
IF(cjob=='MASSB'.or.cjob=='MASSN')WRITE(9,'(a)') " -- NVT algorithm  : massive-NHC"
IF(cjob=='MASSB')WRITE(9,'(a)') " -- NHC integrator  : middle-NHC"
IF(cjob=='MASSN')WRITE(9,'(a)') " -- NHC integrator  : side-NHC"
IF(cjob=='BAOAB'.or.cjob=='OBABO')WRITE(9,'(a)') " -- NVT algorithm  : Langevin"
IF(cjob=='BAOAB')WRITE(9,'(a)') " -- NVT algorithm  : middle-Langevin"
IF(cjob=='OBABO')WRITE(9,'(a)') " -- NVT algorithm  : side-Langevin"
IF(cjob=='Ander'.or.cjob=='oriAn')WRITE(9,'(a)') " -- NVT algorithm  : Andersen"
IF(cjob=='Ander')WRITE(9,'(a)') " -- NVT algorithm  : middle-Andersen"
IF(cjob=='oriAn')WRITE(9,'(a)') " -- NVT algorithm  : end-Andersen"
IF(cjob=='SINR'.or.cjob=='SINRM')WRITE(9,'(a)') " -- NVT algorithm  : SINR"
IF(cjob=='SINR')WRITE(9,'(a)') " -- NVT algorithm  : original-SINR"
IF(cjob=='SINRM')WRITE(9,'(a)') " -- NVT algorithm  : Middle-SINR"
WRITE(9,'(a,i5,i5,a,i2,i5,a1,i2)') " -- Starting Time : " &
 &, v1(1),v1(2),"/",v1(3),v1(5),":",v1(6)
WRITE(9,'(a)') " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"
WRITE(9,'(a,a)') "    NJOB    =         ",cjob
WRITE(9,'(a,i12,a,i12)') "    NATOM   = ",natom,"    NMOL    = ",nmol
WRITE(9,'(a,i12,a,i12)') "    MAXSTEP = ",maxstep,"    NINIT   = ",ninit
WRITE(9,'(a,i12)') "    NOUT    = ",nout
WRITE(9,'(a,i12,a,i12)') "    NPOL    = ",npol,"    NTRJ    = ",ntrj
WRITE(9,'(a,f12.5,a,f12.5)') "    DT(fs)  = ",dt,"    PBOX    = ",pbox
WRITE(9,'(a,f12.5,a,f12.5)') "    CUT1    = ",cut(1),"    CUT2    = ",cut(2)
WRITE(9,'(a,f12.5,a,I12)') "    EWK     = ", ewk, "    NEWH    = ", newh
WRITE(9,'(a,f12.5,a,f12.5)') "    TEMP    = ", rtemp, "    PEXT    = ",pext*kcalmolang2bar
WRITE(9,'(a,i12)') "    NCHAIN   = ", nchain
WRITE(9,'(a,f12.5,a,f12.5,a)') "    VP      = ", vp/0.00284968d0, "    WPT     = ", wpt*1d-8,"(10^8)"
WRITE(9,'(a,f12.5,a,f12.5,a)') "    CGCONV  = ",cgconv,"    EPSRF   = ", epsrf*1d-40,"(10^40)"
WRITE(9,'(a)') " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"

IF(ninit > 0)THEN
  WRITE(9,'(a,a45)') "  MASTER FILE(INOUT) = ",mfile
ENDIF
  WRITE(9,'(a,a45)') "  MASTER FILE(OUT) = ",nmfile

WRITE(9,'(a,a45)') "  OUTPUT FILE     = ",outfile
IF(ntrj > 0)THEN
  WRITE(9,'(a,a45)') "  cor     = ",trjfile
ENDIF
IF(npol > 0)THEN
  WRITE(9,'(a,a45)') "  dip     = ",dipfile
  WRITE(9,'(a,a45)') "  pol     = ",polfile
ENDIF

WRITE(9,'(a)') " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"

fread = .false.

IF(ninit == 1)THEN
  OPEN(11,FILE=mfile,FORM="UNFORMATTED")
  READ(11) cor, ptt, snvt, pnvt, nttime, pbox
  READ(11,IOSTAT=ierr) mutmp(:,:)
  IF(ierr /= 0) mutmp(:,:) = 0d0

  snpt = log(pbox**3)/3d0
  pnpt = 0d0
  snvtm = 0d0
  pnvtm = 0d0
  CLOSE(11)
ELSE IF(ninit == 2)THEN
  OPEN(11,FILE=mfile,FORM="UNFORMATTED")
  READ(11) cor, ptt, snvt, pnvt, nttime, pbox
  READ(11,IOSTAT=ierr) mutmp(:,:)
  IF(ierr /= 0) mutmp(:,:) = 0d0
  
  snvt = 0d0
  snpt = log(pbox**3)/3d0
  pnvt = 0d0
  pnpt = 0d0
  snvtm = 0d0
  pnvtm = 0d0
  CLOSE(11)
ELSE IF(ninit == 3)THEN ! read snpt, pnpt
  OPEN(11,FILE=mfile,FORM="UNFORMATTED")
  READ(11) cor, ptt, snvt, pnvt, nttime, pbox, snpt, pnpt
  READ(11,IOSTAT=ierr) mutmp(:,:)
  IF(ierr /= 0) mutmp(:,:) = 0d0
  READ(11,IOSTAT=ierr) for_f(:)
  fread = .true.
  IF(ierr /= 0) fread = .false.
  CLOSE(11)

  snvtm = 0d0
  pnvtm = 0d0
ELSE IF(ninit == 4)THEN ! read from massive nose-hoover chain
  OPEN(11,FILE=mfile,FORM="UNFORMATTED")
  READ(11) cor, ptt, snvtm, pnvtm, nttime, pbox
  READ(11,IOSTAT=ierr) mutmp(:,:)
  IF(ierr /= 0) mutmp(:,:) = 0d0

  snpt = log(pbox**3)/3d0
  pnpt = 0d0
  CLOSE(11)
ELSE IF(ninit == 5)THEN ! read from SIN(R)
  OPEN(11,FILE=mfile,FORM="UNFORMATTED")
  READ(11) cor, ptt, psin1, psin2, nttime, pbox
  READ(11,IOSTAT=ierr) mutmp(:,:)
  IF(ierr /= 0) mutmp(:,:) = 0d0

  snpt = log(pbox**3)/3d0
  pnpt = 0d0
  snvtm = 0d0
  pnvtm = 0d0
  CLOSE(11)
ELSE
  mutmp(:,:) = 0d0
  snvt = 0d0
  pnvt = 0d0
  snpt = log(pbox**3)/3d0
  pnpt = 0d0
  snvtm = 0d0
  pnvtm = 0d0
  nttime = 0
  n = 0
  
  CALL INIT_COM

  !ninit=0
  CALL SEED(nseed)
  ALLOCATE(rot(3,3,nmol))
  rot(:,:,:) = 0d0
  rot(1,1,:) = 1d0
  rot(2,2,:) = 1d0
  rot(3,3,:) = 1d0

  DO i = 1,nmol
    phi = genrand_real()*pi2
    c = cos(phi)
    s = sin(phi)
    r11= rot(1,1,i)
    r13= rot(1,3,i)
    r31= rot(3,1,i)
    r33= rot(3,3,i)
    rot(1,1,i) = c *r11 - s * r31
    rot(3,1,i) = s *r11 - c * r31
    rot(1,3,i) = c *r13 - s * r33
    rot(3,3,i) = s *r13 + c * r33
  ENDDO

  CALL GET_COR(cor,com,rot)
  DEALLOCATE(rot)
  !ninit=0

ENDIF 


	
IF(ntinit == 0)THEN
  nttime = 0
ENDIF

  EqMD = jobtype(1:3)=="EMD"
  NoneqTHz   = jobtype(1:8)=="NEMD_THz"
  NoneqRAMAN = jobtype(1:10)=="NEMD_RAMAN"
  Noneq    = jobtype(1:4)=="NEMD" !NoneqTHz.or.NoneqRAMAN
  NoneqTHzchk   = jobtype=="NEMD_THz_Chk" &
  & .or.jobtype=="NEMD_THz_Perm_Chk".or.jobtype=="NEMD_THz_Ind_Chk" &
  & .or.jobtype=="NEMD_THz_PERM_Chk".or.jobtype=="NEMD_THz_IND_Chk"
  NoneqRAMANchk = jobtype=="NEMD_RAMAN_Chk" &
  & .or.jobtype=="NEMD_RAMAN_Perm_Chk".or.jobtype=="NEMD_RAMAN_Ind_Chk" &
  & .or.jobtype=="NEMD_RAMAN_PERM_Chk".or.jobtype=="NEMD_RAMAN_IND_Chk" &  
  & .or.jobtype=="NEMD_RAMANxy_Chk"         &
  & .or.jobtype=="NEMD_RAMANxy_Perm_Chk".or.jobtype=="NEMD_RAMANxy_Ind_Chk" &
  & .or.jobtype=="NEMD_RAMANxy_PERM_Chk".or.jobtype=="NEMD_RAMANxy_IND_Chk"
  
  CALL INIT_FORCE(nmol,cut,pbox,ewk,newh,epsrf,fpol)
  tmp = CG_CONV(cgconv)

  IF(fmin == .true.)THEN
  ! Minimize energy
    iter = 1000
    CALL FRPRMN(cor, 3*natom, 1d-5, iter, tmp)
  ENDIF


  for_s(:) = 0d0

  IF(fread == .false.) for_f(:) = 0d0

  IF(ninit==0) THEN
     DO i = 1, 3*natom
        CALL gasdev(eta)
        ptt(i)=sqrt(kbt*mass(i))*eta
     END DO
  END IF
  IF(ninit/=5 .and. (cjob=='SINR'.or.cjob=='SINRM')) THEN
     DO i = 1, 3*natom
        !IF(ptt(i)**2/mass(i)>nchain*kbt*0.81d0) THEN
        !   ptt(i) = sign(sqrt(nchain*kbt*mass(i))*0.9d0,ptt(i))
        !END IF
        DO j = 1, nchain
           !psin1(j,i) = genrand_real()
           CALL gasdev(eta)
           psin1(j,i) = eta
        END DO
        !psin1(:,i) = sqrt(psin1(:,i)/sum(psin1(:,i))*(nchain*kbt-ptt(i)**2/mass(i))*(nchain+1)/nchain/Qsin1)
        !DO j = 1, nchain
        !   eta = genrand_real()
        !   IF(eta<0.5d0) psin1(j,i) = -psin1(j,i)
        !END DO
        DO j = 1, nchain
           CALL gasdev(eta)
           psin2(j,i) = sqrt(kbt/Qsin2)*eta
        END DO
     END DO
  END IF

  CALL GETCOM(cor,com)
  IF(fread == .false.) CALL FORCE_BOND(cor,for_f)
  CALL ENERGY_BOND(cor,ebond)
  IF(fread == .false.) CALL FORCE_LJ_EL(cor,com,for_s,elj,eel,epol,0,mutmp)
  CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
  CALL GETKENG(ptt,keng,ktemp)
  
  !initialize prop
  n = 0
  m = 0    
  CALL INITPROP(prop)
  proptmp = GETPROP()
  prop%pint(1:10) = proptmp%pint(1:10)
  prop%enel = elj
  prop%eel = eel
  prop%epol = epol
  prop%ebond = ebond
  prop%keng = keng
  prop%temp = ktemp
  volume = SETPB(0d0)**3
  prop%pbox = SETPB(0d0)
  prop%density = tmass/(volume*1d-24)

  IF(cjob == 'NVT'.or.cjob == 'BANAB')THEN
    envt = gkt * snvt(1) 
    envt = envt + 0.5d0*pnvt(1)**2*qt
    DO i = 2,nchain
      envt = envt + gk1 *snvt(i) 
      envt = envt + 0.5d0*pnvt(i)**2*qt1
    ENDDO

    CALL OUTPUT(n,nttime,dt,prop,envt)
  ELSE IF(cjob == 'MASSB'.or.cjob == 'MASSN')THEN
    envt = gkt * sum(snvtm(1,:)) 
    envt = envt + 0.5d0*sum(pnvtm(1,:)**2)*qt
    DO i = 2,nchain
      envt = envt + gk1 *sum(snvtm(i,:)) 
      envt = envt + 0.5d0*sum(pnvtm(i,:)**2)*qt1
    ENDDO

    CALL OUTPUT(n,nttime,dt,prop,envt)
  ELSE IF(cjob == 'NPT')THEN
    envt = gkt * snvt(1) 
    envt = envt + 0.5d0*pnvt(1)**2*qt
    DO i = 2,nchain
      envt = envt + gk1 *snvt(i) 
      envt = envt + 0.5d0*pnvt(i)**2*qt1
    ENDDO
    enpt =  0.5d0*pnpt**2*qtp + pext*volume

    CALL OUTPUT(n,nttime,dt,prop,envt,enpt)
  ELSE 

    CALL OUTPUT(n,nttime,dt,prop)
  ENDIF

  tmp = prop%enel + prop%eel + prop%epol + prop%ebond + prop%keng
  WRITE(*,*) "Step: 0"
  WRITE(*,"(a,f15.5,a)") " NVE  =",tmp," kcal/mol"
  WRITE(*,"(a,f15.5,a)") " FORCE=",SQRT(DOT_PRODUCT(for_s,for_s))," kcal/mol/A"
  WRITE(*,*) 
  
  om(1) = 0.78451361047756d0
  om(2) = 0.235573213359357d0
  om(3) = -1.17767998417887d0
  om(4) = 1d0-2d0*(om(1) + om(2) + om(3))
  om(5) = om(3)
  om(6) = om(2)
  om(7) = om(1)
  
  
  IF(NoneqTHz) THEN
    
   IF(jobtype=="NEMD_THz".or.jobtype=="NEMD_THz_Chk")THEN
    NEMDtype = "Total"
   ELSE IF(jobtype(1:13)=="NEMD_THz_Perm".or.jobtype(1:13)=="NEMD_THz_PERM")THEN
    NEMDtype = "Perm"
   ELSE IF(jobtype(1:12)=="NEMD_THz_Ind".or.jobtype(1:12)=="NEMD_THz_IND")THEN
    NEMDtype = "Ind"
   ELSE
    WRITE(*,*) "Error jobtype THz: ",jobtype
    READ(*,*)
    STOP    
   END IF
  
  ELSE IF(NoneqRAMAN) THEN
  
   IF(jobtype=="NEMD_RAMAN".or.jobtype=="NEMD_RAMAN_Chk")THEN
    NEMDtype = "Total"
   ELSE IF(jobtype(1:15)=="NEMD_RAMAN_Perm".or.jobtype(1:15)=="NEMD_RAMAN_PERM")THEN
    NEMDtype = "Perm"
   ELSE IF(jobtype(1:14)=="NEMD_RAMAN_Ind".or.jobtype(1:14)=="NEMD_RAMAN_IND")THEN
    NEMDtype = "Ind"
   ELSE IF(jobtype=="NEMD_RAMANxy".or.jobtype=="NEMD_RAMANxy_Chk")THEN
    NEMDtype = "Totalxy"
   ELSE IF(jobtype(1:17)=="NEMD_RAMANxy_Perm".or.jobtype(1:17)=="NEMD_RAMANxy_PERM")THEN
    NEMDtype = "Permxy"
   ELSE IF(jobtype(1:16)=="NEMD_RAMANxy_Ind".or.jobtype(1:16)=="NEMD_RAMANxy_IND")THEN
    NEMDtype = "Indxy"
   ELSE
    WRITE(*,*) "Error jobtype Raman: ",jobtype
    READ(*,*)
    STOP       
   END IF
  
  ELSE IF(EqMD) THEN
   NEMDtype = "EMD" !No Need
  ELSE IF(jobtype=="ANALYSIS") THEN
   NEMDtype = "ANALYSIS" !No Need
  ELSE 
   WRITE(*,*) "Error jobtype : ",jobtype
   READ(*,*)
   STOP
  END IF
    
  
   OPEN(19,FILE=polfile,FORM="UNFORMATTED")
   OPEN(39,FILE=dipfile,FORM="UNFORMATTED")
   OPEN(59,FILE=chafile,FORM="UNFORMATTED")  
   IF(Noneq)THEN
    OPEN(25,FILE=neqfile,FORM="UNFORMATTED")
   ELSE
    nttime = nttime + 1
   END IF
  
  WRITE(*,*) "LOOP START"

CALL CPU_TIME( ct1 )
  
IF(Noneq)THEN

    !NonEquilibrium Trajectry
  DO m = 1,mstep
	
     READ(25) cor,ptt,pbox
   	
   IF(mod(m,nneq) == 0)THEN
   
   cor_eq0(:) = cor(:)
   ptt_eq0(:) = ptt(:)
   for_eq0(:) = for_s(:)
   
   WRITE(*,*) "mstep",m,"/",mstep
   
   DO l=1,2
   
	IF(l==1)THEN
	 pm = 1d0
	 WRITE(*,*) "p+"
	ELSE
	 pm = -1d0
	 WRITE(*,*) "p-"
	END IF   
    cor(:) = cor_eq0(:)
    ptt(:) = ptt_eq0(:)
	for_s(:) = for_eq0(:)
	
    CALL GETKENG(ptt,keng,ktemp)
	WRITE(*,*) "Initial T (K)",ktemp
	
 	!(t+dh)
	for_f(:) = 0d0
	for_s(:) = 0d0
	CALL GETCOM(cor,com)
	IF(NoneqTHz)THEN
	 CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
	 CALL FORCE_ELEC_THz(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
	ELSE IF(NoneqRAMAN)THEN 
	 CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
	 CALL FORCE_ELEC_RAMAN(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)	
	END IF
	!(t+3dh/2) <= (t+dh/2)
	ptt(:) = ptt(:) + for_s(:)*dt*pm

   	! (t+2dh)
	for_f(:) = 0d0
    for_s(:) = 0d0
    CALL GETCOM(cor,com)
	CALL FORCE_BOND(cor,for_f)
    CALL FORCE_LJ_EL(cor,com,for_s,elj,eel,epol,1,mutmp)

	!(t+2dh = t)
    n = 0
	CALL GET_TRAJECTRY	
	CALL GET_PROPERTY
	WRITE(*,*) "After T (K)",prop%temp

	 DO n = 1,maxstep

	   CALL Solve_Newton(dt)

	   CALL GET_TRAJECTRY
	   
	   IF(mod(n,nout) == 0)THEN
	    CALL GET_PROPERTY
       END IF

	   
     END DO !n
	    
		WRITE(*,*) "Simulation Time(fs)",dble(maxstep)*dt
	    WRITE(*,*) "Final T (K)",prop%temp
		WRITE(*,*) "Final <T>(K)",prop%avgtemp/dble(prop%nsmp)
		
	END DO !l

       CALL GET_TIME
	
	END IF !nneq
	
  END DO !m
  
  
ELSE IF(EqMD)THEN 
 
  !Equilibrium Trajectry
  DO n = 1,maxstep

   IF(cjob /= 'NPT' )THEN

      IF(cjob == 'NVT')THEN

         DO i = 1,7
            CALL NVTCOR(ptt,pnvt,snvt,dt*0.5d0*om(i))
         ENDDO

         CALL Solve_Newton(dt)

         DO i = 1,7
            CALL NVTCOR(ptt,pnvt,snvt,dt*0.5d0*om(i))
         ENDDO

      ELSE IF(cjob == 'MASSN')THEN

         DO j = 1, 3*natom
            DO i = 1,7
               IF(mod(j,9)>=1.and.mod(j,9)<=6)THEN
                  CALL MASSIVE_NVTCOR(ptt(j),pnvtm(:,j),snvtm(:,j),dt*0.5d0*om(i),mass_h)
               ELSE
                  CALL MASSIVE_NVTCOR(ptt(j),pnvtm(:,j),snvtm(:,j),dt*0.5d0*om(i),mass_o)
               END IF
            END DO
         END DO

         CALL Solve_Newton(dt)

         DO j = 1, 3*natom
            DO i = 1,7
               IF(mod(j,9)>=1.and.mod(j,9)<=6)THEN
                  CALL MASSIVE_NVTCOR(ptt(j),pnvtm(:,j),snvtm(:,j),dt*0.5d0*om(i),mass_h)
               ELSE
                  CALL MASSIVE_NVTCOR(ptt(j),pnvtm(:,j),snvtm(:,j),dt*0.5d0*om(i),mass_o)
               END IF
            END DO
         END DO

      ELSE IF(cjob == 'OBABO')THEN

         CALL LangCOR(ptt,dt*0.5d0)

         CALL Solve_Newton(dt)

         CALL LangCOR(ptt,dt*0.5d0)

      ELSE IF(cjob == 'oriAn')THEN

         CALL Solve_Newton(dt)

         CALL Andersen_Thermostat(ptt,dt)

      ELSE IF(cjob == 'BANAB' .or. cjob == 'MASSB' .or. cjob == 'BAOAB' .or. &
      & cjob == 'Ander')THEN

         CALL BATABCOR(dt)

      ELSE IF(cjob == 'SINR') THEN

         CALL SINRCOR(dt)

      ELSE IF(cjob == 'SINRM') THEN

         CALL SINRMCOR(dt)

      END IF
	
   ELSE !NPT

    proptmp = GETPROP()
    prop%pint(1:10) = proptmp%pint(1:10)
    volume = SETPB(0d0)**3
    gnpt = (prop%pint(1)-pext)*(3d0*volume)
    DO i = 1, 7
     CALL NPTCOR(ptt,pnvt,snvt,pnpt,snpt,gnpt,dt*0.5d0*om(i))
    ENDDO

     CALL EVFAST_NPT(dt)

    proptmp = GETPROP()
    prop%pint(1:10) = proptmp%pint(1:10)
    volume = SETPB(0d0)**3
    gnpt = (prop%pint(1)-pext)*(3d0*volume)
    DO i = 1, 7
     CALL NPTCOR(ptt,pnvt,snvt,pnpt,snpt,gnpt,dt*0.5d0*om(i))
    ENDDO

   ENDIF !NPT


   	CALL GET_TRAJECTRY
	
    IF(mod(n,nout) == 0)THEN
     CALL GET_PROPERTY
    ENDIF
	
    nttime = nttime + 1

    IF(feng == .true.)THEN
     sengtmp = prop%enel + prop%eel + prop%epol + prop%ebond + prop%keng
      IF(sengtmp < seng + 1d0 .and. sengtmp > seng -1d0)THEN
        WRITE(9,*) 'ENERGY: ',sengtmp
            EXIT
      ENDIF
    ENDIF

    
   IF(mod(n,10000) == 0)THEN
      pbox = SETPB(0d0)
      OPEN(12,FILE=nmfile,FORM="UNFORMATTED")
    IF(cjob == 'NPT')THEN
      WRITE(12) cor,ptt,snvt,pnvt,nttime,pbox,snpt,pnpt
      WRITE(12) mutmp(:,:)
      WRITE(12) for_f
    ELSE IF(cjob == 'MASSB' .or. cjob == 'MASSN')THEN
      WRITE(12) cor,ptt,snvtm,pnvtm,nttime,pbox
      WRITE(12) mutmp(:,:)
    ELSE IF(cjob == 'SINR'.or.cjob == 'SINRM')THEN
      WRITE(12) cor,ptt,psin1,psin2,nttime,pbox
      WRITE(12) mutmp(:,:)
    ELSE 
      WRITE(12) cor,ptt,snvt,pnvt,nttime,pbox
      WRITE(12) mutmp(:,:)
    ENDIF
      CLOSE(12)
   END IF

   !OUTPUT to Display (time) 
   IF(mod(n,nout)==0)THEN
    CALL GET_TIME
   END IF
   
  ENDDO !maxstep / finished loop 


  !Finished Solving Newton Eq
    pbox = SETPB(0d0)
    OPEN(11,FILE=mfile,FORM="UNFORMATTED")
    IF(cjob == 'NPT')THEN
      WRITE(11) cor,ptt,snvt,pnvt,nttime,pbox,snpt,pnpt
      WRITE(11) mutmp(:,:)
      WRITE(11) for_f
    ELSE IF(cjob == 'MASSB' .or. cjob == 'MASSN')THEN
      WRITE(11) cor,ptt,snvtm,pnvtm,nttime,pbox
      WRITE(11) mutmp(:,:)
    ELSE IF(cjob == 'SINR'.or.cjob == 'SINRM')THEN
      WRITE(11) cor,ptt,psin1,psin2,nttime,pbox
      WRITE(11) mutmp(:,:)
    ELSE 
      WRITE(11) cor,ptt,snvt,pnvt,nttime,pbox
      WRITE(11) mutmp(:,:)
    ENDIF
    CLOSE(11)

 
ELSE IF(jobtype=="ANALYSIS")THEN 

	 DO n = 1,maxstep
	   ntrj = 0
	   IF(fbtrj)THEN
	    READ(15) cor,ptt,pbox
	   ELSE
	    WRITE(*,*) "Error : fbtrj"
		READ(*,*)
		STOP
	   END IF
	   
	    CALL GET_TRAJECTRY	   
	   IF(mod(n,nout) == 0)THEN
	    CALL GET_PROPERTY
	   END IF
	   
	   IF(mod(n,nout) == 0)THEN
	    CALL GET_TIME
	   END IF
     END DO !n

	 
END IF !jobtype
 
 
  !FINISHED
  
  !Check the calculation time
  CALL DATE_AND_TIME(values=v1)
  WRITE(9,*) " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"
  WRITE(9,*) " -- -- CPU Time : "
  WRITE(9,'(i5,a,i5,a,i5,a,i5,a,i5,a)') time_d,"d",time_h,"h",time_m,"m",time_s,"s",time_ms,"ms"
  WRITE(9,*) " -- -- Finish Time : "
  WRITE(9,'(i5,i5,a,i2,i5,a1,i2)') v1(1),v1(2),"/",v1(3),v1(5),":",v1(6)
  WRITE(9,*) " -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"

  CALL FREEING_UP
  WRITE(*,*) 
  WRITE(*,*) "Finished"
  WRITE(*,*) 
  
	
  CONTAINS


  SUBROUTINE GETKENG(ptt,ekn,ktemp)
  IMPLICIT NONE
  REAL(KIND = double),INTENT(IN) :: ptt(3*natom)
  REAL(KIND = double),INTENT(OUT) :: ktemp
  REAL(KIND = double),INTENT(OUT) :: ekn
  REAL(KIND = double) :: keng
  INTEGER :: i,ix,iy,iz
  keng = 0.0d0
  ekn = 0.0d0
  DO i = 1,nmol
    ix = (i - 1)*9 + 1
    iy = ix + 1
    iz = iy + 1
    ekn = ekn + 0.5d0*massi_h*(ptt(ix)**2+ptt(iy)**2+ptt(iz)**2)
    ix = ix + 3
    iy = ix + 1
    iz = iy + 1
    ekn = ekn + 0.5d0*massi_h*(ptt(ix)**2+ptt(iy)**2+ptt(iz)**2)
    ix = ix + 3
    iy = ix + 1
    iz = iy + 1
    ekn = ekn + 0.5d0*massi_o*(ptt(ix)**2+ptt(iy)**2+ptt(iz)**2)
  ENDDO
  ktemp = ekn *2d0 /(dble(3 * natom - 3) * 1.98624d-3)

  END SUBROUTINE


  SUBROUTINE OUTPUT(n,ntt,dt,prop,envt,enpt)
  IMPLICIT NONE
  INTEGER ,INTENT(IN) :: n
  INTEGER ,INTENT(IN) :: ntt
  REAL(KIND = double),INTENT(IN) :: dt
  TYPE(property), INTENT(IN) :: prop
  REAL(KIND = double),INTENT(IN),OPTIONAL :: envt
  REAL(KIND = double),INTENT(IN),OPTIONAL :: enpt
  REAL(KIND = double) :: et
  REAL(KIND = double) :: ep
  REAL(KIND = double) :: nr
  REAL(KIND = double) :: avgeng, avgpot, avgpint,pint
  INTEGER :: i, ix
  REAL(KIND=double) :: ptot(3)
  
  
  nr = dble(prop%nsmp)
  IF(prop%nsmp == 0)THEN
    nr = 1d0
    avgpint = 0d0
  ELSE
    avgpint = (prop%avgpint(1) + &
    &          prop%avgpkint)*kcalmolang2bar
  ENDIF
  
  ep = prop%enel + prop%eel + prop%epol + prop%ebond
  et = ep + prop%keng
  avgpot = prop%avgenel + prop%avgeel + prop%avgepol + prop%avgebond
  avgeng = avgpot + prop%avgkeng
  pint = (prop%pint(1) + 2d0*keng/(3d0*volume))*kcalmolang2bar

!IF(mod(n,10000)==0) THEN
!  write(*,*) 'Pkeng',2d0*keng/(3d0*volume)*kcalmolang2bar
!  write(*,*) 'Pintra',prop%pint(2)*kcalmolang2bar
!  write(*,*) 'Prepulsion',prop%pint(3)*kcalmolang2bar
!  write(*,*) 'Pelec',prop%pint(4)*kcalmolang2bar
!  write(*,*) 'Pelec(ew)',prop%pint(5)*kcalmolang2bar
!  write(*,*) 'PTotal',(2d0*keng/(3d0*volume)+prop%pint(1))*kcalmolang2bar
!END IF
 
 IF(PRESENT(enpt))THEN
  IF(n==0) H0=et+envt+enpt
  WRITE(10,*) dt*dble(n)*0.001d0,dabs(et+envt+enpt-H0)/dabs(H0),(et+envt+enpt)
 ELSE IF(PRESENT(envt))THEN
  IF(n==0) H0=et+envt
  WRITE(10,'(5f)') dt*dble(n)*0.001d0,dabs(et+envt-H0)/dabs(H0),(et+envt),ep,prop%keng
 ELSE IF(cjob=='SINR'.or.cjob=='SINRM')THEN
  IF(n==0) H0=et
  WRITE(10,'(6f)') dt*dble(n)*0.001d0,dabs(et-H0)/dabs(H0),(et),ep,prop%keng,prop%keng*2+dble(nchain)/(nchain+1)*Qsin1*sum(psin1**2)
 ELSE
  IF(n==0) H0=et
  WRITE(10,'(5f)') dt*dble(n)*0.001d0,dabs(et-H0)/dabs(H0),(et),ep,prop%keng
 ENDIF

 
  WRITE(9,*)
  WRITE(9,'(a,I12,a,I10,a)') "---STEP",m,"    (",0,")"
  WRITE(9,'(a,I12,a,I10,a)') "---STEP",n,"    (",ntt,")"
  WRITE(9,'(a,f18.5,a)') "       "  ,dt * n * 0.001d0 , " ps"
  IF(PRESENT(enpt))THEN
    WRITE(9,'(a,f25.15,a)') "--ETOT    =",et+envt+enpt,"  (kcal/mol)"
    WRITE(9,'(a,f12.5,a)') "   -ENPT  =",enpt,"  (kcal/mol)"
    WRITE(9,'(a,f12.5,a)') "   -ENVT  =",envt,"  (kcal/mol)"
    WRITE(9,'(a,f12.5,a,a,f12.5,a)') "   -ENVE  =",et,"  (kcal/mol)" , &
    & "  <ENVE>  =", avgeng/nr , "  (kcal/mol)"
  ELSE IF(PRESENT(envt))THEN
    WRITE(9,'(a,f25.15,a)') "--ETOT    =",et+envt,"  (kcal/mol)"
    WRITE(9,'(a,f12.5,a)') "   -ENVT  =",envt,"  (kcal/mol)"
    WRITE(9,'(a,f12.5,a,a,f12.5,a)') "   -ENVE  =",et,"  (kcal/mol)" , &
    & "  <ENVE>  =", avgeng/nr , "  (kcal/mol)"
  ELSE
    WRITE(9,'(a,f25.15,a)') "--ETOT    =",et,"  (kcal/mol)"
    WRITE(9,'(a,f25.15,a)') " <ETOT>   =",avgeng/nr,"  (kcal/mol)"
  ENDIF
  WRITE(9,*)
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   T     =",prop%temp,'  (K)',&
                                &"         <T>     =",prop%avgtemp/nr ,'  (K)'
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   P     =",pint,'  (bar)',&
                                &"       <P>     =",avgpint/nr ,'  (bar)'
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   Pbox  =",prop%pbox,'  ( A )',&
                                &"       <Pbox>  =",prop%avgpbox/nr ,'  ( A )'
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   EKIN  =",prop%keng,"  (kcal/mol)",&
                                &"  <EKIN>  =",prop%avgkeng/nr, '  (kcal/mol)'
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   EPOT  =",ep,"  (kcal/mol)",&
  &                    "  <EPOT>  =",avgpot/nr ,"  (kcal/mol)"
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   ELJ   =",prop%enel,"  (kcal/mol)",&
                                &"  <ELJ>   =", prop%avgenel/nr,"  (kcal/mol)"
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   ELEC  =",prop%eel,"  (kcal/mol)",&
                                &"  <ELEC>  =", prop%avgeel/nr,"  (kcal/mol)"

  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   EPOL  =",prop%epol,"  (kcal/mol)",&
                                &"  <EPOL>  =",prop%avgepol/nr,"  (kcal/mol)"
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "-   EBOND =",prop%ebond,"  (kcal/mol)",&
                                &"  <EBOND> =",prop%avgebond/nr,"  (kcal/mol)"
  WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- DENSITY =",prop%density,"  (g/cm^3)",&
                            &"   <DENSITY>=",prop%avgdensity/nr,"  (g/cm^3)"
  WRITE(9,'(a,f12.5,a,i5,a,i5)') "  CG iteration : ", CG_ITER() &
  &                  ," max : ", CG_MAXITER(), " min : ", CG_MINITER()


  ptot(1:3) = 0d0
  DO i = 1, natom
    ix = (i-1)*3 + 1
    ptot(1:3) = ptot(1:3) + ptt(ix:ix + 2)
  ENDDO
  WRITE(9,'(a,3ES18.8)')'P ', ptot(1:3)


  IF(EqMD)THEN
  IF(mod(n,nout)==0) THEN
    WRITE(*,*) 
    IF(n==0) WRITE(*,*) "Initial Condition"
    IF(n/=0) WRITE(*,*) "Step: ",n,"/",maxstep
    WRITE(*,"(a,f15.5,a)") " Time =",dt*n*0.001d0," ps"
	WRITE(*,"(a,f15.5,a)") " NVE  =",et," kcal/mol"
	WRITE(*,"(a,f15.5,a15,f15.5)")   " <T>  =",prop%avgtemp/nr," T = ",prop%temp
    WRITE(*,"(a,f15.5,a15,f15.5)")   " <P>  =",avgpint/nr     ," P = ",pint
    WRITE(*,"(a,f15.5,a15,f15.5)")   " <D>  =",prop%avgdensity/nr," D = ",prop%density	
    WRITE(*,"(a,f15.5,a15,f15.5)")   " <L>  =",prop%avgpbox/nr   ," L = ",prop%pbox
	WRITE(*,"(a,f15.5,f15.5,f15.5)") " Ptot =",ptot(1:3)
  END IF
  END IF

  END SUBROUTINE OUTPUT

  SUBROUTINE ZPTT(ptt)
  IMPLICIT NONE
  REAL(KIND=double), INTENT(INOUT) :: ptt(:)
  REAL(KIND=double) :: ptot(3)
  
  ptot(1:3) = 0d0
  DO i = 1, natom
    ix = (i-1)*3 + 1
    ptot(1:3) = ptot(1:3) + ptt(ix:ix + 2)
  ENDDO
  
  WRITE(9,'(a,3ES18.8)')'P ', ptot(1:3)
  WRITE(9,*)
  
  ptot(1:3) = ptot(1:3)/dble(natom)
  DO i = 1, natom
    ix = (i-1)*3 + 1
    ptt(ix:ix + 2) = ptt(ix:ix + 2) - ptot(1:3)
  ENDDO
  
  END SUBROUTINE ZPTT

  SUBROUTINE GETCOM(cor,com)
  IMPLICIT NONE
  REAL(KIND=double),INTENT(IN) :: cor(3*natom)
  REAL(KIND=double),INTENT(OUT) :: com(3*nmol)
  INTEGER :: i , ix , iy , iz , kx , ky , kz , lx , ly , lz ,jx, jy ,jz

   DO i = 1,nmol
    ix = ( i - 1 ) * 3 + 1
    iy = ix + 1
    iz = iy + 1
    jx = 9 * (i - 1) + 1
    jy = jx + 1
    jz = jy + 1
    kx = jx + 3
    ky = kx + 1
    kz = ky + 1
    lx = kx + 3
    ly = lx + 1
    lz = ly + 1

    com(ix) = cor(jx) * mass_h + cor(kx) * mass_h + cor(lx) * mass_o
    com(iy) = cor(jy) * mass_h + cor(ky) * mass_h + cor(ly) * mass_o
    com(iz) = cor(jz) * mass_h + cor(kz) * mass_h + cor(lz) * mass_o
    com(ix) = com(ix) * massi_wat
    com(iy) = com(iy) * massi_wat
    com(iz) = com(iz) * massi_wat
  ENDDO

  END SUBROUTINE GETCOM




  SUBROUTINE NVTCOR(ptt,pnvt,snvt,dtrs)
  REAL(KIND=double),INTENT(INOUT) :: ptt(:)
  REAL(KIND=double),INTENT(INOUT) :: pnvt(:)
  REAL(KIND=double),INTENT(INOUT) :: snvt(:)
  REAL(KIND=double),INTENT(IN) :: dtrs
  REAL(KIND=double) :: pnvt1, etr, ktemp, tmp
  INTEGER :: i

  CALL GETKENG(ptt,etr,ktemp)
  etr = etr * 2d0

  ggt(nchain) = (qt1*pnvt(nchain-1)**2 - gk1)/qt1
  pnvt(nchain) = pnvt(nchain) + 0.5d0 * dtrs * ggt(nchain)
  DO i = 1, nchain - 1


    pnvt(nchain - i) = pnvt(nchain - i) * exp(-0.25d0*dtrs*pnvt(nchain-i+1))

    IF(nchain - i == 1)ggt(1)  =  (etr - gkt)/qt
    IF(nchain - i == 2)ggt(2) = (qt*pnvt(1)**2 - gk1)/qt1
    IF(nchain - i > 2)ggt(nchain-i) = (qt1*pnvt(nchain-i-1)**2 - gk1)/qt1
  
    pnvt(nchain - i) = pnvt(nchain - i) + 0.5d0 * dtrs * ggt(nchain - i)

    pnvt(nchain - i) = pnvt(nchain - i) * exp(-0.25d0*dtrs*pnvt(nchain-i+1))

  ENDDO

  pnvt1  = exp(-dtrs*pnvt(1))
  ptt(:) = ptt(:) * pnvt1

  DO i = 1, nchain
     snvt(i) = snvt(i) + dtrs * pnvt(i)
  ENDDO

  CALL GETKENG(ptt,etr,ktemp)
  etr = etr * 2d0 

  ggt(1)  =  (etr - gkt)/qt

  DO i = 1, nchain - 1

    pnvt(i) = pnvt(i) * exp(-0.25d0*dtrs*pnvt(i+1))

    pnvt(i) = pnvt(i) + 0.5d0 * dtrs * ggt(i)

    pnvt(i) = pnvt(i) * exp(-0.25d0*dtrs*pnvt(i+1))

    IF(i == 1)ggt(2) = (qt*pnvt(1)**2 - gk1)/qt1
    IF(i > 1)ggt(i+1) = (qt1*pnvt(i)**2 - gk1)/qt1
 
  ENDDO

  pnvt(nchain) = pnvt(nchain) + 0.5d0 * dtrs * ggt(nchain)

  END SUBROUTINE NVTCOR
  
  
  SUBROUTINE MASSIVE_NVTCOR(ptt,pnvt,snvt,dtrs,mass)
  REAL(KIND=double),INTENT(INOUT) :: ptt
  REAL(KIND=double),INTENT(INOUT) :: pnvt(:)
  REAL(KIND=double),INTENT(INOUT) :: snvt(:)
  REAL(KIND=double),INTENT(IN) :: dtrs, mass
  REAL(KIND=double) :: etr, pnvt1, ktemp, tmp
  INTEGER :: i

  etr = ptt**2 / mass / 2d0
  etr = etr * 2d0

  ggt(nchain) = (qt1*pnvt(nchain-1)**2 - gk1)/qt1
  pnvt(nchain) = pnvt(nchain) + 0.5d0 * dtrs * ggt(nchain)
  DO i = 1, nchain - 1


    pnvt(nchain - i) = pnvt(nchain - i) * exp(-0.25d0*dtrs*pnvt(nchain-i+1))

    IF(nchain - i == 1)ggt(1)  =  (etr - gkt)/qt
    IF(nchain - i == 2)ggt(2) = (qt*pnvt(1)**2 - gk1)/qt1
    IF(nchain - i > 2)ggt(nchain-i) = (qt1*pnvt(nchain-i-1)**2 - gk1)/qt1
  
    pnvt(nchain - i) = pnvt(nchain - i) + 0.5d0 * dtrs * ggt(nchain - i)

    pnvt(nchain - i) = pnvt(nchain - i) * exp(-0.25d0*dtrs*pnvt(nchain-i+1))

  ENDDO

  pnvt1  = exp(-dtrs*pnvt(1))
  ptt = ptt * pnvt1

  DO i = 1, nchain
     snvt(i) = snvt(i) + dtrs * pnvt(i)
  ENDDO

  etr = ptt**2 / mass / 2d0
  etr = etr * 2d0 

  ggt(1)  =  (etr - gkt)/qt

  DO i = 1, nchain - 1

    pnvt(i) = pnvt(i) * exp(-0.25d0*dtrs*pnvt(i+1))

    pnvt(i) = pnvt(i) + 0.5d0 * dtrs * ggt(i)

    pnvt(i) = pnvt(i) * exp(-0.25d0*dtrs*pnvt(i+1))

    IF(i == 1)ggt(2) = (qt*pnvt(1)**2 - gk1)/qt1
    IF(i > 1)ggt(i+1) = (qt1*pnvt(i)**2 - gk1)/qt1
 
  ENDDO

  pnvt(nchain) = pnvt(nchain) + 0.5d0 * dtrs * ggt(nchain)

  END SUBROUTINE MASSIVE_NVTCOR


  SUBROUTINE NPTCOR(ptt,pnvt,snvt,pnpt,snpt,gnptv,dtrs)
  REAL(KIND=double),INTENT(INOUT) :: ptt(:)
  REAL(KIND=double),INTENT(INOUT) :: pnvt(:)
  REAL(KIND=double),INTENT(INOUT) :: snvt(:)
  REAL(KIND=double),INTENT(INOUT) :: snpt
  REAL(KIND=double),INTENT(INOUT) :: pnpt
  REAL(KIND=double),INTENT(INOUT) :: gnptv
  REAL(KIND=double),INTENT(IN) :: dtrs
  REAL(KIND=double) :: pnvt1, etr, ktemp, tmp
  REAL(KIND=double) :: gnpt
  INTEGER :: i

  CALL GETKENG(ptt,etr,ktemp)
  etr = etr * 2d0

  ggt(nchain) = (qt1*pnvt(nchain-1)**2 - gk1)/qt1
  pnvt(nchain) = pnvt(nchain) + 0.5d0 * dtrs * ggt(nchain)
  DO i = 1, nchain - 1


    pnvt(nchain - i) = pnvt(nchain - i) * exp(-0.25d0*dtrs*pnvt(nchain-i+1))

    IF(nchain - i == 1)ggt(1)  =  (etr + qtp*pnpt**2 - gkt)/qt
    IF(nchain - i == 2)ggt(2) = (qt*pnvt(1)**2 - gk1)/qt1
    IF(nchain - i > 2)ggt(nchain-i) = (qt1*pnvt(nchain-i-1)**2 - gk1)/qt1
  
    pnvt(nchain - i) = pnvt(nchain - i) + 0.5d0 * dtrs * ggt(nchain - i)

    pnvt(nchain - i) = pnvt(nchain - i) * exp(-0.25d0*dtrs*pnvt(nchain-i+1))

  ENDDO

  gnpt = (gnptv + odnf*etr)/qtp
  tmp = EXP(-0.25d0*dtrs*pnvt(1))
  pnpt = pnpt*tmp**2 + 0.5d0*dtrs*tmp*gnpt
! CALL NPT(pnvt,pnpt,gnpt,0.5d0*dtrs)

  pnvt1  = exp(-dtrs*(pnvt(1)+odnf*pnpt))
  ptt(:) = ptt(:) * pnvt1
  DO i = 1, nchain
     snvt(i) = snvt(i) + dtrs * pnvt(i)
  ENDDO
  CALL GETKENG(ptt,etr,ktemp)
  etr = etr * 2d0 

  gnpt = (gnptv + odnf*etr)/qtp
  tmp = EXP(-0.25d0*dtrs*pnvt(1))
  pnpt = pnpt*tmp**2 + 0.5d0*dtrs*tmp*gnpt
! CALL NPT(pnvt,pnpt,gnpt,0.5d0*dtrs)

  ggt(1)  =  (etr + qtp*pnpt**2 - gkt)/qt

  DO i = 1, nchain - 1

    pnvt(i) = pnvt(i) * exp(-0.25d0*dtrs*pnvt(i+1))

    pnvt(i) = pnvt(i) + 0.5d0 * dtrs * ggt(i)

    pnvt(i) = pnvt(i) * exp(-0.25d0*dtrs*pnvt(i+1))

    IF(i == 1)ggt(2) = (qt*pnvt(1)**2 - gk1)/qt1
    IF(i > 1)ggt(i+1) = (qt1*pnvt(i)**2 - gk1)/qt1
 
  ENDDO

  pnvt(nchain) = pnvt(nchain) + 0.5d0 * dtrs * ggt(nchain)



  END SUBROUTINE NPTCOR




  SUBROUTINE GET_COR(cor,com,rot)
  IMPLICIT NONE
  REAL(KIND=double),INTENT(INOUT) :: cor(:)
  REAL(KIND=double),INTENT(INOUT) :: com(:)
  REAL(KIND=double),INTENT(INOUT) :: rot(:,:,:)
  INTEGER :: i , ia , ix , iy ,iz , jx , jy , jz

  DO i = 1,nmol

    ix = ( i - 1 ) * 3 + 1
    iy = ix + 1
    iz = iy + 1


    ia = (i - 1) * 9

!   H1
    jx = ia + 1
    jy = ia + 2
    jz = ia + 3

    cor(jx) = rot(1,1,i)*xh1 +rot(1,2,i)*yh1 + rot(1,3,i)*zh1 +com(ix)
    cor(jy) = rot(2,1,i)*xh1 +rot(2,2,i)*yh1 + rot(2,3,i)*zh1 +com(iy)
    cor(jz) = rot(3,1,i)*xh1 +rot(3,2,i)*yh1 + rot(3,3,i)*zh1 +com(iz)

!   H2
    jx = ia + 4
    jy = ia + 5
    jz = ia + 6


    cor(jx) =  rot(1,1,i)*xh2 + rot(1,2,i)*yh2 + rot(1,3,i)*zh2 + com(ix)
    cor(jy) =  rot(2,1,i)*xh2 + rot(2,2,i)*yh2 + rot(2,3,i)*zh2 + com(iy)
    cor(jz) =  rot(3,1,i)*xh2 + rot(3,2,i)*yh2 + rot(3,3,i)*zh2 + com(iz)

!  O
    jx = ia + 7
    jy = ia + 8
    jz = ia + 9


    cor(jx) =  rot(1,1,i)*xo + rot(1,2,i)*yo + rot(1,3,i)*zo + com(ix)
    cor(jy) =  rot(2,1,i)*xo + rot(2,2,i)*yo + rot(2,3,i)*zo + com(iy)
    cor(jz) =  rot(3,1,i)*xo + rot(3,2,i)*yo + rot(3,3,i)*zo + com(iz)

  ENDDO

  END SUBROUTINE GET_COR


  SUBROUTINE EVFAST(dt)
  IMPLICIT NONE
  REAL(KIND=double),INTENT(IN) :: dt
  REAL(KIND=double) :: dt5
  INTEGER :: i,ix,iy,iz
  
   dt5 = 0.5d0*dt
  
   ptt(:) = ptt(:) + for_f(:)*dt5

   DO i = 1,nmol
    ix = ( i - 1 )*9 + 1
    iy = ix + 1
    iz = ix + 2
    cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dt * massi_h
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dt * massi_h
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dt * massi_o
   ENDDO
  
   for_f(:) = 0d0
   CALL FORCE_BOND(cor,for_f)
   ptt(:) = ptt(:) + for_f(:)*dt5
   
  END SUBROUTINE EVFAST

  
  !MODIFIED STEP
  SUBROUTINE EVFAST_NPT(dt)
  IMPLICIT NONE
  REAL(KIND=double),INTENT(IN) :: dt
  REAL(KIND=double) :: dt5
  REAL(KIND=double), PARAMETER :: e2 = 1d0/6d0
  REAL(KIND=double), PARAMETER :: e4 = e2/20d0
  REAL(KIND=double), PARAMETER :: e6 = e4/42d0
  REAL(KIND=double), PARAMETER :: e8 = e6/72d0
  REAL(KIND=double) :: aa, aa2, arg2, poly, bb
  REAL(KIND=double) :: tmp
  INTEGER :: i,ix,iy,iz
  dt5 = 0.5d0*dt

  ptt(:) = ptt(:) + for_f(:)*dt5


  aa = exp(dt5*pnpt)
  aa2 = aa**2
  arg2 = (pnpt*dt5)**2
  poly = (((e8*arg2 + e6)*arg2+e4)*arg2+e2)*arg2+1d0
  bb = aa*poly*dt

  DO i = 1,nmol
    ix = ( i - 1 )*9 + 1
    iy = ix + 1
    iz = ix + 2
    cor(ix:iz) = cor(ix:iz)*aa2 + ptt(ix:iz) * bb * massi_h 
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    cor(ix:iz) = cor(ix:iz)*aa2 + ptt(ix:iz) * bb * massi_h
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    cor(ix:iz) = cor(ix:iz)*aa2 + ptt(ix:iz) * bb * massi_o
 
  ENDDO

  snpt = snpt + pnpt * dt
  tmp = SETPB(exp(snpt))



  for_f(:) = 0d0
  CALL FORCE_BOND(cor,for_f)
  CALL GETCOM(cor,com)
  CALL FORCE_LJ_EL(cor,com,for_f,elj,eel,epol,1,mutmp)
  ptt(:) = ptt(:) + for_f(:)*dt5
  
  END SUBROUTINE EVFAST_NPT


  

SUBROUTINE FRPRMN(cor, ndim, ftol, iter, fret)
USE PARAM_WATER
  IMPLICIT NONE
  REAL(KIND=double), INTENT(INOUT) :: cor(:)
  INTEGER          , INTENT(IN) :: ndim
  REAL(KIND=double), INTENT(IN) :: ftol
  INTEGER          , INTENT(OUT) :: iter
  REAL(KIND=double), INTENT(OUT) :: fret
  REAL(KIND=double) :: fp, fp1, fp2, fp3, fp4
  REAL(KIND=double) :: gam, dgg, gg
  REAL(KIND=double), PARAMETER :: eps = 1d-10

  ALLOCATE(gvec(ndim))
  ALLOCATE(hvec(ndim))
  ALLOCATE(xivec(ndim))
  ALLOCATE(pcom(ndim))
  ALLOCATE(xicom(ndim))
  ALLOCATE(df(ndim))
  ALLOCATE(xtmp(ndim))

  CALL GETCOM(cor,com)


  CALL ENERGY_LJ_EL(cor,com,fp1,fp2,fp3,mutmp,dipole,dipole_ind,polv,polv_perm)
  CALL ENERGY_BOND(cor,fp4)
 
  fp = fp1 + fp2 + fp3 + fp4
  xivec(:) = 0d0

  CALL FORCE_BOND(cor,xivec)
  CALL FORCE_LJ_EL(cor,com,xivec,fp1,fp2,fp3,0,mutmp)

  gvec(:) = -xivec(:)
  hvec(:) = gvec(:)
  xivec(:) = gvec(:)

  iter = 0
  DO WHILE(.true. )
  iter = iter + 1

    CALL DLINMIN(cor,xivec,ndim,fret)
    IF(mod(iter,1) == 0)THEN
      WRITE(*,*) 'MIN STEP: ', iter
      WRITE(*,*) 'ENERGY  : ',fret
      WRITE(*,*) 'ETOL    : ',fp
    ENDIF

    IF(2d0*ABS(fret - fp) <= ftol*(ABS(fret) + ABS(fp) + eps))EXIT
    
    CALL GETCOM(cor,com)


    CALL ENERGY_LJ_EL(cor,com,fp1,fp2,fp3,mutmp,dipole,dipole_ind,polv,polv_perm)
    CALL ENERGY_BOND(cor,fp4)
    fp = fp1 + fp2 + fp3 + fp4
    xivec(:) = 0d0


    CALL FORCE_BOND(cor,xivec)
    CALL FORCE_LJ_EL(cor,com,xivec,fp1,fp2,fp3,0,mutmp)
    xivec(:) = -xivec(:)

    dgg = 0d0
    gg = 0d0
    gg = DOT_PRODUCT(gvec,gvec)
    dgg = DOT_PRODUCT(gvec+xivec,xivec)
    
    if(gg == 0d0) EXIT
    
    gam = dgg/gg
    gvec(:) = -xivec(:)
    hvec(:) = gvec(:) + gam*hvec(:)
    xivec(:) = hvec(:)

  ENDDO


  DEALLOCATE(pcom)
  DEALLOCATE(xicom)
  DEALLOCATE(gvec)
  DEALLOCATE(hvec)
  DEALLOCATE(xivec)
  DEALLOCATE(df)
  DEALLOCATE(xtmp)



  END SUBROUTINE FRPRMN


    SUBROUTINE DLINMIN(cor, xivec, ndim, fret)
    IMPLICIT NONE
    REAL(KIND=double), INTENT(INOUT) :: cor(:)
    REAL(KIND=double), INTENT(INOUT) :: xivec(:)
    INTEGER          , INTENT(IN)    :: ndim
    REAL(KIND=double), INTENT(INOUT)    :: fret
    REAL(KIND=double) :: xmin, ax, xx, bx, fa, fb, fx
    REAL(KIND=double), PARAMETER :: TOL = 1d-4
    INTEGER :: ncom
  
    ncom = ndim
    pcom(:)  = cor(:)
    xicom(:) = xivec(:)  
    ax = 0d0
    xx = 0.0001d0
    
    CALL MNBRAK(ax,xx,bx,fa,fx,fb)
    CALL DBRENT(ax,xx,bx,TOL,xmin,fret)
    xivec(:) = xivec(:) * xmin
    cor(:) = cor(:) + xivec(:)

    END SUBROUTINE DLINMIN

    SUBROUTINE MNBRAK(ax,bx,cx,fa,fb,fc)
    IMPLICIT NONE
    REAL(KIND=double), INTENT(INOUT) :: ax
    REAL(KIND=double), INTENT(INOUT) :: bx
    REAL(KIND=double), INTENT(INOUT) :: cx
    REAL(KIND=double), INTENT(INOUT) :: fa
    REAL(KIND=double), INTENT(INOUT) :: fb
    REAL(KIND=double), INTENT(INOUT) :: fc
    REAL(KIND=double), PARAMETER :: GOLD = 1.618034d0
    REAL(KIND=double), PARAMETER :: GLIMIT = 100d0
    REAL(KIND=double), PARAMETER :: TINY = 1d-20
    REAL(KIND=double) :: fu, q, r, u, ulim
    REAL(KIND=double) :: dum
  
    fa = F1DIM(ax)
    fb = F1DIM(bx)
    IF( fb > fa)THEN
      dum = ax
      ax = bx
      bx = dum
      dum = fb
      fb = fa
      fa = dum
    ENDIF
    cx = bx+GOLD*(bx - ax)
    fc = F1DIM(cx)



    DO WHILE(fb >= fc)
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx - ((bx-cx)*q-(bx-ax)*r)/(2d0*SIGN(MAX(ABS(q-r),TINY),q-r))
        ulim = bx + GLIMIT*(cx-bx)
        IF((bx-u)*(u-cx) > 0d0)THEN
          fu = F1DIM(u)
          IF(fu < fc)THEN
            ax = bx
            fa = fb
            bx=u
            fb=fu
            RETURN
          ELSE IF(fu > fb)THEN
            cx = u
            fc = fu
            RETURN
          ENDIF
          u = cx + GOLD*(cx-bx)
          fu = F1DIM(u)
        ELSE IF((cx-u)*(u-ulim) > 0d0)THEN
          fu = F1DIM(u)
          IF(fu < fc)THEN
            bx = cx
            cx = u
            u = cx+GOLD*(cx-bx)
            fb = fc
            fc = fu
            fu = F1DIM(u)
          ENDIF
        ELSE IF((u-ulim)*(ulim-cx) >= 0d0)THEN
          u = ulim
          fu = F1DIM(u)
        ELSE
          u = cx+GOLD*(cx-bx)
          fu = F1DIM(u)
        ENDIF
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
    ENDDO
    END SUBROUTINE MNBRAK 

    SUBROUTINE DBRENT(ax, bx, cx, tol, xmin,fret)
    IMPLICIT NONE
    REAL(KIND=double) :: ax
    REAL(KIND=double) :: bx
    REAL(KIND=double) :: cx
    REAL(KIND=double) :: tol
    REAL(KIND=double) :: xmin
    REAL(KIND=double) :: fret
    REAL(KIND=double) :: a, b, d, d1, d2, du, dv, dw, dx, e, fu, fx, fv,fw
    REAL(KIND=double) :: olde, tol1, tol2, u, u1, u2, v, w, x, xm
    INTEGER  :: ITMAX = 100
    INTEGER  :: iter
    REAL(KIND=double) :: ZEPS = 1d-10
    LOGICAL :: ok1, ok2, flag1
    
    a = MIN(ax, cx)
    b = MAX(ax, cx)
    v = bx
    w = v
    x = v
    e = 0d0
    fx = F1DIM(x)
    fv = fx
    fw = fx
    dx = DF1DIM(x)
    dv = dx
    dw = dx
    fu = 0d0  
    u = 0d0
 
    DO iter = 1, ITMAX
      flag1 = .true.
      xm = 0.5d0*(a+b)
      tol1 = tol*ABS(x) + ZEPS
      tol2 = 2d0*tol1
      if(ABS(x-xm) <= (tol2 - 0.5d0*(b-a))) EXIT
      IF(ABS(e) > tol1)THEN

        d1 = 2d0*(b-a)
        d2 = d1
        IF(dw /= dx) d1 = (w-x)*dx/(dx-dw)
        IF(dv /= dx) d2 = (v-x)*dx/(dx-dv)
        u1 = x + d1
        u2 = x + d2
        ok1 = ((a-u1)*(u1-b) > 0) .and. (dx*d1 <= 0d0)
        ok2 = ((a-u2)*(u2-b) > 0) .and. (dx*d2 <= 0d0)
        olde =  e
        e = d
        IF((ok1 .or. ok2))THEN
          IF(ok1 .and. ok2)THEN
            IF(ABS(d1) < ABS(d2))THEN
              d = d1
            ELSE
              d = d2
            ENDIF
          ELSE IF(ok1)THEN
            d = d1
          ELSE 
            d = d2
          ENDIF


          IF(.not. (ABS(d) > ABS(0.5d0*olde)))THEN


            u = x + d
            IF(u-a < tol2 .or. b - u < tol2) d = SIGN(tol1,xm-x)
            flag1  = .false.
          ENDIF
        ENDIF
      ENDIF

        IF(flag1)THEN
          IF(dx >= 0d0)THEN
            e = a - x
          ELSE
            e = b - x
          ENDIF
          d = 0.5d0*e
        ENDIF

        IF(ABS(d) >= tol1)THEN

          u = x + d
          fu = F1DIM(u)
        ELSE

          u = x + SIGN(tol1,d)
          fu = F1DIM(u)
          IF(fu > fx) EXIT
        ENDIF
      
        du = DF1DIM(u)
        IF(fu <= fx)THEN
          IF(u >= x)THEN
            a = x
          ELSE
            b = x
          ENDIF
          v = w 
          fv = fw
          dv = dw
          w = x
          fw = fx
          dw = dx
          x = u
          fx = fu
          dx = du
        ELSE
          IF( u < x)THEN
            a = u
          ELSE
            b = u
          ENDIF
          IF(fu <= fw .or. w == x)THEN
            v = w
            fv =fw
            dv = dw
            w = u 
            fw = fu
            dw = du 
          ELSE IF(fu <= fv .or. v == x .or. v == w)THEN
            v = u
            fv = fu
            dv = du
          ENDIF
        ENDIF
    ENDDO
    IF(iter > 99)THEN
      write(*,*) 'iter',iter
    ENDIF
    xmin = x 
    fret = fx

    END SUBROUTINE DBRENT
      
    



    REAL(KIND=double) FUNCTION F1DIM(x)
    IMPLICIT NONE
    REAL(KIND=double), INTENT(IN) :: x
    REAL(KIND=double) :: fp1, fp2, fp3, fp4

    xtmp(:) = pcom(:) + x*xicom(:)
    CALL GETCOM(xtmp,com)

    CALL ENERGY_LJ_EL(xtmp,com,fp1,fp2,fp3,mutmp,dipole,dipole_ind,polv,polv_perm)
    CALL ENERGY_BOND(xtmp,fp4)
    F1DIM = fp1 + fp2 + fp3 + fp4

 

    END FUNCTION

	
    REAL(KIND=double) FUNCTION DF1DIM(x)
    IMPLICIT NONE
    REAL(KIND=double), INTENT(IN) :: x
    REAL(KIND=double) :: fp1, fp2, fp3, fp4

    xtmp(:) = pcom(:) + x*xicom(:)
    CALL GETCOM(xtmp,com)
    df = 0d0
    CALL FORCE_BOND(xtmp,df)
    CALL FORCE_LJ_EL(xtmp,com,df,fp1,fp2,fp3,0,mutmp)
    df = -df
    DF1DIM = DOT_PRODUCT(df,xicom)

    END FUNCTION


	
	
	
	
SUBROUTINE INIT_COM
IMPLICIT NONE


n = 0

IF(nmol <= 32)THEN
  DO i = 0,1
    DO j = 0,3
      DO k = 0,3
        n = n + 1
		IF(n>nmol) EXIT
        ix = ( n - 1 ) * 3 + 1
        iy = ix + 1
        iz = iy + 1
        com(ix) = pbox0 * dble(i) / 2d0 + 0.01d0
        if(mod(i+1,2)==0)then
        com(iy) = pbox0 * dble(j) / 4d0 + 0.01d0
        else
        com(iy) = pbox0 * dble(j) / 4d0 + pbox0/8d0
        endif
        if(mod(j+1,2)==0)then
        com(ix) = com(ix) + pbox0/8d0
        endif
        if(mod(i+1,2)==0)then
        com(iz) = pbox0 * dble(k) / 4d0 + 0.01d0
        else
        com(iz) = pbox0 * dble(k) / 4d0 + pbox0/8d0
        endif
        if(mod(k+1,2)==0)then
        com(ix) = com(ix) + pbox0/8d0
        endif
      ENDDO
    ENDDO
  ENDDO
  
  
ELSE IF(nmol <= 64)THEN
  DO i = 0,3
    DO j = 0,3
      DO k = 0,3
        n = n + 1
		IF(n>nmol) EXIT
        ix = ( n - 1 ) * 3 + 1
        iy = ix + 1
        iz = iy + 1
        com(ix) = pbox0 * dble(i) / 4d0 + 0.01d0
        if(mod(i+1,2)==0)then
        com(iy) = pbox0 * dble(j) / 4d0 + 0.01d0
        else
        com(iy) = pbox0 * dble(j) / 4d0 + pbox0/8d0
        endif
        if(mod(j+1,2)==0)then
        com(ix) = com(ix) + pbox0/8d0
        endif
        if(mod(i+1,2)==0)then
        com(iz) = pbox0 * dble(k) / 4d0 + 0.01d0
        else
        com(iz) = pbox0 * dble(k) / 4d0 + pbox0/8d0
        endif
        if(mod(k+1,2)==0)then
        com(ix) = com(ix) + pbox0/8d0
        endif
      ENDDO
    ENDDO
  ENDDO
  

ELSE IF(nmol <= 108)THEN
  DO i = 0,2
    DO j = 0,5
      DO k = 0,5
        n = n + 1
		IF(n>nmol) EXIT
        ix = ( n - 1 ) * 3 + 1
        iy = ix + 1
        iz = iy + 1
        IF(mod(j+1,2)==0)THEN
        com(ix) = pbox0 * dble(i) / 3d0 + 0.01d0
        else
        com(ix) = pbox0 * dble(i) / 3d0 + pbox0/6d0
        endif
        if(mod(i+1,2)==0)then
        com(iy) = pbox0 * dble(j) / 6d0 + 0.01d0
        else
        com(iy) = pbox0 * dble(j) / 6d0 + pbox0/12d0
        endif
        if(mod(i+1,2)==0)then
        com(iz) = pbox0 * dble(k) / 6d0 + 0.01d0
        else
        com(iz) = pbox0 * dble(k) / 6d0 + pbox0/12d0
        endif
      ENDDO
    ENDDO
  ENDDO




ELSE IF(nmol <= 216)THEN
  DO i = 0,5
    DO j = 0,5
      DO k = 0,5
        n = n + 1
		IF(n>nmol) EXIT
        ix = ( n - 1 ) * 3 + 1
        iy = ix + 1
        iz = iy + 1
        IF(mod(j+1,2)==0)THEN
        com(ix) = pbox0 * dble(i) / 6d0 + 0.01d0
        else
        com(ix) = pbox0 * dble(i) / 6d0 + pbox0/12d0
        endif
        if(mod(i+1,2)==0)then
        com(iy) = pbox0 * dble(j) / 6d0 + 0.01d0
        else
        com(iy) = pbox0 * dble(j) / 6d0 + pbox0/12d0
        endif
        if(mod(i+1,2)==0)then
        com(iz) = pbox0 * dble(k) / 6d0 + 0.01d0
        else
        com(iz) = pbox0 * dble(k) / 6d0 + pbox0/12d0
        endif
      ENDDO
    ENDDO
  ENDDO
  
ELSE IF(nmol <= 432)THEN
  DO i = 0,8
    DO j = 0,7
      DO k = 0,5
        n = n + 1
		IF(n>nmol) EXIT
        ix = ( n - 1 ) * 3 + 1
        iy = ix + 1
        iz = iy + 1
        IF(mod(j+1,2)==0)THEN
        com(ix) = pbox0 * dble(i) / 9d0 + 0.01d0
        else
        com(ix) = pbox0 * dble(i) / 9d0 + pbox0/18d0
        endif
        if(mod(i+1,2)==0)then
        com(iy) = pbox0 * dble(j) / 8d0 + 0.01d0
        else
        com(iy) = pbox0 * dble(j) / 8d0 + pbox0/16d0
        endif
        if(mod(i+1,2)==0)then
        com(iz) = pbox0 * dble(k) / 6d0 + 0.01d0
        else
        com(iz) = pbox0 * dble(k) / 6d0 + pbox0/12d0
        endif
      ENDDO
    ENDDO
  ENDDO

ENDIF


END SUBROUTINE INIT_COM



SUBROUTINE Solve_Newton(dt)
IMPLICIT NONE
REAL(KIND=double),INTENT(IN) :: dt
INTEGER :: i,j,k


IF(jobtype=="EMD_6th")THEN 

    DO i = 1,7
      CALL EVFAST(dt*0.5d0*om(i))
    ENDDO
	
    for_s(:) = 0d0
    CALL GETCOM(cor,com)
    CALL FORCE_LJ_EL(cor,com,for_s,elj,eel,epol,1,mutmp)
	
    ptt(:) = ptt(:) + for_s(:)*dt
	
    DO i = 1,7
      CALL EVFAST(dt*0.5d0*om(i))
    ENDDO
	
ELSE !V-velet
    
	! (t+h <= t + t+h/2)
	ptt(:) = ptt(:) + for_s(:) * dt * 0.5d0
	
	DO k = 1, nfst
	
	 ptt(:) = ptt(:) + for_f(:) * dtf * 0.5d0 
	 
    DO i = 1,nmol
      ix = ( i - 1 )*9 + 1
      iy = ix + 1
      iz = ix + 2
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf * massi_h
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf * massi_h
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf * massi_o
    ENDDO !i
    
     for_f(:) = 0d0
     CALL FORCE_BOND(cor,for_f)		
     ptt(:) = ptt(:) + for_f(:) * dtf * 0.5d0
	 
	END DO !k
	
	 for_s(:) = 0d0
	 CALL GETCOM(cor,com)
     CALL FORCE_LJ_EL(cor,com,for_s,elj,eel,epol,1,mutmp)
	 IF(NoneqTHzchk)THEN
	  CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
	  CALL FORCE_ELEC_THz(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
	 ELSE IF(NoneqRAMANchk)THEN 
	  CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
	  CALL FORCE_ELEC_RAMAN(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
	 END IF	
	
	 ! (t+h <= t + t+h/2)
	 ptt(:) = ptt(:) + for_s(:) * dt * 0.5d0
	
END IF



END SUBROUTINE Solve_Newton


SUBROUTINE SINRCOR(dt)
IMPLICIT NONE
REAL(KIND=double),INTENT(IN) :: dt
INTEGER :: i,j,k

DO i = 1, natom*3
   CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
END DO
DO i = 1, natom*3
   CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i)+for_s(i)*nfst,dtf*0.5d0)
END DO
cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
DO i = 1, natom*3
   CALL SINR_O(psin2(:,i),dtf)
END DO
cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
IF(nfst>1) THEN
   for_f = 0d0
   CALL FORCE_BOND(cor,for_f)
   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO
   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO
END IF

DO k = 1, nfst - 2

   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO
   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
   DO i = 1, natom*3
      CALL SINR_O(psin2(:,i),dtf)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
   for_f = 0d0
   CALL FORCE_BOND(cor,for_f)		
   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO
   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO

END DO !k

IF(nfst>1) THEN
   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO
   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
   DO i = 1, natom*3
      CALL SINR_O(psin2(:,i),dtf)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
END IF
for_f = 0d0
CALL FORCE_BOND(cor,for_f)		
for_s(:) = 0d0
CALL GETCOM(cor,com)
CALL FORCE_LJ_EL(cor,com,for_s,elj,eel,epol,1,mutmp)
IF(NoneqTHzchk)THEN
   CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
   CALL FORCE_ELEC_THz(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
ELSE IF(NoneqRAMANchk)THEN 
   CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
   CALL FORCE_ELEC_RAMAN(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
END IF	
DO i = 1, natom*3
   CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i)+for_s(i)*nfst,dtf*0.5d0)
END DO
DO i = 1, natom*3
   CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
END DO

END SUBROUTINE SINRCOR


SUBROUTINE SINRMCOR(dt)
IMPLICIT NONE
REAL(KIND=double),INTENT(IN) :: dt
INTEGER :: i,j,k

DO i = 1, natom*3
   CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i)+for_s(i)*nfst,dtf*0.5d0)
END DO
cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
DO i = 1, natom*3
   CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
END DO
DO i = 1, natom*3
   CALL SINR_O(psin2(:,i),dtf)
END DO
DO i = 1, natom*3
   CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
END DO
cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
IF(nfst>1) THEN
   for_f = 0d0
   CALL FORCE_BOND(cor,for_f)
   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO
END IF

DO k = 1, nfst - 2

   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO
   DO i = 1, natom*3
      CALL SINR_O(psin2(:,i),dtf)
   END DO
   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
   for_f = 0d0
   CALL FORCE_BOND(cor,for_f)		
   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO

END DO !k

IF(nfst>1) THEN
   DO i = 1, natom*3
      CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i),dtf*0.5d0)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO
   DO i = 1, natom*3
      CALL SINR_O(psin2(:,i),dtf)
   END DO
   DO i = 1, natom*3
      CALL SINR_N(ptt(i),mass(i),psin1(:,i),psin2(:,i),dtf*0.5d0)
   END DO
   cor(:) = cor(:) + ptt(:) / mass(:) * dtf / 2
END IF
for_f = 0d0
CALL FORCE_BOND(cor,for_f)		
for_s(:) = 0d0
CALL GETCOM(cor,com)
CALL FORCE_LJ_EL(cor,com,for_s,elj,eel,epol,1,mutmp)
IF(NoneqTHzchk)THEN
   CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
   CALL FORCE_ELEC_THz(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
ELSE IF(NoneqRAMANchk)THEN 
   CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
   CALL FORCE_ELEC_RAMAN(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
END IF	
DO i = 1, natom*3
   CALL SINR_V(ptt(i),mass(i),psin1(:,i),for_f(i)+for_s(i)*nfst,dtf*0.5d0)
END DO

END SUBROUTINE SINRMCOR


SUBROUTINE SINR_N(p, mass, v1k, v2k, dtime)
IMPLICIT NONE
REAL(KIND=double), INTENT(INOUT) :: p, v1k(:), v2k(:)
REAL(KIND=double), INTENT(IN) :: mass, dtime
REAL(KIND=double) :: H
INTEGER :: alpha, j
DO alpha = 1, 7
   DO j = 1, nchain
      v2k(j) = v2k(j) + (Qsin1*v1k(j)**2-kbt)/Qsin2*om(alpha)*dtime/2
   END DO
   H = sqrt(nchain*kbt/(p**2/mass+nchain*Qsin1*sum(v1k**2*exp(-2*v2k*om(alpha)*dtime))/(nchain+1)))
   p = p * H
   DO j = 1, nchain
      v1k(j) = v1k(j) * H * exp(-v2k(j)*om(alpha)*dtime)
   END DO
   DO j = 1, nchain
      v2k(j) = v2k(j) + (Qsin1*v1k(j)**2-kbt)/Qsin2*om(alpha)*dtime/2
   END DO
END DO
END SUBROUTINE SINR_N


SUBROUTINE SINR_V(p, mass, v1k, f, dtime)
IMPLICIT NONE
REAL(KIND=double), INTENT(INOUT) :: p, v1k(:)
REAL(KIND=double), INTENT(IN) :: mass, f, dtime
REAL(KIND=double) :: v0, s, sd1
integer :: j
v0 = p / mass
IF(abs(f)*dtime/sqrt(mass*kbt)<1d-5) THEN
   p = p + ((f/mass-f*v0**2/kbt/nchain)*dtime-f**2*(v0*nchain*kbt-mass*v0**3)*dtime**2/((nchain*kbt)**2*mass)&
        -f**3*((nchain*kbt)**2-4*nchain*mass*v0**2*kbt+3*mass**2*v0**4)*dtime**3/(3*(nchain*kbt)**3*mass**2)&
        +f**4*(2*(nchain*kbt)**2*v0-5*nchain*mass*v0**3*kbt+3*mass**2*v0**5)*dtime**4/(3*(nchain*kbt)**4*mass**2))*mass
   DO j = 1, nchain
      v1k(j) = v1k(j) - f*v0*v1k(j)*dtime/kbt/nchain-f**2*(nchain*kbt-2*mass*v0**2)*v1k(j)*dtime**2/(2*(nchain*kbt)**2*mass)&
           +f**3*(5*nchain*v0*kbt-6*mass*v0**3)*v1k(j)*dtime**3/(6*(nchain*kbt)**3*mass)&
           +f**4*(5*(nchain*kbt)**2-28*nchain*mass*v0**2*kbt+24*mass**2*v0**4)*v1k(j)*dtime**4/(24*(nchain*kbt)**4*mass**2)
   END DO
ELSE
   s = 1/sqrt(f**2/mass/nchain/kbt)*sinh(sqrt(f**2/mass/nchain/kbt)*dtime)+f*v0/nchain/kbt/(f**2/mass/nchain/kbt)&
        *(cosh(sqrt(f**2/mass/nchain/kbt)*dtime)-1)
   sd1 = cosh(sqrt(f**2/mass/nchain/kbt)*dtime)+f*v0/nchain/kbt/sqrt(f**2/mass/nchain/kbt)*&
        sinh(sqrt(f**2/mass/nchain/kbt)*dtime)
   p = (p + f * s) / sd1
   v1k = v1k / sd1
END IF
END SUBROUTINE SINR_V


SUBROUTINE SINR_O(v2k,dtime)
REAL(KIND=double),INTENT(INOUT) :: v2k(:)
REAL(KIND=double),INTENT(IN) :: dtime
REAL(KIND=double) :: c1, c2
INTEGER :: j

c1=exp(-gsin*dtime)
c2=sqrt(1-c1**2)

DO j = 1, nchain
   CALL gasdev(eta)
   v2k(j) = c1*v2k(j)+c2*sqrt(kbt/Qsin2)*eta
END DO

END SUBROUTINE SINR_O


SUBROUTINE BATABCOR(dt)
IMPLICIT NONE
REAL(KIND=double),INTENT(IN) :: dt
INTEGER :: i,j,k


! (t+h <= t + t+h/2)
ptt(:) = ptt(:) + for_s(:) * dt * 0.5d0
!WRITE(*,*) 'ptt', 0
!WRITE(*,*) ptt(1:3)

DO k = 1, nfst

   ptt(:) = ptt(:) + for_f(:) * dtf * 0.5d0 
   !WRITE(*,*) 'ptt', k, 1
   !WRITE(*,*) ptt(1:3)

   DO i = 1,nmol
      ix = ( i - 1 )*9 + 1
      iy = ix + 1
      iz = ix + 2
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf / 2 * massi_h
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf / 2 * massi_h
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf / 2 * massi_o
   ENDDO !i
   !WRITE(*,*) 'cor', k, 1
   !WRITE(*,*) cor(1:3)


   IF(cjob == 'BANAB')THEN

      DO i = 1,7
         CALL NVTCOR(ptt,pnvt,snvt,dtf*om(i))
      ENDDO

   ELSE IF(cjob == 'MASSB')THEN

      DO j = 1, 3*natom
         DO i = 1,7
            IF(mod(j,9)>=1.and.mod(j,9)<=6)THEN
               CALL MASSIVE_NVTCOR(ptt(j),pnvtm(:,j),snvtm(:,j),dtf*om(i),mass_h)
            ELSE
               CALL MASSIVE_NVTCOR(ptt(j),pnvtm(:,j),snvtm(:,j),dtf*om(i),mass_o)
            END IF
         END DO
      END DO

   ELSE IF(cjob == 'BAOAB')THEN

      CALL LangCOR(ptt,dtf)

   ELSE IF(cjob == 'Ander')THEN

      CALL Andersen_Thermostat(ptt,dtf)

   END IF
   !WRITE(*,*) 'ptt', k, 2
   !WRITE(*,*) ptt(1:3)

   DO i = 1,nmol
      ix = ( i - 1 )*9 + 1
      iy = ix + 1
      iz = ix + 2
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf / 2 * massi_h
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf / 2 * massi_h
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      cor(ix:iz) = cor(ix:iz) + ptt(ix:iz) * dtf / 2 * massi_o
   ENDDO !i
   !WRITE(*,*) 'cor', k, 2
   !WRITE(*,*) cor(1:3)

   for_f(:) = 0d0
   CALL FORCE_BOND(cor,for_f)		
   ptt(:) = ptt(:) + for_f(:) * dtf * 0.5d0
   !WRITE(*,*) 'ptt', k, 3
   !WRITE(*,*) ptt(1:3)

END DO !k

for_s(:) = 0d0
CALL GETCOM(cor,com)
CALL FORCE_LJ_EL(cor,com,for_s,elj,eel,epol,1,mutmp)
IF(NoneqTHzchk)THEN
   CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
   CALL FORCE_ELEC_THz(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
ELSE IF(NoneqRAMANchk)THEN 
   CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
   CALL FORCE_ELEC_RAMAN(cor,com,for_s,Elec,polv,polv_perm,NEMDtype)
END IF	

! (t+h <= t + t+h/2)
ptt(:) = ptt(:) + for_s(:) * dt * 0.5d0
!WRITE(*,*) 'ptt', 2
!WRITE(*,*) ptt(1:3)


END SUBROUTINE BATABCOR

  
SUBROUTINE LangCOR(ptt,dtime) ! Langevin
REAL(KIND=double),INTENT(INOUT) :: ptt(:)
REAL(KIND=double),INTENT(IN) :: dtime
REAL(KIND=double) :: c1, c2
INTEGER :: i

c1=exp(-lgam*dtime)
c2=sqrt(1-c1**2)

DO i = 1, 3*natom
   CALL gasdev(eta)
   IF(mod(i,9)>=1.and.mod(i,9)<=6)THEN
      ptt(i)=c1*ptt(i)+c2*sqrt(kbt*mass_h)*eta
   ELSE
      ptt(i)=c1*ptt(i)+c2*sqrt(kbt*mass_o)*eta
   END IF
END DO

END SUBROUTINE LangCOR


SUBROUTINE Andersen_Thermostat(ptt,dtime) ! Andersen
REAL(KIND=double),INTENT(INOUT) :: ptt(:)
REAL(KIND=double),INTENT(IN) :: dtime
INTEGER :: i, j

DO i = 1, natom
   eta = genrand_real()

   IF(eta < nu*dtime)THEN
      DO j = 3*i-2, 3*i
         CALL gasdev(eta)
         IF(mod(i,3)>=1.and.mod(i,3)<=2)THEN
            ptt(j)=sqrt(kbt*mass_h)*eta
         ELSE
            ptt(j)=sqrt(kbt*mass_o)*eta
         END IF
      END DO
   END IF

END DO

END SUBROUTINE Andersen_Thermostat


SUBROUTINE GET_PROPERTY
IMPLICIT NONE
INTEGER :: i

      CALL GETCOM(cor,com)

      CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
      CALL ENERGY_BOND(cor,ebond)
      CALL GETKENG(ptt,keng,ktemp)


   IF(NoneqTHzchk)THEN
	IF(NEMDtype=="Total")THEN
	 eel  = eel  - (dipt(3))*Elec*23.06013831069d0
	 epol = epol - (dipit(3))*Elec*23.06013831069d0
	ELSE IF(NEMDtype=="Perm")THEN
	 eel  = eel  - (dipt(3))*Elec*23.06013831069d0
	ELSE IF(NEMDtype=="Ind")THEN
	 epol = epol - (dipit(3))*Elec*23.06013831069d0
	END IF
   ELSE IF(NoneqRAMANchk)THEN
	IF(NEMDtype=="Total")THEN
	 epol = epol - 0.5d0*polt(3)*(Elec*23.06013831069d0)**2 &
	 & /(1d0/1.112650056d-10*(1.60217733d-19)**2/(1.660540d-21*4.184d0)*1d10)
	ELSE IF(NEMDtype=="Perm")THEN
	 epol = epol - 0.5d0*polpt(3)*(Elec*23.06013831069d0)**2 &
	 & /(1d0/1.112650056d-10*(1.60217733d-19)**2/(1.660540d-21*4.184d0)*1d10)
    ELSE IF(NEMDtype=="Ind")THEN
	 epol = epol - 0.5d0*(polt(3)-polpt(3))*(Elec*23.06013831069d0)**2 &
	 & /(1d0/1.112650056d-10*(1.60217733d-19)**2/(1.660540d-21*4.184d0)*1d10)	 
	ELSE IF(NEMDtype=="Totalxy")THEN
	 epol = epol - 0.5d0*polt(4)*(Elec*23.06013831069d0)**2 &
	 & /(1d0/1.112650056d-10*(1.60217733d-19)**2/(1.660540d-21*4.184d0)*1d10)
    ELSE IF(NEMDtype=="Permxy")THEN
	 epol = epol - 0.5d0*polpt(4)*(Elec*23.06013831069d0)**2 &
	 & /(1d0/1.112650056d-10*(1.60217733d-19)**2/(1.660540d-21*4.184d0)*1d10)
    ELSE IF(NEMDtype=="Indxy")THEN
	 epol = epol - 0.5d0*(polt(4)-polpt(4))*(Elec*23.06013831069d0)**2 &
	 & /(1d0/1.112650056d-10*(1.60217733d-19)**2/(1.660540d-21*4.184d0)*1d10)
	END IF
   END IF
	  
	  
      volume = SETPB(0d0)**3
      prop%pbox = SETPB(0d0)
      prop%density = tmass/(volume*1d-24)

      IF(cjob /= 'NPT')THEN
        proptmp = GETPROP()
        prop%pint(1:10) = proptmp%pint(1:10)
      ENDIF
   
      prop%enel = elj
      prop%eel = eel
      prop%epol = epol
      prop%ebond = ebond
      prop%keng = keng
      prop%temp = ktemp
      ntemp = UPDATEAVG(prop)

      IF(cjob == 'NVT'.or.cjob == 'BANAB')THEN
        envt = gkt * snvt(1) 
        envt = envt + 0.5d0*pnvt(1)**2*qt
        DO i = 2,nchain
          envt = envt + gk1 *snvt(i) 
          envt = envt + 0.5d0*pnvt(i)**2*qt1
        ENDDO

        CALL OUTPUT(n,nttime,dt,prop,envt)
     ELSE IF(cjob == 'MASSB' .or. cjob == 'MASSN')THEN
        envt = gkt * sum(snvtm(1,:)) 
        envt = envt + 0.5d0*sum(pnvtm(1,:)**2)*qt
        DO i = 2,nchain
           envt = envt + gk1 *sum(snvtm(i,:)) 
           envt = envt + 0.5d0*sum(pnvtm(i,:)**2)*qt1
        ENDDO

        CALL OUTPUT(n,nttime,dt,prop,envt)
      ELSE IF(cjob == 'NPT')THEN
        envt = gkt * snvt(1) 
        envt = envt + 0.5d0*pnvt(1)**2*qt
        DO i = 2,nchain
          envt = envt + gk1 *snvt(i) 
          envt = envt + 0.5d0*pnvt(i)**2*qt1
        ENDDO
        enpt =  0.5d0*pnpt**2*qtp + pext*volume

        CALL OUTPUT(n,nttime,dt,prop,envt,enpt)

      ELSE 

        CALL OUTPUT(n,nttime,dt,prop)
      ENDIF
  
  
!===========================================================================
!  output charge and dipole 
!===========================================================================

  IF(mod(n,nout)== 0)THEN
 		 					
    CALL CONFOUT(roh,rhh,th)
	rohav = rohav + roh
	rhhav = rhhav + rhh
	thav  = thav  + th
	nconf = nconf + 1
	
	 WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <rOH>     =",rohav/dble(nconf),"   (A)",&
                                &"-  rOH    =",roh, "   (A)"		
	 WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <rHH>     =",rhhav/dble(nconf),"   (A)",&
                                &"-  rHH    =",rhh, "   (A)"		
	 WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <theta>   =",thav/dble(nconf), " (deg)",&
                                &"-  theta  =",th,  " (deg)"		
								
    IF(ndip /=0)THEN
	 WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <qHperm>  =",qHpav/dble(ndip),"   (e)",&
                                &"-  qHperm =",qHperm, "   (e)"				
	 WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <qHind>   =",qHiav/dble(ndip),"   (e)",&
                                &"-  qHind  =",qHind,  "   (e)"
	 WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <qOperm>  =",qOpav/dble(ndip),"   (e)",&
                                &"-  qOperm =",qOperm, "   (e)"		
	 WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <qOind>   =",qOiav/dble(ndip),"   (e)",&
                                &"-  qOind  =",qOind,  "   (e)"
								
     WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <DIPOLE Perm> =",dppav/dble(ndip)/0.20819d0,"  (Debye)",&
                                &"- DIPOLE Perm =",dppinst/0.20819d0,"  (Debye)"
     WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <DIPOLE>      =",dipav/dble(ndip)/0.20819d0,"  (Debye)",&
                                &"- DIPOLE      =",dipinst/0.20819d0,"  (Debye)"
     WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <POLISO Perm> =",popav/dble(ndip),  "    (A^3)",&
                                &"- POLISO Perm =",popinst,  "    (A^3)"
     WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <POLISO>      =",polav/dble(ndip),  "    (A^3)",&
                                &"- POLISO      =",polinst,  "    (A^3)"
     WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <POLANI Perm> =",anpopav/dble(ndip),"    (   )",&
                                &"- POLANI Perm =",anpopinst,"    (   )"
     WRITE(9,'(a,f12.5,a,a,f12.5,a)') "- <POLANI>      =",anpolav/dble(ndip),"    (   )",&
                                &"- POLANI      =",anpolinst,"    (   )"
    ENDIF
	
	 WRITE(9,*)
	 
  ENDIF
  
  
END SUBROUTINE GET_PROPERTY


!Polar	
SUBROUTINE GET_TRAJECTRY
IMPLICIT NONE
INTEGER :: i
INTEGER :: ix,iy,iz

!===================================================
!polarization
!===================================================
	
    IF(n==1 .and. npol /=0)THEN
	  dppav = 0d0
	  dipav = 0d0
	  popav = 0d0
	  polav = 0d0
	  anpopav = 0d0
	  anpolav = 0d0
      ndip  = 0
	  
	  qHpav = 0d0
	  qHiav = 0d0
	  qOpav = 0d0
	  qOiav = 0d0	  
	END IF
	
    IF(npol /= 0 .and. mod(n,npol) == 0)THEN
	
      CALL GETCOM(cor,com)
      CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)
	  
      dipt(:) = 0d0
      dipit(:) = 0d0
	  polpt(:) = 0d0
      polt(:) = 0d0
	  dppinst = 0d0
      dipinst = 0d0
	  polinst = 0d0
      DO i = 1, nmol
        ix = (i-1)*3 + 1
        iz = ix + 2
        dipt(1:3) = dipt(1:3) + dipole(ix:iz)
        dipit(1:3) = dipit(1:3) + dipole_ind(ix:iz)
        dppav = dppav + SQRT(DOT_PRODUCT(dipole(ix:iz),dipole(ix:iz)))/dble(nmol)
        dppinst = dppinst + SQRT(DOT_PRODUCT(dipole(ix:iz),dipole(ix:iz)))/dble(nmol)
		
        dipav = dipav + SQRT(DOT_PRODUCT(dipole(ix:iz)+dipole_ind(ix:iz) ,&
        &               dipole(ix:iz) + dipole_ind(ix:iz)))/dble(nmol)
        dipinst = dipinst + SQRT(DOT_PRODUCT(dipole(ix:iz)+dipole_ind(ix:iz) ,&
        &               dipole(ix:iz) + dipole_ind(ix:iz)))/dble(nmol) 
      ENDDO
	  
      DO i = 1, nmol
        ix = 9*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
        polt(1) = polt(1) + polv(ix,1) + polv(ix + 3,1) + polv(ix + 6,1) 
        polt(2) = polt(2) + polv(iy,2) + polv(iy + 3,2) + polv(iy + 6,2) 
        polt(3) = polt(3) + polv(iz,3) + polv(iz + 3,3) + polv(iz + 6,3) 
        polt(4) = polt(4) + polv(ix,2) + polv(ix + 3,2) + polv(ix + 6,2) 
        polt(5) = polt(5) + polv(iy,3) + polv(iy + 3,3) + polv(iy + 6,3) 
        polt(6) = polt(6) + polv(iz,1) + polv(iz + 3,1) + polv(iz + 6,1) 
        polt(7) = polt(7) + polv(iy,1) + polv(iy + 3,1) + polv(iy + 6,1) 
        polt(8) = polt(8) + polv(iz,2) + polv(iz + 3,2) + polv(iz + 6,2) 
        polt(9) = polt(9) + polv(ix,3) + polv(ix + 3,3) + polv(ix + 6,3) 

        polpt(1) = polpt(1) + polv_perm(ix,1) + polv_perm(ix + 3,1) + polv_perm(ix + 6,1) 
        polpt(2) = polpt(2) + polv_perm(iy,2) + polv_perm(iy + 3,2) + polv_perm(iy + 6,2) 
        polpt(3) = polpt(3) + polv_perm(iz,3) + polv_perm(iz + 3,3) + polv_perm(iz + 6,3) 
        polpt(4) = polpt(4) + polv_perm(ix,2) + polv_perm(ix + 3,2) + polv_perm(ix + 6,2) 
        polpt(5) = polpt(5) + polv_perm(iy,3) + polv_perm(iy + 3,3) + polv_perm(iy + 6,3) 
        polpt(6) = polpt(6) + polv_perm(iz,1) + polv_perm(iz + 3,1) + polv_perm(iz + 6,1) 
        polpt(7) = polpt(7) + polv_perm(iy,1) + polv_perm(iy + 3,1) + polv_perm(iy + 6,1) 
        polpt(8) = polpt(8) + polv_perm(iz,2) + polv_perm(iz + 3,2) + polv_perm(iz + 6,2) 
        polpt(9) = polpt(9) + polv_perm(ix,3) + polv_perm(ix + 3,3) + polv_perm(ix + 6,3) 		
      ENDDO

      DO i = 1, nmol
        j = 2*(i-1) + 1
        ix = 9*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
        jx = ix + 3 
        jy = ix + 4 
        jz = ix + 5 
        kx = ix + 6 
        ky = ix + 7 
        kz = ix + 8 
        polt(1) = polt(1) + polv(9*nmol + j,1)*(cor(ix)-cor(kx)) &
        &                 + polv(9*nmol + j + 1,1)*(cor(jx)-cor(kx))
        polt(2) = polt(2) + polv(9*nmol + j,2)*(cor(iy)-cor(ky)) &
        &                 + polv(9*nmol + j + 1,2)*(cor(jy)-cor(ky))
        polt(3) = polt(3) + polv(9*nmol + j,3)*(cor(iz)-cor(kz)) &
        &                 + polv(9*nmol + j + 1,3)*(cor(jz)-cor(kz))
        polt(7) = polt(7) + polv(9*nmol + j,1)*(cor(iy)-cor(ky)) &
        &                 + polv(9*nmol + j + 1,1)*(cor(jy)-cor(ky))
        polt(8) = polt(8) + polv(9*nmol + j,2)*(cor(iz)-cor(kz)) &
        &                 + polv(9*nmol + j + 1,2)*(cor(jz)-cor(kz))
        polt(9) = polt(9) + polv(9*nmol + j,3)*(cor(ix)-cor(kx)) &
        &                 + polv(9*nmol + j + 1,3)*(cor(jx)-cor(kx))
        polt(4) = polt(4) + polv(9*nmol + j,2)*(cor(ix)-cor(kx)) &
        &                 + polv(9*nmol + j + 1,2)*(cor(jx)-cor(kx))
        polt(5) = polt(5) + polv(9*nmol + j,3)*(cor(iy)-cor(ky)) &
        &                 + polv(9*nmol + j + 1,3)*(cor(jy)-cor(ky))
        polt(6) = polt(6) + polv(9*nmol + j,1)*(cor(iz)-cor(kz)) &
        &                 + polv(9*nmol + j + 1,1)*(cor(jz)-cor(kz))
		
        polpt(1) = polpt(1) + polv_perm(9*nmol + j,1)*(cor(ix)-cor(kx)) &
        &                 + polv_perm(9*nmol + j + 1,1)*(cor(jx)-cor(kx))
        polpt(2) = polpt(2) + polv_perm(9*nmol + j,2)*(cor(iy)-cor(ky)) &
        &                 + polv_perm(9*nmol + j + 1,2)*(cor(jy)-cor(ky))
        polpt(3) = polpt(3) + polv_perm(9*nmol + j,3)*(cor(iz)-cor(kz)) &
        &                 + polv_perm(9*nmol + j + 1,3)*(cor(jz)-cor(kz))
        polpt(7) = polpt(7) + polv_perm(9*nmol + j,1)*(cor(iy)-cor(ky)) &
        &                 + polv_perm(9*nmol + j + 1,1)*(cor(jy)-cor(ky))
        polpt(8) = polpt(8) + polv_perm(9*nmol + j,2)*(cor(iz)-cor(kz)) &
        &                 + polv_perm(9*nmol + j + 1,2)*(cor(jz)-cor(kz))
        polpt(9) = polpt(9) + polv_perm(9*nmol + j,3)*(cor(ix)-cor(kx)) &
        &                 + polv_perm(9*nmol + j + 1,3)*(cor(jx)-cor(kx))
        polpt(4) = polpt(4) + polv_perm(9*nmol + j,2)*(cor(ix)-cor(kx)) &
        &                 + polv_perm(9*nmol + j + 1,2)*(cor(jx)-cor(kx))
        polpt(5) = polpt(5) + polv_perm(9*nmol + j,3)*(cor(iy)-cor(ky)) &
        &                 + polv_perm(9*nmol + j + 1,3)*(cor(jy)-cor(ky))
        polpt(6) = polpt(6) + polv_perm(9*nmol + j,1)*(cor(iz)-cor(kz)) &
        &                 + polv_perm(9*nmol + j + 1,1)*(cor(jz)-cor(kz))		
      ENDDO       

      polt(4) = 0.5d0*(polt(4) + polt(7))
      polt(5) = 0.5d0*(polt(5) + polt(8))
      polt(6) = 0.5d0*(polt(6) + polt(9))
	  
      polpt(4) = 0.5d0*(polpt(4) + polpt(7))
      polpt(5) = 0.5d0*(polpt(5) + polpt(8))
      polpt(6) = 0.5d0*(polpt(6) + polpt(9))
	  
	  polinst = ( polt(1) + polt(2) + polt(3) )/3d0/dble(nmol)
	  polav = polav + polinst

	  popinst = ( polpt(1) + polpt(2) + polpt(3) )/3d0/dble(nmol)
	  popav = popav + popinst
	  
	  anpolinst = (polt(1)/dble(nmol)-polinst)**2 + &
	  & (polt(2)/dble(nmol)-polinst)**2 + (polt(3)/dble(nmol)-polinst)**2
	  anpolinst = SQRT( anpolinst / 6d0 ) / polinst
	  anpolav = anpolav + anpolinst
	   
	  anpopinst = (polpt(1)/dble(nmol)-popinst)**2 + &
	  & (polpt(2)/dble(nmol)-popinst)**2 + (polpt(3)/dble(nmol)-popinst)**2
	  anpopinst = SQRT( anpopinst / 6d0 ) / popinst
	  anpopav = anpopav + anpopinst
	   
	  
      IF(fdpdt)THEN
        edt = dt * 1d-3
        cor_tmp = cor
        ptt_tmp = ptt
		for_stmp = for_s
        for_ftmp = for_f
        
	    CALL Solve_Newton(edt)

        CALL GETCOM(cor,com)
        CALL ENERGY_LJ_EL(cor,com,elj,eel,epol,mutmp,dipole,dipole_ind,polv,polv_perm)

        cor = cor_tmp
        ptt = ptt_tmp
        for_s = for_stmp
        for_f = for_ftmp 

        diptp(:) = 0d0
        dipitp(:) = 0d0
		polptp(:) = 0d0
        poltp(:) = 0d0
        DO i = 1, nmol
          ix = (i-1)*3 + 1
          iz = ix + 2
          diptp(1:3) = diptp(1:3) + dipole(ix:iz)
          dipitp(1:3) = dipitp(1:3) + dipole_ind(ix:iz)
        ENDDO
		
        DO i = 1, nmol
         ix = 9*(i-1) + 1
         iy = ix + 1
         iz = ix + 2
         poltp(1) = poltp(1) + polv(ix,1) + polv(ix + 3,1) + polv(ix + 6,1) 
         poltp(2) = poltp(2) + polv(iy,2) + polv(iy + 3,2) + polv(iy + 6,2) 
         poltp(3) = poltp(3) + polv(iz,3) + polv(iz + 3,3) + polv(iz + 6,3) 
         poltp(4) = poltp(4) + polv(ix,2) + polv(ix + 3,2) + polv(ix + 6,2) 
         poltp(5) = poltp(5) + polv(iy,3) + polv(iy + 3,3) + polv(iy + 6,3) 
         poltp(6) = poltp(6) + polv(iz,1) + polv(iz + 3,1) + polv(iz + 6,1) 
         poltp(7) = poltp(7) + polv(iy,1) + polv(iy + 3,1) + polv(iy + 6,1) 
         poltp(8) = poltp(8) + polv(iz,2) + polv(iz + 3,2) + polv(iz + 6,2) 
         poltp(9) = poltp(9) + polv(ix,3) + polv(ix + 3,3) + polv(ix + 6,3)
		 
         polptp(1) = polptp(1) + polv_perm(ix,1) + polv_perm(ix + 3,1) + polv_perm(ix + 6,1) 
         polptp(2) = polptp(2) + polv_perm(iy,2) + polv_perm(iy + 3,2) + polv_perm(iy + 6,2) 
         polptp(3) = polptp(3) + polv_perm(iz,3) + polv_perm(iz + 3,3) + polv_perm(iz + 6,3) 
         polptp(4) = polptp(4) + polv_perm(ix,2) + polv_perm(ix + 3,2) + polv_perm(ix + 6,2) 
         polptp(5) = polptp(5) + polv_perm(iy,3) + polv_perm(iy + 3,3) + polv_perm(iy + 6,3) 
         polptp(6) = polptp(6) + polv_perm(iz,1) + polv_perm(iz + 3,1) + polv_perm(iz + 6,1) 
         polptp(7) = polptp(7) + polv_perm(iy,1) + polv_perm(iy + 3,1) + polv_perm(iy + 6,1) 
         polptp(8) = polptp(8) + polv_perm(iz,2) + polv_perm(iz + 3,2) + polv_perm(iz + 6,2) 
         polptp(9) = polptp(9) + polv_perm(ix,3) + polv_perm(ix + 3,3) + polv_perm(ix + 6,3)		 
        ENDDO


       DO i = 1, nmol
        j = 2*(i-1) + 1
        ix = 9*(i-1) + 1
        iy = ix + 1
        iz = ix + 2
        jx = ix + 3 
        jy = ix + 4 
        jz = ix + 5 
        kx = ix + 6 
        ky = ix + 7 
        kz = ix + 8 
        poltp(1) = poltp(1) + polv(9*nmol + j,1)*(cor(ix)-cor(kx)) &
        &                 + polv(9*nmol + j + 1,1)*(cor(jx)-cor(kx))
        poltp(2) = poltp(2) + polv(9*nmol + j,2)*(cor(iy)-cor(ky)) &
        &                 + polv(9*nmol + j + 1,2)*(cor(jy)-cor(ky))
        poltp(3) = poltp(3) + polv(9*nmol + j,3)*(cor(iz)-cor(kz)) &
        &                 + polv(9*nmol + j + 1,3)*(cor(jz)-cor(kz))
        poltp(7) = poltp(7) + polv(9*nmol + j,1)*(cor(iy)-cor(ky)) &
        &                 + polv(9*nmol + j + 1,1)*(cor(jy)-cor(ky))
        poltp(8) = poltp(8) + polv(9*nmol + j,2)*(cor(iz)-cor(kz)) &
        &                 + polv(9*nmol + j + 1,2)*(cor(jz)-cor(kz))
        poltp(9) = poltp(9) + polv(9*nmol + j,3)*(cor(ix)-cor(kx)) &
        &                 + polv(9*nmol + j + 1,3)*(cor(jx)-cor(kx))
        poltp(4) = poltp(4) + polv(9*nmol + j,2)*(cor(ix)-cor(kx)) &
        &                 + polv(9*nmol + j + 1,2)*(cor(jx)-cor(kx))
        poltp(5) = poltp(5) + polv(9*nmol + j,3)*(cor(iy)-cor(ky)) &
        &                 + polv(9*nmol + j + 1,3)*(cor(jy)-cor(ky))
        poltp(6) = poltp(6) + polv(9*nmol + j,1)*(cor(iz)-cor(kz)) &
        &                 + polv(9*nmol + j + 1,1)*(cor(jz)-cor(kz))
		
        polptp(1) = polptp(1) + polv_perm(9*nmol + j,1)*(cor(ix)-cor(kx)) &
        &                 + polv_perm(9*nmol + j + 1,1)*(cor(jx)-cor(kx))
        polptp(2) = polptp(2) + polv_perm(9*nmol + j,2)*(cor(iy)-cor(ky)) &
        &                 + polv_perm(9*nmol + j + 1,2)*(cor(jy)-cor(ky))
        polptp(3) = polptp(3) + polv_perm(9*nmol + j,3)*(cor(iz)-cor(kz)) &
        &                 + polv_perm(9*nmol + j + 1,3)*(cor(jz)-cor(kz))
        polptp(7) = polptp(7) + polv_perm(9*nmol + j,1)*(cor(iy)-cor(ky)) &
        &                 + polv_perm(9*nmol + j + 1,1)*(cor(jy)-cor(ky))
        polptp(8) = polptp(8) + polv_perm(9*nmol + j,2)*(cor(iz)-cor(kz)) &
        &                 + polv_perm(9*nmol + j + 1,2)*(cor(jz)-cor(kz))
        polptp(9) = polptp(9) + polv_perm(9*nmol + j,3)*(cor(ix)-cor(kx)) &
        &                 + polv_perm(9*nmol + j + 1,3)*(cor(jx)-cor(kx))
        polptp(4) = polptp(4) + polv_perm(9*nmol + j,2)*(cor(ix)-cor(kx)) &
        &                 + polv_perm(9*nmol + j + 1,2)*(cor(jx)-cor(kx))
        polptp(5) = polptp(5) + polv_perm(9*nmol + j,3)*(cor(iy)-cor(ky)) &
        &                 + polv_perm(9*nmol + j + 1,3)*(cor(jy)-cor(ky))
        polptp(6) = polptp(6) + polv_perm(9*nmol + j,1)*(cor(iz)-cor(kz)) &
        &                 + polv_perm(9*nmol + j + 1,1)*(cor(jz)-cor(kz))		
       ENDDO       
        poltp(4) = 0.5d0*(poltp(4) + poltp(7))
        poltp(5) = 0.5d0*(poltp(5) + poltp(8))
        poltp(6) = 0.5d0*(poltp(6) + poltp(9))

        polptp(4) = 0.5d0*(polptp(4) + polptp(7))
        polptp(5) = 0.5d0*(polptp(5) + polptp(8))
        polptp(6) = 0.5d0*(polptp(6) + polptp(9))
	   
        diptp(:) = (diptp(:) - dipt(:))/edt
        dipitp(:) = (dipitp(:) - dipit(:))/edt
		polptp(:) = (polptp(:) - polpt(:))/edt
        poltp(:) = (poltp(:) - polt(:))/edt
      ENDIF

      IF(fdpdt)THEN
        WRITE(39) dipt(:) , dipit(:)
		WRITE(39) diptp(:), dipitp(:)
		WRITE(19) polpt(1:6),(polt(1:6)-polpt(1:6))
		WRITE(19) polptp(1:6),(poltp(1:6)-polptp(1:6))
      ELSE
        WRITE(39) dipt(:) , dipit(:)
    	WRITE(19) polpt(1:6),(polt(1:6)-polpt(1:6))
      ENDIF
	  
      ndip = ndip + 1
	  
	  
	  
	 !Get Charge
	 CALL MPOUT(cp,ci)
	 qHperm = 0d0
	 qOperm = 0d0
	 qHind  = 0d0
	 qOind  = 0d0
	 
	 DO i=1,nmol
	  ix = 3*(i-1) + 1
	  iz = ix + 2
	  qHperm = qHperm + cp(ix) + cp(ix+1)
	  qOperm = qOperm + cp(iz)
      qHind  = qHind  + ci(ix) + ci(ix+1)
	  qOind  = qOind  + ci(iz)
	 END DO
	 
      qHperm = qHperm/dble(2*nmol)
	  qHind  = qHind/dble(2*nmol)
	  qOperm = qOperm/dble(nmol)
	  qOind  = qOind/dble(nmol)
	  
	  qHpav = qHpav + qHperm
	  qHiav = qHiav + qHind
	  qOpav = qOpav + qOperm
	  qOiav = qOiav + qOind
	  
    ENDIF !npol

	
	IF(ncha /= 0 .and. mod(n,ncha) == 0)THEN
	    CALL MPOUT(cp,ci)
        WRITE(59) cp,ci              ! charge - induced charge
	END IF
	
	
    IF(ntrj > 0 .and. mod(n,ntrj) == 0 )THEN
      cor_tmp = cor
      IF(ftrjmin == .TRUE.)THEN
        iter = 1000
        CALL FRPRMN(cor_tmp, 3*natom, 1d-5, iter, tmp)
        WRITE(9,*) n,'ENGMIN', tmp
      ENDIF

      IF(fbtrj)THEN
        WRITE(15) cor_tmp,ptt,SETPB(0d0)
      ELSE

        WRITE(15,'(3ES18.10)')&
       &(cor_tmp(3*(i-1)+1),cor_tmp(3*(i-1)+2),cor_tmp(3*(i-1)+3)&
       &    ,i=1,natom)
        IF(cjob == 'NPT')THEN
          WRITE(15,'(ES18.10)') SETPB(0d0)
        ENDIF

        WRITE(15,*)
      ENDIF
    ENDIF !ntrj
	
	
  
END SUBROUTINE GET_TRAJECTRY


SUBROUTINE GET_TIME	
IMPLICIT NONE

INTEGER :: nfin
INTEGER :: nnow

 CALL CPU_TIME( ct2 )

 ct=ct2-ct1 !CPU TIME (s)

 time_s=dint(ct)
 time_m=time_s/60
 time_h=time_m/60
 time_d=time_h/24

 time_ms=dint( (ct-time_s)*1d3 )
 time_s=mod(time_s,60)
 time_m=mod(time_m,60)
 time_h=mod(time_h,24)

 WRITE(*,*) "Simulation Time" 
 WRITE(*,"(a,i5,a,i5,a,i5,a,i5,a,i5,a)") &
  & " Time",time_d," d",time_h," h",time_m," m",time_s," s",time_ms," ms" 
 
 IF(EqMD)THEN
  nfin = maxstep
  nnow = n
 ELSE IF(Noneq.or.jobtype=="HESSIAN")THEN
  nfin = mstep
  nnow = m 
 ELSE
  nfin = maxstep
  nnow = n 
 END IF

 ct = ct2 - ct1 !CPU TIME (s)  
 ct = ct * dble( nfin - nnow ) / dble(nnow)
 time_s = dint(ct)
 time_m = time_s/60
 time_h = time_m/60
 time_d = time_h/24

 time_ms=dint( (ct-time_s)*1d3 )
 time_s=mod(time_s,60)
 time_m=mod(time_m,60)
 time_h=mod(time_h,24)

 WRITE(*,*) "Remaining Time"
 WRITE(*,"(a,i5,a,i5,a,i5,a,i5,a,i5,a)") &
  & " Time",time_d," d",time_h," h",time_m," m",time_s," s",time_ms," ms"  
  
 ct = ct2 - ct1 !CPU TIME (s) 
 ct = ct * dble( nfin ) / dble(nnow)
 time_s = dint(ct)
 time_m = time_s/60
 time_h = time_m/60
 time_d = time_h/24

 time_ms=dint( (ct-time_s)*1d3 )
 time_s=mod(time_s,60)
 time_m=mod(time_m,60)
 time_h=mod(time_h,24)

 WRITE(*,*) "Total Time"
 WRITE(*,"(a,i5,a,i5,a,i5,a,i5,a,i5,a)") &
  & " Time",time_d," d",time_h," h",time_m," m",time_s," s",time_ms," ms"  

 
 WRITE(*,*)
			 
END SUBROUTINE GET_TIME	
			 
			 
END PROGRAM WATER_1
