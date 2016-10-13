program cgrid_shallow_water
  !C-grid shallow water model written by Peter Dueben based on F77 code of David Marshall
  
  USE rp_emulator !Use edited emulator for reduced precision to mimic bit flips (see difference_emulator.txt)
  implicit none

!Grid parameters:  
  integer :: nx_f,ny_f,nx_c,ny_c,nx_i,ny_i,nt_f,nt_c  
  parameter(nx_f=181,ny_f=31) !Resolution fine grid
  parameter(nx_c=61,ny_c=11)  !Resolution coarse grid
  parameter(nx_i=181,ny_i=31) !Resolution input file
  parameter(nt_f=3,nt_c=3)    !Time stepping fine/coarse grid
 
!General parameters:
  integer :: i,j,k,l,n_c,n_f,nstop,itestcase, nnudge  !Time intervals for nudging in timesteps
  integer*8 :: k8,l8
  integer :: nsec,nwrite,ninter
  integer :: nxstep,nystep !Number of steps to jump for "on-screen" output
  character*10 :: cinput
  REAL*8 :: Lx,Ly !Domain size 
  REAL*8 :: slip,gp
  real*8 :: time1, time2, time3, time_c, time_f,rdummy
  logical :: lboundary   ! .T. if with boundaries.
  logical :: lrestart    ! .T. if restart
  REAL*8 :: random1, random2, random3
  REAL*8 :: pfault
  INTEGER :: ntime_f(0:nx_f,0:ny_f)
  INTEGER*8 :: ntest1, ntest2
  LOGICAL :: lerror, lbackup, ltopography
  REAL*8 :: hr_c(0:nx_c,0:ny_c,0:1),ur_c(0:nx_c,0:ny_c,0:1),vr_c(0:nx_c,0:ny_c,0:1)
  CHARACTER(len=70) :: output 

!Parameters coarse grid:
  REAL*8 :: h_c(0:nx_c,0:ny_c),u_c(0:nx_c,0:ny_c),v_c(0:nx_c,0:ny_c),taux_c(0:ny_c),tauy_c(0:nx_c) !     layer thickness (h), velocity components (u,v) and wind forcing
  REAL*8 :: ht_c(0:nx_c,0:ny_c),ut_c(0:nx_c,0:ny_c),vt_c(0:nx_c,0:ny_c) !Fields to transform between the coarse and the fine model.
  REAL*8 :: HC_c(0:nx_c,0:ny_c)  !Height of the fluid column
  REAL*8 :: dh_c(0:nx_c,0:ny_c,nt_c),du_c(0:nx_c,0:ny_c,nt_c),dv_c(0:nx_c,0:ny_c,nt_c) !Time increments for AB timestepping
  REAL*8 :: ab_c(nt_c) !AB coefficients
  REAL*8 :: dx_c,dy_c,dt_c,rdx_c,rdy_c
  REAL*8 :: fu_c(0:ny_c),fv_c(0:ny_c) !Coriolis parameter at u and v grid-points respectively
  REAL*8 :: b_c(0:nx_c,0:ny_c) ! Bernoulli potential and relative vorticity 
  REAL*8 :: nudge_c, nuu_c, nuh_c

!Parameters fine grid:
  REAL*8 :: h_f(0:nx_f,0:ny_f),u_f(0:nx_f,0:ny_f),v_f(0:nx_f,0:ny_f),taux_f(0:ny_f),tauy_f(0:nx_f) !     layer thickness (h), velocity components (u,v) and wind forcing
  REAL*8 :: ht_f(0:nx_f,0:ny_f),ut_f(0:nx_f,0:ny_f),vt_f(0:nx_f,0:ny_f) !Fields to transform between the coarse and the fine model.
  REAL*8 :: HC_f(0:nx_f,0:ny_f)  !Height of the fluid column
  REAL*8 :: dh_f(0:nx_f,0:ny_f,nt_f),du_f(0:nx_f,0:ny_f,nt_f),dv_f(0:nx_f,0:ny_f,nt_f) !Time increments for AB timestepping
  REAL*8 :: ab_f(nt_f) !AB coefficients
  REAL*8 :: dx_f,dy_f,dt_f,rdx_f,rdy_f
  REAL*8 :: fu_f(0:ny_f),fv_f(0:ny_f) !Coriolis parameter at u and v grid-points respectively
  REAL*8 :: b_f(0:nx_f,0:ny_f)! Bernoulli potential and relative vorticity 
  REAL*8 :: nudge_f, nuu_f, nuh_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Emulator  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RPE_ACTIVE = .FALSE.
  RPE_STOCHASTIC = .FALSE.
  rpe_fault_rate = 0.0000001_8 !For simulated bit flips
  lbackup = .TRUE.
  ltopography = .FALSE. !Allows reproduction of coarse topography if hardware faults are simulated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nstop =  100000      ! number of timesteps on coarse grid
  nwrite=    1000      ! number of timesteps between file output
  nxstep=4             !Advanced terminal output
  nystep=nxstep*2
  ninter = 1           !Ratio between coarse and fine timestep
  lrestart = .FALSE.   
  cinput = '600000'    ! Input number
  nudge_c =  0.0_8     ! Strength of nudging
  nudge_f =  0.0_8
  nnudge = 1           ! Timestep to update transformation fields
  itestcase = 5        ! Testcase: 1 = Gaussian hill (unstable), 2 = Stommel, 3 = Gaussian hill (stable), 4 = Random fields, 5 = Isolated mountain test
  slip=1._8            ! free-slip (0.) or no-slip (1.)?  
  nsec=0               ! Running output number
! Initialise model fields:
  CALL initialise(itestcase,lboundary,lrestart,Lx,Ly,gp,cinput,nx_i,ny_i,ninter,&
       & nx_f,ny_f,nt_f,nuh_f,nuu_f,dt_f,HC_f,&
       & dx_f,dy_f,rdx_f,rdy_f,ab_f,fu_f,fv_f,taux_f,tauy_f,h_f,dh_f,u_f,du_f,v_f,dv_f, &
       & nx_c,ny_c,nt_c,nuh_c,nuu_c,dt_c,HC_c,&
       & dx_c,dy_c,rdx_c,rdy_c,ab_c,fu_c,fv_c,taux_c,tauy_c,h_c,dh_c,u_c,du_c,v_c,dv_c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    MAIN LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 time_c = 0.0_8
 time_f = 0.0_8
 
 ut_c = 0.0_8
 vt_c = 0.0_8
 ht_c = 0.0_8
 ut_f = 0.0_8
 vt_f = 0.0_8
 ht_f = 0.0_8
 
 ntime_f(:,:) = 1
 ntest1 = 0

 CALL CPU_TIME(time1)

!Initialise backup system:
 IF(lbackup)THEN
    CALL finetocoarse(nx_c,ny_c,hr_c(:,:,0),ur_c(:,:,0),vr_c(:,:,0),&
         &nx_f,ny_f,h_f,u_f,v_f) 
    CALL finetocoarse(nx_c,ny_c,hr_c(:,:,1),ur_c(:,:,1),vr_c(:,:,1),&
         &nx_f,ny_f,h_f,u_f,v_f) 
 END IF


 DO n_c=1,nstop+1

    !Output model fields:
    RPE_ACTIVE = .FALSE.
    IF(mod(n_c,nwrite)==1.or.nwrite==1) CALL output_fields(lboundary,nsec,nxstep,nystep,&
         & u_c,du_c,v_c,dv_c,h_c,dh_c,dx_c,dy_c,nx_c,ny_c,nt_c,&
         & u_f,du_f,v_f,dv_f,h_f,dh_f,dx_f,dy_f,nx_f,ny_f,nt_f)
    
    RPE_ACTIVE = .TRUE.
    !Perform model timestep:
    CALL timeloop(n_c,lboundary,lrestart,gp,slip,&
         & nx_f,ny_f,nt_f,nuh_f,nuu_f,dt_f,&
         & rdx_f,rdy_f,ab_f,fu_f,fv_f,taux_f,tauy_f,&
         & h_f,dh_f,u_f,du_f,v_f,dv_f,ut_f,vt_f,ht_f,HC_f,nudge_f)  
    
    !This will mimic a hardware fault and set large parts of the domain to NAN:
    IF(n_c==9000.or.n_c==19000.or.n_c==49000.or.n_c==99000)THEN
       k8 = HUGE(k8)
       rdummy = TRANSFER(k8,rdummy)
       k = (nx_f-1)/2  !The definition of k and l will define the part of the domain that is set to NAN
       l = (ny_f-1)
       DO i=1,k
          DO j=1,l
             h_f(i,j) = rdummy
             u_f(i,j) = rdummy
             v_f(i,j) = rdummy
             dh_f(i,j,1) = rdummy
             dh_f(i,j,2) = rdummy
             dh_f(i,j,3) = rdummy
             du_f(i,j,1) = rdummy
             du_f(i,j,2) = rdummy
             du_f(i,j,3) = rdummy
             dv_f(i,j,1) = rdummy
             dv_f(i,j,2) = rdummy
             dv_f(i,j,3) = rdummy
          END DO
       END DO
    END IF

    !This will run the backup system:
    IF(lbackup)THEN
       call fault_secure(n_c,nx_c,ny_c,hr_c,ur_c,vr_c,HC_c,nx_f,ny_f,h_f,u_f,v_f,HC_f,ntime_f,ltopography)
       u_c = ur_c(:,:,mod(n_c,2))
       v_c = vr_c(:,:,mod(n_c,2))
       h_c = hr_c(:,:,mod(n_c,2))
    END IF
    
    nsec =nsec + 1
   
 END DO
 
 CALL CPU_TIME(time2)
 write(*,*) 'Time model run: ', time2-time1
 
end program cgrid_shallow_water



subroutine fault_secure(n_c,nx_c,ny_c,h_c,u_c,v_c,HC_c,&
          &nx_f,ny_f,h_f,u_f,v_f,HC_f,ntime_f,ltopography)

  !This is the backup system

  USE rp_emulator
  implicit none

  INTEGER :: nx_c, ny_c, n_c
  REAL*8 :: u_c(0:nx_c,0:ny_c,0:1), v_c(0:nx_c,0:ny_c,0:1), h_c(0:nx_c,0:ny_c,0:1), HC_c(0:nx_c,0:ny_c)

  INTEGER :: nx_f, ny_f
  REAL*8 :: u_f(0:nx_f,0:ny_f), v_f(0:nx_f,0:ny_f), h_f(0:nx_f,0:ny_f), HC_f(0:nx_f,0:ny_f)
  INTEGER :: i,j,k,l
  INTEGER :: ntime_f(0:nx_f,0:ny_f)

  REAL*8 :: gf(3,3),gc(2,2),tc,uc
 
  REAL*8 :: ut_f(0:nx_f,0:ny_f), vt_f(0:nx_f,0:ny_f), ht_f(0:nx_f,0:ny_f), HC2_f(0:nx_f,0:ny_f)
  logical :: lerror
  logical :: hl(0:nx_c,0:ny_c),ul(0:nx_c,0:ny_c),vl(0:nx_c,0:ny_c),ltopography

  character(len=100) filename

  !Flag that defines whether a hardware fault was found for a specific parameter
  DO i=1,nx_c-1
     Do j=1,ny_c-1
        hl(i,j) = .TRUE.
        ul(i,j) = .TRUE.
        vl(i,j) = .TRUE.
     END DO
  END DO

  !This will map the prognostic fields from the model to the backup grid:
  CALL finetocoarse(nx_c,ny_c,h_c(:,:,mod(n_c,2)),u_c(:,:,mod(n_c,2)),v_c(:,:,mod(n_c,2)),&
       &nx_f,ny_f,h_f,u_f,v_f)  
  
  !This will test whether the parameter on the backup grid has changed too much:
  DO i=1,nx_c-1
     Do j=1,ny_c-1
        IF(abs(h_c(i,j,0)-h_c(i,j,1)).gt.0.05_8.or.h_c(i,j,mod(n_c,2)).ne.h_c(i,j,mod(n_c,2)))THEN !Put in flow depebent threshold here for h, second part tests for NAN
           h_c(i,j,mod(n_c,2)) = h_c(i,j,mod(n_c+1,2))
           hl(i,j) = .FALSE.
           hl(i-1,j) = .FALSE.
           hl(i,j-1) = .FALSE.
           hl(i-1,j-1) = .FALSE.
           IF(i==1)THEN
              hl(nx_c-1,j) = .FALSE.
              hl(nx_c-1,j-1) = .FALSE.
           END IF
           IF(j==1)THEN
              hl(i,ny_c-1) = .FALSE.
              hl(i-1,ny_c-1) = .FALSE.
           END IF
           IF(i==1.and.j==1)THEN
              hl(nx_c-1,ny_c-1) = .FALSE.
           END IF
        END IF
        IF(abs(u_c(i,j,0)-u_c(i,j,1)).gt.0.01_8.or.u_c(i,j,mod(n_c,2)).ne.u_c(i,j,mod(n_c,2)))THEN !Put in flow depebent threshold here for u, second part tests for NAN
           u_c(i,j,mod(n_c,2)) = u_c(i,j,mod(n_c+1,2))
           ul(i,j) = .FALSE.
           ul(i-1,j) = .FALSE.
           ul(i,j-1) = .FALSE.
           ul(i-1,j-1) = .FALSE.
           IF(i==1)THEN
              ul(nx_c-1,j) = .FALSE.
              ul(nx_c-1,j-1) = .FALSE.
           END IF
           IF(j==1)THEN
              ul(i,ny_c-1) = .FALSE.
              ul(i-1,ny_c-1) = .FALSE.
           END IF
           IF(i==1.and.j==1)THEN
              ul(nx_c-1,ny_c-1) = .FALSE.
           END IF
        END IF
        IF(abs(v_c(i,j,0)-v_c(i,j,1)).gt.0.01_8.or.v_c(i,j,mod(n_c,2)).ne.v_c(i,j,mod(n_c,2)))THEN !Put in flow depebent threshold here for v, second part tests for NAN
           v_c(i,j,mod(n_c,2)) = v_c(i,j,mod(n_c+1,2))
           vl(i,j) = .FALSE.
           vl(i-1,j) = .FALSE.
           vl(i,j-1) = .FALSE.
           vl(i-1,j-1) = .FALSE.
           IF(i==1)THEN
              vl(nx_c-1,j) = .FALSE.
              vl(nx_c-1,j-1) = .FALSE.
           END IF
           IF(j==1)THEN
              vl(i,ny_c-1) = .FALSE.
              vl(i-1,ny_c-1) = .FALSE.
           END IF
           IF(i==1.and.j==1)THEN
              vl(nx_c-1,ny_c-1) = .FALSE.
           END IF    
        END IF
     END DO
  END DO
  
  !Care for periodicity:
  do j=1,ny_c-1
     u_c(nx_c,j,mod(n_c,2)) = u_c(1,j,mod(n_c,2))
     v_c(nx_c,j,mod(n_c,2)) = v_c(1,j,mod(n_c,2))
     h_c(nx_c,j,mod(n_c,2))= h_c(1,j,mod(n_c,2))
  end do
  do j=1,nx_c-1
     u_c(j,ny_c,mod(n_c,2)) = u_c(j,1,mod(n_c,2))
     v_c(j,ny_c,mod(n_c,2)) = v_c(j,1,mod(n_c,2))
     h_c(j,ny_c,mod(n_c,2))=h_c(j,1,mod(n_c,2))
  end do
  u_c(nx_c,ny_c,mod(n_c,2))=u_c(1,1,mod(n_c,2))
  v_c(nx_c,ny_c,mod(n_c,2))=v_c(1,1,mod(n_c,2))
  h_c(nx_c,ny_c,mod(n_c,2))=h_c(1,1,mod(n_c,2))


  !If hardware fault was detected, test whether the corresponding values on the model grid have resonable values:
  DO i=1,nx_c-1
     Do j=1,ny_c-1
        IF(.not.hl(i,j))THEN 
           gc(1,1)=h_c(i,j,mod(n_c,2))
           gc(2,1)=h_c(i+1,j,mod(n_c,2))
           gc(1,2)=h_c(i,j+1,mod(n_c,2))
           gc(2,2)=h_c(i+1,j+1,mod(n_c,2))
           CALL map_fine_cell(gf,gc) 
           DO k=1,3
              DO l=1,3
                 IF(abs(h_f((i-1)*3+k,(j-1)*3+l)).ge.8.0_8.or.h_f((i-1)*3+k,(j-1)*3+l).ne.h_f((i-1)*3+k,(j-1)*3+l)) h_f((i-1)*3+k,(j-1)*3+l) = gf(k,l) !Put in flow dependent limites for h here
             END DO
           END DO
        END IF
          
        IF(.not.ul(i,j))THEN
           gc(1,1)=u_c(i,j,mod(n_c,2))
           gc(2,1)=u_c(i+1,j,mod(n_c,2))
           gc(1,2)=u_c(i,j+1,mod(n_c,2))
           gc(2,2)=u_c(i+1,j+1,mod(n_c,2))
           CALL map_fine_cell(gf,gc) 
           DO k=1,3
              DO l=1,3
                 IF(u_f((i-1)*3+k-1,(j-1)*3+l).ne.u_f((i-1)*3+k-1,(j-1)*3+l).or.abs(u_f((i-1)*3+k-1,(j-1)*3+l)-10.0_8).ge.2.0_8) u_f((i-1)*3+k-1,(j-1)*3+l) = gf(k,l) !Put in flow dependent limites for u here
              END DO
           END DO
        END IF
 
        IF(.not.vl(i,j))THEN
           gc(1,1)=v_c(i,j,mod(n_c,2))
           gc(2,1)=v_c(i+1,j,mod(n_c,2))
           gc(1,2)=v_c(i,j+1,mod(n_c,2))
           gc(2,2)=v_c(i+1,j+1,mod(n_c,2))
           CALL map_fine_cell(gf,gc) 
           DO k=1,3
              DO l=1,3
                 IF(v_f((i-1)*3+k,(j-1)*3+l-1).ne.v_f((i-1)*3+k,(j-1)*3+l-1).or.abs(v_f((i-1)*3+k,(j-1)*3+l-1)).ge.1.0_8) v_f((i-1)*3+k,(j-1)*3+l-1) = gf(k,l) !Put in flow dependent limites for v here
              END DO
           END DO
        END IF
     END DO
  END DO

  !Take care for periodicity
  DO j=1,ny_f-1
     u_f(nx_f-1,j) = u_f(0,j)
  END DO
  
  DO i=1,nx_f-1
     v_f(i,ny_f-1) = v_f(i,0)
  END DO

  !Change topography if the coarse value should be restored:
  IF(ltopography)THEN
     lerror = .FALSE.
     LOOP1: DO i=1,nx_f-1
        LOOP2: DO j=1,ny_f-1
           IF(HC_f(i,j).ne.HC_f(i,j))lerror = .TRUE.
           EXIT LOOP1
        END DO LOOP2
     END DO LOOP1

     IF(lerror)THEN
        CALL coarsetofine_h(nx_c,ny_c,HC_c,nx_f,ny_f,HC2_f)
        DO i=1,nx_f-1
           DO j=1,ny_f-1
              IF(abs(HC_f(i,j)).gt.1000.0_8.or.HC_f(i,j).ne.HC_f(i,j))THEN
                 HC_f(i,j) = HC2_f(i,j)
              END IF
           END DO
        END DO
     END IF
  END IF

  return
end subroutine fault_secure



!    --------------------------------------------------------------------------


subroutine timeloop(n,lboundary,lrestart,gp,slip,&
       & nx,ny,nt,nuh,nuu,dt,&
       & rdx,rdy,ab,fu,fv,taux,tauy,&
       & h,dh,u,du,v,dv,ut,vt,ht,HC,nudge)


  !Perform model timestep
  
  USE rp_emulator
  implicit none

  INTEGER :: i,j,n,nx,ny,nt,k,l
  logical :: lboundary, lrestart
  REAL*8 :: gp,nuu,nuh,dt,rdx,rdy,slip

  REAL*8 :: h(0:nx,0:ny),u(0:nx,0:ny),v(0:nx,0:ny),taux(0:ny),tauy(0:nx)
  REAL*8 :: ht(0:nx,0:ny),ut(0:nx,0:ny),vt(0:nx,0:ny)
  REAL*8 :: HC(0:nx,0:ny)
  REAL*8 :: dh(0:nx,0:ny,nt),du(0:nx,0:ny,nt),dv(0:nx,0:ny,nt)
  REAL*8 :: ab(nt)
  REAL*8 :: fu(0:ny),fv(0:ny)
  REAL*8 :: b(0:nx,0:ny),zeta(0:nx,0:ny)
  REAL*8 :: nudge
  REAL*8 :: td125,td25,td5,t1,t2

  !Define constants are reals to make sure the emulator is working correctly
  td125=0.125_8
  td25 = 0.25_8
  td5 = 0.5_8
  t1 = 1.0_8
  t2 = 2.0_8

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the grid:
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------...
!  h(1,ny-1)   u(2,ny-1)   h(2,ny-1)   u(3,ny-1)   h(3,ny-1)   u(4,ny-1)       ... u(nx-1,ny-1)   h(nx-1,ny-1)   
!
!--v(1,ny-1)---------------v(1,ny-1)---------------v(1,ny-1)-...                                  v(nx-1,ny-1) 
!
!  h(1,ny-2)   u(2,ny-2)   h(2,ny-2)   u(3,ny-2)   h(3,ny-2)   u(4,ny-2)       ... u(nx-1,ny-2)   h(nx-1,ny-2)   
!.
!.
!.
!----------------------------------------------...
!
!  h(1,2)   u(2,2)   h(2,2)   u(3,2)   h(3,2)   u(4,2) ... u(nx-1,2)   h(nx-1,2)
!
!--v(1,2))-----------v(2,2)------------v(3,2)--...                     v(nx-1,2)
!
!  h(1,1)   u(2,1)   h(2,1)   u(3,1)   h(3,1)   u(4,1) ... u(nx-1,1)   h(nx-1,1)  
!------------------------------------------------------

  
  !Set condition for periodicity:
  IF(.not.lboundary)THEN
     do j=1,ny-1
        u(0,j) = u(nx-1,j)
        u(nx,j) = u(1,j)
        v(0,j) = v(nx-1,j)
        v(nx,j) = v(1,j)
        h(0,j)=h(nx-1,j)
        h(nx,j)=h(1,j)
     end do
     do j=1,nx-1
        u(j,0) = u(j,ny-1)
        u(j,ny) = u(j,1)
        v(j,0) = v(j,ny-1)
        v(j,ny) = v(j,1)
        h(j,0)=h(j,ny-1)
        h(j,ny)=h(j,1)
     end do
     u(0,0)=u(nx-1,ny-1)
     v(0,0)=v(nx-1,ny-1)
     h(0,0)=h(nx-1,ny-1)
     u(nx,ny)=u(1,1)
     v(nx,ny)=v(1,1)
     h(nx,ny)=h(1,1)
     u(0,ny)=u(nx-1,1)
     v(0,ny)=v(nx-1,1)
     h(0,ny)=h(nx-1,1)
     u(nx,0)=u(1,ny-1)
     v(nx,0)=v(1,ny-1)
     h(nx,0)=h(1,ny-1)
  END IF

  ! calculate Bernoulli potential
  ! 1/g*h+0.5(u^2+v^2)
  do j=0,ny-1
     do i=0,nx-1
        b(i,j)=gp*h(i,j)+td125*&
             & ((u(i,j)+u(i+1,j))**t2+(v(i,j)+v(i,j+1))**t2)
     end do
  end do
  
  ! calculate relative vorticity
  ! d_f v - d_y u
  do j=1,ny
     do i=1,nx
        zeta(i,j)=(v(i,j)-v(i-1,j))*rdx-(u(i,j)-u(i,j-1))*rdy
     end do
  end do
  
  !    calculate forcing of u,v, and h
  do j=1,ny-1
     do i=2,nx-1
        du(i,j,3)=du(i,j,2)  !For Adams Bashforth
        du(i,j,2)=du(i,j,1)  !For Adams Bashforth
        du(i,j,1)= nuu*(u(i+1,j)+u(i-1,j)-t2*u(i,j))*rdx**t2+nuu*(u(i,j+1)+u(i,j-1)-t2*u(i,j))*rdy**t2 &  !\nu (dx dx u + dy dy u)
             & +td25*(fu(j)+td5*(zeta(i,j)+zeta(i,j+1)))*(v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1)) & ! +(f+zeta)v
             & -(b(i,j)-b(i-1,j))*rdx + taux(j) + nudge*(ut(i,j)-u(i,j)) !-dx b = dx (g*h+0.5(u^2+v^2))
     end do
  end do
  do j=2,ny-1
     do i=1,nx-1
        dv(i,j,3)=dv(i,j,2)  !For Adams Bashforth
        dv(i,j,2)=dv(i,j,1)  !For Adams Bashforth
        dv(i,j,1)= nuu*(v(i+1,j)+v(i-1,j)-t2*v(i,j))*rdx**t2+nuu*(v(i,j+1)+v(i,j-1)-t2*v(i,j))*rdy**t2 & !\nu (dx dx v + dy dy v)
             & -td25*(fv(j)+td5*(zeta(i,j)+zeta(i+1,j)))*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j)) & ! -(f+zeta)u
             & - (b(i,j)-b(i,j-1))*rdy + tauy(i) + nudge*(vt(i,j)-v(i,j)) !-dx b = dx (g*h+0.5(u^2+v^2))
     end do
  end do
  
  IF(.not.lboundary)THEN
     i=1
     do j=1,ny-1
        du(i,j,3)=du(i,j,2)  !For Adams Bashforth
        du(i,j,2)=du(i,j,1)  !For Adams Bashforth
        du(i,j,1)=nuu*(u(i+1,j)+u(i-1,j)-t2*u(i,j))*rdx**t2+nuu*(u(i,j+1)+u(i,j-1)-t2*u(i,j))*rdy**t2 &  !\nu (dx dx u + dy dy u)
             & +td25*(fu(j)+td5*(zeta(i,j)+zeta(i,j+1)))*(v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1)) & ! +(f+zeta)v
             & -(b(i,j)-b(i-1,j))*rdx + taux(j) + nudge*(ut(i,j)-u(i,j)) !-dx b = dx (g*h+0.5(u^2+v^2))        
     end do
     j=1
     do i=1,nx-1
        dv(i,j,3)=dv(i,j,2)  !For Adams Bashforth
        dv(i,j,2)=dv(i,j,1)  !For Adams Bashforth
        dv(i,j,1)=nuu*(v(i+1,j)+v(i-1,j)-t2*v(i,j))*rdx**t2+nuu*(v(i,j+1)+v(i,j-1)-t2*v(i,j))*rdy**t2 & !\nu (dx dx v + dy dy v)
             & -td25*(fv(j)+td5*(zeta(i,j)+zeta(i+1,j)))*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j)) & ! -(f+zeta)u
             & - (b(i,j)-b(i,j-1))*rdy + tauy(i) + nudge*(vt(i,j)-v(i,j)) !-dx b = dx (g*h+0.5(u^2+v^2))
     end do
  END IF
  
  do j=1,ny-1
     do i=1,nx-1
        dh(i,j,3)=dh(i,j,2)
        dh(i,j,2)=dh(i,j,1)
        dh(i,j,1)= & !nuh*(h(i+1,j)+h(i-1,j)-2.*h(i,j))/dx**2 +nuh*(h(i,j+1)+h(i,j-1)-2.*h(i,j))/dy**2 & !\nu_h (dx dx h + dy dy h)
             & +(td5*(HC(i-1,j)+h(i-1,j)+HC(i,j)+h(i,j)))*u(i,j)*rdx &
             & +(td5*(HC(i,j-1)+h(i,j-1)+HC(i,j)+h(i,j)))*v(i,j)*rdy &
             & -(td5*(HC(i+1,j)+h(i+1,j)+HC(i,j)+h(i,j)))*u(i+1,j)*rdx &
             & -(td5*(HC(i,j+1)+h(i,j+1)+HC(i,j)+h(i,j)))*v(i,j+1)*rdy + nudge*(ht(i,j)-h(i,j)) + nudge*(ht(i,j)-h(i,j))
     end do
  end do
  
  !    step forward for u,v, and h
  do j=1,ny-1
     do i=2,nx-1
        IF(n.lt.3.and.(.not.lrestart))THEN
           u(i,j)=u(i,j)+du(i,j,1)*dt
        ELSE
           u(i,j)=u(i,j)+ab(1)*du(i,j,1)+ab(2)*du(i,j,2)+ab(3)*du(i,j,3)
        END IF
     end do
  end do
  do j=2,ny-1
     do i=1,nx-1
        IF(n.lt.3.and.(.not.lrestart))THEN
           v(i,j)=v(i,j)+dv(i,j,1)*dt
        ELSE
           v(i,j)=v(i,j)+ab(1)*dv(i,j,1)+ab(2)*dv(i,j,2)+ab(3)*dv(i,j,3)
        END IF
     end do
  end do
  
  IF(.not.lboundary)THEN
     i=1
     do j=1,ny-1
        IF(n.lt.3.and.(.not.lrestart))THEN
           u(i,j)=u(i,j)+du(i,j,1)*dt
        ELSE
           u(i,j)=u(i,j)+ab(1)*du(i,j,1)+ab(2)*du(i,j,2)+ab(3)*du(i,j,3)
        END IF
     end do
     j=1
     do i=1,nx-1
        IF(n.lt.3.and.(.not.lrestart))THEN
           v(i,j)=v(i,j)+dv(i,j,1)*dt
        ELSE
           v(i,j)=v(i,j)+ab(1)*dv(i,j,1)+ab(2)*dv(i,j,2)+ab(3)*dv(i,j,3)
        END IF
     end do
  END IF
  
  do j=1,ny-1
     do i=1,nx-1
        IF(n.lt.3.and.(.not.lrestart))THEN
           h(i,j)=h(i,j)+dh(i,j,1)*dt
        ELSE
           h(i,j)=h(i,j)+ab(1)*dh(i,j,1)+ab(2)*dh(i,j,2)+ab(3)*dh(i,j,3)
        END IF
     end do
  end do
  
  ! evaluate dummy array elements from boundary conditions
  IF(lboundary)THEN
     do j=1,ny-1
        v(0,j)=(t1-t2*slip)*v(1,j)
        v(nx,j)=(t1-t2*slip)*v(nx-1,j)
        h(0,j)=h(1,j)
        h(nx,j)=h(nx-1,j)
     end do
     do i=1,nx-1
        u(i,0)=(t1-t2*slip)*u(i,1)
        u(i,ny)=(t1-t2*slip)*u(i,ny-1)
        h(i,0)=h(i,1)
        h(i,ny)=h(i,ny-1)
     end do
  END IF
  
  !Set condition for periodicity:
  IF(.not.lboundary)THEN
     do j=1,ny-1
        u(0,j) = u(nx-1,j)
        u(nx,j) = u(1,j)
        v(0,j) = v(nx-1,j)
        v(nx,j) = v(1,j)
        h(0,j)=h(nx-1,j)
        h(nx,j)=h(1,j)
     end do
     do j=1,nx-1
        u(j,0) = u(j,ny-1)
        u(j,ny) = u(j,1)
        v(j,0) = v(j,ny-1)
        v(j,ny) = v(j,1)
        h(j,0)=h(j,ny-1)
        h(j,ny)=h(j,1)
     end do
     u(0,0)=u(nx-1,ny-1)
     v(0,0)=v(nx-1,ny-1)
     h(0,0)=h(nx-1,ny-1)
     u(nx,ny)=u(1,1)
     v(nx,ny)=v(1,1)
     h(nx,ny)=h(1,1)
     u(0,ny)=u(nx-1,1)
     v(0,ny)=v(nx-1,1)
     h(0,ny)=h(nx-1,1)
     u(nx,0)=u(1,ny-1)
     v(nx,0)=v(1,ny-1)
     h(nx,0)=h(1,ny-1)
  END IF
  
end subroutine timeloop

 
subroutine finetocoarse(nx_c,ny_c,h_c,u_c,v_c,&
     &nx_f,ny_f,h_f,u_f,v_f)
  
  !Perform the mapping between the coarse and the fine grid

  USE rp_emulator
  implicit none
  
  INTEGER :: nx_c, ny_c
  REAL*8 :: u_c(0:nx_c,0:ny_c), v_c(0:nx_c,0:ny_c), h_c(0:nx_c,0:ny_c)

  INTEGER :: nx_f, ny_f
  REAL*8 :: u_f(0:nx_f,0:ny_f), v_f(0:nx_f,0:ny_f), h_f(0:nx_f,0:ny_f)
  INTEGER :: i,j,k,l
  REAL*8 :: t1,t6,t8,t9,t12
  t1 = 1.0_8
  t6 = 6.0_8
  t8 = 8.0_8
  t9 = 9.0_8
  t12 = 12.0_8

  !Height:
  DO j=1,ny_c-1
     DO i=1,nx_c-1
        h_c(i,j) = t1/t6*h_f((i-1)*3+1,(j-1)*3+1)&
             &  + t1/t8*(h_f((i-1)*3+1,(j-1)*3  )+h_f((i-1)*3,(j-1)*3+1)+    h_f((i-1)*3+2,(j-1)*3+1)+h_f((i-1)*3+1,(j-1)*3+2))&
             &  + t1/t12*(h_f((i-1)*3,(j-1)*3  )+h_f((i-1)*3+2,(j-1)*3  )+h_f((i-1)*3,(j-1)*3+2)+h_f((i-1)*3+2,(j-1)*3+2)) 
     END DO
  END DO
  
  !Zonal velocity:
  DO j=1,ny_c-1
     u_c(1,j) = t1/t6*(u_f(nx_f-1,(j-1)*3+1)) &
          &   + t1/t8*(u_f(nx_f-1,(j-1)*3  )+ u_f(nx_f-2,(j-1)*3+1)+u_f(nx_f,(j-1)*3+1)+u_f(nx_f-1,(j-1)*3+2))&
          &   +t1/t12*(u_f(nx_f-2,(j-1)*3  )+ u_f(nx_f,(j-1)*3)+u_f(nx_f-2,(j-1)*3+2)+u_f(nx_f,(j-1)*3+2))
     DO i=2,nx_c-1
        u_c(i,j) =  t1/t6*u_f((i-1)*3  ,(j-1)*3+1) &
             & + t1/t8*(u_f((i-1)*3  ,(j-1)*3  )+u_f((i-1)*3-1,(j-1)*3+1)+u_f((i-1)*3+1,(j-1)*3+1)+u_f((i-1)*3  ,(j-1)*3+2))&
             & +t1/t12*(u_f((i-1)*3-1,(j-1)*3  )+u_f((i-1)*3+1,(j-1)*3)+u_f((i-1)*3-1,(j-1)*3+2)+u_f((i-1)*3+1,(j-1)*3+2))
     END DO
  END DO
  
  !Vertical velocity:
  DO i=1,nx_c-1
     v_c(i,1) = t1/t6*v_f((i-1)*3+1,ny_f-1)&
          & + t1/t8*(v_f((i-1)*3+1,ny_f-2)+v_f((i-1)*3  ,ny_f-1)+v_f((i-1)*3+2,ny_f-1)+v_f((i-1)*3+1,ny_f  ))&
          & +t1/t12*(v_f((i-1)*3  ,ny_f-2)+v_f((i-1)*3+2,ny_f-2)+v_f((i-1)*3  ,ny_f  )+v_f((i-1)*3+2,ny_f  ))
  END DO
  DO j=2,ny_c-1
     DO i=1,nx_c-1
        v_c(i,j) = t1/t6*v_f((i-1)*3+1,(j-1)*3  )&
             & + t1/t8*(v_f((i-1)*3+1,(j-1)*3-1)+v_f((i-1)*3  ,(j-1)*3  )+v_f((i-1)*3+2,(j-1)*3  )+v_f((i-1)*3+1,(j-1)*3+1))&
             & +t1/t12*(v_f((i-1)*3  ,(j-1)*3-1)+v_f((i-1)*3+2,(j-1)*3-1)+v_f((i-1)*3  ,(j-1)*3+1)+v_f((i-1)*3+2,(j-1)*3+1))
     END DO
  END DO

  return
end subroutine finetocoarse


subroutine coarsetofine(nx_c,ny_c,h_c,u_c,v_c,&
     &nx_f,ny_f,h_f,u_f,v_f)

  !Mapping between the coarse and the fine grid:

  USE rp_emulator
  implicit none

  INTEGER :: nx_c, ny_c
  REAL*8 :: u_c(0:nx_c,0:ny_c), v_c(0:nx_c,0:ny_c), h_c(0:nx_c,0:ny_c)

  INTEGER :: nx_f, ny_f
  REAL*8 :: u_f(0:nx_f,0:ny_f), v_f(0:nx_f,0:ny_f), h_f(0:nx_f,0:ny_f)
  INTEGER :: i,j,k,l

  REAL*8 :: gf(3,3),gc(2,2),tc,uc

  !Height:
  DO j=1,ny_c-2
     DO i=1,nx_c-2
        gc(1,1)=h_c(i,j)
        gc(2,1)=h_c(i+1,j)
        gc(1,2)=h_c(i,j+1)
        gc(2,2)=h_c(i+1,j+1)
        CALL map_fine_cell(gf,gc) 
        DO k=1,3
           DO l=1,3
              h_f((i-1)*3+1,(j-1)*3+1) = gf(1,1)
              h_f((i-1)*3+2,(j-1)*3+1) = gf(2,1)
              h_f((i-1)*3+3,(j-1)*3+1) = gf(3,1)
              h_f((i-1)*3+1,(j-1)*3+2) = gf(1,2)
              h_f((i-1)*3+2,(j-1)*3+2) = gf(2,2)
              h_f((i-1)*3+3,(j-1)*3+2) = gf(3,2)
              h_f((i-1)*3+1,(j-1)*3+3) = gf(1,3)
              h_f((i-1)*3+2,(j-1)*3+3) = gf(2,3)
              h_f((i-1)*3+3,(j-1)*3+3) = gf(3,3)
           END DO
        END DO
     END DO
  END DO
  
  DO i=1,nx_c-2
     gc(1,1)=h_c(i,ny_c-1)
     gc(2,1)=h_c(i+1,ny_c-1)
     gc(1,2)=h_c(i,1)
     gc(2,2)=h_c(i+1,1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           h_f((i-1)*3+1,(ny_f-4)+1) = gf(1,1)
           h_f((i-1)*3+2,(ny_f-4)+1) = gf(2,1)
           h_f((i-1)*3+3,(ny_f-4)+1) = gf(3,1)
           h_f((i-1)*3+1,(ny_f-4)+2) = gf(1,2)
           h_f((i-1)*3+2,(ny_f-4)+2) = gf(2,2)
           h_f((i-1)*3+3,(ny_f-4)+2) = gf(3,2)
           h_f((i-1)*3+1,(ny_f-4)+3) = gf(1,3)
           h_f((i-1)*3+2,(ny_f-4)+3) = gf(2,3)
           h_f((i-1)*3+3,(ny_f-4)+3) = gf(3,3)
        END DO
     END DO
  END DO
  
  DO j=1,ny_c-2
     gc(1,1)=h_c(nx_c-1,j)
     gc(2,1)=h_c(1,j)
     gc(1,2)=h_c(nx_c-1,j+1)
     gc(2,2)=h_c(1,j+1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           h_f(nx_f-4+1,(j-1)*3+1) = gf(1,1)
           h_f(nx_f-4+2,(j-1)*3+1) = gf(2,1)
           h_f(nx_f-4+3,(j-1)*3+1) = gf(3,1)
           h_f(nx_f-4+1,(j-1)*3+2) = gf(1,2)
           h_f(nx_f-4+2,(j-1)*3+2) = gf(2,2)
           h_f(nx_f-4+3,(j-1)*3+2) = gf(3,2)
           h_f(nx_f-4+1,(j-1)*3+3) = gf(1,3)
           h_f(nx_f-4+2,(j-1)*3+3) = gf(2,3)
           h_f(nx_f-4+3,(j-1)*3+3) = gf(3,3)
        END DO
     END DO
  END DO
  
  gc(1,1)=h_c(nx_c-1,ny_c-1)
  gc(2,1)=h_c(1,ny_c-1)
  gc(1,2)=h_c(nx_c-1,1)
  gc(2,2)=h_c(1,1)
  CALL map_fine_cell(gf,gc) 
  DO k=1,3
     DO l=1,3
        h_f(nx_f-4+1,ny_f-4+1) = gf(1,1)
        h_f(nx_f-4+2,ny_f-4+1) = gf(2,1)
        h_f(nx_f-4+3,ny_f-4+1) = gf(3,1)
        h_f(nx_f-4+1,ny_f-4+2) = gf(1,2)
        h_f(nx_f-4+2,ny_f-4+2) = gf(2,2)
        h_f(nx_f-4+3,ny_f-4+2) = gf(3,2)
        h_f(nx_f-4+1,ny_f-4+3) = gf(1,3)
        h_f(nx_f-4+2,ny_f-4+3) = gf(2,3)
        h_f(nx_f-4+3,ny_f-4+3) = gf(3,3)
     END DO
  END DO
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Zonal velocity: 
  DO j=1,ny_c-2
     DO i=1,nx_c-2
        gc(1,1)=u_c(i,j)
        gc(2,1)=u_c(i+1,j)
        gc(1,2)=u_c(i,j+1)
        gc(2,2)=u_c(i+1,j+1)
        CALL map_fine_cell(gf,gc) 
        DO k=1,3
           DO l=1,3
              u_f((i-1)*3,(j-1)*3+1)   = gf(1,1)
              u_f((i-1)*3+1,(j-1)*3+1) = gf(2,1)
              u_f((i-1)*3+2,(j-1)*3+1) = gf(3,1)
              u_f((i-1)*3,(j-1)*3+2)   = gf(1,2)
              u_f((i-1)*3+1,(j-1)*3+2) = gf(2,2)
              u_f((i-1)*3+2,(j-1)*3+2) = gf(3,2)
              u_f((i-1)*3,(j-1)*3+3)   = gf(1,3)
              u_f((i-1)*3+1,(j-1)*3+3) = gf(2,3)
              u_f((i-1)*3+2,(j-1)*3+3) = gf(3,3)
           END DO
        END DO
     END DO
  END DO
  
  
  DO i=1,nx_c-2
     gc(1,1)=u_c(i,ny_c-1)
     gc(2,1)=u_c(i+1,ny_c-1)
     gc(1,2)=u_c(i,1)
     gc(2,2)=u_c(i+1,1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           u_f((i-1)*3  ,(ny_f-4)+1) = gf(1,1)
           u_f((i-1)*3+1,(ny_f-4)+1) = gf(2,1)
           u_f((i-1)*3+2,(ny_f-4)+1) = gf(3,1)
           u_f((i-1)*3  ,(ny_f-4)+2) = gf(1,2)
           u_f((i-1)*3+1,(ny_f-4)+2) = gf(2,2)
           u_f((i-1)*3+2,(ny_f-4)+2) = gf(3,2)
           u_f((i-1)*3  ,(ny_f-4)+3) = gf(1,3)
           u_f((i-1)*3+1,(ny_f-4)+3) = gf(2,3)
           u_f((i-1)*3+2,(ny_f-4)+3) = gf(3,3)
        END DO
     END DO
  END DO
  
  DO j=1,ny_c-2
     gc(1,1)=u_c(nx_c-1,j)
     gc(2,1)=u_c(1,j)
     gc(1,2)=u_c(nx_c-1,j+1)
     gc(2,2)=u_c(1,j+1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           u_f(nx_f-4  ,(j-1)*3+1) = gf(1,1)
           u_f(nx_f-4+1,(j-1)*3+1) = gf(2,1)
           u_f(nx_f-4+2,(j-1)*3+1) = gf(3,1)
           u_f(nx_f-4  ,(j-1)*3+2) = gf(1,2)
           u_f(nx_f-4+1,(j-1)*3+2) = gf(2,2)
           u_f(nx_f-4+2,(j-1)*3+2) = gf(3,2)
           u_f(nx_f-4  ,(j-1)*3+3) = gf(1,3)
           u_f(nx_f-4+1,(j-1)*3+3) = gf(2,3)
           u_f(nx_f-4+2,(j-1)*3+3) = gf(3,3)
        END DO
     END DO
  END DO
  
  gc(1,1)=u_c(nx_c-1,ny_c-1)
  gc(2,1)=u_c(1,ny_c-1)
  gc(1,2)=u_c(nx_c-1,1)
  gc(2,2)=u_c(1,1)
  CALL map_fine_cell(gf,gc) 
  DO k=1,3
     DO l=1,3
        u_f(nx_f-4  ,ny_f-4+1) = gf(1,1)
        u_f(nx_f-4+1,ny_f-4+1) = gf(2,1)
        u_f(nx_f-4+2,ny_f-4+1) = gf(3,1)
        u_f(nx_f-4  ,ny_f-4+2) = gf(1,2)
        u_f(nx_f-4+1,ny_f-4+2) = gf(2,2)
        u_f(nx_f-4+2,ny_f-4+2) = gf(3,2)
        u_f(nx_f-4  ,ny_f-4+3) = gf(1,3)
        u_f(nx_f-4+1,ny_f-4+3) = gf(2,3)
        u_f(nx_f-4+2,ny_f-4+3) = gf(3,3)
     END DO
  END DO
  
  DO j=1,ny_f-1
     u_f(nx_f-1,j) = u_f(0,j)
  END DO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Vertical velocity:
  DO j=1,ny_c-2
     DO i=1,nx_c-2
        gc(1,1)=v_c(i,j)
        gc(2,1)=v_c(i+1,j)
        gc(1,2)=v_c(i,j+1)
        gc(2,2)=v_c(i+1,j+1)
        CALL map_fine_cell(gf,gc) 
        DO k=1,3
           DO l=1,3
              v_f((i-1)*3+1,(j-1)*3  ) = gf(1,1)
              v_f((i-1)*3+2,(j-1)*3  ) = gf(2,1)
              v_f((i-1)*3+3,(j-1)*3  ) = gf(3,1)
              v_f((i-1)*3+1,(j-1)*3+1) = gf(1,2)
              v_f((i-1)*3+2,(j-1)*3+1) = gf(2,2)
              v_f((i-1)*3+3,(j-1)*3+1) = gf(3,2)
              v_f((i-1)*3+1,(j-1)*3+2) = gf(1,3)
              v_f((i-1)*3+2,(j-1)*3+2) = gf(2,3)
              v_f((i-1)*3+3,(j-1)*3+2) = gf(3,3)
           END DO
        END DO
     END DO
  END DO
  
  DO i=1,nx_c-2
     gc(1,1)=v_c(i,ny_c-1)
     gc(2,1)=v_c(i+1,ny_c-1)
     gc(1,2)=v_c(i,1)
     gc(2,2)=v_c(i+1,1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           v_f((i-1)*3+1,(ny_f-4)  ) = gf(1,1)
           v_f((i-1)*3+2,(ny_f-4)  ) = gf(2,1)
           v_f((i-1)*3+3,(ny_f-4)  ) = gf(3,1)
           v_f((i-1)*3+1,(ny_f-4)+1) = gf(1,2)
           v_f((i-1)*3+2,(ny_f-4)+1) = gf(2,2)
           v_f((i-1)*3+3,(ny_f-4)+1) = gf(3,2)
           v_f((i-1)*3+1,(ny_f-4)+2) = gf(1,3)
           v_f((i-1)*3+2,(ny_f-4)+2) = gf(2,3)
           v_f((i-1)*3+3,(ny_f-4)+2) = gf(3,3)
        END DO
     END DO
  END DO
  
  DO j=1,ny_c-2
     gc(1,1)=v_c(nx_c-1,j)
     gc(2,1)=v_c(1,j)
     gc(1,2)=v_c(nx_c-1,j+1)
     gc(2,2)=v_c(1,j+1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           v_f(nx_f-4+1,(j-1)*3  ) = gf(1,1)
           v_f(nx_f-4+2,(j-1)*3  ) = gf(2,1)
           v_f(nx_f-4+3,(j-1)*3  ) = gf(3,1)
           v_f(nx_f-4+1,(j-1)*3+1) = gf(1,2)
           v_f(nx_f-4+2,(j-1)*3+1) = gf(2,2)
           v_f(nx_f-4+3,(j-1)*3+1) = gf(3,2)
           v_f(nx_f-4+1,(j-1)*3+2) = gf(1,3)
           v_f(nx_f-4+2,(j-1)*3+2) = gf(2,3)
           v_f(nx_f-4+3,(j-1)*3+2) = gf(3,3)
        END DO
     END DO
  END DO
  
  gc(1,1)=v_c(nx_c-1,ny_c-1)
  gc(2,1)=v_c(1,ny_c-1)
  gc(1,2)=v_c(nx_c-1,1)
  gc(2,2)=v_c(1,1)
  CALL map_fine_cell(gf,gc) 
  DO k=1,3
     DO l=1,3
        v_f(nx_f-4+1,ny_f-4  ) = gf(1,1)
        v_f(nx_f-4+2,ny_f-4  ) = gf(2,1)
        v_f(nx_f-4+3,ny_f-4  ) = gf(3,1)
        v_f(nx_f-4+1,ny_f-4+1) = gf(1,2)
        v_f(nx_f-4+2,ny_f-4+1) = gf(2,2)
        v_f(nx_f-4+3,ny_f-4+1) = gf(3,2)
        v_f(nx_f-4+1,ny_f-4+2) = gf(1,3)
        v_f(nx_f-4+2,ny_f-4+2) = gf(2,3)
        v_f(nx_f-4+3,ny_f-4+2) = gf(3,3)
     END DO
  END DO
  
  DO i=1,nx_f-1
     v_f(i,ny_f-1) = v_f(i,0)
  END DO
  
  return
end subroutine coarsetofine



subroutine map_fine_cell(gf,gc)

  !Map individual cell
  USE rp_emulator
  implicit none

  REAL*8 :: gf(3,3),gc(2,2),tc,uc
  REAL*8 :: t1,t2,t3
  t1 = 1.0_8
  t2 = 2.0_8
  t3 = 3.0_8

!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !This is the grid:
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------!
!  gc(1,2)                                gc(2,2)
!
!  gf(1,3)           gf(2,3)   gf(3,3)  
!
!  gf(1,2)           gf(2,2)   gf(3,2)  
!
!  gc(1,1)/gf(1,1)   gf(2,1)   gf(3,1)    gc(2,1)
!--------------------------------------------------!



  !See numerical recipes 3.6 p. 123
  gf(1,1)=gc(1,1)
  gf(1,2)=t2/t3*gc(1,1)+t1/t3*gc(1,2)
  gf(1,3)=t1/t3*gc(1,1)+t2/t3*gc(1,2)
  gf(2,1)=t2/t3*gc(1,1)+t1/t3*gc(2,1)
  gf(3,1)=t1/t3*gc(1,1)+t2/t3*gc(2,1)
  tc=t1/t3
  uc=t1/t3
  gf(2,2)=(t1-tc)*(t1-uc)*gc(1,1)+tc*(t1-uc)*gc(2,1)+tc*uc*gc(2,2)+(t1-tc)*uc*gc(1,2)
  tc=t1/t3
  uc=t2/t3
  gf(2,3)=(t1-tc)*(t1-uc)*gc(1,1)+tc*(t1-uc)*gc(2,1)+tc*uc*gc(2,2)+(t1-tc)*uc*gc(1,2)
  tc=t2/t3
  uc=t1/t3
  gf(3,2)=(t1-tc)*(t1-uc)*gc(1,1)+tc*(t1-uc)*gc(2,1)+tc*uc*gc(2,2)+(t1-tc)*uc*gc(1,2)
  tc=t2/t3
  uc=t2/t3
  gf(3,3)=(t1-tc)*(t1-uc)*gc(1,1)+tc*(t1-uc)*gc(2,1)+tc*uc*gc(2,2)+(t1-tc)*uc*gc(1,2)

  return
END subroutine map_fine_cell



subroutine coarsetofine_h(nx_c,ny_c,h_c,&
          &nx_f,ny_f,h_f)

  USE rp_emulator
  implicit none

  INTEGER :: nx_c, ny_c
  REAL*8 :: h_c(0:nx_c,0:ny_c)

  INTEGER :: nx_f, ny_f
  REAL*8 :: h_f(0:nx_f,0:ny_f)
  INTEGER :: i,j,k,l

  REAL*8 :: gf(3,3),gc(2,2),tc,uc

  !Height:
  
  DO j=1,ny_c-2
     DO i=1,nx_c-2
        gc(1,1)=h_c(i,j)
        gc(2,1)=h_c(i+1,j)
        gc(1,2)=h_c(i,j+1)
        gc(2,2)=h_c(i+1,j+1)
        CALL map_fine_cell(gf,gc) 
        DO k=1,3
           DO l=1,3
              h_f((i-1)*3+1,(j-1)*3+1) = gf(1,1)
              h_f((i-1)*3+2,(j-1)*3+1) = gf(2,1)
              h_f((i-1)*3+3,(j-1)*3+1) = gf(3,1)
              h_f((i-1)*3+1,(j-1)*3+2) = gf(1,2)
              h_f((i-1)*3+2,(j-1)*3+2) = gf(2,2)
              h_f((i-1)*3+3,(j-1)*3+2) = gf(3,2)
              h_f((i-1)*3+1,(j-1)*3+3) = gf(1,3)
              h_f((i-1)*3+2,(j-1)*3+3) = gf(2,3)
              h_f((i-1)*3+3,(j-1)*3+3) = gf(3,3)
           END DO
        END DO
     END DO
  END DO

  DO i=1,nx_c-2
     gc(1,1)=h_c(i,ny_c-1)
     gc(2,1)=h_c(i+1,ny_c-1)
     gc(1,2)=h_c(i,1)
     gc(2,2)=h_c(i+1,1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           h_f((i-1)*3+1,(ny_f-4)+1) = gf(1,1)
           h_f((i-1)*3+2,(ny_f-4)+1) = gf(2,1)
           h_f((i-1)*3+3,(ny_f-4)+1) = gf(3,1)
           h_f((i-1)*3+1,(ny_f-4)+2) = gf(1,2)
           h_f((i-1)*3+2,(ny_f-4)+2) = gf(2,2)
           h_f((i-1)*3+3,(ny_f-4)+2) = gf(3,2)
           h_f((i-1)*3+1,(ny_f-4)+3) = gf(1,3)
           h_f((i-1)*3+2,(ny_f-4)+3) = gf(2,3)
           h_f((i-1)*3+3,(ny_f-4)+3) = gf(3,3)
        END DO
     END DO
  END DO
  
  DO j=1,ny_c-2
     gc(1,1)=h_c(nx_c-1,j)
     gc(2,1)=h_c(1,j)
     gc(1,2)=h_c(nx_c-1,j+1)
     gc(2,2)=h_c(1,j+1)
     CALL map_fine_cell(gf,gc) 
     DO k=1,3
        DO l=1,3
           h_f(nx_f-4+1,(j-1)*3+1) = gf(1,1)
           h_f(nx_f-4+2,(j-1)*3+1) = gf(2,1)
           h_f(nx_f-4+3,(j-1)*3+1) = gf(3,1)
           h_f(nx_f-4+1,(j-1)*3+2) = gf(1,2)
           h_f(nx_f-4+2,(j-1)*3+2) = gf(2,2)
           h_f(nx_f-4+3,(j-1)*3+2) = gf(3,2)
           h_f(nx_f-4+1,(j-1)*3+3) = gf(1,3)
           h_f(nx_f-4+2,(j-1)*3+3) = gf(2,3)
           h_f(nx_f-4+3,(j-1)*3+3) = gf(3,3)
        END DO
     END DO
  END DO
  
  gc(1,1)=h_c(nx_c-1,ny_c-1)
  gc(2,1)=h_c(1,ny_c-1)
  gc(1,2)=h_c(nx_c-1,1)
  gc(2,2)=h_c(1,1)
  CALL map_fine_cell(gf,gc) 
  DO k=1,3
     DO l=1,3
        h_f(nx_f-4+1,ny_f-4+1) = gf(1,1)
        h_f(nx_f-4+2,ny_f-4+1) = gf(2,1)
        h_f(nx_f-4+3,ny_f-4+1) = gf(3,1)
        h_f(nx_f-4+1,ny_f-4+2) = gf(1,2)
        h_f(nx_f-4+2,ny_f-4+2) = gf(2,2)
        h_f(nx_f-4+3,ny_f-4+2) = gf(3,2)
        h_f(nx_f-4+1,ny_f-4+3) = gf(1,3)
        h_f(nx_f-4+2,ny_f-4+3) = gf(2,3)
        h_f(nx_f-4+3,ny_f-4+3) = gf(3,3)
     END DO
  END DO
  
END subroutine COARSETOFINE_H
   

subroutine initialise(itestcase,lboundary,lrestart,Lx,Ly,gp,cinput,nx_i,ny_i,ninter,&
     & nx_f,ny_f,nt_f,nuh_f,nuu_f,dt_f,HC_f,&
     & dx_f,dy_f,rdx_f,rdy_f,ab_f,fu_f,fv_f,taux_f,tauy_f,h_f,dh_f,u_f,du_f,v_f,dv_f, &
     & nx_c,ny_c,nt_c,nuh_c,nuu_c,dt_c,HC_c,&
     & dx_c,dy_c,rdx_c,rdy_c,ab_c,fu_c,fv_c,taux_c,tauy_c,h_c,dh_c,u_c,du_c,v_c,dv_c)

!Initialise testcases

  USE rp_emulator
  implicit none

  LOGICAL :: lboundary 
  LOGICAL :: lrestart
  INTEGER :: itestcase
  REAL*8 :: Lx,Ly,h0,f0,beta,gp
  character*10 :: cinput  
  INTEGER :: nx_i,ny_i,ninter

  INTEGER :: nx_f,ny_f,nt_f
  REAL*8 :: nuh_f,nuu_f,dt_f,dx_f,dy_f,rdx_f,rdy_f,ab_f(nt_f),fu_f(0:ny_f),fv_f(0:ny_f),taux_f(0:ny_f),tauy_f(0:nx_f)
  REAL*8 :: h_f(0:nx_f,0:ny_f),dh_f(0:nx_f,0:ny_f,nt_f),u_f(0:nx_f,0:ny_f),du_f(0:nx_f,0:ny_f,nt_f),&
       & v_f(0:nx_f,0:ny_f),dv_f(0:nx_f,0:ny_f,nt_f)
  REAL*8 :: HC_f(0:nx_f,0:ny_f)

  INTEGER :: nx_c,ny_c,nt_c
  REAL*8 :: nuh_c,nuu_c,dt_c,dx_c,dy_c,rdx_c,rdy_c,ab_c(nt_c),fu_c(0:ny_c),fv_c(0:ny_c),taux_c(0:ny_c),tauy_c(0:nx_c)
  REAL*8 :: h_c(0:nx_c,0:ny_c),dh_c(0:nx_c,0:ny_c,nt_c),u_c(0:nx_c,0:ny_c),du_c(0:nx_c,0:ny_c,nt_c),&
       & v_c(0:nx_c,0:ny_c),dv_c(0:nx_c,0:ny_c,nt_c)
  REAL*8 :: HC_c(0:nx_c,0:ny_c)

  REAL*8 :: pi
  REAL*8 :: rdummy1, rdummy2,sigmax,sigmay
  INTEGER :: i,j
  REAL*8 :: h_i(nx_i,ny_i),u_i(nx_i,ny_i),v_i(nx_i,ny_i)

  real*8 :: vec(5)
  real*8 :: real_dummy

  IF(itestcase==1.or.itestcase==3.or.itestcase==4.or.itestcase==5)THEN
     lboundary = .FALSE.
  ELSE
     lboundary = .TRUE.
  END IF

  pi=3.14159265358979_8
  gp=9.81_8
  nuh_f=0._8
  nuh_c=0._8

  IF(itestcase==1.or.itestcase==3)THEN
     Lx=5000000.0_8
     Ly=5000000.0_8
     nuu_c=0.0_8
     h0=5000.0_8
     dt_c=10.0_8
     f0=4.46e-5_8
     beta=0.0_8
  ELSE IF(itestcase==2)THEN
     Lx=3480000.0_8
     Ly=3480000.0_8
     nuu_c=470.23_8   
     h0=500.0_8
     dt_c=25.0_8
     f0=4.46e-5_8
     beta=2.e-11_8
  ELSE IF(itestcase==4)THEN
     Lx=10000000.0_8
     Ly=10000000.0_8
     nuu_c=0.0_8
     h0=400.0_8
     dt_c=25.0_8
     f0=1.0e-4_8
     beta=0.0_8
  ELSE IF(itestcase==5)THEN
     Lx= 1200000.0_8
     Ly=  200000.0_8
     nuu_c=0.0_8
     h0=400.0_8
     dt_c=2.0_8
     f0=0.0_8
     beta=0.0_8 
  ELSE
     write(*,*) 'Wrong test case!'
     STOP
  END IF
  dx_f=Lx/real(nx_f-1,8)
  dy_f=Ly/real(ny_f-1,8)
  dx_c=dx_f*3.0_8
  dy_c=dy_f*3.0_8
  rdx_f=1.0_8/dx_f
  rdy_f=1.0_8/dy_f
  rdx_c=1.0_8/dx_c
  rdy_c=1.0_8/dy_c
  dt_f = dt_c
  nuu_f = nuu_c

! Adams Bashforth parameters
  ab_c(1)=(23._8/12._8)*dt_c
  ab_c(2)=-(16._8/12._8)*dt_c
  ab_c(3)=(5._8/12._8)*dt_c
  ab_f(1)=(23._8/12._8)*dt_f
  ab_f(2)=-(16._8/12._8)*dt_f
  ab_f(3)=(5._8/12._8)*dt_f

! Define Coriolis parameter and initial u,v,h
  do j=0,ny_c
     fu_c(j)=f0+beta*Ly*(real(j,8)-0.5_8)/real(ny_c-1,8)
     fv_c(j)=f0+beta*Ly*real(j-1,8)/real(ny_c-1,8)
  end do
  do j=0,ny_f
     fu_f(j)=f0+beta*Ly*(real(j,8)-0.5_8)/real(ny_f-1,8)
     fv_f(j)=f0+beta*Ly*real(j-1,8)/real(ny_f-1,8)
  end do


! Define the wind forcing:
  IF(itestcase==1.or.itestcase==3.or.itestcase==4.or.itestcase==5)THEN
     DO i=0,ny_c-1        
        taux_c(i) = 0.0_8
     END DO
     DO i=0,ny_f-1        
        taux_f(i) = 0.0_8
     END DO
  ELSE IF(itestcase==2)THEN
     DO i=0,ny_c-1        
        taux_c(i) = 0.12_8*(cos(2.0_8*pi*((real(i,8)-0.5_8)*Ly/real(ny_c-1,8)-0.5_8*Ly)/Ly)&
             & +2.0_8*sin(pi*((real(i,8)-0.5_8)*Ly/real(ny_c-1,8)-0.5_8*Ly)/Ly))/(999.8_8*h0)
        taux_f(i) = 0.12_8*(cos(2.0_8*pi*((real(i,8)-0.5_8)*Ly/real(ny_f-1,8)-0.5_8*Ly)/Ly)&
             & +2.0_8*sin(pi*((real(i,8)-0.5_8)*Ly/real(ny_f-1,8)-0.5_8*Ly)/Ly))/(999.8_8*h0)
     END DO
  END IF

  IF(itestcase.le.4)THEN
     DO i=0,nx_c
        DO j=0,ny_c
           HC_c(i,j) = h0
        END DO
     END DO
     DO i=0,nx_f
        DO j=0,ny_f
           HC_f(i,j) = h0
        END DO
     END DO
  ELSE IF(itestcase==5)THEN
     sigmax = 3.0_8*Ly/20.0_8
     sigmay = 3.0_8*Ly/20.0_8
     DO i=1,nx_c-1
        DO j=1,ny_c-1
           rdummy1 = (Lx*real(i-1,8)/real(nx_c-1,8)-Lx/8.0_8)*(Lx*real(i-1,8)/real(nx_c-1,8)-Lx/8.0_8)/(sigmax*sigmax)
           rdummy2 = (Ly*real(j-1,8)/real(ny_c-1,8)-Ly/2.0_8)*(Ly*real(j-1,8)/real(ny_c-1,8)-Ly/2.0_8)/(sigmay*sigmay)
           HC_c(i,j)=h0-100.0_8*exp(-rdummy1-rdummy2)
        END DO
     END DO
     DO i=1,nx_f-1
        DO j=1,ny_f-1
           rdummy1 = (Lx*real(i-1,8)/real(nx_f-1,8)-Lx/8.0_8)*(Lx*real(i-1,8)/real(nx_f-1,8)-Lx/8.0_8)/(sigmax*sigmax)
           rdummy2 = (Ly*real(j-1,8)/real(ny_f-1,8)-Ly/2.0_8)*(Ly*real(j-1,8)/real(ny_f-1,8)-Ly/2.0_8)/(sigmay*sigmay)
           HC_f(i,j)=h0-100.0_8*exp(-rdummy1-rdummy2)
        END DO
     END DO
     do j=1,ny_c-1
        HC_c(0,j)=HC_c(nx_c-1,j)
        HC_c(nx_c,j)=HC_c(1,j)
     end do
     do j=1,nx_c-1
        HC_c(j,0)=HC_c(j,ny_c-1)
        HC_c(j,ny_c)=HC_c(j,1)
     end do
     HC_c(0,0)=HC_c(nx_c-1,ny_c-1)
     HC_c(nx_c,ny_c)=HC_c(1,1)
     HC_c(0,ny_c)=HC_c(nx_c-1,1)
     HC_c(nx_c,0)=HC_c(1,ny_c-1)

    do j=1,ny_f-1
        HC_f(0,j)=HC_f(nx_f-1,j)
        HC_f(nx_f,j)=HC_f(1,j)
     end do
     do j=1,nx_f-1
        HC_f(j,0)=HC_f(j,ny_f-1)
        HC_f(j,ny_f)=HC_f(j,1)
     end do
     HC_f(0,0)=HC_f(nx_f-1,ny_f-1)
     HC_f(nx_f,ny_f)=HC_f(1,1)
     HC_f(0,ny_f)=HC_f(nx_f-1,1)
     HC_f(nx_f,0)=HC_f(1,ny_f-1)
  END IF

  DO i=0,nx_c
     tauy_c(i)= 0.0_8
  END DO
  DO i=0,nx_f
     tauy_f(i)= 0.0_8
  END DO
  
  !Define fields:
  IF(.not.lrestart)THEN
     IF(itestcase==1.or.itestcase==3)THEN
        sigmax = 3.0_8*Lx/20.0_8
        sigmay = 3.0_8*Ly/20.0_8
        DO j=1,ny_c-1
           DO i=1,nx_c-1
              rdummy1 = (Lx*real(i-1,8)/real(nx_c-1,8)-Lx/2.0_8)*(Lx*real(i-1,8)/real(nx_c-1,8)-Lx/2.0_8)/(sigmax*sigmax)
              rdummy2 = (Ly*real(j-1,8)/real(ny_c-1,8)-Ly/2.0_8)*(Ly*real(j-1,8)/real(ny_c-1,8)-Ly/2.0_8)/(sigmay*sigmay)

              h_c(i,j)=100.0_8*exp(-rdummy1-rdummy2)

              IF(itestcase==3)THEN
                 rdummy1 = (Lx*real(i-1.5_8,8)/real(nx_c-1,8)-Lx/2.0_8)*(Lx*real(i-1.5_8,8)&
                      &/real(nx_c-1,8)-Lx/2.0_8)/(sigmax*sigmax)
                 rdummy2 = (Ly*real(j-1,8)/real(ny_c-1,8)-Ly/2.0_8)*(Ly*real(j-1,8)&
                      &/real(ny_c-1,8)-Ly/2.0_8)/(sigmay*sigmay)        
                 u_c(i,j)= +(100.0_8*gp*2.0_8*(Ly*real(j-1,8)/real(ny_c-1,8)-Ly/2.0_8))&
                      &/(6.147E-5_8*3.0_8*Ly/20.0_8*3.0_8*Ly/20.0_8)*exp(-rdummy1-rdummy2)
                 
                 rdummy1 = (Lx*real(i-1,8)/real(nx_c-1,8)-Lx/2.0_8)*(Lx*real(i-1,8)&
                      &/real(nx_c-1,8)-Lx/2.0_8)/(sigmax*sigmax)
                 rdummy2 = (Ly*(real(j,8)-1.5_8)/real(ny_c-1,8)-Ly/2.0_8)*(Ly*(real(j,8)-1.5_8)&
                      &/real(ny_c-1,8)-Ly/2.0_8)/(sigmay*sigmay)              
                 v_c(i,j)= -(100.0_8*gp*2.0_8*(Lx*real(i-1,8)/real(nx_c-1,8)-Lx/2.0_8))&
                      &/(6.147E-5*3.0_8*Lx/20.0_8*3.0_8*Lx/20.0_8)*exp(-rdummy1-rdummy2)
              ELSE
                 u_c=0._8 
                 v_c=0._8 
              END IF
           END DO
        END DO
        DO j=1,ny_f-1
           DO i=1,nx_f-1
              rdummy1 = (Lx*real(i-1,8)/real(nx_f-1,8)-Lx/2.0_8)*(Lx*real(i-1,8)/real(nx_f-1,8)-Lx/2.0_8)/(sigmax*sigmax)
              rdummy2 = (Ly*real(j-1,8)/real(ny_f-1,8)-Ly/2.0_8)*(Ly*real(j-1,8)/real(ny_f-1,8)-Ly/2.0_8)/(sigmay*sigmay)
              h_f(i,j)=100.0_8*exp(-rdummy1-rdummy2)

              IF(itestcase==3)THEN
                 rdummy1 = (Lx*real(i-1.5_8,8)/real(nx_f-1,8)-Lx/2.0_8)*(Lx*real(i-1.5_8,8)&
                      &/real(nx_f-1,8)-Lx/2.0_8)/(sigmax*sigmax)
                 rdummy2 = (Ly*real(j-1,8)/real(ny_f-1,8)-Ly/2.0_8)*(Ly*real(j-1,8)&
                      &/real(ny_f-1,8)-Ly/2.0_8)/(sigmay*sigmay)        
                 u_f(i,j)= +(100.0_8*gp*2.0_8*(Ly*real(j-1,8)/real(ny_f-1,8)-Ly/2.0_8))&
                      &/(6.147E-5_8*3.0_8*Ly/20.0_8*3.0_8*Ly/20.0_8)*exp(-rdummy1-rdummy2)
                 
                 rdummy1 = (Lx*real(i-1,8)/real(nx_f-1,8)-Lx/2.0_8)*(Lx*real(i-1,8)&
                      &/real(nx_f-1,8)-Lx/2.0_8)/(sigmax*sigmax)
                 rdummy2 = (Ly*(real(j,8)-1.5_8)/real(ny_f-1,8)-Ly/2.0_8)*(Ly*(real(j,8)-1.5_8)&
                      &/real(ny_f-1,8)-Ly/2.0_8)/(sigmay*sigmay)              
                 v_f(i,j)= -(100.0_8*gp*2.0_8*(Lx*real(i-1,8)/real(nx_f-1,8)-Lx/2.0_8))&
                      &/(6.147E-5*3.0_8*Lx/20.0_8*3.0_8*Lx/20.0_8)*exp(-rdummy1-rdummy2)
              ELSE
                 u_f=0._8 
                 v_f=0._8 
              END IF
           END DO
        END DO
     ELSE IF(itestcase==2)THEN
        h_c= 0.0_8
        h_f= 0.0_8
        u_c=0._8
        u_f=0._8           
        v_c=0._8
        v_f=0._8           
     ELSE IF(itestcase==4)THEN
        DO j=1,ny_f-1
           DO i=1,nx_f-1
              CALL random_number(real_dummy)
              h_f(i,j)= 100.0_8*(real_dummy-0.5_8)
              CALL random_number(real_dummy)
              u_f(i,j)= real_dummy-0.5_8
              CALL random_number(real_dummy)
              v_f(i,j)= real_dummy-0.5_8
           END DO
        END DO
     ELSE IF(itestcase==5)THEN
        h_c= 0.0_8
        h_f= 0.0_8
        u_c=10.0_8
        u_f=10.0_8           
        v_c=0._8
        v_f=0._8   
     END IF
  ELSE 
     open(12,file='./../improve_fine/long_ref/Output/ufine.'//TRIM(cinput), STATUS='OLD', ACTION='read')
     open(13,file='./../improve_fine/long_ref/Output/vfine.'//TRIM(cinput), STATUS='OLD', ACTION='read')
     open(14,file='./../improve_fine/long_ref/Output/hfine.'//TRIM(cinput), STATUS='OLD', ACTION='read')
     
     IF(nx_f.gt.nx_i.or.ny_f.gt.ny_i)THEN
        write(*,*) 'Model input needs to be on a larger grid!'
        STOP
     END IF
     
     IF(nx_f==nx_i.and.ny_f==ny_i)THEN
        do j=1,ny_i-1
           do i=1,nx_i-1
              read(14,*) vec
              h_f(i,j) = vec(3)
              dh_f(i,j,1:(nt_f-1)) = vec(4:)
           end do
        end do
        
        IF(lboundary)THEN
           do j=1,ny_i-1
              do i=1,nx_i
                 read(12,*) vec
                 u_f(i,j) = vec(3)
                 du_f(i,j,1:(nt_f-1)) = vec(4:)
              end do
           end do
           do j=1,ny_i
              do i=1,nx_i-1
                 read(13,*) vec
                 v_f(i,j) = vec(3)
                 dv_f(i,j,1:(nt_f-1)) = vec(4:)
              end do
           end do
        ELSE
           do j=1,ny_i-1
              do i=1,nx_i-1
                 read(12,*) vec
                 u_f(i,j) = vec(3)
                 du_f(i,j,1:(nt_f-1)) = vec(4:)
              end do
           end do
           do j=1,ny_i-1
              do i=1,nx_i-1
                 read(13,*) vec
                 v_f(i,j) = vec(3)
                 dv_f(i,j,1:(nt_f-1)) = vec(4:)
              end do
           end do
        END IF
        CLOSE(12)
        CLOSE(13)
        CLOSE(14)
        CALL finetocoarse(nx_c,ny_c,h_c,u_c,v_c,&
          &nx_f,ny_f,h_f,u_f,v_f)
     ELSE IF((nx_f-1)==(nx_i-1)/3.and.(ny_f-1)==(ny_i-1)/3)THEN
        do j=1,ny_i-1
           do i=1,nx_i-1
              read(12,*) vec
              u_i(i,j) = vec(3)
              read(13,*) vec
              v_i(i,j) = vec(3)
              read(14,*) vec
              h_i(i,j) = vec(3)
           end do
        end do
        CALL finetocoarse(nx_f,ny_f,h_f,u_f,v_f,&
          &nx_i,ny_i,h_i,u_i,v_i)
        CALL finetocoarse(nx_c,ny_c,h_c,u_c,v_c,&
          &nx_f,ny_f,h_f,u_f,v_f)
        IF(lboundary)THEN
           write(*,*) 'Boundary restart not implemented yet...'
           STOP
        END IF
     END IF
  END IF
  dh_c=0._8
  dh_f=0._8
  du_c=0._8
  du_f=0._8
  dv_c=0._8
  dv_f=0._8
  
     !Set condition for periodicity:
     IF(.not.lboundary)THEN
        do j=1,ny_f-1
           u_f(0,j) = u_f(nx_f-1,j)
           u_f(nx_f,j) = u_f(1,j)
           v_f(0,j) = v_f(nx_f-1,j)
           v_f(nx_f,j) = v_f(1,j)
           h_f(0,j)=h_f(nx_f-1,j)
           h_f(nx_f,j)=h_f(1,j)
        end do
        do j=1,nx_f-1
           u_f(j,0) = u_f(j,ny_f-1)
           u_f(j,ny_f) = u_f(j,1)
           v_f(j,0) = v_f(j,ny_f-1)
           v_f(j,ny_f) = v_f(j,1)
           h_f(j,0)=h_f(j,ny_f-1)
           h_f(j,ny_f)=h_f(j,1)
        end do
        u_f(0,0)=u_f(nx_f-1,ny_f-1)
        v_f(0,0)=v_f(nx_f-1,ny_f-1)
        h_f(0,0)=h_f(nx_f-1,ny_f-1)
        u_f(nx_f,ny_f)=u_f(1,1)
        v_f(nx_f,ny_f)=v_f(1,1)
        h_f(nx_f,ny_f)=h_f(1,1)
        u_f(0,ny_f)=u_f(nx_f-1,1)
        v_f(0,ny_f)=v_f(nx_f-1,1)
        h_f(0,ny_f)=h_f(nx_f-1,1)
        u_f(nx_f,0)=u_f(1,ny_f-1)
        v_f(nx_f,0)=v_f(1,ny_f-1)
        h_f(nx_f,0)=h_f(1,ny_f-1)
     END IF
  return
end subroutine initialise


subroutine output_fields(lboundary,nsec,nxstep,nystep,&
     & u_c,du_c,v_c,dv_c,h_c,dh_c,dx_c,dy_c,nx_c,ny_c,nt_c,&
     & u_f,du_f,v_f,dv_f,h_f,dh_f,dx_f,dy_f,nx_f,ny_f,nt_f)

!Write model output:

  USE rp_emulator
  implicit none

  LOGICAL :: lboundary
  INTEGER :: nsec,nxstep,nystep

  INTEGER :: nx_c, ny_c, nt_c
  REAL*8 :: u_c(0:nx_c,0:ny_c),v_c(0:nx_c,0:ny_c),h_c(0:nx_c,0:ny_c),dx_c,dy_c
  REAL*8 :: dh_c(0:nx_c,0:ny_c,nt_c),du_c(0:nx_c,0:ny_c,nt_c),dv_c(0:nx_c,0:ny_c,nt_c)
  
  INTEGER :: nx_f, ny_f, nt_f
  REAL*8 :: u_f(0:nx_f,0:ny_f),v_f(0:nx_f,0:ny_f),h_f(0:nx_f,0:ny_f),dx_f,dy_f
  REAL*8 :: dh_f(0:nx_f,0:ny_f,nt_f),du_f(0:nx_f,0:ny_f,nt_f),dv_f(0:nx_f,0:ny_f,nt_f)

  character*9 :: num

  !Write to screen:
  print *,' '
  print *,'Time elapsed: ',nsec,' seconds'
  print *,' '
  call contour(h_f,nx_f,ny_f,nxstep,nystep,'h')

  !Write to file: 

  IF (nsec.le.9) write(num,'(I1)') nsec
  IF (nsec.ge.10        .and.nsec.le.99) write(num,'(I2)') nsec
  IF (nsec.ge.100       .and.nsec.le.999) write(num,'(I3)') nsec
  IF (nsec.ge.1000      .and.nsec.le.9999) write(num,'(I4)') nsec
  IF (nsec.ge.10000     .and.nsec.le.99999) write(num,'(I5)') nsec
  IF (nsec.ge.100000    .and.nsec.le.999999) write(num,'(I6)') nsec
  IF (nsec.ge.1000000   .and.nsec.le.9999999) write(num,'(I7)') nsec
  IF (nsec.ge.10000000  .and.nsec.le.99999999) write(num,'(I8)') nsec
  IF (nsec.ge.100000000 .and.nsec.le.999999999) write(num,'(I9)') nsec
  IF (nsec.ge.1000000000)THEN 
     write(*,*) 'Output number too big!!!'
     STOP
  END IF

  !Write coarse field:
  call write_hdata_file(h_c,dh_c,nt_c,nx_c,ny_c,dx_c,dy_c,'./Output/hcoarse.'//num)
  call write_udata_file_period(u_c,du_c,nt_c,nx_c,ny_c,dx_c,dy_c,'./Output/ucoarse.'//num)
  call write_vdata_file_period(v_c,dv_c,nt_c,nx_c,ny_c,dx_c,dy_c,'./Output/vcoarse.'//num)           
  
  !Write fine field:
  call write_hdata_file(h_f,dh_f,nt_f,nx_f,ny_f,dx_f,dy_f,'./Output/hfine.'//num)
  call write_udata_file_period(u_f,du_f,nt_f,nx_f,ny_f,dx_f,dy_f,'./Output/ufine.'//num)
  call write_vdata_file_period(v_f,dv_f,nt_f,nx_f,ny_f,dx_f,dy_f,'./Output/vfine.'//num)           

END SUBROUTINE output_fields


subroutine contour(c,nx,ny,nxstep,nystep,name)

! Produce simple contour plot on screen
  USE rp_emulator
  implicit none

  REAL*8 :: cmin,cmax
  integer nx,ny,nxstep,nystep
  integer i,j
  integer nc(nx-1,ny-1)
  REAL*8 :: c(0:nx,0:ny)
  REAL*8 :: dummy
  character(len=*) name
  real*8 :: rcmin,rcmax

!    determine maximum and minimum values and write to screen

  cmax=maxval(c(:,:))
  cmin=minval(c(:,:))
  rcmin =cmin
  rcmax= cmax
  print *,name,'   max:',rcmax,'; min:',rcmin
!
!    convert values to numbers between 0 and 9

  do j=1,ny-1
     do i=1,nx-1
        dummy = 9.999_8*(c(i,j)-cmin)/(cmax-cmin+1.e-9_8)
        nc(i,j)=int(dummy)
     end do
  end do

!    write values to screen to produce contour plot

  do j=ny-1,1,-nystep
     write(6,100)(nc(i,j),i=1,nx-1,nxstep)
  end do
100 format(1x,400i1)
  
  return
end subroutine contour

!    --------------------------------------------------------------------------

subroutine write_hdata_file(data,ddata,nt_f,nx,ny,dx,dy,filename)

!    function: write data array to file

  USE rp_emulator
  implicit none
  
  integer nx,ny,nt_f
  REAL*8 :: data(0:nx,0:ny),ddata(0:nx,0:ny,nt_f),dx,dy
  character(len=*) filename
  integer i,j
  real*8 :: rvec(nt_f+2)

  print *,filename
  open(unit=9,file=filename,status='unknown')      
  
  do j=1,ny-1
     do i=1,nx-1
        rvec(1) = real(i-1,8)*dx
        rvec(2) = real(j-1,8)*dy
        rvec(3) = data(i,j)
        rvec(4:) = ddata(i,j,1:(nt_f-1))
        write(9,*) rvec
     end do
  end do
  
  close(9)

  return
end subroutine write_hdata_file

!    --------------------------------------------------------------------------

subroutine write_udata_file(data,ddata,nt_f,nx,ny,dx,dy,filename)

!    function: write data array to file

  USE rp_emulator
  implicit none

  integer nx,ny,nt_f
  REAL*8 :: data(0:nx,0:ny),ddata(0:nx,0:ny,nt_f),dx,dy
  character(len=*) filename
  integer i,j
  real*8 :: rvec(nt_f+2)

  print *,filename
  open(unit=9,file=filename,status='unknown')      

  do j=1,ny-1
     do i=1,nx
        rvec(1) = (real(i,8)-0.5_8)*dx
        rvec(2) = real(j-1,8)*dy
        rvec(3) = data(i,j)
        rvec(4:) = ddata(i,j,1:(nt_f-1))
        write(9,*) rvec
     end do
  end do

  close(9)

  return
end subroutine write_udata_file

!    --------------------------------------------------------------------------

subroutine write_udata_file_period(data,ddata,nt_f,nx,ny,dx,dy,filename)

!    function: write data array to file

  USE rp_emulator
  implicit none

  integer nx,ny,nt_f
  REAL*8 :: data(0:nx,0:ny),ddata(0:nx,0:ny,nt_f),dx,dy
  character(len=*) filename
  integer i,j
  real*8 :: rvec(nt_f+2)

  print *,filename
  open(unit=9,file=filename,status='unknown')      

  do j=1,ny-1
     do i=1,nx-1
        rvec(1) = (real(i,8)-1.5_8)*dx
        rvec(2) = real(j-1,8)*dy
        rvec(3) = data(i,j)
        rvec(4:) = ddata(i,j,1:(nt_f-1))
        write(9,*) rvec
     end do
  end do

  close(9)

  return
end subroutine write_udata_file_period

!    --------------------------------------------------------------------------

subroutine write_vdata_file(data,ddata,nt_f,nx,ny,dx,dy,filename)

!    function: write data array to file

  USE rp_emulator
  implicit none

  integer nx,ny,nt_f
  REAL*8 :: data(0:nx,0:ny),ddata(0:nx,0:ny,nt_f),dx,dy
  character(len=*) filename
  integer i,j
  real*8 :: rvec(nt_f+2)

  print *,filename
  open(unit=9,file=filename,status='unknown')      

  do j=1,ny
     do i=1,nx-1
        rvec(1) = real(i-1,8)*dx
        rvec(2) = (real(j,8)-0.5_8)*dy
        rvec(3) = data(i,j)
        rvec(4:) = ddata(i,j,1:(nt_f-1))
        write(9,*) rvec
     end do
  end do

  close(9)

  return
end subroutine write_vdata_file
!
!    --------------------------------------------------------------------------

subroutine write_vdata_file_period(data,ddata,nt_f,nx,ny,dx,dy,filename)

!    function: write data array to file

  USE rp_emulator
  implicit none

  integer nx,ny,nt_f
  REAL*8 :: data(0:nx,0:ny),ddata(0:nx,0:ny,nt_f),dx,dy
  character(len=*) filename
  integer i,j
  real*8 :: rvec(nt_f+2)

  print *,filename
  open(unit=9,file=filename,status='unknown')      

  do j=1,ny-1
     do i=1,nx-1
        rvec(1) = real(i-1,8)*dx
        rvec(2) = (real(j,8)-1.5_8)*dy
        rvec(3) = data(i,j)
        rvec(4:) = ddata(i,j,1:(nt_f-1))
        write(9,*) rvec
     end do
  end do

  close(9)

  return
end subroutine write_vdata_file_period
