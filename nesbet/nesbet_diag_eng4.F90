program diagonalization
implicit none
integer,parameter :: wp=8
integer,parameter :: maxchar=100
integer,parameter :: ifp=100
integer,parameter :: ilp=20
integer,parameter :: iout=6
character(maxchar) :: buff
integer :: nloop
real(wp),allocatable :: hmat(:,:)
real(wp),allocatable :: hmat_mu_mu(:,:)
real(wp),allocatable :: smat(:,:)
real(wp),allocatable :: smat_mu_mu(:,:)
real(wp),allocatable :: inp_mat(:,:)
real(wp),allocatable :: dot_s(:)
real(wp),allocatable :: ccurr(:)
real(wp),allocatable :: sigma(:)
real(wp),allocatable :: delta_ccurr(:)
real(wp),allocatable :: test_mat_h(:,:)
real(wp),allocatable :: test_mat_s(:,:)
real(wp) :: ecurr
real(wp) :: xsum
real(wp) :: delta_d
real(wp) :: delta_e
real(wp) :: Dsum
real(wp) :: xnorm
integer :: nci
integer :: i,j

  !write(*,*) 'ENTER nci'
  !read(*,*) nci

  if(allocated(test_mat_h)) deallocate(test_mat_h)

  open(unit=ifp,file='test_mat_h.inp',action='read')
    read(ifp,'(A)') buff
    read(ifp,*) nci
    allocate(test_mat_h(nci,nci))
    read(ifp,'(A)') buff
      do i=1,nci
        read(ifp,*) (test_mat_h(i,j),j=1,nci)
      end do
    close(ilp)

    write(*,*) 'test_mat_h =',test_mat_h
  
  if(allocated(test_mat_s)) deallocate(test_mat_s)

  open(unit=ifp,file='test_mat_s.inp',action='read')
    read(ifp,'(A)') buff
    read(ifp,*) nci
    allocate(test_mat_s(nci,nci))
    read(ifp,'(A)') buff
      do i=1,nci
        read(ifp,*) (test_mat_s(i,j),j=1,nci)
      end do
    close(ilp)

    write(*,*) 'test_mat_s =',test_mat_s
  
  if(allocated(hmat)) deallocate(hmat)
  if(allocated(hmat_mu_mu)) deallocate(hmat_mu_mu)
  if(allocated(smat)) deallocate(smat)
  if(allocated(smat_mu_mu)) deallocate(smat_mu_mu)
  if(allocated(inp_mat)) deallocate(inp_mat)
  if(allocated(dot_s)) deallocate(dot_s)
  if(allocated(ccurr)) deallocate(ccurr)
  if(allocated(sigma)) deallocate(sigma)
  if(allocated(delta_ccurr)) deallocate(delta_ccurr)
  allocate(hmat(nci,nci))
  allocate(hmat_mu_mu(nci,nci))
  allocate(smat(nci,nci))
  allocate(smat_mu_mu(nci,nci))
  allocate(inp_mat(nci,nci))
  allocate(dot_s(nci))
  allocate(ccurr(nci))
  allocate(sigma(nci))
  allocate(delta_ccurr(nci))
  
  !write(*,*) 'ENTER',nci,'x',nci,'symmetric matrix'
  !  do i=1,nci
  !  do j=1,nci
  !    read(*,*) hmat(i,j)
  !  end do
  !  end do
!  write(*,*) 'hmat =',hmat
!  write(*,*) 'hmat(1,1) =',hmat(1,1)
!  write(*,*) 'hmat(1,2) =',hmat(1,2)
!  write(*,*) 'hmat(2,1) =',hmat(2,1)
!  write(*,*) 'hmat(2,2) =',hmat(2,2)
 
  !write(*,*) 'ENTER smat'
  !  do i=1,nci
  !  do j=1,nci
  !    read(*,*) smat(i,j)
  !  end do 
  !  end do 

!  write(*,*) 'smat =',smat
!  write(*,*) 'smat(1,1) =',smat(1,1)
!  write(*,*) 'smat(1,2) =',smat(1,2)
!  write(*,*) 'smat(2,1) =',smat(2,1)
!  write(*,*) 'smat(2,2) =',smat(2,2)

  write(*,*) 'ENTER ccurr'
    do i=1,nci  
      read(*,*) ccurr(i)
    end do
!  write(*,*) 'ccurr(1) =', ccurr(1)
!  write(*,*) 'ccurr(2) =', ccurr(2)
!  do i=1,5 
!  call calc_eng(nci,ccurr,hmat,smat,Dsum,ecurr)
 
!  write(*,*)
!  write(*,*) '************************This is the initial ecurr ****************************'
!  write(*,*) 'ecurr =',ecurr
!  write(*,*) '************************This is the initial ecurr ****************************'
!  write(*,*)

  hmat = test_mat_h
  smat = test_mat_s

  write(*,*) 'hmat =',hmat
  write(*,*) 'smat =',smat
  
!  write(*,*) 'smat(1,1) =',smat(1,1)
!  write(*,*) 'smat(1,2) =',smat(1,2)
!  write(*,*) 'smat(2,1) =',smat(2,1)
!  write(*,*) 'smat(2,2) =',smat(2,2)

  xnorm = 1.0e0_wp/sqrt(dot_product(ccurr,ccurr))
  ccurr(:) = xnorm*ccurr(:)
  write(*,*) "enter nloop"
  read(*,*) nloop
  do i=1,nloop
    write(*,*) '=================='
    write(*,*) 'ecurr =',ecurr
    write(*,*) '=================='
    call nesbet_diag(nci,ccurr,hmat,smat,ecurr)
    xnorm = 1.0e0_wp/sqrt(dot_product(ccurr,ccurr))
    ccurr(:) = xnorm*ccurr(:)
    write(*,*) '******************************************'
    write(*,*) 'ccurr =',ccurr
    write(*,*) 'root =',ecurr
    write(*,*) '******************************************'
  end do

  if(allocated(hmat)) deallocate(hmat)
  if(allocated(hmat_mu_mu)) deallocate(hmat_mu_mu)
  if(allocated(smat)) deallocate(smat)
  if(allocated(smat_mu_mu)) deallocate(smat_mu_mu)
  if(allocated(inp_mat)) deallocate(inp_mat)
  if(allocated(dot_s)) deallocate(dot_s)
  if(allocated(ccurr)) deallocate(ccurr)
  if(allocated(sigma)) deallocate(sigma)
  if(allocated(delta_ccurr)) deallocate(delta_ccurr)
  if(allocated(test_mat_h)) deallocate(test_mat_h)
  if(allocated(test_mat_s)) deallocate(test_mat_s)

contains

!subroutine calc_eng(nci,ccurr,hmat,smat,Dsum,ecurr)
!implicit none
!integer,parameter :: wp=8
!integer,intent(in) :: nci
!real(wp),intent(inout) :: ccurr(nci)
!real(wp),intent(inout) :: hmat(nci,nci)
!real(wp),intent(inout) :: smat(nci,nci)
!real(wp),intent(out) :: Dsum
!real(wp),intent(out) :: ecurr
!real(wp) :: Nsum 
  
!  Nsum = 0.0e0_wp
!  Dsum = 0.0e0_wp
!  
!  do i=1,nci
!  do j=1,nci
!    Nsum = Nsum + ( hmat(i,j) * ccurr(i) * ccurr(j) )
!    Dsum = Dsum + ( smat(i,j) * ccurr(i) * ccurr(j) )
!  end do 
!  end do 
!  write(*,*) 'Dsum =',Dsum
!  ecurr = ( Nsum / Dsum )
!end subroutine calc_eng

subroutine nesbet_diag(nci,ccurr,hmat,smat,ecurr)
implicit none
integer,parameter :: wp=8
integer,intent(in) :: nci
real(wp),intent(inout) :: ccurr(nci)
real(wp),intent(in) :: hmat(nci,nci)
real(wp),intent(in) :: smat(nci,nci)
real(wp),intent(out) :: ecurr
real(wp) :: hmat_mu_mu
real(wp) :: smat_mu_mu 
real(wp) :: dot_s(nci)
real(wp) :: sigma(nci) 
real(wp) :: delta_ccurr(nci)
real(wp) :: Dsum
real(wp) :: Nsum 
real(wp) :: xsum
real(wp) :: ysum
real(wp) :: delta_d
real(wp) :: delta_e
integer :: mu
integer :: i,j

  
!  write(*,*) 'nci =',nci
  write(*,*)
  write(*,*)
  
  Nsum = 0.0e0_wp
  Dsum = 0.0e0_wp

  do i=1,nci
  do j=1,nci
    Nsum = Nsum + ( hmat(i,j) * ccurr(i) * ccurr(j) )
    Dsum = Dsum + ( smat(i,j) * ccurr(i) * ccurr(j) )
    write(*,*) '--dsum--',dsum
  end do 
  end do 
  ecurr = ( Nsum / Dsum )
  
!  write(*,*) 'Nsum =',Nsum
!  write(*,*) 'Dsum =',Dsum
!  write(*,*) 'eccur 1 =',ecurr
  
!  do 
    mu_loop: do mu=1,nci
      xsum = 0.0e0_wp
      ysum = 0.0e0_wp
      hmat_mu_mu = hmat(mu,mu)
      smat_mu_mu = smat(mu,mu)
      write(*,*) 'hmat_mu_mu =', hmat_mu_mu
      write(*,*) 'smat_mu_mu =', smat_mu_mu
      trace_loop: do i=1,nci
        xsum = xsum + ( (hmat(mu,i)*ccurr(i)) - (ecurr*smat(mu,i)*ccurr(i)) ) 
        ysum = ysum + ( smat(mu,i)*ccurr(i) )
        write(*,*) 'xsum =',xsum
        write(*,*) 'ysum =',ysum
      end do trace_loop
      sigma = xsum
      dot_s = ysum
      write(*,*) 'sigma =', sigma
      write(*,*) 'dot_s =', dot_s

      if( sigma(mu) == 0.0e0_wp ) then 
        delta_ccurr(mu) = 0.0e0_wp
      else
        delta_ccurr(mu) = ( (sigma(mu)) / ((ecurr * smat_mu_mu) - (hmat_mu_mu)) ) 
      end if

      write(*,*) 'delta_ccurr =', delta_ccurr(mu)
      write(*,*) "sigma",  (sigma(mu)) 
      write(*,*) "deno", ((ecurr * smat_mu_mu) - (hmat_mu_mu))  
      write(*,*) "ecurr smumu hmumu", ecurr,smat_mu_mu,hmat_mu_mu
      delta_d = ( (2.0e0_wp * dot_s(mu)) + (smat_mu_mu * delta_ccurr(mu)) )
!      write(*,*) 'delta_d =', delta_d
      delta_e = sigma(mu) * delta_ccurr(mu) / ( Dsum + delta_d )
!      write(*,*) 'delta_e =',delta_e
      write(*,*) 'Dsum =',Dsum
      ecurr = ecurr + delta_e
      ccurr(mu) = ccurr(mu) + delta_ccurr(mu)
      write(*,*) '***************'
      write(*,*) 'ecurr =', ecurr
    end do mu_loop
!    if( delta_e <= 10e-6 ) stop
!  end do

end subroutine nesbet_diag
end program diagonalization
