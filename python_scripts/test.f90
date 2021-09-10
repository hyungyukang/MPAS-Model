program test
implicit none
integer,parameter :: nx = 121,ncpu = 60
real(kind=8),dimension(nx,nx) :: matp
integer :: i,j,k,n,rank
character(100) :: iname

do n = 1,ncpu

  write(iname,('(A,i3.3,A)')) 'output/matp_',n-1,'.csv'
  print*, trim(iname)

  open(11,file=trim(iname),action='read')
    do j = 1,nx
      read(11,*) matp(j,:)
    end do
  close(11)

end do


print*, matp(1,:)


end program test
