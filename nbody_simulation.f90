program nbody_simulation
 implicit none
 real,allocatable ::a(:,:),r(:,:),dr(:),v(:,:),m(:),orbiters(:,:)
 integer :: N,t,dt,ios,i,k

!opens and reads the input file containing information about the distances, velocities and masses of the planets
 open(unit=8,file='input2.dat',status='old',iostat=ios)

 	if(ios .ne. 0) then
 	 print*,"Error opening the input file."
     stop
 	else
     read(8,*) N
 	 read(8,*) k
 	 read(8,*) dt
 	 read(8,*) t
     print*,t
 	 allocate(orbiters(7,N),r(3,N),a(3,N),dr(3),v(3,N),m(N))
	 do i=1,N
      read(8,*) orbiters(:,i)
   	 end do
    end if
 close(8)

 r(1:3,:)=orbiters(1:3,:)
 v(1:3,:)=orbiters(4:6,:)
 m(:)=orbiters(7,:)

!opens the output file to record the generated data
open(unit=25,file='output5.dat',status='old',iostat=ios)
 if(ios .ne. 0) then
   print*,"Error opening the output file."
   stop
 end if
 
 call verletint(t,dt,r,v,m,N,k)  !the subroutine verletint is defined below
close(25)

contains

subroutine verletint(t,dt,r,v,m,N,k)
 integer,intent(in) :: N,k,t,dt
 real,intent(inout) :: r(3,N),v(3,N),m(N)
 real :: dr(3),r_mag,a(3,N),anew(3,N), m2AU
 real, parameter :: G=6.6738D-11
 integer :: counter,i,j,steps
 steps=t/dt
 counter=0
 m2AU=6.68458712e-12

 a=0
 anew=0

 do while(counter .le. steps)
  do i=1,N
   do j=1,N
    if(i .le. j) cycle   !ignores its own gravitational effect
    dr(1:3)=r(1:3,j)-r(1:3,i)
    r_mag=sqrt(dr(1)**2+dr(2)**2+dr(3)**2)

    anew(1:3,i)=anew(1:3,i)+G*m(j)*dr(1:3)/r_mag**3  !new gravitational acceleration
    
   end do
   v(1:3,i)=v(1:3,i)+0.5*(a(1:3,i)+anew(1:3,i))*dt
   a(1:3,i)=anew(1:3,i)
   r(1:3,i)=r(1:3,i)+v(1:3,i)*dt+0.5*a(1:3,i)*dt**2

   anew=0

  end do

   if(mod(counter,k) .eq. 0) then 
    do i=1,N
     write(25,*) r(1:3,i)*m2AU
    end do
   end if

!output
   print 10,'# objects:',N
   print 10,'# steps:',counter
   print'(a)','The x,y-coordinates of each body in the Solar System in astronomical units'
   print 11,'Sun      	:',r(1:2,1)*m2AU
   print 11,'Mercury  	:',r(1:2,2)*m2AU
   print 11,'Venus  	:',r(1:2,3)*m2AU
   print 11,'Earth  	:',r(1:2,4)*m2AU
   print 11,'Mars     :',r(1:2,5)*m2AU
   !print 11,'Jupiter  :',r(1:2,6)*m2AU
   !print 11,'Saturn :',r(1:2,7)*m2AU
   !print 11,'Uranus   :',r(1:2,8)*m2AU
   !print 11,'Neptune :',r(1:2,9)*m2AU
   print'(/)'

!if needed to extend to 8 planets, remove the exclamation from the last four entries above

   10 format (a,1x,i0)
   11 format (a,1x,f15.10,f15.10)
  counter=counter+1
 end do

end subroutine verletint
end program nbody_simulation