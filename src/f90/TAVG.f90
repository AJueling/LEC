program TAVG
implicit none

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Program takes a list of binary files and averages over them. This is
!  used to create yearly averages from monthly average output of the
!  POP model.
!  by Andre Juling 2015
!
!  program needs the following input:
!  1. number of files to be averaged
!  2. dimensions of unformatted binary file =  (imt*jmt) * number of records
!  3. name of the output file
!  4. names of files to be averaged
!
!  a good webiste with information about Fortran file i/o:
!  http://cs.ubishops.ca/ljensen/fortran/fortran.htm#files
!
!  currently with the 40Gb output files, it takes about 5 min/file on Cartesius
!
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

integer                           :: imt, jmt, nrecs, length, i, j, k
double precision                  :: frac
character*120                     :: bin_file, new_file
real, dimension(:,:), allocatable :: NEW, RECORD

write (*,*) 'Enter number of files to be averaged:'
read  (*,*) length
frac = 1.0/real(length)
write (*,*) 'number of files:', length, ' -> frac =', frac
write (*,*) 'Enter dimensions of the binary files:'
read  (*,*) imt, jmt, nrecs
write(*,*) imt,jmt
allocate( NEW(imt,jmt), RECORD(imt,jmt) )

write (*,*) 'Enter the name of the new file:'
read  (*,'(a120)') new_file

open(1,file=new_file,access='direct',form='unformatted',recl=imt*jmt,status='unknown')

! loop over files to be averaged
do i = 1,length

  write(*,*) 'File number: ', i
  read(*,'(a120)') bin_file
  write(*,*) i, trim(bin_file)
  open(2,file=bin_file,access='direct',form='unformatted',recl=imt*jmt,status='old',action='read')

  ! loop over records 
  do k = 1,nrecs

    if ( i==1 ) then
      NEW(:,:) = 0.0
    else
      read  (1,rec=k) NEW
    endif

    read  (2,rec=k) RECORD
    NEW = NEW + frac*RECORD
    write (1,rec=k) NEW

  enddo
  close (2)
enddo

close (1)

end program TAVG
