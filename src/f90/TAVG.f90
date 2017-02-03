!
   program TAVG_BIN

!===============================================================================
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
!===============================================================================
   implicit none

!===============================================================================
!  variables
!===============================================================================

   integer                           :: imt, jmt, nrecs, length
   integer                           :: i,j,k,rec_length
   double precision                  :: frac
   character*120                     :: bin_file, new_file
   real, dimension(:,:), allocatable :: TAVG, RECORD, ZRECORD

   write (*,*) 'Enter number of files to be averaged:'
   read  (*,*) length
   frac = 1.0/real(length)
   write (*,*) 'number of files:', length, ' -> frac =', frac

   write (*,*) 'Enter dimensions of the binary files:'
   read  (*,*) imt, jmt, nrecs

   allocate( TAVG(imt,jmt), RECORD(imt,jmt), ZRECORD(imt,jmt) )
   do j = 1,jmt
   do i = 1,imt
     ZRECORD(i,j) = 0.0
   enddo
   enddo

   write (*,*) 'Enter the name of the new file:'
   read  (*,'(a120)') new_file

   inquire ( iolength = rec_length ) TAVG
   open(1,file=new_file,access='direct',form='unformatted',recl=rec_length,status='unknown')

   ! loop over files to be averaged
   do i = 1,length

     write (*,*) 'File number: ', i
     read  (*,'(a120)') bin_file

     open  (2,file=bin_file,access='direct',form='unformatted',recl=rec_length,status='old',action='read')

     ! loop over records 
     do k = 1,nrecs

     if ( i==1 ) then
     write (1,rec=k) ZRECORD
     endif

     read  (1,rec=k) TAVG
     read  (2,rec=k) RECORD

     TAVG = TAVG + frac*RECORD

     write (1,rec=k) TAVG

     enddo

     close (2)
   enddo

   close (1)


   end program
