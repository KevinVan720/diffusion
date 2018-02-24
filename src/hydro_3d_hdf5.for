!**********************************************************
      Subroutine readHydroFiles_initial_3D(H5HydroFilename_in, 
     &                                  bufferSize_in)

      Use HDF5
      Implicit none

      Character(Len=15) :: H5HydroFilename_in
      Character(Len=15) :: hydroFileH5name ! Filename
      Character(Len=8) :: groupEventname = "/Event" ! group name

      Common /fileInfo3D/ hydroFileH5name, groupEventname

      Integer(HID_T) :: file_id   ! file identifier
      Integer(HID_T) :: group_id  ! group identifier

      Integer :: bufferSize_in
      Integer :: error

      hydroFileH5name = H5HydroFilename_in
      
      ! initialize FORTRAN interface
      call h5open_f(error)

      ! open hydro file
      call h5fopen_f (hydroFileH5name, H5F_ACC_RDWR_F, file_id, error)

      ! open a group
      call h5gopen_f (file_id, groupEventname, group_id, error)

      ! read Attribute for group "Event"
      call readHydroGridInfo_3D(group_id)
      call printHydroGridInfo_3D()

      ! read datasets from the file
      call readHydroinfoBuffered_initialization_3D(bufferSize_in)
      call readHydroinfoBuffered_total_3D(group_id)

      ! close the groups
      call h5gclose_f(group_id, error)

      ! close the file
      call h5fclose_f(file_id, error)

      ! close FORTRAN interface
      call h5close_f(error)
      end

!----------------------------------------------------------------------


!*********************************************************************
      Subroutine readHydroGridInfo_3D(group_id)
      use HDF5
      Implicit None

      character(Len=15) :: hydroFileH5name
      character(len=8) :: groupEventname

      Common /fileInfo3D/ hydroFileH5name, groupEventname

      Integer:: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Integer:: hydroGrid_ZL, hydroGrid_ZH
      Integer:: hydroGrid_NX, hydroGrid_NY, hydroGrid_NZ   ! here z is actually eta
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_xmax, hydroGrid_ymax, hydroGrid_zmax
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes

      Integer(HID_T) :: group_id
      Integer :: error

      call readH5Attribute_int3D(group_id, "nx", hydroGrid_NX)
      call readH5Attribute_int3D(group_id, "ny", hydroGrid_NY)
      call readH5Attribute_int3D(group_id, "nz", hydroGrid_NZ)
      
      hydroGrid_XL = - int(hydroGrid_NX/2)
      hydroGrid_XH = int(hydroGrid_NX/2)
      hydroGrid_YL = -int(hydroGrid_NY/2)
      hydroGrid_YH = int(hydroGrid_NY/2)
      hydroGrid_ZL = - int(hydroGrid_NZ/2)
      hydroGrid_ZH = int(hydroGrid_NZ/2)
    
      call readH5Attribute_double3D(group_id,"xmin", hydroGrid_xmin)
      call readH5Attribute_double3D(group_id,"ymin", hydroGrid_ymin)
      call readH5Attribute_double3D(group_id,"zmin", hydroGrid_zmin)
      call readH5Attribute_double3D(group_id,"xmax", hydroGrid_xmax)
      call readH5Attribute_double3D(group_id,"ymax", hydroGrid_ymax)
      call readH5Attribute_double3D(group_id,"zmax", hydroGrid_zmax)
      call readH5Attribute_double3D(group_id,"dx", hydroGrid_dx)
      call readH5Attribute_double3D(group_id,"dy", hydroGrid_dy)
      call readH5Attribute_double3D(group_id,"dz", hydroGrid_dz)
      call readH5Attribute_double3D(group_id,"tau0", hydroGrid_tau0)
      call readH5Attribute_double3D(group_id,"dtau", hydroGrid_dtau)

      call h5gn_members_f(group_id, groupEventname, 
     &              hydroGrid_numOfframes, error)
      
      hydroGrid_taumax=hydroGrid_tau0
     &              +(hydroGrid_numOfframes-1)*hydroGrid_dtau

!debug
!      write(6,*) hydroGrid_NX,hydroGrid_NY,hydroGrid_NZ,
!     &  hydroGrid_XL,hydroGrid_YL,hydroGrid_ZL,
!     &  hydroGrid_XH,hydroGrid_YH,hydroGrid_ZH,
!     &  hydroGrid_xmin,hydroGrid_ymin,hydroGrid_zmin,
!     &  hydroGrid_xmax,hydroGrid_ymax,hydroGrid_zmax,
!     &  hydroGrid_dx,hydroGrid_dy,hydroGrid_dz,
!     &  hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax


      end
!---------------------------------------------------------------------

!********************************************************************
      Subroutine printHydroGridInfo_3D()
      Implicit none

      character(Len=15) :: hydroFileH5name
      character(len=8) :: groupEventname
      Common /fileInfo3D/ hydroFileH5name, groupEventname

      Integer:: hydroGrid_XL, hydroGrid_YL, hydroGrid_ZL   ! here z is actually eta
      Integer:: hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes


      write(*,'(A)') "----------------------------------------------"
      write(*,'(A)') "------------hydro grid info-------------------"
      write(*,'(A)') "----------------------------------------------"
      write(*,'(A, A)')"Filename : ", hydroFileH5name
      write(*,'(A,I5, I5, F6.3,A)') "XL,XH = ", hydroGrid_XL,
     &             hydroGrid_XH, hydroGrid_dx, "fm"
      write(*,'(A,I5, I5, F6.3,A)') "YL,YH = ", hydroGrid_YL,
     &             hydroGrid_YH, hydroGrid_dy, "fm"
      write(*,'(A,I5, I5, F6.3, A)') "etaL, etaH = ", hydroGrid_ZL,
     &             hydroGrid_ZH, hydroGrid_dz, "fm"
      write(*,'(A, F5.3, A)')"Tau0 = ",hydroGrid_tau0,"fm/c"
      write(*,'(A, F5.3, A)')"dTau = ",hydroGrid_dTau,"fm/c"
      write(*,'(A, I5)') "number of Frames = ",hydroGrid_numOfframes
      write(*,'(A, F7.3, A)')"Tau_max =", hydroGrid_Taumax,"fm/c"
      write(*,'(A)')"--------------------------------------------------"

      end
!-----------------------------------------------------------------------



!***********************************************************************
      Subroutine readH5Attribute_int3D(group_id, aname, avalue)
      Use HDF5
      Implicit none

      CHARACTER(LEN=*) :: aname       ! Attribute name
      Integer :: avalue      ! Attribute value
      INTEGER(HID_T) :: group_id      ! Group identifier

      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      Integer(HID_T) :: atype_id

      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank

      INTEGER     ::   error ! Error flag

      ! open an attribute
      Call h5aopen_name_f(group_id, aname, attr_id, error)
      ! read an attribute
      Call h5aread_f(attr_id, H5T_NATIVE_INTEGER, avalue, adims, error)
      ! close an attribute
      Call h5aclose_f(attr_id, error)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readH5Attribute_double3D(group_id, aname, avalue)
      Use HDF5
      Implicit none

      CHARACTER(LEN=*) :: aname       ! Attribute name
      Double precision :: avalue      ! Attribute value
      INTEGER(HID_T) :: group_id      ! Group identifier

      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      Integer(HID_T) :: atype_id

      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank

      INTEGER     ::   error ! Error flag

      ! open an attribute
      Call h5aopen_name_f(group_id, aname, attr_id, error)
      ! read an attribute
      Call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, avalue, adims, error)
      ! close an attribute
      Call h5aclose_f(attr_id, error)

      end
!-----------------------------------------------------------------------



!***********************************************************************
      Subroutine readHydroinfoBuffered_initialization_3D(bufferSize_in)
      Implicit None

      Integer:: hydroGrid_XL, hydroGrid_YL, hydroGrid_ZL   ! here z is actually eta
      Integer:: hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes


      Integer:: bufferSize
      Integer:: bufferSize_in

      ! the last index is the buffer "layer" index that goes from 1 to buffersize
      Double Precision, Pointer::
     &  TM(:,:,:,:), vxM(:,:,:,:), vyM(:,:,:,:), vzM(:,:,:,:)

      Common /bufferedData3D/ bufferSize,
     &  TM, vxM, vyM, vzM

      bufferSize = bufferSize_in
      Allocate(TM(hydroGrid_XL:hydroGrid_XH,hydroGrid_YL:hydroGrid_YH,
     &            hydroGrid_ZL:hydroGrid_ZH, 1:bufferSize))
      Allocate(vxM(hydroGrid_XL:hydroGrid_XH,hydroGrid_YL:hydroGrid_YH,
     &            hydroGrid_ZL:hydroGrid_ZH, 1:bufferSize))
      Allocate(vyM(hydroGrid_XL:hydroGrid_XH,hydroGrid_YL:hydroGrid_YH,
     &            hydroGrid_ZL:hydroGrid_ZH, 1:bufferSize))
      Allocate(vzM(hydroGrid_XL:hydroGrid_XH,hydroGrid_YL:hydroGrid_YH,
     &           hydroGrid_ZL:hydroGrid_ZH, 1:bufferSize))


      end
!-----------------------------------------------------------------------------



!****************************************************************************
      Subroutine readHydroinfoBuffered_total_3D(group_id)
      use HDF5
      Implicit None

      character(Len=15) :: hydroFileH5name
      character(len=8) :: groupEventname
      Common /fileInfo3D/ hydroFileH5name, groupEventname

      character(len=10) :: frameName  !group frame name
      character(len=4) :: frame_id_string 
    
      Integer(HID_T) :: group_id
      Integer(HID_T) :: groupFrame_id
      Integer :: error

      Integer:: hydroGrid_XL, hydroGrid_YL, hydroGrid_ZL   ! here z is actually eta
      Integer:: hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes


 
      Double Precision, Dimension(hydroGrid_XL:hydroGrid_XH,
     &     hydroGrid_YL:hydroGrid_YH, hydroGrid_ZL:hydroGrid_ZH,1:1)
     &  :: Temp, Vx, Vy, Vz


      Double Precision :: Time
      Integer:: bufferSize

      Double Precision, Pointer::
     &  TM(:,:,:,:), vxM(:,:,:,:), vyM(:,:,:,:), vzM(:,:,:,:)

      Common /bufferedData3D/ bufferSize,
     &  TM, vxM, vyM, vzM

       Integer :: J
       Integer :: xidx


       if (bufferSize.lt.hydroGrid_numOfframes) then
         write(*,*) "BufferSize if too small, increase it to at least",
     &      hydroGrid_numOfframes
         stop
       endif

       Do J=1, hydroGrid_numOfframes
         write(unit=frame_id_string, fmt='(I3.3)') J-1
         frameName = "Frame_" // frame_id_string
         call h5gopen_f(group_id, frameName, groupFrame_id, error)

         call readH5Dataset_double3D(groupFrame_id, 'Temp', Temp)
         call readH5Dataset_double3D(groupFrame_id, 'vx', Vx)
         call readH5Dataset_double3D(groupFrame_id, 'vy', Vy)
         call readH5Dataset_double3D(groupFrame_id, 'vz', Vz)

         ! write to the buffer
         TM(:,:,:,J) = Temp(:,:,:,1)
         vxM(:,:,:,J) = Vx(:,:,:,1)
         vyM(:,:,:,J) = Vy(:,:,:,1)
         vzM(:,:,:,J) = Vz(:,:,:,1)
        
        ! close the groups
        call h5gclose_f(groupFrame_id, error)
      enddo

      End
!-------------------------------------------------------------------------



!************************************************************************
      Subroutine readH5Dataset_double3D(group_id, datasetName,dset_data)
      Use HDF5
      Implicit None

      Character(Len=*) :: datasetName
      Integer(HID_T) :: group_id
      Integer(HID_T) :: dset_id

      Integer:: hydroGrid_XL, hydroGrid_YL, hydroGrid_ZL   ! here z is actually eta
      Integer:: hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes


      Integer(HSIZE_T), Dimension(3) :: data_dims
      Integer(HSIZE_T), Dimension(3) :: data_dims_Cstyle

      Double precision, Dimension(hydroGrid_XL:hydroGrid_XH,
     & hydroGrid_YL:hydroGrid_YH,hydroGrid_ZL:hydroGrid_ZH,1:1)
     &  :: dset_data
      Double precision, Dimension(hydroGrid_ZL:hydroGrid_ZH,
     & hydroGrid_YL:hydroGrid_YH,hydroGrid_XL:hydroGrid_XH,1:1)
     &  :: dset_data_Cstyle
      Integer :: error
      Integer :: i, j, k

      ! read in data matrix assuming in c style
      ! need to perform transport to convert it into fortran stle
      data_dims(1) = hydroGrid_XH - hydroGrid_XL +1
      data_dims(2) = hydroGrid_YH - hydroGrid_YL + 1
      data_dims(3) = hydroGrid_ZH - hydroGrid_ZL +1

      data_dims_Cstyle(1) = data_dims(3)
      data_dims_Cstyle(2) = data_dims(2)
      data_dims_Cstyle(3) = data_dims(1)

      ! open an existing dataset
      call h5dopen_f(group_id, datasetName, dset_id, error)

      ! read the dataset
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset_data_Cstyle,
     &      data_dims_Cstyle, error)
      
      do i = hydroGrid_XL, hydroGrid_XH, 1
        do j = hydroGrid_YL, hydroGrid_YH, 1
          do k = hydroGrid_ZL, hydroGrid_ZH, 1
            dset_data(i,j,k,1) = dset_data_Cstyle(k,j,i,1)    !!!!!     ATTENTTION (CHECK THIS)
                                                              
          enddo
        enddo
      enddo

      ! close the dataset
      call h5dclose_f(dset_id, error)
      end
!---------------------------------------------------------------------------------



!********************************************************************************
      Subroutine readHydroInfoYingru_3D(t,x,y,z,Temp,vx,vy,vz,ctl)
! t, x, y, z: input (cartesian corrdinate)
! Temp, vx, vy, vz, ctl: output
! ctl: 0-success, 1-x/y/z out of range, 2-tau>maxTau,2-tau<tau0

      Implicit None
      double precision:: t,x,y,z,Temp,vx,vy,vz
      Integer:: ctl
      Double precision:: tau, eta, bfactor,gamma

      Integer:: hydroGrid_XL, hydroGrid_YL, hydroGrid_ZL   ! here z is actually eta
      Integer:: hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes

      Temp=0d0
      vx=0d0
      vy=0d0
      vz=0d0

      ctl = 0
      if((t*t-z*z).lt.0) then
        write(6,*) "Warning! t**2-z**2<0 ..."
        tau=0
        eta=0d0
      else
        tau =sqrt(t*t-z*z)
        eta=0.5*Log((t+z)/(t-z))
      endif
      
      if(tau.lt.hydroGrid_tau0) then
        ctl=3
        return
      endif

      if(tau.gt.hydroGrid_taumax) then
        ctl=2
        return 
      endif

      ! assume symmetric grid
      if(abs(x).gt.abs(hydroGrid_xmin)
     & .or. abs(y).gt.abs(hydroGrid_ymin)
     & .or. abs(eta).gt.abs(hydroGrid_zmin)) then
        ctl = 1
        return
      endif
       
!debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       write(6,*) "tracking-hydro-3d-1: ", tau,eta,Temp,vx,vy,vz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      call readHydroinfoBuffered_ideal_3D(tau,x,y,eta,Temp,vx,vy,vz)
!!! debug
!!! this is a 3D hydro but mimic 2D behavior!
!      if (abs(eta) .lt. 2d0) then
!        vz = z/(t+1D-30) * 1.3
!      else
!       vz = z/(t+1D-30)
!      endif
!      gamma = 1d0/(sqrt(1d0-vz*vz) + 1D-30)
!      vx = vx/gamma
!      vy = vy/gamma


!### not really correct...
      if (abs(vz).gt. (abs(z)+0.2)/(t+1D-30)) then
        vz = sign((abs(z)+0.2)/(t+1D-30), z)
      endif
   

!!! a real 3D hydro
      bfactor = vx**2 + vy**2 + vz**2
      if (bfactor.gt.1.0) then
        write(6,*) "fluid velocity exceeds 1!"
        vx = vx/sqrt(bfactor)
        vy = vy/sqrt(bfactor)
        vz = vz/sqrt(bfactor)
      endif
       
!debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       write(6,*) "tracking-hydro-3d: ", tau,eta,Temp,vx,vy,vz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end subroutine
!-----------------------------------------------------------------------------





!***********************************************************************
      Subroutine readHydroinfoBuffered_ideal_3D(tau,
     &                          x,y,eta,Temp,vx,vy,vz)
      Implicit None

      Double Precision :: tau,x,y,eta,Temp,vx,vy,vz

      Integer:: hydroGrid_XL, hydroGrid_YL, hydroGrid_ZL   ! here z is actually eta
      Integer:: hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes

      Double Precision, Dimension(1:2,1:2,1:2) :: Temp1,VxB1,VyB1,VzB1
      Double Precision, Dimension(1:2,1:2,1:2) :: Temp2,VxB2,VyB2,VzB2

      Integer :: tauI ! tau in tau0+dtau*tauI, tau0+dtau*(tauI+1))
      Double Precision:: tauInc ! tau = tau0+dtau*(tauI+tauInc)

      Integer:: xi, yi, zi
      Double Precision:: xInc, yInc, zInc

      Double Precision :: var1 ! temporary variables

      if (tau<hydroGrid_tau0.and.tau>hydroGrid_tau0-1d-15) then
        tau = hydroGrid_tau0
      endif

      var1 = (tau-hydroGrid_tau0)/hydroGrid_dtau
      tauI = floor(var1)
      tauInc = var1 - tauI

      var1 = (x)/hydroGrid_dx
      xi = floor(var1)
      xInc = var1 - xi

      var1 = (y)/hydroGrid_dy
      yi = floor(var1)
      yInc = var1 - yi

      var1 = (eta)/hydroGrid_dz
      zi = floor(var1)
      zInc = var1 - zi

!!!debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write(6,*) "tracking-hydro-read: ", tau, tauI, tauInc
!     &  ,x, xi, xInc, y, yi, yInc, eta, zi, zInc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call readHydroBlockBufferedOrdered_ideal_3D
     &  (tauI+1,xi,yi,zi,Temp1,VxB1,VyB1,VzB1)

      call readHydroBlockBufferedOrdered_ideal_3D
     &  (tauI+2,xi,yi,zi,Temp2,VxB2,VyB2,VzB2)


     
      call cubeInterp4D(tauInc,xInc,yInc,zInc,Temp,
     & Temp1(1,1,1),Temp1(2,1,1),Temp1(1,2,1),Temp1(1,1,2),
     & Temp1(1,2,2),Temp1(2,1,2),Temp1(2,2,1),Temp1(2,2,2),
     & Temp2(1,1,1),Temp2(2,1,1),Temp2(1,2,1),Temp2(1,1,2),
     & Temp2(1,2,2),Temp2(2,1,2),Temp2(2,2,1),Temp2(2,2,2))

      call cubeInterp4D(tauInc,xInc,yInc,zInc,vx,
     & VxB1(1,1,1),VxB1(2,1,1),VxB1(1,2,1),VxB1(1,1,2),
     & VxB1(1,2,2),VxB1(2,1,2),VxB1(2,2,1),VxB1(2,2,2),
     & VxB2(1,1,1),VxB2(2,1,1),VxB2(1,2,1),VxB2(1,1,2),
     & VxB2(1,2,2),VxB2(2,1,2),VxB2(2,2,1),VxB2(2,2,2))

       call cubeInterp4D(tauInc,xInc,yInc,zInc,vy,
     & VyB1(1,1,1),VyB1(2,1,1),VyB1(1,2,1),VyB1(1,1,2),
     & VyB1(1,2,2),VyB1(2,1,2),VyB1(2,2,1),VyB1(2,2,2),
     & VyB2(1,1,1),VyB2(2,1,1),VyB2(1,2,1),VyB2(1,1,2),
     & VyB2(1,2,2),VyB2(2,1,2),VyB2(2,2,1),VyB2(2,2,2))

     
      call cubeInterp4D(tauInc,xInc,yInc,zInc,vz,
     & VzB1(1,1,1),VzB1(2,1,1),VzB1(1,2,1),VzB1(1,1,2),
     & VzB1(1,2,2),VzB1(2,1,2),VzB1(2,2,1),VzB1(2,2,2),
     & VzB2(1,1,1),VzB2(2,1,1),VzB2(1,2,1),VzB2(1,1,2),
     & VzB2(1,2,2),VzB2(2,1,2),VzB2(2,2,1),VzB2(2,2,2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if (Temp.ge.0.75) then
!        write(6,*) "debug!!..",tau,tauI,tauInc,x,xi,xInc,
!     &   y,yi,yInc,eta,zi,zInc,Temp
!      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write(6,*) "debug: readHydroInfo(T,vx,vy,vz): ", Temp,vx,vy,vz
      end subroutine
!--------------------------------------------------------------------




!******************************************************************
      Subroutine readHydroBlockBufferedOrdered_ideal_3D
     &  (idxTau,idxX,idxY,idxZ,Temp23,VxB23,VyB23,VzB23)
! read only the 2x3 block from the buffer
      Implicit None
      Integer :: idxTau, idxX, idxY, idxZ
      
      Integer:: hydroGrid_XL, hydroGrid_YL, hydroGrid_ZL   ! here z is actually eta
      Integer:: hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH
      Double precision:: hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin
      Double precision:: hydroGrid_dx, hydroGrid_dy, hydroGrid_dz
      Double precision:: hydroGrid_tau0,hydroGrid_dtau,hydroGrid_taumax
      Integer:: hydroGrid_numOfframes

      Common /hydroInfo3D/ hydroGrid_XL, hydroGrid_YL,hydroGrid_ZL,
     &          hydroGrid_XH, hydroGrid_YH, hydroGrid_ZH,
     &          hydroGrid_xmin, hydroGrid_ymin, hydroGrid_zmin,
     &          hydroGrid_dx, hydroGrid_dy, hydroGrid_dz,
     &          hydroGrid_tau0, hydroGrid_dtau, hydroGrid_taumax,
     &          hydroGrid_numOfframes

      Double precision, Dimension(hydroGrid_XL:hydroGrid_XH, 
     & hydroGrid_YL:hydroGrid_YH, hydroGrid_ZL:hydroGrid_ZH,
     & 1:1) :: Temp, VxB, VyB, VzB

      Double precision, Dimension(1:2,1:2,1:2) ::
     & Temp23, VxB23, VyB23, VzB23

      Integer:: bufferSize

      Double Precision, Pointer::
     & TM(:,:,:,:), vxM(:,:,:,:), vyM(:,:,:,:), vzM(:,:,:,:)

      Common /bufferedData3D/ bufferSize,
     &  TM, vxM, vyM, vzM

      If (idxTau.lt.0 .or. idxTau.gt.bufferSize 
     & .or. idxTau > hydroGrid_numOfframes) then
        Temp23(:,:,:) = 0d0
        VxB23(:,:,:) = 0d0
        VyB23(:,:,:) = 0d0
        VzB23(:,:,:) = 0d0
        Return
      End if
      
      If (idxX.lt.hydroGrid_XL .or. idxX.gt.hydroGrid_XH) then
        Temp23(:,:,:) = 0d0
        VxB23(:,:,:) = 0d0
        VyB23(:,:,:) = 0d0
        VzB23(:,:,:) = 0d0
        Return
      end if

      If (idxY.lt.hydroGrid_YL .or. idxY.gt.hydroGrid_YH) then
        Temp23(:,:,:) = 0d0
        VxB23(:,:,:) = 0d0
        VyB23(:,:,:) = 0d0
        VzB23(:,:,:) = 0d0
        Return
      end if


       If (idxZ.lt.hydroGrid_ZL .or. idxZ.gt.hydroGrid_ZH) then
        Temp23(:,:,:) = 0d0
        VxB23(:,:,:) = 0d0
        VyB23(:,:,:) = 0d0
        VzB23(:,:,:) = 0d0
        Return
      end if

       ! read it from buffer
       Temp23(:,:,:) = TM(idxX:idxX+1,idxY:idxY+1,idxZ:idxZ+1,idxTau)
       VxB23(:,:,:) = vxM(idxX:idxX+1,idxY:idxY+1,idxZ:idxZ+1,idxTau)
       VyB23(:,:,:) = vyM(idxX:idxX+1,idxY:idxY+1,idxZ:idxZ+1,idxTau)
       VzB23(:,:,:) = vzM(idxX:idxX+1,idxY:idxY+1,idxZ:idxZ+1,idxTau)

!debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write(6,*) "debug-read-from-buffer: ", idxX, idxY, idxZ, idxTau
!     & , TM(idxX,idxY,idxZ,idxTau)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       End Subroutine
!--------------------------------------------------------------------------



!******************************************************************
      Subroutine cubeInterp4D(t,x,y,z,Atxyz,
     & A0000,A0100,A0010,A0001,A0011,A0101,A0110,A0111,
     & A1000,A1100,A1010,A1001,A1011,A1101,A1110,A1111)
! perform 4D interpolate
! it is actually twice 3d interpolation and then 2d-linear cube
      Implicit None
      Double Precision:: t,x,y,z,Atxyz
      Double Precision:: A0000,A0100,A0010,A0001,A0011,
     &             A0101,A0110,A0111,A1000,A1100,A1010,
     &             A1001,A1011,A1101,A1110,A1111

      Double Precision:: Axyz1, Axyz2

      call cubeInterp3D(x,y,z,Axyz1, 
     &   A0000,A0100,A0010,A0001,A0011,A0101,A0110,A0111)

      call cubeInterp3D(x,y,z,Axyz2,
     &   A1000,A1100,A1010,A1001,A1011,A1101,A1110,A1111)

      Atxyz = t*(Axyz2 - Axyz1) + Axyz1

!!!! debug!!!!!!!!!!!!!!!!!!!!      
!      write(6,*) "debug-hydro...",A0000,A0100,A0010,A0001,A0011,
!     &             A0101,A0110,A0111,A1000,A1100,A1010,
!     &             A1001,A1011,A1101,A1110,A1111

!      write(6,*) "debug-hydro...",t, x,y,z,Axyz1, Axyz2, Atxyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end
!------------------------------------------------------------------------


!*******************************************************************
      Subroutine cubeInterp3D(x,y,z,Axyz,
     &     A000,A100,A010,A001,A011,A101,A110,A111)
! perform 3D interpolateion
      Implicit none
      Double Precision :: x,y,z,Axyz
      Double Precision :: A000,A100,A010,A110,A001,A101,A011,A111

      
      Axyz=A000*(1-x)*(1-y)*(1-z) + A100*x*(1-y)*(1-z)
     &      +A010*(1-x)*y*(1-z) + A001*(1-x)*(1-y)*z
     &      +A011*(1-x)*y*z + A101*x*(1-y)*z
     &      +A110*x*y*(1-z) + A111*x*y*z


!!!debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Axyz = Max(A000,A100,A010,A001,A110,A101,A011,A111)
        !Axyz = A000

!       if(A000.eq.0 .or.A100.eq.0 .or. A010.eq.0 .or. A001.eq.1
!     &  .or. A011.eq.0 .or.A101.eq.0 .or.A110.eq.0 .or.A111.eq.1) then
!      write(6,*) "debug: interpolate! ", x, y, z,
!     & A000, A010, A100, A110, A001, A011, A101, A111
!       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end
!--------------------------------------------------------------------

