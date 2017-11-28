cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine read_xy_init(iret)

      implicit none
      include 'df_coms.f'
      include 'ucoms.f'

      integer numPart,iret,index_xy
      double precision mass,sigma_HQ,dum_rx,dum_ry
      double precision rlu
      character*80 file30

      numXY=0

 2201 continue

!      read(unit=5,fmt=*,err=2204,end=2202) dum_rx,dum_ry
! Yingru (Read from file)
      call getenv('ftn30', file30)
      If (file30(1:4) .NE. '    ') Then
        OPEN(UNIT=30, FILE=file30, STATUS='OLD', FORM='FORMATTED')
      End If
     
      read(unit=30,fmt=*,err=2204,end=2202) dum_rx, dum_ry
!      write(6, *) dum_rx, dum_ry
! end of Yingru modify
      numXY=numXY+1
      initX(numXY)=dum_rx
      initY(numXY)=dum_ry

      if(numXY.lt.mxpart) then
         goto 2201
      else
         write(6,*) "XY table is full: ",numXY
         write(6,*) "End reading in initial positions."
         iret=1
         goto 2203
      endif

 2202 continue
      iret=1
      write(6,*) 'EOF reached in the table of initial positions.'
      write(6,*) 'Number of XY positions: ',numXY
      write(6,*) "Finish reading in initial positions."

 2203 continue 

      if(exp_setup.eq.1) then ! LHC 2.76~TeV
         sigma_pptot = 61.3564d0
         Ap = 208
         Zp = 82
         At = 208
         Zt = 82
         ecm = 2760d0
      elseif(exp_setup.eq.2) then ! RHIC 200~GeV
         sigma_pptot = 41.9357d0
         Ap = 197
         Zp = 79
         At = 197
         Zt = 79
         ecm = 200d0
      else
         write(*,*) "Unexpected experimental setup ..."
         write(*,*) "Terminating ..."
         stop
      endif

      if(reweight.ne.1) then
         sigma_ctot=sigma_c0*2d0*eta_cut
         sigma_btot=sigma_b0*2d0*eta_cut
      endif

      if(iflav.eq.4) then ! c quark
         sigma_HQ=sigma_ctot
      elseif(iflav.eq.5) then ! b quark
         sigma_HQ=sigma_btot
      else
         write(*,*) "Unexpected id of heavy quark ..."
         write(*,*) "Terminating ..."
         stop
      endif

      if(num_binary.eq.0) then 
         npart=numXY
      elseif(num_binary.lt.0) then
         npart=-num_binary
      else
         npart=int(sigma_HQ/sigma_pptot*num_binary)*2
      endif     

c generate initial table
      
      entry reSampleXY3

      numPart=0
      if(iflav.eq.4) then ! c quark
         mass=cMass
      elseif(iflav.eq.5) then ! b quark
         mass=bMass
      endif

      do while (numPart.lt.npart)
         numPart = numPart+1
         index_xy = int(rlu(0)*numXY)+1
         if(index_xy.gt.numXY) index_xy = numXY
         ityp(numPart) = iflav
         rx(numPart) = initX(index_xy)
         ry(numPart) = initY(index_xy)
         rz(numPart) = 0d0
         r0(numPart) = 0d6
         fmass(numPart) = mass
         call pQCDwt(ityp(numPart),px(numPart),py(numPart),pz(numPart),
     &           p0(numPart),fmass(numPart),pweight(numPart))

c now its anti-particle (if consider pair production)
         if(corr_flag.eq.1) then
            ityp(numPart+1) = -ityp(numPart)
            px(numPart+1) = -px(numPart)
            py(numPart+1) = -py(numPart)
            pz(numPart+1) = -pz(numPart)
            p0(numPart+1) = p0(numPart)
            fmass(numPart+1) = fmass(numPart)
            rx(numPart+1) =  rx(numPart)
            ry(numPart+1) =  ry(numPart)
            rz(numPart+1) =  rz(numPart)
            r0(numPart+1) =  r0(numPart)
            pweight(numPart+1) = pweight(numPart)
            numPart=numPart+1
         endif

      enddo
      return

 2204 continue
      write(6,*) 'READ-ERROR in weight table'
      write(6,*) 'terminating ...'
      stop
      return
      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine read_xyz_init(iret)

      implicit none
      include 'df_coms.f'
      include 'ucoms.f'

      integer numPart,iret,index_xyz
      double precision mass,sigma_HQ,dum_rx,dum_ry, dum_rz
      double precision rlu
      character*80 file30

      numXY=0

 2201 continue

!      read(unit=5,fmt=*,err=2204,end=2202) dum_rx,dum_ry
! Yingru (Read from file)
      call getenv('ftn30', file30)
      If (file30(1:4) .NE. '    ') Then
        OPEN(UNIT=30, FILE=file30, STATUS='OLD', FORM='FORMATTED')
      End If
     
      read(unit=30,fmt=*,err=2204,end=2202) dum_rx, dum_ry
!      write(6, *) dum_rx, dum_ry
! end of Yingru modify
      numXY=numXY+1
      initX(numXY)=dum_rx
      initY(numXY)=dum_ry
      initZ(numXY)=dum_rz

      if(numXY.lt.mxpart) then
         goto 2201
      else
         write(6,*) "XY table is full: ",numXY
         write(6,*) "End reading in initial positions."
         iret=1
         goto 2203
      endif

 2202 continue
      iret=1
      write(6,*) 'EOF reached in the table of initial positions.'
      write(6,*) 'Number of XY positions: ',numXY
      write(6,*) "Finish reading in initial positions."

 2203 continue 

      if(exp_setup.eq.1) then ! LHC 2.76~TeV
         sigma_pptot = 61.3564d0
         Ap = 208
         Zp = 82
         At = 208
         Zt = 82
         ecm = 2760d0
      elseif(exp_setup.eq.2) then ! RHIC 200~GeV
         sigma_pptot = 41.9357d0
         Ap = 197
         Zp = 79
         At = 197
         Zt = 79
         ecm = 200d0
      else
         write(*,*) "Unexpected experimental setup ..."
         write(*,*) "Terminating ..."
         stop
      endif

      if(reweight.ne.1) then
         sigma_ctot=sigma_c0*2d0*eta_cut
         sigma_btot=sigma_b0*2d0*eta_cut
      endif

      if(iflav.eq.4) then ! c quark
         sigma_HQ=sigma_ctot
      elseif(iflav.eq.5) then ! b quark
         sigma_HQ=sigma_btot
      else
         write(*,*) "Unexpected id of heavy quark ..."
         write(*,*) "Terminating ..."
         stop
      endif

      if(num_binary.eq.0) then 
         npart=numXY
      elseif(num_binary.lt.0) then
         npart=-num_binary
      else
         npart=int(sigma_HQ/sigma_pptot*num_binary)*2
      endif     

c generate initial table
      
      entry reSampleXY4

      numPart=0
      if(iflav.eq.4) then ! c quark
         mass=cMass
      elseif(iflav.eq.5) then ! b quark
         mass=bMass
      endif

      do while (numPart.lt.npart)
         numPart = numPart+1
         index_xyz = int(rlu(0)*numXY)+1
         if(index_xyz.gt.numXY) index_xyz = numXY
         ityp(numPart) = iflav
         rx(numPart) = initX(index_xyz)
         ry(numPart) = initY(index_xyz)
         rz(numPart) = initZ(index_xyz)
         r0(numPart) = 0d6
         fmass(numPart) = mass
         call pQCDwt(ityp(numPart),px(numPart),py(numPart),pz(numPart),
     &           p0(numPart),fmass(numPart),pweight(numPart))

c now its anti-particle (if consider pair production)
         if(corr_flag.eq.1) then
            ityp(numPart+1) = -ityp(numPart)
            px(numPart+1) = -px(numPart)
            py(numPart+1) = -py(numPart)
            pz(numPart+1) = -pz(numPart)
            p0(numPart+1) = p0(numPart)
            fmass(numPart+1) = fmass(numPart)
            rx(numPart+1) =  rx(numPart)
            ry(numPart+1) =  ry(numPart)
            rz(numPart+1) =  rz(numPart)
            r0(numPart+1) =  r0(numPart)
            pweight(numPart+1) = pweight(numPart)
            numPart=numPart+1
         endif

      enddo
      return

 2204 continue
      write(6,*) 'READ-ERROR in weight table'
      write(6,*) 'terminating ...'
      stop
      return
      end

