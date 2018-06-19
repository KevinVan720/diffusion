cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_PHSD_init(iret)

      implicit none
      include 'ucoms.f'
      include 'df_coms.f'

      integer numPart, iret, i
      double precision dum_rx, dum_ry, dum_rz, dum_r0
      double precision dum_px, dum_py, dum_pz, dum_E, dum_ipT
      double precision mass
      character*80 file30
   
      call getenv('ftn30', file30)

      if (file30(1:4).ne.'    ') then
        open(unit=30, file=file30, status='old')
      endif
    
      numXY = 0
 
 1101 continue

      numXY = numXY+1
      read(30, *, err=1199, end=1103) dum_rx, dum_ry, dum_rz,
     &     dum_r0, dum_px, dum_py, dum_pz, dum_E, dum_ipT
      initX(numXY) = dum_rx
      initY(numXY) = dum_ry
      initZ(numXY) = dum_rz
      initZ0(numXY) = dum_r0
      initPX(numXY) = dum_px
      initPY(numXY) = dum_py
      initPZ(numXY) = dum_pz
      initE0(numXY) = sqrt(dum_px**2+dum_py**2+dum_pz**2+cMass**2)
      initIPT(numXY) = dum_ipT

      if(numXY .lt.mxpart) then
        goto 1101
      else
        write(6, *) "Initial table is full: ", numXY
        goto 1103
      endif


 1103 continue
      iret = 1

      numXY = numXY -1  ! substract the extra 1 which is added in last step

      if(num_binary.eq.0) then ! use numXY as npart
        npart = numXY
      elseif (num_binary.lt.0) then
        npart = -num_binary
      endif

      write(6,*) 'EOF reached in table of initial positions.'
      write(6,*) 'number of PHSD initial charm: ', npart
        
      numPart = 0
      if(iflav.eq.4) then ! cquark
        mass = cMass
      elseif(iflav.eq.5) then ! bquark
        mass = bMass
      endif

c let us now do not run random ic
      do while (numPart.lt.npart)
        numPart = numPart + 1
        ityp(numPart) = iflav
        rx(numPart) = initX(numPart)
        ry(numPart) = initY(numPart)
        rz(numPart) = initZ(numPart)
        r0(numPart) = initZ0(numPart)
        fmass(numPart) = mass
        px(numPart) = initPX(numPart)
        py(numPart) = initPY(numPart)
        pz(numPart) = initPZ(numPart)
        p0(numPart) = initE0(numPart)
        pweight(numPart) = initIPT(numPart)
      enddo
      return
     
 1199 continue
      write(6,*) "Read-error in table"
      write(6,*) 'terminating ...'
      stop 
      return

      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  read static medium initial condition
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_static_init(iret)
      implicit none
      include 'df_coms.f'
      include 'ucoms.f'

      integer numPart, iret, index_xy
      double precision mass, dum_rx, dum_ry
      double precision rlu

      write(6,*) "read in static init"
      numXY = 0
 2011 continue
      numXY = numXY + 1
      initX(numXY) = 0d0
      initY(numXY) = 0d0

      if (numXY.lt.mxpart) then
        goto 2011
      else
        write(6,*) "read_static_init: XY table is full", numXY
        iret = 1
        goto 2012
      endif

 2012 continue
      
      if (num_binary.eq.0) then
        npart = numXY
      elseif (num_binary .lt. 0) then
        npart = -num_binary
      endif

c generate initial table
      entry reSampleXY4

      numPart=0
      if (iflav.eq.4) then
        mass = cMass
      elseif (iflav.eq.5) then
        mass = bMass
      endif

      do while (numPart.lt.npart)
        numPart = numPart + 1
        ityp(numPart) = iflav
        rx(numPart) = 0d0
        ry(numPart) = 0d0
        rx(numPart) = 0d0
        r0(numPart) = 0d0
        fmass(numPart) = mass

! senario one: same initial momentum for all the particle
      px(numPart) = 0d0
      py(numPart) = 0d0
      pz(numPart) = p_static
      p0(numPart) = sqrt(pz(numPart)**2 + fmass(numPart) **2)

! senario two: random initial momentum
!      call pQCDwt(ityp(numPart), px(numPart), py(numPart), pz(numPart),
!     &     p0(numPart), fmass(numPart), pweight(numPart))

      enddo

      return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       read PHSD diffusion coefficients
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_PHSD_qhat
      implicit none
      include 'df_coms.f'

      integer i, j
      character*10 dummy_c
      double precision dummy_f
   
      !debug
      !open(unit=16, file='Duke-collision_diffusion.dat',status='old')
      !open(unit=16, file='PHSD_diffusion.dat', status='old')
      !open(unit=16, file='Catania_pQCD_diffusion.dat', status='old')
      !open(unit=16, file='Catania_QPM_diffusion.dat', status='old')
      open(unit=16, file='LBT_diffusion.dat', status='old')
      !open(unit=16, file='Nantes_diffusion.dat', status='old')
      !!!! remember to change the grid size as well!!!
      read(16,*) dummy_c, dummy_c, dummy_c, dummy_c, dummy_c
      do 3011 i=1, PHSD_nT
        do 3011 j=1, PHSD_nE
          if ((i.eq.1).and.(j.eq.1)) then
            read(16,*) gamma_TL, gamma_EL, PHSD_A(i,j), PHSD_BL(i,j),
     &                 PHSD_BT(i,j)
          else if ((i.eq.PHSD_nT) .and. (j.eq.PHSD_nE)) then
            read(16, *) gamma_TH, gamma_EH, PHSD_A(i,j), PHSD_BL(i,j),
     &                  PHSD_BT(i,j)
          else
            read(16,*) dummy_f, dummy_f, PHSD_A(i,j), PHSD_BL(i,j),
     &                 PHSD_BT(i,j)
          endif

 3011 continue

      gamma_dT = (gamma_TH - gamma_TL)/(PHSD_nT-1)
      gamma_dE = (gamma_EH - gamma_EL)/(PHSD_nE-1)

      !write(6,*) "read in tabled diffusion: ", gamma_TL, gamma_TH,
      !&           gamma_EL, gamma_EH, gamma_dT, gamma_dE
 
      close(16)
      return

      end subroutine



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        return PHSD qhat
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PHSD_qhat(drag,kappaL,
     &                              kappaT,temperature,momentum)
      implicit none
      include 'df_coms.f'

      
      integer i, j
      double precision delta_x, delta_y, resI, resJ
      double precision drag,kappaL,kappaT,temperature,momentum
      double precision mass, energ
    
      if (iflav.eq.4) then
        mass = cMass
      else if (iflav.eq.5) then
        mass = bMass
      endif

      !!debug
      !temperature = 1.038961*0.154
      !momentum = 1.207958


      resI = (temperature/0.154 - gamma_TL)/gamma_dT + 1
      resJ = (momentum - gamma_EL)/gamma_dE + 1
      i = int(resI)
      j = int(resJ)
      delta_x = resI - i
      delta_y = resJ - j

      if (resI.lt.1) then
        i=1
        delta_x=0
      endif

      if (resI.gt.PHSD_nT) then
        i=PHSD_nT
        delta_x=0
      endif

      if (resJ.lt.1) then
        j=1
        delta_y=0
      endif

      if (resJ.gt.PHSD_nE) then
        j=PHSD_nE
        delta_y=0
      endif

c interpolate gamma_table
      drag =  PHSD_A(i,j)*(1-delta_x)*(1-delta_y)
     &      + PHSD_A(i+1,j)*delta_x*(1-delta_y)
     &      + PHSD_A(i,j+1)*(1-delta_x)*delta_y
     &      + PHSD_A(i+1,j+1)*delta_x*delta_y

      kappaL = PHSD_BL(i,j)*(1-delta_x)*(1-delta_y)
     &       + PHSD_BL(i+1,j)*delta_x*(1-delta_y)
     &       + PHSD_BL(i,j+1)*(1-delta_x)*delta_y
     &       + PHSD_BL(i+1,j+1)*delta_x*delta_y

      kappaT = PHSD_BT(i,j)*(1-delta_x)*(1-delta_y)
     &       + PHSD_BT(i+1,j)*delta_x*(1-delta_y)
     &       + PHSD_BT(i,j+1)*(1-delta_x)*delta_y
     &       + PHSD_BT(i+1,j+1)*delta_x*delta_y
    
      !corrections!
      energ = sqrt(momentum**2 + mass**2)

      ! debug
      !write(6,*) i,j,delta_x, delta_y, 
      !&      temperature/0.154, momentum, drag, kappaL, kappaT

      drag = 2*temperature* drag*energ/momentum/inv_fm_to_GeV
      kappaL = kappaL/inv_fm_to_GeV
      kappaT = kappaT/inv_fm_to_GeV

      end subroutine


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           read PHSD qhat (and impose Einstein relation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_PHSD_qhat_ER(drag,kappaL,kappaT,T,p)
      implicit none
      include 'df_coms.f'

      integer i, j
      double precision delta_x, delta_y
      double precision drag, kappaL, kappaT, T, p, energ
      double precision kappaL_
      double precision mass

      if (iflav.eq.4) then
        mass = cMass
      else if (iflav.eq.5) then
        mass = bMass
      endif

      i = (T/0.154 - gamma_TL)/gamma_dT
      j = (p - gamma_EL)/gamma_dE
      delta_x = i - int(i)
      delta_y = j - int(j)
      i = int(i) + 1
      j = int(j) + 1

      if (i.lt.1) then
        i = 1
        delta_x = 0
      endif

      if (i.gt.PHSD_nT) then
        i = PHSD_nT
        delta_x = 0
      endif

      if (j.lt.1) then
        j=1
        delta_y=0
      endif

      if (j.gt.PHSD_nE) then
        j=PHSD_nE
        delta_y=0
      endif

cc interpolate gamma_table
      kappaL = PHSD_BL(i,j)*(1-delta_x)*(1-delta_y)
     &       + PHSD_BL(i+1,j)*delta_x*(1-delta_y)
     &       + PHSD_BL(i,j+1)*(1-delta_x)*delta_y
     &       + PHSD_BL(i+1,j+1)*delta_x*delta_y

      if (j.eq.1) then
        kappaL_ = kappaL
      else
        kappaL_ = PHSD_BL(i,j-1)*(1-delta_x)*(1-delta_y)
     &       + PHSD_BL(i+1,j-1)*delta_x*(1-delta_y)
     &       + PHSD_BL(i,j+1-1)*(1-delta_x)*delta_y
     &       + PHSD_BL(i+1,j+1-1)*delta_x*delta_y
      endif

      kappaT = PHSD_BT(i,j)*(1-delta_x)*(1-delta_y)
     &       + PHSD_BT(i+1,j)*delta_x*(1-delta_y)
     &       + PHSD_BT(i,j+1)*(1-delta_x)*delta_y
     &       + PHSD_BT(i+1,j+1)*delta_x*delta_y


      energ = sqrt(p**2 + mass**2)

cc for ito (pre-point):
cc     Gamma = BL/(2*E*T) - (BL-BT)/p^2 - d.BL/d.p^2
!       drag = kappaL/(2*T*energ) - (kappaL - kappaT)/p**2
!     &    -(kappaL - kappaL_)/(gamma_dE * 2*p)

cc for post-point (not a very good choice)
c!!     Gamma = BL/(2*T*E) - 1./p^2 (sqrt(BL) - sqrt(BT))^2
!      drag= kappaL/(2*T*energ)-1d0/p**2*(sqrt(kappaL)-sqrt(kappaT))**2 
    
c for isotropic case (determine kT and kL from drag)
      drag = kappaT/(2*T*energ)
      kappaL = kappaT
      


      drag = drag*energ * 2*T/inv_fm_to_GeV
      kappaL = kappaL/inv_fm_to_GeV
      kappaT = kappaT/inv_fm_to_GeV
      
      end subroutine
