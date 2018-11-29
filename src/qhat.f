      subroutine read_qhat_table

      implicit none
      include 'df_coms.f'

      integer i,j
      character*10 dummy_c
      integer dummy_n
      double precision dummy_f

      if(iflav.eq.4) then
         open(unit=15,file='gamma-table_charm.dat', status='old')
      elseif(iflav.eq.5) then
         open(unit=15,file='gamma-table_bottom.dat', status='old')
      else
         write(*,*) "un-recognized heavy quark ID ..."
         stop
      endif

      read(15,*) dummy_c, dummy_c, dummy_c
      do 2096 j=1, gamma_nE
        do 2096 i=1, gamma_nT
            if ((i.eq.1) .and. (j.eq.1)) then
                read(15,*) gamma_TL, gamma_EL,qhat_over_T3(i,j)
            else if ((i.eq. gamma_nT) .and. (j.eq.gamma_nE)) then
                read(15,*) gamma_TH, gamma_EH, qhat_over_T3(i,j)
            else
                read(15,*, err=2098,end=2097)
     &      dummy_f, dummy_f, qhat_over_T3(i,j)
            endif
    
 2096 continue
        
      gamma_dT=(gamma_TH-gamma_TL)/(gamma_nT-1)
      gamma_dE=(gamma_EH-gamma_EL)/(gamma_nE-1)

      close(15)
      return

 2097 continue
      write(6,*) 'ERROR: EOF reached in qhat table'
      write(6,*) 'terminating ...'
      stop
      return

 2098 continue
      write(6,*) 'READ-ERROR in qhat table'
      write(6,*) 'terminating ...'
      stop
      return

      end subroutine


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_qhat(qhat_value, temperature, energy)
c a subroutine that reads in qhat_over_T3 table and return qhat value
      implicit none
      include 'df_coms.f'

      integer i, j
      double precision delta_x, delta_y
      double precision qhat_value, temperature, energy
     
      i = (temperature - gamma_TL)/gamma_dT 
      j = (energy - gamma_EL)/gamma_dE
      delta_x = i - int(i)
      delta_y = j - int(j)
      i = int(i) + 1
      j = int(j) + 1


      if (i.lt.1) then
        i=1
        delta_x=0
      endif

      if (i.gt.gamma_nT) then
        i=gamma_nT
        delta_x=0
      endif

      if (j.lt.1) then
        j=1
        delta_y=0
      endif

      if (j.gt.gamma_nE) then 
        j=gamma_nE
        delta_y=0
      endif

c interpolate gamma_table
      qhat_value = qhat_over_T3(i,j)*(1-delta_x)*(1-delta_y) 
     &        + qhat_over_T3(i+1,j)*delta_x*(1-delta_y)
     &        + qhat_over_T3(i,j+1)*(1-delta_x)*delta_y
     &        + qhat_over_T3(i+1,j+1)*delta_x*delta_y
      end subroutine




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   read all the diffusion coefficients (anisotropic version
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine read_kappaLT_table
       implicit none
       include 'df_coms.f'

       integer i, j
       character*10 dummy_c
       double precision dummy_f

       open(unit=16, file='gamma-table_charm.dat', status='old')
       read(16, *) dummy_c, dummy_c, dummy_c, dummy_c, dummy_c
     
       do 3011 j=1, gamma_nE
         do 3011 i=1, gamma_nT
           if ((i.eq.1) .and. (j.eq.1)) then
             read(16,*) gamma_TL, gamma_EL, kappaA(i,j), kappaBL(i,j)
     &                          , kappaBT(i,j)
           else if ((i.eq.gamma_nT) .and. (j.eq.gamma_nE)) then
             read(16,*) gamma_TH, gamma_EH, kappaA(i,j), kappaBL(i,j)
     &                          , kappaBT(i,j)
           else
             read(16,*,err=3098,end=3097) dummy_f,dummy_f,kappaA(i,j)
     &                          , kappaBL(i,j), kappaBT(i,j)
           endif

        
 3011 continue
    
       gamma_dT = (gamma_TH - gamma_TL)/(gamma_nT -1)
       gamma_dE = (gamma_EH - gamma_EL)/(gamma_nE -1)
       close(16)
       return 

 3097 continue
       write(6,*)  "ERROR: EOF reached in kappa table"
       write(6,*)  "teminating ...."
       return 

 3098 continue
       write(6,*)   "READ-ERROR in kappa table"
       write(6,*)   "teminating ...."
       stop
       return 

       end subroutine

 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_kappaLT(drag,kappaL,kappaT,temperature,energy)
c a subroutine that reads in qhat_over_T3 table and return qhat value
      implicit none
      include 'df_coms.f'

      integer i, j
      double precision delta_x, delta_y, resI, resJ
      double precision drag,kappaL,kappaT, temperature, energy
      double precision momentum, mass
    
      if (iflav .eq. 4) then
        mass = cMass
      else if (iflav .eq. 5) then
        mass = bMass
      endif

      resI = (temperature - gamma_TL)/gamma_dT + 1
      resJ = (energy - gamma_EL)/gamma_dE + 1
      i = int(resI)
      j = int(resJ)
      delta_x = resI - i
      delta_y = resJ - j


      if (resI.lt.1) then
        i=1
        delta_x=0
      endif

      if (resI.gt.gamma_nT) then
        i=gamma_nT
        delta_x=0
      endif

      if (resJ.lt.1) then
        j=1
        delta_y=0
      endif

      if (resJ.gt.gamma_nE) then 
        j=gamma_nE
        delta_y=0
      endif

c interpolate gamma_table
      drag = kappaA(i,j)*(1-delta_x)*(1-delta_y) 
     &        + kappaA(i+1,j)*delta_x*(1-delta_y)
     &        + kappaA(i,j+1)*(1-delta_x)*delta_y
     &        + kappaA(i+1,j+1)*delta_x*delta_y

      kappaL = kappaBL(i,j)*(1-delta_x)*(1-delta_y) 
     &        + kappaBL(i+1,j)*delta_x*(1-delta_y)
     &        + kappaBL(i,j+1)*(1-delta_x)*delta_y
     &        + kappaBL(i+1,j+1)*delta_x*delta_y

      kappaT = kappaBT(i,j)*(1-delta_x)*(1-delta_y) 
     &        + kappaBT(i+1,j)*delta_x*(1-delta_y)
     &        + kappaBT(i,j+1)*(1-delta_x)*delta_y
     &        + kappaBT(i+1,j+1)*delta_x*delta_y

    
      if (kappaL .lt. 0) then
        write(6,*) kappaL,kappaBL(i,j), kappaBL(i+1,j), kappaBL(i,j+1),
     &          kappaBL(i+1,j+1), delta_x, delta_y
      endif

      ! convert the unit
      momentum = sqrt(energy**2 - mass**2)
      drag = 2*temperature*drag*energy/momentum/inv_fm_to_GeV
      kappaL = kappaL/inv_fm_to_GeV
      kappaT = kappaT/inv_fm_to_GeV


      end subroutine


          
