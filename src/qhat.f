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


