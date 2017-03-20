      subroutine read_qhat_table

      implicit none
      include 'df_coms.f'

      integer i,j
      character*10 dummy_c
      integer dummy_n
      double precision dummy_f
      
      open(unit=10,file='gamma-table.dat',
     &        status='old',form='formatted')

      i=7
      do while(i.le.(gamma_nT+6))
         j=-2
         do while(j.le.gamma_np)
            read(unit=10,fmt=101,err=2098,end=2097)
     &      dummy_n,dummy_n,dummy_f,dummy_f,dummy_f,qhat_over_T3(i,j)
c            write(6,*) i,j,qhat_over_T3(i,j)
            j=j+1
         enddo
         qhat_over_T3(i,gamma_np+1)=0d0
         i=i+1
      enddo

 101  format(I4,2X,I4,4(2X,D12.6))

      close(10)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_qhat(qhat_value,temperature,momentum)

      implicit none
      include 'df_coms.f'

      integer i,j
      double precision delta_x,delta_y
      double precision qhat_value,temperature,momentum
      double precision qhat_j1,qhat_j2

      if(log(momentum).lt.0d0) then
         j=int(log(momentum))-1
      else 
         j=int(log(momentum))
      endif

      i=int(temperature/delta_Te)

      if(i.lt.7.or.i.gt.gamma_nT+6) then
         write(6,*) "Error: T out of range for qhat!"
         write(6,*) "Terminating ..."
         stop
      endif

      if(j.lt.-2) j=-2
      if(j.gt.7) j=7

c interpolate in the x (temperature) direction

      delta_x=temperature-delta_Te*i
      qhat_j1=qhat_over_T3(i,j)+delta_x/delta_Te*
     &        (qhat_over_T3(i+1,j)-qhat_over_T3(i,j))
      qhat_j2=qhat_over_T3(i,j+1)+delta_x/delta_Te*
     &        (qhat_over_T3(i+1,j+1)-qhat_over_T3(i,j+1))

c interpolate in the y (momentum) direction

      delta_y=log(momentum)-j
      if(delta_y.lt.0d0) delta_y=0d0

      if(j.eq.7) then
         qhat_value=qhat_j1
      else
         qhat_value=qhat_j1+delta_y*(qhat_j2-qhat_j1)
      endif

      qhat_value=qhat_value*temperature**3

      end subroutine

           

