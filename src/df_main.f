cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     df_main.f version 1.0
c
c     Steffen A. Bass 9/2005
c
c     This program reads in the (Lagrangian) 3D Hydro evolution
c     and a set of (test)particles (OSCAR format), which then 
c     propate through the hydro-medium according to a Langevin Equation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program df_main

      implicit none
      include 'ucoms.f'
      include 'df_coms.f'


      DATA MRLU/19780503,0,0,97,33,0/

      integer nevent,iret,i,pevent,j,rand_seed,wt_table,sampleCount
      integer oflag,hflag,ioflag,geoflag
      double precision deltat,energ
      double precision rlu,p_length,utheta,phi

      rand_seed=-1
      call sseed(rand_seed)

      write(6,*) '*****************************************************'
      write(6,*) '* HYDRO+PCM+LANGEVIN 1.0 - Asakawa - Bass - Mueller *'
      write(6,*) '*****************************************************'


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call parameterRead
      
c     "static" can be: 
c                       1 -- Static Medium;
c                       2 -- OSU hydro;
c                       3 -- PHSD Medium

      if(static.ne.1.and.static.ne.2.and.static.ne.3) then
         write(6,*) "Inappropriate choice for static ..."
         write(6,*) "Terminating ..."
         stop
      endif

c     "flag_rad" can be:  0 -- Only collsional energy loss;
c                         1 -- Only radiative energy loss;
c                         2 -- Both collsional and radiative.

      if(flag_rad.ne.0.and.flag_rad.ne.1.and.flag_rad.ne.2) then
         write(6,*) "Inappropriate choice for flag_rad ..."
         write(6,*) "Terminating ..."
         stop
      endif

      if(reweight.ne.0.and.reweight.ne.1.and.reweight.ne.2) then
         write(6,*) "Inappropriate choice for reweight ..."
         write(6,*) "Terminating ..."
         stop
      endif

c     open hydro file  ccccccccccccccccccccccccccccccccccccccccccccc
c     read in OSU hydro
      if(static.eq.2.and.out_skip.ne.0) then
         call readHydroFiles_initialEZ("JetData.h5")
      endif
    
      if(static.eq.3.and.out_skip.ne.0) then
         call readPHSD_MediumFile()
      endif


c     read weights of pT distributions of initial heavy quark
      if(reweight.eq.1) then
         wt_table=0
         call read_wt_table(wt_table)
         if(wt_table.ne.0) then
            write(6,*) "error occurs during reading weight table ..."
            write(6,*) "terminating ..."
            stop
         endif
      endif

c     read table for qhat if it depends on T and p
      if(qhat_TP.eq.1 .or. qhat_TP.eq.5) call read_qhat_table
      if(qhat_TP.eq.6) call read_PHSD_qhat

c     read table for gluon radiation
      if(flag_rad.ne.0) call read_radTable
      sampleCount=NUMSAMP
 

c     read in initial particles
      if(HQ_input.eq.3) then ! read from an x-y position list
         if(reweight.eq.0) then
            write(*,*) "Please set reweight to 1 or 2"
            stop
         endif
         call read_xy_init(iret)
      elseif(HQ_input.eq.4) then ! read from PHSD initial condition
         call read_PHSD_init(iret)
      elseif(HQ_input.eq.0) then ! read in static initial condition
         call read_static_init(iret)
      else
         write(*,*) "Wrong value for HQ_input"
         stop
      endif

      if(iret.eq.0) then
         write(6,*) '*** end of particle event-input ***'
         write(6,*) ' dEsum = ',dEsum/des_cntr
         stop
      elseif(iret.lt.0) then
         write(6,*) '**** error during read in initial HQ ***'
         write(6,*) '**** last event was skipped, terminating...***'
         stop
      endif

      if(npart.gt.nmax) then
         write(6,*)'warning: too many particles in event ',nevent
         write(6,*)'terminating...'
         stop
      endif

      if(ref_frame.eq.1) then
         reffram = 'CMS'
      elseif(ref_frame.eq.2) then
         reffram = 'LRF'
      else
         write(6,*) 'Improper choice of reference frame: ref_frame = ',
     .               ref_frame
         write(6,*) 'terminating ...'
         stop
      endif

c     loop over samples cccccccccccccccccccccccccccccccccccccccccccccccc

      pevent=0

 3    continue

      dEsum=0d0
      des_cntr=0
      flag_stop=0
      pevent=pevent+1

      write(6,*)
      write(6,*)
      write(6,*) sampleCount," SAMPLES LEFT ..."
      write(6,*)'-> reading PCM event # ',pevent

c     initialize particles' positions
      if(HQ_input.eq.3.and.sampleCount.ne.NUMSAMP) then
         call reSampleXY3
      endif


cccccccc now filter out particles used for calculation cccccccccccccccc
      npt=0
c debug
      write(6,*) '->  selecting particles'
      write(6,*) npart
      do 10 i=1,npart
c select flavor and tune particle phase space distribution for different
c purposes
         if(iabs(ityp(i)).eq.iflav) then
            npt=npt+1   
            do 11 j=1,evsamp
               pid(npt,j)=ityp(i)
               p_px(npt,j)=px(i)
               p_py(npt,j)=py(i)
               p_pz(npt,j)=pz(i)
               p_mass(npt,j)=fmass(i)
               p_p0(npt,j)=sqrt(px(i)**2+py(i)**2+pz(i)**2+fmass(i)**2)
               p_r0(npt,j)=r0(i)
               p_rx(npt,j)=rx(i)
               p_ry(npt,j)=ry(i)
               p_rz(npt,j)=rz(i)
               if(reweight.ne.0) p_wt(npt,j)=pweight(i)

c sample initial charms as Moore and Teaney did (with my weight method)
               if (Moore.eq.1) call MTCharm(ityp(i),p_px(npt,j),
     &              p_py(npt,j),p_pz(npt,j),p_p0(npt,j),p_mass(npt,j),
     &              p_wt(npt,j))

c sample initial HQ pT distribution with input weight table
               if (HQ_input.eq.1.and.reweight.ne.0) 
     &              call pQCDwt(ityp(i),p_px(npt,j),
     &              p_py(npt,j),p_pz(npt,j),p_p0(npt,j),p_mass(npt,j),
     &              p_wt(npt,j))
               
c special initialization for calculating correlation
               if (HQ_input.eq.1.and.corr_flag.eq.1.and.mod(npt,2).eq.0)
     &            then
                  pid(npt,j)=-pid(npt-1,j)
                  p_px(npt,j)=-p_px(npt-1,j)
                  p_py(npt,j)=-p_py(npt-1,j)
                  p_pz(npt,j)=-p_pz(npt-1,j)
                  p_p0(npt,j)=p_p0(npt-1,j)
                  p_mass(npt,j)=p_mass(npt-1,j)
                  p_rx(npt,j)=p_rx(npt-1,j)
                  p_ry(npt,j)=p_ry(npt-1,j)
                  p_rz(npt,j)=p_rz(npt-1,j)
                  p_r0(npt,j)=p_r0(npt-1,j)
                  p_wt(npt,j)=p_wt(npt-1,j)
               endif

               p_ipx(npt,j)=p_px(npt,j)
               p_ipy(npt,j)=p_py(npt,j)
               p_ipz(npt,j)=p_pz(npt,j)
               p_ip0(npt,j)=p_p0(npt,j)
               p_ipT(npt,j)=sqrt(p_px(npt,j)**2+p_py(npt,j)**2)

c initialize cell velocity of hydro medium
               c_vx(npt,j)=0d0
               c_vy(npt,j)=0d0
               c_vz(npt,j)=0d0

 11         continue
         endif
 10   continue ! end re-distribute particles in the phase space

c now synchronize particles to 1st time-step

      if(static.eq.2) initt = 0.6d0 ! OSU hydro starts at 0.6~fm/c

      do 20 i=1,npt
         do 21 j=1,evsamp
            deltat=initt-p_r0(i,j)
c note: back-propagation is possible here for the particle
c       freeze-out time being later than the hydro initial
c       time! Model works best if there is a clear separation
c       of time-scales between PCM and Hydro...

            energ = p_p0(i,j)
            p_r0(i,j)  = initt
            p_rx(i,j)  = p_rx(i,j) + p_px(i,j)/energ*deltat
            p_ry(i,j)  = p_ry(i,j) + p_py(i,j)/energ*deltat
c            p_rz(i,j)  = p_rz(i,j) + p_pz(i,j)/energ*deltat
            p_rz(i,j)  = 0d0
            if(abs(p_rz(i,j)).gt.initt) then
               p_rz(i,j)=sign(initt-1d-10,p_rz(i,j))
            endif
            p_reta(i,j)=0.5d0*log((p_r0(i,j)+p_rz(i,j))/
     &                            (p_r0(i,j)-p_rz(i,j)))

            time_lim(i,j)=initt   ! record time of last interaction
 21      continue
 20   continue

cdebug
      write(6,*) '-> ', npt,' particles initialized for calculation'
      write(6,*) '-> ', evsamp,' events will be run in parallel'

      if(out_skip.eq.0) then ! write out init. distri. and end program
         call output(pevent)
         write(6,*) "Write out initial distribution only."
         stop
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for gluon radiation, initialization

      if(flag_rad.ne.0) then

         delta_HQener=HQener_max/HQener_gn
         delta_tg=t_max/t_gn
         delta_temp=(temp_max-temp_min)/temp_gn
         cprob_gt1=0
         ctemp_gtmax=0
         ct_gtmax=0
         cHQE_gtmax=0
         num_gluon=0
         count_Ecut=0
         Etot_gluon=0d0

         OPEN(30,FILE='radiation.log',ACCESS='APPEND',STATUS='UNKNOWN')

         time_tg=0.d0
         i=1
         do while(i.lt.temp_gn+1)
            j=1
            do while(j.lt.HQener_gn+1)
               dNg_over_dt(1,i,j)=0.d0
               max_dNgfnc(1,i,j)=0.d0
               j=j+1
            enddo
            i=i+1
         enddo

         
c initialize fgluon(i,j): whether the (i,j) HQ has emitted a gluon
         do 40 i=1,npt
            do 41 j=1,evsamp
               fgluon(i,j)=0
               init_int(i,j)=0
               Tinteval_lrf(i,j)=0d0
               t_init(i,j)=0.6d0
               previous_kT(i,j)=1000d0
c variables related to the last radiated gluon -- for tau cut purpose
               lastG_tau(i,j)=0d0
               lastG_px(i,j)=0d0
               lastG_py(i,j)=0d0
               lastG_pz(i,j)=0d0
 41         continue
 40      continue
         
         write(6,*) "Intialization of radiation succeeds :)"

      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(static.eq.1.or.static.eq.2) then
         ntsteps=tsteps_cut  ! total time steps for Langevin evolution
         tau=0.6d0-0.1d0
      endif

      do 22 tstep=1,ntsteps        
         tau=tau+0.1d0
         
c do physics for each time step -- Langevin evolution
         do 32 lstep=1,nlang
            if(flag_stop.ge.npt*evsamp) then
               write(6,*) "Diffusion ends at time step: ", 
     &                    (tstep-1)*nlang+lstep-1
               goto 23
            endif

            call exec_dif
 32      continue

c output in time-step loop
c insert output condition
      if(out_skip.ne.1.and.mod(tstep,out_skip).eq.0) then
c      or check an exact time-step
c      if(tstep.eq.174) then

c  transform to local rest frame if necessary (determine T from p)
         if((static.ne.1).and.(ref_frame.eq.2)) then
           do 34 i=1,npt
             do 35 j=1,evsamp
               call rotbos(0d0,0d0,-c_vx(i,j),-c_vy(i,j),
     &           -c_vz(i,j),p_px(i,j),p_py(i,j),p_pz(i,j),
     &           p_p0(i,j)) 
 35          continue
 34        continue
         endif

         call output(pevent)

c  transform back to the lab frame
         if((static.ne.1).and.(ref_frame.eq.2)) then
           do 36 i=1,npt
             do 37 j=1,evsamp
               call rotbos(0d0,0d0,c_vx(i,j),c_vy(i,j),c_vz(i,j),
     &              p_px(i,j),p_py(i,j),p_pz(i,j),p_p0(i,j)) 
 37          continue
 36        continue
         endif
            
      endif
c     close timeteps cccccccccccccccccccccccccccccccccccccccccccccccc

 22   continue
 23   continue


c free-stream particles back to their space-time of last interaction
      do 24 i=1,npt
         do 25 j=1,evsamp
            deltat = p_r0(i,j)-time_lim(i,j)
            energ = p_p0(i,j)
            p_r0(i,j) = time_lim(i,j)
            p_rx(i,j) = p_rx(i,j) - p_px(i,j)/energ*deltat
            p_ry(i,j) = p_ry(i,j) - p_py(i,j)/energ*deltat
            p_rz(i,j) = p_rz(i,j) - p_pz(i,j)/energ*deltat
            p_reta(i,j) = 0.5d0*log((p_r0(i,j)+p_rz(i,j))/
     &                    (p_r0(i,j)-p_rz(i,j)))
 25      continue
 24   continue

c check whether the last gluon is able to form, if not, add its 
c momentum back to the heavy quark (only a rough manipulation)
      do 26 i=1,npt
         do 27 j=1,evsamp
            if(Tinteval_lrf(i,j).lt.lastG_tau(i,j)) then
               p_px(i,j) = p_px(i,j) + lastG_px(i,j)
               p_py(i,j) = p_py(i,j) + lastG_py(i,j)
               p_pz(i,j) = p_pz(i,j) + lastG_pz(i,j)
               p_p0(i,j) = sqrt(p_px(i,j)**2+p_py(i,j)**2+
     &                     p_pz(i,j)**2+p_mass(i,j)**2)
               energ = p_p0(i,j)
            endif
 27      continue
 26   continue

c if out_skip is 1, output the last time-step only
      if(out_skip.eq.1) call output(pevent)
      
      if(flag_rad.ne.0) then
         write(30,*) "times for probability to exceed 1:",cprob_gt1
         write(30,*) "times for time to exceed t_max:",ct_gtmax
         write(30,*) "times for temperature to exceed temp_max:",
     &               ctemp_gtmax
         write(30,*) "times for HQenergy to exceed HQener_max:",
     &               cHQE_gtmax
         write(30,*) "total number of gluons radiated:",num_gluon
         write(30,*) "average gluon number:",1d0*num_gluon/npt/evsamp
         write(30,*) "number of cutted gluons: ",count_Ecut
         write(30,*) "average gluon energy: ",Etot_gluon/num_gluon
         close(30)
      endif

      sampleCount = sampleCount-1
      if(sampleCount.gt.0) goto 3

      write(6,*) "PROGRAM ENDS SUCCESSFULLY :)"

      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SUBROUTINES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c most important subroutine -- do Langevin evolution

      subroutine exec_dif

      implicit none
      include 'df_coms.f'
      include 'ucoms.f'

      double precision tau_p,deltat,energ,betax,betay,betaz,tmp
      integer ix,iy,iz,i,j,k,iflag,count_sample
      double precision rlu,random(2),theta_gluon,dNg_over_dxdydt
      double precision kpx_gluon,kpy_gluon,kperp_gluon,klong_gluon
      double precision dkpx_gluon,dkpy_gluon,width_gluon,tau_f
      double precision dkT_gluon,plength
c      double precision kT_dum
c      double precision rgau_mean,rgau_width(3),rgau,kappa_d,kappa_t,
c     &                 kappa_l
      double precision rgau_mean, rgau_width(3), rgau
      double precision xi_gauss(3),T,v_x,v_y,v_z,v2,e_old,xi(3)
      double precision theta,phi,rr,dum
      double precision lim_low,lim_high,lim_int,max_Ng,deltat_lrf
c dummy for OSU hydro at the moment
      double precision edensity0,sdensity0
      integer ctl_OSU
      double precision KTfactor,KPfactor,KPTfactor

      double precision dum_D2piT, p0_lab
      rgau_mean=0d0
      flag_stop=0

c     set deltat
      deltat = (tmax-initt)/(ntsteps*nlang)
      if(static.eq.1.or.static.eq.2) deltat = 0.1d0/nlang
      tau_p=tau+deltat*(lstep-1)

      if(static.eq.1) then
         T=T_static
         initt=0.6d0
         if(static_cool.eq.1) then
            T=T_static*(tau_p/initt)**(-1d0/3d0)
         endif
         !write(6,*) '    infinite matter, temperature: ',T
      endif 

      do 10 i=1,npt
         do 11 j=1,evsamp
            iflag=1
c use OSU hydro
            if(static.eq.2) then

c This is the key call that reads hydro info (e,s,T,vx,vy) at a given (tau,x,y)-tuple
! use hdf5 file instead
                call readHydroInfoYingru(p_r0(i,j),p_rx(i,j),p_ry(i,j),
     &    p_rz(i,j),edensity0,sdensity0,T,betax,betay,betaz,ctl_OSU)

               if(temp_cut.ne.1) then
                  write(6,*) "There must be T cut for OSU hydro!"
                  write(6,*) "Terminating ..."
                  stop
               endif
               
               if(ctl_OSU.ne.0.or.T.lt.Tcut_critical) then
                  iflag=0
                  flag_stop=flag_stop+1
               endif
            endif ! OSU hydro


c do particle diffusion/propagation here...
            if(static.eq.1.or.iflag.ne.0) then

c update the time of the last interaction
              time_lim(i,j) = p_r0(i,j)+deltat 
              if(static.ne.1) then
                p0_lab = p_p0(i,j)
                call rotbos(0d0,0d0,-betax,-betay,-betaz,
     &          p_px(i,j),p_py(i,j),p_pz(i,j),p_p0(i,j))
! change the deltat in cell frame to lrf
                deltat_lrf = p_p0(i,j)/p0_lab * deltat
              else
                deltat_lrf=deltat
              endif
               
              energ = sqrt(p_px(i,j)**2+p_py(i,j)**2+p_pz(i,j)**2
     &                      +p_mass(i,j)**2)
              e_old=energ

              if(init_int(i,j).eq.0) then
                init_int(i,j)=1           ! start first interaction
                Tinteval_lrf(i,j)=0d0     ! for gluon radiation
              endif

c set transport coefficient if necessary
            if(qhat_TP.eq.0) then
              qhat = 4d0*alpha*T**3/C_F  ! a constant D2piT
            else if (qhat_TP.eq.5) then
              call get_qhat(qhat,T,energ)  !qhat_over_T3
              D2piT=25.1327d0/qhat
              dum_D2piT=qhatMin*(1d0+qhatSlope*(T/Tcut_critical-1d0))
              plength=sqrt(p_px(i,j)**2+p_py(i,j)**2+p_pz(i,j)**2)
!param18
!D2piT=(qhatPower**2*pT)**2/(1+(qhatPower**2*pT)**2) * pQCD +
!1/(1+(qhatPower**2*pT)**2) * linear
              plength=sqrt(p_px(i,j)**2+p_py(i,j)**2+ p_pz(i,j)**2)
              D2piT=(1d0-1d0/(1d0+((qhatPower**2)*plength)**2))*D2piT 
     &        + 1d0/(1d0+((qhatPower**2)*plength)**2)*dum_D2piT

              if (D2piT .lt. 0.1d0) then
                D2piT = 0.1d0
              endif

              alpha=6.2832d0/D2piT
              qhat=4d0*alpha/C_F*T**3
            else if (qhat_TP.eq.6) then
              plength=sqrt(p_px(i,j)**2+p_py(i,j)**2+p_pz(i,j)**2)
              call get_PHSD_qhat(kappa_d,kappa_l,kappa_t,T,plength)
              D2piT=25.13727d0/(2*kappa_t/ T**3)
              alpha=6.2832d0/D2piT
              qhat=4d0*alpha/C_F*T**3

!              !write(6,*) "read in PHSD diffusion coefficient: " ,T
              write(6,*) T
     &                   ,plength,kappa_d,kappa_l,kappa_t 
            else
              write(6,*) "Wrong value for qhat_TP."
              write(6,*) "Terminating ..."
            endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for gluon radiation

c whether multiple gluon radiation is allowed
               
            if(flag_rad.ne.0) then
              HQenergy=energ
              HQmass=p_mass(i,j)
              temp_med=T
              Tinteval_lrf(i,j)=Tinteval_lrf(i,j)+deltat_lrf
              time_gluon=Tinteval_lrf(i,j)
c if tau_p or T is larger than the corresponding maximum value in 
c the dNg_over_dt table, report error
c                  if(time-time_init.gt.t_max) then
c                     write(30,*) 'time exceeds t_max'
c                     write(30,*) time,time_init,temp_med,HQenergy
c                     ct_gtmax=ct_gtmax+1
c                     time=t_max+time_init
c                  endif

              if(time_gluon.gt.t_max) then
                write(30,*) 'accumulated time exceeds t_max'
                write(30,*) tau_p,time_gluon,temp_med,HQenergy
                ct_gtmax=ct_gtmax+1
                time_gluon=t_max
              endif

              if(temp_med.gt.temp_max) then
                write(30,*) 'temperature exceeds temp_max'
                write(30,*) tau_p,temp_med,HQenergy
                ctemp_gtmax=ctemp_gtmax+1
                temp_med=temp_max
              endif

              if(HQenergy.gt.HQener_max) then
                write(30,*) 'HQenergy exceeds HQener_max'
                write(30,*) tau_p,temp_med,HQenergy
                cHQE_gtmax=cHQE_gtmax+1
                HQenergy=HQener_max
              endif

              if(temp_med.lt.temp_min) then
                write(6,*) 'temperature drops below temp_min'
                write(6,*) 'terminating ...'
                stop
              endif

              time_num=int(time_gluon/delta_tg+0.5d0)+1
              temp_num=int((temp_med-temp_min)/delta_temp+0.5d0)
              HQenergy_num=int(HQenergy/delta_HQener+0.5d0)
              delta_Ng=dNg_over_dt(time_num,temp_num,HQenergy_num)
     &                *6d0/D2piT
              max_Ng=max_dNgfnc(time_num,temp_num,HQenergy_num)
     &                *6d0/D2piT

!             if(qhat_TP.eq.1) delta_Ng=delta_Ng*KTfactor
! by yingru, second parameterization: (let us not to rescale alpha_s this time)

              if(delta_Ng*deltat_lrf.gt.1.d0) then
                write(30,*) 'gluon emission probability exceeds 1'
                write(30,*) tau_p,time_gluon,temp_med,HQenergy,
     &                      delta_Ng*deltat_lrf
                cprob_gt1=cprob_gt1+1
              endif

c check whether gluon emission is allowed for HQ at this energy
              lim_low = PI*temp_med/HQenergy
              lim_high = 1.d0
              lim_int = lim_high-lim_low

              if(lim_int.gt.EPS.and.
     &           rlu(0).lt.(delta_Ng*deltat_lrf-EPS).and.
     &           2d0*HQenergy*(HQenergy-PI*temp_med).gt.HQmass**2
     &           +EPS) then
c a gluon might be emitted
c random(1) corresponds to x, and random(2) corresponds to y
                random(1)=lim_low+lim_int*rlu(0)
                random(2)=rlu(0)
                do while(tau_f(random(1),random(2)).lt.1d0/PI/temp_med)
                  random(1)=lim_low+lim_int*rlu(0)
                  random(2)=rlu(0)
                enddo
 
                count_sample=0
                do while(max_Ng*rlu(0).gt.
     &                  dNg_over_dxdydt(random(1),random(2)))
                  count_sample=count_sample+1
c debug
                  if(count_sample.gt.1e+5) then
                    write(6,*) "give up loop at point 1 ..."
c                   write(6,*) tau_p,time_gluon,temp_med,HQenergy,delta_Ng
                    kpx_gluon=0d0
                    kpy_gluon=0d0
                    klong_gluon=0d0
                    goto 3312
                  endif
 
                  random(1)=lim_low+lim_int*rlu(0)
                  random(2)=rlu(0)
                  do while(tau_f(random(1),random(2)).lt.
     &              1d0/PI/temp_med)   
                    random(1)=lim_low+lim_int*rlu(0)
                    random(2)=rlu(0)
                  enddo
 
                enddo

                theta_gluon=2d0*PI*rlu(0)
                kperp_gluon=random(1)*random(2)*HQenergy
                kpx_gluon=kperp_gluon*cos(theta_gluon)
                kpy_gluon=kperp_gluon*sin(theta_gluon)
                klong_gluon=random(1)*HQenergy*
     &               sqrt(1d0-random(2)*random(2))

c throw away collinear gluon or gluon that disobays kT ordering
c             if(kperp_gluon.gt.previous_kT(i,j)) then
c               if(kperp_gluon.lt.sqrt(2d0/3d0)*PI*temp_med) then
c                 kpx_gluon=0d0
c                 kpy_gluon=0d0
c                 klong_gluon=0d0
c                 count_Ecut=count_Ecut+1
c                 goto 3312
c               endif

c               kT_dum = kperp_gluon

c for checking purpose, four ways to remove the momentum broadening of gluons

                if(narrow.eq.1) then
                  kpx_gluon=0.d0
                  kpy_gluon=0.d0
                endif

                if(narrow.eq.3) then
                  width_gluon=sqrt(N_c*qhat*
     &                  tau_f(random(1),random(2))/2.d0)
 3310 continue
                  dkpx_gluon=rgau(0,0.d0,width_gluon)
                  dkpy_gluon=rgau(0,0.d0,width_gluon)
                  dkT_gluon=sqrt(dkpx_gluon**2+dkpy_gluon**2)

                  if(dkT_gluon.gt.kperp_gluon) goto 3310

                  kperp_gluon=sqrt(kperp_gluon**2-dkT_gluon**2)
                  theta_gluon=2*PI*rlu(0)
                  kpx_gluon=kperp_gluon*cos(theta_gluon)
                  kpy_gluon=kperp_gluon*sin(theta_gluon)     
                endif

                if(narrow.eq.4) then
                  count_sample=0
                  width_gluon=sqrt(N_c*qhat*
     &              tau_f(random(1),random(2))/2.d0)
 3311 continue
                  dkpx_gluon=rgau(0,0.d0,width_gluon)
                  dkpy_gluon=rgau(0,0.d0,width_gluon)
                  dkT_gluon=sqrt(dkpx_gluon**2+dkpy_gluon**2)
                  if(dkT_gluon.gt.kperp_gluon.or.
     &               dkpx_gluon**2+dkpy_gluon**2
     &               -2*dkpx_gluon*kpx_gluon-2*dkpy_gluon*kpy_gluon
     &               .gt.0.d0) then
                    count_sample=count_sample+1
c debug
                  if(count_sample.gt.1e+6) then
                     write(6,*) "give up loop at point2: " 
                     write(6,*) kpx_gluon,kpy_gluon,width_gluon
                     goto 3313
                   endif

                   goto 3311
                 endif

                 kpx_gluon=kpx_gluon-dkpx_gluon
                 kpy_gluon=kpy_gluon-dkpy_gluon
               endif

 3313 continue
               call getang(p_px(i,j),p_py(i,j),p_pz(i,j),
     &                    theta,phi,rr)
c rotate gluon momentum:
               call rotbos(theta,phi,0d0,0d0,0d0,
     &              kpx_gluon,kpy_gluon,klong_gluon,dum)

c gluon energy cannot exceed HQ kinetic energy
               if(kpx_gluon**2+kpy_gluon**2+klong_gluon**2.gt.
     &            (HQenergy-HQmass)**2) then
                 kpx_gluon=0.d0
                 kpy_gluon=0.d0
                 klong_gluon=0.d0
                 count_Ecut=count_Ecut+1
               else
c this gluon is valid
                 fgluon(i,j)=fgluon(i,j)+1
                 num_gluon=num_gluon+1
c                t_init(i,j)=tau_p
                 Tinteval_lrf(i,j)=0d0
c                previous_kT(i,j)=kT_dum
                 Etot_gluon=Etot_gluon+random(1)*HQenergy
                 lastG_px(i,j)=kpx_gluon
                 lastG_py(i,j)=kpy_gluon
                 lastG_pz(i,j)=klong_gluon
                 lastG_tau(i,j)=tau_f(random(1),random(2))*inv_fm_to_GeV
               endif

             else
c no gluon emission in this time step 
               kpx_gluon=0.d0
               kpy_gluon=0.d0
               klong_gluon=0.d0
             endif

           else
             kpx_gluon=0.d0
             kpy_gluon=0.d0
             klong_gluon=0.d0
           endif

 3312 continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c the following part is for collisional energy loss

           if(flag_rad.ne.1) then
             v_x = p_px(i,j)/energ 
             v_y = p_py(i,j)/energ 
             v_z = p_pz(i,j)/energ

             v2 = v_x**2 + v_y**2 + v_z**2

c debug
             v2 = min(0.9999999,v2)

c     get paramters for noise calculation:
             if(qhat_TP.ne.6) then
               kappa_d=2*alpha*T**3/inv_fm_to_GeV ! alpha is a dimensionless factor
               kappa_l=2*alpha*T**3/inv_fm_to_GeV ! alpha is a dimensionless factor
               kappa_t=2*alpha*T**3/inv_fm_to_GeV ! alpha is a dimensionless factor
             endif
               
             if(pdep_flag.eq.1) then
               kappa_t = kappa_t/sqrt(1d0-v2)
               kappa_l = kappa_l/(1d0-v2)
               kappa_d = 2*T* kappa_d*(1d0/(2d0*T*(1d0-v2)) 
     &                           - (1d0+v2-sqrt(1d0-v2))/
     &                           (p_mass(i,j)*v2*sqrt(1d0-v2)))
             endif


             rgau_width(1)=sqrt(kappa_t/deltat_lrf)
             rgau_width(2)=sqrt(kappa_t/deltat_lrf)
             rgau_width(3)=sqrt(kappa_l/deltat_lrf)

c This generation of noise is originally at the end of the subroutine
c generate new noise
             do k=1,3
               xi_gauss(k)=rgau(0,rgau_mean,rgau_width(k))
             enddo

c get angles
             call getang(p_px(i,j),p_py(i,j),p_pz(i,j),theta,phi,rr)

c rotate xi_gauss to particle axis:
             call rotbos(theta,phi,0d0,0d0,0d0,
     &            xi_gauss(1),xi_gauss(2),xi_gauss(3),dum)

c update noise

             xi_x(i,j) = w*xi_gauss(1) + (1-w)*xi_x(i,j)
             xi_y(i,j) = w*xi_gauss(2) + (1-w)*xi_y(i,j)
             xi_z(i,j) = w*xi_gauss(3) + (1-w)*xi_z(i,j)
c end of generation of noise

c special case: initialize noise for zeroth time-step:
             if((tstep.eq.1).and.(lstep.eq.1)) then
               rgau_width(1)=sqrt(kappa_t*w/((2d0-w)*deltat_lrf))
               rgau_width(2)=sqrt(kappa_t*w/((2d0-w)*deltat_lrf))
               rgau_width(3)=sqrt(kappa_l*w/((2d0-w)*deltat_lrf))

c these have to be oriented along the particle axis of motion
               xi(1)=rgau(0,rgau_mean,rgau_width(1))
               xi(2)=rgau(0,rgau_mean,rgau_width(2))
               xi(3)=rgau(0,rgau_mean,rgau_width(3))

c rotate xi's to particle axis:

c     get angles
               call getang(p_px(i,j),p_py(i,j),p_pz(i,j),
     &              theta,phi,rr)

c rotate xi's:
               call rotbos(theta,phi,0d0,0d0,0d0,
     &              xi(1),xi(2),xi(3),dum)

               xi_x(i,j)=xi(1)
               xi_y(i,j)=xi(2)
               xi_z(i,j)=xi(3)

c initialize v_bar's and kv/T at zero time-step

               p_vbx(i,j)=(2d0-w)/(2d0*w) * v_x * kappa_d/T
               p_vby(i,j)=(2d0-w)/(2d0*w) * v_y * kappa_d/T
               p_vbz(i,j)=(2d0-w)/(2d0*w) * v_z * kappa_d/T

               p_kvtx(i,j)= kappa_d*v_x/T
               p_kvty(i,j)= kappa_d*v_y/T
               p_kvtz(i,j)= kappa_d*v_z/T
c end initialization
             else ! for normal time steps
c update v^bar's:
               p_vbx(i,j) = 0.5d0*p_px(i,j)/p_p0(i,j) *kappa_d/T
     &                 + (1d0-w)*(0.5d0*p_kvtx(i,j) + p_vbx(i,j))
               p_vby(i,j) = 0.5d0*p_py(i,j)/p_p0(i,j) *kappa_d/T
     &                 + (1d0-w)*(0.5d0*p_kvty(i,j) + p_vby(i,j))
               p_vbz(i,j) = 0.5d0*p_pz(i,j)/p_p0(i,j) *kappa_d/T
     &                 + (1d0-w)*(0.5d0*p_kvtz(i,j) + p_vbz(i,j))
c store v*kappa/T, to serve as "old" value in next time-step
               p_kvtx(i,j) = p_px(i,j)/p_p0(i,j) * kappa_d/T
               p_kvty(i,j) = p_py(i,j)/p_p0(i,j) * kappa_d/T
               p_kvtz(i,j) = p_pz(i,j)/p_p0(i,j) * kappa_d/T

             endif ! end calculation of drag term

           endif              

c update momenta (combine radiation and collision):

           if(flag_rad.eq.1) then ! only radiative energy loss
             p_px(i,j)  = p_px(i,j) - kpx_gluon 
             p_py(i,j)  = p_py(i,j) - kpy_gluon
             p_pz(i,j)  = p_pz(i,j) - klong_gluon

           else             ! combine collisional and radiative
             p_px(i,j) = p_px(i,j) - kpx_gluon - p_vbx(i,j)*deltat_lrf
     &                   + xi_x(i,j)*deltat_lrf
             p_py(i,j) = p_py(i,j) - kpy_gluon - p_vby(i,j)*deltat_lrf
     &                   + xi_y(i,j)*deltat_lrf
             p_pz(i,j) = p_pz(i,j) - klong_gluon - p_vbz(i,j)*deltat_lrf
     &                   + xi_z(i,j)*deltat_lrf
           endif

c update energy
           p_p0(i,j)  = sqrt(p_px(i,j)**2+p_py(i,j)**2+p_pz(i,j)**2
     &                +p_mass(i,j)**2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c reset momentum and energy to the initial value
c for debug purpose
c$$$
c$$$               p_px(i,j)=p_ipx(i,j)
c$$$               p_py(i,j)=p_ipy(i,j)
c$$$               p_pz(i,j)=p_ipz(i,j)
c$$$               p_p0(i,j)=p_ip0(i,j)
c$$$               
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           dEsum=dEsum+(p_p0(i,j)-e_old)**2
           des_cntr=des_cntr+1
c boost back into computational frame

           if(static.ne.1) then
             call rotbos(0d0,0d0,betax,betay,betaz,
     &       p_px(i,j),p_py(i,j),p_pz(i,j),p_p0(i,j))
           endif            ! static.ne.1
c record cell information, the last interaction will be kept
           Thydro(i,j)=T
           edensity(i,j)=edensity0
           sdensity(i,j)=sdensity0
           if(static.ne.1) then
             c_vx(i,j)=betax
             c_vy(i,j)=betay
             c_vz(i,j)=betaz
           endif
         endif               ! particle in medium
      
c propagate particle:
         energ = sqrt(p_px(i,j)**2+p_py(i,j)**2+p_pz(i,j)**2
     &            +p_mass(i,j)**2)
         p_p0(i,j)= energ
         p_r0(i,j)  = tau_p+deltat
         p_rx(i,j)  = p_rx(i,j) + p_px(i,j)/energ*deltat
         p_ry(i,j)  = p_ry(i,j) + p_py(i,j)/energ*deltat
         p_rz(i,j)  = p_rz(i,j) + p_pz(i,j)/energ*deltat


c debug
         if(abs(p_rz(i,j)).gt.tau_p+deltat) then
           write(6,*) ' tachyon scaled! ',p_rz(i,j),
     &                  tau_p,p_pz(i,j)/energ
          p_rz(i,j)=sign(tau_p+deltat-1d-10,p_rz(i,j))
          write(6,*) ' new rz ',p_rz(i,j) 
        endif
        p_reta(i,j)=0.5d0*log((p_r0(i,j)+p_rz(i,j))/
     &              (p_r0(i,j)-p_rz(i,j)))

 11      continue
 10   continue

c      write(6,*)    '    time in dNg_over_dt table   : ', time_tg

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine output(pevent)

      implicit none

      include 'df_coms.f'
      include 'ucoms.f'

      integer pevent,i,j,count_out
      character*77 file20
      double precision vers

      IF(pevent.EQ.1) THEN
          VERS=1D0
          file20='      '
          call getenv('ftn20',file20)
          if (file20(1:4).ne.'    ') then
             OPEN(20,FILE=file20,STATUS='unknown',FORM='FORMATTED')
          else
             OPEN(20,FILE='HPLOSC.DAT',STATUS='UNKNOWN')
          endif

c avoid writing out the head from the second timestep          
          if (out_skip.le.1.or.tstep.eq.(out_skip-6)) then
             WRITE(20,2000)
             WRITE(20,2100) VERS,AP,ZP,AT,ZT,reffram,ECM,0
          endif
       
      ENDIF

      do 10 j=1,evsamp
          if(out_skip.le.1) then
             WRITE(20,2300) evsamp*(pevent-1)+j,npt,bimp,0d0,
     &            1,1,evsamp
          else
             WRITE(20,2300) evsamp*(pevent-1)+j,npt,bimp,0d0,
     &            int(ntsteps/out_skip),int(tstep/out_skip+1),evsamp
          endif

          write(6,*) ' ->  writing output for event# ',
     &                 evsamp*(pevent-1)+j

          do 2401 i=1,npt
              if(corr_flag.eq.1) then
                 count_out=int(i/2d0+0.6d0)
              else
                 count_out=i
              endif
c for e-by-e study, rotate the final momentum by the (par) plane angle
c flow velocity is also rotated
c position vector should also be rotated but not useful yet
              if(out_skip.le.1.and.ebe_flag.eq.1) then
                 call rot2D(-par_plane_phi,p_px(i,j),p_py(i,j))
                 call rot2D(-par_plane_phi,c_vx(i,j),c_vy(i,j))
              endif ! end rotation
              WRITE(20,2420) count_out,pid(i,j),
     &                       p_px(i,j),p_py(i,j),p_pz(i,j),
     &                       p_p0(i,j),p_mass(i,j),p_rx(i,j),p_ry(i,j),
     &                       p_rz(i,j),p_r0(i,j),Thydro(i,j),c_vx(i,j),
     &                       c_vy(i,j),c_vz(i,j),edensity(i,j),
     &                       p_ip0(i,j),p_ipT(i,j),p_wt(i,j)
 2401     continue
 10   continue

 2000 FORMAT('OSC1997A'/'final_id_p_x')
 2100 FORMAT('HPL/ABM ',5X,F5.3,2X,'(',I3,',',I6,')+(',I3,',',I6,')',
     &2X,A4,D10.4,2X,I8)
 2300 FORMAT(I10,2X,I10,2X,F8.3,2X,F8.3,2x,i4,2x,i4,2x,i7)
c sab, oscar modifications
 2420 FORMAT(I10,2X,I10,17(2X,D12.6))

      return

      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE rotbos(THE,PHI,BEX,BEY,BEZ,p1,p2,p3,p4)
c
c INPUT: the,phi,bex,bey,bez,p
c OUTPUT: p
c 1)rotate 4-vector p according to the and phi 2/ boost 4-vector p
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      double precision P(4),BEX,BEY,BEZ,GA,BEP,GABEP,rot(3,3),the,phi
      double precision p1,p2,p3,p4,bb2
      integer j

      IF(THE**2+PHI**2.GT.1d-10) THEN
C...ROTATE
        ROT(1,1)=COS(THE)*COS(PHI)
        ROT(1,2)=-SIN(PHI)
        ROT(1,3)=SIN(THE)*COS(PHI)
        ROT(2,1)=COS(THE)*SIN(PHI)
        ROT(2,2)=COS(PHI)
        ROT(2,3)=SIN(THE)*SIN(PHI)
        ROT(3,1)=-SIN(THE)
        ROT(3,2)=0.
        ROT(3,3)=COS(THE)
        DO 108 J=1,3
 108       P(J)=ROT(J,1)*P1+ROT(J,2)*P2+ROT(J,3)*P3
        p(4)=p4
      else
        p(1)=p1
        p(2)=p2
        p(3)=p3
        p(4)=p4
      ENDIF

      bb2=BEX**2+BEY**2+BEZ**2
      IF(bb2.GT.1d-10) THEN
C...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
        GA=1D0/DSQRT(1D0-bb2)
        BEP=BEX*P(1)+BEY*P(2)+BEZ*P(3)
        GABEP=GA*(GA*BEP/(1D0+GA)+P(4))
        P(1)=P(1)+GABEP*BEX
        P(2)=P(2)+GABEP*BEY
        P(3)=P(3)+GABEP*BEZ
        P(4)=GA*(P(4)+BEP)
      ENDIF

      p1=p(1)
      p2=p(2)
      p3=p(3)
      p4=p(4)

      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine getang(x,y,z,th,ph,r)
c
c gives spherical coordinates of cartesian 3-vector $(x,y,z)$
c
c input : 3-vector x,y,z
c output: angles {\tt th}($\vartheta$), {\tt ph}($\varphi$) and radius {\tt r}
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      double precision x,y,z,th,ph,cut,r
      parameter (cut=1d-9)
      if(abs(x).lt.cut.and.abs(y).lt.cut) then
         ph=0d0
      else
         ph=datan2(y,x)
      endif
      r=sqrt(x*x+y*y+z*z)
      th=dacos(z/max(r,cut))
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Sample initial charms as Moore and Teaney did

      subroutine MTCharm(par_id,MTpx,MTpy,MTpz,MTp0,MTmass,pT2_weight)
      
      implicit none

      include 'df_coms.f'
      include 'ucoms.f'

      integer par_id
      double precision rlu,pT,lambda,ialpha,phi
      double precision pT_max,pT2,rapidity
      double precision MTpx,MTpy,MTpz,MTp0,MTmass,pT2_weight

      pT_max = 30.d0
 
      if(abs(par_id).eq.4) then
         lambda=2.104d0
         ialpha=3.901d0
      elseif(abs(par_id).eq.5) then
         lambda=7.502d0
         ialpha=4.858d0
      else
         write(6,*) "particle information not available ... "
         write(6,*) "terminating ..."
         stop
      endif

      pT2 = pT_max*pT_max*rlu(0)
      pT2_weight = 100000.d0/(pT2+lambda*lambda)**ialpha

c      pT = 30.d0
c      pT2 = pT*pT

      pT = sqrt(pT2)
      phi = rlu(0)*2.0d0*PI
      MTpx = pT*cos(phi)
      MTpy = pT*sin(phi)

      rapidity = 2.d0*rlu(0)-1.0d0
      MTpz = sqrt(pT2+MTmass*MTmass)*sinh(rapidity)
      MTp0 = sqrt(pT2+MTpz*MTpz+MTmass*MTmass)
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine read_wt_table(wt_flag)

      implicit none
      include 'df_coms.f'
      include 'ucoms.f'

      double precision dummy_pT1,dummy_pT2,wtsum1,wtsum2
      integer wt_flag,i
      character*77 file11,file12

      wtsum1=0d0
      wtsum2=0d0

      file11='      '
      call getenv('ftn11',file11)
      if (file11(1:4).ne.'    ') then
         OPEN(UNIT=11,FILE=file11,STATUS='OLD',FORM='FORMATTED')
      else
         OPEN(UNIT=11,FILE='pTc_wt.dat',STATUS='OLD',FORM='FORMATTED')
      endif

      file12='      '
      call getenv('ftn12',file12)
      if (file12(1:4).ne.'    ') then
         OPEN(UNIT=12,FILE=file12,STATUS='OLD',FORM='FORMATTED')
      else
         OPEN(UNIT=12,FILE='pTb_wt.dat',STATUS='OLD',FORM='FORMATTED')
      endif

      i=1
      do while(i.le.wt_num)
         read(unit=11,fmt=*,err=2198,end=2197) dummy_pT1,cpT_wt(i)
         read(unit=12,fmt=*,err=2198,end=2197) dummy_pT2,bpT_wt(i)
         cpT_wt(i)=dummy_pT1*cpT_wt(i)
         bpT_wt(i)=dummy_pT2*bpT_wt(i)
         wtsum1=wtsum1+cpT_wt(i)
         wtsum2=wtsum2+bpT_wt(i)
         i=i+1
      enddo 

c this was a d\sigma/d^2p_T/d\eta table
c we change it to dN/dp table normalized to 1
c by default, the bin size is 0.5~GeV, might be changed due to request

c calculate total cross section of HQ production
      sigma_ctot=wtsum1*2d0*PI*wt_int*2d0*eta_cut
      sigma_btot=wtsum2*2d0*PI*wt_int*2d0*eta_cut

      write(23,*) "# ",wtsum1,wtsum2
      write(23,*) "# ",sigma_ctot,sigma_btot

      i=1
      do while(i.le.wt_num)
         cpT_wt(i)=cpT_wt(i)/wtsum1/wt_int
         bpT_wt(i)=bpT_wt(i)/wtsum2/wt_int
         write(23,*) wt_int*i,cpT_wt(i),bpT_wt(i)
         i=i+1
      enddo 

      wt_flag=0
      write(6,*) "weight table has been read in successfully :)"
      return

 2197 continue
      wt_flag=1
      write(6,*) 'ERROR: EOF reached in weight table'
      write(6,*) 'terminating ...'
      return

 2198  continue
       wt_flag=-1
       write(6,*) 'READ-ERROR in weight table'
       write(6,*) 'terminating ...'
       return

       end

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
         r0(numPart) = 0d0
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c rotation in a 2 dimensional plane

      SUBROUTINE rot2D(theta,x,y)

      implicit none
      double precision theta,x,y
      double precision x0,y0

      if(theta.gt.1d-6) then
         x0 = x
         y0 = y
         x = x0*cos(theta)-y0*sin(theta)
         y = x0*sin(theta)+y0*cos(theta)
      endif

      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine pQCDwt(par_id,MTpx,MTpy,MTpz,MTp0,MTmass,pT_weight)
      
      implicit none

      include 'df_coms.f'
      include 'ucoms.f'

      integer par_id,wt_i
      double precision rlu,pT,phi
      double precision pT2,rapidity,pT_len
      double precision MTpx,MTpy,MTpz,MTp0,MTmass,pT_weight

      pT_len = pT_init_max-pT_init_min

      pT=pT_init_min+pT_len*rlu(0)

      if(reweight.eq.2) then ! record initial pT
         pT_weight=pT_len
      elseif(reweight.eq.1) then ! calculate weight from initial spectrum
         wt_i=int((pT-wt_Tab_min)/wt_int+1.5d0)
         if(wt_i.lt.1.or.wt_i.gt.wt_num) then ! debug
            pT_weight=0d0
            write(6,*) "Warning: initial pT is out of weight range!"
         elseif(abs(par_id).eq.4) then
            pT_weight=cpT_wt(wt_i)*pT_len
         elseif(abs(par_id).eq.5) then
            pT_weight=bpT_wt(wt_i)*pT_len
         else
            write(6,*) "particle information not available ... "
            write(6,*) "terminating ..."
            stop
         endif
      else
         write(6,*) "Inappropriate value for reweight ..."
         write(6,*) "terminating ..."
         stop
      endif

      pT2=pT*pT 
      phi = rlu(0)*2.0d0*PI
      MTpx = pT*cos(phi)
      MTpy = pT*sin(phi)

      rapidity = 2d0*eta_cut*rlu(0)-eta_cut
      MTpz = sqrt(pT2+MTmass*MTmass)*sinh(rapidity)
      MTp0 = sqrt(pT2+MTpz*MTpz+MTmass*MTmass)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
