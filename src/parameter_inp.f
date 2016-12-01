      subroutine parameterRead

      implicit none

      character*10 flag
      character*78 inputstr
      integer open_status
      include 'df_coms.f'
      include 'ucoms.f'

c set default values
      NUMSAMP=1
      flag_rad=2 ! 0-collsional; 1-radiative; 2-both.

      iflav=4     ! flavor of parton to select
      exp_setup=1 ! 1 for LHC 2.76~TeV and 2 for RHIC 200~GeV
      out_skip=0  ! if 1, only output the last time step
      reweight=2  ! 0: no re-weighting
                  ! 1: re-sample pT distribution according to wt table
                  ! 2: use uniform distribution
                  ! if 0, no evolution, output init. condition
      corr_flag=0 ! 1 for calculation of correlation
      HQ_input=1  ! 1:read in oscar file, 2:read in OSU file
                  ! 3:read in a list of x-y positions only
      ebe_flag=2  ! 0:smooth initial condition, 1:e-by-e study
                  ! 2:no rotation -- plain 
                  ! (use 2 when HQ_input is not 2)
      qhat_TP=0   ! 0: constant qhat/T^3
                  ! 1: read a table of qhat/T^3 as functios of T,p

c diffusion/noise parameters:
      D2piT=6d0
      KFactor=1d0
      
      KPamp=7d0
      KPsig=5d0
      KTamp=2d0
      KTsig=0.05d0
      preKT=0.7d0

      qhatMin= 0.2d0
      qhatSlope = 2.0d0
      qhatPower = 0d0
      preP = 1d0     
      qhatC = 0d0
 
      alpha=6.2832d0/D2piT
      eta_cut=0.5d0   ! rapidity cut, symmetric inteval
      pT_init_min=0.25d0
      pT_init_max=70.25d0

      TLow=0.165d0
      THigh=0.4d0
      qhatLow=5.5d0
      qhatHigh=5.5d0

c other paramters:
      static=2 ! 0-Chiho's hydro, 1-Static, 2-OSU hydro
      tsteps_cut=140
      static_cool=0
      T_static=0.30d0
      temp_cut=1 ! 0 for no cut, 1 for cut
      Tcut_critical=0.165d0
      wt_num=140 ! only meaningful for reweight=1
      wt_int=0.5d0 ! bin size of the weight table
      wt_Tab_min=0.5d0

      num_binary=0 ! 0: use number of xy points as npart
                   ! <0: use -num_binary as npart
                   ! >0: calculate npart with cross section
      sigma_c0=0.68432d0
      sigma_b0=0.056263d0

      OPEN(UNIT=15,FILE="parameters_df.dat",
     &     STATUS='OLD',FORM='FORMATTED',IOSTAT=open_status)

      if (open_status.ne.0) then
         write(6,*) "Error occurs during reading parameters_df.dat."
         write(6,*) "Use default parameters."
         goto 2
      endif

c read input lines
 1    continue
      read(15,99) flag,inputstr
 99   format(1A10,1A78)

c # : treat line as a comment
      if(flag(1:1).eq.'#') goto 1
      if(flag(1:4).eq."    ") goto 1
c xx: treat line as end of input marker
      if(flag(1:2).eq.'xx') goto 2

      SELECT CASE (flag)
         CASE ("NUMSAMP...")
            read(inputstr,fmt=*,err=88,end=88) NUMSAMP
         CASE ("iflav.....")
            read(inputstr,fmt=*,err=88,end=88) iflav
         CASE ("out_skip..")
            read(inputstr,fmt=*,err=88,end=88) out_skip
         CASE ("reweight..")
            read(inputstr,fmt=*,err=88,end=88) reweight
         CASE ("eta_cut...")
            read(inputstr,fmt=*,err=88,end=88) eta_cut
         CASE ("qhat_TP...")
            read(inputstr,fmt=*,err=88,end=88) qhat_TP
         CASE ("corr_flag.")
            read(inputstr,fmt=*,err=88,end=88) corr_flag
         CASE ("HQ_input..")
            read(inputstr,fmt=*,err=88,end=88) HQ_input
         CASE ("ebe_flag..")
            read(inputstr,fmt=*,err=88,end=88) ebe_flag
         CASE ("exp_setup.")
            read(inputstr,fmt=*,err=88,end=88) exp_setup
         CASE ("D2piT.....")
            read(inputstr,fmt=*,err=88,end=88) D2piT
            alpha=6.2832d0/D2piT
         CASE ("TLow......")
            read(inputstr,fmt=*,err=88,end=88) TLow
         CASE ("THigh.....")
            read(inputstr,fmt=*,err=88,end=88) THigh
         CASE ("qhatLow...")
            read(inputstr,fmt=*,err=88,end=88) qhatLow
         CASE ("qhatHigh..")
            read(inputstr,fmt=*,err=88,end=88) qhatHigh
         CASE ("KFactor...")
            read(inputstr,fmt=*,err=88,end=88) KFactor
         CASE ("KPamp.....")
            read(inputstr,fmt=*,err=88,end=88) KPamp
         CASE ("KPsig.....")
            read(inputstr,fmt=*,err=88,end=88) KPsig
         CASE ("KTamp.....")
            read(inputstr,fmt=*,err=88,end=88) KTamp
         CASE ("KTsig.....")
            read(inputstr,fmt=*,err=88,end=88) KTsig
         CASE ("preKT.....")
            read(inputstr,fmt=*,err=88,end=88) preKT
         CASE ("qhatMin...")
            read(inputstr,fmt=*,err=88,end=88) qhatMin
         CASE ("preP......")
            read(inputstr,fmt=*,err=88,end=88) preP
         CASE ("qhatC.....")
            read(inputstr,fmt=*,err=88,end=88) qhatC
         CASE ("qhatSlope.")
            read(inputstr,fmt=*,err=88,end=88) qhatSlope
         CASE ("qhatPower.")
            read(inputstr,fmt=*,err=88,end=88) qhatPower
         CASE ("pT_min....")
            read(inputstr,fmt=*,err=88,end=88) pT_init_min
         CASE ("pT_max....")
            read(inputstr,fmt=*,err=88,end=88) pT_init_max
         CASE ("flag_rad..")
            read(inputstr,fmt=*,err=88,end=88) flag_rad
         CASE ("static....")
            read(inputstr,fmt=*,err=88,end=88) static
         CASE ("tsteps_cut")
            read(inputstr,fmt=*,err=88,end=88) tsteps_cut
         CASE ("stat_cool.")
            read(inputstr,fmt=*,err=88,end=88) static_cool
         CASE ("T_static..")
            read(inputstr,fmt=*,err=88,end=88) T_static
         CASE ("temp_cut..")
            read(inputstr,fmt=*,err=88,end=88) temp_cut
         CASE ("Tcut......")
            read(inputstr,fmt=*,err=88,end=88) Tcut_critical
         CASE ("wt_num....")
            read(inputstr,fmt=*,err=88,end=88) wt_num
         CASE ("wt_int....")
            read(inputstr,fmt=*,err=88,end=88) wt_int   
         CASE ("wt_Tab_min")
            read(inputstr,fmt=*,err=88,end=88) wt_Tab_min
         CASE ("sigma_c0..")
            read(inputstr,fmt=*,err=88,end=88) sigma_c0
         CASE ("sigma_b0..")
            read(inputstr,fmt=*,err=88,end=88) sigma_b0
         CASE ("num_binary")
            read(inputstr,fmt=*,err=88,end=88) num_binary
         CASE DEFAULT
            write(6,*) flag,"NOT a valid parameter!"
            write(6,*) "Terminating ..."
            stop
      END SELECT


      goto 1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 2    continue


      call readInputFromCML

      KTsig=KTsig*Tcut_critical
      qhatA=(qhatLow-qhatHigh)/(TLow-THigh)
      qhatB=(TLow*qhatHigh-THigh*qhatLow)/(TLow-THigh)

      if(open_status.eq.0)
     &   write(6,*) "Parameters for Langevin has been read in."
      if(.TRUE.) then
         write(6,*) "flag_rad: ",flag_rad
         write(6,*) "iflav: ",iflav
         write(6,*) "KPamp: ",KPamp
         write(6,*) "KPsig: ",KPsig
         write(6,*) "KTamp: ",KTamp
         write(6,*) "KTsig: ",KTsig
         write(6,*) "preKT: ",preKT
         write(6,*) "qhatMin: ", qhatMin
         write(6,*) "qhatSlope: ", qhatSlope
         write(6,*) "qhatPower: ", qhatPower
         write(6,*) "preP: ", preP
         write(6,*) "qhatC: ", qhatC
         write(6,*) "D(2piT): ",D2piT
         write(6,*) "alpha: ", alpha
         write(6,*) "qhatA: ",qhatA
         write(6,*) "qhatB: ",qhatB
         write(6,*) "HQ_input: ",HQ_input
         write(6,*) "NUMSAMP: ",NUMSAMP
         write(6,*) "reweight: ",reweight
         write(6,*) "eta_cut: ",eta_cut
         write(6,*) "pT_min: ",pT_init_min
         write(6,*) "pT_max: ",pT_init_max
         write(6,*) "wt_num: ",wt_num
         write(6,*) "wt_int: ",wt_int
         write(6,*) "wt_Tab_min: ",wt_Tab_min
         write(6,*) "out_skip: ",out_skip
         write(6,*) "exp_setup: ",exp_setup
         write(6,*) "ebe_flag: ",ebe_flag
         write(6,*) "corr_flag: ",corr_flag
         write(6,*) "static: ",static
         write(6,*) "static_cool: ",static_cool
         write(6,*) "T_static: ",T_static
         write(6,*) "tsteps_cut: ",tsteps_cut 
         write(6,*) "temp_cut: ",temp_cut
         write(6,*) "Tcut_critical: ",Tcut_critical
         write(6,*) "sigma_c0: ",sigma_c0
         write(6,*) "sigma_b0: ",sigma_b0
         write(6,*) "num_binary: ",num_binary
      endif 
      return

c error-exit
 88   write(6,*) 'Syntax-error in the parameters_df.dat ...'
      write(6,*) 'Terminating ...'
      stop
      end
      
      
! Yingru
      Subroutine readInputFromCML()
!     Purpose:
!     Read inputs from command line
      Implicit None
      Include 'df_coms.f'
      Include 'ucoms.f'
      
      Integer QNum, ArgIndex   !QNum is the total number of arguments, ArgIndex gives the index to the one currently reading
      
      Character*60 :: buffer
      Character*20 :: varName
!      Integer IResult
      Double Precision DResult     
      
!      write(6, *) "here reading from command line" 
      QNum = iargc()
      
      Do ArgIndex = 1, QNum
        call getarg(ArgIndex, buffer)
!        write(6, *) buffer
        call processAssignment(buffer,"=", varName, DResult)
        
        if (varName .EQ. "kpamp") KPamp = DResult  ! the amplitude for kP
        if (varName .EQ. "kpsig") KPsig = DResult  ! the sigma for kP       
        if (varName .EQ. "ktamp") KTamp = DResult  ! the amplitude for kT
        if (varName .EQ. "ktsig") KTsig = DResult  ! the sigma for kT
        if (varName .EQ. "prekt") preKT = DResult  ! the sigma for preKT
        if (varName .EQ. "qhatmin") qhatMin = DResult   !for linear parameterization
        if (varName .EQ. "qhatslope") qhatSlope = DResult  ! for linear parameterization
        if (varName .EQ. "qhatpower") qhatPower = DResult   ! for linear parameterization
        if (varName .EQ. "prep") preP = DResult
        if (varName .EQ. "qhatc") qhatC = DResult
        if (varName .EQ. "d2pit") then
                D2piT = DResult
                alpha = 6.2832d0/D2piT
        endif

      End Do  ! ArgIndex
      
      End Subroutine
      
      
      Subroutine processAssignment(string, separator, 
     &                            varName, DResult)
!     This subroutine process a string assignment
!     First it seprate string into LHS and RHS according to separator.
!     The the LHS is coverted into variable using only lower case letters,
!     The RHS is converted into numerical values

      Implicit None
      
      Character (*) :: string, varName
      Character*60 :: LHS, RHS
      Character separator
      Double Precision DResult
      
      Integer break_here, I, cha
      
      varName = ""
      break_here = index(string, separator)
      LHS = adjustl(string(:break_here-1))
      RHS = adjustl(string(break_here+1:))
      
      ! convert LHS to lower case:
      Do I = 1, len_trim(LHS)
        cha = ichar(LHS(I:I))
        If (cha>=65 .AND. cha<90) Then
          varName(I:I) = char(cha+32)
        Else
          varName(I:I) = LHS(I:I)
        Endif
      EndDo
      
!      write(6, *) varName 
      ! convert RHS to numerics(here only take double precision):
      Read(RHS, fmt='(f15.8)') DResult
      
      End Subroutine
      
