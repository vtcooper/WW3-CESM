!/ ------------------------------------------------------------------- /
      MODULE WAV_COMP_MCT
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     A generic cpl7 interface for WAVEWATCH III
!     using input fields from cpl7.
!
!  2. Method :
!
!     MCT component for the actual wave model (W3WAVE).
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!       NHMAX   I.P.  Maximum number of homogeneous fields.
!
!       NDSO    Int.  General output unit number (shell only).
!       NDSE    Int.  Error output unit number (shell only).
!       NDST    Int.  Test output unit number (shell only).
!       FLH     L.A.  Flags for homogeneous fields.
!       NH      I.A.  Number of times for homogeneous fields.
!       THO     I.A.  Times of homogeneous fields.
!       TIME0   I.A.  Starting time.  
!       TIMEN   I.A.  Ending time.    
!     ----------------------------------------------------------------
!
!       NDS, NTRACE, ..., see W3WAVE
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set nummber of data structures
!      W3SETG    Subr.   Id.    Point to data structure.
!      W3NDAT    Subr. W3WDATMD Set nummber of data structures
!      W3SETW    Subr.   Id.    Point to data structure.
!      W3NMOD    Subr. W3ADATMD Set nummber of data structures
!      W3NAUX    Subr.   Id.    Point to data structure.
!      W3NOUT    Subr. W3ODATMD Set nummber of data structures
!      W3SETO    Subr.   Id.    Point to data structure.
!      W3NINP    Subr. W3IDATMD Set nummber of data structures
!      W3SETI    Subr.   Id.    Point to data structure.
!
!      STME21    Subr. W3TIMEMD Print date and time readable.
!      DSEC21    Func.   Id.    Difference between times.
!      TICK21    Subr.   Id.    Increment time.
!
!      W3INIT    Subr. W3INITMD Wave model initialization.
!      W3WAVE    Subr. W3WAVEMD Wave model.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!     - Checks on I-O.
!     - Check on time interval.
!
!  7. Remarks :
!
!     - A rigourous input check is made in W3INIT.
!     - See W3WDAS for documentation on the set-up of the data
!       assimilation.
!
!  8. Structure :
!
!     ----------------------------------------------------------------
!
!     wav_comp_init 
!
!        0.   Set up data structures.                ( W3NMOD, etc. )
!        1.   I-O setup.
!          a  For shell.
!          b  For WAVEWATCH III.
!          c  Local parameters.
!        2.   Define input fields
!        3.   Set time frame.
!        4.   Define output 
!          a  Loop over types, do
!        +--------------------------------------------------------+
!        | b    Process standard line                             |
!        | c    If type 1: fields of mean wave parameters         |
!        | d    If type 2: point output                           |
!        | e    If type 3: track output                           |
!        | f    If type 4: restart files                          |
!        | g    If type 5: boundary output                        |
!        | h    If type 6: separated wave fields                  |
!        +--------------------------------------------------------+
!        5.   Initialzations
!          a  Wave model.                              ( W3INIT )
!          d  Set field times.
!
!     wav_comp_run 
!
!        7.   Run model with input
!             Do until end time is reached
!        +--------------------------------------------------------+
!        | a  Determine next time interval and input fields.      |
!        |   1  Preparation                                       |
!        |      Loop over input fields                            |
!        | +------------------------------------------------------|
!        | | 2  Check if update is needed                         |
!        | | 4  Update next ending time                           |
!        | +------------------------------------------------------|
!        | b  Run wave model.                          ( W3WAVE ) |
!        | c  If requested, data assimilation.         ( W3WDAS ) |
!        | d  Final output if needed.                  ( W3WAVE ) |
!        | e  Check time                                          |
!        +--------------------------------------------------------+
!
!     wav_comp_fin
!
!     ----------------------------------------------------------------
!
!  9. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!       !/MPI   Id.
!
!       !/LLG   Spherical grid.
!       !/XYG   Cartesian grid.
!       !/MGW   Moving grid wind correction.
!       !/MGP   Moving grid propagation correction.
!
!       !/T     Enable test output.
!       !/O7    Echo input homogeneous fields.
!
!       !/NCO   NCEP NCO modifications for operational implementation.
!
!       !/F90   Timer function included for F90.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      use w3gdatmd
      use w3wdatmd, only: time, w3ndat, w3dimw, w3setw, UST, USTDIR
      use w3adatmd
      use w3idatmd, only: flags, w3seti, w3ninp 
      use w3odatmd, only: w3nout, w3seto, naproc, iaproc, napout, naperr,             &
                          nogrd, idout, fnmpre, iostyp
!/
      use w3initmd
      use w3wavemd
!/
      use w3iopomd
      use w3timemd
      use w3cesmmd, only : casename

      use esmf_mod
      use mct_mod 
      use seq_flds_mod
      use shr_sys_mod      , only : shr_sys_flush 
      use shr_kind_mod     , only : in=>shr_kind_in, r8=>shr_kind_r8, &
                                    cs=>shr_kind_cs, cl=>shr_kind_cl
      use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
      use seq_timemgr_mod  , only : seq_timemgr_eclockgetdata
      use seq_infodata_mod , only : seq_infodata_type, seq_infodata_getdata, seq_infodata_putdata, &
                                    seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                    seq_infodata_start_type_brnch
      use shr_file_mod     , only : shr_file_setlogunit, shr_file_setloglevel, &
                                    shr_file_getlogunit, shr_file_getloglevel, &
                                    shr_file_getunit, shr_file_setio
!
      implicit none
!
      public :: wav_init_mct
      public :: wav_run_mct
      public :: wav_final_mct

      private

      private :: wav_setgsmap_mct
      private :: wav_domain_mct

      integer,save :: index_x2w_Sa_u
      integer,save :: index_x2w_Sa_v
      integer,save :: index_x2w_Sa_tbot
      integer,save :: index_x2w_So_u
      integer,save :: index_x2w_So_v
      integer,save :: index_x2w_So_t
      integer,save :: index_x2w_Si_ifrac

      integer,save :: index_w2x_Sw_langnum
      integer,save :: index_w2x_Fw_taux
      integer,save :: index_w2x_Fw_tauy

      include "mpif.h"
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CONTAINS
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    SUBROUTINE WAV_INIT_MCT( EClock, cdata, x2w, w2x, NLFilename )
      !/
      !/ ------------------------------------------------------------------- /
      !/ Parameter list
      !/
      TYPE(ESMF_CLOCK), INTENT(IN)    :: ECLOCK
      TYPE(SEQ_CDATA) , INTENT(INOUT) :: CDATA
      TYPE(MCT_AVECT) , INTENT(INOUT) :: X2W, W2X
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: NLFILENAME ! NAMELIST FILENAME
      !/
      !/ ------------------------------------------------------------------- /
      !/ Local PARAMETER statements
      !/
      INTEGER, PARAMETER  :: NHMAX =    200
      !/
      !/ ------------------------------------------------------------------- /
      !/ Local parameters
      !/
      integer :: n
      integer :: compid
      integer :: mpi_comm
      integer :: lsize
      integer :: shrlogunit, shrloglev  
      integer :: hh,mm,ss
      integer :: dtime_sync        ! integer timestep size
      integer :: start_ymd         ! start date (yyyymmdd)
      integer :: start_tod         ! start time of day (sec)
      integer :: stdout

      character(CL)            :: starttype
      type(mct_gsmap), pointer :: gsmap
      type(mct_ggrid), pointer :: dom
      type(seq_infodata_type), pointer :: infodata   ! input init object

      integer             :: ndso, ndse, nds(13), ntrace(2), time0(2), &
                             timen(2), odat(30), nh(4), iprt(6)
      integer             :: i,j,npts
      integer             :: ierr, ierr_mpi
      real                :: a(nhmax,4)
      real, allocatable   :: x(:), y(:)
      logical             :: flgrd(nogrd), prtfrm, flt 
!
      character(len=3)    :: idstr(8), idtst
      character(len=10)   :: pn
      character(len=13)   :: idflds(8)
      character(len=20)   :: strng
      character(len=23)   :: dtme21
      character(len=30)   :: idotyp(6)
      character(len=10), allocatable :: pnames(:)
      !/
      !/ ------------------------------------------------------------------- /
      !/
      DATA IDFLDS / 'water levels ' , 'currents     ' ,               &
                    'winds        ' , 'ice fields   ' ,               &
                    'mean param.  ' , '1D spectra   ' ,               &
                    '2D spectra   ' , 'moving grid  ' /
      DATA IDOTYP / 'Fields of mean wave parameters' ,                &
                    'Point output                  ' ,                &
                    'Track point output            ' ,                &
                    'Restart files                 ' ,                &
                    'Nesting data                  ' ,                &
                    'Partitioned wave field data   ' /
      DATA IDSTR  / 'LEV', 'CUR', 'WND', 'ICE', 'DT0', 'DT1', 'DT2',  &
                    'MOV' /
      !
      !--------------------------------------------------------------------
      ! Set up data structures
      !--------------------------------------------------------------------
      !
      call w3nmod ( 1, 6, 6 )
      call w3ndat (    6, 6 )
      call w3naux (    6, 6 )
      call w3nout (    6, 6 )
      call w3ninp (    6, 6 )
      !
      call w3setg ( 1, 6, 6 )
      call w3setw ( 1, 6, 6 )
      call w3seta ( 1, 6, 6 )
      call w3seto ( 1, 6, 6 )
      call w3seti ( 1, 6, 6 )
      !
      !--------------------------------------------------------------------
      ! Initialize mpi
      !--------------------------------------------------------------------
      !
      call seq_cdata_setptrs(cdata, id=compid, mpicom=mpi_comm, &
           gsmap=gsmap, dom=dom, infodata=infodata)

      call mpi_comm_size(mpi_comm, naproc, ierr)
      call mpi_comm_rank(mpi_comm, iaproc, ierr)
      iaproc = iaproc + 1
      !
      !--------------------------------------------------------------------
      ! IO set-up
      !--------------------------------------------------------------------
      !
      napout = 1
      naperr = 1
      !
      ! 1.b For WAVEWATCH III (See W3INIT) ??? ask adrean if i am missing something
      !
      ! The following units are referenced in module w3initmd
      ! NDS(1) ! OUTPUT LOG: General output unit number ("log file") (NDS0)
      ! NDS(2) ! OUTPUT LOG: Error output unit number (NDSE)
      ! NDS(3) ! OUTPUT LOG: Test output unit number (NDST)
      ! NDS(4) ! OUTPUT LOG: Unit for 'direct' output (SCREEN) 
      !
      ! NDS(5) ! INPUT: mod_def.ww3 file (model definition) unit number
      ! NDS(9) ! INPUT: unit for read in boundary conditions (based on FLBPI) 
      !
      ! The following units are referenced in module w3wavemd for output
      ! NDS( 6) ! OUTPUT DATA: restart(N).ww3 file (model restart) unit number
      ! NDS( 7) ! unit for output for FLOUT(1) flag
      ! NDS( 8) ! unit for output for FLOUT(2) flag
      ! NDS(11) ! unit for output for FLOUT(3) flag
      ! NDS(12) ! unit for output for FLOUT(3) flag
      ! NDS(13) ! unit for output for FLOUT(6) flag
      ! NDS(10) ! unit for output for FLOUT(5) flag
      !
      if (iaproc .eq. napout) then
         stdout = shr_file_getunit()
         call shr_file_setio('wav_modelio.nml',stdout)
      else
         stdout = 6
      endif

      nds( 1) = stdout
      nds( 2) = stdout
      nds( 3) = stdout
      nds( 4) = stdout
      nds( 5) = shr_file_getunit()  
      nds( 6) = shr_file_getunit()
      nds( 7) = shr_file_getunit()
      nds( 8) = shr_file_getunit()
      nds( 9) = shr_file_getunit()
      nds(10) = shr_file_getunit()
      nds(11) = shr_file_getunit()
      nds(12) = shr_file_getunit()
      nds(13) = shr_file_getunit()
      !
      ndso      =  stdout
      ndse      =  stdout
      ntrace(1) =  nds(3)
      ntrace(2) =  10
      !
      ! Redirect share output to wav log
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (stdout)

      if ( iaproc .eq. napout ) write (ndso,900)
      call shr_sys_flush(stdout)
      !
      !--------------------------------------------------------------------
      ! Define input fields
      !--------------------------------------------------------------------
      !
      ! NOTE: no longer read NDSI (ww3_shel.inp) - input is hard-wired below
      flags(1:8) = .false.

      !--------------------------------------------------------------------
      ! Set time frame
      !--------------------------------------------------------------------
      !
      ! TIME0 = from ESMF clock
      ! NOTE - are not setting TIMEN here

      if ( iaproc .eq. napout ) write (ndso,930)
      call shr_sys_flush(stdout)
      !
      call seq_timemgr_EClockGetData(EClock, &
           start_ymd=start_ymd, start_tod=start_tod)

      hh = start_tod/3600
      mm = (start_tod - (hh * 3600))/60
      ss = start_tod - (hh*3600) - (mm*60) 

      time0(1) = start_ymd
      time0(2) = hh*10000 + mm*100 + ss
      call stme21 ( time0 , dtme21 )
      if ( iaproc .eq. napout ) write (ndso,931) dtme21
      call shr_sys_flush(stdout)
      time = time0
      !
      !--------------------------------------------------------------------
      ! Define output type and fields
      !--------------------------------------------------------------------
      !
      iostyp = 1
      write (ndso,940) 'no dedicated output process, any file system '
      call shr_sys_flush(stdout)
      !
      ! TODO - need to enable restart files in run
      ! Actually will need a new restart flag - since all of the ODAT
      ! should be set to 0 - since they are initializated in w3initmd
      ! ODAT    I.A.   I   Output data, five parameters per output type
      !                          1 YYYMMDD for first output.
      !                          2 HHMMSS for first output.
      !                          3 Output interval in seconds.
      !                          4 YYYMMDD for last output.
      !                          5 HHMMSS for last output.
      !                     1-5  Data for OTYPE = 1; gridded fields.
      !                     6-10 Id.  for OTYPE = 2; point output.
      !                    11-15 Id.  for OTYPE = 3; track point output.
      !                    16-20 Id.  for OTYPE = 4; restart files.
      !                    21-25 Id.  for OTYPE = 5; boundary data.
      ! FLGRD   L.A.   I   Flags for gridded output.
      ! NPT     Int.   I   Number of output points
      ! X/YPT   R.A.   I   Coordinates of output points.
      ! PNAMES  C.A.   I   Output point names.
      !
      npts = 0
      allocate ( x(1), y(1), pnames(1) )
      pnames(1) = ' '

      do j=1, 6
         odat(5*(j-1)+3) = 0
         !DEBUG
         ! Hardwire gridded output for now 
         odat(1) = 10101
         odat(2) = 0
         odat(3) = 86400
         odat(4) = 99990101
         odat(5) = 0
         !DEBUG
      end do
                                                                          
      ! Output Type 1: fields of mean wave parameters gridded output

      flgrd( 1) = .true. !   1. depth (m)                              
      flgrd( 2) = .true. !   2. mean current vel (vec, m/s)            
      flgrd( 3) = .true. !   3. mean wind vel (vec, m/s)               
      flgrd( 4) = .true. !   4. air-sea temp diff (deg C)              
      flgrd( 5) = .true. !   5. skin friction vel (scalar, m/s)        
      flgrd( 6) = .true. !   6. significant wave height (m)            
      flgrd( 7) = .true. !   7. mean wave length (m)                   
      flgrd( 8) = .true. !   8. mean wave period (Tn1, s)              
      flgrd( 9) = .true. !   9. mean wave dir (deg: met conv)          
      flgrd(10) = .true. !  10. mean dir spread (deg: )                
      flgrd(11) = .true. !  11. peak freq (Hz)                         
      flgrd(12) = .true. !  12. peak dir (deg: )                       
      flgrd(13) = .true. !  13. peak freq of wind-sea part             
      flgrd(14) = .true. !  14. wind-sea dir (deg: met conv)           
      flgrd(15) = .false.!  15. wave height of partitions               
      flgrd(16) = .false.!  16. peak periods of partitions             
      flgrd(17) = .false.!  17. peak wave length of partitions         
      flgrd(18) = .false.!  18. mean dir of partitions                 
      flgrd(19) = .false.!  19. dir spread of partitions               
      flgrd(20) = .false.!  20. wind-sea frac of partitions             
      flgrd(21) = .false.!  21. wind-sea frac of entire spec           
      flgrd(22) = .false.!  22. number of partitions                   
      flgrd(23) = .true. !  23. average time step (s)                  
      flgrd(24) = .true. !  24. cut-off freq (Hz)                      
      flgrd(25) = .true. !  25. ice concentration (frac)               
      flgrd(26) = .true. !  26. water level (m?)                       
      flgrd(27) = .false.!  27. near-bottom rms exclusion amp          
      flgrd(28) = .false.!  28. near-bottom rms orbital vel            
      flgrd(29) = .false.!  29. radiation stresses                     
      flgrd(30) = .false.!  30. user defined (1)                       
      flgrd(31) = .false.!  31. user defined (2)                       

      if ( iaproc .eq. napout ) then
         flt = .true.
         do i=1, nogrd
            if ( flgrd(i) ) then
               if ( flt ) then
                  write (ndso,1945) idout(i)
                  flt    = .false.
               else
                  write (ndso,1946) idout(i)
               end if
            end if
         end do
         if ( flt ) write (ndso,1945) 'no fields defined'
      end if
      call shr_sys_flush(stdout)
      !
      !--------------------------------------------------------------------
      ! Wave model initializations
      !--------------------------------------------------------------------
      !
      if ( iaproc .eq. napout ) write (ndso,950)
      if ( iaproc .eq. napout ) write (ndso,951) 'wave model ...'
      call shr_sys_flush(stdout)

      call seq_infodata_GetData(infodata,case_name=casename)

      call w3init ( 1, 'ww3', nds, ntrace, odat, flgrd, npts, x, y,   &
           pnames, iprt, prtfrm, mpi_comm )
      call shr_sys_flush(stdout)

      ! overwrite dt values with variables from coupler
      ! is this a problem with any things being set in w3init?
      call seq_timemgr_eclockgetdata(eclock, dtime=dtime_sync )      
      dtmax  = real(dtime_sync)
      dtcfl  = real(dtime_sync) / 3. !todo ??? ask adrean
      dtcfli = real(dtime_sync)      !todo ??? ask adrean
      dtmin  = real(dtime_sync) / 12 !todo ??? ask adrean

      call mpi_barrier ( mpi_comm, ierr_mpi )
      !
      !--------------------------------------------------------------------
      ! cpl7/mct initialization
      !--------------------------------------------------------------------
      !
      ! initialize mct gsmap

      call wav_setgsmap_mct(mpi_comm, compid, gsmap)
      lsize = mct_gsmap_lsize(gsmap, mpi_comm)
      write(stdout,*)'lsize= ',lsize
      call shr_sys_flush(stdout)

      ! initialize mct domain

      call wav_domain_mct(lsize, gsmap, dom)
      
      ! set flags in infodata
      ! wav_prognostic is set to .false. for debugging purposes only
      call seq_infodata_putdata(infodata, wav_present=.true., &
           wav_prognostic=.true., wav_nx=nx, wav_ny=ny)

      ! initialize mct attribute vectors
      
      call mct_avect_init(w2x, rlist=seq_flds_w2x_fields, lsize=lsize)
      call mct_avect_zero(w2x)
      
      call mct_avect_init(x2w, rlist=seq_flds_x2w_fields, lsize=lsize)
      call mct_avect_zero(x2w)

      index_x2w_Sa_u       = mct_avect_indexRA(x2w,'Sa_u'      ,perrWith='quiet')
      index_x2w_Sa_v       = mct_avect_indexRA(x2w,'Sa_v'      ,perrWith='quiet')
      index_x2w_Sa_tbot    = mct_avect_indexRA(x2w,'Sa_tbot'   ,perrWith='quiet')
      index_x2w_So_u       = mct_avect_indexRA(x2w,'So_u'      ,perrWith='quiet')
      index_x2w_So_v       = mct_avect_indexRA(x2w,'So_v'      ,perrWith='quiet')
      index_x2w_So_t       = mct_avect_indexRA(x2w,'So_t'      ,perrWith='quiet')
      index_x2w_Si_ifrac   = mct_avect_indexRA(x2w,'Si_ifrac'  ,perrWith='quiet')

      index_w2x_Sw_langnum = mct_avect_indexRA(w2x,'Sw_langnum',perrWith='quiet')
      index_w2x_Fw_taux    = mct_avect_indexRA(w2x,'Fw_taux'   ,perrWith='quiet')
      index_w2x_Fw_tauy    = mct_avect_indexRA(w2x,'Fw_tauy'   ,perrWith='quiet')
      
      ! add call to gptl timer

      ! end redirection of share output to wav log

      call shr_sys_flush(stdout)
      call shr_file_setlogunit (shrlogunit)
      call shr_file_setloglevel(shrloglev)


900   FORMAT (/15X,'      *** WAVEWATCH III Program shell ***      '/ &
               15X,'==============================================='/)
901   FORMAT ( '  Comment character is ''',A,''''/)

930   FORMAT (/'  Time interval : '/                                  &
           ' --------------------------------------------------')
931   FORMAT ( '       Starting time : ',A)

940   FORMAT (/'  Output requests : '/                                &
           ' --------------------------------------------------'/ &
           '       ',A)

950   FORMAT (/'  Initializations :'/                                 &
           ' --------------------------------------------------')
951   FORMAT ( '       ',A)
1945  FORMAT ( '            Fields   : ',A)
1946  FORMAT ( '                       ',A)

    END SUBROUTINE WAV_INIT_MCT
    
!=====================================================================
!=====================================================================
!=====================================================================

    SUBROUTINE WAV_RUN_MCT(EClock, cdata_w, x2w_w, w2x_w)

      ! Parameters
      !
      type(ESMF_Clock)            ,intent(in)    :: EClock
      type(seq_cdata)             ,intent(inout) :: cdata_w
      type(mct_aVect)             ,intent(inout) :: x2w_w
      type(mct_aVect)             ,intent(inout) :: w2x_w
      !
      !/
      !/ ------------------------------------------------------------------- /
      !/ Local parameters
      !/
      integer :: time0(2), timen(2), ttime(2), ttt(2), ierr, i, j, ix, iy
      integer :: ymd              ! current year-month-day
      integer :: tod              ! current time of day (sec)
      integer :: hh,mm,ss
      integer :: n,jsea,isea
      integer :: ierr_mpi
      real    :: uwind, vwind ! This will not work with real*8???
      real    :: tbot, tocn 
      !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! 7.  Model with input
      ! 7.a Determine next time interval and input fields
      ! 7.a.1 Preparation
      ! 7.a.3 Update time and fields / data
      !
      ! W3ADATM (type WDATA)
      !
      !      Name      Type  Scope    Description
      !     ----------------------------------------------------------------
      !      TIME      I.A.  Public   Valid time for spectra.
      !      TLEV      I.A.  Public   Valid time for water levels.
      !      TICE      I.A.  Public   Valid time for ice.
      !      VA        R.A.  Public   Storage array for spectra.
      !      WLV       R.A.  Public   Water levels.
      !      ICE       R.A.  Public   Ice coverage.
      !      UST       R.A.  Public   Friction velocity (absolute).
      !      USTDIR    R.A.  Public   Friction velocity direction.
      !     ----------------------------------------------------------------
      !
      ! W3ADATM (type WADATA)
      !
      !     Fields of mean wave parameters:
      !
      !      Name      Type  Scope    Description
      !     ----------------------------------------------------------------
      !      DW        R.A.  Public   Water depths.
      !      UA        R.A.  Public   Absolute wind speeds.
      !      UD        R.A.  Public   Absolute wind direction.
      !      U10       R.A.  Public   Wind speed used.
      !      U10D      R.A.  Public   Wind direction used.
      !      AS        R.A.  Public   Stability parameter.
      !      CX/Y      R.A.  Public   Current components.
      !      EMN       R.A.  Public   Mean energy.
      !      FMN       R.A.  Public   Mean frequency.
      !      WNM       R.A.  Public   Mean wavenumber.
      !      AMX       R.A.  Public   Spectral maximum.
      !      CDS       R.A.  Public   Drag coefficient.
      !      Z0S       R.A.  Public   Roughness parameter.
      !      HS        R.A.  Public   Wave Height.
      !      WLM       R.A.  Public   Mean wave length.
      !      TMN       R.A.  Public   Mean wave period.
      !      THM       R.A.  Public   Mean wave direction.
      !      THS       R.A.  Public   Mean directional spread.
      !      FP0       R.A.  Public   Peak frequency.
      !      THP0      R.A.  Public   Peak direction.
      !      FP1       R.A.  Public   Wind sea peak frequency.
      !      THP1      R.A.  Public   Wind sea peak direction.
      !      DTDYN     R.A.  Public   Mean dynamic time step (raw).
      !      FCUT      R.A.  Public   Cut-off frequency for tail.
      !      ABA       R.A.  Public   Near-bottom rms wave ex. apmplitude.
      !      ABD       R.A.  Public   Corresponding direction.
      !      UBA       R.A.  Public   Near-bottom rms wave velocity.
      !      UBD       R.A.  Public   Corresponding direction.
      !      Sxx       R.A.  Public   Radiation stresses.
      !      DDDx      R.A.  Public   Spatial derivatives of the depth.
      !      DCxDx     R.A.  Public   Spatial dirivatives of the current.
      !
      ! Need to have the following loop from n = 1, isea not over 2d
      ! But need to map to 2d - since that is what is stored here?
      ! But is this what is used in the computation?
      ! do not call w3uwnd

!!$      do isea=1, nsea
!!$         uwind = x2w_w%rattr(index_x2w_sa_u,isea)         
!!$         vwind = x2w_w%rattr(index_x2w_sa_v,isea)
!!$         tbot  = x2w_w%rattr(index_x2w_sa_tbot,isea)
!!$         tocn  = x2w_w%rattr(index_x2w_so_t   ,isea)
!!$         ua(isea) = sqrt ( uwind**2 + vwind**2 )
!!$         if ( ua(isea) .gt. 1.e-7) then
!!$            ud(isea) = mod ( tpi+atan2(vwind,uwind) , tpi )
!!$         else
!!$            ud(isea) = 0.
!!$         end if
!!$         as(isea)   = tbot - tocn
!!$         u10 (isea) = max ( ua(isea) , 0.001 )
!!$         u10d(isea) = ud(isea)
!!$      end do
      U10    = 0.01
      U10D   = 0.
      UST    = 0.05
      USTDIR = 0.05
      !
      ! 7.b Run the wave model for the given interval
      !
      call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )

      hh = tod/3600
      mm = (tod - (hh * 3600))/60
      ss = tod - (hh*3600) - (mm*60) 

      timen(1) = ymd
      timen(2) = hh*10000 + mm*100 + ss

      call seq_timemgr_EClockGetData( EClock, prev_ymd=ymd, prev_tod=tod )

      hh = tod/3600
      mm = (tod - (hh * 3600))/60
      ss = tod - (hh*3600) - (mm*60) 

      time0(1) = ymd
      time0(2) = hh*10000 + mm*100 + ss

      time = time0

      write(6,*)'time0= ',time0
      write(6,*)'timen= ',timen
      call w3wave ( 1, timen )

!!$      copy ww3 data to coupling datatype
!!$      do isea=1, nsea
!!$         w2x_w%rattr(index_w2x_fw_taux,isea) = ??
!!$         w2x_w%rattr(index_w2x_fw_tauy,isea) = ??
!!$         w2x_w%rattr(index_w2x_sw_langnum,isea) = ??
!!$      enddo

      !
      ! TODO Put in gptl timer calls
      !
      ! Formats
      !
    END SUBROUTINE WAV_RUN_MCT
    
!=====================================================================
!=====================================================================
!=====================================================================

    SUBROUTINE WAV_FINAL_MCT
      ! do nothing now
    END SUBROUTINE WAV_FINAL_MCT

!=====================================================================
!=====================================================================
!=====================================================================

    subroutine wav_setgsmap_mct(mpi_comm, compid, gsmap)

      !/ ------------------------------------------------------------------- /
      !use w3gdatmd, only: nx, ny, nseal
      !use w3odatmd, only: naproc
      implicit none
      !/
      !/ ------------------------------------------------------------------- /
      !/ parameter list
      !/
      integer        , intent(in)  :: mpi_comm
      integer        , intent(in)  :: compid
      type(mct_gsmap), intent(out) :: gsmap
      !/
      !/ ------------------------------------------------------------------- /
      !/ local parameters
      !/
      integer, allocatable :: gindex(:)
      integer :: n,jsea,isea,ix,iy
      integer :: ier   
      ! -------------------------------------------------------------------- /
      
      allocate(gindex(nseal))
      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2) 
         gindex(jsea) = ix + (iy-1)*nx 
      end do
      call mct_gsmap_init( gsmap, gindex, mpi_comm, compid, nseal, nx*ny)
      deallocate(gindex)

    end subroutine wav_setgsmap_mct

!=====================================================================
!=====================================================================
!=====================================================================

    subroutine wav_domain_mct(lsize, gsmap, dom)

      integer        , intent(in)   :: lsize
      type(mct_gsmap), intent(in)   :: gsmap
      type(mct_ggrid), intent(inout):: dom  

      integer  :: n,i,ix,iy,isea,jsea   ! indices	
      real(r8) :: lon, lat, mask
      real(r8), pointer  :: data(:)     ! temporary
      integer , pointer  :: idata(:)    ! temporary
      real(r8), parameter:: radtodeg = 180.0_r8/shr_const_pi

      !
      ! initialize mct domain type
      ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
      ! note that in addition land carries around landfrac for the purposes of domain checking
      ! 
      call mct_ggrid_init( ggrid=dom, coordchars=trim(seq_flds_dom_coord), &
           otherchars=trim(seq_flds_dom_other), lsize=lsize )
      !
      ! allocate memory
      !
      allocate(data(lsize))
      !
      ! determine global gridpoint number attribute, globgridnum, which is set automatically by mct
      !
      call mct_gsMap_orderedPoints(gsmap, iaproc-1, idata)
      call mct_gGrid_importIattr(dom,'GlobGridNum',idata,lsize)
      !
      ! determine domain (numbering scheme is: west to east and south to north to south pole)
      ! initialize attribute vector with special value
      !
      data(:) = -9999.0_r8 
      call mct_ggrid_importrattr(dom,"lat"  ,data,lsize) 
      call mct_ggrid_importrattr(dom,"lon"  ,data,lsize) 
      call mct_ggrid_importrattr(dom,"area" ,data,lsize) 
      call mct_ggrid_importrattr(dom,"aream",data,lsize) 
      data(:) = 0.0_r8     
      call mct_ggrid_importrattr(dom,"mask" ,data,lsize) 
      !
      ! fill in correct values for domain components
      ! note aream will be filled in in the atm-lnd mapper
      ! sx, sy  real  i  grid increments (deg.).        
      !

      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2) 
         lon = x0 + real(ix-1)*sx 
         data(jsea) = lon
         !write(6,*)' jsea= ',jsea,' lon is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"lon",data,lsize) 

      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2) 
         lat = y0 + real(iy-1)*sy 
         data(jsea) = lat
         !write(6,*)' jsea= ',jsea,' lat is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"lat",data,lsize) 

      do jsea = 1,nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2) 
         lat = y0 + real(iy-1)*sy 
         !data(i) = !??? todo - calculate area
      end do
      !call mct_ggrid_importrattr(dom,"area",data,lsize) -todo fill this in

      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2) 
         mask = mapsta(iy,ix)
         data(jsea) = mask
         !write(6,*)' jsea= ',jsea,' mask is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"mask",data,lsize) 
      call mct_ggrid_importrattr(dom,"frac",data,lsize) 

      n = mct_aVect_lSize(dom%data)

      deallocate(data)
      deallocate(idata)

    end subroutine wav_domain_mct

  END MODULE WAV_COMP_MCT



