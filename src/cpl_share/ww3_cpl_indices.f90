module ww3_cpl_indices
  
  use seq_flds_mod
  use mct_mod

  implicit none

  SAVE
  public                               ! By default make data private

  integer :: index_x2w_Sa_u     
  integer :: index_x2w_Sa_v     
  integer :: index_x2w_Sa_tbot  
  integer :: index_x2w_So_t     

  integer :: index_x2a_Faxx_fco2_lnd   ! co2 flux from land   
  integer :: index_x2a_Faxx_fco2_ocn   ! co2 flux from ocean  
  integer :: index_x2a_Faxx_fdms_ocn   ! dms flux from ocean

contains

  subroutine ww3_cpl_indices_set( )

    type(mct_aVect) :: w2x      ! temporary
    type(mct_aVect) :: x2w      ! temporary

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2w, rList=seq_flds_x2w_fields, lsize=1)
    call mct_aVect_init(w2x, rList=seq_flds_w2x_fields, lsize=1)

    index_x2w_Sa_u     = mct_avect_indexra(x2w,'Sa_u')
    index_x2w_Sa_v     = mct_avect_indexra(x2w,'Sa_v')
    index_x2w_Sa_tbot  = mct_avect_indexra(x2w,'Sa_tbot')
    index_x2w_So_t     = mct_avect_indexra(x2w,'So_t')

    call mct_aVect_clean(x2w)
    call mct_aVect_clean(w2x)

  end subroutine ww3_cpl_indices_set

end module ww3_cpl_indices
