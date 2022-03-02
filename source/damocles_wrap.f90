!----------------------------------------------------------------------------!
!  this subroutine is only called if the python mcmc wrapper is being used.  !
!  the logical lg_mcmc is set to true and input parameters passed from the   !
!  python routine.                                                           !
!----------------------------------------------------------------------------!
subroutine run_damocles_wrap()

    use globals
    use input
    use class_dust
    use initialise
    use vector_functions
    use driver

    implicit none

    

!f2py   intent(in) mcmc_v_min
!f2py   intent(in) mcmc_v_max
!f2py   intent(in) mcmc_rho_index
!f2py   intent(in) mcmc_mdust
!f2py   intent(in) mcmc_grain_size
!f2py   intent(in) mcmc_doublet_ratio
!f2py   intent(in) n

!f2py   intent(out) mcmc_mod
!f2py   depend(mcmc_v_max) mcmc_mod
!f2py   depend(n) mcmc_mod

   
    call run_damocles()

end subroutine
