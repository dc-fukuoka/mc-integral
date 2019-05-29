! integral gaussian function exp(-x^2)
! exact solution is pi^(1/2)

module mc_module
  use omp_lib
  implicit none
  integer,parameter :: dp = kind(1.0d0)

  interface mc
     module procedure init0 ! no argument
     module procedure init1 ! one argument
  end interface mc
  
  type mc
     integer :: npts
     real(dp) :: ends
     real(dp),allocatable,dimension(:) :: array
     real(dp) :: ans
   contains
     procedure,nopass :: init0, init1
     generic :: init => init0, init1
     procedure :: integrate
     procedure :: show_ans
     final :: fini
  end type mc

contains
  function init0() result(this)
    implicit none
    type(mc) :: this
    integer :: npts

    npts = 1024
    
    this%npts = npts
    this%ends = 1.0d5
    this%ans  = 0.0d0
    allocate(this%array(this%npts))
  end function init0

  function init1(npts) result(this)
    implicit none
    integer,intent(in) :: npts
    type(mc) :: this

    this%npts = npts
    this%ends = 1.0d5
    this%ans  = 0.0d0
    allocate(this%array(this%npts))
  end function init1

  subroutine fini(this)
    implicit none
    type(mc),intent(inout) :: this

    if (allocated(this%array)) deallocate(this%array)
  end subroutine fini

  pure function gauss(x) result(res)
    implicit none
    real(dp),intent(in) :: x
    real(dp) :: res

    res = exp(-1.0d0*x**2)
  end function gauss

  subroutine gen_rand(len, a, val_min, val_max, seed)
    use mkl_vsl_type
    use mkl_vsl
    implicit none
    
    integer,intent(in) :: len
    real(dp),dimension(len),intent(out) :: a
    real(dp),intent(in) :: val_min, val_max
    integer,intent(in) :: seed
    integer::ierr
    integer::brng, method
    type(vsl_stream_state)::stream

    brng   = vsl_brng_mt19937
    method = vsl_rng_method_uniform_std_accurate

    ierr = vslnewstream(stream, brng, seed)
    ierr = vdrnguniform(method, stream, len, a, val_min, val_max)
    ierr = vsldeletestream(stream)
  end subroutine gen_rand

  subroutine integrate(this)
    implicit none
    class(mc),intent(inout) :: this
    real(dp) :: val_min, val_max
    real(dp) :: ans
    integer :: seed
    integer :: i

    val_min = -this%ends
    val_max = this%ends
    seed    = 7777

    call gen_rand(this%npts, this%array, val_min, val_max, seed)

    ans = 0.0d0
    !$omp parallel do reduction(+:ans)
    do i = 1, this%npts
       ans = ans + (gauss(this%array(i)))/this%npts
    end do
    ans = 2*this%ends*ans
    this%ans = ans
  end subroutine integrate

  subroutine show_ans(this)
    implicit none
    class(mc),intent(inout) :: this

    write(6, *) "Monte Carlo integration:", this%ans
    write(6, *) "Exact solution:         ", sqrt(4.0d0*atan(1.0d0))
  end subroutine show_ans
end module mc_module

program main
  use mc_module
  implicit none
  
  type(mc) :: a
  character(len=32) :: argv1
  integer :: npts

  npts = 1024*1024*512
  if (command_argument_count() .ge. 1) then
     call get_command_argument(1, argv1)
     read(argv1, *) npts
  end if
  write(6, *) "npts:", npts

  a = mc(npts)
  write(6, '(2(a, 1pe14.5))') "integrate gaussian function from", -a%ends, " to ", a%ends
  call a%integrate
  call a%show_ans
  
  stop
end program main
