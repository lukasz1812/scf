!The Program was written during the Quantum Chemistry II internship
!at the University of Bonn
!Written by Łukasz Wantoch, Cologne 2020

module scf_main
    !> Include standard Fortran environment for IO
    use iso_fortran_env, only : output_unit, error_unit

    ! ------------------------------------------------------------------------
    !> library functions provided by your lab assistents:

    !> interface to LAPACK's double precision symmetric eigenvalue solver (dspev)
    !  examples:
    !  call solve_spev(mat, eigval, eigvec)
    use linear_algebra, only : solve_spev

    !> expansion of slater-functions into contracted gaussians,
    !  coefficients and primitive exponents are taken from R.F. Stewart, JCP, 1970
    !  example:
    !  call expand_slater(zeta, alpha, coeff)
    use slater, only : expand_slater

    !> calculates one-electron integrals and two-electron integrals over
    !  spherical gaussians (s-functions). One-electron quanities supported
    !  are overlap, kinetic energy and nuclear attraction integrals.
    !  Two-electron integrals are provided in chemist notation.
    !  examples:
    !  call oneint(xyz, chrg, r_a, r_b, alp, bet, ca, ca, s, t, v)
    !  call twoint(r_a, r_b, r_c, r_d, alp, bet, gam, del, ca, cb, cc, cd, g)
    use integrals, only : oneint, twoint

    !> prints a matrix quantity to screen
    !  examples:
    !  call write_vector(vec, name='vector')
    !  call write_matrix(mat, name='matrix')
    !  call write_matrix(mat, name='packed matrix')
    use print_matrix, only : write_vector, write_matrix

    !> other tools that may help you jump ahead with I/O-heavy tasks
    !  example:
    !  call read_line(input, line)
    use io_tools, only : read_line

    !> Always declare everything explicitly
    implicit none

    !> All subroutines within this module are not exported, except for scf_prog
    !  which is the entry point to your program
    private
    public :: scf_prog

    !> Selecting double precision real number
    integer, parameter :: wp = selected_real_kind(15)




contains


!> This is the entry point to your program, do not modify the dummy arguments
!  without adjusting the call in lib/prog.f90
subroutine scf_prog(input)

    !> Always declare everything explicitly
    implicit none

    !>Prooving repeating of Program
    integer :: quit

    !> IO unit bound to the input file___________"
    integer, intent(in) :: input

    !> System specific data
    !> Number of atoms
    integer :: nat

    !> Number of electrons
    integer :: nel


    !> Atom coordinates of the system, all distances in bohr
    real(wp), allocatable :: xyz(:,:)

    !> Nuclear charges
    real(wp), allocatable :: chrg(:)

    !>Exponents of primitive Gaussian used for Slater Approximation
    real(wp), allocatable :: alpha(:)
    real(wp), allocatable :: allalpha(:)

    !>Coefficients of primitive Gaussian used for Slater Approximation
    real(wp), allocatable :: coeff(:)
    real(wp), allocatable :: allcoeff(:)

    !> Number of basis functions
    integer :: nbf

    !>Number of Gaussian functions used for Slater Approximation
    integer :: ng

    !>number of slater functionss for each Atom
    integer, allocatable :: number_Functions(:)

    !> Slater exponents of basis functions
    real(wp),allocatable :: zeta(:)

    !>one electron integrals input
    real(wp), allocatable :: alphaa(:)
    real(wp), allocatable :: alphab(:)
    real(wp), allocatable :: coeffa(:)
    real(wp), allocatable ::coeffb(:)
    real(wp), allocatable :: alphac(:)
    real(wp), allocatable :: alphad(:)
    real(wp), allocatable :: coeffc(:)
    real(wp), allocatable ::coeffd(:)

    !>Overlap Matrix
    real(wp), allocatable :: smatrix(:,:)

    !>kinetic energy Matrix
    real(wp), allocatable :: tmatrix(:,:)
    real(wp), allocatable :: pack_tmatrix(:)

    !>Nuclear Attraction Matrix
    real(wp), allocatable :: vmatrix(:,:)
    real(wp), allocatable :: pack_vmatrix(:)


    !>Symmetric orthonormalizer packed Matrix
  !  real(wp), allocatable :: orthonormalizer(:,:)


    !>Fock packed Matrix
    real(wp), allocatable :: Fmatrix(:,:)
    real(wp), allocatable :: pack_fmatrix(:)

    !>Fock packed Matrix
    real(wp), allocatable :: hcore(:,:)




    !>Fock packed Matrix
    real(wp), allocatable :: Gmatrix(:,:)

    !>initial coefficients matrix
    real(wp), allocatable :: icmatrix(:,:)

    !>final coefficients matrix
    real(wp), allocatable :: fcmatrix(:,:)

    !>initial density matrix Matrix
    real(wp), allocatable :: ipmatrix(:,:)

    !>initial density matrix Matrix
    real(wp), allocatable :: fpmatrix(:,:)
    !>Coefficints paked matrix
    real(wp), allocatable :: pack_Cmatrix(:)


    real(wp), allocatable :: seigen_func(:)
    real(wp), allocatable :: v(:,:)

    integer :: length
    integer:: choose



    !> Nuclei repulsion energy
    real(wp) :: nnrep
    real(wp) :: nnrep_part

    !> Hartree-Fock energy
    real(wp) :: escf
    real(wp) :: iescf
    real(wp) :: fescf

    !>SCF quit paramter
    real(wp) :: delta




    !>Name of the input file
    character :: Filename*100

    !>counters
    integer:: i_Atom,i, j, k, l, m, n, lk, ji, jilk

        !Printing title
        write(*,*)"________________________________________________________________________________"
        Write(*,*)"|#######                                       ___________               #######|"
        Write(*,*)"|######                ||          ||        ||                           ######|"
        Write(*,*)"|#####                 ||          ||        ||                            #####|"
        Write(*,*)"|####                  ||__________||        ||__________                   ####|"
        Write(*,*)"|###=                  ||          ||        ||                             =###|"
        Write(*,*)"|##==                  ||          ||        ||                             ==##|"
        Write(*,*)"|#====                 ||          ||artree  ||ock                         ====#|"
        write(*,*)"_________________________________________________________________________________"



write(*,*)
write(*,*)"       ----------------------Reading an input file-----------------------"
write(*,*)

!Users Input
Write (*,*) "Give the name of input file (*.in) you want to calculate "
read (*,"(a)") filename
write(*,*)"_________________________________________________________________________"

!Reading Data from file
  call input_reader(Filename, nat, nel, nbf, chrg, xyz, number_Functions, zeta)


!Calculating nuclei repulsion Energy
  call Nuclei_Rep(nat, xyz, chrg, nnrep)


!!///////////////////////  Check of start variables ////////////////////////////

    write(*,*)
    Write(*,*) "Total Number of Atoms:    ", nat
    write(*,*) "Total Number of Electrons:", nel
    write(*,*) "Total Number of Basissets:", nbf
    write(*,*) "Nuclei repulsion Energy /[Hartree]", nnrep
    write(*,*)

    !Check of Matrices
    !Print of position Matrix
    Write(*,*)"+++++++++++++++++++++++++++ Position Matrix [Bohr] +++++++++++++++++++++++++++++++"
    write(*,*)
    write(*,*)"    #---------[x]------------------------[y]-----------------------[z]-------"

    !Print Matrix in table form
    i_Atom=1

    do while(i_Atom<=nat)
      write(*,'(i6)',advance='no') i_Atom
      write(*,*) xyz(i_Atom,1:3)
      write(*,*)
      i_Atom=i_Atom+1
    end do
    write(*,*)

    !Print charge Matrix
    Write(*,*)"++++++++++++++++++++++   Charge Matrix [el. Charge] ++++++++++++++++++++++++++++++"
    write(*,*)

    !Print Matrix in table form
      do i_Atom=1,nat
        write(*,'(i6)',advance='no') i_Atom
        write(*,"(f15.2)") chrg(i_Atom)
      end do

    write(*,*)

    !Print charge Matrix
    Write(*,*)"++++++++++++++++++++++   Exponents of ζ ++++++++++++++++++++++++++++++"
    write(*,*)

    !Print Matrix in table form
    i=1

    do i_Atom=1,nat
      write(*,*)
      write(*,'(i6)',advance='no') i_Atom

      do while(i<=(sum(number_functions(1:i_Atom))))
        write(*,'(f15.8)',advance='no') zeta(i)
        i=i+1
      end do

    end do
!================================= Check End ===========================================

    !Users input
    write(*,*)
      write (*,*) "Give the Value x of the Basis STO-xG"
      read(*,*) ng

!++++++++++++++++++++++++++ Slater expansion ++++++++++++++++++++++++++++++++++++++++++++

      !Allocating memory for Slater Expansion
      allocate(alpha(ng))
      allocate(coeff(ng))
      allocate(allalpha(ng*nbf))
      allocate(allcoeff(ng*nbf))

    !  allocate(hcore(nbf,nbf))




      !>Calculating the primitive gaussian functions
      do i=1,nbf
        call expand_slater(ng, zeta(i),alpha, coeff)
        allalpha(i*ng-(ng-1):i*ng)=alpha
        allcoeff(i*ng-(ng-1):i*ng)=coeff
      end do
!======================= End slater Expansion ===================================

!//////////////////////// Expamsion Check //////////////////////////////////////
      !Printing all gaussian Functions as a table
      write(*,*)
      write(*,*)"__________________________________________________________________"
      write(*,*)"     coefficients         |", "         Exponents"
      write(*,*)"__________________________________________________________________"

      do i=1, nbf*ng
        write(*,*) allalpha(i),"|", allcoeff(i)
      end do

        write(*,*)"__________________________________________________________________"
!======================= End  Expansion Check ===================================



!Allocate memory for one electron integral calculations input
allocate(alphaa(ng))
allocate(alphab(ng))
allocate(coeffa(ng))
allocate(coeffb(ng))
allocate(alphac(ng))
allocate(alphad(ng))
allocate(coeffc(ng))
allocate(coeffd(ng))

!Allocate memory for one electron integral calculations output
allocate(smatrix(nbf,nbf))
allocate(tmatrix(nbf,nbf))
allocate(vmatrix(nbf,nbf))

write(*,*)
write(*,*)"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
write(*,*)"|                       -Method Menu-                           |"
write(*,*)"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
write(*,*)
write(*,*)"1 Restricted HF (closed shell systems)"
write(*,*)"      (minimal basis sets STO-NG)"
write(*,*)"2 Unrestricted HF (opened shell systems)"
write(*,*)"      (minimal basis sets STO-NG)"
write(*,*)"3 Geometry optimization"
write(*,*)"              (using RHF)"

read(*,*)choose
select case (choose)
case (1)
1   call restricted_HF(ng,nbf,nel,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore,nnrep,nat, number_Functions)
 !case (2)
  !  call unrestricted_HF(ng,nbf,nel,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore,nnrep)
  !case (3)
!     call geo_opt(ng,nbf,xyz,nel,chrg,allalpha,allcoeff,nat,number_functions, nnrep)
end select


   !deallocate(coeffa,coeffb,alphaa, alphab)

end subroutine scf_prog


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine input_reader(Filename, nat, nel, nbf, chrg, xyz, number_Functions, zeta)

    !definition of local variables
    integer, parameter :: wp = selected_real_kind(15)
    character Filename*100
    integer :: n_max
    integer :: io
    integer i, i_Atom, i_Basis, i_Dim, n_Basis
    integer nat, nel, nbf, Dimensions
    real(wp):: read_vec(5)
    integer, allocatable :: number_Functions(:)
    real(wp), allocatable :: xyz(:,:), chrg(:), zeta(:)

    !open input file
    open(file=fileName,newunit=io)
    !Input of start Variables
    read(io,*) nat, nel, nbf

    !Setting Atomcounter with value 1
     i_Atom=1

    !Setting x,y,z Dimensions
    Dimensions=3

    !Set position xyzs Matrix for all atoms
    allocate(xyz(nat,Dimensions))

    !Set chrg matrix for all nuclei
    allocate(chrg(nat))


    !Set zeta matrix of all basis functions
    allocate(number_Functions(nat))
    allocate(zeta(nbf))



    !First Loop for all read xyzs
      do while (i_Atom<=nat)

      read (io,*) read_vec
        i_dim=1

        !Filling Position Matrix
        do while(i_Dim<4)
            xyz(i_Atom,i_dim)=read_vec(i_dim)
            i_Dim=i_dim+1
        end do

        !Filling chrg Matrix
        chrg(i_Atom)=read_vec(4)

        number_Functions(i_Atom)=read_vec(5)
        !Filling number of primitive Gaussian Matrix
        i_Basis=1

        !skipping all Gaussian zeta
          do while (i_Basis<=read_vec(5))
            read(io,*)
            i_basis=i_basis+1
          end do

          i_Atom=i_Atom+1
        end do
        close(io)
!-------------------End of first reading------------------------------


!-----------------Start of second reading---------------------------------
      open(file=Filename, newunit=io)
      !Skipping first textline
      read(io,*)

      !setting all counters
      n_max=0
      i_Basis=1
      i_Atom=1

      do while(i_Atom<=nat)

        !Reading 5 el.lines
        read(io,*) read_vec
        n_max=n_Max+Number_Functions(i_Atom)
        !Reading all zeta of prim. Gaussian
        do while(i_basis<=n_Max)
        read(io,*)zeta(i_basis)
        i_basis=i_Basis+1

      end do
      i_Atom=i_Atom+1

    end do


end subroutine input_reader

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine Nuclei_Rep(nat, xyz, chrg, nnrep)

!definition of local variables
  integer, parameter :: wp = selected_real_kind(15)
  integer :: i, j
  integer nat
  intrinsic :: Sum
  real(wp):: read_vec(5), length(3)
  real(wp) distance, nnrep, nnrep_part
  real(wp), allocatable :: xyz(:,:), chrg(:)


  nnrep=0
  !Case of 1 Atom
  if (nat==1) then
    nnrep=0
  else
    !Counting over all Nuclei pairs
  i=1
     do while(i<nat)

       j=i+1
        do while (j<=nat)

          !Calculating all distances between the nuclei paires
          length=xyz(i,1:3)-xyz(j,1:3)
          distance=sqrt(sum(length**2))
          !Calculating nuclei repulsion between 2 Atoms
          nnrep_part=chrg(i)*chrg(j)/(distance)

          !Adding all Repulsion energies
          nnrep=nnrep+nnrep_part
          j=j+1
        end do
        i=i+1
      end do
    end if



end subroutine nuclei_Rep

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine restricted_HF(ng,nbf,nel,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore, nnrep,nat, number_Functions)
  integer:: ng, nbf,nel, i, nocc,nat, number_functions(:)
  !Number of electrons in orbital


  real(wp) :: delta, fescf, nnrep
  real(wp),allocatable:: gmatrix(:,:) ,ipmatrix(:,:), fpmatrix(:,:),fmatrix(:,:), hcore(:,:)
  real(wp), allocatable :: chrg(:), xyz(:,:), allalpha(:),allcoeff(:)
  real(WP), allocatable ::smatrix(:,:), tmatrix(:,:), vmatrix(:,:)
  real(wp), allocatable :: pack_smatrix(:), pack_tmatrix(:), pack_vmatrix(:)
  !>Symmetric orthonormalizer packed Matrix
  real(wp), allocatable :: orthonormalizer(:,:)
  real(wp), allocatable:: fcmatrix(:,:)
  !>Calculation time variables
  real  ::start, finish, startscf, finishscf

  real(wp), allocatable :: proove(:,:)
   !Array Matrix for convergence check
  real(wp), allocatable :: array(:,:)

  !start measure of calculation time
  call cpu_time(start)

  nocc=2

  !Subroutine for calculation of 1 electron integrals
  call one_int(nel,ng,nbf,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore)


 !Allocating memmory for some matrices
  allocate(gmatrix(nbf,nbf))
  allocate(ipmatrix(nbf,nbf))
  allocate(fpmatrix(nbf,nbf))
  allocate(fmatrix(nbf,nbf), fcmatrix(nbf,nbf))
  allocate(hcore(nbf,nbf))
  allocate(orthonormalizer(nbf,nbf))

  !Calculating H_Core
  hcore=tmatrix+vmatrix

!////////////////////////   Integrals Check    //////////////////////////////////////

!Print of all one electron matrices
  write(*,*)
  write(*,*) "++++++++++++++++++++++  Overlap Matrix  ++++++++++++++++++++++++++++++"
  call write_matrix(smatrix)
  write(*,*)
  write(*,*) "++++++++++++++++++++++  Kinetic Matrix  ++++++++++++++++++++++++++++++"
  call write_matrix(tmatrix)
  write(*,*)
  write(*,*) "+++++++++++++++++++  Nuc. Attraction Matrix  ++++++++++++++++++++++++"
  call write_matrix(vmatrix)
  write(*,*)
  write(*,*) "++++++++++++++++++++++  Hcore Matrix  ++++++++++++++++++++++++++++++"
  call write_matrix(hcore)
!======================= End one electron integrals check ==========================================


!+++++++++++++++++++ Symmetry check and packin of all symmetric matrices +++++++++++++++++++++

  call packing(nbf,smatrix,pack_smatrix)

  call packing(nbf,tmatrix,pack_tmatrix)

  call packing(nbf,vmatrix,pack_vmatrix)

!========================== end of packing =========================================


!+++++++++++++++++++++++ Calculation of symmetric orthonormalizer +++++++++++++++++++++++

    call sym_orthonormalizer(pack_smatrix,orthonormalizer)


!========================== End of sym. orth. calculation ==================================



!/////////////////////////////// Print of Orthonormalizer ///////////////////////////////
write(*,*)
write(*,*) "++++++++++++++++++++++  Sym. Orthonormalizer  ++++++++++++++++++++++++++++++"
call write_Matrix(orthonormalizer)

allocate(proove(nbf,nbf))
proove=matmul(smatrix,orthonormalizer)
proove=matmul(transpose(orthonormalizer),proove)
write(*,*)
write(*,*) "++++++++++++++++++ Proove of  Sym. Orthonormalizer  ++++++++++++++++++++++++++++++"
call write_matrix(proove)
deallocate(proove)
!================================ End of print ==================================================



!+++++++++++++++++++++++++ Initial guess of density Matrix ++++++++++++++++++++++++++
ipmatrix=0
gmatrix=0
write(*,*)
write(*,*) "++++++++++++++++++++++  Density Matrix  ++++++++++++++++++++++++++++++"
call write_matrix(ipmatrix)

!========================== End of initial guess ==================================



!+++++++++++++++++++++++++++ Initial Fock Matrix +++++++++++++++++++++++++++
fmatrix=hcore+gmatrix


!=============================End of new Fock Matrix ===============================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SCF PROCEDURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cpu_time(startscf)
delta=16
i=1
allocate(array(nbf,nbf))
do while(delta>0.00001)


if (i<=26) then

   call iteration(i, nbf,ng,xyz,allcoeff,allalpha,ipmatrix,gmatrix,fMatrix,orthonormalizer,fcmatrix,nel,fpmatrix,vmatrix,tmatrix,fescf, nnrep, hcore, nocc)
    array=fpmatrix-ipmatrix
    delta=sqrt(sum(array**2)/4)
    write(*,*) Delta, "Delta"
    ipmatrix=0
    ipmatrix=fpmatrix
    i=i+1
 else
        write(*,*)"!!!! System could not be converged !!!!"
        exit
  endif

end do
write(*,*)
write(*,*)"#################################################################################"
write(*,*)"#                         -- End of SCF-Procedure --                            #"
write(*,*)"#################################################################################"

call cpu_time(finishscf)
deallocate(array)
fescf=fescf+nnrep
call cpu_time(finish)
write(*,*)
write(*,*)"Converged after [s]: ",finishscf-startscf
write(*,*)
Write(*,*) "final SCF energy",fescf
write(*,*)
 call Mulliken(nbf,ipmatrix,smatrix,nat,chrg, number_functions)
write(*,*)"Total calculation time [s]: ",finish-start



deallocate(gmatrix,ipmatrix,fpmatrix,fmatrix,fcmatrix,hcore, orthonormalizer)

end subroutine restricted_HF


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 subroutine one_int(nel,ng,nbf,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore)

   real(wp), allocatable :: xyz(:,:),chrg(:), allalpha(:), allcoeff(:)
   real(wp),allocatable :: smatrix(:,:), tmatrix(:,:), vmatrix(:,:), hcore(:,:)
   integer :: i, j, nbf, ng, nel
   !Loop over all slater functions
      do i=1,nbf
        do j=1,nbf

            !One electron integral calculation
          call oneint(transpose(xyz),chrg,xyz(i,1:3),xyz(j,1:3),allalpha(ng*(i-1)+1:ng*i),allalpha(ng*(j-1)+1:ng*j),allcoeff(ng*(i-1)+1:ng*i),allcoeff(ng*(j-1)+1:ng*j),smatrix(i,j),tmatrix(i,j), vmatrix(i,j))

       end do
     end do




 end subroutine one_int

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine packing(ndim,matrix,pack_smatrix)
  !definition of variables
  integer, parameter :: wp = selected_real_kind(15)
  integer, intent(in) :: ndim
  integer :: i,j,k,  check
  real(wp), allocatable,intent(in) :: matrix (:,:)
  real(wp), allocatable :: pack_smatrix(:)

  allocate(pack_smatrix(ndim*(1+ndim)/2))


!  Symmetry Check

      j=1
      k=1

        do while(j<=ndim)

          i=1

          do while(i<=j)!loops over all matrix elements

            pack_sMatrix(k)=matrix(j,i)
            k=k+1
            i=i+1

          end do

          j=j+1

        end do



end subroutine packing


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine new_Fock(hcore,gmatrix,fmatrix)

  !Declaration of local variables
  real(wp) :: hcore(:,:), gmatrix(:,:)
  real(wp) :: fmatrix(:,:)

  !calculaton of new Fock Matrix
  fmatrix=0
  fmatrix=hcore+gmatrix

  !Printing the Fock Matrix
  write(*,*)
  write(*,*) "++++++++++++++++++++++  Fock Matrix  ++++++++++++++++++++++++++++++"
  call write_matrix(fmatrix)

end subroutine new_Fock

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine sym_orthonormalizer(pack_smatrix,orthonormalizer)

    real(wp), intent(inout) :: pack_smatrix(:)
    real(wp), allocatable :: orthonormalizer(:,:)
    real(wp), allocatable :: Seigen_func(:)
    real(wp), allocatable :: v(:,:)
    real(wp), allocatable :: eigen_smatrix(:,:)
    real(wp), allocatable :: eigen_smatrixinv(:)
    integer :: i, length, j



      length=nint(sqrt(8*real(size(pack_smatrix), wp)+1)-1)/2

      allocate(Seigen_func(length))
      allocate(v(length,length))
      allocate(eigen_smatrixinv(length))
      allocate(eigen_smatrix(length, length))
      call solve_spev(pack_smatrix, Seigen_func, v, i)
      if (i /= 0) then
        error stop "Could not diagonlize matrix"
      endif


        eigen_smatrixinv=(1/sqrt(Seigen_func))



        do i=1,Length

          do j=1, Length

              if (i==j) then
               eigen_smatrix(i,j)=eigen_smatrixinv(i)
            else
              eigen_smatrix(i,j)=0
            end if

            end do

          end do
          orthonormalizer=matmul(v, eigen_smatrix)
          orthonormalizer=matmul(orthonormalizer,transpose(v))


end subroutine sym_orthonormalizer
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


subroutine iteration(i, nbf,ng,xyz,allcoeff,allalpha,ipmatrix,gmatrix,fMatrix,orthonormalizer,fcmatrix,nel,fpmatrix,vmatrix,tmatrix,fescf, nnrep, hcore, nocc)

  integer :: nbf,ng,nel, i, nocc
  real(wp), allocatable :: xyz(:,:)
  real(wp), allocatable :: gmatrix(:,:), ipmatrix(:,:),allalpha(:),allcoeff(:)
  real(wp), allocatable :: fmatrix(:,:), fcmatrix(:,:), orthonormalizer(:,:)
  real(wp), allocatable :: vMatrix(:,:),tMatrix(:,:),fpMatrix(:,:), hcore(:,:)
  real(wp) :: fescf, nnrep
  fescf=0
        write(*,*)
        write(*,*)
        write(*,*)"=================================================================="
        write(*,*)
        write(*,*)"                   Iteration step",i
        write(*,*)" ================================================================="
        write(*,*)
        write(*,*)
  call new_gmatrix(nbf,ng,xyz,allcoeff, allalpha,ipmatrix, gmatrix)

  call new_Fock(hcore,gmatrix,fmatrix)

  call iHF(nel, nbf,fmatrix,hcore,ipmatrix,fescf)

  call new_coefficients(fmatrix,orthonormalizer,nbf,fcmatrix)
  nel=nel/2
  call  density(nocc,nbf, nel,fcmatrix, fpmatrix)
  nel=nel*2


write(*,*) "Electronic energy", fescf


end subroutine iteration

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine new_gmatrix(nbf,ng,xyz,allcoeff, allalpha,ipmatrix, gmatrix)

  !Local variable
  integer :: i,j,k,l,m,n, ji, lk, jilk, nbf, ng,max, jl, ik, jlik
  real(wp), allocatable :: coeffa(:),coeffb(:),coeffc(:),coeffd(:)
  real(wp), allocatable :: alphaa(:),alphab(:),alphac(:),alphad(:)
  real(wp), allocatable :: tei(:)
  real(wp), allocatable :: gmatrix(:,:), ipmatrix(:,:), xyz(:,:),allalpha(:),allcoeff(:)      !

        !allocate Memory for the two electrons integral input
        allocate(alphaa(ng))
        allocate(coeffa(ng))
        allocate(alphab(ng))
        allocate(coeffb(ng))
        allocate(alphac(ng))
        allocate(coeffc(ng))
        allocate(alphad(ng))
        allocate(coeffd(ng))

        ji=(nbf*(nbf-1)/2+nbf)
        jilk=ji*(ji-1)/2+ji
        allocate(tei(jilk))




!counters set
      i=1
      j=1
      k=1
      l=1
      m=0
      n=0

      gmatrix=0

  do while(n<=m)

      do while(l<=nbf)

        do while(k<=l)

          i=1


            do while(i<=nbf)

                j=1

                  do while(j<=i)
                    if(l>k) then
                    lk=l*(l-1)/2+k
                  ELSE
                    lk=k*(k-1)/2+l
                  endif

                  if(j>i)then
                    ji=j*(j-1)/2+i
                  else
                    ji=i*(i-1)/2+j
                  endif
                  if (ji>lk) then
                    jilk=ji*(ji-1)/2+lk
                  ELSE
                    jilk=lk*(lk-1)/2+ji
                  end if


                    !two electron integrals calculation
                    coeffa=allcoeff(ng*(j-1)+1:ng*j)
                    alphaa=allalpha(ng*(j-1)+1:ng*j)
                    coeffb=allcoeff(ng*(i-1)+1:ng*i)
                    alphab=allalpha(ng*(i-1)+1:ng*i)
                    coeffc=allcoeff(ng*(l-1)+1:ng*l)
                    alphac=allalpha(ng*(l-1)+1:ng*l)
                    coeffd=allcoeff(ng*(k-1)+1:ng*k)
                    alphad=allalpha(ng*(k-1)+1:ng*k)

                    !Two electron integrals subroutine
                    call twoint(xyz(j,1:3), xyz(i,1:3), xyz(l,1:3), xyz(k,1:3), alphaa, alphab, alphac, alphad, coeffa, coeffb, coeffc, coeffd, tei(jilk))



                    j=j+1
                    m=i*j
                    n=l*k

                  end do
                  i=i+1
                end do
                k=k+1
              end do
              k=1

              l=l+1

            end do

          n=n+1
        end do

call write_Matrix(tei)


do j=1,nbf
do i=1,nbf
  do k=1,nbf
    do l=1, nbf

      if(l>k) then
      lk=l*(l-1)/2+k
    ELSE
      lk=k*(k-1)/2+l
    endif

    if(j>i)then
      ji=j*(j-1)/2+i
    else
      ji=i*(i-1)/2+j
    endif
    if (ji>lk) then
      jilk=ji*(ji-1)/2+lk
    ELSE
      jilk=lk*(lk-1)/2+ji
    end if

    if(i>k) then
    ik=i*(i-1)/2+k
  ELSE
    ik=k*(k-1)/2+i
  endif

  if(j>l)then
    jl=j*(j-1)/2+l
  else
    jl=l*(l-1)/2+j
  endif
  if (jl>ik) then
    jlik=jl*(jl-1)/2+ik
  ELSE
    jlik=ik*(ik-1)/2+jl
  end if

                        gmatrix(j,i)=gmatrix(j,i)+ipmatrix(l,k)*(tei(jilk)-0.5*tei(jlik))

    end do
  end do

end do
end do


! Printing of G-matrix
write(*,*)
write(*,*) "++++++++++++++++++++++  G array  ++++++++++++++++++++++++++++++"
call write_Matrix(Gmatrix)





end subroutine new_gmatrix
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@q

subroutine new_coefficients(fmatrix,orthonormalizer,nbf,icmatrix)

    integer :: nbf
    real(wp), allocatable :: fmatrix(:,:), icmatrix(:,:), feigen_val(:), v(:,:), pack_fmatrix(:), orthonormalizer(:,:)


      !Calculation of F’ Matrix
      icmatrix=matmul(transpose(orthonormalizer),fmatrix)
      icmatrix=matmul(icmatrix,orthonormalizer)
      write(*,*)
      write(*,*) "++++++++++++++++++++++  F' Matrix  ++++++++++++++++++++++++++++++"
      call  write_Matrix(icmatrix)

      !Packing f'matrix
      call packing(nbf,icmatrix,pack_fmatrix)

      !Calculating C'Matrix
      allocate(Feigen_val(nbf))
      allocate(v(nbf,nbf))
      call solve_spev(pack_fmatrix, feigen_val, v)
      write(*,*)
      write(*,*) "++++++++++++++++++++++  c' array  ++++++++++++++++++++++++++++++"
      call write_matrix(v)
      write(*,*)
      write(*,*) "++++++++++++++++++++++  E Matrix  ++++++++++++++++++++++++++++++"
      write(*,*)
      write(*,*)feigen_val



      !transforming c’ into C
      icmatrix=matmul(orthonormalizer,v)
      write(*,*)
      write(*,*) "++++++++++++++++++++++  C array  ++++++++++++++++++++++++++++++"
      call  write_Matrix(icmatrix)

end subroutine new_coefficients

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@q

subroutine iHF(nel,nbf,fmatrix,hcore,ipmatrix,iescf)

    !Declaration of local variables
    real(wp), allocatable :: fMatrix(:,:),ipMatrix(:,:), hcore(:,:)
    integer :: i,j, nbf,nel
    real(wp) :: iescf
    real(wp) :: HFMatrix(nbf,nbf)
    iescf=0

    HFMatrix=hcore+fmatrix

    HFmatrix=matmul(HFmatrix,ipmatrix)

        do i=1,nbf

          do j=1,nbf

            if(i==j) then

              Iescf=iescf+0.5*HFmatrix(i,j)
            endif

          end do

        end do
        write(*,*)
    write(*,*)"Electronic energy=", IESCF


end subroutine iHF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine density( nocc, nbf, nel,icmatrix, ipmatrix)



!Declarataion of local variables
integer :: nocc, nel
integer :: nbf, i, j, k
real(wp), allocatable :: icmatrix(:,:), ipmatrix(:,:), occ_M(:,:)
allocate(occ_M(nbf,nbf))


occ_m=0
iPmatrix=0
do i=1,nel
  do j=1, nel

    if(i==j) then
      occ_m(i,j)=nocc

    end if

  end do
end do
call write_Matrix(occ_m, "OCcupation Matrix")
ipmatrix=matmul(occ_M,transpose(icmatrix))
ipmatrix=(matmul(icmatrix,ipmatrix))
call write_Matrix(ipmatrix, "Density Matrix")


 end subroutine density
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


 subroutine Mulliken(nbf,ipmatrix,smatrix,nat,chrg, number_functions)

   integer :: u,v, nbf,i, nat, initial, final
   real(wp),allocatable ::ipmatrix(:,:),Smatrix(:,:), mlkn(:,:)
   real(wp), allocatable:: Mulliken_chrg(:), chrg(:)
   integer :: number_functions(:)
   real ::nel
   allocate(mlkn(nbf,nbf))

   mlkn=matmul(ipmatrix,smatrix)


   nel=0
  !////////////////////////////////Check if everythings orrectly
   do u=1, nbf
     do v=1,nbf
       if(u==v) then

       nel=nel+Mlkn(u,v)
       end if
     end do
   end do

   write (*,*) "number electrons:", nel
   initial=1
   final=0
   Mulliken_chrg=chrg


  !Mulliken Population Calculation
  do i=1,nat

     initial=initial+number_Functions(i-1)

     final=final+number_Functions(i)

     u=initial

     do while (u<=final)

mlkn(u,u)=sngl(mlkn(u,u))
mulliken_chrg(i)=sngl(mulliken_chrg(i))
       mulliken_chrg(i)=mulliken_chrg(i)-mlkn(u,u)
       u=u+1

     end do

   end do
   write(*,*)
   write(*,*)"      +++  Mulliken Population  +++"
do i=1, Nat
  write(*,*)"Atom", i
  write(*,*)"Charge", mulliken_chrg(i)
end do

 end subroutine mulliken

 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module scf_main
