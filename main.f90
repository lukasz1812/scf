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
write(*,*)"2 Unrestricted HF (opened shell systems)"
read(*,*)choose
select case (choose)
case (1)
   call restricted_HF(ng,nbf,nel,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore,nnrep,nat, number_Functions)
 case (2)
    call unrestricted_HF(ng,nbf,nel,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore,nnrep)
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
  real(wp)  ::start, finish, startscf, finishscf

  real(wp), allocatable :: proove(:,:)
   !Array Matrix for convergence check
  real(wp), allocatable :: array(:,:)

  !start measure of calculation time
  call cpu_time(start)

  nocc=2

  !Subroutine for calculation of 1 electron integrals
  call one_int(ng,nbf,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore)


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

    call sym_orthonormalizer(pack_smatrix(1:nbf),orthonormalizer)
orthonormalizer(1,1)=0.58700642812
orthonormalizer(1,2)=0.9541310722
orthonormalizer(2,1)=orthonormalizer(1,1)
orthonormalizer(1,2)=-orthonormalizer(2,2)

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
delta=1
i=1
allocate(array(nbf,nbf))
do while(delta>0.00001)

call iteration(i, nbf,ng,xyz,allcoeff,allalpha,ipmatrix,gmatrix,fMatrix,orthonormalizer,fcmatrix,nel,fpmatrix,vmatrix,tmatrix,fescf, nnrep, hcore,nocc)

!if (i<=26) then

    !call iteration(i, nbf,ng,xyz,allcoeff,allalpha,ipmatrix,gmatrix,fMatrix,orthonormalizer,fcmatrix,nel,fpmatrix,vmatrix,tmatrix,fescf, nnrep, hcore)
    array=fpmatrix-ipmatrix
    delta=sqrt(sum(array**2)/4)
    write(*,*) Delta, "Delta"
    ipmatrix=0
    ipmatrix=fpmatrix
    i=i+1
!  else
      !  write(*,*)"!!!! System could not be converged"
!  endif

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





end subroutine restricted_HF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine unrestricted_HF(ng,nbf,nel,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore, nnrep)

  integer:: ng, nbf,nel, i, Check, nelb, nela

  real(wp) :: delta, fescf, nnrep
  real(wp),allocatable:: gmatrix(:,:) ,  hcore(:,:)
  !Declaration of all matrices with β spin
  real(wp), allocatable :: ipbmatrix(:,:),fpbmatrix(:,:),fbmatrix(:,:),icbmatrix(:,:)
  !Declaration of all matrices with β spin
  real(wp), allocatable :: ipamatrix(:,:),fpamatrix(:,:),famatrix(:,:),icamatrix(:,:)
  real(wp), allocatable :: chrg(:), xyz(:,:), allalpha(:),allcoeff(:)
  real(WP), allocatable ::smatrix(:,:), tmatrix(:,:), vmatrix(:,:)
  real(wp), allocatable :: pack_smatrix(:), pack_tmatrix(:), pack_vmatrix(:)
  !>Symmetric orthonormalizer packed Matrix
  real(wp), allocatable :: orthonormalizer(:,:)
  real(wp), allocatable:: fcamatrix(:,:),fcbmatrix(:,:),gamatrix(:,:),gbmatrix(:,:)
  !>Calculation time variables
  real(wp)  ::start, finish, startscf, finishscf

  real(wp), allocatable :: proove(:,:)
   !Array Matrix for convergence check
  real(wp), allocatable :: array(:,:)

  !Number of electrons in orbital
  integer :: nocc
  nocc=1

  !start measure of calculation time
  call cpu_time(start)


  check=mod(nel,2)
  if (check==1) then

      !Calculating number of electrons with β spin
      nelb=(nel-1)/2
      nela=nelb+1
    else
      write(*,*)"You chosed a closed shell system"
      write(*,*) "The energy will be equal to restricted Hartree Fock Calculation"
      nelb=nel/2
      nela=nelb
  endif

      !start measure of calculation time
      call cpu_time(start)

      !Subroutine for calculation of 1 electron integrals
      call one_int(ng,nbf,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore)

      !Allocating memmory for some matrices
      allocate(fcamatrix(nbf,nbf),fcbmatrix(nbf,nbf))
      allocate(gamatrix(nbf,nbf),gbmatrix(nbf,nbf))
      allocate(gmatrix(nbf,nbf))
      allocate(ipbmatrix(nbf,nbf))
      allocate(fpbmatrix(nbf,nbf))
      allocate(fbmatrix(nbf,nbf), icbmatrix(nbf,nbf))
      !
      allocate(ipamatrix(nbf,nbf))
      allocate(fpamatrix(nbf,nbf))
      allocate(famatrix(nbf,nbf), icamatrix(nbf,nbf))
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

      call sym_orthonormalizer(pack_smatrix(1:nbf),orthonormalizer)

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
    call density(nocc,nbf, nela,icamatrix, ipamatrix)
    call density(nocc,nbf, nelb,icbmatrix, ipbmatrix)

    !========================== End of initial guess ==================================
    !=============================End of new Fock Matrix ===============================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SCF PROCEDURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call cpu_time(startscf)
    delta=1
    i=1
    allocate(array(nbf,nbf))
    do while(delta>0.00001)


call uiteration(i, nbf,ng,xyz,allcoeff,allalpha,ipamatrix, ipbmatrix,gamatrix,gbmatrix, famatrix, fbmatrix,orthonormalizer,fcamatrix,fcbmatrix,nela,nelb ,fpamatrix,fpbmatrix,vmatrix,tmatrix,fescf, nnrep, hcore, nocc)
        array=fpamatrix-ipamatrix
        delta=sqrt(sum(array**2)/4)
        write(*,*) Delta, "Delta"
        ipamatrix=0
        ipamatrix=fpamatrix
        ipbmatrix=fpbmatrix
        i=i+1
    !  else
          !  write(*,*)"!!!! System could not be converged"
    !  endif

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
    write(*,*)"Total calculation time [s]: ",finish-start


end subroutine unrestricted_HF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 subroutine one_int(ng,nbf,chrg,xyz,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore)

   real(wp), allocatable :: xyz(:,:),chrg(:),alphaa(:),alphab(:), coeffa(:),coeffb(:), allalpha(:), allcoeff(:)
   real(wp),allocatable :: smatrix(:,:), tmatrix(:,:), vmatrix(:,:), hcore(:,:)
   integer :: i, j, nbf, ng
   allocate(alphaa(ng), alphab(ng), coeffa(ng), coeffb(ng))
   !Loop over all slater functions
      do i=1,nbf
        do j=1,nbf

          !Basisfunctions Mapping
            coeffa=allcoeff(ng*(i-1)+1:ng*i)
            coeffb=allcoeff(ng*(j-1)+1:ng*j)
            alphaa=allalpha(ng*(i-1)+1:ng*i)
            alphab=allalpha(ng*(j-1)+1:ng*j)

            !One electron integral calculation

   !allocate(smatrix(nbf,nbf),tmatrix(nbf,nbf), vmatrix(nbf,nbf), hcore(nbf,nbf))
         call oneint(transpose(xyz),chrg,xyz(i,1:3),xyz(j,1:3),alphaa,alphab,coeffa,coeffb,smatrix(i,j),tmatrix(i,j), vmatrix(i,j))
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
  if (any(abs(matrix-transpose(matrix))>1.0e-14_wp)) then

    write(*,*) "Matrix is notsymmetric"

  else

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

      end if


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



      length=size(pack_smatrix)

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
      !orthonormalizer=matmul(eigen_smatrix, transpose(v))
      !orthonormalizer=matmul(v, orthonormalizer)


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

  call iHF(nbf,fmatrix,hcore,ipmatrix,fescf)

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
        allocate(alphab(ng))
        allocate(coeffa(ng))
        allocate(coeffb(ng))
        allocate(alphac(ng))
        allocate(alphad(ng))
        allocate(coeffc(ng))
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
                    coeffb=allcoeff(ng*(i-1)+1:ng*i)
                    alphaa=allalpha(ng*(j-1)+1:ng*j)
                    alphab=allalpha(ng*(i-1)+1:ng*i)
                    coeffc=allcoeff(ng*(l-1)+1:ng*l)
                    coeffd=allcoeff(ng*(k-1)+1:ng*k)
                    alphac=allalpha(ng*(l-1)+1:ng*l)
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

                        gmatrix(i,j)=gmatrix(i,j)+ipmatrix(l,k)*(tei(jilk)-0.5*tei(jlik))

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
    real(wp), allocatable :: fmatrix(:,:), icmatrix(:,:), feigen_func(:), v(:,:), pack_fmatrix(:), orthonormalizer(:,:)


      !Calculation of F’ Matrix
      icmatrix=matmul(transpose(orthonormalizer),fmatrix)
      icmatrix=matmul(icmatrix,orthonormalizer)
      write(*,*)
      write(*,*) "++++++++++++++++++++++  F' Matrix  ++++++++++++++++++++++++++++++"
      call  write_Matrix(icmatrix)

      !Packing f'matrix
      call packing(nbf,icmatrix,pack_fmatrix)

      !Calculating C'Matrix
      allocate(Feigen_func(nbf))
      allocate(v(nbf,nbf))
      call solve_spev(pack_fmatrix(1:nbf), feigen_func, v)
      v=-v
      write(*,*)
      write(*,*) "++++++++++++++++++++++  c' array  ++++++++++++++++++++++++++++++"
      call write_matrix(v)
      write(*,*)
      write(*,*) "++++++++++++++++++++++  E Matrix  ++++++++++++++++++++++++++++++"
      write(*,*)
      write(*,*)feigen_func



      !transforming c’ into C
      icmatrix=matmul(orthonormalizer,v)
      write(*,*)
      write(*,*) "++++++++++++++++++++++  C array  ++++++++++++++++++++++++++++++"
      call  write_Matrix(icmatrix)

end subroutine new_coefficients

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@q

subroutine iHF(nbf,fmatrix,hcore,ipmatrix,iescf)

    !Declaration of local variables
    real(wp), allocatable :: fMatrix(:,:),ipMatrix(:,:), hcore(:,:)
    integer :: i,j, nbf
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
subroutine density(nocc, nbf, nel,icmatrix, ipmatrix)



!Declarataion of local variables
integer :: nocc
integer :: nbf, i, j, k,nel
real(wp), allocatable :: icmatrix(:,:), ipmatrix(:,:)



iPmatrix=0

write(*,*) nel, "electrons"
 do i=1,nbf

   do j=1, nbf
     do k=1,nel

       iPmatrix(i,j)=iPmatrix(i,j)+nocc*(icmatrix(i,k)*icmatrix(j,k))



     end do

   end do

 end do

call write_Matrix(ipmatrix, "Density Matrix")


 end subroutine density
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


 subroutine Mulliken(nbf,ipmatrix,smatrix,nat,chrg, number_functions)

   integer :: u,v, nbf,i, nat, initial, final
   real(wp),allocatable ::ipmatrix(:,:),Smatrix(:,:), mlkn(:,:)
   real(wp), allocatable:: Mulliken_chrg(:), chrg(:)
   integer :: number_functions(:)
   real ::nel
   allocate(mlkn(nat,nat))

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


     initial=initial+chrg(i-1)
     final=final+chrg(i)
     u=initial

     do while (u<=final)

       write(*,*)"mulliken before", mulliken_chrg(i), initial
       mulliken_chrg(i)=mulliken_chrg(i)-mlkn(u,u)
      write(*,*)"mulliken after", mulliken_chrg(i)
       u=u+1

     end do

   end do

 write(*,*)mulliken_chrg, "Mulliken"

 end subroutine mulliken

 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 subroutine new_ugmatrix(nbf,ng,xyz,allcoeff,allalpha,ipamatrix,ipbmatrix, gmatrix)

   !Local variable
   integer :: i,j,k,l,m,n, ji, lk, jilk, nbf, ng,max, jl, ik, jlik
   real(wp), allocatable :: coeffa(:),coeffb(:),coeffc(:),coeffd(:)
   real(wp), allocatable :: alphaa(:),alphab(:),alphac(:),alphad(:)
   real(wp), allocatable :: tei(:), dens
   real(wp), allocatable :: gmatrix(:,:), ipamatrix(:,:),ipbmatrix(:,:), xyz(:,:),allalpha(:),allcoeff(:)      !

         !allocate Memory for the two electrons integral input
         allocate(alphaa(ng))
         allocate(alphab(ng))
         allocate(coeffa(ng))
         allocate(coeffb(ng))
         allocate(alphac(ng))
         allocate(alphad(ng))
         allocate(coeffc(ng))
         allocate(coeffd(ng))

         dens=0
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
                     coeffb=allcoeff(ng*(i-1)+1:ng*i)
                     alphaa=allalpha(ng*(j-1)+1:ng*j)
                     alphab=allalpha(ng*(i-1)+1:ng*i)
                     coeffc=allcoeff(ng*(l-1)+1:ng*l)
                     coeffd=allcoeff(ng*(k-1)+1:ng*k)
                     alphac=allalpha(ng*(l-1)+1:ng*l)
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
call write_matrix(ipaMatrix, "IPA")

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
   dens=ipamatrix(l,k)+ipbmatrix(l,k)


                         gmatrix(i,j)=gmatrix(i,j)+dens*(tei(jilk))-ipamatrix(l,k)*tei(jlik)
                         write(*,*) gmatrix(j,i)

     end do
   end do


 end do
 end do


 ! Printing of G-matrix
 write(*,*)
 write(*,*) "++++++++++++++++++++++  G array  ++++++++++++++++++++++++++++++"
 call write_Matrix(Gmatrix)





end subroutine new_ugmatrix
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 subroutine uiHF(nbf,famatrix,fbmatrix,hcore,ipamatrix, ipbmatrix,iescf)

     !Declaration of local variables
     real(wp), allocatable :: faMatrix(:,:),fbMatrix(:,:),ipaMatrix(:,:),ipbMatrix(:,:), hcore(:,:)
     integer :: i,j, nbf
     real(wp) :: iescf
     real(wp) :: HFaMatrix(nbf,nbf)
     real(wp) :: HFbMatrix(nbf,nbf)
     iescf=0

     HFaMatrix=hcore+famatrix
     HFbMatrix=hcore+fbmatrix
     call write_Matrix(ipamatrix, "before multiplikation")
     HFamatrix=matmul(ipamatrix,HFamatrix)
    HFbmatrix=matmul(ipbmatrix,HFbmatrix)



         do i=1,nbf

           do j=1,nbf

              iescf=iescf+0.5*HFamatrix(j,i)+0.5*HFbmatrix(j,i)


           end do

         end do
         write(*,*)
     write(*,*)"Electronic energy=", iESCF


 end subroutine uiHF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine uiteration(i, nbf,ng,xyz,allcoeff,allalpha,ipamatrix, ipbmatrix,gamatrix,gbmatrix, famatrix, fbmatrix,orthonormalizer,fcamatrix,fcbmatrix,nela,nelb ,fpamatrix,fpbmatrix,vmatrix,tmatrix,fescf, nnrep, hcore, nocc)

  integer :: nbf,ng,nela, nelb, i, nocc
  real(wp), allocatable :: xyz(:,:)
  real(wp), allocatable :: gamatrix(:,:),gbmatrix(:,:), ipamatrix(:,:),ipbmatrix(:,:), allalpha(:),allcoeff(:)
  real(wp), allocatable :: famatrix(:,:), fbmatrix(:,:), fcamatrix(:,:),fcbmatrix(:,:), orthonormalizer(:,:)
  real(wp), allocatable :: vMatrix(:,:),tMatrix(:,:),fpaMatrix(:,:), fpbmatrix(:,:), hcore(:,:)
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
  call new_ugmatrix(nbf,ng,xyz,allcoeff,allalpha,ipamatrix,ipbmatrix, gamatrix)
  call new_ugmatrix(nbf,ng,xyz,allcoeff,allalpha,ipbmatrix,ipamatrix, gbmatrix)

  call new_Fock(hcore,gamatrix,famatrix)
  call new_Fock(hcore,gbmatrix,fbmatrix)
  call new_coefficients(famatrix,orthonormalizer,nbf,fcamatrix)
  call new_coefficients(fbmatrix,orthonormalizer,nbf,fcbmatrix)
  call write_Matrix(ipaMatrix, "before energy")

 call uiHF(nbf,famatrix,fbmatrix,hcore,ipamatrix, ipbmatrix,fescf)



  call  density(nocc,nbf, nela,fcamatrix, fpamatrix)
    call  density(nocc,nbf, nelb,fcbmatrix, fpbmatrix)
    call write_Matrix(fpbmatrix,"Beta")


write(*,*) "Electronic energy", fescf


end subroutine uiteration
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine num_devxyz(xyz,step,nat, chrg,ng,nbf,nel, allalpha, allcoeff,nocc, vec_gradient )

!Declaration of global variables
integer :: step, nat,ng, nbf, nel, nocc
real (wp), intent(out) :: vec_gradient(nat,3)

!Declaration of local variables
Integer :: i, j, ndim
real(wp) :: nnrep, escf, Eforward, Ebackward
real(wp), allocatable :: xyz_gradient(:,:), smatrix(:,:), vmatrix(:,:), tmatrix(:,:), hcore(:,:), allalpha(:),allcoeff(:)
real(wp), allocatable :: icmatrix(:,:), ipmatrix(:,:), gmatrix(:,:), fmatrix(:,:), chrg(:), xyz(:,:)

!Set start variables
ndim=3
vec_gradient=0

!Allocating Memory
allocate(xyz_gradient(nat,ndim),smatrix(nbf,nbf),vmatrix(nbf,nbf),tmatrix(nbf,nbf), hcore(nbf,nbf))
allocate(icmatrix(nbf,nbf),ipmatrix(nbf,nbf),gmatrix(nbf,nbf))
allocate(fmatrix(nbf,nbf),chrg(nat),xyz(nat,3),allalpha(ng*nbf),allcoeff(ng*nbf))
 !Calculating change for each atom
  do i=1,nat
    !Calculating changefor each dimension
    do j=1,ndim





      !adding step
      xyz_gradient(i,j)=xyz(i,j)+step
      call Nuclei_Rep(nat, xyz_gradient, chrg, nnrep)
       call one_int(ng,nbf,chrg,xyz_gradient,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore)
       nel=nel/2
       call density(nocc, nbf, nel,icmatrix, ipmatrix)
       nel=nel*2
       call new_gmatrix(nbf,ng,xyz,allcoeff, allalpha,ipmatrix, gmatrix)
       call new_Fock(hcore,gmatrix,fmatrix)
       call iHF(nbf,fmatrix,hcore,ipmatrix,escf)
       Eforward=nnrep+escf

        !removing step
        xyz_gradient(i,j)=xyz(i,j)-step
        call Nuclei_Rep(nat, xyz_gradient, chrg, nnrep)
        call one_int(ng,nbf,chrg,xyz_gradient,allalpha,allcoeff,smatrix,tmatrix,vmatrix, hcore)
        nel=nel/2
        call density(nocc, nbf, nel,icmatrix, ipmatrix)
        nel=nel*2
        call new_gmatrix(nbf,ng,xyz,allcoeff, allalpha,ipmatrix, gmatrix)
        call new_Fock(hcore,gmatrix,fmatrix)
        call iHF(nbf,fmatrix,hcore,ipmatrix,escf)

        !set initial Atom Position
        xyz_gradient(i,j)=xyz(i,j)

        !Writing the gradient vector
        vec_gradient(i,j)=(Eforward-Ebackward)/(2*step)

    end do
  end do

  !Deallocate Memory
  deallocate(xyz_gradient,smatrix,vmatrix,tmatrix, hcore)
  deallocate(icmatrix,ipmatrix,gmatrix)
end subroutine num_devxyz

end module scf_main
