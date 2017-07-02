! This unit defines all specific routines for a eight-atomic molecule of ethane-type

module mol_c2h6
  use accuracy
  use moltype

  implicit none

  public ML_coordinate_transform_C2H6, ML_b0_C2H6, ML_symmetry_transformation_C2H6, ML_rotsymmetry_C2H6
  private

  integer(ik), parameter :: verbose = 3 ! Verbosity level

  !--------------------------
  !      ZMAT_4BETA_1TAU
  !--------------------------
  !
  !ZMAT
  !C   0  0  0  0  12.00000000
  !C   1  0  0  0  12.00000000
  !H   1  2  0  0   1.00782505
  !H   1  2  3  2   1.00782505
  !H   1  2  3 -2   1.00782505
  !H   2  1  3  2   1.00782505
  !H   2  1  5  2   1.00782505
  !H   2  1  5 -2   1.00782505
  !end
  ! 1  r_12        0         1.52579576
  ! 2  r_31        0         1.09074923
  ! 3  r_41        0         1.09074923
  ! 4  r_51        0         1.09074923
  ! 5  r_62        0         1.09074923
  ! 6  r_72        0         1.09074923
  ! 7  r_82        0         1.09074923
  ! 8   alpha_312   0         111.0 DEG
  ! 9   alpha_412   0         111.0 DEG
  ! 10  alpha_512   0         111.0 DEG
  ! 11  alpha_621   0         111.0 DEG
  ! 12  alpha_721   0         111.0 DEG
  ! 13  alpha_821   0         111.0 DEG
  ! 14  thet12      0         120.00 DEG
  ! 15  thet13      0         120.00 DEG
  ! 16  thet45      0         120.00 DEG
  ! 17  thet64      0         120.00 DEG
  ! 18  t14         0         180.00 DEG
  !end




  contains


  function ML_coordinate_transform_C2H6(src,ndst,direct) result (dst)
    !
    ! Transformtation from Z-matrix to TROVE coords
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct
    !
    real(ark),dimension(ndst) :: dst
    integer(ik) :: nsrc
    real(ark) :: tau4213,tau5124,tau6213,tau5126
    real(ark) :: tau1, tau2, b1_zmat, b2_zmat, tau1_zmat, dtau, db1, db2
    !

    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C2H6/start'
    !
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_coordinate_transform_C2H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_coordinate_transform_C2H6 error: bad coordinate type'
      !
    case('ZMAT_4BETA_1TAU','C2H6_4BETA_1TAU')
      ! ORDER CHANGED HERE AS WANT TBAR TO BE 18TH COORDINATE BUT NEEDS
      !  TO BE 16TH COORDINATE IN Z-MAT
      !
      if (direct) then ! transform from Z-matrix coords to TROVE coords

        dst(1:13) = src(1:13)-molec%local_eq(1:13)

        dst(14)  = (3.D0*src(14) - 2.0*pi)/(6.D0**0.5)
        dst(15)  = ( 2.0*src(15) - 2.0*pi + src(14))/(2.0**0.5)
        dst(16) = (3.0*src(18) - 2.0*pi)/(6.0**0.5)
        dst(17) = ( 2*src(17) - 2.0*pi + src(18))/(2.0**0.5)
        dst(18) = ( ( 3.D0*src(16) + src(14) - src(15) + src(17) - src(18))/3.D0 ) - molec%local_eq(16)
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:13) = src(1:13)+molec%local_eq(1:13)

        dst(14) = ((6.D0**0.5D0)*src(14) + 2.0*pi)/3.D0
        dst(15) = ( (2.D0**0.5D0)*src(15) - src(14) + 2.0*pi)/2.D0
        dst(18) = ((6.D0**0.5D0)*src(17) + 2.0_rk*pi)/3.D0
        dst(17) = ( (2.D0**0.5D0)*src(16) - dst(18) + 2.0*pi)/2.D0
        dst(16) = ( (3.D0*src(18) - dst(14) + dst(15) - dst(16) + dst(17) )/3.D0  ) + molec%local_eq(16)
        !
      endif
      !
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C2H6/end'
    !
  end function ML_coordinate_transform_C2H6



  subroutine ML_b0_C2H6(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
    !
    ! Reference structure
    !
    integer(ik),intent(in) :: Npoints, Natoms
    real(ark),intent(out) :: b0(Natoms,3,0:Npoints)
    real(ark),intent(inout),optional :: rho_i(0:Npoints)
    real(ark),intent(out),optional :: rho_ref
    real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
    !
    real(ark) :: a0(molec%Natoms,3),CM_shift,tau,alpha0
    integer(ik) :: i, n, iatom
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C2H6/start'
    !
    a0 = 0.0_ark
    !
      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = -molec%req(1)*0.5_ark
      !
      a0(2,1) = 0.0_ark
      a0(2,2) = 0.0_ark
      a0(2,3) = molec%req(1)*0.5_ark
      !
      a0(3,1) = molec%req(2)*SIN(PI - molec%alphaeq(1))
      a0(3,2) = 0.0
      a0(3,3) = -molec%req(2)*COS(PI - molec%alphaeq(1)) - 0.5d0*molec%req(1)
      !
      a0(4,1) = molec%req(3)*SIN(PI - molec%alphaeq(2))*COS(molec%taueq(1))
      a0(4,2) = -molec%req(3)*SIN(PI - molec%alphaeq(2))*SIN(molec%taueq(1)) 
      a0(4,3) = -molec%req(3)*cos(pi - molec%alphaeq(2))- 0.5*molec%req(1)
      !
      a0(5,1) = molec%req(4)*SIN(PI - molec%alphaeq(3))*COS(-molec%taueq(2))
      a0(5,2) = -molec%req(4)*SIN(PI - molec%alphaeq(3))*SIN(-molec%taueq(2))
      a0(5,3) = -molec%req(4)*cos(pi - molec%alphaeq(3))- 0.5*molec%req(1)
      !
      a0(6,1) = molec%req(5)*SIN(PI - molec%alphaeq(4))*COS(molec%taueq(3))
      a0(6,2) = molec%req(5)*SIN(PI - molec%alphaeq(4))*SIN(molec%taueq(3))
      a0(6,3) = molec%req(5)*cos(pi - molec%alphaeq(4)) + 0.5*molec%req(1)
      !
      a0(7,1) = molec%req(6)*SIN(PI - molec%alphaeq(5))*COS(molec%taueq(3)+ molec%taueq(4) )
      a0(7,2) = molec%req(6)*SIN(PI - molec%alphaeq(5))*SIN(molec%taueq(3) + molec%taueq(4) )
      a0(7,3) = molec%req(6)*cos(pi - molec%alphaeq(5)) + 0.5*molec%req(1) 
     !
      a0(8,1) = molec%req(7)*SIN(PI - molec%alphaeq(5))*COS(molec%taueq(3) - molec%taueq(5) )
      a0(8,2) = molec%req(7)*SIN(PI - molec%alphaeq(6))*SIN(molec%taueq(3) - molec%taueq(5) )
      a0(8,3) = molec%req(7)*cos(pi - molec%alphaeq(6)) + 0.5*molec%req(1)  
     !
    !
    
    do n=1, 3
      CM_shift = sum(a0(:,n)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
      a0(:,n) = a0(:,n) - CM_shift
    enddo
    !
    if (verbose>=3) then
      do iatom=1, Natoms
        write(out, '(1x,a,1x,3(1x,es16.8))') 'H', a0(iatom,1:3)
      enddo
    endif
    !
    b0(:,:,0) = a0(:,:)
    !
    if (Npoints/=0) then
    stop
    end if
       !
!        call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i))
        !
!        do n = 1,3
!          CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
!          b0(:,n,i) = b0(:,n,i) - CM_shift
!        enddo
        !
!        if (verbose>=5) then
!          write(out, '(1x,a,1x,i6,100(1x,es16.8))') 'b0', i, b0(:,:,i)
!        endif
        !
!      enddo
      !
!    endif
    rho_ref = 0
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C2H6/end'
    !
  end subroutine ML_b0_C2H6



  subroutine ML_symmetry_transformation_C2H6(ioper, nmodes, src, dst)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer(ik), intent(in)  ::  ioper
    integer(ik), intent(in)  ::  nmodes
    real(ark), intent(in)    ::  src(1:nmodes)
    real(ark), intent(out)   ::  dst(1:nmodes)
    real(ark) :: a,b
    !
    a = 0.5d0
    b = 0.5d0*sqrt(3.0d0)
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C2H6/start'
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_symmetry_transformation_C2H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_symmetry_transformation_C2H6 error: bad coordinate type'
      !
    case('ZMAT_4BETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop
        !
      case('D3D(M)')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H6 error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:18) = src(1:18)
          !
        case (2) ! !C3+/(123)(465)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(2)
          dst(4) = src(3)
          dst(5) = src(6)
          dst(6) = src(7)
          dst(7) = src(5)
          dst(8) = src(10)
          dst(9) = src(8)
          dst(10) = src(9)
          dst(11) = src(12)
          dst(12) = src(13)
          dst(13) = src(11)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) = -b*src(14) - a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) = -b*src(16) - a*src(17)
          dst(18) = src(18)
          !
        case (3) !C3-/(123)(465)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(4)
          dst(4) = src(2)
          dst(5) = src(7)
          dst(6) = src(5)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(10)
          dst(10) = src(8)
          dst(11) = src(13)
          dst(12) = src(11)
          dst(13) = src(12)
          dst(14) = -a*src(14) - b*src(15)
          dst(15) =  b*src(14) - a*src(15) 
          dst(16) = -a*src(16) - b*src(17)
          dst(17) =  b*src(16) - a*src(17)
          dst(18) = src(18)
          !
        case (4) ! C2/(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(5)
          dst(4) = src(6)
          dst(5) = src(3)
          dst(6) = src(4)
          dst(7) = src(2)
          dst(8) = src(13)
          dst(9) = src(11)
          dst(10) = src(12)
          dst(11) = src(9)
          dst(12) = src(10)
          dst(13) = src(8)
          dst(14) = src(16)
          dst(15) = -src(17) 
          dst(16) = src(14)
          dst(17) = -src(15)
          dst(18) = src(18)
          !
        case (5) !C2'/(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(6)
          dst(3) = src(7)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(2)
          dst(7) = src(3)
          dst(8) = src(12)
          dst(9) = src(13)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(8)
          dst(13) = src(9)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) =  b*src(16) + a*src(17) 
          dst(16) = -a*src(14) - b*src(16)
          dst(17) = -b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (6) ! C2''(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(6)
          dst(4) = src(7)
          dst(5) = src(2)
          dst(6) = src(3)
          dst(7) = src(4)
          dst(8) = src(11)
          dst(9) = src(12)
          dst(10) = src(13)
          dst(11) = src(8)
          dst(12) = src(9)
          dst(13) = src(10)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) =  b*src(16) + a*src(17) 
          dst(16) = -a*src(14) + b*src(16)
          dst(17) = b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (7) ! i/(14)(26)(35)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(7)
          dst(4) = src(6)
          dst(5) = src(2)
          dst(6) = src(4)
          dst(7) = src(3)
          dst(8) = src(11)
          dst(9) = src(13)
          dst(10) = src(12)
          dst(11) = src(8)
          dst(12) = src(10)
          dst(13) = src(9)
          dst(14) = src(16)
          dst(15) = src(17) 
          dst(16) = src(14)
          dst(17) = src(15)
          dst(18) = -src(18)
          !
        case (8) ! S6/(163425)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(6)
          dst(3) = src(5)
          dst(4) = src(7)
          dst(5) = src(4)
          dst(6) = src(3)
          dst(7) = src(2)
          dst(8) = src(12)
          dst(9) = src(11)
          dst(10) = src(13)
          dst(11) = src(10)
          dst(12) = src(9)
          dst(13) = src(8)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) =  -b*src(16) - a*src(17) 
          dst(16) = -a*src(14) + b*src(16)
          dst(17) = -b*src(14) - a*src(15)
          dst(18) = -src(18)
          !
        case (9) !S6'/(14)(26)(35)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(6)
          dst(4) = src(5)
          dst(5) = src(3)
          dst(6) = src(2)
          dst(7) = src(4)
          dst(8) = src(13)
          dst(9) = src(12)
          dst(10) = src(11)
          dst(11) = src(9)
          dst(12) = src(8)
          dst(13) = src(10)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) =  b*src(16) - a*src(17) 
          dst(16) = -a*src(14) - b*src(16)
          dst(17) =  b*src(14) - a*src(15)
          dst(18) = -src(18)
          !
        case (10) !sigmad/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(3)
          dst(4) = src(2)
          dst(5) = src(6)
          dst(6) = src(5)
          dst(7) = src(7)
          dst(8) = src(10)
          dst(9) = src(9)
          dst(10) = src(8)
          dst(11) = src(12)
          dst(12) = src(11)
          dst(13) = src(13)
          dst(14) = -a*src(14) - b*src(15)
          dst(15) =  -b*src(14) + a*src(15) 
          dst(16) = -a*src(16) - b*src(17)
          dst(17) =  -b*src(16) + a*src(17)
          dst(18) = -src(18)
          !
        case (11) !sigmad'/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(2)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(5)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(8)
          dst(9) = src(10)
          dst(10) = src(9)
          dst(11) = src(11)
          dst(12) = src(13)
          dst(13) = src(12)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) =  b*src(14) + a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) =  b*src(16) + a*src(17)
          dst(18) = -src(18)
          !
        case (12) !sigmad''/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(4)
          dst(5) = src(7)
          dst(6) = src(6)
          dst(7) = src(5)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(13)
          dst(12) = src(12)
          dst(13) = src(11)
          dst(14) = src(14)
          dst(15) = -src(15) 
          dst(16) =  src(16)
          dst(17) =  -src(17)
          dst(18) = -src(18)
        end select
        !
        !
      end select
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C2H6/end'
    !
  end subroutine ML_symmetry_transformation_C2H6



  subroutine ML_rotsymmetry_C2H6(J,K,tau,gamma,ideg)
    !
    ! Symmetry transformation rules of TROVE rotational functions
    !
    integer(ik),intent(in) :: J,K,tau
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C2H6/start'
    !
    select case(trim(molec%coords_transform))
      !
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_rotsymmetry_C2H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_rotsymmetry_C2H6 error: bad coordinate type'
      !
      !
    case('ZMAT_2BETA_1TAU','C2H6_2BETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C2H6 error: bad symmetry type'
        !
      case('D2H(M)')
        !
        gamma = 0
        ideg = 1
        !
       ! if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
       ! if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
       ! if (mod(K+2,2)/=0.and.tau==1) gamma = 7 ! B3g
       ! if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 7 ! B1g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 3 ! B3g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
!       CHANGING TO TRY TO GET TO WORK: BARRY
        !
      case('D2H')
        !
        gamma = 0
        ideg = 1
        !
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 7 ! B3g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 3 ! B1g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
        !
      end select
      !
      !
    case('R_ALPHA_4TAU')
      !
      select case(trim(molec%symmetry))
        !
     case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C2H6 error: bad symmetry type'
        !
      case('D2H','D2H(M)')
        !
        gamma = 0
        ideg = 1
        !
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 5 ! B2g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 7 ! B3g
        !
      case('D2H(S)')
        !
        gamma = 0
        ideg = 1
        !
        !if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        !if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
        !if (mod(K+2,2)/=0.and.tau==1) gamma = 6 ! B2u
        !if (mod(K+2,2)/=0.and.tau==0) gamma = 8 ! B3u
        !
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 7 ! B3g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
        !
      end select
      !
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C2H6/end'
    !
  end subroutine ML_rotsymmetry_C2H6


end module mol_c2h6
