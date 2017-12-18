! ============================================================================
!
! ***************************************************************************
!
Module gauss_data
    !   ------------------
    !
    implicit none
    save
    !     
    integer:: iset
    real*8:: gset
    !
    data iset /0/
    !
end Module gauss_data
!       
! ===========================================================================

! =====================================================================================
! A C-program for MT19937, with initialization improved 2002/1/26.
! Coded by Takuji Nishimura and Makoto Matsumoto.

! Code converted to Fortran 95 by Josi Rui Faustino de Sousa
! Date: 2002-02-01

! Before using, initialize the state by using init_genrand(seed)
! or init_by_array(init_key, key_length).

! This library is free software.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

! Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura.
! Any feedback is very welcome.
! http://www.math.keio.ac.jp/matumoto/emt.html
! email: matumoto@math.keio.ac.jp

module mt19937

    implicit none

    intrinsic :: bit_size

    private
    public  :: init_genrand, init_by_array
    public  :: genrand_int32, genrand_int31
    public  :: genrand_real1, genrand_real2, genrand_real3, genrand_res53
    public  :: genrand_real, init_genrand_real, standard_normal

    integer,  parameter  :: intg = selected_int_kind( 9 )
    integer,  parameter  :: long = selected_int_kind( 18 )
    integer,  parameter  :: flot = selected_real_kind( 6, 37 )
    integer,  parameter  :: dobl = selected_real_kind( 15, 307 )

    integer,  public, parameter :: wi = intg
    integer,  public, parameter :: wl = long
    integer,  public, parameter :: wr = dobl

    ! Period parameters
    integer( kind = wi ), parameter :: n = 624_wi
    integer( kind = wi ), parameter :: m = 397_wi
    integer( kind = wi ), parameter :: hbs = bit_size( n ) / 2_wi
    integer( kind = wi ), parameter :: qbs = hbs / 2_wi
    integer( kind = wi ), parameter :: tbs = 3_wi * qbs

    integer( kind = wi )  :: mt(n)          ! the array for the state vector
    logical( kind = wi )  :: mtinit = .false._wi ! means mt[N] is not initialized
    integer( kind = wi )  :: mti = n + 1_wi ! mti==N+1 means mt[N] is not initialized

    interface genrand_real
        module procedure genrand_real3 
        end interface

    contains

        ! init_genrand_real was added by C. Predescu to help construct a simulacrum of 
        ! a parallel rng. It initializes the mt19937 rng with seeds constructed with 
        ! the help of the native rng. Different seeds are generated for different id's. 
        ! The subroutine must be called after the MPI environment has been initialized.
        ! We use the maximal number of 624 integers of 32 bits each for initialization,
        ! so that to maximize the Hamming distance for the states corresponding to
        ! different id's and, consequently, increase the independence of the generated
        ! parallel streams.  
        !
        subroutine init_genrand_real( id )

            integer, intent(in)  :: id 

            integer    ::  k
            real(wr)      ::  aa(624)
            integer(wi)   ::  init_seed(624)

            if(id < 0) then
                do k = 0, -2 * id
                    call random_number( aa )
                end do
            else
                do k = 0, 2 * id + 1
                    call random_number( aa )
                end do
            end if

            init_seed = int( huge(1_wi) * (2.0_wr * aa - 1.0_wr), kind = wi) 

            call init_by_array( init_seed )

        end subroutine init_genrand_real

        ! Added by C. Predescu. It returns a standard normal variable computed via the polar
        ! form of the Box-Muller method. Although the method generates two uncorrelated Gauss
        ! variables, we only use one. 

        function standard_normal( ) result(res)      
            real(dobl) :: res
            real(dobl) :: w

            do
                res = 2.0_dobl * genrand_real() - 1.0_dobl
                w   = 2.0_dobl * genrand_real() - 1.0_dobl
                w   = res**2 + w**2
                if(w < 1.0_dobl) exit
            end do
            w = sqrt(-2.0_dobl * log(w) / w)
            res = res * w
        end function standard_normal

        elemental function uiadd( a, b ) result( c )

            implicit none

            intrinsic :: ibits, ior, ishft

            integer( kind = wi ), intent( in )  :: a, b

            integer( kind = wi )  :: c

            integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

            a1 = ibits( a, 0, hbs )
            a2 = ibits( a, hbs, hbs )
            b1 = ibits( b, 0, hbs )
            b2 = ibits( b, hbs, hbs )
            s1 = a1 + b1
            s2 = a2 + b2 + ibits( s1, hbs, hbs )
            c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

        end function uiadd

        elemental function uisub( a, b ) result( c )

            implicit none

            intrinsic :: ibits, ior, ishft

            integer( kind = wi ), intent( in )  :: a, b

            integer( kind = wi )  :: c

            integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

            a1 = ibits( a, 0, hbs )
            a2 = ibits( a, hbs, hbs )
            b1 = ibits( b, 0, hbs )
            b2 = ibits( b, hbs, hbs )
            s1 = a1 - b1
            s2 = a2 - b2 + ibits( s1, hbs, hbs )
            c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

        end function uisub

        elemental  function uimlt( a, b ) result( c )

            implicit none

            intrinsic :: ibits, ior, ishft

            integer( kind = wi ), intent( in )  :: a, b

            integer( kind = wi )  :: c

            integer( kind = wi )  :: a0, a1, a2, a3
            integer( kind = wi )  :: b0, b1, b2, b3
            integer( kind = wi )  :: p0, p1, p2, p3

            a0 = ibits( a, 0, qbs )
            a1 = ibits( a, qbs, qbs )
            a2 = ibits( a, hbs, qbs )
            a3 = ibits( a, tbs, qbs )
            b0 = ibits( b, 0, qbs )
            b1 = ibits( b, qbs, qbs )
            b2 = ibits( b, hbs, qbs )
            b3 = ibits( b, tbs, qbs )
            p0 = a0 * b0
            p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
            p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
            p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
            c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
            c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
            c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )

        end function uimlt

        ! initializes mt[N] with a seed
        subroutine init_genrand( s )

            implicit none

            intrinsic :: iand, ishft, ieor, ibits

            integer( kind = wi ), intent( in )  :: s

            integer( kind = wi )  :: i, mult_a

            data mult_a /z'6C078965'/

            mtinit = .true._wi
            mt(1) = ibits( s, 0, 32 )
            do i = 2, n, 1
                mt(i) = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
                mt(i) = uimlt( mt(i), mult_a )
                mt(i) = uiadd( mt(i), uisub( i, 1_wi ) )
                ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
                ! In the previous versions, MSBs of the seed affect
                ! only MSBs of the array mt[].
                ! 2002/01/09 modified by Makoto Matsumoto
                mt(i) = ibits( mt(i), 0, 32 )
                ! for >32 bit machines
            end do

        end subroutine init_genrand

        ! initialize by an array with array-length
        ! init_key is the array for initializing keys
        ! key_length is its length
        subroutine init_by_array( init_key )

            implicit none

            intrinsic :: iand, ishft, ieor

            integer( kind = wi ), intent( in )  :: init_key(:)

            integer( kind = wi )  :: i, j, k, tp, key_length
            integer( kind = wi )  :: seed_d, mult_a, mult_b, msb1_d

            data seed_d /z'12BD6AA'/
            data mult_a /z'19660D'/
            data mult_b /z'5D588B65'/
            data msb1_d /z'80000000'/

            key_length = size( init_key, dim = 1 )
            call init_genrand( seed_d )
            i = 2_wi
            j = 1_wi
            do k = max( n, key_length ), 1, -1
                tp = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
                tp = uimlt( tp, mult_a )
                mt(i) = ieor( mt(i), tp )
                mt(i) = uiadd( mt(i), uiadd( init_key(j), uisub( j, 1_wi ) ) ) ! non linear
                mt(i) = ibits( mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
                i = i + 1_wi
                j = j + 1_wi
                if ( i > n ) then
                    mt(1) = mt(n)
                    i = 2_wi
                end if
                if ( j > key_length) j = 1_wi
            end do
            do k = n-1, 1, -1
                tp = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
                tp = uimlt( tp, mult_b )
                mt(i) = ieor( mt(i), tp )
                mt(i) = uisub( mt(i), uisub( i, 1_wi ) ) ! non linear
                mt(i) = ibits( mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
                i = i + 1_wi
                if ( i > n ) then
                    mt(1) = mt(n)
                    i = 2_wi
                end if
            end do
            mt(1) = msb1_d ! MSB is 1; assuring non-zero initial array 
        end subroutine init_by_array

        ! generates a random number on [0,0xffffffff]-interval
        function genrand_int32( ) result( y )

            implicit none

            intrinsic :: iand, ishft, ior, ieor, btest, ibset, mvbits

            integer( kind = wi )  :: y

            integer( kind = wi )  :: kk
            integer( kind = wi )  :: seed_d, matrix_a, matrix_b, temper_a, temper_b

            data seed_d   /z'5489'/
            data matrix_a /z'9908B0DF'/
            data matrix_b /z'0'/
            data temper_a /z'9D2C5680'/
            data temper_b /z'EFC60000'/

            if ( mti > n ) then ! generate N words at one time
                if ( .not. mtinit ) call init_genrand( seed_d ) ! if init_genrand() has not been called,
                ! a default initial seed is used
                do kk = 1, n-m, 1
                    y = ibits( mt(kk+1), 0, 31 )
                    call mvbits( mt(kk), 31, 1, y, 31 )
                    if ( btest( y, 0 ) ) then
                        mt(kk) = ieor( ieor( mt(kk+m), ishft( y, -1 ) ), matrix_a )
                    else
                        mt(kk) = ieor( ieor( mt(kk+m), ishft( y, -1 ) ), matrix_b )
                    end if
                end do
                do kk = n-m+1, n-1, 1
                    y = ibits( mt(kk+1), 0, 31 )
                    call mvbits( mt(kk), 31, 1, y, 31 )
                    if ( btest( y, 0 ) ) then
                        mt(kk) = ieor( ieor( mt(kk+m-n), ishft( y, -1 ) ), matrix_a )
                    else
                        mt(kk) = ieor( ieor( mt(kk+m-n), ishft( y, -1 ) ), matrix_b )
                    end if
                end do
                y = ibits( mt(1), 0, 31 )
                call mvbits( mt(n), 31, 1, y, 31 )
                if ( btest( y, 0 ) ) then
                    mt(kk) = ieor( ieor( mt(m), ishft( y, -1 ) ), matrix_a )
                else
                    mt(kk) = ieor( ieor( mt(m), ishft( y, -1 ) ), matrix_b )
                end if
                mti = 1_wi
            end if
            y = mt(mti)
            mti = mti + 1_wi
            ! Tempering
            y = ieor( y, ishft( y, -11) )
            y = ieor( y, iand( ishft( y, 7 ), temper_a ) )
            y = ieor( y, iand( ishft( y, 15 ), temper_b ) )
            y = ieor( y, ishft( y, -18 ) )

        end function genrand_int32

        ! generates a random number on [0,0x7fffffff]-interval
        function genrand_int31( ) result( i )

            implicit none

            intrinsic :: ishft

            integer( kind = wi )  :: i

            i = ishft( genrand_int32( ), -1 )

        end function genrand_int31

        ! generates a random number on [0,1]-real-interval
        function genrand_real1( ) result( r )

            implicit none

            real( kind = wr )  :: r

            integer( kind = wi )  :: a, a1, a0

            a = genrand_int32( )
            a0 = ibits( a, 0, hbs )
            a1 = ibits( a, hbs, hbs )
            r = real( a0, kind = wr ) / 4294967295.0_wr
            r = real( a1, kind = wr ) * ( 65536.0_wr / 4294967295.0_wr ) + r
            ! divided by 2^32-1

        end function genrand_real1

        ! generates a random number on [0,1)-real-interval
        function genrand_real2( ) result( r )

            implicit none

            intrinsic :: ibits

            real( kind = wr )  :: r

            integer( kind = wi )  :: a, a1, a0

            a = genrand_int32( )
            a0 = ibits( a, 0, hbs )
            a1 = ibits( a, hbs, hbs )
            r = real( a0, kind = wr ) / 4294967296.0_wr
            r = real( a1, kind = wr ) / 65536.0_wr + r
            ! divided by 2^32

        end function genrand_real2

        ! generates a random number on (0,1)-real-interval
        function genrand_real3( ) result( r )

            implicit none

            real( kind = wr )  :: r

            integer( kind = wi )  :: a, a1, a0

            a = genrand_int32( )
            a0 = ibits( a, 0, hbs )
            a1 = ibits( a, hbs, hbs )
            r = ( real( a0, kind = wr ) + 0.5_wr ) / 4294967296.0_wr
            r = real( a1, kind = wr ) / 65536.0_wr + r
            ! divided by 2^32

        end function genrand_real3

        ! generates a random number on [0,1) with 53-bit resolution
        function genrand_res53( )  result( r )

            implicit none

            intrinsic :: ishft

            real( kind = wr )  :: r

            integer( kind = wi )  :: a, a0, a1
            integer( kind = wi )  :: b, b0, b1

            a = ishft( genrand_int32( ), -5 )
            a0 = ibits( a, 0, hbs )
            a1 = ibits( a, hbs, hbs )
            b = ishft( genrand_int32( ), -6 )
            b0 = ibits( b, 0, hbs )
            b1 = ibits( b, hbs, hbs )
            r = real( a1, kind = wr ) / 2048.0_wr
            r = real( a0, kind = wr ) / 134217728.0_wr + r
            r = real( b1, kind = wr ) / 137438953472.0_wr + r
            r = real( b0, kind = wr ) / 9007199254740992.0_wr + r

        end function genrand_res53
        ! These real versions are due to Isaku Wada, 2002/01/09 added

    end module mt19937
    !
    ! =============================================================================
    !
    ! ------------------------------------------------------------------------------
    !
    Subroutine gasdev (r)
        !      
        use gauss_data
        use mt19937 
        implicit none
        !      
        real*8:: r, v1, v2, rsq, fac
        real*8:: randnm
        !
        if (iset .eq. 0)then
            randnm = genrand_real()
            v1 = 2.0d0*randnm - 1.0d0
            randnm = genrand_real()
            v2 = 2.0d0*randnm - 1.0d0
            rsq = v1**2 + v2**2

            do while ((rsq .ge. 1.0d0) .or. (rsq .eq. 0.0d0))
                randnm = genrand_real()
                v1 = 2.0d0*randnm - 1.0d0
                randnm = genrand_real()
                v2 = 2.0d0*randnm - 1.0d0
                rsq = v1**2 + v2**2     
            end do

            fac = sqrt(-2.0d0*log(rsq)/rsq)
            gset = v1*fac
            r = v2*fac
            iset = 1
        else
            r = gset
            iset = 0
        end if

        return
    end Subroutine gasdev
    !
    ! ================================================================================
    !




    module constant
        implicit none
        real*8,parameter ::mass=3.35d-26  !kg
        real*8,parameter ::womega=1.d0
        real*8,parameter ::convert_fs=1d-15  !s to fs
        real*8,parameter ::dt=20.d0*convert_fs
        real*8,parameter ::lgama=0.1d0/convert_fs
        integer, parameter :: nne=13
        integer,parameter :: ndim=3
        real*8, parameter :: Temperature=14.0d0  ! K
        real*8,parameter ::hbar=1.05457266d-34
        real*8, parameter :: kb=1.3806488d-23
        real*8:: beta=1d0/(kb*Temperature)
        real*8,parameter:: ESI_to_au=43.60d-19
        !...............
        !Lenard Jones
        real*8, parameter :: elj=35.6d0
        real*8, parameter :: dlj=2.749d-10
        real*8, parameter :: tstep=dlj/5



    endmodule constant






    program main
        use mt19937
        use gauss_data
        use constant
        implicit none
        real*8 ::x(nne,ndim),E,force(nne,ndim),d2v,E_tot,p(nne,ndim)
        real*8,parameter ::ntime=5d-10
        integer,parameter ::eq_step=int(2000000.d0*convert_fs/dt)
        integer ::nsample=1000!0
        integer,parameter ::nstep=int(200000000.d0*convert_fs/dt)
        integer ::i,j,n,m
        integer,parameter ::nblock=20
        real*8 ::sum_energy,Ek,Ep,random,ran
        real*8 ::E_block(nblock),error_bar,E0,average_E,step_energy(nstep),sum_E(nstep)
        integer ::sysclc



        call SYSTEM_CLOCK( sysclc )
        call init_genrand_real(mod(sysclc,1000000))

        open(13,file='energy.dat')


        sum_E=0.d0
        do j=1,1!nsample
            
            sum_energy=0.d0
            call LJ_init(x)
          !  write(*,*) x(:,1)
            do n=1, nne
                do m=1, 3
                    call gasdev(ran)
                    p(n,m)=ran*sqrt(mass/beta)
                end do
            end do

            call forcelj(force,x)
            do i=1,eq_step
                call run(x,p,force,1)
                call forcelj(force,x)
                call run(x,p,force,2)
            enddo
            
            
            do i=1,nstep
                call run(x,p,force,1)
                call forcelj(force,x)
                call run(x,p,force,2)

                call energy(x,p,Ep,Ek,E_tot)
                !     write(12,*) E_tot
                write(13,*) Ep,E_tot
          !      if(i.eq.1) E0=E
          !      average_E=average_E+E
          !      step_energy(i)= E
            enddo

            average_E=average_E/nstep
            E0=E0-average_E
            step_energy=step_energy-average_E

            do i=1,nstep
                sum_E(i)=sum_E(i)+E0*step_energy(i)
            enddo

        enddo



close(13)


        !  close(12)

        !
        !
    endprogram main



    subroutine energy(x,p,Ep,Ek,E_tot)
        use constant
        implicit none
        real*8 ,intent(in) ::x(nne,ndim),p(nne,ndim)
        real*8 ,intent(out) ::E_tot
        real*8 ::Ek,Ep
        integer ::i,j


        call vf(Ep, x)
        Ek=0.d0
        do i=1,ndim
            do j=1,nne
                Ek=Ek+p(j,i)**2*0.5d0/mass
            enddo
        enddo

        E_tot=Ek+EP

        return
    endsubroutine energy



    subroutine run(x,p,force,flag)
        use mt19937
        use gauss_data
        use constant
        implicit none
        integer,intent(in) ::flag
        real*8,intent(inout) ::x(nne,ndim),p(nne,ndim)
        real*8,intent(in) ::force(nne,ndim)
        real*8 ::ran 
        real*8 ::random,sigma,c1k,c2k
        integer ::i,j

        if(flag.eq.1) then
            p = p - 0.50d0 * force * dt
            x = x + 0.50d0 * p * dt / mass

            c1k = exp( - dt * lgama )
            c2k = sqrt( 1.0d0 - c1k**2 )

            do i=1,ndim
                do j=1,nne
                    call gasdev(ran)
                    p(j,i) = c1k * p(j,i) + c2k * sqrt( mass/beta ) * ran
                enddo
            enddo


            x = x + 0.50d0 * p * dt / mass
        elseif(flag.eq.2) then

            p = p - 0.50d0 * force * dt

        endif


        return
    endsubroutine run


    subroutine vf(v ,x)
        use constant
        implicit none
        !Psssed variables
        !real*8 :: epot, xlj(nne,3)
        real*8 :: v
        real*8 :: x(nne,3)
        !Local variables
        real*8,external :: vlj, vc
        integer :: i,j,k,l
        real*8 :: rcm(3), temp(3)

        !x(:,:,1)=xlj(:,:)

        rcm(:)=sum(x(:,:),1)/nne
        v=0.0d0
        !do k=1,1
        do i=1, nne
            v=v+vc(x(i,:),rcm(:))
            if(i==1)cycle
            do j=1, i-1
                !if(i==j)cycle
                v=v+vlj(x(i,:),x(j,:))
            end do
        end do
        !end do
        v=v*kb
        !epot=v(1)

        return
    end subroutine

    subroutine forcelj(f,x)
        use constant
        implicit none
        !Psssed variables
        !real*8 :: flj(nne,3), xlj(nne,3)
        real*8 :: f(nne,3)
        real*8 :: x(nne,3)
        !Local variables
        integer :: i,j,k,l
        real*8 :: rcm(3), temp(3)
        real*8 :: rc=2*dlj

        !x(:,:,1)=xlj(:,:)

        rcm(:)=sum(x(:,:),1)/nne
        f=0.0d0
        !do k=1,1
        do i=1, nne
            temp=0.0d0
            do l=1,nne
                temp=temp+elj*10*(sum((x(l,:)-rcm(:))**2))**9*2*(x(l,:)-rcm(:))*&
                    (-1.0d0/nne)/(rc)**20
            end do
            f(i,:)=f(i,:)+elj*10*(sum((x(i,:)-rcm(:))**2))**9*2*(x(i,:)-&
                rcm(:))/(rc)**20+temp
            if(i==1)cycle
            do j=1, i-1
                !if(i==j)cycle
                f(i,:)=f(i,:)+4*elj*(dlj**12/(sum((x(i,:)-x(j,:))**2))**7*&
                    (-12*(x(i,:)-x(j,:)))-dlj**6/(sum((x(i,:)-x(j,:))**2))**4*&
                    (-6*(x(i,:)-x(j,:))))
                f(j,:)=f(j,:)+4*elj*(dlj**12/(sum((x(i,:)-x(j,:))**2))**7*&
                    (-12*(x(j,:)-x(i,:)))-dlj**6/(sum((x(i,:)-x(j,:))**2))**4*&
                    (-6*(x(j,:)-x(i,:))))
            end do
        end do
        !end do
        !   f=-f*kb
        f=f*kb

        !flj(:,:)=f(:,:,1)

        return
    end subroutine

    function vc(x,r)
        use constant
        implicit none
        !Passed variables
        real*8 :: x(3), r(3)
        !Local variables
        real*8 :: vc
        real*8 :: rc=2*dlj

        vc=elj*((sum((x-r)**2))/rc**2)**10

        return
    end function

    function vlj(x1,x2)
        use constant
        implicit none
        !Passed variables
        real*8 :: x1(3),x2(3)
        !Locl variables
        real*8 :: vlj

        vlj=4*elj*(dlj**12/(sum((x1-x2)**2))**6-dlj**6/(sum((x1-x2)**2))**3)

        return
    end function

    subroutine LJ_init(x)
        use constant
        implicit none
        !Passed variables
        !real*8 :: xlj(nne,3)
        real*8 :: x(nne,3),u(nne,3)
        !Local variables
        real*8 :: fie=(1+dsqrt(5.0d0))/2
        real*8 :: rmin=(2.0d0)**(1.0d0/6)*dlj/2
        integer :: i,j,k
        integer :: LDA = 1
        character(len=1) :: TRANS = 'N'

        !do i=1,1
        x(1,:)=(/0.0d0,0.0d0,0.0d0/)
        x(2,:)=(/0.0d0,rmin,fie*rmin/)
        x(3,:)=(/0.0d0,-rmin,fie*rmin/)
        x(4,:)=(/0.0d0,rmin,-fie*rmin/)
        x(5,:)=(/0.0d0,-rmin,-fie*rmin/)
        x(6,:)=(/rmin,fie*rmin,0.0d0/)
        x(7,:)=(/-rmin,fie*rmin,0.0d0/)
        x(8,:)=(/rmin,-fie*rmin,0.0d0/)
        x(9,:)=(/-rmin,-fie*rmin,0.0d0/)
        x(10,:)=(/fie*rmin,0.0d0,rmin/)
        x(11,:)=(/fie*rmin,0.0d0,-rmin/)
        x(12,:)=(/-fie*rmin,0.0d0,rmin/)
        x(13,:)=(/-fie*rmin,0.0d0,-rmin/)
        !end do

        !xlj(:,:)=x(:,:,1)

        return
    end subroutine LJ_init






