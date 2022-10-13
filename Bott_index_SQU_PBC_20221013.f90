!calculate Bott index of square lattice with periodic boundary condition when pair potential is given(not self-consist solution)
!This program does NOT produce the correct results.
program Bott_index_SQU_OBC
    use F95_LAPACK

    implicit none
    integer :: i, unit_write_result
    real, parameter :: pi=4*atan(1.)
    !physical quantity
    integer, parameter :: size = 20
    real, parameter :: chemical_potential = -4.0, magnetic_field = 0.5, spin_orbit_coupling = 0.5
    real, parameter :: pair_potential = 0.2
    real :: Bott_index
    real :: eigenvalues(4*size**2)
    complex :: Hamiltonian(4*size**2, 4*size**2)
    !to use lapack
    integer :: info
    complex :: work(2*4*size**2 - 1)
    real :: rwork(3*4*size**2-2)

    print *, "mu =", chemical_potential

    !make Hamiltonian
    call make_BdG_Hamiltonian(size, Hamiltonian, chemical_potential, magnetic_field, spin_orbit_coupling,&
                              pair_potential)
    
    !diagonalize
    call CHEEV('V', 'L', 4*size**2, Hamiltonian, 4*size**2, eigenvalues, work, 2*4*size**2 - 1, rwork, info)

    !calculate Bott index
    call get_Bott_index(size, Bott_index, Hamiltonian)

    !output data
    !quasiparticle energy level
    open(newunit=unit_write_result, file="eigvals.txt")
        do i = 1, 4*size**2
            write(unit_write_result, *) eigenvalues(i)
        end do
    close(unit_write_result)
    !Bott index
    print *, "Bott index:", Bott_index
contains
    subroutine make_BdG_Hamiltonian(L, H, mu, h_z, V, Delta)
        implicit none
        integer :: n, L
        real :: mu, h_z, V, Delta
        complex :: H(4*L**2, 4*L**2)

        n = L**2
        H = 0

        !chemical potential and Zeeman term
        do i = 1, n
            !particle
            !spin up
            H(i, i) = -mu + h_z
            !spin down
            H(i+n, i+n) = -mu - h_z
            !hole
            !spin up
            H(i+2*n, i+2*n) = -H(i, i)
            !spin down
            H(i+3*n, i+3*n) = -H(i+n, i+n)
        end do

        !hopping
        !along x-axis
        do i = 1, n
            !edge
            if ( modulo(i, L) == 0 ) then
                !particle
                !spin up
                H(i, i-L+1) = -1
                H(i-L+1, i) = H(i, i-L+1)
                !spin down
                H(i+n, i-L+1+n) = -1
                H(i-L+1+n, i+n) = H(i+n, i-L+1+n)
                !hole
                !spin up
                H(i+2*n, i-L+1+2*n) = 1
                H(i-L+1+2*n, i+2*n) = H(i+2*n, i-L+1+2*n)
                !spin down
                H(i+3*n, i-L+1+3*n) = 1
                H(i-L+1+3*n, i+3*n) = H(i+3*n, i-L+1+3*n)
            !bulk
            else
                !particle
                !spin up
                H(i, i+1) = -1
                H(i+1, i) = H(i, i+1)
                !spin down
                H(i+n, i+1+n) = -1
                H(i+1+n, i+n) = H(i+n, i+1+n)
                !hole
                !spin up
                H(i+2*n, i+1+2*n) = 1
                H(i+1+2*n, i+2*n) = H(i+2*n, i+1+2*n)
                !spin down
                H(i+3*n, i+1+3*n) = 1
                H(i+1+3*n, i+3*n) = H(i+3*n, i+1+3*n)
            end if
        end do
        !along y-axis
        do i = 1, n
            !edge
            if ( i > n - L ) then
                !particle
                !spin up
                H(i, mod(i+L, n)) = -1
                H(mod(i+L, n), i) = H(i, mod(i+L, n))
                !spin down
                H(i+n, mod(i+L, n)+n) = -1
                H(mod(i+L, n)+n, i+n) = H(i+n, mod(i+L, n)+n)
                !hole
                !spin up
                H(i+2*n, mod(i+L, n)+2*n) = 1
                H(mod(i+L, n)+2*n, i+2*n) = H(i+2*n, mod(i+L, n)+2*n)
                !spin down
                H(i+3*n, mod(i+L, n)+3*n) = 1
                H(mod(i+L, n)+3*n, i+3*n) = H(i+3*n, mod(i+L, n)+3*n)
            !bulk
            else
                !particle
                !spin up
                H(i, i+L) = -1
                H(i+L, i) = H(i, i+L)
                !spin down
                H(i+n, i+L+n) = -1
                H(i+L+n, i+n) = H(i+n, i+L+n)
                !hole
                !spin up
                H(i+2*n, i+L+2*n) = 1
                H(i+L+2*n, i+2*n) = H(i+2*n, i+L+2*n)
                !spin down
                H(i+3*n, i+L+3*n) = 1
                H(i+L+3*n, i+3*n) = H(i+3*n, i+L+3*n)
            end if
        end do

        !spin-orbit coupling
        !along x-axis
        do i = 1, n
            !edge
            if ( modulo(i, L) == 0 ) then
                !particle
                !hat{x}=(-1,0)
                H(i, i-L+1+n) = V
                H(i-L+1+n, i) = H(i, i-L+1+n)
                !hat{x}=(1,0)
                H(i+n, i-L+1) = -V
                H(i-L+1, i+n) = H(i+n, i-L+1)
                !hole
                !hat{x}=(-1,0)
                H(i+2*n, i-L+1+3*n) = -H(i, i-L+1+n)
                H(i-L+1+3*n, i+2*n) = H(i+2*n, i-L+1+3*n)
                !hat{x}=(1,0)
                H(i+3*n, i-L+1+2*n) = -H(i+n, i-L+1)
                H(i-L+1+2*n, i+3*n) = H(i+3*n, i-L+1+2*n)
            !bulk
            else
                !particle
                !hat{x}=(-1,0)
                H(i, i+1+n) = V
                H(i+1+n, i) = H(i, i+1+n)
                !hat{x}=(1,0)
                H(i+n, i+1) = -V
                H(i+1, i+n) = H(i+n, i+1)
                !hole
                !hat{x}=(-1,0)
                H(i+2*n, i+1+3*n) = -H(i, i+1+n)
                H(i+1+3*n, i+2*n) = H(i+2*n, i+1+3*n)
                !hat{x}=(1,0)
                H(i+3*n, i+1+2*n) = -H(i+n, i+1)
                H(i+1+2*n, i+3*n) = H(i+3*n, i+1+2*n)
            end if
        end do
        !along y-axis
        do i = 1, n
            !edge
            if ( i > n - L ) then
                !particle
                !hat{y} = (0,-1)
                H(i, mod(i+L, n)+n) = cmplx(0, -V)
                H(mod(i+L, n)+n, i) = conjg(H(i, mod(i+L, n)+n))
                !hat{y} = (0,1)
                H(i+n, mod(i+L, n)) = cmplx(0, V)
                H(mod(i+L, n), i+n) = conjg(H(i+n, mod(i+L, n)))
                !hole
                !hat{y} = (0,-1)
                H(i+2*n, mod(i+L, n)+3*n) = -conjg(H(i, mod(i+L, n)+n))
                H(mod(i+L, n)+3*n, i+2*n) = conjg(H(i+2*n, mod(i+L, n)+3*n))
                !hat{y} = (0,1)
                H(i+3*n, mod(i+L, n)+2*n) = -conjg(H(i+n, mod(i+L, n)))
                H(mod(i+L, n)+2*n, i+3*n) = conjg(H(i+3*n, mod(i+L, n)+2*n))
            !bulk
            else
                !particle
                !hat{y} = (0,-1)
                H(i, i+L+n) = cmplx(0, -V)
                H(i+L+n, i) = conjg(H(i, i+L+n))
                !hat{y} = (0,1)
                H(i+n, i+L) = cmplx(0, V)
                H(i+L, i+n) = conjg(H(i+n, i+L))
                !hole
                !hat{y} = (0,-1)
                H(i+2*n, i+L+3*n) = -conjg(H(i, i+L+n))
                H(i+L+3*n, i+2*n) = conjg(H(i+2*n, i+L+3*n))
                !hat{y} = (0,1)
                H(i+3*n, i+L+2*n) = -conjg(H(i+n, i+L))
                H(i+L+2*n, i+3*n) = conjg(H(i+3*n, i+L+2*n))
            end if
        end do

        !pair potential
        do i = 1, n
            !particle
            H(i, i+3*n) = Delta
            H(i+3*n, i) = conjg(H(i, i+3*n))
            !hole
            H(i+n, i+2*n) = -Delta
            H(i+2*n, i+n) = conjg(H(i+n, i+2*n))
        end do
    end subroutine make_BdG_Hamiltonian

    subroutine get_Bott_index(L, B, eigm)
        implicit none
        integer :: n, L
        real :: B
        real :: rwork_g(2*4*L**2)
        complex :: X(4*L**2, 4*L**2), Y(4*L**2, 4*L**2)
        complex :: eigenvalues_g(4*L**2), VL(1, 4*L**2), VR(1, 4*L**2), work_g(2*4*L**2)
        complex :: eigm(4*L**2, 4*L**2), P(4*L**2, 4*L**2), Q(4*L**2, 4*L**2)
        complex :: U_X(4*L**2, 4*L**2), U_Y(4*L**2, 4*L**2), Uprod(4*L**2, 4*L**2)

        n = L**2
        !make rescaled position matrix
        X = 0
        Y = 0
        do i = 1, n
            !particle
            !spin up
            X(i, i) = exp(cmplx(0, 2*pi*real(mod(i-1, L))/(L-1)))
            Y(i, i) = exp(cmplx(0, 2*pi*(real(ceiling(real(i)/L)) - 1)/(L-1)))
            !spin down
            X(i+n, i+n) = X(i, i)
            Y(i+n, i+n) = Y(i, i)
            
            !hole
            !spin up
            X(i+2*n, i+2*n) = X(i, i)
            Y(i+2*n, i+2*n) = Y(i, i)
            !spin down
            X(i+3*n, i+3*n) = X(i, i)
            Y(i+3*n, i+3*n) = Y(i, i)
        end do
        
        !transform X, Y into quasiparticle basis
        X = matmul(transpose(conjg(eigm)), matmul(X, eigm))
        Y = matmul(transpose(conjg(eigm)), matmul(Y, eigm))

        !make P, Q matrices
        P = 0
        Q = 0
        do i = 1, 2*n
            P(i, i) = 1
            Q(4*n-i+1, 4*n-i+1) = 1
        end do

        !make projected position operators
        U_X = matmul(P, matmul(X, P)) + Q
        U_Y = matmul(P, matmul(Y, P)) + Q

        !get eigenvalues of the product
        Uprod = matmul(U_Y, matmul(U_X, matmul(transpose(conjg(U_Y)), transpose(conjg(U_X)))))
        call CGEEV('N', 'N', 4*n, Uprod, 4*n, eigenvalues_g, VL, 1, VR, 1, work_g, 2*4*n, rwork_g, info)

        ! get Bott index
        B = 0
        do i = 1, 4*n
            B = B + aimag(log(eigenvalues_g(i))/abs(eigenvalues_g(i)))
        end do
        B = B/(2*pi)
    end subroutine get_Bott_index
end program Bott_index_SQU_OBC