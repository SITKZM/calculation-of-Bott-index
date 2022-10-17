!gfortran Bott_index_SQU_PBC.f90 -I/usr/local/include -llapack95 -llapack -lblas -o Bott_index_SQU_PBC.out
!./Bott_index_SQU_PBC.out
program Bott_index_SQU_PBC
    !use F95_LAPACK

    implicit none
    !parameter
    integer, parameter :: N_x = 20, N_y = 20, N = N_x*N_y
    real, parameter :: mu = -4.0, h_z = 0.5, V = 0.5, Delta = 0.2, pi = 4*atan(1.)
    !variable
    integer :: INFO
    complex :: Hamiltonian(4*N, 4*N), WORK(2*4*N-1)
    real :: quasiparticle_energy(4*N), RWORK(3*4*N-2), Bott_index

    call make_BdG_Hamiltonian(Hamiltonian)
    call CHEEV('V', 'U', 4*N, Hamiltonian, 4*N, quasiparticle_energy, WORK, 2*4*N-1, RWORK, INFO)
    call get_Bott_index(Bott_index, Hamiltonian)

    print *, Bott_index
    contains
        subroutine make_BdG_Hamiltonian(H)
            complex :: H(4*N, 4*N)
            integer :: i

            H = 0

            ! diagonal elements
            do i = 1, N
                !particle
                !spin up
                H(2*i-1,2*i-1) = -mu + h_z
                !spin down
                H(2*i, 2*i) = -mu - h_z
                !hole
                !spin up
                H(2*i-1+2*N, 2*i-1+2*N) = -H(2*i-1,2*i-1)
                !spin down
                H(2*i+2*N, 2*i+2*N) = -H(2*i, 2*i)
            end do

            ! hopping elements
            do i = 1, N
                !along the x-axis
                !edge
                if (modulo(i,N_x) == 0) then
                    !particle
                    !spin up
                    H(2*i-1, 2*(i+1-N_x)-1) = -1
                    H(2*(i+1-N_x)-1, 2*i-1) = H(2*i-1, 2*(i+1-N_x)-1)
                    !spin down
                    H(2*i, 2*(i+1-N_x)) = -1
                    H(2*(i+1-N_x), 2*i) = H(2*i, 2*(i+1-N_x))
                    !hole
                    !spin up
                    H(2*i-1+2*N, 2*(i+1-N_x)-1+2*N) = 1
                    H(2*(i+1-N_x)-1+2*N, 2*i-1+2*N) = H(2*i-1+2*N, 2*(i+1-N_x)-1+2*N)
                    !spin down
                    H(2*i+2*N, 2*(i+1-N_x)+2*N) = 1
                    H(2*(i+1-N_x)+2*N, 2*i+2*N) = H(2*i+2*N, 2*(i+1-N_x)+2*N)
                !bulk
                else
                    !particle
                    !spin up
                    H(2*i-1, 2*(i+1)-1) = -1
                    H(2*(i+1)-1, 2*i-1) = H(2*i-1, 2*(i+1)-1)
                    !spin down
                    H(2*i, 2*(i+1)) = -1
                    H(2*(i+1), 2*i) = H(2*i, 2*(i+1))
                    !hole
                    !spin up
                    H(2*i-1+2*N, 2*(i+1)-1+2*N) = 1
                    H(2*(i+1)-1+2*N, 2*i-1+2*N) = H(2*i-1+2*N, 2*(i+1)-1+2*N)
                    !spin down
                    H(2*i+2*N, 2*(i+1)+2*N) = 1
                    H(2*(i+1)+2*N, 2*i+2*N) = H(2*i+2*N, 2*(i+1)+2*N)
                end if

                !along the y-axis
                !edge
                if (i > (N_y - 1)*N_x) then
                    !particle
                    !spin up
                    H(2*i-1, 2*(modulo(i+1+N_x, N)-1)-1) = -1
                    H(2*(modulo(i+1+N_x, N)-1)-1, 2*i-1) = H(2*i-1, 2*(modulo(i+1+N_x, N)-1)-1)
                    !spin down
                    H(2*i, 2*(modulo(i+1+N_x, N)-1)) = -1
                    H(2*(modulo(i+1+N_x, N)-1), 2*i) = H(2*i, 2*(modulo(i+1+N_x, N)-1))
                    !hole
                    !spin up
                    H(2*i-1+2*N, 2*(modulo(i+1+N_x, N)-1)-1+2*N) = 1
                    H(2*(modulo(i+1+N_x, N)-1)-1+2*N, 2*i-1+2*N) = H(2*i-1+2*N, 2*(modulo(i+1+N_x, N)-1)-1+2*N)
                    !spin down
                    H(2*i+2*N, 2*(modulo(i+1+N_x, N)-1)+2*N) = 1
                    H(2*(modulo(i+1+N_x, N)-1)+2*N, 2*i+2*N) = H(2*i+2*N, 2*(modulo(i+1+N_x, N)-1)+2*N)
                else
                    !particle
                    !spin up
                    H(2*i-1, 2*(i+N_x)-1) = -1
                    H(2*(i+N_x)-1, 2*i-1) = H(2*i-1, 2*(i+N_x)-1)
                    !spin down
                    H(2*i, 2*(i+N_x)) = -1
                    H(2*(i+N_x), 2*i) = H(2*i, 2*(i+N_x))
                    !hole
                    !spin up
                    H(2*i-1+2*N, 2*(i+N_x)-1+2*N) = 1
                    H(2*(i+N_x)-1+2*N, 2*i-1+2*N) = H(2*i-1+2*N, 2*(i+N_x)-1+2*N)
                    !spin down
                    H(2*i+2*N, 2*(i+N_x)+2*N) = 1
                    H(2*(i+N_x)+2*N, 2*i+2*N) = H(2*i+2*N, 2*(i+N_x)+2*N)
                end if
            end do

            !spin-orbit coupling
            do i = 1, N
                !along the x-axis
                !edge
                if (modulo(i, N_x) == 0) then
                    !particle
                    !hat{x} = (-1,0)
                    H(2*i-1, 2*(i+1-N_x)) = V
                    H(2*(i+1-N_x), 2*i-1) = H(2*i-1, 2*(i+1-N_x))
                    !hat{x} = (1,0)
                    H(2*i, 2*(i+1-N_x)-1) = -V
                    H(2*(i+1-N_x)-1, 2*i) = H(2*i, 2*(i+1-N_x)-1)
                    !hole
                    !hat{x} = (-1,0)
                    H(2*i-1+2*N, 2*(i+1-N_x)+2*N) = -V
                    H(2*(i+1-N_x)+2*N, 2*i-1+2*N) = H(2*i-1+2*N, 2*(i+1-N_x)+2*N)
                    !hat{x} = (1,0)
                    H(2*i+2*N, 2*(i+1-N_x)-1+2*N) = V
                    H(2*(i+1-N_x)-1+2*N, 2*i+2*N) = H(2*i+2*N, 2*(i+1-N_x)-1+2*N)
                !bulk
                else
                    !particle
                    !hat{x} = (-1,0)
                    H(2*i-1, 2*(i+1)) = V
                    H(2*(i+1), 2*i-1) = H(2*i-1, 2*(i+1))
                    !hat{x} = (1,0)
                    H(2*i, 2*(i+1)-1) = -V
                    H(2*(i+1)-1, 2*i) = H(2*i, 2*(i+1)-1)
                    !hole
                    !hat{x} = (-1,0)
                    H(2*i-1+2*N, 2*(i+1)+2*N) = -V
                    H(2*(i+1)+2*N, 2*i-1+2*N) = H(2*i-1+2*N, 2*(i+1)+2*N)
                    !hat{x} = (1,0)
                    H(2*i+2*N, 2*(i+1)-1+2*N) = V
                    H(2*(i+1)-1+2*N, 2*i+2*N) = H(2*i+2*N, 2*(i+1)-1+2*N)
                end if

                !along the y-axis
                !edge
                if (i > (N_y - 1)*N_x) then
                    !particle
                    !hat{y} = (0,-1)
                    H(2*i-1, 2*(modulo(i+1+N_x, N)-1)) = cmplx(0, V)
                    H(2*(modulo(i+1+N_x, N)-1), 2*i-1) = conjg(H(2*i-1, 2*(modulo(i+1+N_x, N)-1)))
                    !hat{y} = (0,1)
                    H(2*(modulo(i+1+N_x, N)-1)-1, 2*i) = cmplx(0, V)
                    H(2*i, 2*(modulo(i+1+N_x, N)-1)-1) = conjg(H(2*(modulo(i+1+N_x, N)-1)-1, 2*i))
                    !hole
                    !hat{y} = (0,-1)
                    H(2*i-1+2*N, 2*(modulo(i+1+N_x, N)-1)+2*N) = -conjg(cmplx(0, V))
                    H(2*(modulo(i+1+N_x, N)-1)+2*N, 2*i-1+2*N) = conjg(H(2*i-1+2*N, 2*(modulo(i+1+N_x, N)-1)+2*N))
                    !hat{y} = (0,1)
                    H(2*(modulo(i+1+N_x, N)-1)-1+2*N, 2*i+2*N) = -conjg(cmplx(0, V))
                    H(2*i+2*N, 2*(modulo(i+1+N_x, N)-1)-1+2*N) = conjg(H(2*(modulo(i+1+N_x, N)-1)-1+2*N, 2*i+2*N))
                else
                    !particle
                    !hat{y} = (0,-1)
                    H(2*i-1, 2*(i+N_x)) = -cmplx(0, V)
                    H(2*(i+N_x), 2*i-1) = conjg(H(2*i-1, 2*(i+N_x)))
                    !hat{y} = (0,1)
                    H(2*(i+N_x)-1, 2*i) = cmplx(0, V)
                    H(2*i, 2*(i+N_x)-1) = conjg(H(2*(i+N_x)-1, 2*i))
                    !hole
                    !hat{y} = (0,-1)
                    H(2*i-1+2*N, 2*(i+N_x)+2*N) = conjg(cmplx(0, V))
                    H(2*(i+N_x)+2*N, 2*i-1+2*N) = conjg(H(2*i-1+2*N, 2*(i+N_x)+2*N))
                    !hat{y} = (0,1)
                    H(2*(i+N_x)-1+2*N, 2*i+2*N) = -conjg(cmplx(0, V))
                    H(2*i+2*N, 2*(i+N_x)-1+2*N) = conjg(H(2*(i+N_x)-1+2*N, 2*i+2*N))
                end if
            end do    

            !pair potential
            do i = 1, N
                !particle
                H(2*i, 2*i-1+2*N) = Delta
                H(2*i-1+2*N, 2*i) = conjg(H(2*i, 2*i-1+2*N))
                !hole
                H(2*i-1, 2*i+2*N) = -Delta
                H(2*i+2*N, 2*i-1) = conjg(H(2*i-1, 2*i+2*N))
            end do
        end subroutine make_BdG_Hamiltonian

        subroutine get_Bott_index(B, eigm)
            integer :: i
            real :: B, RWORK_B(2*4*N)
            complex :: eigm(4*N, 4*N), eigval(4*N)
            complex :: X(4*N, 4*N), Y(4*N, 4*N), U_X(4*N, 4*N), U_Y(4*N, 4*N), U_prod(4*N, 4*N)
            complex :: P(4*N, 4*N), Q(4*N, 4*N)
            complex :: VL(1, 4*N), VR(1, 4*N), WORK_B(2*4*N)
            !position operators
            X = 0
            Y = 0
            do i = 1, N
                !particle
                !spin up
                X(2*i-1, 2*i-1) = exp(cmplx(0, 2*pi*real(modulo(i-1, N_x))/N_x))
                Y(2*i-1, 2*i-1) = exp(cmplx(0, 2*pi*real(floor(real(i-1)/N_x))/N_y))
                !spin down
                X(2*i, 2*i) = exp(cmplx(0, 2*pi*real(modulo(i-1, N_x))/N_x))
                Y(2*i, 2*i) = exp(cmplx(0, 2*pi*real(floor(real(i-1)/N_x))/N_y))
                !hole
                !spin up
                X(2*i-1+2*N, 2*i-1+2*N) = exp(cmplx(0, 2*pi*real(modulo(i-1, N_x))/N_x))
                Y(2*i-1+2*N, 2*i-1+2*N) = exp(cmplx(0, 2*pi*real(floor(real(i-1)/N_x))/N_y))
                !spin down
                X(2*i+2*N, 2*i+2*N) = exp(cmplx(0, 2*pi*real(modulo(i-1, N_x))/N_x))
                Y(2*i+2*N, 2*i+2*N) = exp(cmplx(0, 2*pi*real(floor(real(i-1)/N_x))/N_y))
            end do

            !projection operators
            P = 0
            Q = 0
            do i = 1, 2*N
                P(i,i) = 1
                Q(4*N+1-i, 4*N+1-i) = 1
            end do

            !projected position operators
            U_X = matmul(P, matmul(transpose(conjg(eigm)), matmul(X, matmul(eigm, P)))) + Q
            U_Y = matmul(P, matmul(transpose(conjg(eigm)), matmul(Y, matmul(eigm, P)))) + Q
            U_prod = matmul(U_Y, matmul(U_X, matmul(transpose(conjg(U_Y)), transpose(conjg(U_X)))))

            !Bott index
            B = 0
            call CGEEV('N', 'N', 4*N, U_prod, 4*N, eigval, VL, 1, VR, 1, WORK_B, 2*4*N, RWORK_B, INFO)
            do i = 1, 4*N
                B = B + aimag(log(eigval(i)/abs(eigval(i))))
            end do
            B = B/(2*pi)
        end subroutine get_Bott_index
end program Bott_index_SQU_PBC