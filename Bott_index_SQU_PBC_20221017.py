from numpy.lib.type_check import imag
import numpy as np
import numpy.linalg as LA

def make_BdG_Hamiltonian(N_x, N_y, mu, h_z, V, Delta):
    N = int(N_x * N_y)
    Hamiltonian = np.zeros((4*N, 4*N), dtype=complex)


    # diagonal elements
    for i in range(N):
        #particle
        #spin up
        Hamiltonian[2*i,2*i] = -mu + h_z
        #spin down
        Hamiltonian[2*i+1, 2*i+1] = -mu - h_z
        #hole
        #spin up
        Hamiltonian[2*i+2*N, 2*i+2*N] = -Hamiltonian[2*i,2*i]
        #spin down
        Hamiltonian[2*i+1+2*N, 2*i+1+2*N] = -Hamiltonian[2*i+1, 2*i+1]


    # hopping elements
    for i in range(N):
        #along the x-axis
        #edge
        if (i+1) % N_x == 0:
            #particle
            #spin up
            Hamiltonian[2*i, 2*(i+1-N_x)] = -1
            Hamiltonian[2*(i+1-N_x), 2*i] = Hamiltonian[2*i, 2*(i+1-N_x)]
            #spin down
            Hamiltonian[2*i+1, 2*(i+1-N_x)+1] = -1
            Hamiltonian[2*(i+1-N_x)+1, 2*i+1] = Hamiltonian[2*i+1, 2*(i+1-N_x)+1]
            #hole
            #spin up
            Hamiltonian[2*i+2*N, 2*(i+1-N_x)+2*N] = 1
            Hamiltonian[2*(i+1-N_x)+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*(i+1-N_x)+2*N]
            #spin down
            Hamiltonian[2*i+1+2*N, 2*(i+1-N_x)+1+2*N] = 1
            Hamiltonian[2*(i+1-N_x)+1+2*N, 2*i+1+2*N] = Hamiltonian[2*i+1+2*N, 2*(i+1-N_x)+1+2*N]
        #bulk
        else:
            #particle
            #spin up
            Hamiltonian[2*i, 2*(i+1)] = -1
            Hamiltonian[2*(i+1), 2*i] = Hamiltonian[2*i, 2*(i+1)]
            #spin down
            Hamiltonian[2*i+1, 2*(i+1)+1] = -1
            Hamiltonian[2*(i+1)+1, 2*i+1] = Hamiltonian[2*i+1, 2*(i+1)+1]
            #hole
            #spin up
            Hamiltonian[2*i+2*N, 2*(i+1)+2*N] = 1
            Hamiltonian[2*(i+1)+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*(i+1)+2*N]
            #spin down
            Hamiltonian[2*i+1+2*N, 2*(i+1)+1+2*N] = 1
            Hamiltonian[2*(i+1)+1+2*N, 2*i+1+2*N] = Hamiltonian[2*i+1+2*N, 2*(i+1)+1+2*N]
        
        #along the y-axis
        #edge
        if (i+1) > (N_y - 1)*N_x:
            #particle
            #spin up
            Hamiltonian[2*i, 2*((i+1+N_x)%N-1)] = -1
            Hamiltonian[2*((i+1+N_x)%N-1), 2*i] = Hamiltonian[2*i, 2*((i+1+N_x)%N-1)]
            #spin down
            Hamiltonian[2*i+1, 2*((i+1+N_x)%N-1)+1] = -1
            Hamiltonian[2*((i+1+N_x)%N-1)+1, 2*i+1] = Hamiltonian[2*i+1, 2*((i+1+N_x)%N-1)+1]
            #hole
            #spin up
            Hamiltonian[2*i+2*N, 2*((i+1+N_x)%N-1)+2*N] = 1
            Hamiltonian[2*((i+1+N_x)%N-1)+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*((i+1+N_x)%N-1)+2*N]
            #spin down
            Hamiltonian[2*i+1+2*N, 2*((i+1+N_x)%N-1)+1+2*N] = 1
            Hamiltonian[2*((i+1+N_x)%N-1)+1+2*N, 2*i+1+2*N] = Hamiltonian[2*i+1+2*N, 2*((i+1+N_x)%N-1)+1+2*N]
        else:
            #particle
            #spin up
            Hamiltonian[2*i, 2*(i+N_x)] = -1
            Hamiltonian[2*(i+N_x), 2*i] = Hamiltonian[2*i, 2*(i+N_x)]
            #spin down
            Hamiltonian[2*i+1, 2*(i+N_x)+1] = -1
            Hamiltonian[2*(i+N_x)+1, 2*i+1] = Hamiltonian[2*i+1, 2*(i+N_x)+1]
            #hole
            #spin up
            Hamiltonian[2*i+2*N, 2*(i+N_x)+2*N] = 1
            Hamiltonian[2*(i+N_x)+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*(i+N_x)+2*N]
            #spin down
            Hamiltonian[2*i+1+2*N, 2*(i+N_x)+1+2*N] = 1
            Hamiltonian[2*(i+N_x)+1+2*N, 2*i+1+2*N] = Hamiltonian[2*i+1+2*N, 2*(i+N_x)+1+2*N]


    #spin-orbit coupling
    for i in range(N):
        #along the x-axis
        #edge
        if (i+1) % N_x == 0:
            #particle
            #hat{x} = (-1,0)
            Hamiltonian[2*i, 2*(i+1-N_x)+1] = V
            Hamiltonian[2*(i+1-N_x)+1, 2*i] = Hamiltonian[2*i, 2*(i+1-N_x)+1]
            #hat{x} = (1,0)
            Hamiltonian[2*i+1, 2*(i+1-N_x)] = -V
            Hamiltonian[2*(i+1-N_x), 2*i+1] = Hamiltonian[2*i+1, 2*(i+1-N_x)]
            #hole
            #hat{x} = (-1,0)
            Hamiltonian[2*i+2*N, 2*(i+1-N_x)+1+2*N] = -V
            Hamiltonian[2*(i+1-N_x)+1+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*(i+1-N_x)+1+2*N]
            #hat{x} = (1,0)
            Hamiltonian[2*i+1+2*N, 2*(i+1-N_x)+2*N] = V
            Hamiltonian[2*(i+1-N_x)+2*N, 2*i+1+2*N] = Hamiltonian[2*i+1+2*N, 2*(i+1-N_x)+2*N]
        #bulk
        else:
            #particle
            #hat{x} = (-1,0)
            Hamiltonian[2*i, 2*(i+1)+1] = V
            Hamiltonian[2*(i+1)+1, 2*i] = Hamiltonian[2*i, 2*(i+1)+1]
            #hat{x} = (1,0)
            Hamiltonian[2*i+1, 2*(i+1)] = -V
            Hamiltonian[2*(i+1), 2*i+1] = Hamiltonian[2*i+1, 2*(i+1)]
            #hole
            #hat{x} = (-1,0)
            Hamiltonian[2*i+2*N, 2*(i+1)+1+2*N] = -V
            Hamiltonian[2*(i+1)+1+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*(i+1)+1+2*N]
            #hat{x} = (1,0)
            Hamiltonian[2*i+1+2*N, 2*(i+1)+2*N] = V
            Hamiltonian[2*(i+1)+2*N, 2*i+1+2*N] = Hamiltonian[2*i+1+2*N, 2*(i+1)+2*N]

        #along the y-axis
        #edge
        if (i+1) > (N_y - 1)*N_x:
            #particle
            #hat{y} = (0,-1)
            Hamiltonian[2*i, 2*((i+1+N_x)%N-1)+1] = V*1j
            Hamiltonian[2*((i+1+N_x)%N-1)+1, 2*i] = Hamiltonian[2*i, 2*((i+1+N_x)%N-1)+1].conjugate()
            #hat{y} = (0,1)
            Hamiltonian[2*((i+1+N_x)%N-1), 2*i+1] = V*1j
            Hamiltonian[2*i+1, 2*((i+1+N_x)%N-1)] = Hamiltonian[2*((i+1+N_x)%N-1), 2*i+1].conjugate()
            #hole
            #hat{y} = (0,-1)
            Hamiltonian[2*i+2*N, 2*((i+1+N_x)%N-1)+1+2*N] = -V*1j.conjugate()
            Hamiltonian[2*((i+1+N_x)%N-1)+1+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*((i+1+N_x)%N-1)+1+2*N].conjugate()
            #hat{y} = (0,1)
            Hamiltonian[2*((i+1+N_x)%N-1)+2*N, 2*i+1+2*N] = -V*1j.conjugate()
            Hamiltonian[2*i+1+2*N, 2*((i+1+N_x)%N-1)+2*N] = Hamiltonian[2*((i+1+N_x)%N-1)+2*N, 2*i+1+2*N].conjugate()
        else:
            #particle
            #hat{y} = (0,-1)
            Hamiltonian[2*i, 2*(i+N_x)+1] = -V*1j
            Hamiltonian[2*(i+N_x)+1, 2*i] = Hamiltonian[2*i, 2*(i+N_x)+1].conjugate()
            #hat{y} = (0,1)
            Hamiltonian[2*(i+N_x), 2*i+1] = V*1j
            Hamiltonian[2*i+1, 2*(i+N_x)] = Hamiltonian[2*(i+N_x), 2*i+1].conjugate()
            #hole
            #hat{y} = (0,-1)
            Hamiltonian[2*i+2*N, 2*(i+N_x)+1+2*N] = V*1j.conjugate()
            Hamiltonian[2*(i+N_x)+1+2*N, 2*i+2*N] = Hamiltonian[2*i+2*N, 2*(i+N_x)+1+2*N].conjugate()
            #hat{y} = (0,1)
            Hamiltonian[2*(i+N_x)+2*N, 2*i+1+2*N] = -V*1j.conjugate()
            Hamiltonian[2*i+1+2*N, 2*(i+N_x)+2*N] = Hamiltonian[2*(i+N_x)+2*N, 2*i+1+2*N].conjugate()
        

    #pair potential
    for i in range(N):
        #particle
        Hamiltonian[2*i+1, 2*i+2*N] = Delta
        Hamiltonian[2*i+2*N, 2*i+1] = Hamiltonian[2*i+1, 2*i+2*N].conjugate()
        #hole
        Hamiltonian[2*i, 2*i+1+2*N] = -Delta
        Hamiltonian[2*i+1+2*N, 2*i] = Hamiltonian[2*i, 2*i+1+2*N].conjugate()

    return Hamiltonian

def get_Bott_index(N_x, N_y, eigm):
    N = N_x*N_y
    #position operators
    X = np.zeros((4*N, 4*N), dtype=complex)
    Y = np.zeros((4*N, 4*N), dtype=complex)
    for i in range(N):
        #particle
        #spin up
        X[2*i, 2*i] = (i % N_x)/N_x
        Y[2*i, 2*i] = (np.floor(i/N_x))/N_y
        #spin down
        X[2*i+1, 2*i+1] = (i % N_x)/N_x
        Y[2*i+1, 2*i+1] = (np.floor(i/N_x))/N_y
        #hole
        #spin up
        X[2*i+2*N, 2*i+2*N] = (i % N_x)/N_x
        Y[2*i+2*N, 2*i+2*N] = (np.floor(i/N_x))/N_y
        #spin down
        X[2*i+1+2*N, 2*i+1+2*N] = (i % N_x)/N_x
        Y[2*i+1+2*N, 2*i+1+2*N] = (np.floor(i/N_x))/N_y

    for i in range(4*N):
        X[i,i] = np.exp(2*np.pi*(X[i,i])*1j)
        Y[i,i] = np.exp(2*np.pi*(Y[i,i])*1j)

    #projection operators
    P = np.zeros((4*N, 4*N), dtype=complex)
    Q = np.zeros((4*N, 4*N), dtype=complex)
    for i in range(2*N):
        P[i,i] = 1
        Q[4*N-1-i, 4*N-1-i] = 1

    #projected position operators
    U_X = P @ np.conjugate(eigm.T) @ X @ eigm @ P + Q
    U_Y = P @ np.conjugate(eigm.T) @ Y @ eigm @ P + Q
    U_prod = U_Y @ U_X @ np.conjugate(U_Y.T) @ np.conjugate(U_X.T)

    #Bott index
    B = 0
    eigval, _ = LA.eig(U_prod)
    for i in range(4*N):
        B = B + imag(np.log(eigval[i]/abs(eigval[i])))

    B = B/(2*np.pi)
    return B

def main(N_x, N_y, mu, h_z, V, Delta):
    Hamiltonian = make_BdG_Hamiltonian(N_x, N_y, mu, h_z, V, Delta)

    _, eigvec = LA.eigh(Hamiltonian)
    
    B = get_Bott_index(N_x, N_y, eigvec)

    return B
