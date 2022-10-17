# failed-calculation-of-Bott-index
I made a fortran program to calculate Bott index of (topological) superconductor. But it is not correct.


Here is what we are trying to do.

We consider the $L \times L$ square lattice system with periodic boundary condition.

Hamiltonian of superconductivity with perpendicular magnetic field, h_z, and Rashba spin-orbit coupling, V:

$\mathcal{H} = -t\sum_{\langle i,j\rangle,\sigma} c_{i\sigma}^\dagger c_{j\sigma} - \sum_{i} (\mu + h_z\sigma) c_{i\sigma}^\dagger c_{i\sigma} + \Delta\sum_i (c_{i\uparrow}^{\dagger}c_{i\downarrow}^\dagger + \mathrm{H.c.}) -V \sum_i [ (c_{i-\hat{x}\downarrow}^\dagger c_{i\uparrow} - c_{i+\hat{x}\downarrow}^\dagger c_{i\uparrow}) + i(c_{i-\hat{y}\downarrow}^\dagger c_{i\uparrow} - c_{i+\hat{y}\downarrow}^\dagger c_{i\uparrow}) + \mathrm{H.c.} ]$.

We can transform this Hamiltonian into BdG(Bogoliubov-de Gennens) form, using commutation relation of creation and annhilation operators:

$\mathcal{H} = ({}^t\boldsymbol{c}^\dagger \ \ {}^t\boldsymbol{c}) \mathcal{H}_{\mathcal{BdG}} {}^t({}^t\boldsymbol{c} \ \ {}^t\boldsymbol{c}^\dagger) + \mathrm{const.}$

Diagonalizing the "BdG Hamiltonian", $\mathcal{H}_{\mathrm{BdG}}$, we can get eigenenergies, occupied and excitation energy levels of Bogoliubov quasiparticles, and eigenvectors.

Then, we make these matrices:

$U_X = Pe^{i2\pi X}P + Q, U_Y = Pe^{i2\pi Y}P + Q$,

where $P$ is the projector to the occupied energies, $Q$ is the projector to the excitation energies, and $X(Y)$ is the $x(y)$ position operator rescaled to $[0,1]$.

Bott index, $B$, is defined

$B = \frac{1}{2\pi} \mathrm{Imtrln}(U_Y U_X U_Y^\dagger U_X^\dagger)$.

The matrix logarithm is calculated using this lemma:

$\mathrm{trln}(U_Y U_X U_Y^\dagger U_X^\dagger) = \sum_i \mathrm{ln}\frac{\lambda_i}{|\lambda_i|}$,

where $\lambda_i$ are the eigenvalues of $U_Y U_X U_Y^\dagger U_X^\dagger$.

This index $B$ should be integer. But this program produces real numbers.

### (Added on October 16, 2022)

I realized there might be a mistake in the spin-orbit part in the subroutine, make_BdG_Hamiltonian.

### (Added on October 17, 2022)

I fixed the mistakes and turned it into python code and it worked correctly.

I'll make its Fortran version for larger matrices.
