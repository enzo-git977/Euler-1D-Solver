-- INSTRUCTIONS --------------------------------------------------------

COMPILE :
make

RUN / 
./solve

CLEAN DATA : 
make clean

-- GENERAL INFO ABOUT THE CODE  --------------------------------------------------------

First-order Roe: 
++ It is extremely robust and provides a monotone solution (no oscillations)
-- Suffers from heavy numerical diffusion, and without a fix, it can produce non-physical solutions

Harten Entropy Fix (delta_star) : smooths the absolute value of the eigenvalues ($\lambda$) of the Roe matrix when they are close to zero.
++ better accuracy of the rarefaction waves, mandatory for testcase 3
-- introduce numerical diffusion

MUSCL Reconstruction: Achieving Second-Order or even third-order Precision
Parameter PHI: This parameter controls the centering of the reconstruction. 
Phi = -1 (fully upwind), which is the most stable choice for supersonic flows, Phi = 1/3 ; third order, Phi = 1 2nd order centered (can be non stable)

Parameter $b$ (Compression): This factor multiplies the neighboring slopes. A value near 1 ensures stability, 
while values closer to 2 compress shocks to be very sharp, though they may introduce small dispersive oscillations behind the front.

Limiter : 
In smooth regions, ϕ -> 1 which gives second-order accuracy
Near shocks, ϕ -> 0 which falls back to first-order Roe


-- Robustess test (cas= 3)  --------------------------------------------------------
The third testcase is designed to assess the numerical robustness of the Riemann solver
It consists of two symmetric rarefaction waves moving away from each other, creating a near-vacuum state at the center of the domain ($x=0.5$).
Positivity Preservation: This is a "crash test" for the Roe solver, which is not naturally positivity-preserving.
The strong expansion can lead to non-physical negative pressures at the expansion fan's center.
Entropy Fix Requirement: STANDARD ROE SCHEME WILL FAIL HERE. 
A significant Harten-type entropy fix (delta_star) is mandatory to add enough numerical diffusion to maintain a positive pressure.
INDICATION : delta_star = 20, cfl = 0.1