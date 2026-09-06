# Independent audit of the conserved-amplitude orbit

**Status: INDEPENDENT ANALYTIC AND EXACT SOURCE AUDIT — PASS.**
Root independently audited [the proof](long_frontier2_sep06_duality.md),
[source](../../04-computation/long_frontier2_sep06_duality.py), and
[output](long_frontier2_sep06_duality.out). The degree-support hypothesis
was made explicit during review; formulas and source bytes were unchanged.

The regular Phi and Psi rows use rising factorials, retaining all interior
amplitudes and genuine endpoint coefficients. Direct cancellation of
consecutive products proves both parameter steps for every height. The
doubled integral uses u^(2z)(1-u) and argument u^2 t, with strictly positive
denominator factors. Its absence of an inverse channel comes from the full
regular count fibre, not deletion of the old carry.

For a family supported in degrees0 through2H, each coefficient evolves
by the same nonzero scalar as its genuine counterpart. Every ratio is
therefore constant on the integer orbit, and the converse is coefficientwise.
The support assumption is now explicit; higher-degree terms would require
an additional description. The real-parameter H1 formula is an explicit
extension, not an interpolation uniquely implied by the discrete recurrence.

At H1 the discriminant is A_z^2(lambda^2-R_z/(10A_z))/9. Since
A_z-R_z=8z+10 and R_z/A_z tends to one, its exact global threshold is
lambda>=1/sqrt10. At equality the discriminant is strictly positive for
every finite parameter. Positive constant and middle coefficients then
make both roots negative. Below the threshold the discriminant is
eventually negative, including along integer parameters.

Root independently reconstructed the ell_z sign evaluation, gap,
consecutive difference, lambda8/25 transition polynomial and paired
residual using SymPy, separately from the producer's Fraction polynomials.
The positive gap13/40-ell_z tends to zero. Thus strict same-root negativity
for every z>=1 is equivalent to lambda>=13/40, also including equality.
Neither boundary threshold is confused with a finite double root or zero.

The polynomial z^2-69z-112 has one positive root between70 and71 and one
negative root, proving the complete integer sign transition. The exact
attenuation10981/34320 creates the named same-zero model while remaining
above the global reality threshold. The displayed supports and masses are
actual first-return addresses, but their attenuated response coefficients
are models. Changing a Laurent phase cannot create that normalized
interior attenuation. No genuine Laurent counterexample follows.

The paired factor divides the characteristic coefficient for every
attenuation. Its residual has nonnegative coefficients with positive
constant and linear terms exactly in the stated calibrated region. This
is a complete H1 condition, not an assertion for higher heights. Conserved
ratios explain why more parameter steps cannot remove the defect.

The literal multinomial controls and h3,r2 complementary reflection retain
their full normalizing factors. The inherited singular h1,x=-1 example
correctly separates generic quotient -1/90720 from the raw specialized
row t^2/720. No argument identifies the zero block with a nonzero
complementary root.

Root replayed all93 always-active gates in normal and optimized Python.
Both raw byte streams agree with the frozen output. The source imports
no repository producer. Five additional independent rational identities
reconstruct the key sign, threshold, transition and residual formulas.

    source SHA256 e3d9f21f61d1a5999f02ff83f9e3f6cd565e23da094278dc0f2cfbfd6297f274
    output SHA256 d3427c911cffd4e919577ae715caef7e1ebb12091c4176892924880bd58f6020

Reproduction:

    python3 04-computation/long_frontier2_sep06_duality.py
    python3 -O 04-computation/long_frontier2_sep06_duality.py

The conserved-amplitude theorem, global model obstruction and sharp H1
thresholds pass. Actual all-height separation and higher-height amplitude
regions remain OPEN.
