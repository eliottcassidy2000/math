# Score-scale grading and exact random-tournament controls

Status: PROVED finite/coefficientwise identities and elementary bounds. No
pointwise Hamiltonian asymptotic follows merely from these controls.

## 1. The score scale is a genuine expansion parameter

For a connected quotient of a collection of directed runs, let e=E−V and let
ell denote its number of degree-one vertices. Set M=S/n and D=max_i|(S1)_i|.
Summing each leaf label gives an entry of ±M1, bounded byD/n. The remaining
E−ell edges have magnitude1/n, and there are V−ell free labels. Therefore

    |hom(quotient,M)| <= D^ell n^(-e-ell).

The isolated single-edge graph contributes zero because1^T M1=0; loops also
vanish. For all other connected shapes the leaf vertices can be summed
independently. This statement is before multiplying by the finite combinatorial
Möbius/kernel coefficient.

For D=1 this gives weight w=e+ell. For D=O(n^alpha), the corresponding bound is
n^(-e-(1-alpha)ell), with a shape-dependent constant. In particular:

- path trees: exponent1−2alpha;
- a cycle with one terminal path: exponent1−alpha;
- vacuum two-loop cores: exponent1.

At the square-root score scale, path trees cease to be a small1/n correction.
This is a power-counting fact, not a claim that every upper bound is attained,
and not a uniform summation theorem. It explains why the near-regular expansion
must not be transplanted unchanged to arbitrary random tournaments.

## 2. Uniform independent orientations: exact expectations

Let every unordered edge be oriented independently with probability1/2. Every
nonempty signed linear forest is a product of distinct independent mean-zero
edge signs. Hence E N_L=0 for L>0 and

    E Phi_T(z)=1,       E kappa(T)=1

exactly for every n. The latter also follows by summing the probability2^(-(n−1))
for each of the n! vertex orderings. This is an average identity, not pointwise
concentration and not a regular-score ensemble identity.

Writing d=S1 and W=||d||^2, distinct two-edge monomials are orthogonal under
this measure. Expansion gives

    E W=n(n−1),
    Var W=2n(n−1)(n−2).

Consequently W/n² tends to1 in mean square. The exact second forest coefficient
N2=n(n−1)−W has mean0 and variance2n(n−1)(n−2). By contrast every almost-regular
tournament has W=n and [z²]logPhi=1−2/n.

## 3. The spectral normalizer stays nontrivial

For any tournament let the eigenvalues of M=S/n be ix_j. Put b=(n−1)/(2n).
The eigenpair/Frobenius bound gives0≤x_j²≤b and sum_j x_j²=2b. Since
atanh(u)/u is increasing,

    (n−1)/n <= log rho = sum_j atanh(x_j²)
                    <= 2atanh(b),

so

    exp((n−1)/n) <= rho <= (3n−1)/(n+1) < 3.

Under independent orientations,

    E tr(M^4) = (n−1)(2n−3)/n^3.

Indeed (S²)_ii=−(n−1), while E[(S²)_ij²]=n−2 for i≠j. For0≤u≤1/2,
atanh(u)−u ≤ u³/[3(1−u²)] ≤ (2/9)u². Thus

    0 <= E[log rho]−(n−1)/n
       <= (2/9)(n−1)(2n−3)/n^3.

It follows that log rho→1 inL1 and rho→e inL1, since rho is uniformly bounded
by3. Yet E kappa=1 exactly. Therefore a spectrum-only approximation kappa≈rho
cannot hold with vanishing mean absolute error over the unconstrained ensemble:

    E|kappa−rho| >= |E kappa−E rho| -> e−1.

This is a proved obstruction to extending the regular predictor without score/
endpoint corrections. It does NOT prove kappa concentrates at1, does NOT prove
E[kappa/rho]→1/e, and does NOT assume independence of kappa and rho.

## 4. Interpretation and limits

The almost-regular regime suppresses degree-one ends by1/n. A random tournament
has aggregate score energy of ordern² instead ofn, and its nonzero tree level
is essential. The coincidence of a formal tree penalty near−1 with a spectral
log gain near+1 is an informative heuristic cancellation; the average identities
above are the precise statements proved here. No Collatz or lattice implication
is inferred from this shared quadratic-energy language.
