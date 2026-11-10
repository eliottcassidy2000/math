# Heat covariance sequences and connected replica kernels

**Session 249, 10 November 2026. Status: PROVED IN NOTE** for the identities
and recurrence below, conditional only on the explicitly inherited forest
identity. **FINITE-EXACT** for the new bounded computations. No canon ID or
literature-priority claim. **OPEN:** an all-order higher-moment recurrence or
positivity theorem sufficient for Chan–Pak.

## 0. Inheritance and the useful change of object

Closest inherited mechanisms, accepted here rather than re-audited:

- [THM-4761, forest/Walsh chaos](../../../../01-canon/theorems/THM-4761-random-tournaments-the-forest-identity-is-the-walsh-chaos-decomposition-of-kappa.md)
  supplies the expansion of the normalized Hamiltonian-path count into even
  linear forests.
- [THM-5127, heat flow](../../../../01-canon/theorems/THM-5127-heat-flow-reading-of-the-w-polynomial-tilted-ensemble-mean-and-the-mod-four-quadratic-form-law-of-a-vertex-cube.md)
  reads arc noise probabilistically.
- [S248 support spectrum](../../session248_recursive_symmetry_20261110/bridges/PROOF.md),
  sections 2 and 5, and [oscillator](../../session248_recursive_symmetry_20261110/OSCILLATOR.md)
  give the commuting vertex-support coordinate and birth/death ladders.
- [The October run-kernel compiler](../../oeis_transfer_atlas_20261014h.md),
  section 3, already pulls back the factorial-series ODE through a rational
  kernel. The ODE method and its seven-step unweighted specialization are
  **inherited**, not newly discovered here.

Canonical hostile: the unweighted no-unit-descent recurrence is not the
weighted second-moment recurrence. Corrected near miss: a scalar second
moment is not the Hamiltonian-path spectrum. Least-used sidecar: the number
of vertices in a Walsh forest, in addition to its number of edges.

Concept board: even-path components; support-size oscillator; rational
block kernels; connected replica overlaps; exact coefficient positivity.
The anchor is an across-size covariance recurrence, the niche is subcubic
three-replica graphs, and the wildcard is a finite-moment support certificate.

A000568 counts **unlabelled tournament types**. The probabilities in this note
are uniform over **labelled tournaments**, equivalently class weights
`n!/(2^binom(n,2)|Aut(T)|)`. The new covariance sequences are neither A000568
nor a recurrence for that census. A fixed-size flip-time recurrence is a third,
different object, separated in section 4.

## 1. An exactly solvable two-parameter covariance array

For `n >= 1`, write `kappa_n(T)=2^(n-1) H(T)/n!`. Let `N_rho` independently
reverse each arc with probability `(1-rho)/2`. Define `K_z` by retaining each
vertex independently with probability `z`, keeping all arcs between retained
vertices, and fairly resampling every other arc. Initially `0<=rho,z<=1`.
Both operators are diagonal in Walsh characters:

    N_rho K_z chi_F = rho^e(F) z^v(F) chi_F.                 (1)

This is the binomial mixture of S248's fixed-cardinality retention kernels.
It is also the time-t vertex refresh semigroup when `z=exp(-t)` and every
vertex refreshes at rate one. Put `y=rho^2` and

    F_n(y,z) = n! E[kappa_n(T) (N_sqrt(y) K_z kappa_n)(T)].  (2)

For formal generating purposes define `F_0=1` separately. Equation (2) is
only claimed for `n>=1`; there is no empty-tournament normalization trick.
The probabilistic definition extends uniquely to a polynomial in `y,z`.

**Theorem 1.** With `k=n-2j`,

    F_n(y,z) = n! + sum_(j>=1, 1<=p<=min(j,k))
                  k! binom(k,p) binom(j-1,p-1) 2^p y^j z^(2j+p).       (3)

In particular every coefficient is a nonnegative integer, and

    sum_(n>=0) F_n(y,z) x^n = sum_(k>=0) k! f(x;y,z)^k,
    f(x;y,z) = x + 2 y z^3 x^3 / (1-y z^2 x^2).                     (4)

These are **formal ordinary generating series**, not convergent analytic
series near a positive x. At `y=z=1`, (4) is exactly the inherited W_n kernel.

**Proof.** In the inherited forest expansion, a forest with `2j` edges and
`p` nontrivial components uses `2j+p` vertices. Its Walsh coefficient in
kappa is `2^p epsilon_F/(n)_(2j)`, where `epsilon_F=+/-1`. Distinct edge sets
are orthogonal, so the covariance is a sum of squared coefficients times
`y^j z^(2j+p)`. The number of such forests is

    (n)_(2j+p) binom(j-1,p-1)/(2^p p!).

Multiplication by `n!` gives (3), since
`n!(n)_(2j+p)/(n)_(2j)^2 = k!^2/(k-p)!`.
For (4), view the k components, including singleton vertices, as ordered
odd blocks. A nontrivial block of size `2a+1` carries weight
`2 y^a z^(2a+1)`. Choose p nontrivial blocks among k and compose j into p
positive parts. The additional `k!` is exactly the factor in (3). ∎

For example,

    F_6 = 720 + 192 y z^3 + 8 y^2 z^5 + 8 y^2 z^6,
    F_8 = 40320 + 8640 y z^3 + 192 y^2 z^5 + 576 y^2 z^6
           + 8 y^3 z^7 + 16 y^3 z^8.

Thus the equal-edge-degree modes `(e,v)=(4,5)` and `(4,6)`, indistinguishable
by arc heat alone, have separate exact coefficient sequences. Differentiating
in z at z=1 records support-weighted energies; differentiating in y records
half-edge-degree energies. This is a new use of the S248 sidecar, rather than
a new proof of its spectrum.

## 2. A seven-step recurrence across tournament size

Set `a=y z^2`, `b=(2z-1)a`, and `Y(x)=sum F_n x^n`. Then

    d=1-a x^2,  h=x+b x^3,  R=1+(a+3b)x^2-ab x^4,
    d h^2 Y' + (h-d) R Y + d R = 0.                       (5)

Indeed the factorial series `G(t)=sum k!t^k` satisfies
`t^2 G'+(t-1)G+1=0`, and `Y=G(h/d)` with `(h/d)'=R/d^2`.
This is the inherited rational-kernel compiler, now applied to the new
support-refined covariance kernel.

**Corollary 2.** With negative-index terms zero,

    F_n = n F_(n-1) - 3b F_(n-2)
        + ((2b-a)n+4a-2b) F_(n-3) + (a^2+4ab) F_(n-4)
        + ((b^2-2ab)n-2b^2+10ab) F_(n-5)
        - a^2 b F_(n-6) - ab^2(n-6) F_(n-7) + c_n,        (6)

where `c_n=[x^n](1+3b x^2-(a^2+4ab)x^4+a^2b x^6)`.
The recurrence is homogeneous from n=7. At `y=z=1` it specializes exactly
to the October atlas recurrence, including its forcing terms. At y=0 or z=0
it reduces to `F_n=n!`. At z=1/2, b=0, several terms disappear.

For each fixed rational pair `(y,z)`, this computes the first N values in
O(N) rational arithmetic operations and constant auxiliary streaming storage.
No minimal-order or bit-complexity claim is made. Computing all polynomial
coefficients has larger output size and is not an O(N) assertion.

## 2a. Coherent birth gives a second, simpler sequence compiler

There is a useful alternative to normalizing each chaos back to H. Couple
`T_1,T_2,...` by adjoining independent fair incidences. Let `N_(2j)(T_n)`
be the unnormalized signed even-forest sum of THM-4761, and set

    Z_n(t)=sum_(j>=0) N_(2j)(T_n) t^(2j),   N_0=1.

Every forest monomial touching a new vertex averages to zero. Therefore
`E[Z_m(t)|T_n]=Z_n(t)` for m>=n: this is a polynomial-valued birth martingale,
a direct use of S248's coherent unnormalized embedding coordinates. For
`a=(st)^2`, its two-time covariance depends only on the earlier time:

    E[Z_n(s) Z_m(t)] = B_n(a),                 m>=n,
    sum_(n>=0) B_n(a) x^n/n! = exp(x+2a x^3/(1-a x^2)).               (7a)

The proof squares the forest coefficients before the H-normalization:

    B_n(a)=1+sum_(j>=1,p>=1,2j+p<=n)
               2^p (n)_(2j+p) binom(j-1,p-1) a^j / p!.

This is the same connected even-path kernel, now with **exponential assembly**
instead of the factorial-series composition. Its EGF obeys

    (1-2a x^2+a^2 x^4) B'(x)=(1+4a x^2-a^2 x^4)B(x),

hence, with negative subscripts zero and B_0=1,

    B_(n+1)=B_n+2a(n)_2 B_(n-1)+4a(n)_2 B_(n-2)
                 -a^2(n)_4(B_(n-3)+B_(n-4)).                     (7b)

At a=1 its first terms, beginning at n=0, are
`1,1,1,13,49,361,3121,25621,296353,3176209,42700321`.
These are martingale-energy values, **not** moments of H without additional
factorial extraction and **not** tournament class counts. For each fixed j,
`E[N_(2j)(T_n)^2]` is a polynomial in n of degree at most 3j, so its
`(3j+1)`-st forward difference vanishes identically. This provides another
honest across-size recurrence, organized by fixed Fourier degree.

## 3. Higher moments: replace a path kernel by connected overlap graphs

The second moment works because two Walsh forests must coincide. For higher
moments the correct generalization is a **connected colored overlap kernel**,
not a second copy of the same rational function.

Let `F_1,...,F_r` be even-path linear forests on a common labelled vertex set;
empty layers are allowed. Require every undirected edge to occur in an even
number of layers. Write `e_i=|E(F_i)|`, `p_i` for its nontrivial components,
and give the tuple weight

    w(F_1,...,F_r) = product_i 2^p_i epsilon_(F_i).        (7)

Because each edge occurs evenly, this weight is invariant under relabelling:
all changes of the antisymmetric reference orientation cancel. It is also
multiplicative over connected components of the union graph.

Define `C_r(u;t_1,...,t_r)` by summing `w t_1^e_1...t_r^e_r u^v/v!`
over tuples on `[v]` whose union is connected. Include the isolated singleton
with all layers empty, so its first term is u. Let `A_r=exp(C_r)`.
For the coefficient of u^n define the linear extraction map

    L_n(t_1^e_1...t_r^e_r) = product_i (n-e_i)!.

**Theorem 3 (exact replica compiler).** For `n>=1`,

    E[H_n^r] = n! 2^(-r(n-1)) [u^n] L_n(exp(C_r)).       (8)

**Proof.** Expand each H by the inherited forest formula. Averaging the
product of Walsh characters kills exactly the tuples with an edge of odd
multiplicity. In the surviving terms, `n!/(n)_e_i=(n-e_i)!`.
The labelled exponential formula assembles the union's connected components,
with the sign retained; `[u^n]exp(C_r)` is the labelled tuple sum divided by
n!. This proves (8). Every layer is a forest, so `e_i<=n-1` and all factorials
are defined. ∎

For r=2, the two forests coincide and

    C_2 = u + 2 sum_(j>=1) u^(2j+1)(t_1 t_2)^(2j).

After the factorial extraction, this recovers the ordinary odd-block kernel
of section 1. So (8) is a genuine extension, not an unrelated moment notation.

### Three replicas expose a subcubic, signed graph problem

For r=3, each union edge belongs to exactly two layers and has one of the
colors `12,13,23`. Every union vertex has degree at most three. At a trivalent
vertex all three pair-colors occur once. The edges of each layer i must form
an even-path forest, rather than a cycle or an arbitrary bounded-degree graph.
This is the exact map from the moment problem to colored subcubic graphs.

The first terms are

    C_3 = u
        + 2u^3(t_1^2 t_2^2+t_1^2 t_3^2+t_2^2 t_3^2)
        + 8(u^3-u^4)t_1^2 t_2^2 t_3^2 + O(u^5).          (9)

The `+8u^3` term is the triangle with all three pair-colors. The `-8u^4`
term is the three-star with all three pair-colors. In a cyclic reference
orientation of the triangle each path product gives a positive total sign;
in an outward reference orientation of the star each two-edge path contributes
-1, giving a negative triple product. There are six pair-colorings in either
case, and four choices of the labelled star center; division by labelled
vertex factorials gives the displayed coefficients.

**Hostile control:** connected replica coefficients are not all positive,
even though the final moments are. Discarding their signs loses the true
moment recurrence. At n=3 the normalized third central moment is `2/9`; at
n=4 triangle and star contributions cancel and it is 0. At n=5 it is
`11/900`. These are computed directly in the runnable example, not asserted
from a fitted extrapolation. The output also gives every u^5 coefficient.

**OPEN next step:** obtain a finite kernel/recurrence for the signed subcubic
classes in C_3, first for tree union graphs, then with controlled cycle excess.
For the tree sector, a rooted edge interface can retain, for each layer, the
parity of the unfinished path and its orientation state. That supplies a
concrete finite-state grammar target, not a claimed completed grammar. A
finite-state result for the tree sector alone would not enumerate all of C_3.
No rationality or D-finiteness of the complete C_3 is claimed.

## 4. Fixed-size time recurrences are a separate, easier consequence

At a fixed n with `N=binom(n,2)>0`, let Q resample one uniformly chosen arc.
A degree-e Walsh character has eigenvalue `1-e/N`. Since H is unchanged by
reversing all arcs and has degree at most n-1, H^r has only even degrees,
at most `min(r(n-1),N)`. Therefore, for any starting T,

    a_k = (Q^k H^r)(T)

satisfies a constant-coefficient recurrence annihilated by

    product_(0<=j<=floor(min(r(n-1),N)/2))
                 (E - (1-2j/N)),                         (10)

where E shifts k forward. Some factors can be absent. Equivalently its OGF
in k is rational with denominator dividing the product of
`1-(1-2j/N)s`. A zero eigenvalue causes only a finite initial contribution.
This recurrence is in the **number of refresh steps k with n frozen**.
It gives no recurrence in n for A000568 or for E[H_n^r].
The two-parameter construction in sections 1–2 is what supplies a genuine
across-n recurrence for the second-moment covariance.

## 5. What would actually certify a Chan–Pak coefficient?

Let `Y_n=(H_n-1)/2`, an integer in `[0,M_n]`, where the safe bound is
`M_n=floor((n!-1)/2)`. Let

    beta_(n,r)=E binom(Y_n,r),
    P_n(q)=E q^Y_n = sum_(r=0)^M_n beta_(n,r)(q-1)^r.

Each beta is an explicit rational linear combination of the first r ordinary
H moments, so (8) is an exact input interface. For `0<=j<=M_n`,

    P(Y_n=j) = sum_(r=j)^M_n (-1)^(r-j) binom(r,j) beta_(n,r).         (11)

This gives an exact necessary-and-sufficient coefficient/support test. There
is also a finite-moment **lower bound**, which can succeed without recovering
the whole distribution. For every odd positive L, define beta beyond M_n as
zero and put

    B_(n,j,L) = sum_(r=j)^(j+L) (-1)^(r-j) binom(r,j) beta_(n,r).
    B_(n,j,L) <= P(Y_n=j).                                           (12)

Thus `B_(n,j,L)>0` certifies a tournament with `H=2j+1`.
For even L the inequality reverses. The proof is pointwise: for an integer m,

    sum_(r=j)^(j+L) (-1)^(r-j) binom(r,j) binom(m,r)
      = binom(m,j) sum_(s=0)^L (-1)^s binom(m-j,s).

When m>j and L<m-j the last sum is `(-1)^L binom(m-j-1,L)`;
otherwise it is the exact indicator contribution. This proves both bounds.
The smallest odd truncation is `beta_j-(j+1)beta_(j+1)`.

Since a sufficiently long odd truncation equals (11), for each fixed j,

    2j+1 occurs  <=>  there exist n and odd L with B_(n,j,L)>0.        (13)

Chan–Pak's existence assertion is precisely (13) for every j except 3 and 10;
the impossible values 7 and 21 are inherited exclusions, not proved here.
The missing step is a **uniform** construction/control of the required
higher-moment kernels or a positivity estimate for these alternating sums.
Second moments, near-Gaussian shape, and mod-4 balance do not discharge it.
The safe n! bound is intentionally crude; no claim of efficient all-j
support recovery is hidden in the equivalence.

## 6. Reproduction and scope ledger

Run from repository root:

    python3 05-knowledge/results/sequence_gf_20261110/heat/compile_heat_covariance.py

The standard-library program writes [output.json](output.json). Its checks
use explicit exceptions and remain active under Python -O. It checks:

- Formula (3) against recurrence (6) for six rational parameter pairs through
  n=40, including reset boundaries and a formal nonprobabilistic pair `(2,3)`.
- Formula (3) against a direct Hamiltonian-path DP plus Walsh transform on
  every labelled tournament of orders 2–5, with edge degree and vertex support
  retained separately. No class-census weights are guessed.
- Moment-to-support inversion and odd/even truncated Bonferroni inequalities
  against the same exact labelled histograms, including full-truncation equality.
- The coherent birth-energy recurrence (7b) against its finite forest sum
  through n=40 at a=1 and a=1/4.
- Connected C_3 coefficients through five vertices by two forest masks and
  their XOR third mask, keeping the sign and connectedness condition.

This is a cheap check of the newly derived objects, not an inherited-results
audit. The map to S248 preserves spectral covariance and both noise degrees;
it loses individual types, automorphism groups, and higher moments. The
replica kernel restores the higher-moment sidecar. The last coefficient test
restores exact support, at a currently uncontrolled computational cost.
