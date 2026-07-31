---
id: THM-3003
title: "Antipodal circuit rigidity and the multipole spread criterion"
status: PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428
depends_on:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3001-newton-circuit-reversal-involution-and-two-end-curvature-law
related:
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return
script: 04-computation/gmc_antipodal_circuit_rigidity_and_multipole_spread_criterion_thm3003.py
output: 05-knowledge/results/gmc_antipodal_circuit_rigidity_and_multipole_spread_criterion_thm3003.out
script_sha256: 594ab6c509ff42b61516ff47526e1811c33cf433e400b4d95d467a8932a9697e
output_sha256: 9571dfcf261bd1ed51e3675071101cb4217bc3de1c1c4167727422292d6fd202
hash_basis: LF-normalized bytes
---

# THM-3003 -- antipodal circuit rigidity and the multipole spread criterion

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.**

Three layers about one involution.  Notation is THM-3000 section 1 throughout;
`c_k=log(R_k/R_(k-1))=-Delta^3(log h)_(k-2)` is the **circuit**.

## 1. Rigidity: the circuit detects reciprocal symmetry exactly

**Theorem.**  Let `N` have degree `d` and all coefficients positive.  Then

    R_k=R_(d-k) for every 1<=k<=d-1
      <=> {r_i}={mu/r_i} as multisets, with mu=e_d^(2/d),
      <=> the empirical measure of log r is symmetric about its mean.  (1)

THM-3001 section 5 proved `(<=)`.  Here is `(=>)`, so (1) is an **iff** and the
Newton circuit is a *complete* detector of the reversal symmetry.

*Proof.*  `R_k=exp(-Delta^2(log h)_(k-1))`, so the hypothesis says the sequence
`phi_j=Delta^2(log h)_j`, `0<=j<=d-2`, is a **palindrome**.  Since
`log h_k=log e_k-log binom(d,k)` and

    Delta^2(log binom(d,.))_j=log[binom(d,j+2)binom(d,j)/binom(d,j+1)^2]

is already palindromic in `j` (binomials are symmetric), the hypothesis is
equivalent to `Delta^2(log e)` being a palindrome.  Two sequences with equal
second differences differ by an affine function of the index, so applying this
to `log e_k` and `log e_(d-k)` gives constants `A,B` with

    log e_k-log e_(d-k)=A+Bk,  i.e.  e_k=C mu^k e_(d-k),
    C=e^A,  mu=e^B.                                     (2)

At `k=0`, `1=C e_d`, so `C=1/e_d`; at `k=d`, `e_d=C mu^d`, so `mu^d=e_d^2`.
With `Q(t)=sum_k e_k t^k=prod_i(1+r_i t)`, (2) reads

    t^d Q(1/t)=e_d Q(t/mu).                              (3)

Both sides are degree-`d` polynomials.  The left vanishes exactly at `t=-r_i`,
the right exactly at `t=-mu/r_i`, so the multisets agree: `{r_i}={mu/r_i}`.
Taking logarithms, `log r` is symmetric about `(1/2)log mu=(1/d)sum log r_i`.
QED.

Equivalence of the two natural phrasings is elementary and used below: the
third-difference antipalindrome `c_k=-c_(d+1-k)` and the second-difference
palindrome `R_k=R_(d-k)` are the **same** condition, because
`Delta(phi)_j=-Delta(phi)_(n-1-j)` forces `phi_j-phi_(n-j)` constant and
odd, hence zero.

**Control.**  An exhaustive search over all multisets of size `d=4,5,6` drawn
from a ten-element rational pool (including reciprocal pairs and scalings) finds
**zero** ratio-palindromic multisets that are not reciprocal-closed.  A separate
`d=6` numerical solve of `c_2+c_5=0, c_3+c_4=0` from `400` random starts returned
`399` solutions, **all** log-symmetric.

## 2. Reversal is the antipodal map; Brouwer and Borsuk--Ulam

Work in centered log-root coordinates `ell_i=log r_i-(1/d)sum_j log r_j`, which
live in `V={ell: sum ell_i=0}`.  Because `R_k` is invariant under root scaling,
`R` is a function of `ell` alone, and by THM-3001 (2) reversal acts by

    ell -> -ell:   the ANTIPODAL MAP on V.                (4)

The circuit map `c:V->R^(d-2)` is equivariant,

    c(-ell)=-reverse(c(ell)),                             (5)

verified to `10^(-15)` for `d=5,6,7,8,12`.  Splitting `c=Sc+Ac` into
index-symmetric and index-antisymmetric parts, (5) says

    Sc is ODD in ell,   Ac is EVEN in ell.                (6)

**(2a) Balanced-locus crossing (IVT/Brouwer form).**  For any
reversal-symmetric weight `w` (`w_k=w_(d+1-k)`), `Phi_w(ell)=sum_k w_k c_k(ell)`
is odd by (6).  Hence in any **connected** reversal-closed class `H`, every path
from `N` to `N*` contains a point with `Phi_w=0`.  With `w=1`,
`Phi_1=log(R_(d-1)/R_1)`, so:

> every path from a polynomial to its reversal crosses the **balanced locus**
> `R_(d-1)=R_1`.

This is strictly stronger than THM-3001 section 2: it needs **no** "H implies
no-return" hypothesis and it produces an explicit witness inside `H`.  If in
addition `H` forces `c>=0` and `w>0`, the witness has `c=0` identically, hence
by (1) has log-symmetric roots and constant ratio sequence.

The straight path `ell(t)=(1-2t)ell` is a *degenerate* witness (it passes through
`ell=0`, all roots equal).  The output therefore uses **great-circle** paths
`ell(t)=cos(pi t)ell+sin(pi t)ell_perp` of exactly constant norm: for
`d=6,9,14` the endpoints have opposite `Phi_1` and the crossing occurs in the
interior, away from any degeneracy (`d=9` has three crossings).

**(2b) Essentiality (Borsuk--Ulam form).**  Restricted to the direction sphere
`S(V)=S^(d-2)`, `Sc` is a continuous **odd** map into `R^(ceil((d-2)/2))`.  Since
`ceil((d-2)/2)<=d-2` for `d>=4`, Borsuk--Ulam gives a zero, and the zero set has
dimension at least `floor((d-2)/2)`.  By section 1 that zero set **is** the
log-symmetric locus -- whose dimension in `S(V)` is `floor(d/2)-1`, equal to
`floor((d-2)/2)`.  So the two agree exactly:

> the log-symmetric locus is precisely the Borsuk--Ulam zero set of an odd map,
> hence **topologically essential**: no equivariant deformation removes it.

The failure of global no-return on log-symmetric profiles is therefore a stable
topological fact, not a feature of one convenient family.

## 3. The multipole dictionary and the spread criterion

The jet expansion of THM-3000 (2),

    log(N(n)/(a_d n^d))=sum_(j>=1)(ell_j/j)n^(-j),  ell_j=(-1)^(j-1)p_j,

is literally the **two-dimensional multipole expansion** of the potential
`sum_m log(n+r_m)` about the origin, with unit charges at the roots.  Hence the
dictionary:

| Newton-circuit object | potential-theory object |
|---|---|
| jet `ell_j` | multipole moment of the root charge |
| `u=ell_1` | dipole (centre of charge) |
| `x,z,w` | normalized quadrupole, octupole, hexadecapole |
| jets additive over `N=FG` | superposition / M2M |
| THM-2997 (21) wall subtraction | multipole subtraction |
| reversal `N->N*` | Kelvin transform `r->1/r` |
| THM-3001 two-end law | multipole/local (far/near field) duality |
| cumulant form of the curvature | translation (M2M) gauge |
| bounded normalized jets | FMM well-separatedness |

Additivity `p_j(FG)=p_j(F)+p_j(G)` is checked exactly; it is why the repo's
wall-stripping bookkeeping works at all, and it identifies THM-2997's `w_1,w_2,w_3`
as the **wall's** multipole moments.

**Spread criterion (new sufficient condition).**  Let

    kappa=max_i|r_i|/(p_1/d)                              (7)

be the root **spread ratio** (max modulus over mean).  From `|p_j|<=d max|r|^j`,

    m_j/m_1^j<=kappa^j,                                   (8)

so THM-3000's graded hypothesis `m_j/m_1^j=o(d^(j-3))` for all `3<=j<=k+1`
follows from the single geometric bound

    kappa=o(d^(1-3/(k+1))).                               (9)

The binding index is always `j=4`, so edge `3` needs only `kappa=o(d^(1/4))`,
edge `4` needs `kappa=o(d^(2/5))`, and so on.  **(9) replaces a jet-by-jet
invoice by one number**, and it is exactly the FMM aspect-ratio condition.
Verified against the exact bound (8) on uniform, two-cluster and heavy-tail
families at `d=60`.

**First-gap consequence.**  THM-2997 (24) gives `d/M->62/3` and `u/M^2->131/12`,
so the mean root tends to `0.528226 M`.  Then (9) at `k=3` reads

    max|root of N_M|=o(1.1263 M^(5/4)).                   (10)

> A root-modulus bound `|r|=o(M^(5/4))` on the wall-stripped core CLOSES the
> third edge.  No `P_4`, no new Macaulay chart, no new wall.

This supersedes THM-3000 section 7's `|p_4|=O(M^5)` phrasing as the cheapest
route, because (10) is a statement about root location rather than about a
resultant jet, and root-location bounds (Cauchy, Kojima, Fujiwara) are available
from coefficient data alone.

## 4. Boundaries and what is not claimed

- Section 1 is exact for every `d` and needs only positivity of the
  coefficients; it does **not** need real-rootedness.
- Section 2a needs `H` connected and the path to lie in `H`; it says nothing
  about disconnected classes.  Section 2b is a statement about the zero set of
  `Sc`, not about `c` itself: `Ac` is unconstrained there.
- (9) is **sufficient, not necessary**.  It is lossy: (8) uses the crude
  `|p_j|<=d max|r|^j`, which is tight only when the roots are concentrated at
  the maximum modulus.  A family with one heavy root has `kappa` large while its
  jets may still be tame; do not read a failure of (9) as a failure of THM-3000.
- Nothing here proves the THM-3001 section 6 classifier, which still needs a
  discrete-unimodality interior theorem, nor GMC(2), ULC, or any all-width family
  conclusion; THM-2989/2997's continuation wall invoice is untouched.
- (10) is a *conditional reduction*: no root-modulus bound for the core is proved
  here.

## 5. Reproduction

    python3 04-computation/gmc_antipodal_circuit_rigidity_and_multipole_spread_criterion_thm3003.py

Three layers, all reporting `True`: rigidity (with the exhaustive rational
control and the `e_k e_d=mu^k e_(d-k)` identity check), antipodal equivariance
with non-degenerate great-circle crossings and the Borsuk--Ulam dimension table,
and the multipole dictionary with the exact spread bound (8).
