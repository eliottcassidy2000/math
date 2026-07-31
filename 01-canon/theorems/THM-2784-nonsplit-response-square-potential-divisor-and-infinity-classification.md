---
id: THM-2784
title: "Nonsplit response square-potential divisor and infinity classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On a genuine
  quadratic deck U^2=V, an anti-invariant response
  R'=kappa/U is equivalent to the base square-potential equation
  V(F')^2=4kappa^2F for F=R^2.  Its complete finite/infinity divisor ledger
  forbids double roots of V, forces every squarefree V to be linear, and
  places every non-pure-power survivor in an explicit clean three-value
  permutation passport; [F]=[V] retains the quadratic cover.  This is a
  uniform response-layer sieve, not chart entry, full nonsplit closure,
  JC(2), or DC(2).
source: jc-chart-entry/root-nonsplit-square-potential-2026-07-28
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2158-quartic-quadratic-deck-parity-and-exact-finite-pole-criterion
  - THM-2189-nonsplit-quartic-deck-forces-the-remaining-pole-congruence
  - THM-2202-uniform-all-degree-quartic-pole-closure
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
related:
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2217-square-prefix-pole-alternative-and-odd-leading-degree-terminal-wall
  - THM-2778-all-degree-complete-chosen-sheet-split-exact-prefix-closure
  - THM-2781-terminal-tail-perfect-power-rigidity-and-response-count
script: 04-computation/jc2_nonsplit_response_square_potential_thm2784.py
output: 05-knowledge/results/jc2_nonsplit_response_square_potential_thm2784.out
script_sha256: 7fcbe78730b20ffe06dcc16bf640600ca5810631789e9b6d74eda113f7b053d1
output_sha256: 07922e330700c495d0f2f541727f5a4383b17dd37c03151ea6658c111f5b5fee
independent_script: 04-computation/jc2_nonsplit_response_square_potential_hostile_audit_thm2784.py
independent_output: 05-knowledge/results/jc2_nonsplit_response_square_potential_hostile_audit_thm2784.out
independent_script_sha256: 77f32b028ac8276c1fa292202c3e82ab02c3618a49207eb81b62f58a4362c136
independent_output_sha256: 3e6fad51f117d612ca8de5a08d0e62f22ecc693b2eb7fd7e276ad622b1667b49
hash_basis: LF-normalized bytes
---

# THM-2784 -- nonsplit response square-potential classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

**Exact scope:** the response layer of the polynomial exact-square-prefix,
genuine nonsplit terminal quartic chart.  The theorem is uniform in the
reduced Faber degree.  It does not derive that chart from an arbitrary planar
Keller pair, close the surviving linear-leading-coefficient case, or prove
`JC(2)`/`DC(2)`.

## 1. Inheritance and the newly retained coordinate

The closest proved mechanism is the all-degree nonsplit response identity in
THM-2214, Section 4 (and its degree-14/18 descendants):

```text
Phi_Q=0,            Psi_Q in C,            R_Q'=kappa/U,
U^2=V,              kappa in C*.
```

THM-2202 supplies polynomiality of the quadratic approximate root, and
THM-2230 makes the full Faber representative an exact target-shear gauge.
THM-2217 is the corrected near miss: its first-tooth divisibility gives only
`R deg(q)>=deg(V)` when `deg(V)` is odd, not terminal descent.

The least-used sidecar is the quadratic deck itself.  The split argument sees
a rational primitive on the base.  On the genuine nonsplit deck, the response
is instead an anti-invariant primitive.  Squaring that primitive does not lose
the differential equation: it produces one exact rational function on the
base.

The new coordinate is

```text
F=R_Q^2 in C(x).
```

It simultaneously records the anti-invariant de Rham class and the complete
finite ramification ledger of the response.

## 2. Function-field theorem

Let `C` be algebraically closed of characteristic zero,

```text
K=C(x),             V in C[x]\{0},
L=K(U),             U^2=V,
```

and assume `V` is nonsquare in `K`.  Let `sigma` be the nontrivial deck
involution.  Fix `kappa in C*`.

Suppose `R in L` satisfies

```text
sigma(R)=-R,                  R'=kappa/U.             (1)
```

Then `F=R^2` belongs to `K`, is nonconstant, and satisfies

```text
boxed: V(F')^2=4 kappa^2 F.                           (2)
```

Conversely, if a nonconstant `F in K` satisfies `(2)`, then

```text
R=U F'/(2kappa)                                      (3)
```

is anti-invariant, has `R^2=F`, and satisfies `(1)`.  Thus `(1)` and `(2)`
are equivalent, not merely one-way necessary conditions.

### Proof

Anti-invariance gives `F=R^2 in K`.  Differentiating and using `(1)` gives

```text
F'=2RR'=2kappa R/U.
```

Squaring and using `R^2=F`, `U^2=V` proves `(2)`.

Conversely, `(2)` and `(3)` give

```text
R^2=U^2(F')^2/(4kappa^2)=F.
```

Because `F` is nonconstant in characteristic zero, `F'` is nonzero.  From
`2RR'=F'` and `(3)`,

```text
R'=F'/(2R)=kappa/U.
```

This also proves that `(3)` chooses the correct sign.  QED.

The equation retains the deck squareclass exactly.  Put

```text
G=F'/(2kappa) in K*.                                  (3a)
```

Then `(2)` is equivalent to

```text
F=VG^2,                    2VG'+V'G=2kappa.           (3b)
```

Indeed the first identity is `(2)`, and differentiating it gives the second
after using `F'=2kappa G`; conversely the second makes
`(VG^2)'=2kappa G`.  In particular

```text
[F]=[V] in K*/K*^2,             K(sqrt(F))=K(U).      (3c)
```

Thus squaring `R` loses its sign but not the quadratic cover.  Formula
`(3b)` is also a linear rational differential equation for `G`; it is often
the cheapest algebraic test before constructing the passport of `F`.

## 3. Exact finite-place ramification ledger

Fix a finite point `a`.  Write

```text
m=ord_a(V)>=0,                 n=ord_a(F).
```

If `n!=0`, characteristic zero gives

```text
ord_a(F')=n-1.
```

Taking valuations in `(2)` yields

```text
m+2(n-1)=n,                 hence                 n=2-m.   (4)
```

If `n=0`, then `F` is a DVR unit and `F'` is regular.  Equation `(2)` forces

```text
m+2 ord_a(F')=0.                                  (5)
```

Consequently:

1. if `m=1`, then `F` has a simple zero at `a`;
2. `m=2` is impossible: `(4)` would demand the excluded value `n=0`, while
   the unit case `(5)` has left valuation at least two and right valuation
   zero;
3. if `m>=3`, then `F` has a pole of exact order `m-2`;
4. if `m=0`, every zero of `F` has order exactly two, `F` has no pole, and
   at every finite point where `F` is a unit its derivative is also a unit.

In rational-map language, every finite ramification point lies over zero or
infinity.  Simple roots of `V` are the simple zero fibre of `F`; a root of
`V` of multiplicity `s+2` is a pole of `F` of order `s`; and the remaining
finite zero fibre consists of double zeros.  Infinity may carry an additional
branch value, so this is not a classification of all possible three-branch-
value maps.

The immediate all-degree no-go is:

```text
any finite root of V of multiplicity exactly two
    => no nonsplit terminal Keller response.          (6)
```

### Global infinity ledger

The local data have a sharp completion.  Let:

```text
s = number of simple roots of V,
h = number of roots of V having multiplicity at least three,
e = number of double zeros of F away from V=0,
P = total finite pole degree of F,
Z = total finite zero degree of F.
```

All counts are geometric, over `C`.  If the high root multiplicities are
`m_i`, then

```text
P=sum_i(m_i-2),             Z=s+2e,
deg(V)=s+P+2h.                                      (7)
```

First suppose `n_infty=ord_infty(F)!=0`.  With `t=1/x`,

```text
ord_infty(F')=n_infty+1.
```

Taking the infinity valuation in `(2)` and using
`n_infty=P-Z` gives

```text
n_infty=deg(V)-2,             s+h+e=1.               (8)
```

Thus the genuine nonsplit, nonconstant possibilities in this chamber have
exactly one finite marked point: either `V` is linear, or `V` is a pure
odd power of one linear factor of multiplicity at least three.  (An even
pure power is split over `C`.)  The odd-power case is response-feasible; for
`V=(x-a)^m`, `F` is a scalar multiple of `(x-a)^(2-m)`.

Now suppose `ord_infty(F)=0`, so `F` tends to a nonzero finite value
`F_infty`.  Put

```text
r=ord_infty(F-F_infty)>=1.
```

Then `ord_infty(F')=r+1`, and `(2)` gives

```text
P=Z,
deg(V)=2r+2,
r=s+e+h-1.                                           (9)
```

So every non-pure-power survivor lies in the balanced chamber: `deg(V)` is
even, the finite zero and pole degrees of `F` agree, and infinity is the
only possible third branch value.  This is a ramification ledger, not an
overclassification of the resulting rational maps.

In fact the balanced chamber has an exact genus-zero Belyi passport.  Put
`N=P=Z` and `p_i=m_i-2`.  After sending the three distinguished values to
`0,infinity,1`, the three cycle partitions are

```text
over 0:          2^e 1^s,
over infinity:   (p_1,...,p_h),
over 1:          (r,1^(N-r)).                        (10)
```

Here `2^e 1^s` is the cycle type of an involution or the identity, and the
distinguished length-`r` cycle over `1` is the only possibly nontrivial
cycle.  The boundary `N=r=1` therefore has identity inertia over both `0`
and `1`; no nonexistent nontrivial cycle is being asserted.  Moreover

```text
N-r=e-h+1>=0,                  hence h<=e+1,
p_i=m_i-2<=N for every high root.                     (10a)
```

So extra high-multiplicity poles beyond the first require compensating
double zeros on the zero fibre.

They satisfy

```text
e+(N-h)+(r-1)=2N-2,
```

exactly the genus-zero Riemann--Hurwitz equation.  Equivalently the finite
combinatorial carrier is a transitive product-one permutation triple whose
zero permutation is an involution or identity and whose third permutation
has at most one nontrivial cycle.  This is a clean dessin/passport, not a tournament:
transpositions and the long cycle are typed local inertia, and orienting them
head-to-head would discard their orders.

## 4. Squarefree leading coefficient theorem

Assume now that `V` is squarefree and nonconstant.  The ledger shows that
`F` has no finite pole, hence `F` is a polynomial.  Put

```text
D=deg(V)>=1,                   n=deg(F)>=1.
```

There is no leading-term cancellation in `F'` in characteristic zero.
Comparing degrees in `(2)` gives

```text
D+2(n-1)=n,                   D=2-n.                 (11)
```

The only possibility is

```text
D=n=1.                                              (12)
```

More precisely, if `V=cx+b`, then `(2)` and the common simple zero force

```text
F=(4kappa^2/c^2)V,             R=(2kappa/c)U.        (13)
```

Therefore:

> **All-degree squarefree nonsplit corollary.**  A polynomial
> exact-square-prefix, genuine nonsplit terminal quartic Keller chart whose
> quadratic-root leading coefficient `V` is squarefree of degree at least
> two is empty, for every reduced Faber degree.  If `V` is squarefree and a
> response exists, then `V` is linear and the third response has the unique
> form `(13)`.

The linear case is a sharp boundary, not an emptiness result.  Indeed, for
`V=c(x-a)`,

```text
dx/U=(2/c)dU,
```

so `(13)` really solves the response equation.  Neither THM-2180, THM-2189,
THM-2202, THM-2214, nor THM-2217 closes this linear branch uniformly.
THM-2217 only adds `q!=0` and `R deg(q)>=1` in odd leading degree.  Existing
finite-degree closures still apply in their stated degrees.

### Geometric explanation

Equation `(1)` says on the smooth projective normalization of

```text
Y: U^2=V(x)
```

that

```text
dR=kappa dx/U.                                      (14)
```

For squarefree `deg(V)=2`, the two points at infinity have residues

```text
Res(dx/U)=-1/sqrt(lc(V)), +1/sqrt(lc(V)),
```

both nonzero, while an exact differential has zero residues.  For squarefree
degree at least three, `dx/U` is a nonzero holomorphic differential: at the
odd-degree infinity its order is `2g-2`, and at either even-degree infinity
its order is `g-1`.  A nonzero holomorphic differential cannot be the
differential of a rational function on a complete curve.  This explains
`(11)` geometrically; the polynomial degree count is the cheaper proof.

## 5. Transfer to the terminal quartic chart

In the polynomial exact-square-prefix chart, write

```text
P=H^2+L,
H=Vz^2+Bz+C_0,                 L=Az+E,
U^2=V,
```

and pass to the genuine quadratic deck:

```text
w=Uz+B/(2U),
P=w^4+2dw^2+qw+(d^2-s),        q=A/U.
```

The base polynomial mate `Q` is deck-invariant.  In the full Faber gauge,
the distinct leading fibre degrees force all odd Faber coefficients to
vanish.  For every remaining even row, the deck sends the first and third
Laurent coefficients to their negatives and fixes the second.  Hence

```text
Phi_Q is anti-invariant,
Psi_Q is invariant,
R_Q is anti-invariant.
```

The constant field is `C`, so `Phi_Q=0`, `Psi_Q in C`.  Comparing the three
coefficients in THM-2129's Hamiltonian identity gives

```text
R_Q'=kappa/U.
```

Thus Sections 2--4 apply uniformly, without fixing the reduced degree.  The
argument consumes only the exact response identity and deck character; it is
orthogonal to the degree-specific spectral-curve discriminant atlases.

## 6. Sharp hostile: squarefree part is not enough

It would be false to replace “`V` is squarefree” by “the squarefree part of
`V` has large degree.”  For every `n>=1`, put

```text
F_n=4(1-x^(-n)),
V_n=x^(n+2)(x^n-1)/n^2.
```

Direct differentiation gives

```text
V_n(F_n')^2=4F_n.
```

The polynomial `V_n` is nonsquare and its squarefree radical has degree
`n+1`, yet its repeated root at zero supplies the poles needed to make the
deck differential exact.  This family is the canonical hostile against a
genus-only argument.  The faithful invariant is the full multiplicity
divisor, equivalently the square potential `F`, not the squarefree deck
alone.

## 7. What this says about chart entry

THM-2778 closes the complete chosen-sheet split polynomial exact-prefix
family.  The smallest orthogonal obstruction inside the inherited quartic
terminal branch is now explicit:

```text
split sheet:
    base rational primitive and one-pole localization;

nonsplit sheet:
    anti-invariant primitive
      <=> base square potential V(F')^2=4kappa^2F
      <=> exact deck differential [dx/U]=0.
```

The quotient from the nonsplit deck to the base forgets precisely this
square-potential/de Rham coordinate.  A higher-degree nonsplit attack should
first classify polynomial `V=4kappa^2F/(F')^2` by the rational map `F`, then
intersect those ramification types with the two flux equations.  Repeating
the spectral discriminant atlas without this source-side sieve ignores a
degree-independent obstruction.

This does not solve the earlier chart-entry problem: an arbitrary planar
Keller pair is not known to admit this quartic source-fibre/exact-prefix
normal form.  Nor does it treat a nonpolynomial exact prefix.

## 8. Exact hostile controls

Artifacts:

```text
04-computation/jc2_nonsplit_response_square_potential_thm2784.py
05-knowledge/results/jc2_nonsplit_response_square_potential_thm2784.out
```

The assertion-free exact referee checks:

- three sharp positive potentials, including the linear boundary and the
  two-root multiplicity `(1,3)` boundary;
- all 23 polynomial potentials in a deterministic three-point rational
  monomial atlas;
- 56 exact finite-root valuation rows;
- 32 exact global infinity-ledger rows in both balanced and unbalanced
  chambers;
- the double-root prohibition and squarefree degree wall;
- six high-squarefree-part hostiles, reaching radical degree seven;
- the degree-two residues and every odd/even infinity order through genus
  eight; and
- normal/optimized byte identity against the stored transcript.

Reproduce with

```bash
python3 04-computation/jc2_nonsplit_response_square_potential_thm2784.py
python3 -O 04-computation/jc2_nonsplit_response_square_potential_thm2784.py
python3 04-computation/jc2_nonsplit_response_square_potential_hostile_audit_thm2784.py
python3 -O 04-computation/jc2_nonsplit_response_square_potential_hostile_audit_thm2784.py
```

The companion contains no Python `assert` node and reports `316` explicit
gates.  Normal, optimized, and stored outputs agree byte-for-byte.

The independent path starts from the linearized equation `(3b)`, rather
than the primary potential constructor.  It checks both descent directions,
`275` integer local cases, `1568` infinity ledgers, `96` balanced and `8`
unbalanced integer passports, `48` exact rational-map solutions/hostiles,
the squareclass identity, and the Faber parity transfer through reduced
degree `30`.  Its `3185` gates caught the identity-inertia boundary in the
passport wording and verify the repaired statement with zero `assert`
nodes; normal, optimized, and stored transcripts agree byte-for-byte.

```text
script_sha256 = 7fcbe78730b20ffe06dcc16bf640600ca5810631789e9b6d74eda113f7b053d1
output_sha256 = 07922e330700c495d0f2f541727f5a4383b17dd37c03151ea6658c111f5b5fee
independent_script_sha256 = 77f32b028ac8276c1fa292202c3e82ab02c3618a49207eb81b62f58a4362c136
independent_output_sha256 = 3e6fad51f117d612ca8de5a08d0e62f22ecc693b2eb7fd7e276ad622b1667b49
hash_basis    = LF-normalized bytes
```

Reproduction:

```bash
python3 .scratch/jc_chart_entry/nonsplit_response_square_potential_referee.py
python3 -O .scratch/jc_chart_entry/nonsplit_response_square_potential_referee.py
cmp .scratch/jc_chart_entry/nonsplit_response_square_potential_referee.out \
    <(python3 .scratch/jc_chart_entry/nonsplit_response_square_potential_referee.py)
```

Hashes over current LF bytes:

```text
script_sha256 = 05f2469928494b3d7f47c1ab5f9f40467f32c3aff6660b88ffbb71022b447f23
output_sha256 = 07922e330700c495d0f2f541727f5a4383b17dd37c03151ea6658c111f5b5fee
```

The finite atlas is a hostile referee only.  The quantifiers come from the
function-field identity, DVR proof, and polynomial degree comparison.
