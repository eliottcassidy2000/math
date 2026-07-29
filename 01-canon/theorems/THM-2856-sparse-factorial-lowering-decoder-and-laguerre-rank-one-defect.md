---
id: THM-2856
title: "Sparse factorial lowering decoder and Laguerre rank-one defect"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On every
  known k-term factorial support, the scalar
  mean together with k-1 external derivative-multiplier observations
  L(s^i H') is a rising-factorial Vandermonde decoder.  The bank is
  dimension-minimal on the mean-zero coefficient hyperplane and is
  triangularly equivalent to the ordinary multiplier bank.  The bare
  chain lowerings L(H^r H') add only the boundary value H(0) to scalar
  moments and vanish on the positive-minimum-support null branch.  In the
  optimal Laguerre quotient, canonical-representative differentiation
  obeys an exact rank-one corrected Weyl relation; every intrinsic
  derivation vanishes, and the correction is exactly the highest
  inverse-Hankel coefficient selector.  These are external observability
  statements, not ordinary SFC(4), HYP-8765, or a new Gaussian-moment
  proof.
source: root/gmc-lowering-observable-2026-07-28
depends_on:
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
related:
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2842-ordered-positive-cone-vandermonde-multiplier-observability
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
script: 04-computation/gmc_sparse_lowering_decoder_thm2856.py
output: 05-knowledge/results/gmc_sparse_lowering_decoder_thm2856.out
script_sha256: 359c467af728973af83b67eafc9f2bfa699b120a4a400e7297d06fb1675ceaa2
output_sha256: 7e19e90f3c27c5b4dec0b81d1d752b38f97da51adc323194dc33bfb299e14c98
hash_basis: LF-normalized bytes
---

# THM-2856 -- sparse factorial lowering decoder and Laguerre rank-one defect

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `K` be a characteristic-zero field and let

```text
L:K[s] -> K,                         L(s^n)=n!.          (1)
```

THM-2815 gives the optimal finite Laguerre carrier and exact
inverse-Hankel coefficient selectors.  THM-2842 identifies multiplier
access with radial-variance jets and proves an ordered-positive-cone
Vandermonde theorem.  THM-2848 exposes the integration-by-parts signal
`L(Z^rZ')` but leaves its access open.

The present theorem separates three questions which those results place
next to one another.

1. On a **known sparse monomial support**, a finite external lowering bank
   recovers every coefficient, uniformly in the support height.
2. The unmultiplied chain lowerings from THM-2848 do not supply that bank:
   modulo scalar moments they contain only the endpoint `H(0)`.
3. The optimal Laguerre quotient admits a canonical representative
   lowering, but not an intrinsic derivation.  Its failure is one exact
   rank-one boundary selector.

## 1. The sparse rising-factorial decoder

Fix distinct nonnegative integers

```text
0<=a_0<a_1<...<a_(k-1),              f_a=s^a/a!,       (2)
H=sum_(j=0)^(k-1)c_j f_(a_j).                           (3)
```

Write

```text
x^(overline 0)=1,
x^(overline i)=x(x+1)...(x+i-1),       i>=1.           (4)
```

Define the external observations

```text
O_0(H)=L(H),
O_i(H)=L(s^i H'),                       1<=i<k.         (5)
```

Since

```text
f_a'=s^(a-1)/(a-1)!                    for a>=1
```

and `f_0'=0`, one has, uniformly also at `a=0`,

```text
O_i(f_a)=a^(overline i),                1<=i<k,
O_0(f_a)=1.                                                (6)
```

Therefore

```text
O_i(H)=sum_j c_j a_j^(overline i),      0<=i<k.         (7)
```

The polynomials

```text
1, x, x(x+1), ..., x^(overline (k-1))
```

are monic of respective degrees `0,1,...,k-1`.  Changing from them to
the monomial basis is unit triangular, so

```text
det_(0<=i,j<k)[a_j^(overline i)]
 =product_(0<=r<t<k)(a_t-a_r)
 !=0.                                                     (8)
```

Thus `(5)` is an exact coefficient decoder on every known `k`-term
factorial support.  No positivity, adjacency, ordering of cone blocks,
or height cutoff is used; only distinctness of the exponent labels is
needed.

## 2. Sharp four-slot observability and dimension

The mean-zero coefficient space

```text
W={H as in (3):L(H)=sum_j c_j=0}                        (9)
```

has dimension `k-1`.  Equation `(8)` says that

```text
H |-> (O_1(H),...,O_(k-1)(H))                          (10)
```

is injective on `W`.  Conversely, no bank of fewer than `k-1` scalar
linear observations can be injective on a `(k-1)`-dimensional vector
space.  Hence `(10)` is dimension-minimal among linear external
observation banks.

For four slots this gives the exact uniform statement

```text
L(H)=L(sH')=L(s^2H')=L(s^3H')=0
                    implies H=0                         (11)
```

for every known support `a_0<a_1<a_2<a_3`, at arbitrary height.  Once
`L(H)=0` is already known, precisely three further linear scalar
observations suffice, and three are necessary for uniform linear
coefficient recovery.

This is finite **external observability**, not ordinary SFC(4).  The
ordinary moment hypothesis supplies `L(H^m)`, not the three new values in
`(11)`.  Also, `(8)` assumes the exponent support is known; it does not
solve a simultaneous sparse-support identification problem.

## 3. Lowering and multiplier banks are the same sidecar

Put

```text
M_i(H)=L(s^iH),                         i>=0.            (12)
```

For `i>=1`, integration by parts has zero boundary term because of the
factor `s^i`, and gives

```text
O_i(H)=M_i(H)-i M_(i-1)(H),                              (13)
M_i(H)=i M_(i-1)(H)+O_i(H).                              (14)
```

Together with `O_0=M_0`, equations `(13)--(14)` are mutually inverse
unit-triangular changes of observation coordinates.  Thus the sparse
lowering bank is exactly as informative as

```text
L(H),L(sH),...,L(s^(k-1)H),                            (15)
```

but it is not extra information for free.  It is the derivative form of
the multiplier sidecar already isolated by THM-2815.

This result is distinct from THM-2842.  That theorem proves strict
Vandermonde orientation on ordered positive adjacent-difference cones
and identifies its multipliers with variance jets.  Equation `(8)`
instead treats arbitrary known sparse monomial supports, with no cone
sign or transport-order hypothesis, and uses derivative-multiplier
observations.

## 4. The bare chain lowerings are endpoint-redundant

The THM-2848 integration-by-parts identity is complete.  For every
polynomial `H` and every `r>=0`,

```text
B_r(H):=L(H^rH')
 =[L(H^(r+1))-H(0)^(r+1)]/(r+1).                       (16)
```

Consequently the bank

```text
B_0(H),...,B_(m-1)(H)                                 (17)
```

contains no information beyond the first `m` scalar moments and the one
boundary coordinate `H(0)`.  At a first-window four-moment null point,

```text
B_r(H)=-H(0)^(r+1)/(r+1),                0<=r<=3.      (18)
```

If the minimum support exponent is positive, then `H(0)=0`, and all four
bare chain lowerings vanish.  This boundary is realized before the
quartic exit: THM-2846 gives a nonzero factorial four-slot polynomial on
support `{1,2,3,4}` with

```text
L(H)=L(H^2)=L(H^3)=0,              H(0)=0,
```

so `(16)` makes `B_0=B_1=B_2=0` while `H!=0`.  Its fourth moment is
nonzero, as it must be; this is a hostile to promoting the bare chain
bank into a first-three-moment detector, not a counterexample to SFC(4).

The distinction is sharp:

```text
bare chain access       L(H^rH')       -> endpoint only modulo moments;
multiplied lowering     L(s^iH')       -> full sparse coefficient decoder.
```

## 5. The finite quotient pays one rank-one Weyl defect

Let

```text
q(s)=s^D+q_(D-1)s^(D-1)+...+q_0
```

be monic over `K`, and set

```text
A=K[s]/(q),                    u=[s].                  (19)
```

Every class has a unique representative of degree below `D`.  Relative
to the pointed generator `u`, define

```text
partial_q(sum_(j=0)^(D-1)b_j u^j)
       =sum_(j=1)^(D-1)j b_j u^(j-1),                 (20)
tau(sum_(j=0)^(D-1)b_j u^j)=b_(D-1),                  (21)
M_u(f)=uf.                                             (22)
```

Here `partial_q` differentiates the canonical representative; it is
not claimed to be a derivation of the quotient algebra.

On `1,u,...,u^(D-2)`,

```text
(partial_q M_u-M_u partial_q)(u^j)=u^j.                (23)
```

On the last basis vector, use

```text
u^D=-sum_(j=0)^(D-1)q_j u^j
```

to obtain

```text
(partial_q M_u-M_u partial_q)(u^(D-1))
 =u^(D-1)-q'(u).                                      (24)
```

Equivalently, as an exact operator identity,

```text
[partial_q,M_u]
 =I-q'(u) tensor tau,                                 (25)
```

where `(q'(u) tensor tau)(f)=tau(f)q'(u)`.  In
characteristic zero `q'(u)!=0`, so the correction has rank exactly one.
Its trace is `D`, because the top coefficient of `q'` is `D`; hence the
two terms on the right of `(25)` have equal trace.  This is the minimal
finite-dimensional correction to the Weyl relation
`[partial_s,M_s]=I`.

The sign and side of the correction are load-bearing.  Equation `(25)`
uses

```text
[partial_q,M_u]=partial_q M_u-M_u partial_q,
```

not the opposite commutator.

## 6. The optimal Laguerre quotient has no intrinsic lowering

Now specialize to the optimal carrier

```text
A_D=K[s]/(ell_D).                                     (26)
```

THM-2815 proves that `ell_D` has `D` simple positive roots.  Hence

```text
gcd(ell_D,ell_D')=1
```

and `ell_D'(u)` is a unit in `A_D`.  Every `K`-derivation

```text
delta:A_D->A_D
```

therefore vanishes.  Indeed,

```text
0=delta(ell_D(u))=ell_D'(u)delta(u)
```

forces `delta(u)=0`, and `A_D=K[u]` then forces `delta=0`.  In
particular no intrinsic quotient derivation can satisfy `delta(u)=1`.

The same obstruction is visible without the derivation axiom.  The zero
class has representatives `0` and `ell_D`, but ordinary differentiation
sends them to

```text
0,                         ell_D',
```

whose quotient classes differ by a unit.  Thus `partial_D` in `(20)` is
well-defined only because the unique degree-below-`D` representative was
chosen.  It is pointed-representative data, not an operation descending
from `K[s]`.

The rank-one term in `(25)` is exactly THM-2815's highest coefficient
selector.  Let

```text
phi_(D-1)^(D-1)
```

be its inverse-Hankel selector.  For the canonical representative
`f` of degree below `D`,

```text
tau(f)=[s^(D-1)]f
      =L(f phi_(D-1)^(D-1))
      =Lambda_D([f][phi_(D-1)^(D-1)]).                (27)
```

The product has degree at most `2D-2`, inside the exact Laguerre horizon.
So `(25)` does not remove the access problem: it localizes it to the
single hardest coefficient functional, which still requires an inserted
multiplier observation.

Conversely, the quotient augmented by canonical-representative lowering
does realize the sparse decoder.  If `a_(k-1)=a_max` and
`D=a_max+1`, then `H` is its own canonical representative and, for
`1<=i<k`,

```text
O_i(H)=Lambda_D(u^i partial_D[H]).                     (28)
```

Indeed

```text
deg(s^iH')<=a_max+k-2<=2D-3,
```

because `k` distinct nonnegative exponents imply `k<=a_max+1`.
Thus `(28)` remains inside the exact quadrature horizon.  The added
operator gives finite sparse observability; the multiplicative carrier
and scalar readout alone do not.

## 7. Connection contract and scope

| item | exact content |
|---|---|
| source | a polynomial on one known `k`-term factorial support |
| external map | `H -> (L(H),L(sH'),...,L(s^(k-1)H'))` |
| preserved predicate | every coefficient, by the rising-factorial determinant `(8)` |
| equivalent sidecar | the multiplier bank `L(s^iH)` through `(13)--(14)` |
| bare-chain loss | `L(H^rH')` retains only `H(0)` modulo scalar moments |
| finite carrier | `A_D=K[s]/(ell_D)` with its scalar readout |
| missing coordinate | the representative lowering / top selector `tau` |
| exact obstruction | `[partial_D,M_u]=I-ell_D'(u) tensor tau`; every intrinsic derivation is zero |
| sharp hostiles | THM-2846 for bare-chain blindness; the equivalent representatives `0` and `ell_D` have inequivalent derivatives |

This theorem does not prove ordinary SFC(4), a shifted-window theorem,
HYP-8765, arbitrary Wick-channel separation, or a new Gaussian Moment
result.  It proves exactly what finite lowering access would buy and why
the optimal multiplicative carrier cannot manufacture that access
internally.

## 8. Exact companion

Run

```text
python 04-computation/gmc_sparse_lowering_decoder_thm2856.py
python -O 04-computation/gmc_sparse_lowering_decoder_thm2856.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_sparse_lowering_decoder_thm2856.out.
```

The standard-library companion uses only integers and
`fractions.Fraction`.  It verifies:

1. `(8)` on all `4,095` supports of sizes `1,...,6` in
   `{0,...,12}`;
2. rank three of the mean-zero four-slot lowering bank on all
   `binom(16,4)=1,820` supports in `{0,...,15}`;
3. direct factorial integration, the triangular identities, and exact
   coefficient recovery on the hostile-height support `{0,3,8,15}`;
4. `(16)` on `18` exact polynomial/power controls and the nonzero
   `6s-s^3` mean-zero, boundary-zero, first-lowering-blind hostile;
5. the signed commutator `(25)`, rank-one defect, trace cancellation,
   squarefreeness, and representative hostile for `D=1,...,14`; and
6. the highest inverse-Hankel selector on `105` monomial cases.

Every truth-bearing gate is an explicit exception.  There is no Python
`assert`, floating-point decision, optional CAS, or scratch dependency.
The universal statements are proved above rather than inferred from the
finite controls.

An independent hostile audit rederived the rising-factorial
Vandermonde determinant and its sign, the dimension-minimality
statement, the triangular multiplier/lowering equivalence, the endpoint
redundancy of the bare chain, and the sign and rank of the corrected
commutator.  It also checked the Laguerre squarefree/unit argument,
vanishing of intrinsic derivations, the highest-coefficient selector,
all degree horizons, and the stated non-implications.  Independent
normal and optimized runs byte-match the stored transcript and both
declared LF hashes.

**QED.**
