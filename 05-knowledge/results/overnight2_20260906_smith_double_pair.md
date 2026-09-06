# The ternary four-node double-pair Smith law

**Status: PROVED by a finite symbolic minor certificate + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.** The certificate handles all
positive integer scales and depths, with arbitrary unit lifts. The bounded
integer checks are additional controls, not the source of the universal
quantifiers. No theorem ID, external novelty, or priority claim is made.

Root accepted the weighted-minor proof and attaining witnesses. The
[independent audit](overnight2_20260906_smith_double_pair_audit.md)
reconstructed all 923 minors using a third polynomial-domain engine,
checked all 14,501 coefficient bounds and all thirteen attaining
factorizations, and audited the sorting, boundary, symmetry and precision
claims. The universal proof passed without a mathematical correction.

## Inheritance and the precise missing case

The closest proved mechanism is the weighted derivative-bank competition in
[the first odd-prime mixed cluster](overnight_20260906_smith_mixed_cluster.md),
which covers residue occupancy `(2,1,1)` at prime three. The present occupancy
is `(2,2)`, and is outside that theorem. Its proof retains the complete
observer lattice when two depths compete.

The earlier boundary is
[THM-4419 / twojet-prime-wall-precision-and-dyadic-triple-smith-law](../../01-canon/theorems/THM-4419-twojet-prime-wall-precision-and-dyadic-triple-smith-law.md).
The canonical hostile is its failure of an unchanged Smith prefix after
adjoining a node. The corrected near miss is
[THM-4010 / confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall](../../01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md):
the full determinant does not determine the intermediate invariants. The
least-used sidecar here is the weighted support of every intermediate minor.

The incoming
[arbitrary-three-node note](three-node-smith-overnight-hexagon-sep05.md) is
was PROVISIONAL at inheritance; incoming commit `4d1ad2a39` subsequently
promoted it as
[THM-4429 / arbitrary-three-node-two-jet-Smith-form](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md).
It suggested asking whether the distance tree still suffices for four
nodes; **no part of the present proof depends on that three-node result**.

The live concept board is: the outer scale; the two inner depths; derivative
saturation at prime three; weighted minor supports; and metric information
versus determinant mass. The initial unit-lift test was positive, but a
straight continuation of the shallow-depth spectrum failed. That failure
led to the three exact switching boundaries below.

## 1. Universal formula

Let `e,d,r>=1`, and let `a,b` be arbitrary integers prime to three. Set

```text
A=3^d*a,       B=3^r*b,
x=3^e*(0,1,A,1+B),
s=min(d,r),    t=max(d,r).
```

The four nodes are distinct: within each pair the difference is nonzero,
and between pairs it is a unit after dividing by `3^e`. On the fixed module
of polynomials of degree below eight, observe the value and first Hasse
derivative at every node. This is a square integral eight-by-eight map.

Its two smallest three-adic Smith exponents are zero. Write `D_k` for the
valuation of the gcd of the `k`-minors of the six-by-six matrix remaining
after clearing the first value and derivative, with `D_0=0`. Then

| k | D_k |
|---|---|
| 1 | `e` |
| 2 | `min(4e,3e+s+1)` |
| 3 | `min(7e+s,6e+s+t+1)` |
| 4 | `min(11e+s+t,12e+4s)` |
| 5 | `17e+4s+t` |
| 6 | `24e+4s+4t` |

Thus the full spectrum is `(0,0,D_1,D_2-D_1,...,D_6-D_5)`.
An explicit form is obtained by putting

```text
lambda=min(e,s+1),
mu=min(e,t+1),
nu=min(t,e+3s):

(0,0,
 e,
 2e+lambda,
 3e+s+mu-lambda,
 5e+nu-mu,
 6e+3s+t-nu,
 7e+3t).                                             (1)
```

The list is already increasing after the initial equal zeros. For example,
`lambda<=mu<=e` and `nu>=s`. The gap between the sixth and fifth entries is
at least `nu-s+lambda>0`. The next gap is
`e+3s+t-2nu+mu>=mu>0`, by either branch of `nu`; the final gap is
`e+2t-3s+nu>=e>0`. The earlier gaps are immediate.

In particular **the complete pairwise three-adic distance tree suffices**
for this cluster. No further unit digit is needed. The four cross-pair
differences have valuation `e`, and the two within-pair differences have
valuations `e+s,e+t`. Conversely every four-node configuration with that
distance pattern reduces to this form over `Z_(3)` by translation and a
unit change of coordinate. The polynomial proof works for `a,b` in
`Z_(3)^*` as well: only their unit valuations are used.

The largest exponent is `7e+3t`. For every `N>=1`, observations modulo
`3^(N+7e+3t)` recover all eight source coefficients modulo `3^N`, and this
uniform precision loss is sharp, by passing to Smith coordinates. This is
a statement about the declared degree box and full two-jet map.

## 2. The residual matrix and exact weight

The columns of degrees zero and one, followed by integral row elimination,
clear the value and derivative at zero. The remaining columns are degrees
`q=2,...,7`, and the residual rows are

```text
0: V(1),  1: D(1),  2: V(A),  3: D(A),
4: V(1+B), 5: D(1+B).
```

Here `V(u)_q=u^q` and `D(u)_q=q*u^(q-1)`. If selected columns have degrees
`Q` and selected rows contain `j` derivative rows, their actual minor is

```text
3^(e*w) * P(A,B),       w=sum(Q)-j,                 (2)
```

where `P` is the corresponding residual polynomial determinant. This is
an exact factorization, not a scaling equivalence that discards integral
information. Every residual entry is divisible by three, so the two cleared
unit factors are indeed the first two Smith factors.

For a nonzero term `c*A^i*B^j` of `P`, its valuation after (2) is

```text
F(e,d,r)=w*e+i*d+j*r+v_3(c).                         (3)
```

Cancellation can only increase a minor's valuation. A matching upper
bound requires an actual minor whose least-valuation term cannot cancel.
Both sides are certified below.

## 3. A finite symbolic lower-bound certificate

There are exactly `sum_(k=1)^6 binom(6,k)^2=923` nonempty minors of the
residual matrix. Their complete expansions contain 14,501 nonzero
coefficients. The following table gives a sufficient list of affine forms,
written as tuples `(w,i,j,c)` representing `w*e+i*d+j*r+c`:

| Rank | Forms |
|---|---|
| 1 | `(1,0,0,0)` |
| 2 | `(4,0,0,0)`, `(3,1,0,1)`, `(3,0,1,1)` |
| 3 | `(7,1,0,0)`, `(7,0,1,0)`, `(6,1,1,1)` |
| 4 | `(11,1,1,0)`, `(12,4,0,0)`, `(12,0,4,0)` |
| 5 | `(17,4,1,0)`, `(17,1,4,0)` |
| 6 | `(24,4,4,0)` |

**Finite polynomial lemma.** Every nonzero term of every rank-`k` minor is
bounded below, for all `e,d,r>=1`, by at least one form in its table row.

Here is a fully specified finite certificate for the lemma. Replace a form
`f=(w,i,j,c)` by

```text
T(f)=(w,i,j,w+i+j+c).
```

For every coefficient term of every minor, the accompanying script checks
that some listed form `g` satisfies `T(g)<=T(f)` coordinatewise. This implies
`g(e,d,r)<=f(e,d,r)` because

```text
f-g = (w_f-w_g)(e-1)+(i_f-i_g)(d-1)+(j_f-j_g)(r-1)
      +(w_f+i_f+j_f+c_f)-(w_g+i_g+j_g+c_g) >= 0.     (4)
```

The check therefore covers unbounded integer depths; it is not a sample of
their values. Summing polynomial terms gives the desired lower bound for
each minor, and taking the gcd gives the same bound for the determinantal
divisor.

For reproducibility without trusting one determinant implementation, all
923 coefficient dictionaries are constructed twice:

1. Laplace expansion using cached minors and exact integer polynomial
   addition/multiplication;
2. an independent literal permutation expansion. For a column permutation
   `q_l` assigned to selected rows, its term is
   `sign * product_(derivative rows) q_l * A^u * (1+B)^v`, where `u` and `v`
   are the sums of `q_l-order_l` over the `A` and `1+B` rows respectively.
   The binomial theorem expands this term directly.

Every coefficient dictionary must agree exactly, and every coefficient
must satisfy (4). These finite algebraic gates establish the polynomial
lemma and carry its universal implication. No floating-point optimization,
feasibility solver, rational-function denominator clearing, or guessed
convex hull is used.

## 4. Thirteen explicit attaining minors

The lower bounds are exact because each listed form comes from the minor
in this table. Row labels refer to Section 2; columns give actual degrees.
The final column is its residual determinant `P`, with an irrelevant sign
permitted. All displayed factors other than the powers of `A,B` and the
explicit constants are units modulo three when `A,B` are divisible by
three.

| Rank | Rows | Degrees | w | P(A,B) |
|---|---|---|---|---|
| 1 | `1` | `2` | 1 | `2` |
| 2 | `0,1` | `2,3` | 4 | `1` |
| 2 | `1,3` | `2,3` | 3 | `6A(A-1)` |
| 2 | `1,5` | `2,3` | 3 | `6B(B+1)` |
| 3 | `0,1,3` | `2,3,4` | 7 | `2A(A-1)(2A-1)` |
| 3 | `0,1,5` | `2,3,4` | 7 | `2B(B+1)(2B+1)` |
| 3 | `1,3,5` | `2,3,4` | 6 | `24AB(A-1)(B+1)(A-B-1)` |
| 4 | `0,1,3,5` | `2,3,4,5` | 11 | `2AB(A-1)(B+1)(A-B-1)(10AB+5A-5B-2)` |
| 4 | `0,1,2,3` | `2,3,4,5` | 12 | `A^4(A-1)^4` |
| 4 | `0,1,4,5` | `2,3,4,5` | 12 | `B^4(B+1)^4` |
| 5 | `0,1,2,3,5` | `2,3,4,5,6` | 17 | `2A^4 B(A-1)^4(B+1)(A-B-1)(2AB+A-3B^2-4B-1)` |
| 5 | `0,1,3,4,5` | `2,3,4,5,6` | 17 | `2AB^4(A-1)(B+1)^4(A-B-1)(3A^2-2AB-4A+B+1)` |
| 6 | `0,1,2,3,4,5` | `2,3,4,5,6,7` | 24 | `A^4 B^4(A-1)^4(B+1)^4(A-B-1)^4` |

Thus the valuations of these actual minors are precisely the forms in
Section 3. The script also verifies a stronger direct condition: the
selected term has strictly smaller valuation than every other term for
every positive depth, so its unit coefficient cannot cancel for any lifts.
Taking the minimum in each row gives the six `D_k` formulas and proves
(1). The last row is also the full confluent Vandermonde determinant.

The three switches have concrete meanings. At rank two, the same-node
value/derivative witness competes with a cheaper-scale derivative pair
carrying the factor `6A` or `6B`. At rank three, the three-derivative bank
has the factor `24AB`, and competes with retaining a value row. At rank
four, retaining both shallow-pair values costs `4s` in collision depth,
while retaining three derivative rows costs `s+t`. This last competition
produces the new boundary `t=e+3s`.

## 5. Hostiles, consequences, and the next missing coordinate

For `e=1,s=1,t=5`, the exact spectrum is

```text
(0,0,1,3,4,8,10,22).
```

The straight shallow-depth continuation, even after sorting, would give
`(0,0,1,3,4,9,9,22)`. Its first failed implication is that the middle factors
can be continued independently after one pair becomes much deeper. The
repair is the rank-four comparison `min(11e+s+t,12e+4s)`.

At `e=4,d=r=1`, the exact list is `(0,0,4,10,13,19,27,31)`; the fourth
exponent is already below `3e`. This is a separate scale/saturation switch,
controlled by `s+1=e`. Both examples retain the full observer.

The determinant-only hostile is equally small: at `e=1`, depth splits
`(1,3)` and `(2,2)` both have total determinant valuation 40, but spectra

```text
(0,0,1,3,4,7,9,16),
(0,0,1,3,5,6,12,13).
```

So total collision depth is insufficient, whereas the full two-depth tree
is sufficient here. This is the strongest surviving metric statement, not
a general theorem for all four-node configurations or all primes.

At the separate boundary `e=0`, the two pairs lie in distinct residue
classes. Integral CRT gives the sorted union
`(0,0,d,3d)` and `(0,0,r,3r)`, and substitution into (1) agrees. This is
an independent boundary control; the lower-bound proof above assumes
`e>=1`. CRT cannot be used to glue pair Smith lists across a nonunit outer
scale.

The source-to-target map is now explicit: the complete observer maps to
all weighted coefficient supports of its minors, then to their lower
envelope. The first step preserves every determinantal ideal; the envelope
keeps the ideal valuations. Unit digits are discarded only after actual
unit-leading witness minors prove that cancellation cannot change those
valuations. The needed sidecar is the witness for each competing minimum.
The cheap decisive test was the same-distance/different-unit comparison;
the positive result was followed to a universal certificate.

The next four-node ternary shape is occupancy `(3,1)`, where the inner
triple has its own two-scale tree. A second divided derivative row may meet
the characteristic-three loss in the third ordinary derivative. Whether
its entire metric tree still suffices is OPEN; this report does not infer
it from the `(2,2)` result.

## 6. Exact reproduction and evidence manifest

Run from the repository root:

```text
python3 -B 04-computation/overnight2_20260906_smith_double_pair.py
python3 -B -O 04-computation/overnight2_20260906_smith_double_pair.py
```

The [self-contained source](../../04-computation/overnight2_20260906_smith_double_pair.py)
and [preserved output](overnight2_20260906_smith_double_pair.out) record:

- all 923 symbolic minors, 14,501 nonzero coefficient checks, and 13
  universally attaining witness forms;
- 6,400 complete modular Smith lists: `e,d,r=1,...,4` and every pair from
  the ten units `1,2,4,5,7,8,10,11,-1,-2`;
- 400 independent rational DVR lists: `e=1,...,4`, `d,r=1,...,10`, with
  signed lifts, translation, and a unit change of scale;
- every minor of the original eight-by-eight matrix in four exact controls,
  including the unequal-depth and nonunit-scale hostiles;
- 36 independent unit-separated CRT controls at `e=0`.

The modular precision is one more than the known full determinant valuation,
so no Smith exponent can be hidden by truncation. Normal and optimized
Python executions produce identical output, with 1,359,304 explicit gates.
There are no optimization-disabled assertions or floating-point constants.

```text
source SHA256:
bd2b5d2a9e080dd7a5422892ef5a22a4798282c7ad593376712d03e1aa186a6d
output SHA256:
78fbfe52e3c8e7c51eda03b8b449beaffbe23469848ecc34d2b6c9a478a16fdd
complete symbolic-minor semantic SHA256:
b6800af922d7e9bc60b05550de834c1c0fc06940c67b09bc1e347db82a1beff7
unit-lift semantic SHA256:
56d974addfdd5c7439e1c314925d8470d54ed2737e0d8107156d5c2677df212a
```
