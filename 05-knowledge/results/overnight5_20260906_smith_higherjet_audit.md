# Independent audit of arbitrary-jet precision and the dyadic unit boundary

**PASS: written proof and bounded independent exact controls.** This is an
additional audit of incoming **PROVED** THM-4443, not a new discovery claim.
Read-only source snapshot: `058a8ded98cfa25fd90f061efe5d54903c9b7379`:

* `01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md`;
* `05-knowledge/results/hermite-higher-jet-unit-boundary-overnight-hexagon-sep05.md`;
* its inherited THM-4439 terminal-cluster theorem and complete proof.

The observer comprises Hasse jets of orders `0,...,m_i-1` at every node,
on degree less than `sum m_i`. In the uniform three-jet boundary, the newly
added observation is order two. No ordinary-derivative normalization is
substituted. The all-node arbitrary-multiplicity denominator formula and
the full dyadic three-node uniform-three-jet law both pass this audit.

## Exact inverse and sharpness

Write `Q_i(X)=product_(j!=i)(X-x_j)^(m_j)` and
`a_(i,l)=[T^l]/Q_i(x_i+T)` for `0<=l<m_i`. For datum `(i,r)`, the cardinal
polynomial

```
(X-x_i)^r Q_i(X) sum_(l=0)^(m_i-r-1) a_(i,l)(X-x_i)^l
```

has the required local jet and vanishes to the prescribed orders at every
other node. Its degree is at most `sum m_i-1`. Its top coefficient is
exactly `a_(i,m_i-r-1)` because Q_i is monic and all earlier summands have
lower degree. Thus every proposed denominator is attained as a literal
inverse entry. This proves both inequalities in the largest-integer-factor
claim, including arbitrary unequal multiplicities and zero Taylor entries.
The order-zero coefficient prevents a negative largest prime exponent.
Smith changes then give the sharp precision statement for every output
modulus; a coordinate with largest exponent supplies the one-less-digit
hostile. A single node has an integral translated monomial basis and loss
zero, as stated.

Under common dilation by p^e, a coefficient `(i,l)` has denominator cost
increased by `(sum_(j!=i)m_j+l)e`. Therefore taking the maximum first at one
shallow scale would discard information. The incoming proof correctly
retains the upper envelope over all orders. For uniform multiplicity k,
the repeated-reciprocal complete homogeneous coefficient formula and the
two-node binomial special case follow immediately from the geometric
series expansion; their multiplicity factors are correct.

## Complete dyadic law

At a three-jet node with differences u,v, direct reciprocal expansion gives

```
a0=1/(uv)^3,
a1=-3(u+v)/(uv)^4,
a2=3(2u^2+3uv+2v^2)/(uv)^5.
```

For metric depths `(e,e,e+d)`, d>=2, the closest-pair quadratic numerator
has term valuations `2e+1, 2e+d, 2e+2d+1`; the first is uniquely least.
Its resulting loss `8e+5d-1` dominates the other pair orders and all
outsider orders. No unit-dependent cancellation is possible in this range.

When d=1, choose closest endpoints x,y and outsider z and set
`t=2(z-x)/(y-x)`, an odd two-adic unit. Removing the common numerator
power at the two endpoints leaves

```
P(t)=t^2+3t+4,   P(2-t)=t^2-7t+14.
```

The proof of their *minimum* valuation is exact for arbitrary two-adic
units. For t=3 mod4 it is one. Otherwise write t=1+4a; the polynomials
are `4(2+5a+4a^2)` and `4(2-5a+4a^2)`. Odd a gives minimum two. For a=2b
they become `8(1+5b+8b^2)` and `8(1-5b+8b^2)`; even b gives three. For
odd b both parentheses are even, and their difference has valuation one,
so their minimum valuation is exactly one and the original minimum is
four. This establishes all four residue classes without assuming that an
individual endpoint coefficient has bounded valuation.

The resulting `Gamma=3,2,1,0` corresponds respectively to
`t=3 mod4, 5 mod8, 1 mod16, 9 mod16`. The outsider order-two candidate is
exactly 8e, so it does not exceed the pair maximum since Gamma>=0. The
pair order-one candidate is 7e+4. Therefore

```
L2=max(7e+4,8e+Gamma(t)).
```

Exchanging endpoints replaces t by 2-t and preserves each class; affine
unit changes preserve t. The shallow maximum can mask the unit branch.
The proof appropriately claims a sufficient intrinsic sidecar, not its
necessity at every scale. All four branches become distinguishable at
e>=4 directly from the displayed maximum; this is an algebraic consequence,
not an additional census result.

## What fails in the terminal-metric principle

THM-4439 applies only to the complete *uniform two-jet* observer. Its
metric-only worst inverse entry arises from order-zero denominators and
the first reciprocal sum `2q_x`. At a complete terminal cluster, a
low-degree polynomial prevents simultaneous cancellation beyond one
digit. For dyadic terminal pairs the internal reciprocal sum is itself
a unit; the factor two causes the digit loss.

The third jet introduces a quadratic reciprocal coefficient. Its two
endpoint numerators can simultaneously lose one, two, three, or four
digits depending on the unit residue. Maximizing over the complete
terminal pair therefore does not remove the unit coordinate. This is the
first failed implication. It does **not** refute the two-jet theorem, and
the three-node law does not claim that a nonterminal node must become the
maximizer: here the terminal pair still attains the maximum. The failure
is the universal terminal cancellation budget after changing the observer.

For the isometric uniform triples `2^e(0,1,2)` and `2^e(1,0,3)`, both
largest losses agree at e=0,1 and differ for every e>=2. At e=2 the
independently reconstructed full dyadic partitions are

```
(0,0,0,2,6,8,11,18,18),
(0,0,0,2,6,8,11,17,19).
```

Nonuniformity already breaks the uniform two-jet conclusion without an
order-two observation. The weighted sets `2^e(0,2,1)` and `2^e(1,3,0)`
with multiplicities `(2,2,1)` have the same labelled metric, but their
weighted reciprocal sums have different simultaneous cancellations. At
e=2 the full partitions are `(0,0,4,7,9)` and `(0,0,4,8,8)`. The local
inverse formula gives the stated all-scale losses
`max(3e+2,4e+1)` and `max(3e+2,4e)`. The remaining order-zero candidate
from the outsider is retained, so the proof does not accidentally discard
the branch that can dominate at greater depth.

Minimality is correctly scoped: k=1 and k=2 are metric-only at all node
counts, and one/two-node uniform observations are metric-only for every k.
Thus `(uniform k,n)=(3,3)` is the first possible multiplicity/node-count
boundary. No smallest-diameter or general higher-node classification is
proved or inferred.

## Connection to earlier Smith blindness

The fourth-round terminal audit retained the ternary isometric four-node
pair `(0,9,27,81)` and `(0,9,54,81)`: the largest exponent is 22 for both,
while intermediate factors differ. That was a loss of *full partition*
information under the metric quotient. THM-4443 changes the observable:
the extra Taylor order or the unequal multiplicity makes units visible
even to the largest inverse denominator. These are compatible failures at
different targets, rather than contradictory metric principles.

Source-target-loss-next-test contract: start with multiplicity-labelled
integer nodes; map to their p-adic distance tree; target sharp coefficient
recovery precision. This quotient preserves the complete uniform two-jet
target by THM-4439 but destroys higher reciprocal Taylor coefficients.
Retain observer multiplicities and the coefficient valuations with their
dilation slopes; in the first dyadic uniform case, Gamma is an explicit
finite replacement. The cheapest next test in any larger observer class
is to compare the entire dilation envelopes of isometric sets, not only
their largest exponent at scale zero. This audit does not pursue that
larger classification.

## Independent exact controls and freeze

The companion program imports no incoming producer, referee, or repository
mathematical routine. It constructs literal Hasse matrices, inverts them
by rational Gaussian elimination, reconstructs reciprocal coefficients
using polynomial convolution and formal inversion, and compares the
*largest integer* Smith factor using SymPy's integer Smith routine. It
checks every inverse product, every top-row denominator witness, and the
full Hasse determinant `product_(i<j)|x_i-x_j|^(m_i*m_j)`.

There are **77 exact matrix cases**: all eight odd residue classes at
e=0,2,5; 18 deeper-pair cases; four signed unit-affine controls; both
uniform and heterogeneous hostile families at four scales; singleton,
two-node, signed heterogeneous, and complete two-jet controls. The old
two-jet metric twin control remains equal, preventing the audit from
accidentally switching observers. The finite matrix bank supports the
written proof; it is not an extrapolation premise for the unbounded law.

```
python -B 04-computation/overnight5_20260906_smith_higherjet_audit.py
python -B -O 04-computation/overnight5_20260906_smith_higherjet_audit.py
```

Normal and optimized LF outputs agree: **6,806 optimization-live gates PASS**.
Frozen hashes:

```
source 0a63553b63fa8a9152dc2d6a5c59323787723aa1ce5ae4f961b1f900e1cf1dbb
output 2f7f5c42ef443f44fc504e1d0713140c7e528901b29a6865632cc4d7095453e4
incoming canon 0e29e23680794e5aef6a6c557ff919e0291ef909a899af1c68acc6b617576346
incoming proof c56c9546a62d907ed4ffbf1dda91ac80548321d40119daee8259020428b64369
```

No incoming files, repository navigation, theorem namespaces, or Git state
were changed by this audit. It uses primary repository proofs and makes
no external literature or discovery-priority claim.

**Filing note.** Root filed this completed audit after checkpoint `07b2d91b2`.
Local paths were made repository-relative and output line endings normalized;
normal and optimized verification was rerun with matching output. The proof
and finite universes are unchanged.
