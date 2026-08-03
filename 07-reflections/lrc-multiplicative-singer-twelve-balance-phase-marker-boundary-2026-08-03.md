# Multiplicative Singer twelve-balance and the phase-marker boundary

**Status:** VERIFIED-EXACT DISCOVERY CANDIDATE; independently suitable for
auditing the promoted THM-3255 package, but not itself a theorem reservation or
proved physical LRC handoff.

## Inheritance pass

- **Closest proved mechanism.** THM-3246 supplies the exact all-dilation
  `168`-owner Hodge word and positive numerator quadratics. THM-3252 proves
  full charged cyclicity after arranging the completed Hodge word as a
  `13 by 13` additive coefficient matrix. THM-3253 proves that the actual
  positive masses give a nonsingular additive matrix for every integer
  dilation and all `8,064` Singer gauges.
- **Canonical hostile.** THM-3246's twelve negative seam owners have the
  correct cardinality for a tempting line or phase object but meet every
  Singer vector line in at most one point. Cardinality is not geometry.
- **Corrected near miss.** MISTAKE-166 forbids moving an LRC claim through
  Singer-looking parameters without specifying the carrier and forgotten
  fibre. Here the additive plane matrix and multiplicative regular
  representation are deliberately kept distinct.
- **Least-used relevant sidecar.** THM-2839 says prime-to-`p` total mass is a
  unit on a `p`-power cyclic orbit. The orbit here has order
  `168=2^3*3*7`, so positivity and nonzero total mass do not imply a group-ring
  unit. The exact failure is instead detected by the `C_12` norm quotient.
- **Meta-patterns used.** “The same representation is not the same carrier,”
  “Pay algebraic rank, then audit the physical action,” and “Controlled
  forgetting and unlabeled quotients require a sidecar.”

## Exact candidate statement

Let `q=(q_j)` be THM-3246's Hodge word, and let

```text
P_g(X)=sum_(j=0)^167 N_(g,j) X^j,
Q(X)=sum_(j=0)^167 q_j X^j,
S_12(X)=1+X+...+X^11.
```

For every integer `g>=1`, in `Q[X]`,

```text
gcd(P_g(X),X^168-1)=gcd(Q(X),X^168-1)=S_12(X).       (1)
```

Equivalently, both cyclic convolution matrices on `C_168` have exact rank

```text
168-deg(S_12)=157.                                    (2)
```

Their zero characters have exact orders

```text
2,3,4,6,12,                                           (3)
```

whose Euler degrees total `1+2+2+2+4=11`. There are no other zeros, and
the statement is invariant under every Singer gauge

```text
j |-> b+a*j,      a in (Z/168Z)^*, b in Z/168Z.       (4)
```

This cleanly refutes the bold but false reframe “additive Singer-plane
nonsingularity implies multiplicative `C_168` cyclicity.” The same positive
word has additive rank `13` in THM-3253's matrix and multiplicative rank
`157<168` here.

## Why twelve: uniform norm pushforward

In THM-3234's deterministic field, `Norm(alpha)=6` has order `12`. Thus

```text
F_169^* -> F_13^*,      alpha^j |-> 6^j               (5)
```

is exactly exponent reduction modulo `12`; each norm fibre has `14` points.
The exact class sums are independent of the residue `r`:

```text
sum_(j=r mod 12) q_j       = 1/296352,
sum_(j=r mod 12) N_(g,j)   = 120960g^2-528g+2.         (6)
```

Hence the pushforward of either word to `C_12` is uniform. Every nontrivial
multiplicative character of `F_13^*` vanishes, giving the factor `S_12`.
Gauge `(a,b)` merely permutes norm phases by

```text
r |-> a*r+b mod 12.                                    (7)
```

The negative seam set

```text
{0,1,2,3,4,5,162,163,164,165,166,167}                 (8)
```

is a complete transversal modulo `12`. Thus every Hodge phase slice has
exactly thirteen positive coefficients and one negative coefficient.

## The stronger phase-disintegration reframe

The rank defect is created by summing phases, not carried by an individual
phase. For every residue `r mod12`, both restricted words

```text
Q_r(X)=sum_(j=r mod12) q_j X^j,
P_(g,r)(X)=sum_(j=r mod12) N_(g,j) X^j                 (9)
```

satisfy

```text
gcd(Q_r,X^168-1)=1,
gcd(P_(g,r),X^168-1)=1        for every integer g>=1. (10)
```

So each of the twelve Hodge slices and each fourteen-point **positive mass
slice** already spans the full `168`-dimensional multiplicative regular
module. Cancellation between twelve individually cyclic objects creates the
eleven missing norm characters in their sum.

This full multiplicative rank does not transfer back to the additive matrix.
Each norm conic has active row/column counts `(8,7)` or `(7,8)`, alternately,
so any `13 by 13` matrix supported on one phase has additive rank at most
seven. The whole word reverses the behavior:

```text
one norm phase: multiplicative rank 168, additive rank <=7;
sum of 12 phases: multiplicative rank 157, additive rank 13.             (11)
```

This is the sharp representation-switch boundary exposed by the computation.

## Minimal phase marker and positivity tax

For a norm residue `r`, let

```text
e_r(j)=1_(j=r mod12),
m_r=e_r-(1/12)1.                                      (12)
```

The centered marker `m_r` has Fourier support exactly the missing eleven
characters, so its orbit is the exact `11`-dimensional direct complement to
the rank-`157` word. It has zero sum and must be signed: after scaling by
`12`, its signs are `14` positive and `154` negative.

The Boolean marker `e_r` is nonnegative and has orbit dimension `12`; it
contains the eleven missing characters plus the trivial character. Therefore

```text
span(C_168.w)+span(C_168.e_r)=Q[C_168],                (13)
```

with exactly a one-dimensional overlap. This is optimal among nonnegative
markers: every nonzero nonnegative vector has positive total mass, hence a
nonzero trivial Fourier component, and a marker filling all eleven missing
characters must consequently have orbit dimension at least `12`.

## Exact proof certificate

The companion verifies `(6)` directly from all `168` exact coefficients. For
each cyclotomic order dividing `168`, it reduces the symbolic quadratic
`P_g` modulo `Phi_d` and takes the gcd in `Z[g]` of all remainder
coordinates. The forced-zero orders are exactly `(3)`. On every remaining
order the coordinate gcd is a unit, except

```text
d=8,24: 504g-1,
d=1:    60480g^2-264g+1.                              (14)
```

The first has no integer root and the second is positive for `g>=1`. The same
test on every phase restriction gives no forced factor: quotient characters
see the positive fibre total, orders `8,24` again require `g=1/504`, and all
other coordinate gcds are units. Direct hostile controls cover
`g=1,2,3,17,18,57,58,169`, six nontrivial gauges, and all `8,064` quotient
gauge permutations.

## Typed connection and loss ledger

```text
source:       owner weights on F_169^*
map A:        zero-complete and arrange on additive F_13^2
consumer A:   THM-3250/3252/3253 multiplicity matrix
map M:        identify F_169^* with C_168 and take multiplicative translates
consumer M:   circulant/group-ring cyclicity
quotient:     norm pushforward C_168 -> C_12
preserved:    total mass, cyclic gauge covariance, norm-fibre totals
destroyed:    additive row/column placement, completion point, endpoints,
              target/ancestry, and physical current semantics
sidecar:      one labelled norm phase, centered or Boolean as appropriate
hostile:      full additive rank coexists with multiplicative defect 11;
              one multiplicatively full phase has additive rank at most 7
```

No tournament is manufactured: the intrinsic object is a cyclic group ring
with a norm quotient, not a binary orientation.

## Next physical test

The cheapest decisive test is no longer another determinant. Attach the
owner label `j mod12` to THM-3246's **actual endpoint/target/ancestry record**
and ask whether one positive norm-phase slice is closed under the physical
placement map. In exact terms, construct or refute a commutative square

```text
positive owner mass slice  --physical placement--> endpoint/target current
          | norm phase                              | retained phase
          v                                         v
       one C_12 fibre  -------------------------> lawful target fibre.
```

The first hostile should use the seam transversal `(8)`: every norm phase
contains exactly one negative Hodge owner but fourteen positive mass owners.
If endpoint placement mixes two norm phases, record the first mixed pair and
the lost endpoint coordinate; if one phase survives, test its full
multiplicative orbit against the actual target rather than another abstract
module. Until that square exists, equations `(1)--(14)` remove an algebraic
rank gate but do not exclude an LRC row or decrement `LRC(14)`.

## Reproduction

```text
python3 04-computation/lrc_multiplicative_singer_twelve_balance_discovery_20260803.py
python3 -O 04-computation/lrc_multiplicative_singer_twelve_balance_discovery_20260803.py
```

The stored transcript is
`05-knowledge/results/lrc_multiplicative_singer_twelve_balance_discovery_20260803.out`.
