---
id: THM-3316
title: "Prime right-boundary interpolation forces micro-staircase scalar rigidity"
status: >
  PROVED + INDEPENDENTLY FINITE-EXACT AUDITED.  At every prime modulus p,
  a residue vector covering only the p^2 right-boundary affine tests is
  already a scalar ramp.  Consequently the full micro-staircase blockers are
  exactly the p scalar ramps.  This is a finite residue/cell theorem, not a
  speed-to-residue lift and not an LRC theorem.
source: root/creative-jacobian-lrc/2026-08-03
depends_on:
  - THM-363-lrc-scalar-gauge-reindexing
  - THM-364-lrc-scalar-ramp-cell-blocking
related:
  - HYP-1823-lrc-scalar-gauge-quotient
scripts:
  - 04-computation/lrc_prime_boundary_rigidity_micro_boundary_moments_agent_20260803.py
  - 04-computation/lrc_microstaircase_n13_prime_unit_orbit_micro_n13_rigidity_agent.py
outputs:
  - 05-knowledge/results/lrc_prime_boundary_rigidity_micro_boundary_moments_agent_20260803.out
  - 05-knowledge/results/lrc_microstaircase_n13_prime_unit_orbit_micro_n13_rigidity_agent.out
---

# THM-3316 -- prime right-boundary interpolation forces scalar rigidity

**PROVED + INDEPENDENTLY FINITE-EXACT AUDITED.**

## Statement and scope

Let `p` be prime and let

```text
v=(v_1,...,v_(p-1)) in F_p^(p-1).                       (1)
```

For `a,s in F_p`, call `(a,s)` right-boundary blocked by `v` if some
`1<=i<=p-1` satisfies

```text
a i+s v_i in {0,-1}.                                   (2)
```

If all `p^2` pairs `(a,s)` are blocked, then there is one `m in F_p` such
that

```text
v_i=m i                         for every i.             (3)
```

Consequently, in the full open-cell micro-staircase system at modulus `p`,
the full blockers are exactly the `p` scalar ramps `(mi)_(i=1)^(p-1)`.

The implication from a full blocker to `(2)` is literal.  Choose
`0<epsilon<1/(p(p-1))` immediately to the right of `alpha=a/p`; its open-cell
pattern is

```text
floor(p {i alpha})=a i mod p.                            (4)
```

This theorem classifies a finite residue/cell model.  It neither realizes an
arbitrary speed tuple as `v` nor preserves endpoint, gcd, cell-width, or
physical-time data.  It proves no case of the Lonely Runner Conjecture.

## Proof

Subtract the scalar ramp `v_1 i` and write

```text
w_i=v_i-v_1 i,                 w_1=0.                   (5)
```

The right-boundary predicate is preserved: the expression for `w` at `(a,s)`
is the expression for `v` at `(a-sv_1,s)`.  Suppose for contradiction that
`w` is nonzero.  Put

```text
S={l:w_l!=0},       h=|S|,       q_l=w_l^(-1).           (6)
```

Since `w_1=0`, one has `1<=h<=p-2`.

Fix `j in S` and take `a=-j^(-1)`.  A zero coordinate cannot block: for
`i notin S`, the value `-i/j` is neither zero nor minus one, because the
second equality would force `i=j in S`.  A supported coordinate `l` blocks
precisely at

```text
s=lq_l/j             or             s=(l-j)q_l/j.       (7)
```

After the bijective rescaling `y=js`, condition `(2)` says that, for every
`j in S`, the `2h` values

```text
{lq_l,(l-j)q_l:l in S}                                  (8)
```

cover `F_p`.

Let

```text
X={lq_l:l in S},
A(Y)=product_(l in S)(l-Yw_l),
L=product_(l in S)l.                                    (9)
```

Fix `y notin X`.  The first alternative in `(8)` is then unavailable.  For
every `j in S`, coverage supplies some `l in S` with

```text
y=(l-j)q_l,                   hence j=l-yw_l.            (10)
```

The `h`-element multiset `{l-yw_l:l in S}` therefore contains all `h`
distinct members of `S`; it equals `S`.  Taking products gives

```text
A(y)=L                         for y notin X.             (11)
```

For `y in X`, one factor in `A(y)` is zero, and conversely, so

```text
A(y)=0                         for y in X.                (12)
```

Write `r=|X|`.  Every label in `S` is nonzero, hence `L!=0`, and
`1<=r<=h<p`.  Equations `(11)--(12)` imply

```text
sum_(y in F_p) A(y)=(p-r)L=-rL !=0.                     (13)
```

But `deg A=h<=p-2`.  The sum over `F_p` of every monomial of degree at most
`p-2`, including the constant monomial, is zero.  Therefore the left side of
`(13)` is zero, a contradiction.  Thus `w=0` and `(3)` follows.

Conversely, THM-364 proves that every scalar ramp blocks every shift on every
open cell.  This completes the classification.  For `p=2`, normalization
already leaves no nonzero coordinate and the same conclusion is immediate.

## Constructive witness

The proof finds a missed boundary cell.  If `w` is normalized and nonzero,
then some `y notin X` has `A(y)!=L`; otherwise the same finite-field sum gives
the contradiction `(13)`.  For that `y`, some `j in S` fails `(10)`.  Then

```text
a=-j^(-1),                  s=yj^(-1)                   (14)
```

is an explicit unblocked right-boundary candidate.

## Equality and failure boundary

Primality is load-bearing: it supplies inverses and makes the nonzero integer
`|X|<p` nonzero in the coefficient field.  At composite modulus `14`, the
normalized vector

```text
(0,1,0,...,0)                                             (15)
```

covers all `196` analogous right-boundary tests while missing interior
micro-staircase cells.  Thus boundary sufficiency cannot be transferred to
the LRC(14) target without CRT/valuation and interior-cell sidecars.

Over a general finite field, the same calculation forces low-degree power
sums of `X` to vanish; already the zeroth moment says that `|X|` is
characteristic-torsion.  The prime-field contradiction is exactly the case
where this leakage is impossible.

## Verification record

The independent companion reconstructs `(4)` for every prime through `31`,
exhausts all normalized vectors at `p=3,5,7`, checks a constructive witness
for every nonzero vector in those universes, and tests dense deterministic
banks through `p=19`.  It verifies `5,780,507` scalar-gauge candidate
equivalences, audits the weaker critical-support moment route through `p=101`,
and retains `(15)` as a composite hostile.  Normal and optimized runs equal
the frozen transcript.

```text
script sha256 (LF)  774288510f30ff36cc0b2f863f0f9d3af809c0ff958d8fc035b311ff7dea5cd8
output sha256 (LF)  2eeff12d06d8588f94641877d1900388bfe2605214228f6bbc1b5ce85a31b83b
```

The separate `p=13` four-bit SAT classifier gives an orthogonal exact audit:
after scalar gauge and `F_13^*` first-support normalization, CaDiCaL and
MapleSAT independently make all eleven nonscalar branches UNSAT.  It is no
longer the strongest proof source, but it checks the formerly missing case by
a disjoint mechanism.

## Dependencies and consumers

THM-363 supplies the full-cell gauge interpretation and THM-364 the positive
scalar-ramp direction.  The converse above is self-contained at the
right-boundary level.  The next legitimate consumer is a prime-grid
speed/residue lift retaining width and endpoint data; absent that map, this
theorem is not an LRC result.
