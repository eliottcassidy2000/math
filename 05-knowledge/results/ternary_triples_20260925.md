# A primitive-triple filtration whose denominator drops by three

**Status: PROVED elementary conjugacy, contraction, inverse guards, and scope
obstructions; FINITE-EXACT independent controls.** Integer Collatz convergence
remains **OPEN**. The useful new object is an oriented primitive-triple graph
with an intrinsic contracting region and an exact integer boundary. It is not
an identification of Collatz with the unguarded three-branch Berggren tree.

## Inheritance and concept board

The closest proved mechanism is the odd-root chart in
[the Pythagorean semicircle note, Theorem 1](collatz_mod6_20260917_pythagorean_semicircle.md):
`(s,t) -> (st,(s^2-t^2)/2,(s^2+t^2)/2)`, with odd coprime roots s,t.
Its Theorem 8 already identifies a different halving operation: Gaussian
angle halving of a primitive triple with square hypotenuse. We do not identify
that operation with arithmetic division by two.

The canonical hostile is the same note's loss of the rational half-angle
coordinate, together with the signed Collatz cycles at 5 and 17. The corrected
near miss is treating a ternary address as a Collatz edge. The current proved
[THM-4057, CW18--CW25](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md)
already gives exact Calkin--Wilf/Berggren ordinal transducers and their guards;
we use none of those address maps as a Collatz convergence premise. The
least-used sidecar here is the sign of the even leg, which selects a rational
number rather than its reciprocal.

The classical Pythagorean parameterization is also recorded in the primary
paper [Janičková--Csókási, section 2](https://arxiv.org/html/2304.05230).
All parameterization and dynamical claims used below are proved here. That
paper is background for the classical chart, not a source for this Collatz
lift, contraction, or root claim. The repository's
[Pythagorean source sidecar](../reference/CORE-PAPERS-PYTHAGOREAN.md) supplies
the wider literature route.

| Live concept | Exact coordinate | What it resolves or does not resolve |
|---|---|---|
| Primitive triple | two signed square gaps | recovers both coprime roots |
| Division by two | valuation of an actual linear numerator | retains the full halving count |
| Division by three | primitive common factor of numerator and denominator | decreases an intrinsic geometric level |
| Integer boundary | hypotenuse minus signed even leg equals 1 | recovers ordinary odd Collatz |
| Guarded join | two actual forward clocks and a smaller certified input | certificates transport through the chart |
| Reflected sign | positive numerator guard | preserves cycles and a rational zero boundary |

Anchor: an exact structural lift with a decreasing rank. Niche: the invariant
prime-to-three denominator that determines which boundary is reachable.
Wildcard: an explicit infinite family of root certificates beyond every
fixed geometric level, and noninteger fixed points that exclude an oversized
claim.

## 1. The exact objects and inverse

Let R be the positive rational numbers whose reduced numerator s and
denominator t are both odd. Let P consist of triples `(A,B,C)` satisfying

    A>0, C>0, A odd, B even, A^2+B^2=C^2, gcd(A,|B|)=1.

B is signed. Include `(1,0,1)`, the sole degenerate primitive point. Define

    Phi(s/t) = (st, (s^2-t^2)/2, (s^2+t^2)/2).             (T1)

This is a bijection R -> P. Its inverse is intrinsic:

    s=sqrt(C+B),  t=sqrt(C-B),  x=s/t.                    (T2)

To prove it, for odd coprime s,t the displayed coordinates are integral,
primitive, and satisfy the Pythagorean equation. Every common odd prime
dividing st and `(s^2-t^2)/2` would divide both s and t; the odd leg excludes
a common factor two. Conversely, `C+B` and `C-B` are positive odd coprime
integers with product `A^2`: their gcd divides `2C` and `2B`, and primitivity
gives `gcd(C,B)=1`. Each is therefore a square. Their positive square roots
are odd and coprime, and their product is A. This gives (T2) and uniqueness.
When B=0, primitivity forces A=C=1.

For B>0 this is an ordinary positive primitive Pythagorean triple. Negative
B means the same triangle with the two square roots interchanged. Thus
forgetting the sign identifies x with `1/x`. It is not a valid quotient of
the dynamics below.

## 2. The primitive splitter is exactly one factor of three

For `sigma in {+1,-1}`, define the accelerated rational operation

    U_sigma(x)=(3x+sigma)/2^v2(3x+sigma).                  (T3)

The valuation of a nonzero rational with odd denominator is the valuation
of its reduced numerator. On the positive sheet require `3x+sigma>0`.
For plus this is automatic. For minus it is the real guard `x>1/3`.

For a reduced odd spinor `(s,t)`, put

    k=v2(3s+sigma t),  g=gcd(3s+sigma t,t)=gcd(3,t),
    (s',t')=((3s+sigma t)/(2^k g), t/g).                 (T4)

The equality for g is exact, not an upper bound: since `gcd(s,t)=1`,
`gcd(3s+sigma t,t)=gcd(3s,t)=gcd(3,t)`, which is either 1 or 3 even if t
contains a high power of 3. The numerator is even, so k>=1. After the
displayed cancellations the roots are again positive odd coprime integers.
Consequently (T4), transported by Phi, is an isomorphism of functional
graphs with (T3) on its legal positive domain. At t=1 it is precisely the
ordinary accelerated odd Collatz map, with the chosen sign.

This is a concrete separation of the two operations. First form the actual
linear numerator `3s+sigma t`; divide it by two until it is odd; cancel the
common factor three from both coordinates when present. The two
cancellations commute because three is odd. The ternary operation here is
primitive cancellation, not a count of three children. Arithmetic halving
acts on one spinor coordinate, not on all three triangle coordinates and
not by Gaussian angle halving.

The full triple has enough information to perform the operation: recover
the square roots in (T2). In particular `g=3` iff `9|(C-B)`. The companion
[main bridge](ternary_bridge_20260925.md) gives an entirely intrinsic
coordinate formula for the same edge.

## 3. The power-of-nine gap family contracts to the integer boundary

Write

    t=3^r q,  3 does not divide q.

Every legal edge obeys

    r' = max(r-1,0),             q'=q.                   (T5)

Thus the prime-to-three part q is an invariant. The family

    P_3={P in P: C-B=9^r for some integer r>=0}           (T6)

is exactly the image of the positive odd rationals with reduced denominator
a power of three. Its level is determined by the gap; its level-zero
boundary `C-B=1` consists of

    Phi(n)=(n,(n^2-1)/2,(n^2+1)/2), n positive odd.

For plus, a point of level r reaches this integer boundary after exactly r
accelerated edges. No convergence assumption is needed for that entry.
More strongly, every step while `3|t` has an intrinsic decreasing integer
rank: the hypotenuse C.

When `3|t`, put `d=t/3`. Equation (T4) becomes

    s'=(s+sigma d)/2^k,       t'=d,       k>=1.           (T7)

For plus,

    2C' <= (s+t/3)^2/4+t^2/9 < (s^2+t^2)/3=2C/3,       (T8)

because the difference in the strict inequality is

    [3(s-t)^2+4t^2]/36 > 0.

Hence **C'<C/3** on every positive plus edge with `3|t`. Equality is
impossible. This also gives the simpler root-sum bound
`s'+t' <= (s+t)/2`. The contraction holds for any odd t divisible by three,
not only powers of three; when q>1 it stops at that different invariant
denominator boundary.

For a legal positive minus edge, `s>t/3` and

    2C' <= (s-t/3)^2/4+t^2/9 < (s^2+t^2)/4,

because the difference is `st/6+t^2/9>0`. Thus **C'<C/4** for minus while
the positive guard holds. It does not assert that all positive rational
points have a next minus edge.

There is now an exact equivalence, with its proof boundary exposed:

> Every positive odd integer reaches 1 under U_+ iff every point of P_3
> reaches `(1,0,1)` under the lifted plus operation.

The forward implication uses exactly r contracting steps to reach an
integer and then the assumed integer certificate. The converse restricts
to level zero. This is a reduction to the integer boundary, not a proof
that the boundary is exhausted. It avoids substituting a large coverage
count or a ternary-tree address for an actual Collatz certificate.

## 4. Inverse edges make the ternary guard visible

For a positive target `u/v in R` and an integer k>=1, the inverse candidate
with this exact halving count is

    x=(2^k u-sigma v)/(3v),    require x>0.               (T9)

Reduce this fraction. It has odd numerator and denominator, and direct
substitution gives `3x+sigma=2^k(u/v)`, so the count k is exact. Its forward
edge is therefore certified without searching an orbit.

On P_3, if v is divisible by three, `3` does not divide u, so the numerator
in (T9) is not divisible by three. The denominator becomes exactly `3v`:
every reverse step outside the integer boundary raises the level by one.

If v=1, the parent is an integer **iff**

    2^k u=sigma mod3.                                   (T10)

If that condition fails, the parent has denominator exactly three. For
plus: if `u=1 mod3`, integer parents require even k; if `u=2 mod3`, they
require odd k; if `3|u`, there are none. For minus, swap the two parity
requirements; multiples of three still have no integer parent. These are
the three exact residue cases. The ambient rational inverse operation is
always available once positivity holds, but the integer Collatz inverse
is only its guarded boundary part.

Thus an unguarded inverse word can be realized as a rational/triple route
without being a route on positive integers. Its increasing denominator
level records exactly where that interpretation was lost. Once such a
reverse route has left the integer boundary it cannot return to it by
further reverse steps.

## 5. Orientation, prime support, and finite-bit obstructions

The positive-triangle example

    (117,44,125) -> (3,-4,5) -> (1,0,1)                 (T11)

has spinors `(13,9)->(1,3)->(1,1)`. It crosses the diagonal s=t after the
first step. Discarding the sign of B substitutes the different rational 3
for `1/3`: `(3,+4,5)` goes to `(5,12,13)`, whereas `(3,-4,5)` goes directly
to `(1,0,1)`. Therefore unsigned triangles cannot determine this map.

The invariant q also gives an infinite hostile to extending the root claim
to all primitive triples. For every k>=3,

    x=1/(2^k-3)                  satisfies U_+(x)=x.      (T12)

The denominator is prime to three. Its primitive triple is
`(t,(1-t^2)/2,(1+t^2)/2)`, `t=2^k-3`, and lies outside P_3. The first is
`(5,-12,13)`. In contrast its reciprocal `(5,+12,13)` reaches the root.
This explains both the selected gap family and the necessary orientation.

Modulo three, exactly one leg of every point in P is divisible by three,
and its hypotenuse is not. In the roots, `3|A` iff `3|s` or `3|t`; when
neither root is divisible by three, `3|B` since both root squares are 1
modulo three. The two roots cannot both be divisible by three. Merely
knowing which leg carries three does not distinguish numerator from
denominator; the oriented square gaps do.

The usual dyadic seam is also exact. Both odd squares are 1 modulo eight,
so `4|B`, and for s!=t,

    v2(|B|)=v2(|s-t|)+v2(s+t)-1.                         (T13)

One of `s-t,s+t` has valuation one and the other at least two. This leg
valuation is not the Collatz halving count `v2(3s+sigma t)`.

No fixed number of low bits of the triple determines that full count. For
any h>=1, take t=9 and `s=2^k-3` with k>=h+3. Then
`3s+t=3*2^k`, so the count is k, while all coordinates of Phi(s/t) have the
same residues modulo `2^h` as k varies. The division by two in (T1) is why
we retain more than h input bits in this witness. This excludes a local
residue table that outputs the entire count; it does not exclude a
transducer with an unbounded carry counter or access to the full triple.

## 6. A counter family and the sign controls

The filtration yields a simple infinite family of complete plus root
certificates. For r>=1 put

    t_r=3^r,  s_r=(2*8^(r-1)+3^r)/5.                    (T14)

The numerator is divisible by five because 8 and 3 agree modulo five.
The quotient is positive, odd, and prime to three. We have s_1=1 and
`s_r=8s_(r-1)-3^(r-1)` for r>=2. Thus its exact accelerated exponent word
is r-1 repetitions of 3, followed by 1:

    (s_r,t_r) -> ... -> (s_1,t_1)=(1,3) -> (1,1).

This gives root certificates at arbitrarily high intrinsic levels without
forward search. These starting points are rational, not new certificates
for arbitrary integer inputs. Inverse construction changes the starting
point; it cannot be used to choose the future of a prescribed integer.

Minus supplies both kinds of hostile that the chart must preserve. Starting
with a positive primitive triangle,

    (783,56,785) -> (45,-28,53) -> (3,-4,5),
    (29,27)     -> (5,9)       -> (1,3),                 (T15)

the next rational numerator is zero. The positive-sheet guard rejects it;
zero has no odd accelerated normalization. At the integer boundary,
the distinct minus cycles with odd minima 1,5,17 survive unchanged. Hence
off-boundary geometric contraction cannot decide what happens on the
boundary, and it does not make the sign irrelevant.

## 7. Certificate interface and exact controls

The source is the legal positive odd-rational graph, and the target is the
oriented primitive-triple graph P. Phi and (T2) are explicit inverses.
They preserve every actual step, its sign, exact halving count, and common
future. On the integer boundary n grows strictly with
`C=(n^2+1)/2`, so the guarded join relation
`0<y<x`, `U^a(x)=U^b(y)` from
[the common-future decoder](creative_decoder_20260925.md) transports exactly
with the same clocks and smaller hypotenuse. A known boundary certificate
can be prepended by the finitely many contracting off-boundary edges.

The destroyed information in an unsigned quotient is reciprocal
orientation; the destroyed information in a residue quotient is the full
halving count; the destroyed information in an unguarded inverse tree is
integer-domain membership. The required data are respectively signed B,
the full square roots or an exact carry counter, and the gap level.

Reproduce with

    python 04-computation/experiments/ternary_triples_20260925.py
    python -O 04-computation/experiments/ternary_triples_20260925.py

Both agree with [the frozen output](ternary_triples_20260925.out). All checks
remain active under `-O`. The explicit universes are:

- 6,636 odd coprime pairs `s<=255,t<=127`, including 12,716 legal signed
  edges compared with an independent rational arithmetic implementation;
  all 3,127 edges with denominator divisible by three pass strict rank.
- An independent, Euclid-free enumeration of all 95 oriented primitive
  triples with `C<=300`, including the degenerate point.
- 638 power-of-nine-gap starts with `r=0..6,s<=255`: every plus entry clock
  is exact. All 14,984 positive signed inverse candidates with `k=1..12`
  pass their exact exponent and ternary-level guards.
- The explicit orientation crossings, noninteger fixed points for
  `k=3..20`, counter family `r=1..20`, and low-bit hostiles `h=1..32`.
- The minus zero boundary and all three indicated integer cycle controls.

The infinite statements are consequences of the proofs above, not
extrapolations of these finite universes. The contraction-to-boundary
problem is resolved for P_3; generation of root certificates for every
point on the integer boundary remains OPEN.

Independent proof audit: the parallel Berggren lane checked sections 2--6,
including the primitive gcd, both strict contraction differences, inverse
phases, counter recurrence, orientation and sign controls, with no issues.
