# A whole arithmetic cell escapes every bounded-depth smaller-child portfolio

2026-10-04. **PROVED:** the mixed congruence obstruction for every forward
and inverse depth bound, its exact density, and the scoped family-cover
consequence. **PROVED:** an inherited depth-growing rule supplies valid
smaller dependencies inside each obstructed cell. **FINITE-EXACT:** the
independent controls below. Universal grounded Collatz coverage remains
**OPEN**. No literature-priority claim.

Artifacts: [program](../../04-computation/experiments/collatz_partition_cover_20261004.py)
and [saved output](collatz_partition_cover_20261004.out).

## 1. Inheritance and the kind of cover being considered

The closest mechanisms are the guarded common-future family lift and
bounded-inverse completeness argument in
[the inverse-shield note](collatz_join_shields_and_lifts_20261004.md), and
the distinction between a decreasing family and its completed members in
[the frontier family compiler](frontier_family_compiler_20261004.md).
The elementary smaller-predecessor guards in
[Proposition3 of the connectivity note](collatz_connectivity_from_rigidity_20261001.md#4-route-3-the-minimal-counterexamples-two-sided-trap-realised-on-the-minus-sheet)
are also retained, including `n=4 mod9`; the literal identity is checked
again below rather than importing that note's broader unaudited claims.
[THM-4512, coefficient-descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
keeps exact and coarse valuation guards separate; it does not assert
universal entry into the contracting classes. The
[binary/ternary fusion](collatz_binary_ternary_guard_fusion_20261004.md)
combines selected families but retains a necessary remainder.

Canonical hostile: `27 ->41 ->31` does not descend below the original27.
Corrected near miss: a class whose members all have a smaller obligation
is not thereby a class whose members all have completed root proofs.
The least-used sidecar is the pair of depths: source-to-join and
child-to-join. Both survive any partition into residue classes.

Anchor: a sound portfolio of family obligations. Niche: find a whole
arithmetic cell that every bounded-depth portfolio misses. Wildcard: test
whether growing word templates can enter that very cell. The concept board
is **original source / smaller child / exact words / two depth bounds /
mixed congruence / grounded suffix**. The result concerns the existence of
actual smaller joins, not just the failure of a particular decoder.

## 2. What a family cover must prove

Write `U(n)=oddpart(3n+1)` on positive odd integers, with `U(1)=1` for
arithmetic bookkeeping. A checked dependency at an original source n is

    U^r(n)=U^s(m),       0<m<n,                        (1)

with actual finite valuation words retained on both sides. Empty words
are permitted. This encompasses direct descent (`s=0`), the fixed
half-child, ternary sibling constructions, and variable-child suffix cuts
when their exact guards and original-source inequality hold.

There are two distinct ways to use (1).

* With a supplied first-hit certificate for m, retain its suffix at the
  common endpoint and prepend the source word. If an arithmetic word
  includes root loops, trim them before exporting a first-hit certificate.
* In a complete induction proof, establish a covering collection of such
  dependencies for every positive odd n>1. A least non-home n would have a
  smaller home child m, a contradiction. If the cover has an effective
  selector with finite checks, recursive calls on the strictly smaller
  positive odd m construct a certificate and cannot form a circular chain.

The second use requires the entire cover, not membership in just one paid
family. For example every positive `n=5 mod8` descends in one odd step,
but `29 ->11` leaves that same family. Its descent proof alone does not
supply the missing home proof for11. Likewise, a checked smaller child of
a later larger frontier need not be below the original request.

If a finite residue cover overlaps, disjoint refinement by a common modulus
can retain each applicable dependency and its sidecars. Refinement preserves
valid guards; it cannot create a missing join. The theorem below places a
precise obstruction on every uniform pair of depth bounds, even if the
guards are more general than arithmetic progressions.

## 3. The mixed congruence obstruction

**PROVED.** Fix nonnegative integers R,S. Let n>1 be a positive odd integer
satisfying

    n=-1 mod2^(R+1),       n=0 mod3^S.                (2)

Then there do not exist a positive odd m<n and integers

    0<=r<=R,  0<=s<=S

such that `U^r(n)=U^s(m)`.

All inverse valuation exponents are unrestricted. There is no source-size
threshold, no coefficient-only approximation, no assumption that n or m
reaches1, and no finite search extrapolation in this statement.

The binary guard fixes the first R actual valuations to one:

    U^r(n)=(3/2)^r(n+1)-1.                            (3)

If s<=r, each inverse map `G_a(x)=(2^a x-1)/3`, a>=1, is at least
`G_1(x)` on positive inputs. Composing s such inverse steps from (3)
therefore gives

    m >= (3/2)^(r-s)(n+1)-1 >= n,

which contradicts m<n. It remains to rule out s>r.

Let the actual child word have length s, total valuation D>=s, and carry B:

    F_v(X)=(3^s X+B)/2^D,
    B=sum_(i=0)^(s-1) 3^(s-1-i) 2^(a_1+...+a_i).

Inverting its common endpoint with (3) gives the affine expression

    m = lambda*n+c,
    lambda=3^(r-s)2^(D-r),
    c=lambda-(2^D+B)/3^s.                             (4)

Since D>=s>r and `3^S` divides n, the number `lambda*n` is an integer.
The actual child m is an integer, so c is an integer too. Formally c is
the inverse image of

    z_r=(3/2)^r-1=F_(1^r)(0)                          (5)

under the same child word. The following sign lemma is the missing fact.

**PROVED inverse sign lemma.** For r>=1, an integer inverse image of z_r
under a positive valuation word of length s>r must be strictly positive.
For r=0, no nonempty positive valuation word has an integer inverse
image of z_0=0.

To prove it, read the inverse maps from (5). Initially the reduced
denominator is exactly `2^r`. If an inverse state ever has a factor3 in
its reduced denominator, that factor cannot later disappear. Indeed, for
`x=b/(2^d3^h)` with h>0 and 3 not dividing b, the numerator of
`G_a(x)` remains a3-unit while the denominator acquires another3.
An integer final result thus forces every intermediate state to be dyadic.

Suppose a nonintegral intermediate has dyadic denominator exponent d>0
and value at least `z_d=(3/2)^d-1`. The next exponent a>=1 decreases its
dyadic denominator exponent to `max(d-a,0)`. Since the current value is
nonnegative,

    G_a(x)>=G_1(z_d)=z_(d-1)>=z_max(d-a,0).

This invariant starts at (5), so within at most r inverse steps the
denominator has cleared and the integer value is nonnegative. If that
integer were zero, another inverse would introduce the irreversible
denominator3. Because s>r, another step is required. Therefore the first
integer is positive; every remaining integral inverse of a positive
integer is positive. This proves c>0. Starting from z_0=0, the first
inverse is already `-1/3`, which proves the r=0 case.

Finally each a_i>=1 gives

    B>=3^s-2^s,       2^D>=2^s,
    (2^D+B)/3^s>=1.                                  (6)

If lambda<1, (4)--(6) force c<0, contradicting the lemma. If lambda>=1,
the same lemma gives `m=lambda*n+c>n`, again impossible. This completes
the proof, including arbitrarily large inverse exponents. The empty/empty
case is just m=n and was already excluded.

For instance,27 belongs to the cell R=S=1. The larger inverse child581 is
perfectly legal: `27 --(1)-->41 <--(4,3)--581`. Its intercept at formal
source zero is5 and its slope is64/3. The theorem excludes a *smaller*
child, not all incoming paths or all common futures.

## 4. Exact family-cover consequence

The two moduli in (2) are coprime. Hence (2) is one complete odd residue
class modulo `2^(R+1)3^S`, with exact natural density

    1/(2^(R+1)3^S) among all integers,
    1/(2^R3^S) among odd integers.                    (7)

Excluding the possible root1 changes neither density. These are densities
of a specified obstruction to bounded-depth joins, not densities of
nonconvergence or unresolved inputs under every method.

Some diagonal examples are:

| R=S | Blind source class | Relative density among odd integers |
|---:|---|---:|
|1|3 mod12|1/6|
|2|63 mod72|1/36|
|3|351 mod432|1/216|
|4|1215 mod2592|1/1296|
|5|1215 mod15552|1/7776|
|6|16767 mod93312|1/46656|

Consequently, no portfolio whose source-to-join depths are uniformly
bounded by R and child-to-join depths by S can furnish a smaller-child
dependency for every positive odd source. This holds even for infinitely
many candidate rules, arbitrary modular filters and unbounded individual
valuation exponents, provided the two numbers of odd steps stay bounded.
A finite bank of fixed word pairs is a special case. Splitting their
domains into finer integer families cannot overcome the obstruction.

This is not a prohibition on a finite *program* with recursive or
parameter-dependent word lengths. Nor does it forbid a direct supplied
root certificate, a proof outside the smaller-child model, or further
actual exploration beyond the two bounds. Finite explicit completed
exceptions can remove only finitely many points from the positive-density
cell. The moduli in (2) are convenient sufficient guards; no minimality
claim is made for their precision.

## 5. Positive signals beyond the precise obstruction

The ternary guard in (2) matters. For any R>=0, put

    n=3*2^(R+1)-1,       m=2^(R+2)-1.

Then n has R initial valuation-one steps, but m<n and `U(m)=n`. More
generally the mixed guard `n=2 mod3` gives the smaller inverse source
`m=(2n-1)/3`. Thus an all-ones binary prefix alone cannot rule out every
mixed-modulus smaller-child rule.

The inherited second guard gives another concrete division of labor:

    n=4 mod9,       m=(8n-5)/9<n,       U_(1,2)(m)=n.

For positive odd n in that class, m is a positive odd integer and the two
valuations are exactly1,2. For `n=2^p-1` with `p=5 mod6`, periodicity of
powers of2 modulo9 gives precisely this guard. It therefore supplies a
smaller dependency for an infinite family with arbitrarily long initial
one-runs; at p=5 this is `27 --(1,2)-->31`. The odd exponent classes
`p=1,3 mod6` have source residues1,7 modulo9 and are not covered by that
particular rule. This is inherited applicability, not a new proof that
every member of any of these families reaches1.

Conversely the inherited reset-at-least-three switch enters the very
cell (2) when its word length is allowed to grow. For R>=1 choose the
exact source word

    left=(1^R,3),       right=(1^(R-1),2,1).

For each n in that exact dyadic cylinder,

    U_left(n)=U_right((n-1)/2),       0<(n-1)/2<n.      (8)

This follows from the retained guarded reset identity; direct affine
composition gives equal lengths and cost difference one. Intersect its
source cylinder with `n=0 mod3^S` using CRT. The resulting infinite
progression is a subfamily of (2), and every member has the checked
smaller dependency (8). Its source depth is R+1, honestly beyond the
excluded bound. This construction proves that the obstructed cell is not
an intrinsic non-home class. Its members still need their child proofs
or participation in a complete well-founded cover.

For example at R=S=1, the subfamily starts at51:

    51 --(1,3)-->29 <--(2,1)--25.

The source and child need two steps each, whereas (2) excludes only joins
with at most one step on either side. No fixed bounded table can replace
all such growing-length cases merely by adding more residue filters.

## 6. Exact computation and its independent completeness control

The program tests R,S from0 through8:81 cells and243 positive sources,
using parameters0,1,7 in each cell and replacing root1 by the next member
when necessary. It explores all forward endpoints through R and all
smaller inverse ancestors through S. No cap is imposed on inverse
valuation exponents. Instead, if h inverse steps remain and the present
positive odd value is x, any eventual child m<n requires

    2^h(x+1)<3^h(n+1).                                (9)

Every inverse map is at least G_1, and iterating G_1 proves (9). This gives
an exact finite upper bound on the next exponent. If the inequality fails,
no continuation within the remaining depth can be smaller. The same
completeness mechanism is inherited from the inverse-shield work.

For all75 odd sources3 through151 and four bound pairs `(0,3),(1,2),
(2,3),(3,4)`, a separate literal enumeration of every positive odd label
m<n gives exactly the same join sets:300 independent controls, including
all depths and any root padding present in those finite tests.

A second universe consists of all65535 positive valuation words with
total cost at most16. A separate carry sum verifies the recurrence;
48559 contracting source-all-ones/child-word pairs have a negative,
nonintegral formal-zero intercept, as required. The1024 smaller rational
inverse controls explicitly replay the inverse maps. These finite checks
support the proof; they are not its all-depth justification.

The depth-growing repair has288 exact controls across R=1..8, S=0..8,
and parameters0,1,17,10^30. The dyadic-only hostile, the larger ancestor581,
the nonclosed descending family `29 ->11`, and the unpaid original-source
path `27 ->41 ->31` are explicit controls. Four malformed depth inputs
are rejected. The inherited mod9 predecessor is also replayed at all ten
exponents `p=5 mod6` from5 through59; the remaining odd exponents3 through61
are checked to miss that guard. All checks survive Python optimization.

Reproduce:

    python -X utf8 -B 04-computation/experiments/collatz_partition_cover_20261004.py
    python -O -X utf8 -B 04-computation/experiments/collatz_partition_cover_20261004.py

Normal and optimized stdout agree after LF normalization; SHA256:
`e79548beadcc14cbe15ce258902885a1075a43c9ad63ccca3a834b9f2005c8c6`.

Connection contract: source = an original positive odd integer and a
portfolio of exact common-future dependencies; target = a family-cover
proof or a grounded certificate; map = guarded word composition and
strict descent to the chosen child; preserved = actual source and common
future; lost by residue labels alone = child identity, two word depths,
original-source inequality and root proof. The new obstruction isolates
the depth coordinate that a universal smaller-child portfolio cannot
bound uniformly. Proving that an effective unbounded construction always
finds a grounded dependency remains open.
