# Reversed sibling guards and a mirrored run family

**PROVED:** the retained signed receipts pay their original source on their
exact guards; reversing the last-two run operation gives an infinite family.
**FINITE-EXACT:** the declared proposal box and its new parameter cells.
**CONDITIONAL:** each emitted child still needs its own ROOT certificate.
**OPEN:** universal parameter coverage and universal Collatz completion.

## 1. Inheritance and the omitted orientation

Use the actual parameter chart from
[child compression](collatz_child_compression_synthesis_20261007d.md):

    E(t)=924745897+2^32 t, n(t)=2^E(t)-1, t>=0.

The closest proved mechanism is the
[signed sibling ladder](collatz_complement_ladders_20261007c.md). The new
[positive parameter bank](collatz_parameter_cover_20261007e.md) searched only
positive sibling gaps. Its complement therefore deserves a test with the
same source, cost scale and word order, but the opposite ladder orientation.
The hostile is a formally matching pair whose adjusted final valuation is
nonpositive. The corrected near miss is to omit that terminal reserve.
The least-used sidecar is the source-side sibling gap with its actual cost.

The live board is source ownership, signed gap, native head, terminal
reserve, prefix-free coverage and child obligation. The anchor is the
uncovered parameter set, the niche is reversed affine symmetry, and the
wildcard is a reflected version of the positive run surgery. It gives an
actual operation, rather than an analogy based only on the number of bits.

## 2. The signed rule and its guard

For an ordered positive valuation word u, write
`F_u(x)=(P_u x+B_u)/Q_u`, with `P_u=3^length(u)` and `Q_u=2^sum(u)`.
Let D>0, r<0, `length(v)=length(u)+D` and `sum(v)=sum(u)-2r`.
The exact carrier test is

    F_v(y)=S^r(F_u(x)) whenever x+1=3^D(y+1), S(z)=4z+1.

Negative powers denote the rational inverse affine map; integrality is not
assumed from this identity. It is provided by the guard below.

Set `A=sum(u)`, `R=-2r`, `C=A+R`. The native source cell is

    x = -(3B_u+2^A)/(3P_u) modulo 2^(C+1).             (1)

All denominators inverted here are odd. This cell fixes the entire actual
head u and requires the next source valuation c to satisfy `c>=R+1`.
Thus the child terminal `c+2r` is positive, and

    F_(c+2r)(F_v(y)) = F_c(F_u(x)).                    (2)

For a supplied source `n=2^(e+1)h-1`, h positive odd, the initial run gives
`x=2*3^e*h-1`. The smaller source `n'=(n+1)/2^D-1` has initial run e-D
and run endpoint y. A sufficient uniform cutoff retained by the compiler is

    e >= max(2A, D+2C, D+1).                          (3)

Both initial runs and both heads then remain strictly above their respective
starting sources after the initial expansion. The elementary estimate is
`(3/2)^e/2^A >= (9/8)^A>1` when e>=2A; the child uses e-D and C.
All prefix carries are positive. Odd endpoint integrality forces every
intermediate valuation to be native. There is no earlier ROOT crossing in
either head. Adding the actual terminal in (2) gives an authenticated common
future with `0<n'<n`.

On the Mersenne chart h=1, solve `(x+1)/2=3^(E-1)` in its exact power-three
subgroup. Formula (1) gives a unique exponent class modulo `2^(C-2)` and
hence, for every retained rule, one parameter class of precision `C-34`.
The compiler keeps the cutoff from (3). Every saved finite rule has cutoff
zero in this very large starting exponent chart.

## 3. Concrete complement fills

The proposal universe is the parameters 0..1023 missing the frozen positive
bank, deletions D=1..32, source-head depths 0..128 and retained precision1024.
There are136 successful negative-gap proposals. Their saved ordered words
are recompiled, rather than trusting stored guard labels. Each successful
proposal proves an infinite arithmetic class, not just its seed.

Examples:

| Parameter guard | Deletion D | Sibling gap | Source depth | Effective guard cost C |
|---|---:|---:|---:|---:|
| t=393 modulo2^10 |4|-1|18|44|
| t=9 modulo2^14 |4|-1|21|48|
| t=6 modulo2^455 |4|-1|236|489|

The last row is a separately declared deeper target, not part of the
128-depth census. The known positive-head search misses t=6; this actual
source admits a reversed ladder at depth236. Its tiny density is not
reported as broad coverage.

These137 signed cells add exactly the rational mass stored in
[the output](collatz_signed_parameter_fill_20261007e.out), approximately
0.001507608827920, to the frozen positive bank. The combined finite binary
mass is approximately0.199776235881211. Counting uses a prefix antichain;
the original labelled rules and all their child obligations remain saved.

For the t=9 rule, deleting its two-bit terminal reserve admits the literal
source `351379561727494872721350587442095838424203263`. Its head is native,
but its actual next valuation is1, giving adjusted child valuation -1.
The production consumer rejects it. Endpoint carrier equality alone would
have accepted the wrong domain.

## 4. The mirrored all-run operation

The positive bank's run surgery comes from the common fixed point -1 of
`F_1(x)=(3x+1)/2` and `S(F_2(x))=3x+2`; these affine maps commute.
Reflect the operation across a gap-minus-one identity. If v=p,2, set

    u_k=u,1^k,       v_k=p,1^k,2,       k>=0.          (4)

Then `F_(v_k)(y)=S^-1(F_(u_k)(x))` with the same deletion D. Indeed
`S F_2 F_1^k = F_1^k S F_2`. The source/partner length difference and cost
difference are unchanged. The operation retains the ordered carrier and
the marked terminal reserve; it does not commute arbitrary Collatz words.

The t=9 row has precisely this form. Its k-th native parameter cell has
precision14+k. Different k specify different maximal runs of ones after
the same head u, each ending at a reset of at least3, so the cells are
pairwise disjoint. Their total parameter mass is

    sum_(k>=0) 2^(-14-k) = 1/8192.

The k=0 cell is already in the finite signed bank. All k>=1 lie in a raw
prefix cell disjoint from both finite banks, certified by exact prefix
comparison. Thus the reflected family adds exactly **1/16384**.
The all-run family in the positive bank begins after the old15-letter head
with valuation5, while this reflected family begins there with valuation8;
the two infinite families are disjoint.

The generic compiler's conservative height cuts remain part of pointwise
acceptance. They do not alter the stated natural density: truncate at k<K,
where only finitely many points are cut; the entire remaining union lies
in the uncompleted prefix u,1^K, whose parameter cylinder has mass
`2^(-12-K)`. Its integer count through X is at most that mass times X plus1.
Let X tend to infinity and then K to infinity. The same shrinking-prefix
argument permits fusion with a finite ternary bank by CRT.

Connection contract: source = a marked signed receipt; target = an
unbounded run-indexed bank; map = (4); preserved predicate = exact common
future and strict child payment; lost by a bare gap/type = native phase,
terminal reserve and child identity; sidecars = (1)--(3) and the ROOT port;
cheapest test = exact carrier composition plus a literal missing-reserve
hostile. This is lossless reuse of a small operation in a second orientation.

## 5. Reproduction and scope

Run `python -B 04-computation/experiments/collatz_signed_parameter_fill_20261007e.py`
and repeat with `-O -B`. Add `--rediscover` to repeat the declared finite
proposal census and the separately frozen t=6 target. The
[script](../../04-computation/experiments/collatz_signed_parameter_fill_20261007e.py),
[raw heads](collatz_signed_parameter_fill_20261007e.json) and
[output](collatz_signed_parameter_fill_20261007e.out) retain the evidence.
Independent modular powers check both source/child words and final valuation
adjustment; modest literal sources test the first-hit receipt consumer.
Wrong precision, a changed partner, wrong exact types and missing reserve
are rejected. The signed guard/core proof was independently peer audited.

The remaining task is an authenticated rule on the exact residual parameter
set and a mechanism grounding emitted children. Repeating a marked type or
improving its coverage fraction alone supplies neither of these obligations.
