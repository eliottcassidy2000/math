# Sharp fibre margin at a signed-pair virtual wall

**Status: PROVED elementary corollary + FINITE-EXACT + INDEPENDENTLY AUDITED,
after the explicit correction below.** This separate addendum leaves the frozen source and output of
`overnight9_20260906_lrc_virtual_pair_wall` unchanged.

Retain its actual three-sheet model, and define the tail-fibre clearance

    H_T(y)=max_(j=0,1,2) min_(w in T) ||w(y+j)/3||.

For any signed pair frequency d=|u+-v|/3 supplied by two ternary-unit tails,

    H_T(y) >= | ||dy|| - 1/3 | / 2.                  (1)

This is the optimal continuous lower profile on the entire distance domain
0<=delta=||dy||<=1/2. It is uniform over all three distinct positive
ternary-unit tails, all phases, and each of their three signed pairs.

**Proof.** Put alpha=|delta-1/3|/2, so 0<=alpha<=1/6. If alpha=0 the
statement is trivial. At strict danger threshold alpha, an
order-three grid has at most one bad point: its open danger arc has length
2alpha<=1/3. If both pair tails have different owners, their nearest tooth
integers have the same nonzero combined residue used in the primary proof,
but the signed error after division by three is now strictly below 2alpha.
The circle-distance function is one-Lipschitz in both directions, so full
coverage would force |delta-1/3|<2alpha=|delta-1/3|, a contradiction. Some
sheet therefore has every tail clearance at least alpha, proving (1).
At alpha=1/6 a half-integer effective phase is inactive for the strict
mask, so the nearest-integer step still has no tie. The zero-alpha boundary
requires no owner argument.

For lower-half sharpness let A<B, A=B=2 modulo three, and
take T=(A,B,A+B), y=1/(A+B), d=(B-A)/3. Its literal sheet minima are

    ( A/[3(A+B)], A/[3(A+B)], 0 ),

so delta=(B-A)/(3(A+B)) lies in (0,1/3) and equality holds in (1).
Primitive choices are retained by gcd(A,B)=1;
the ratio A/(A+B) is dense in (0,1/2) before this restriction. Dividing by
the common gcd produces primitive tails, rescaling y and permuting the
three sheet labels; both H_T(y) and ||dy|| are preserved.

For upper-half sharpness let A=C=2 modulo three, C>=4A>0, and take
T=(A,C,A+C), y=1/C, d=(2A+C)/3. This d is supplied by the signed sum of
the first and third tails. Writing r=A/C gives the literal sheet minima

    (r/3, 0, r/3),   delta=(1+2r)/3 in (1/3,1/2].

Thus equality also holds throughout a dense set of the upper half.
Indeed the ratios (3i+2)/(3j+2) are dense in (0,1/4); for either family,
dividing the tails by their common gcd g and replacing y by gy preserves
delta and permutes the three sheet values because 3 does not divide g.
The two dense equality sets and their limits rule out every continuous
lower profile that is everywhere at least (1) and strictly larger anywhere.
This does not claim attainment at every individual real delta. The upper
endpoint is attained by primitive T=(1,4,5), y=1/4, d=2, H=1/12.

In particular T=(11,17,28), d=2 and y=1/28 give

    ||dy||=1/14,  H_T(y)=11/84.

Thus the sharp uniform tail clearance at the primary virtual wall is
**11/84**, substantially above the body threshold 1/14. This is an actual
max-min over a common sheet, not a sum of marginal tail measures.

At the original strict danger threshold 1/14, full spoil implies H_T(y)<1/14.
Hence the signed-pair necessary band strengthens to

    4/21 < ||dy|| < 10/21.                           (2)

Both inequalities are strict, and both constants are sharp. The primary
report proves lower-wall sharpness. For the upper wall take

    T=(9N-1,42N-1,51N-2), y=1/(42N-1), d=20N-1,
    N>=1, N!=5 modulo eleven.

The tails are primitive ternary units: their common gcd divides eleven,
with eleven occurring exactly in the excluded residue. Their sheet margin
is (9N-1)/(3(42N-1))<1/14, so they are genuinely fully spoiled, while

    10/21-||dy||=11/[21(42N-1)] -> 0.

**Correction lineage.** The pre-audit addendum proved the valid lower bound
max(0,(1/3-delta)/2), but incorrectly called it globally sharp: its equality
family only covered delta<1/3. The two-sided Lipschitz inequality above
provides a strictly stronger bound for delta>1/3. For instance
T=(2,11,13), y=1/11, d=5 has delta=5/11 and H=2/33, attaining the new
positive upper-half bound where the old bound was zero. This corrects the
global sharpness assertion, not the original lower inequality or the
11/84 virtual-wall constant. The original report, source, and outputs are
preserved outside the repository with `_pre_audit` suffixes, for historical
audit only; they must not be filed as current proof artifacts.

## Arbitrary coprime h: a strict family extension

Keep the primary d>=31, gcd(d,42)=1, b=3d+1, c=3d+2 and the eight base
body frequencies {d,2d,3d,4d,14,14b,14c,1}. Now let h be any integer larger
than all eight, coprime to their lcm L=42dbc, and put C=base union {h,2h}.
The six coprime virtual blocker counts prove that some y in V_d^+ is
body-safe, without specifying its address.

Every such survivor has all body clearances except d strictly above
1/14. For a coprime residual frequency z, equality on
y=(14k+1)/(14d) would imply

    z(14k+1)=+-d modulo 14d.

Reducing modulo d and using gcd(z,d)=1 forces d to divide 14k+1. Then
y is an integer multiple of 1/14, strictly blocked by body frequency 14,
a contradiction. The other d multiples have clearances 2/14,3/14,4/14.
Frequency d is consequently the unique owner at a positive left endpoint.
Formula (1) supplies a sheet with strict tail margin at least 11/84.
Continuity, or the explicit interval

    eta=min(1/(14d),
       min_(z in C,z!=d)(||zy||-1/14)/(2z),
       3(H_T(y)-1/14)/(2 max T)) > 0,

preserves that sheet and pays every body coordinate to the right of y.
Its interior points are strict witnesses for the full thirteen-speed row.
Thus the primary unbounded closure extends from h=1 mod L to **every
sufficiently large h coprime to L**, with no prescribed h residue and no
common body dilation. Requiring h>91^6 retains the same necessary numerical
cross-height filters. This still makes no assertion of actual decoder entry.

The companion standard-library verifier evaluates the literal max-min
fibre on 56 tail triples and every rational phase of reduced denominator
at most 24, checks both exact equality families and the sharp owner-band
sequences, and reconstructs strict
component witnesses for four arbitrary coprime h values that are not
1 modulo L. The universal proof is the threshold-parameter argument above.

The independent proof and separate piecewise/interval verifier are
`overnight9_20260906_lrc_virtual_pair_wall_audit.md` and its companion source.
They also record the scope overlap with incoming THM-4448: its uniform cone
closes the explicit h=1+mL subfamily, while the present arbitrary-coprime-h
wall criterion applies even to some h just above the base maximum.

**Frozen reproduction:** 36,019 active gates; normal and optimized output
are byte-identical LF. Source SHA256:
`62ee68aff4af17b0d3951e7a13bd93c3ccd5245c22e66d600b7973036bd65bd0`.
Output SHA256:
`9c05b9b64e0e93ccd025722aa26f46614b5bf447a589ca945c2815b3d91f6603`.

**Filing:** root integrated these audited artifacts in the ninth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
