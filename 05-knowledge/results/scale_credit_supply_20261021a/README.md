# Scale phases, linear parent supply and the early-arrival boundary

2026-10-21. **PROVED:** the complete fixed-parent annulus fibre, its
three-phase signed interaction law, two linear-size native parent banks,
and their source-owned finite-horizon identity. **FINITE-EXACT:** all342
parents in those banks at N=8192 have authenticated first-arrival receipts;
121 of159 residual intervals are paid by this bank, and a disjoint38-owner
frontier supplement completes them. **OPEN:** a uniform estimate for the
remaining R', the153N annulus bound, and universal Collatz convergence.

## 1. Inheritance and the unchanged source

This continues [the94 external credits](../residual_external_credit_20261020k/README.md)
and [the exact interval flow](../annulus_residual_flow_20261020k/README.md).
The inverse formula and full braid clock are inherited from
[THM-4690](../../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md)
and [native phase lifts](../session18_native_phase_lifts/README.md).
The six-band charge law and source norm153N are inherited; this packet
does not rebrand the inverse map or negative-term regrouping as new energy.

The named source is a positive odd primitive integer (not5mod8). P is
the canonical primitive map; its edge clock includes every actual
B-strip, B(x)=4x+1. It is not a count of fully odd U-steps. V is the
same627-element grounded frontier, with max V=63665. All intervals use
actual first-arrival P-depth, and charges retain the full clock sign.
The joint mod15/mod35 register comes from the literal integer source.

Closest mechanism: same-parent lifts64s+35/49. Canonical hostile:
negative charge at one scale may be positive at the next. Corrected
near miss: fixed old parents provide only bounded one-step supply, and
large eventual supply need not be available at an early cut. Least-used
sidecar: the actual arrival event of the old parent. Board: inverse
fibre / scale phase / source bands / parent arrival / temporal capacity /
unchanged remainder R.

## 2. A complete bounded fibre, not a linear source supply

For a in{1,2}, every native inverse is

    s(u,a,j)=[2^a(3u+1)4^j-(2^a+3)]/9,  j>=0.       (1)

Integrality selects exactly one j modulo3: 4 generates{1,4,7} modulo9,
and the required unit belongs to that subgroup. The positive integral
inverse has the stated native branch; when j=0 this follows directly
from (3s+1)/2^a=u, and when j>=1 it is3mod8 or1mod8 respectively.
The edge is P(s)=u, with clock a+2j. Advancing j by3 gives

    s' =64s+c_a,   c_1=35, c_2=49.                 (2)

**Exact annulus fibre.** For every N>=1 and primitive u<4N, each a gives
exactly one s in[4N,256N). The first legal j is at most2 and gives
0<s<256N. If it is below4N, advance(2) until it first enters. Its
preceding integer r satisfies r<=4N-1, so64r+c_a<256N. Once inside,
its next successor is at least256N+c_a and lies outside. The affine
offset is essential; an interval-ratio argument alone would be incomplete.

Consequently a fixed bank of K parents has at most2K immediate inverse
addresses in any annulus. For all627 frontier parents this is at most1254,
and exactly1254 before sign filtering when N>=15917. One must still
exclude s in V when declaring its first-arrival depth to be1.

All charges have Q at most3. Indeed L_r=U_r-U_(r+1): a B-chain ending
in an upper term telescopes to U_r, and an all-lower chain to a difference
of two U-vectors. Thus |beta(v,w)|<=6, and the absolute interaction
budget of any fixed K-parent immediate bank is at most12K. Such a bank
cannot supply a positive multiple of N as N grows. This does not prove
that the actual missing interaction requires that much credit. Nor does
it bound longer-depth predecessors of V or all pairs sharing a terminal.

## 3. All three scale phases must be retained

For new s, its charge A_(64N)(s) depends only on s mod3 and the six bands
[4N,8N),...,[128N,256N). The values are respectively

    0, U_r, -L_(r+2), U_r, L_r, U_r.               (3)

At N'=64N the unique child of the same(u,a) is s'=64s+c_a. Integer
band boundaries are preserved because0<c_a<64. Write R(x,y)=(-y,x-y).
Then

    A_(64N')(s')=R^(c_a) A_(64N)(s),
    A_(64N')(u)=A_(64N)(u).                         (4)

The edge clock increases by6. Hence the signed old/new interaction is
three-periodic in the scale index and its three consecutive values sum
to zero, since I+R+R²=0. An initially negative value cannot remain
strictly negative at every consecutive64-scale. The later three-scale
return preserves a chosen negative credit, but does not multiply its
number. The complete1254-slot phase table at N=16384 is retained in the
receipt and determines its whole64-scale ray, not all integer scales.

The full-denominator modular reader evaluates(1) modulo9m before dividing
by9. Along(2), the joint mod105 register has period dividing6: the second
iterate is translation65c_a modulo105 and its third repetition is zero.
No independent mod15/mod35 choices are substituted for this actual source.

## 4. Two uniformly signed banks with growing parent sets

Both following families use native first branch a=2. Set c=ceil(3N/2).

| Family | Original parent range and residue | New source | Edge clock | Actual pair interaction after arrival |
|---|---|---|---:|---:|
| A | c<=u<2N; u mod72 in{9,27,63} | (16u+3)/3 | 4 | -3 |
| B | N<=u<c; u mod72 in{17,35,71} | (64u+19)/3 | 6 | -6 |

The residues enforce primitive odd u and the required integrality. All
parents lie in[N,2N), so A_N(u)=L_(u mod3). In A, u=0mod9,
8N<s<16N and s=1mod3; the charges are L_0=(1,2) and U_1=(-1,-1).
In B, u=8mod9,16N<s<32N and s=0mod3; the charges are L_2 and-L_2.
The even edge clocks preserve the tail sign in the polarization, giving
-3 and-6 exactly. This proves the sign without relying on a finite test.

For each interval[lo,hi) and residue set T, the exact count is

    sum_(r in T) [floor((hi-1-r)/72)-floor((lo-1-r)/72)].  (5)

Each family has between3floor(N/144) and3ceil(N/144) members. When144
divides N, each has exactly N/48, so their eventual signed magnitude is
3N/16. This is potential supplied interaction, not a horizon-independent
available budget. For each fixed residue progression u=r+72t, the joint
(u,s) register modulo105 has exact parameter period35. That parameter
is different from the scale index in section3.

For N>=7959, every new source here exceeds max V. If u has first-arrival
terminal t, clock k and depth h, its child has the same terminal, clock
k+4 or k+6, and depth h+1. At a finite horizon d the exact selected sum is

    C_N(d)=-3#{A parents arrived by d-1}
           -6#{B parents arrived by d-1}.            (6)

This identity is unconditional: if u never enters V, neither does its
child, and that pair contributes zero at every finite horizon. No
convergence assumption is hidden in(6). A production use of a nonzero
term requires the supplied actual parent receipt. In particular when
N>63665 none of these parents is in V, so C_N(1)=0 despite the linear
source count. The missing quantitative premise is timely arrival.

Following a fixed pair through section3 gives phase values(-3,0,3) for
A and(-6,3,3) for B. Rebuilding the two parent bands at each new N is a
different operation from lifting the same old parent. These are not
interchangeable sources or clocks.

## 5. One finite bulk/early-boundary completion

At N8192, each family has171 members. All342 parent and child routes are
authenticated to the original V, retaining terminal, both clocks/depths,
unsigned charges and mod15/mod35. Their eventual capacity is513 units of3;
409 units are available by depth16, but only8 at depth1. The latest child
arrival is64. These are actual finite receipts, not a universal entry law.

The inherited interval bank had16983 completed positives and159 residual
ones. Keep the former certificate untouched. Using only the new linear
bank, the exact interval primal completes121 of159. The weak-dual cut
{1} proves the remaining38 are unavoidable:46 residual intervals cover
depth1, against8 units then present. Later abundance cannot pay that cut.

Enumerate immediate inverse slots of the admitted old frontier parents,
retain beta=-3, and exclude every already-owned pair. The first38 such
owners have actual depth1. Together with the342 bulk owners they pay
all159 residual intervals. Every original pair is unique. A beta=-6
owner supplies two explicitly marked3-unit slots, not two copies of its
physical interaction. The saved temporal allocation uses each slot at
most once at any horizon and is checked on all257 cuts.

For C the selected old/new sum, J the entire old sibling interaction and
R every other original pair term, the exact ledger remains

    R'=R-C, J'=J+C, I=R+J=R'+J',
    J'<=U+C<=0,  R'+U+C=R+U.                       (7)

All external owners are old/new, hence distinct from the original
all-new sibling addresses. Shared endpoints are permitted; repeated pair
addresses are not. Interactions between these owners and every other
source remain in R. Scalar matching need not preserve terminal, but it
never moves a charge vector between terminals. No full8192 R or energy
census is performed. The proved inequality I<=R' alone is not a numerical
153N estimate; the additional condition R'<=M remains to be established.
The finite completion improves the ownership explanation;
it cannot improve physical I by recycling its negative terms.

## 6. Validation and interfaces

`linear_pairs(N,routes)` exports source-authenticated dict rows.
`receipt.json.case8192.linear_rows` is the342-parent bank;
`frontier_supplement` is the38 early owners. Literal ROOT suffixes may
be supplied by another compiler; this packet uses the authenticated
first-arrival/grounded-frontier interface and does not claim an all-family
ROOT result. Given these already authenticated rows, the capacity-only
helper `slots` rejects repeated physical pair addresses before
splitting a beta=-6 capacity into its two marked units.

Small controls completely enumerate inverse fibres for N1..16 and check
both families for N1..288. The fixed627-parent table retains all six
tested scale phases, full-denominator registers and the zero three-phase
sum. Independent charge formulas, dependency hashes, literal type guards,
source reconstruction, exact capacity primal/dual, and every ownership
segment are checked. No external solver or unresolved orbit assumption
is used.

```text
python -B 05-knowledge/results/scale_credit_supply_20261021a/check.py --verify
python -B -O 05-knowledge/results/scale_credit_supply_20261021a/check.py --verify
```

Generated receipt:157146 checks; verification adds one saved comparison.
LF SHA256: `5b9b010d19c5b1d38f380d1e31567d44750f5277e786a5d6b6c446a7785ff64a`.

Independent peer proof/API and optimized saved-receipt audit: **PASS**.
The affine-offset fibre, three-phase sign law, linear counts, literal
parent/child arrival seams and121+38 ownership completion were checked.
The all-N statement retains the parent-arrival event; no full8192
remainder estimate or universal153N conclusion is inferred.
