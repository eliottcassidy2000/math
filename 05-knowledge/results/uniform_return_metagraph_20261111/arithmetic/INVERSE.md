# Exact inversion of the positive five-core charts

**PROVED scoped algorithm / FINITE-EXACT controls.** No all-target coverage,
source-return closure, optimal tournament order, or new primality result is
asserted. This algorithm improves inversion of an inherited grammar.

## 1. Inherited objects and the coordinate retained

Use the positive observer from
[the five-core packet](../../chanpak_five_core_bridges_20261110/cover/positive_raw_signatures.json)
and its signed independent compiler from
[Session251](../../session251_creative_completion_20261111/ternary/compiler.py).
Their exact normalized flank response is

    Phi(A,B)=h+d(g(A)+g(B)-2)+LAB+iA+oB+s(A mod2)+e(B mod2),

where g(0)=g(1)=1, g(2q+1)=g(q), and
g(2q)=g(q)+q+(q mod2). For the actual ninety cuts,
h>=d>0, L>d, i,o,s,e>=0 and s,e<=d. Here L=d+K(1), with K(1)>0.
The two signed and positive
representations are inherited proved identities, not inferred from testing.
The old direct quadratic scan can invert these charts, but its practical
cost conceals large target witnesses. The retained coordinate here is the
entire pair of binary addresses, not merely the represented odd integer.

The live concepts are source fibres, positive observers, dyadic addresses,
short-flank bounds, and the distinction between target coverage and return.
The closest proved mechanism is the positive seven-coordinate chart. The
hostile examples are the axis failure at2501 and the first-three-core
short-flank minimum29 at8820901. The discarded shortcut is to extrapolate a
small fixed flank from low targets. The less-used sidecar is the binary
recursion for g; it supplies an exact inverse without such an extrapolation.

## 2. The one-coordinate inverse

Induction gives 1<=g(B)<=B+1 for every B>=0; for odd positive B the
stronger g(B)<=B holds. In the even induction, an even q uses
g(q)<=q+1, whereas an odd q uses g(q)<=q; either gives g(2q)<=2q+1.
The case B=0 is separate. For positive d,c and 0<=e<=d, consider

    d g(B)+cB+e(B mod2)=y.                                  (1)

There are no solutions if y<d. Otherwise every solution lies in

    max(0,ceil((y-d-e)/(c+d))) <= B <= floor((y-d)/c).         (2)

Split B by parity. The complete two branches are

    even B=2q:   (y,d,c,e) -> (y,d,2c+d,d),
    odd  B=2q+1: (y,d,c,e) -> (y-c-e,d,2c,0).               (3)

Both preserve d>0, c>0, 0<=e<=d. The implementation directly tests an
interval of at most nine integers, and otherwise recurses by(3). Every
returned integer satisfies(1), and every possible integer is in exactly
one branch. These facts prove soundness and completeness by induction on
the decreasing upper bound in(2), with direct testing as the base case.
The parity split itself is inherited, for example in
[Session220 section5](../../session220_binary_envelopes_and_selectors_20261107/selectors/PROOF.md#5-inherited-curvature-obstruction-specialized-to-f2).
The contribution here is its quantitative interval stopping rule and
composition with the positive short-flank bound, with complete controls.

Here is a useful deterministic complexity bound, in arithmetic operations.
The width of(2), before clipping at zero, is at most

    yd/c^2 + e/c.

If c>=sqrt(yd), then c>=d whenever y>=d, and this width is at most2.
At recursion depth j, c has at least doubled j times, while y has not
increased. Therefore the recursion has depth
O(1+log(max(1,sqrt(yd)/c))) and at most
O(1+sqrt(yd)/c) nodes. Each leaf tests at most nine candidates and each
g evaluation takes O(log(y+2)) elementary integer steps. This bound
includes all returned solutions; it is not an average-case heuristic.
The integers have O(log(y+d+c+2)) bits; no unit-cost bit-complexity claim
is made.

## 3. The free axis and the complete two-coordinate inverse

The c=0,e=0 case has infinitely many binary presentations because
g(2q+1)=g(q). Remove trailing binary ones repeatedly. The representative
is even, including0 for the strings consisting only of ones; 0 and1
are consequently the same normalized axis case. The remaining equation

    d g(2q)=d g(q)+dq+d(q mod2)=y

is(1) with c=e=d. The routine returns every even representative, and does
not incorrectly claim to return all infinitely many unnormalized strings.
The actual chart coefficient guards ensure that a free axis has exactly
this form. An arbitrary nonnegative seven-tuple without that guard is not
silently treated as an actual tournament chart.
Nor can the opposite bounded-flank scan reintroduce a positive free-axis
address: if B<=K on such an axis, then
LB^2<=N-d=d(g(B)-1)<=dB, impossible for B>=1 because L>d.

For N>=h and A,B>0, positivity gives LAB<=N-h. Thus

    min(A,B)<=K=floor(sqrt((N-h)/L)).                        (4)

For each A=0,...,K set

    y=N-h-d(g(A)-2)-iA-s(A mod2),  c=LA+o,

and solve(1) in B. Repeat after interchanging the two flanks and their
coefficients, then remove identical pairs. Axes are included by A=0 in
one orientation. Equation(4) proves this is the complete normalized
chart inverse. An optional max_short parameter restricts this to all
solutions with min(A,B)<=max_short; absence in that restricted probe is
never promoted to absence from the complete chart.

For an actual chart, the subtracted base is nonnegative: g(A)>=1 and
h>=d. Hence every recursive y<=N. Summing the node bounds over A>=1,
where c>=LA, gives

    O(K + (sqrt(Nd)/L) log(K+1))

nodes, plus the two axis terms. For fixed chart coefficients this is
O(sqrt(N) log(N+2)) nodes and O(sqrt(N) log^2(N+2)) elementary integer
steps including g evaluation and sorting the output pairs. Finite fixed-bank inversion has the same
asymptotic bound. This is an exact target-dependent size cap, not a proof
that every odd target has a witness.

## 4. Independent finite controls and the reset interface

`chart_inverse.py` checks42240 one-coordinate cases by a direct integer
scan. For all90 actual cuts and all128 odd targets through255, it agrees
with the old signed complete generator on all4402 normalized pairs. It
also checks large bounded-short-flank witnesses with the signed formula.
The receipt is `chart_inverse.json`.

`ray_inverse.py` obtains and independently checks witnesses for
102*2^k+19, k=1,...,25, using short-flank cap16. The last target is
3422552083; one witness is the S4 cut2 chart A1,B89224797. The saved
addresses and exact counters are `ternary_plus_four.json`. These are25
literal grammar witnesses. They provide neither a fixed prefix recurrence
nor a blue/black path from a prescribed source, and no uniform return is
inferred. The parent reset lane must still retain its marked source and
legal next operations.

For the useful C3 cut1 chart, the positive observer is(5,3,12,6,3,2,1).
A target102D+19 forces A odd and B even modulo3. Writing A=2a+1,B=2b,
the exact inverse equation reduces to

    g(a)+g(b)+16ab+4a+11b+(b mod2)=34D+4.

For fixed b this is precisely the one-coordinate inverse with d=1,
c=16b+4,e=0 and target34D+4-g(b)-11b-(b mod2). Thus a changing
short flank can be searched without scanning the other large address.
Whether a bounded grammar of these changing addresses covers every k
remains OPEN.

## 5. Concurrent inheritance and independent scope

The newly integrated
[infinite-birth packet](../../infinite_birth_families_20261111/arithmetic.py)
independently uses the Omega>=3 source bound and the same exact
A135244/A002496 boundary. Its source sieve through10^6 certifies target
prefix10000 and its quadratic sieve reaches root100000. Our larger
prefixes, tighter integer cap, consecutive-prime maximum family and Lucas
certificates complement those checks; the shared elementary bridge is
credited as a concurrent result, with the earlier sharp-bound mechanism
already present in the inherited sessions.

That packet also supplies32 fixed odd-root classes, odd part1,...,63,
with actual tournament seeds for eventual square-plus-one rays. This is
an infinite-family result for finitely many root classes, whereas our2047
short-flank minima are finite complete-chart results. Neither implies
unrestricted all-root coverage. Distinguishing these quantifiers keeps
the positive infinite families useful without losing the remaining gap.
