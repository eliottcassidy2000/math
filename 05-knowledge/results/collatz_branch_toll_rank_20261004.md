# What the four-step barrier measures, and a rank that permits larger children

2026-10-04. **PROVED:** the sharp branch-delay hierarchy, commutator/guard
identity, proper signed graph rank, complete local-minimum classification,
ranked family transport, and the conditional controller criterion below.
**FINITE-EXACT:** the explicitly declared audits and frozen-source search.
**OPEN:** cancellation of every positive critical state and universal
Collatz convergence. No literature-priority claim.

The useful change of objective is to rank **basin-preserving proof moves**,
including moves to larger integers. This produces a genuine well-founded
reduction to an explicit critical set. It does not make that critical set
disappear. The four-step inverse obstruction also extends to a sharp
all-height hierarchy whose delays equal the sibling construction's ternary
precision requirements.

## 1. Inheritance and the changed questions

This session pulled commits `10b5f45b53` and `23fe3869e1`. The
[incoming supplied-state synthesis](supplied_state_extension_limits_and_repairs_20261004.md)
already certifies all239 historical requests from ROOT1 in its finite
experiment. The239-source list below is a frozen comparison universe,
not a claim that those sources remain uncertified. The incoming
[half-child obstruction](half_child_extension_obstruction_20261004.md)
proves that the exact fixed half-child ansatz fails at7 at every length.
Its slope-spectrum theorem and the
[child-cut law](supplied_state_join_extensions_20261004.md) are inherited.

A later clean checkpoint integrated `b58b73e38c`, `891b067700`, and
`f29c414ec2`: the descent-class correction, whole-cell depth obstruction,
and [composed dependency families](collatz_partitioned_completion_20261004.md).
Section8 applies the incoming construction to this note's critical set;
the final audit independently verifies its exact all-height guards.

Closest mechanisms: [the reset shield and family lift](collatz_join_shields_and_lifts_20261004.md),
[the binary/ternary selector](collatz_binary_ternary_guard_fusion_20261004.md),
and the incoming rooted proof graph. Canonical hostiles: 3,7,27,703;
the negative cycles; and a cyclic collection of unsupported obligations.
Corrected near miss: increasing a decoder's length cannot repair an
invariant obstruction. Least-used sidecars: root-centered cofactor,
inverse branch excess, and a well-founded rank distinct from numeric order.

Anchor: a proper rank for proof moves. Niche: explain the branch delay and
its exact ternary coordinate. Wildcard: recover the old Berggren odometer
and interpret critical states as a graph-contraction problem. The live board
is **original source / branch toll / ternary guard / oriented carry /
root distance / proof rank / grounded component**.

The relevant recovered canon is
[THM-3756, odd-square ordinal Berggren descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md).
It supplies a model of descent with explicit domains; it does not supply
a Collatz conjugacy. The older [ternary Berggren note](ternary_berggren_20260925.md)
already proves the odometer isometry used in section3, and excludes
Berggren comparability for consecutive nondegenerate plus-Collatz edges.
The connection here retains that obstruction.

## 2. The entire sharp hierarchy behind q+4

Write `U_epsilon(n)=oddpart(3n+epsilon)` on positive odd integers, for
epsilon=+1 or -1. Suppose

    n=2^q*t-epsilon, q>=2, t positive odd,
    3^q*t-epsilon=2 mod4,
    J=(3^q*t-epsilon)/2.

The canonical reverse path from J has q exponents `(2,1,...,1)` and ends
at n. Its d-th vertex is `2^d*3^(q-d)*t-epsilon`.

**PROVED branch-toll theorem.** Suppose another reverse path first differs
at one of these q positions, using an exponent larger by2k, k>=1. If it
ends at a positive integer smaller than n, its depth is at least

    q+ell(k),   ell(k)=min{ell>=1:3^ell>2^(ell+2k)}.       (1)

This has no cap on subsequent inverse exponents. The sequence begins

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ell(k) | 4 | 7 | 11 | 14 | 18 | 21 | 24 | 28 |

Thus q+4 is the first member of an exact hierarchy, not a standalone
numerical coincidence. The first two coefficient gaps are17 and139.

**Proof, including the carry.** Let the first changed position be j<=q,
the total depth be q+ell, and the new endpoint be z. The inverse letter is
`E_a(x)=(2^a*x-epsilon)/3`. Increasing its exponent by2k changes its output
to `4^k*E_a(x)+epsilon*(4^k-1)/3`. Every subsequent inverse letter is at
least the exponent-one letter on positive inputs. Therefore

    z+epsilon >= C*(n+epsilon)
                   -(2*epsilon/3)*(4^k-1)*(2/3)^(q+ell-j),
    C=4^k*(2/3)^ell.                                  (2)

If C>=1, the minus-sheet correction is positive, so z>=n. On the plus
sheet its magnitude is strictly less than2C/3, giving

    z-n > (C-1)*(n+1)-2C/3 >= -2/3.

Since z-n is an integer, again z>=n. A smaller endpoint requires C<1,
which is exactly (1). The distance2/3 between the affine fixed points
-1 and -1/3 is the carry correction that makes this argument precise.

**PROVED sharpness for every q>=2,k>=1 on the plus sheet.** Let d=ell(k).
Choose n simultaneously in the reset-two dyadic class and the sibling
ternary class:

    n=2^q*(3 if q even else1)-1 mod2^(q+2),
    n=-(4^k+2)/(3*4^k) mod3^d.

CRT supplies infinitely many positive odd n. Set

    S(n)=4n+1,
    m=2^d*(S^k(n)+1)/3^d-1.

The inherited sharp sibling criterion gives 0<m<n. The actual words

    n: (1^(q-1),2),
    m: (1^d,1+2k,1^(q-2),2)

meet at J, with reverse depth q+d and first branch excess2k at position q.
The m path cannot pass through n: that would give n a first reset different
from its specified2, or meet a strictly larger proper prefix of its path.
Thus this is a genuine bypass and attains the bound.

For q=2,k=1 the least source in this construction is283 and the child223.
For q=3,k=2 they are4647 and4351. The earlier examples are successive
instances of the same sharp theorem. The saved audit constructs84 infinite
families for q=2..8,k=1..12; the proof covers all q,k.

## 3. The same toll has a ternary and a noncommutative meaning

Let `E(x)=(2x-1)/3` and `S(x)=4x+1`. Then

    E^ell(S^k(n))-S^k(E^ell(n))
      =2*(4^k-1)*(3^ell-2^ell)/3^(ell+1).              (3)

For k>=1, the ternary valuation of this defect is `v3(k)-ell`.
Consequently the two orders can both be integral only when3^ell divides k.
This condition is also sufficient for their integrality domains to agree:
the first order requires n=R(k) mod3^ell, and the second n=-1 mod3^ell,
where

    R(k)=-(4^k+2)/(3*4^k),
    v3(R(k)-R(0))=v3(k).

At k=ell=1 the defect is2/3: the two operation orders cannot both be legal
integer paths. This is an exact obstruction to a commuting square, not
an invocation of Eckmann-Hilton to erase the carries.

The recovered Berggren height coordinate is

    H_2(j)=2*(4^j-1)/3.

Extended to ternary arguments, the current address is exactly

    R(k)=-1-H_2(-k),    R(k+1)=(R(k)-1)/4.             (4)

Thus it is a reflected inverse odometer. One sibling step costs a factor4
in real size and translates its ternary address by one. The number ell(k)
is simultaneously the minimum inverse contraction needed to pay for that
branch and the ternary precision needed to realize its integer endpoint.
The finite-address isometry does not supply an arbitrarily prescribed
integer endpoint; its height and integrality guards remain essential.

## 4. An integral energy with a genuinely well-founded order

Now use the signed odd map `U(n)=oddpart(3n+1)` on nonzero odd integers.
For n!=1 write

    n-1=2^K*t,   K>=1, t odd and signed,
    E(n)=3^K*t^2,    R(n)=(E(n),K),
    R(1)=(0,0).                                      (5)

Order these pairs lexicographically. They lie in `N x N`, so strict
descent is well founded. Moreover

    E(n)^2 = 9^K*t^4 >= 8^K*abs(t)^3 = abs(n-1)^3.    (6)

The rank has finite sublevel sets, with an explicit integer cube-root
bound. It is defined from the input, not from an unknown stopping time.

**PROVED neutral energy and consumed precision.** On every actual
valuation-two edge,

    U(n)-1=3*(n-1)/4,
    (t,K) -> (3t,K-2),
    R(U(n))=(E(n),K-2)<R(n).                          (7)

The energy is the multiplicative quantity that stays fixed when two binary
digits are exchanged for a factor3 in the odd cofactor. Precision supplies
the strictly decreasing second coordinate. This resembles a conserved
quantity plus a consumed resource; it is an exact arithmetic identity.

The energy also allows increases in the integer:

    171 ->257 ->193 ->145 ->109 ->41 ->31,
    E: 21675,6561,6561,6561,6561,675,675.

At equal energies K strictly falls. The last source31 is a critical state.
No root certificate is inferred merely from this path. This repairs the
incoming171/257 circularity at the level of induction order:257 really is
lower in this independently defined rank.

## 5. Exact classification of every remaining local minimum

Neighbors mean the actual forward neighbor U(n) and **every** odd inverse
neighbor `(2^a*n-1)/3`, a>=1 when integral. They preserve the eventual basin
in either direction. We are orienting proof moves on this undirected graph,
not claiming that the original forward trajectory decreases in R.

**PROVED critical-set theorem.** A signed odd integer has no neighbor of
strictly smaller rank iff it is1, -1, or satisfies

    n=3 mod4,   n!=2 mod3,   n!=11 mod32.              (8)

Equivalently, apart from the two designated fixed points, the critical
classes modulo96 are

    3,7,15,19,27,31,39,51,55,63,67,79,87,91.

The positive critical domain has odd-relative natural density7/24.
It is a set of local minima for this specified rank, not a set of actual
counterexamples or the repository's best necessary domain for one.

**Proof: all inverse exponents.** Treat n=1 separately. For n!=1 and an
inverse neighbor m, the possible
exponents have fixed parity modulo2. At a=1 one has
`E(m)=(n-2)^2/3` and K(m)=1. At a=2 one has E(m)=E(n) and K(m)=K(n)+2,
so this never decreases rank. At every a>=3,

    K(m)=2,   E(m)=(2^(a-2)*n-1)^2>E(n)

for n!=1; the root case is immediate. Thus only the actual exponent-one
predecessor could lower rank. For n=3 mod4 this predecessor lowers rank
exactly when n=2 mod3, except for the self-edge at -1.

**Proof: the forward edge.** If n=1 mod8, (7) applies. If n=5 mod8,
K(n)=2 and the next valuation is at least3; the numerical contraction
implies strict energy decrease, on both signs. The remaining case is
n=3 mod4, so K(n)=1 and U(n)=(3n+1)/2. Put b=K(U(n)). Comparing squares
gives

    E(U(n))<E(n)
      iff 3^(b-1)*(3n-1)^2 < 4^b*(n-1)^2.            (9)

For positive n, b<=3 cannot work. For b>=4 it works except at n=11;
the b=4 comparison reduces to `13n^2-350n+229>0`, true for n>=27,
and every other positive member of11 mod32 is at least43. The exception11
has a smaller inverse-one neighbor7. For negative n, b>=4 always works;
b<=3 never gives strict rank descent (the equality cases -1 and -5 do
not lower the second coordinate). The b>=4 condition is exactly11 mod32.
Combining the two neighbor directions proves (8).

This second appearance of four uses `4^4>3^5`, with gap13. The branch
shield uses `3^4>2^6`, with gap17. They are different balance inequalities,
not the same theorem simply because both thresholds equal four.

**PROVED reduction.** Repeatedly choose a lower-rank neighbor whenever one
exists. Lexicographic well-foundedness forces a finite path to (8),1 or -1.
Each edge preserves basin and sign. Thus positive Collatz is equivalent to
grounding all positive critical states. No convergence assumption entered
the reduction. In particular the rank cannot certify an unsupported cycle.

## 6. A larger-child family connecting223 and233

The [existing exact common future](triplet_crt_223_233_20261004.md) has words

    223 --(1,1,1,1,3)-->425,
    233 --(2,1,1,1,2,3,1,1,2,1)-->425.

Its exact family lift is

    n_t=223+62208t,   m_t=233+65536t, t>=0.            (10)

The child is larger at every height: `m_t-n_t=10+3328t>0`. Its ordinary
lift slope is256/243>1, so the earlier nonexpanding-slope filter rejects it.
But K(n_t)=1 and K(m_t)=3, and

    3*(m_t-1)<4*(n_t-1)

at every height. Therefore **R(m_t)<R(n_t) for every t>=0**. A supplied
child certificate transports across the same actual common future.
This is a valid all-height induction step in a different well-founded
order, not a claim that all children in the family are already grounded.

More generally, refine a join's exact periods by a common power of2 until
K(n_t)=K_n and K(m_t)=K_m stay fixed. If its seed has strictly lower child
energy and its periods P,Q satisfy

    Q^2*3^K_m*4^K_n <= P^2*3^K_n*4^K_m,              (11)

then the entire family has lower child energy. Compare the positive linear
square roots of the energies: the initial value is smaller and the slope
is no greater. Equal seed energies are also allowed when K_m<K_n, with
the corresponding non-strict energy comparison. Exact guard refinement
is necessary; keeping only the seed's K values is insufficient.

**Scope of the gain.** The incoming child-cut operation can replace233 by
its successor175 and produce the broader ordinary-smaller family
`223+20736t <-175+16384t`. Thus (10) establishes a new permitted induction
order, not new coverage beyond every earlier family operation.

There is a general explanation. At a positive critical source N, K(N)=1.
If R(m)<R(N), perform valuation-two steps on m until its K is1 or2.
At K=1, energy comparison gives the resulting integer below N. At K=2,
one further actual step has valuation at least3 and gives an integer below N.
If m reaches1 earlier, stop there. This is a finite, explicit normalization,
not a request to run an unknown orbit until descent. The rank can recognize
a payable future obligation before that normalization has been expanded.

## 7. Signed folds and the old triplet

The square in (5) forgets an orientation:

    R(n)=R(2-n).

It pairs the positive critical states3,7,19 with the known negative-cycle
minima -1,-5,-17. This is a rank symmetry, **not** a conjugacy of U.
Indeed, if a is the actual valuation at n, the valuation at2-n is a when
a=1 or2; valuation3 is exchanged with a valuation at least4. At a=1,

    2-U(2-n)=U(n)-2,

whereas at a=2 the two sides without the final -2 are equal. The oriented
carry is precisely what the energy fold loses. The signed controls retain
the negative cycles rather than accidentally ruling them out.

On residues1 mod6 in modulus18 or36, this reflection has another exact
role: `h*(2-h)=1` in the residue ring. It therefore acts as group inversion
on those six-spaced multiplicative subgroups. In particular it fixes1 and
exchanges7 and13 modulo18. On the three odd classes modulo6 it fixes1
and exchanges3 with5. This is a concrete link to the earlier triplets;
it preserves that quotient's group law, not arbitrary Collatz trajectories.

### The branch gap and the negative-cycle denominator are the same invariant

The [boundary compiler](collatz_boundary_compiler_20261004.md) and the older
[coincidence atlas, section4](collatz_coincidence_atlas_20260927.md) already
identify the -17 cycle's seven odd steps, eleven halvings, and raw gap139.
The new branch theorem gives those counts another exact role. Put
`r=ell(k), A=ell(k)+2k`. The positive gap

    Delta_k=3^r-2^A

is both the strict margin that first allows the branch to be paid and the
raw denominator of the fixed point of every word with these counts.
For an ordered word w, write `F_w(n)=(3^r*n+B_w)/2^A`; its anchor is
`rho_w=-B_w/Delta_k`. The ordering determines B_w and hence whether the
anchor is an integer and which signed basin it represents.

**FINITE-EXACT, exhaustive at three specified pairs:**

| k | r,A | Delta_k | All positive compositions of A into r letters | Integral anchors |
|---:|---|---:|---:|---|
| 1 | 4,6 | 17 | 10 | -5,-7: the two rotations of the repeated word12 |
| 2 | 7,11 | 139 | 210 | The seven rotations of the -17 cycle |
| 3 | 11,17 | 46075=25*19*97 | 8008 | None |

The audit enumerates every ordered positive composition, computes its
exact rational anchor, and replays every integral result with actual
valuations. Thus the first two matches are exact, but the third barrier
does not introduce an integer cycle at those counts. This is not an
assertion excluding cycles at other lengths or costs.

There is an all-word explanation for another apparent prime-generating
pattern. Repeating a word doubles its counts and multiplies both its gap
and its carry by `3^r+2^A`; the reduced anchor stays unchanged. Repeating
the word12 gives `3^4-2^6=17` while still fixing -5. Repeating the -17 word
gives `3^14-2^22=139*5*7*11^2` while still fixing -17. New factors of the
raw gap need not represent new dynamics. The pair `(B_w,Delta_w)` has a
projective quotient, the anchor; retain the word as a sidecar because
that quotient discards the consumed binary fuel and the exact guards.
This supplies a concrete interpretation of recursive number appearances:
some encode new addresses, while others encode repeated presentations
of the same fixed point.

## 8. A more ambitious well-founded controller, with explicit obligations

The graph reduction suggests a different research target: **cancel the
remaining critical states by guarded common-future paths whose final
dependency has lower R**. The path may cross a large numerical or energy
barrier. The original source and its rank stay fixed until the whole path
is checked. Merely descending from a larger intermediate state is not enough.

A complementary way to represent a long path is a finite controller with
rational anchors. Here is a sufficient criterion that admits growing loops.

**PROVED conditional controller criterion.** A finite control graph has
exact guarded, nonempty forward-word edges. For every strongly connected component,
assign a rational anchor rho_s at each mode such that every internal edge
F_e sends rho_s exactly to rho_t. Denominators are odd. Require that no
nonterminal actual input equal its mode's anchor. Then an internal edge
of halving cost A_e consumes exactly A_e units of the nonnegative fuel

    K_s(n)=v2(n-rho_s).

Order components by a decreasing topological index h. The pair `(h,K_s)`
strictly decreases on every edge: internal edges decrease fuel, and exits
decrease h. Acyclic singleton modes need no anchor and may use fuel zero.
If the exact guards cover every admitted nonterminal input and all terminal
exits have supplied basin certificates, the controller terminates and
transports those certificates to every admitted input.

The proof is `F_e(n)-rho_t=(3^r/2^A_e)*(n-rho_s)` and finite-component
acyclicity. It uses no guessed stopping time. A simple cycle has the anchor
`rho=B/(2^A-3^r)`. For the growing word112 it is -19/11, and each repeat
consumes four binary digits of11n+19 while numerical size increases.
The one-letter growth loop uses -1. The -17 shadow supplies another chart.

This gives a concrete proposed carrier: original integer and rank;
current guarded word; controller mode; rational anchor and remaining
precision; alternative child identities; and grounded proof provenance.
Attach such a controller to a critical family and require every exit to
pay the **original** rank. The combined lexicographic tuple
`(E(original),K(original),component-index,anchor-fuel)` then decreases
through internal work and dependency replacement. Universal existence and
guard coverage of these controllers remain OPEN.

Two exact obstructions prevent this from being a formal restatement of
success. First, two internal cycles with different fixed points cannot use
the same single anchor at their common mode. Their branching component
must be refined or supplied with another proved ranking mechanism. Second,
an unweighted minimum of normalized anchor energies does not solve that
problem. For a word of length r,cost A, center b/d, let

    d*n-b=2^K*t,    H=3^(rK)*abs(t)^A.

H is conserved by the word and K drops by A. But with a common degree and
denominator, every negative fixed-point chart has energy greater than
the positive-root chart for n>1: its distance from n is larger and its
coefficient expands. The naive minimum collapses to the rank already
studied above. Phase-dependent guarded controllers, rather than an
unjustified minimum over incompatible charts, are the proposed repair.

A second cheap hostile rules out a broad finite-phase shortcut. For any
modulus M and positive t, `n=4Mt-1` has U(n)=6Mt-1; both have residue -1
modM while the integer grows. A rank `a_r*n+b_r` with every a_r positive
therefore cannot decrease on every odd edge when r is just n modM.
Adding fixed prime phases does not repair that rank class. An unbounded
anchor-precision coordinate is materially different from such a phase.

**An incoming infinite cancellation of critical states.** The concurrent
[partitioned-completion construction](collatz_partitioned_completion_20261004.md)
supplies the exact family

    n_t=27727075633746555+79062194724345216t,
    h_t=25270565367447551+72057594037927936t, t>=0,
    n_t --(1,2,1)-->J_t<--(1^32,2,19)-- h_t.

Every source is27 mod96 and hence belongs to our critical set. Both source
and child are3 mod4, so both K values equal1. The positive differences
between the constants and between the periods prove `0<h_t<n_t`, hence
`R(h_t)<R(n_t)` for the entire family. This is a concrete cancellation
rule inside the residual problem. It still requires a child certificate
or an induction covering every possible dependency. Our audit preserves
the incoming periods; it makes no new claim about their maximality.

The incoming whole-cell theorem says that, for any fixed forward/inverse
depths R,S, the entire cell `n=-1 mod2^(R+1), n=0 mod3^S` has no smaller
common-future child within those bounds. It targets ordinary integer
descent. Permitting a larger ranked child does not contradict it: the
explicit subsequent normalization has a length controlled by the child's
unbounded K. Nor does the theorem exclude an unbounded-word controller.
The two results explain why retaining precision, rather than merely
enlarging a fixed-depth table, is a substantive change of representation.

This connects to established termination research without importing a
solution. Yolcu, Aaronson and Heule construct mixed binary/ternary rewriting
systems equivalent to Collatz termination and use matrix interpretations
to prove restricted variants. Their result does not prove Collatz; our
anchor criterion and graph-rank reduction above are independently proved
scoped constructions. [Primary paper, sections2--4](https://emreyolcu.com/research/rewriting-collatz.pdf).

## 9. Finite experiments and what they actually improve

The signed branch audit checks125 reset-two sources between3 and501,
every smaller positive odd candidate, and all forward depths through
q+ell(6), on both signs. Sixteen noncanonical smaller hits obey the general
delay law. The84 constructed sharp families are replayed at parameters
0,1,17 and one million. The commutator is checked independently by rational
composition, including failed integrality controls.

The rank classification is checked against **all** lower-rank inverse
neighbors for2002 signed odd inputs in[-2001,2001]. The properness bound
makes that search complete without an exponent cap. All10002 signed odd
inputs in[-10001,10001] reduce to the stated critical set. This finite
audit includes1351 magnitude-increasing edges; it does not assert that
each reduction ends at one of the known cycles.

On the frozen239 positive sources, the following search keeps six forward
steps, permits a lower-ranked child, and requires the all-height test (11):

| Reverse depth | Ranked-family hits | Numerically larger children |
|---:|---:|---:|
| 2 | 3 | 3 |
| 4 | 3 | 3 |
| 6 | 3 | 3 |
| 8 | 10 | 3 |

At reverse depth8, all seven inherited ordinary-rank hits survive. The
three additional sources are303,3067,5191. Their larger children1025,5825,
9857 normalize explicitly to61,1229,4159 using words22224,223,222.
They were already within the earlier eight-forward-step comparison. This
is earlier symbolic recognition at a fixed observation window, not new
global coverage or an improvement over every existing macro bank.

The inverse search uses (6) to bound all possible lower-ranked children,
then the inherited complete inverse depth bound. It is independently
checked by forward enumeration on2880 small source/endpoint/depth/cost
combinations. Failed candidates remain unproved; no child home certificate
is invented by the family search.

Each saved ranked family also includes all three exact coefficients of
`E(n_t)-E(m_t)`. They are nonnegative, with a strict constant term or the
specified strictly smaller K in the equality case. This is a directly
checkable polynomial certificate for every parameter, beyond sampled replays.
An independent symbolic replay checks each exact valuation on a whole
affine family: for a word letter a, its numerator constant is2^a mod2^(a+1)
and its numerator period is0 mod2^(a+1). It then divides both coefficients
by2^a. This certifies the common future for every parameter, including the
incoming critical family, without relying on a sample of heights.

## 10. Connection contracts and next decisive tests

| Connection | Exact map or retained object | Preserved predicate | Loss / next test |
|---|---|---|---|
| Branch toll to ternary precision | k -> ell(k), with exact CRT guard | Smaller integer bypass | Carries and feasibility; all-k proof and signed controls |
| Berggren height to sibling address | R(k)=-1-H_2(-k) | Ternary distance and adding-machine action | Ordinary integer height; no PPT/Collatz conjugacy |
| Root distance to proof rank | n -> (3^K*t^2,K) | Well-foundedness and basin along checked edges | Orientation and unresolved critical classes |
| Graph contraction to critical problem | Repeated lower-rank actual neighbors | Connected component / eventual basin | Original directed chronology; critical cancellations still owed |
| Mod18 triplet to rank reflection | n ->2-n | Energy and subgroup inversion | Oriented carry; explicit failure to commute at valuation one |
| Growing word to finite controller fuel | n ->v2(n-rho_s) | Termination while guards and anchors transport | Guard coverage and incompatible cycle anchors |

The next mathematical test is to construct controllers on the remaining
critical classes with complete guards and source-relative rank payment.
Test3,7,27,703 first; then test the signed cycles and the simultaneous
presence of incompatible rational anchors. A finite table of successful
controllers is useful only in its actual covered domain. The long-range
goal is universal critical cancellation, not a larger claim from the
already-complete239-source experiment.

## Reproduction

Run `python3 -B 04-computation/experiments/collatz_branch_toll_rank_20261004.py`
and repeat with `python3 -O -B`. The
[script](../../04-computation/experiments/collatz_branch_toll_rank_20261004.py),
[JSON](collatz_branch_toll_rank_20261004.json), and
[stdout](collatz_branch_toll_rank_20261004.out) contain the universes,
exact families, source-preserving search and explicit checks. Normal and
optimized runs are compared byte-for-byte. The global theorems above are
proved algebraically; the finite audits do not replace their quantifiers.
