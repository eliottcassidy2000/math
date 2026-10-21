# Sharp arrival-order capacity and native temporal block completion

2026-10-20. **PROVED:** the exact twelve-cut subset envelope, sharp source-norm
capacity 27 for blocks safe under every arrival order, primitive-circuit
prefix safety, and source-owned interval completion. **FINITE-EXACT:** the
retained native N=4096 packets and completion bank, and the N=7 failure.
**OPEN:** a uniform 153N interaction bound and universal Collatz entry.
There is no optimal packing or historical-priority claim.

Files: [checker](check.py), [complete deterministic receipt](receipt.json).

## 1. Inheritance, original sources and the missing coordinate

The exact ledger is inherited from
[native annuli](../session18_native_annulus/README.md): for primitive positive
sources s in [4N,256N), retain the enlarged charge A_(64N)(s), actual first
arrival at the fixed 627-point frontier V, full canonical-P braid clock t_s,
and depth d_s. These are not literal fully odd U-word clocks. Write

    v_s=(-1)^t_s A_(64N)(s),  Q(a,b)=a²-ab+b²,
    beta(v,w)=Q(v+w)-Q(v)-Q(w).

Compatible source registers s mod15 and s mod35 jointly encode mod105;
they are not 525 independent states. The checker authenticates the original
source, both registers, charge, terminal and full clock before using a row.

At each finite horizon d, the established identities are

    Delta E=A_d+I_d,    A_d+M_d=153N,
    Delta E<=153N iff I_d<=M_d.                         (1)

The closest mechanisms are [chronological opposite matching and its six-source
hostile](../owned_pair_blocks_20261020i/README.md), [primitive positive
circuits](../eisenstein_circuit_completion_20261018g.md), and [native sibling
triples](../sibling_annulus_guards_20261020i/README.md). The correction to
timeless cancellation is that unfinished members still contribute. Three
individually safe opposite pairs can simultaneously expose three equal
vectors: their energy is 9Q(v), above their total norm 6Q(v).

Board: original source / signed alphabet / actual arrival interval / quadratic
coherence / unique pair address / unspent moving norm. The anchor is a
stronger residual guard. The niche is an exact finite-alphabet envelope.
The wildcard is larger blocks whose temporal support, rather than mere final
cancellation, proves their payment. No sunflower theorem is needed here.

## 2. Exactly twelve cuts replace all arrived subsets

The signed nonzero alphabet, in angular order, is

    (1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),
    (-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1).

It alternates norm1 and norm3. For physical integer multiplicities c_i>=0,
put T=sum c_i Q(v_i), and

    H_i=sum_{j=i}^{i+5} c_(j mod12) v_(j mod12),
    K=max_i Q(H_i).                                    (2)

**Exact envelope theorem.** K is the maximum of Q over all submultisets of
the packet. Consequently every possible arrival order obeys packet
energy<=T at every horizon if and only if K<=T.

Proof. Choose a maximizing subset with sum w. Removing any included copy v
gives beta(w,v)>=Q(v)>0. Adding any excluded copy gives
beta(w,v)<=-Q(v)<0. Thus no type is partly included: it cannot satisfy both
inequalities. The chosen types are precisely the populated directions in
the open half-plane beta(w,v)>0. No populated direction lies on its
boundary. Perturb the boundary through absent directions if necessary;
one of the twelve consecutive six-direction cuts gives exactly the same
sum. Every such cut is itself an admissible subset. The empty packet has
K=T=0. Opposite directions and zero-count boundary ties cause no exception.

For the native application, all packet members must have the same actual
terminal. Its original norm T is assigned once from 153N. A prefix with
arrived norm A_B and vector H_B obeys Q(H_B)<=T precisely when its internal
interaction Q(H_B)-A_B is at most its own moving norm T-A_B. This spends no
future negative interaction; it uses already owned source norm. It does
not authorize removing a partially arrived zero packet.

## 3. Sharp capacity 27; all primitive circuits are safe

**Sharp bound.** For every such packet,

    K>=T²/27.                                         (3)

Therefore a nonempty block safe under every arrival order has T<=27.
No zero-sum or circuit-free premise is used.

Proof. Regard (2) as a convex function of real nonnegative counts: it is a
maximum of positive semidefinite quadratic forms. The integral rotation
R(a,b)=(a-b,a) permutes the cuts and the alphabet. Average the count vector
over its six rotations. This cannot increase K and preserves T. In the
average each unit direction has multiplicity u, each norm3 direction
multiplicity w. This fractional average is only a convexity proof device;
the source packet still consists of physical integer members. A cut sums
to (2u+2w,2u+4w), so

    K>=4u²+12uw+12w²,    T=6u+18w,
    27K-T²>=72u²+108uw>=0.

The constant and the safe upper endpoint are both attained. On the six
norm3 directions

    (2,1),(1,2),(-1,1),(-2,-1),(-1,-2),(1,-1)

take multiplicities (2,1,2,1,2,1). The two alternating equilateral triples
show total vector zero. Every cut has norm27, so this nine-member packet
has T=K=27. Merely repeating any safe block need not preserve safety.

The inherited six opposite and forty primitive three-type circuits are
all safe under every internal arrival order. For a three-type circuit
sum c_i v_i=0, convexity on the rectangular count box shows a maximum at a
corner. Nonzero corner sums are exactly +/-c_i v_i. Therefore

    K=max_i c_i² Q(v_i).                               (4)

The exhaustive determinant classifier gives the following exact table.
The same coefficient shape can have different direction norms.

| Sorted multiplicities | Number | T | K |
|---|---:|---:|---:|
| 1,1,1 | 2 | 3 | 1 |
| 1,1,1 | 6 | 5 | 3 |
| 1,1,1 | 2 | 9 | 3 |
| 1,1,2 | 12 | 6 | 4 |
| 1,1,3 | 6 | 9 | 9 |
| 1,2,3 | 12 | 12 | 12 |

For an opposite pair T=2q,K=q. Repeating a fixed circuit k times gives
T_k=kT and K_k=k²K, so the exact safe repetition bound is k<=floor(T/K).
These statements include four-, five- and six-member circuits, without
assuming that repeated types are simultaneous physical arrivals.

## 4. Actual clocks permit a much larger native block

The nine-member sharp template is realized at N=4096, terminal1, by sources

    92051,96971,103555,72729,92049,276145,107623,122731,77667.

The receipt retains each full source row, not just these residues. Its
all-order envelope is27; its actual prefix energy is at most3.

There is also a much larger clock-dependent block at terminal711. Start
with the inherited canonical rank-matching of the +(2,1) and -(2,1)
classes. Discard synchronous pairs for this control. Sort exposure
intervals [min(d_x,d_y),max(d_x,d_y)) by finishing depth and keep disjoint
intervals. In the fixed N4096 bank this selects 56 pairs, or112 sources.

    total original norm T=336;
    maximum actual prefix energy=3;
    maximum over arbitrary arrival orders K=9408;
    final arrival depth=82; final packet vector=0.

Proof of the actual bound: at any horizon at most one pair is unfinished;
all other selected pairs are either absent or have summed to zero. The
remaining vector is zero or +/- (2,1). This is a reusable certificate for
arbitrarily many disjoint exposure intervals when such actual intervals
are supplied. It is not a claim that every native terminal supplies them.
The finite bank is a selection from the already authenticated source
universe, not a new convergence census. Different controls may reuse a
source; they are not presented as one simultaneously spendable partition.

The distinction from (3) is essential: this block needs its original
clock order. Permuting arrivals can expose56 equal vectors and attain
9408. This is a concrete reason to enlarge blocks through timing rather
than through unconditional cancellation alone.

**Actual complete-interaction payment.** Insert this one block last,
after every other old and new source. At terminal711 let E_d be the sum
of every arrived source outside the block, with its original sign and
clock. The full insertion bill is

    delta_B(d)=Q(B_d)+beta(E_d,B_d).

For all finite horizons in the complete393,216-source population,
delta_B<=156<336. The maximum156 occurs exactly at depth33, where
B=(2,1), E=(51,39), and beta(E,B)=153. The minimum is-147 at depth62.
The receipt retains every horizon and the hash of the complete source
population. Beyond the last arrival these rows are constant.

This has a reusable exact directional guard. For k disjoint-exposure
opposite pairs of vector v, original budget2kQ(v), an active block has
B_d=sigma*v with sigma in{+1,-1}. Its full insertion pays if and only if

    sigma*beta(E_d,v)<=(2k-1)Q(v).                      (5a)

When B=0 it pays automatically. For this block v=(2,1), beta(E,v)=3a
for E=(a,b), so the exact guard is sigma*a<=111. The actual maximum
sigma*a is51. The orthogonal coordinate b is unrestricted by this
directional guard; no blanket bound on Q(E) is needed.

The general ordered-insertion theorem is equally explicit. Partition
all new source addresses into blocks, assign each its original norm T_j,
and fix one insertion order. At a horizon, E_j is the old vector plus
the earlier blocks at the same terminal. Then

    delta_j=Q(B_j)+beta(E_j,B_j),
    sum_j delta_j=Delta E,      sum_j T_j=153N.         (5b)

Thus delta_j<=T_j for every j and every horizon suffices for the exact
153N bound. This follows by telescoping Q(E_j+B_j)-Q(E_j). The present
control authenticates one block in last position; it does not establish
the bills for all preceding blocks. Several independently successful
"last block" tests cannot all claim the same last-position credit.

**Infinite source-preserving subpacket lifts.** The56 pairs split into24
with both sources on branch a=1,10 with both on branch a=2, and22 mixed
pairs. Retain either same-branch class. With c_1=35,c_2=49, set

    N_k=64^k*4096,
    s_k=L_a^k(s)=64^k*s+c_a*(64^k-1)/63,   k>=0.

These are the inherited [native phase lifts](../session18_native_phase_lifts/README.md),
now restricted to an annulus packet with one common branch. Every member
is outside V; all its positive lifts exceed max(V)=63665. The exact
parent is unchanged, the depth and terminal are unchanged, and its full
clock increases by6k. No convergence assertion is needed beyond the
supplied original routes.

The band is also unchanged. Its integer endpoints are even; an original
odd source is at least one integer inside. Multiplication by64 and
addition of c_a<64 therefore preserve the corresponding band of the
scaled annulus. The band-charge formula uses the same sequence of lower
and upper palettes, while s mod3 changes by c_a. As both palettes obey
R(v_r)=v_(r+1), with R(a,b)=(-b,a-b),

    A_(64N_k)(s_k)=R^(c_a*k) A_(64N_0)(s).

The even clock increment preserves the sign, so every signed vector in
the subpacket undergoes the same norm-preserving rotation. All exposure
intervals and all prefix energies persist for every k. Thus there are
two explicit infinite families of actual native zero packets:

| Branch | Members | Original norm | Actual maximum | All-order maximum |
|---|---:|---:|---:|---:|
| 1 | 48 | 144 | 3 | 1728 |
| 2 | 20 | 60 | 3 | 300 |

For arbitrary modulus m>=1 the exact symbolic source reader is

    (((63s+c_a)*64^k-c_a) mod(63m))/63.                (5c)

The numerator representative is divisible by63. Keeping modulus63m
before division is necessary, including when m shares a factor with63.
The joint mod105 rows have period dividing6: modulo105,
L_a²(s)=s+65c_a, and3*65c_a=0. Their mod15 and mod35 projections remain
compatible. The receipt gives all six rows and independently checks
literal k=0,...,6; the all-k conclusion comes from the identities above.
External populations change with k. These internal safety families do
not inherit the complete-interaction bound156: each later insertion
still needs the directional premise (5a) for its rotated vector.

Even at k=0 the two subpackets fail their own complete insertion budgets.
Branch1 costs156 at depth33 against144. Branch2 costs135 at depth61
against60, where B=(-2,-1) and E=(-44,-85). These values follow directly
from receipt `full_native_control.last_insertion_rows[33]` and `[61]`,
and the exposure intervals under `same_branch_lifts`.
All56 intervals are mutually disjoint, so when one component is active
every other component has vector zero; its external vector is exactly
the whole packet's stored E. Pooling all112 original source budgets
therefore succeeds where these two separately budgeted components fail.
No energy or source norm has been copied in that pooling.

## 5. Completing unpaid interaction intervals, with no reused credit

The inherited N4096 sibling bank contains21,152 disjoint source triples of
pairs. Each pair is (x,2x+1), x=1 mod8; actual frontier guards force the
same terminal/depth and a clock difference1. Its interaction is exactly
-beta(A_(64N)(x),A_(64N)(2x+1)), born at that actual depth. A consecutive
three-pair block has weights -3,0,+3, or all zero. Hence its full selected
debit is either zero, -3 on a half-open depth interval, or +3 on one.
The sign counts are8906 negative intervals,3394 identically zero blocks,
and8852 positive intervals. The zero class includes simultaneous -3/+3.

**Temporal completion theorem.** Let one negative block have interval
[a,b). Select positive intervals [c_i,e_i) contained in [a,b), mutually
disjoint as half-open intervals. Their combined selected debit is

    -3*1_[a,b)(d)+3*sum_i 1_[c_i,e_i)(d)<=0             (5)

at every finite horizon. This can cover many successive positive
intervals using one negative interval, because its amplitude is never
spent twice at the same horizon. Overlapping positive intervals fail
the guard; e.g. [2,4) and[3,4) against[1,5) leave debit+3 at depth3.

The deterministic bounded experiment processes negative containers by
decreasing length and greedily takes the earliest-ending available
disjoint positive intervals. It completes6088 positive blocks using2927
negative containers, leaving2764 positive blocks. No global matching
optimum is claimed. One retained negative interval[9,110) covers35
positive intervals:36 native triples and216 distinct source members.
The complete certificate gives every block address; the representative
also stores every source row and compatible mod15/mod35 registers.
Reproduction reauthenticates all selected pair terms from their original
sources and actual P routes, and recomputes the full bank hash.

This improves a precise residual upper bound. Write I_d=R_d+J_d, where
J is exactly the selected sibling bank and R retains every other pair
address, including all inter-block interactions. Dropping the separate
nonpositive blocks formerly gave

    I_d <= R_d + all_positive_interval_debit(d).

Drop the new nonpositive groups from (5), and separately the unused
negative blocks. The stronger result is

    I_d <= R_d + unmatched_positive_interval_debit(d). (6)

The improvement is exactly the completed positive debit, pointwise
nonnegative and reaching7116 at depth13 in this fixed bank. The remaining
positive debit reaches7677 at depth18. This does not alter the
physical I_d. A matched negative container is no longer independently
available; counting it again would spend the same pair terms twice.
Selected pair addresses, not just source names, govern this operation.
No source norm is withdrawn in (5) or (6).

More explicitly, R_d is the sum of all original old/new and unordered
new/new interactions not among the selected sibling pair addresses,
with their actual same-terminal indicator and later-arrival birth depth.
It includes every inter-block term. If U_all(d) and U_left(d) count the
active positive intervals before and after completion, the two bounds
are R_d+3U_all(d) and R_d+3U_left(d). The finite full-source replay checks
both against I_d at every horizon. At N4096 their excess over M_d is at
most-201674, attained after all arrivals; the final I=-201674,
new/new energy433053, old/new cross term-8039, and original new norm626688
agree with the inherited ledger. No new global N census was performed.

The completion theorem and the implication

    R_d+3U_left(d)<=M_d  ==>  Delta E<=153N             (6a)

hold for every N>=1 and every supplied source-authenticated disjoint
family satisfying (5). Existence of sufficiently many such families,
or a uniform bound on R_d, is not a premise established by this experiment.

## 6. The retained global obligation

For a source-disjoint partition into blocks with vectors H_B(d), the new
energy contains Q(H_B) for each block, but also every beta(H_B,H_C) and
the old-core terms beta(W,sum H_B). Thus standalone block safety is not
global safety. Once the original norms T_B are assigned once, the exact
remaining requirement is

    sum_{B<C} beta(H_B,H_C)+beta(W,sum_B H_B)
       <= sum_B [T_B-Q(H_B)].                          (7)

For multiple terminals sum only within each terminal. Residual sources
can be retained as singleton blocks, so (7) is an exact criterion for
the full partition, not an omitted error term. Equivalently, (6) pays
the unit-rate annulus increment if its explicitly retained right side
is at most M_d. Neither premise is proved uniformly here.

The N7,d8 hostile remains literal: arriving new norm1054, moving norm17,
new/new interaction32, old/new interaction-2. Thus I=30>17, deficit13.
Nonpositive old/new cross terms alone do not suffice. This small unit-rate
failure does not refute a differently scoped large-N or larger-rate
conjecture. The inherited N4096 global control also remains unchanged.

[Synchronous three-phase lifts](../session18_native_phase_lifts/README.md)
give another valid block mechanism. They compose with this one only
through retained source or pair-address ownership and the mixed old/new
baseline. A mixed packet's old norm cannot be copied into the new153N
budget. All statements here are at Nyquist clock signs, not all frequencies.

## 7. Validation and exact next target

The checker compares all4096 binary count vectors and all3^12 subset
choices against the twelve-cut formula; independently enumerates every
primitive determinant circuit and its full count box; checks the sharp
nine-member packet; reauthenticates the two native source packets and
all matched sibling terms; and rejects wrong full clocks, duplicate
sources, Boolean/float impostors and overlapping credit claims. Receipt
comparison is LF-normalized. Ordinary and optimized commands are:

```text
python -B 05-knowledge/results/annulus_block_completion_20261020j/check.py --verify
python -B -O 05-knowledge/results/annulus_block_completion_20261020j/check.py --verify
```

The reproducible body has1,615,263 checks; verification adds one saved
receipt comparison. LF SHA256:
`b85cfd8a95962c827d9e21d1418a8fe7f4086c6cca6203f042f547d32ad3b609`.

The next obligation is now explicit: control R_d plus the unmatched
positive intervals in (6), or exhibit a source-owned partition satisfying
(7). Enlarging an all-order-safe packet beyond norm27 is impossible;
the112-source control demonstrates the correct temporal escape from
that obstruction. These are guarded certificate improvements, not a
universal annulus or Collatz proof.

Independent peer proof/API/receipt audit: **PASS**. The reviewer rebuilt
all257 complete-insertion rows from every original source using the
independent charge engine, checked both subpacket budget failures, and
verified additional lift indices0,1,2,7,19 and moduli sharing factors
with63. Ordinary/optimized saved comparison also passed. Files frozen.
