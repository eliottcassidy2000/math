# Nontransitive centers give exact repair windows

A stronger replacement for the transitive mirror center: a seven-vertex
center has a width-four window in every matched transitive frame. More
importantly, a six-vertex width-three window survives **independent**
left and right addresses. Thus its baseline has two free quantities, not
the previous single mirror square.

## Interface and baseline

Mark a subset of a tournament B; the added apex beats exactly its marks.
Write h=H(B), n=H(B+apex), and q for the number of ordered covers of B by
two paths, allowing empty slots. Let s and t count weighted nonempty
paths starting marked or ending unmarked, respectively, with weight
H(complement). Let e,f be their signed weighted counts, sign positive
for odd path order and negative for even order.

Place a transitive word on each side of B and orient all cross arcs left
to right before adding the common apex. Encode the left marks by binary x
in the transitive order (MSB first); encode the right marks by reversing
and complementing the bits of positive y. Put

    p(x) = 2 sum_{i<j} 2^max(0,j-i-2) w_i(1-w_j).

Then the exact count of the whole tournament is

    qxy + t x + s y + f [x odd] + e [y odd] + h(p(x)+p(y)) + n.   (1)

For a direct explanation, use the positive odd-cycle packing formula for
Hamilton path counts (as in the inherited CUT_STATE.md). Contributions
whose apex cycle uses only one block, along with the no-apex
contributions, sum to n+h(p(x)+p(y)). Crosses of the left/center seam
supply x t+[x odd] f; crosses of the center/right seam supply s y+e[y odd].
For a crossing apex cycle meeting both exterior blocks, the middle has
an optional path: split its contributions by the parity of that path.
The signed sum of the optional path counts is zero because B is
nonempty, so allowing either parity produces the cross term qxy.
This is the weighted ordinal transfer, with ordinary counts of both
exterior blocks equal to one. Choosing x,y>=1 makes the graph strong:
the apex beats the first vertex of the left exterior and loses to the
last vertex of the right exterior; these give routes to and from all
vertices.

A full collision of (h,q,s,t,e,f) lets n vary in (1) without affecting
anything else. A mirror collision needs only (h,q,s+t,e+f) when x=y.

## Two literal banks

An orientation code uses lexicographic unordered pairs (i,j), i<j, and
bit 1 means i beats j. Mark masks have bit i set exactly when the apex
beats vertex i (opposite binary display order from the exterior address).

* At center order 6, take (code,mask)=(8850,54),(8881,54),(8881,53).
  They all have (h,q,s,t,e,f)=(45,432,258,129,30,15), with n=135,137,139.
  Thus the row offsets are exactly 0,2,4 for **every pair** x,y>=1.
  At x=y=1 the counts are 999,1001,1003; x=2,y=3 gives
  3879,3881,3883, and exchanging the addresses gives 3735,3737,3739.
* At center order 7, take code 528449 and masks 99,66,93,94.
  The common mirror packet is (h,q,s+t,e+f)=(27,564,537,27);
  n=151,153,155,157. The exact quadruple is therefore
  564x²+537x+27[x odd]+54p(x)+151+{0,2,4,6}. At x=1 this gives
  1279,1281,1283,1285. The normalized unit steps persist at every
  matched transitive exterior. The quartet is **not** four consecutive odds
  when x≠y: the first two variants have (s,t,e,f)=(360,177,20,7)
  and (174,363,10,17), so their difference is
  2−186(y−x)−10([y odd]−[x odd]), of magnitude at least 174 if x≠y.

The candidate banks were located by grouping all center/mask packets
through order eight by their full or mirror interface. The literal
witnesses and formulas, not completeness of the candidate catalogue,
are the basis of the theorem. Compared to the older three-odd mirror,
the first new bank removes the enforced equality x=y and the second
enlarges the exact window. Neither places that window at every target.

## Exact inversion within the first bank

To invert (1), choose which of x,y is smaller. Their product is at most
N/432 because all the terms are nonnegative and p>=0. Thus the smaller
address ranges from 1 to floor(sqrt(N/432)). Fixing this address and
j in {0,1,2} leaves a single nonnegative bit-function equation

    A u + 45 p(u) + b [u odd] = R,

where A>0 and b>=0. Moreover j is determined before the search:
j=2N modulo 3 because every other term is a multiple of 3. Modulo 9
one also obtains (N−135−2j)/3 ≡ 2(floor(x/2)−floor(y/2)) modulo 3. Unravel the low bit using
p(2k)=p(k)+k+[k odd] and p(2k+1)=p(k).
The two reduced states (slope,bias,target) are
(2A+45,45,R) and (2A,0,R−A−b), respectively. Both double the slope;
a positive u is impossible if R<A, and u=0 occurs exactly when R=0.
This gives an exact target-first solver for the bank with no approximate
rounding assumption. The solver is not a proof that some branch must
exist. Replay finds only 378 represented targets up to 20,000.

## Replay

    python selector.py
    python replay.py HOST_CODES
    python full_trial.py HOST_CODES

HOST_CODES is the two-column order/orientation-code catalogue used in
the pair-line work, or any desired center list. selector.py verifies
the literal examples and cross-checks inversion versus brute family
generation. The mirror and full catalogues share the same definitions
but intentionally use different keys.

## Why swapping two ordinal blocks is not a full-port repair

For nonempty marked A,B, in their ordinal concatenation AB the transfer
gives

    s_AB = q_B s_A + h_A s_B;       t_AB = h_B t_A + q_A t_B;
    e_AB = h_A e_B;                 f_AB = h_B f_A;
    n_AB = h_B n_A + h_A n_B - h_A h_B + s_A t_B + e_A f_B.

The signed coordinates reset to the terminal blocks because the
nonempty optional path parity sum in the other block is zero.
Put d_i=q_i−h_i>0. If AB and BA have equal full interfaces, equality of
s and t forces (s_A/d_A,t_A/d_A)=(s_B/d_B,t_B/d_B), and equality of e
and f forces (e_A/h_A,f_A/h_A)=(e_B/h_B,f_B/h_B). Hence both cross
products commute: s_A t_B=s_B t_A and e_A f_B=e_B f_A. In particular
n_AB=n_BA. There is no hidden unit intercept from merely commuting two
blocks while preserving the state seen by **independent** exteriors.

This does not rule out matched mirrors, where only sums of state
coordinates are compared, nor a compensating change of masks or a
longer word permutation that is not a single locally interface-
preserving exchange.
