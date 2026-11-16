# Three physical states make a canonical radix-four carry

Use the transitive two-block seam formula (rederived in the replay):
a forward right mask a and reverse left mask b give
H=M·P+C, where M=2^t,
P=(2^k−b)(2^r+a), and
C=1+rev_k(b)−rev_r(a)−J(a)−J(b).
The interaction J sums K(d)=2^(d−1)+[d=1] over marked pairs
at distance d. These are literal path counts on a transitive tape
with the marked head and a sentinel. For clarity the construction
is: all tape arcs point forward, all tape vertices beat the
sentinel, the sentinel beats the head, and the head loses to
precisely the marked tape vertices. First and last tape vertices
are unmarked, so every tournament below is strong.

For the seam identity, first group paths by predecessor and successor
of the head. A backward jump of length d has 2^(d−1) ways of
assigning vertices between the two runs; accounting for forward
adjacent jumps changes this to K(d) and gives total
1+Σ K(j−i), over unmarked i preceding marked j (including
the marked sentinel). Expand the departures from the empty-mask
baseline 1+2^(k+t+r). The linear contribution is the difference
of the two reversals plus 2^t(2^k a−2^r b). The interaction of
the left and right masks is 2^t ab, separately from J(a)+J(b),
since no across-seam pair is adjacent. Collecting the growing
terms gives the product P above.

**General three-state theorem.** Fix r≥5, D=2^r, a positive multiple
a of 4 with a<D/2−1, b=D−a, and let p=2^v be the lowest set bit of
a. Take widths k=r+1, r and spacer t=2, so M=4. Let m=D+a and
C denote the remainder at (a,b). Change the actual head contacts
to run through the three rows:

| Step | Right mask | Left mask | High product P | Remainder |
|--:|--:|--:|--:|--:|
| 0 | a+1 | b+1 | m²−1 | C |
| 1 | a−1 | b−1 | m²−1 | C+2 |
| 2 | a | b | m² | C |

These give consecutive odd counts. The first edge changes only the
low readout while keeping the literal high product; the second resets
the low readout and advances the literal high product by one. Nothing
has been added to the arithmetic readout by convention. In particular,
if a=3·2^v with 2≤v≤r−3, then C=1, and the three exact counts are

$$
4(m^2-1)+1,\quad 4(m^2-1)+3,\quad 4m^2+1.
$$

This is an unshifted canonical radix-four rollover. The first
instance r=5,v=2 uses masks (13,21)→(11,19)→(12,20) and counts
7741→7743→7745. The spine and sentinel do not change. The flips
affect v tape marks in each block on the first step.
On the second step v+1 arcs in each block change.
The simple choice v=2 gives fixed costs four and six at every height.

## Proof of the full mask table

The factor identity follows by using 2^k=2D: the two factors at
( a±1, b±1 ) are m∓1 and m±1. The low bits of a,b vanish
below v and their bit v is set; above v their binary strings
are complementary through bit r−1. Those upper marks have aggregate
binary weight D−2p.

Adding 1 introduces only bit zero in each block. The total interaction
increase is (a+b)/2=D/2, with no adjacency bonus because v≥2.
The reversal difference increases by D−D/2=D/2 as well. Hence
the remainder is unchanged.

Subtracting 1 replaces the two isolated lowest bits v by runs
of v ones at positions 0,…,v−1. In each block the reversal
increases by the same geometric sum, with twice the coefficient
on the left. Thus rev_k(b)−rev_r(a) increases by
D−3D/(2p). The new internal-run energies total 2(p−2).
For each upper bit j>v, replacing its old partner bit v by the
run increases its energy by 2^j(1−3/(2p)), minus one if j=v+1
(the lost adjacency bonus). Exactly one of the complementary
upper masks has that bit, because a<D/2 and b>D/2. Hence the
total J increase is
`2(p−2)+(D−2p)(1−3/(2p))−1 = D−3D/(2p)−2`.
The remainder increases by exactly two.

To see C=1 for a=3p, write S=D/(4p). Then a has exactly bits
v,v+1, b has bit v and all bits v+2 through r−1, and
rev_r(a)=3S, rev_(r+1)(b)=6S−2, J(a)=2, J(b)=3S−4.
The last sum is (S−2) inside the upper run plus (2S−2) from
its contacts with bit v. Substitution gives C=1.

## How often these are true radix resets

If r≥v+3, then C≡1 mod4 precisely when bit v+2 of a is zero.
Only distances 1,2 contribute to J mod4, each by 2. On the
upper interval of L=r−v−1 bits, the paired masks are opposite;
the indicator of being in the same mask is congruent mod2
to 1+α_i+α_j. The graph on those upper bits with edges at distance
one or two has odd-degree positions the second from each end,
with multiplicity (which cancels for L=3; for L=2 these are
the two endpoints). The shared low bit v contributes one
contact at each allowed distance. Thus the combined pair count
mod2 is 1+α_(v+2)+α_(r−2). The reversal contribution mod4 is
2−2α_(r−2), giving C≡1−2α_(v+2).

In particular fix v=2 and use a≡4 or 12 mod32 with 0<a<D/2
for r≥6. At least D/32 distinct canonical rollovers occur at
each such r, using the same constant edit budgets. Distinctness
follows from the leading term: a changes by at least 8 (or more),
so 4(D+a)² changes by at least 64D, while the digital
remainders have range width <5D. This yields Θ(√X) distinct
rollover windows up to X in this fixed-v family (counts are
bounded above and below by constant multiples of D², and
sum the geometric number of choices over r).

## Two physical routes to the same carry

At the lower high-product state with masks (a−1,b−1), there is
an alternate +2 move: mark both of the hitherto unmarked spacer
vertices and change nothing else. Write L=2D−b_outer and
R=D+a_outer for the current outer factors. The contribution of
a single contact j=0,1, before inserting the other, is
L·2^j−R·2^(1−j), minus one if there is an adjacent old mark
across the outer boundary. This follows from the signed linear
term minus the two geometric cross-contact sums. Here L−R=2
and both boundary marks are on, so adding both yields exactly
3·2−2−K(1)=2. It lands at the very same integer as the high-product
reset in the table. In fact these two tournaments are converses:
reverse the full spine order (including the sentinel) and complement
all its head bits, as explained in VARIABLE_RADIX.md.
In the balanced high state (a,b), both boundary marks are off and
L=R; marking the same two spacer vertices instead contributes −2.
So slack at the actual spacer has a sign that depends on the
two outer ports. This is not an argument that they are interchangeable
for subsequent extensions.

## Exact limit of this result

This closes one low-digit cycle, not a full sequence of consecutive
high addresses. The high reset always arrives at a square m²;
the next interval may start at an unrelated m'²−1. If m=4q
with q an odd prime and r≥5, the next state with raw high product
m² and remainder C+2 is impossible at that same cut. Either
factor having odd part 1 would be a power of two at least D,
but v₂(m²)=4. Thus both factors must have odd part q. Each allowed
factor interval meets a fixed odd part's power-of-two orbit in
at most one point, and m is already in both intervals. So the
same-cut square product has just its original state. This is
conditional on q prime and does not exclude changing the cut,
a different nontransitive sideboard, or a different high product
with a signed remainder. It gives a target for the next search.
