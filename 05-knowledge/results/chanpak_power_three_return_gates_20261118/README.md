# The odd centers 81 through 2187, and an unbounded plus lift

This extends the [odd-center-27 construction](../chanpak_odd_center27_20261118/README.md).
The next pair at **81** is proved, as are both signs at **243, 729 and
2187**, at every dyadic width m≥1. One fixed marked seed supplies
**all plus counts 3^j·2^m+1 for j≥3 and m≥3 by induced growth**.
There is no matching unbounded minus mechanism here, and no claim
of universal Chan–Pak coverage.

## 1. One scalar-read formula, several exact centers

Use a marked core (D,q), where q is the mask of vertices beaten by
the apex v. Write h=H(D), κ=H(D+v). Let R count ordered covers P,Q of D
by two paths, allowing Q to be empty, with P nonempty and ending at a
vertex that beats v. Thus W=h+R also allows P to be empty. Let R⁻ be
the same count as R weighted by (-1)^|P|.

Precede D by a transitive chain L, with every L vertex beating D.
The apex tag on L is the binary expansion of a positive odd a
followed by t zeroes. The apex's contacts on D are q. Define F(0)=1,
F(2b+1)=F(b), F(2b)=F(b)+b+[b odd] for positive 2b. Then for t≥1,

    H = a W 2^t + κ - h(a-F(a)).

For t=0 the same expression has additional correction -h-R⁻.
Proof: the [Grinberg–Stanley formula](https://arxiv.org/abs/2307.05569)
counts odd-cycle packings, each cycle weighted 2. All cycles that
touch L also touch v. Packings avoiding L contribute κ. A v-cycle
wholly on L contributes h(F(a2^t)-1) in total: the signed binary
recurrence F counts those apex cycles plus the empty packing. A
crossing v-cycle starts at a winning tag on L, chooses a subset of
later L vertices and traverses a nonempty path P in D ending at an
in-contact. Each start has a free slot on L for t≥1, so the parity
restriction takes exactly half the subsets, canceling the cycle
weight 2. The crossing total is (a2^t)R. The recurrence gives
F(a2^t)=a2^t+F(a)-a+1. At t=0 just the final winning tag has no
free vertex; it changes the crossing total by -R⁻, while the
formal F expression was one too high, giving -h as stated.
The one-apex recurrence also follows by appending a final losing
tag: a start with d intervening vertices contributes 2^d+[d=0].

Here are literal fixtures. Rows are outgoing bit masks and q indexes
the same zero-based labels. The number e is defined by aW=3^j·2^e;
choose t=m-e. Each used family is strong: D can reach the apex from
every vertex, and the apex reaches all D via L.

| Center | Sign | a | D rows | q | h | W | κ | e |
|---:|:-:|---:|---|---:|---:|---:|---:|---:|
| 81 | − | 3 | 2,24,3,21,37,15,63 | 47 | 25 | 108 | 49 | 2 |
| 81 | + | 3 | 2,4,9,3,15,31 | 4 | 5 | 108 | 11 | 2 |
| 243 | − | 3 | 50,32,3,7,14,28 | 47 | 37 | 81 | 73 | 0 |
| 243 | + | 3 | 30,100,96,102,110,1,33 | 100 | 47 | 324 | 95 | 2 |
| 729 | − | 3 | 18,44,89,65,106,77,3 | 95 | 109 | 243 | 217 | 0 |
| 729 | + | 3 | 60,113,74,34,12,20,57 | 63 | 109 | 243 | 219 | 0 |
| 2187 | − | 9 | 38,4,104,3,79,26,43 | 23 | 81 | 486 | 323 | 1 |
| 2187 | + | 9 | 20,109,120,1,42,9,57 | 11 | 59 | 486 | 237 | 1 |

For 81 both have R⁻=-h and so cover m=2 on the very same templates.
The remaining finite boundaries 161,163,487,973,4373 have literal
rows in the replay. For 4375 use 35×5³ by ordinal gluing, supplied
with all factor witnesses. An additional shorter 81 minus fixture
with a=9 and h=15,W=144,κ=59 gives the same order m+7 as
the plus fixture for every m≥5. Adjacent
signs here need not be single-arc or blue metagraph mates.

## 2. Return capacity is a bipartite matching number

Let a distinguished vertex r of any tournament be replaced by a
transitive chain of m twins. Set ν equal to the maximum matching in
the bipartite graph of arcs x→y from out-neighbors x of r to
in-neighbors y of r.

For Hamilton paths the maximum possible number of visits to the
whole twin chain, once m is large enough, is exactly 1+ν.
For a cover by q labelled, possibly empty paths, it is q+ν.
To prove the upper bound, compress each maximal chain visit into a
copy of r. Every extra visit within a slot encloses a return excursion
out of r and back. Each excursion uses an out→in crossing, and
excursions are disjoint in the exterior. For the lower bound, use
matched edges as triangles r→x→y→r; concatenate the excursions in
one slot. Insert remaining exterior vertices by the usual tournament
path-insertion argument, which still works in a walk with repeated r:
an unused vertex can go first, last, or at a change from beating it
to being beaten by it. Add a chain visit in every other slot.

For a skeleton with j visits each onto assignment of the m chain
vertices gives exactly one way to fill the visits. Thus the count is
Σ_j A_j·j!S(m,j) with A_(q+ν)>0. In particular the exact leading
exponential base is q+ν. When ν=1, Hamilton growth has exactly the
dyadic form H=H₁+A₂(2^m−2). The return graph then has a single
common endpoint for all crossings. For the marked W-read below,
the same upper bound holds, and the lower bound is obtained by
putting excursions in the second slot and a lone chain visit in
the first, at an admissible end. This is why two readings of the
same physical edit can have exponential bases 2 and 3.

## 3. Select one pure mode in each reading: all plus powers of three

Use the plus core at 27, D=(2,20,9,3,13), q=13. Append k
transitive twins of the marked apex to the deletion core in place
of that apex. The mark itself is a new last twin, so every old
twin beats it. Call the resulting marked profile h_k,R_k,κ_k.
The exterior crossing graph has one common endpoint, vertex 1.
Hence with k old twins, there are at most two clone visits in h_k
and at most three in W_k. The latter is a two-slot read whose first
slot must end at an in-contact or be empty; each degree is still
in the same onto basis because all clones beat the new mark.
From exact values

    k       0    1    2    3
    h_k    11   23   47   95
    W_k    36  108  324  972

the one-slot and marked two-slot reads for k≥1 are uniquely
determined (triangular first values). Explicitly their onto
coefficients are (23,12) and (108,108,36). Therefore, including k=0,

    h_k=12·2^k−1,    W_k=36·3^k,    κ_k−2h_k=h_(k+1)−2h_k=1.

Take a=3 and t=m−2 in Section 1. This gives exactly 3^(k+3)·2^m+1
for all k≥0,m≥3. All hosts are strong: the marked part without L is a
transitive-twin expansion of the strong D+v, and v reaches the
left side via its first vertex. To see induced growth in k, insert a
new unmarked twin immediately before the unchanged marked vertex,
with all its cross-arcs to L pointing from L to the new twin. Combined with deleting one trailing zero of the flank,
this physically lifts the positive signed Collatz pair endpoint
from center 3^j to 3^(j+1), at width one lower (starting width ≥4).
This gives the plus half unboundedly but does not propagate both
mates under the same edit.

The old minus core at 27 has the same unit defect -1 under this
operation, but its marked read is

    W_k=24·3^k+16·2^k−4.

The coefficients follow by its same single crossing gate and the
four reads 36,100,276,772. The relative error shrinks, but the
absolute lower-mode error grows. This is an exact failure mechanism:
sharing the dominant base is weaker than having the same scalar
transfer.

## 4. Why a neutral extension cannot manufacture the odd multiplier

Suppose an induced extension of a marked core preserves κ exactly.
Each new vertex must be a singleton SCC of the full core with v:
a vertex in a nonsingleton SCC is in a directed triangle, adding
a positive odd-cycle packing not present before. Hence h is
also preserved. A new singleton before the SCC of v can join
either of the two slots in the W-read without changing validity
(the first slot's failed out-contact, if any, occurs later).
A singleton after the SCC of v must join slot Q: placing it in P
would make P end at an out-contact, since no later vertex beats v.
Thus W scales by exactly 2^s with s the number of new singletons
before v's SCC. In particular no κ-preserving induced extension
can triple W. The successful clone step simultaneously moves κ
and h by the dyadic rule, and only their **defect** remains fixed.

## 5. The unifying conjecture and a numerical consequence

A sharply testable strengthening of the two-sided power-of-three
ladder: find a minus seed and a repeatable induced move with W→3W
and κ−2h=-1 (or another bounded-state construction producing that
same read). The gate theorem controls the candidate bases; it
does not force cancellation of all lower bases. An exact search
should track the full low-mode signature, not just H near a target.
The completed plus half shows such cancellation is physically
possible. Neither direction alone settles the remaining odd prime
or 7p absolute targets of Chan–Pak.

There is an exact mod-7 byproduct: 2=3² mod 7, with 3 of order 6.
Thus at odd j, only the plus sign of 3^j·2^m±1 is divisible by 7
and m=(3−j)/2 mod 3. At even j, only the minus sign is divisible
and m=−j/2 mod 3. In particular the orientation flips from 27 to
81, where the first minus endpoint is 161=7·23. Each eligible
width class has the usual one-digit 7-adic lift under +3 in m
because 8=1+7. There is no assertion about infinitely many prime
quotients.

Run `python3 verify.py` for literal boundary fixtures, signed core
profiles, clone reads, and the first terms of all eight rails.
