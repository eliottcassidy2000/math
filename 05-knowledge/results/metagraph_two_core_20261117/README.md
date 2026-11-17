# Two-core metalines: recognizing the seam and escaping mirror parity

This extends the [three-coordinate tail decoder](../metagraph_tail_tomography_20261117/README.md) to a movable right core. It supplies an explicit infinite asymmetric family in both of the sevenfold quotient classes missed by mirrors, with literal prime examples, and gives an exact score certificate for recognizing an existing seam. It does not prove prime infinitude in the families or a universal selector.

## Recognition and reading

Take two nonempty framed cores A and B. The first frame ends at v and the second starts at u. Between the cores orient all arcs B→A except v→u. This is a literal framed metaline: reverse the nonpath arcs within the cores, take their ordinal sum in A,B order, and cross the line along the concatenated frame.

Use terminal profile (a,b,c) on (A,v): H(A), covers by one or two ordered paths with the first ending at v, and full paths starting at v. Use (d,e,f) similarly on (B^op,u). The count of the resulting tournament is

    H = ad + be − cf.

Without the exceptional arc, a full path runs through B and then A, except for the forbidden transition u→v, eliminating cf paths. With the arc, the two marked cover choices are independent: the portions before their distinguished runs or after them may be absent. This gives be. Sizes at least two ensure the concatenated frame closes to a cycle, makes the low-outdegree set exactly A, and splits the automorphism group into the two port stabilizers. A's scores are all below |A|, since the port v already has a predecessor in its frame; each score of B is at least |A|, since u has a successor.

Conversely in a strong tournament, with increasing scores s_i, let the k lowest-score vertices make a candidate set A (including any necessary choice across a tie). Then

    sum_(i≤k) s_i − binom(k,2) = number of arcs A→V\A.

An interior cut of excess one, with both sides of size at least two, is exactly a two-core seam. Its sole outgoing arc must occur in every spanning directed cycle, so the cycle traverses each core in one interval and gives compatible frames. The degree separation above then shows there could not in fact have been a tie across that cut. Thus all such cuts are canonical, nested score prefixes. More generally, reversing their sole outgoing arc destroys strong connectivity. For a fixed arc, precisely the cuts it opens become boundaries of the ordered SCCs after reversal, with the cut interpreted as a final segment of the resulting SCC order. In the reverse direction, every such final segment has only that outgoing arc in the original strong tournament. This describes the excess-one cuts as data from critical individual arcs, not just a score coincidence. This version of recognition allows singleton cuts if their actual sets are retained; score ties need not be unique there.

## Mirroring and its true arithmetic obstruction

Taking B=A^op and u=v, with the reversed reflected frame, makes the line blue and gives H=a²+b²−c². Rédei parity implies a odd and b≡c mod 2 (also from the singleton-tail decoder), so every such count is 1 mod 4. For H=7p, this forces p≡3 mod 4. Moreover if c=0 and 7 divides H, then 7 divides both a and b (−1 is a nonsquare mod 7) and H/7 can only be prime if it equals 7 and (a/7,b/7) is a unit-axis pair. Both a,b are positive here, so even that exception is impossible.

There is no such barrier for a nonmirrored seam. Writing b=c+2E and e=f+2F gives

    H = ad + 2(Ef+Fc) + 4EF.

The mod-4 change is the mixed contact Ef+Fc, not a larger copy of the same mirror. Below is a single switch of port that controls it.

## A two-parameter family with different ports

Let R_j have j+2 vertices ordered 0,…,j+1 and every arc increasing except j+1→0. For j≥1, its H is 2^j+1: either no reversed arc is used (natural path) or a freely chosen subset of the j interior vertices is placed before that arc, with the other vertices after.

At terminal port 0 the profile is (2^j+1,2^j+1,1). There is exactly one path starting at 0, namely the natural path. Any two-cover ending its first run at 0 concatenates: if the terminal run is not the singleton 0 then its predecessor is j+1 and the second run starts at an interior vertex; if the first run is the singleton 0, the second is the natural path on the rest and starts at 1. Splitting a Hamiltonian path at 0 inverts concatenation.

At terminal port 1 the profile is (2^j+1,2^j+2,2^(j−1)). A path starting at 1 must use j+1→0 and it freely chooses which other interior vertices appear before that arc. A two-cover ending its first run at 1 fails to concatenate exactly when the second starts at 0, the only vertex beating 1. In this case the first run avoids 0 and has to be the singleton 1, because the induced order on the other vertices is transitive with 1 first. There is exactly one second run starting at 0 (natural order with 1 removed). Both ports are eligible by the spanning cycle.

Put the port 1 profile on the left using R_j, and port 0 on the right using R_k^op. The exact count is

    N(j,k) = (2^(j+1)+3)(2^k+1) − 2^(j−1).

For every positive j divisible by 3 and k≡2 mod 3 this count is a multiple of 7, since it reduces to 5·2^k+1 mod 7. For k≥5:
- j=3 gives N=19·2^k+15, with N/7≡1 mod 8. k=5 gives 7·89; k=2 (no need for the mod-8 claim) gives 7·13.
- j=6 gives N=131·2^k+99, with N/7≡5 mod 8. k=5 gives 7·613.

Thus a fixed change of port plus an independently chosen right size gives both classes excluded by mirroring, with actual primes, not just formal residues. This is an exact infinite family of physical counts in the 7-divisible channel, not a claim that the quotient is always prime. Its sparse fixed-core rays cannot prove that every prescribed p is hit. The meaningful next inverse is to let an excess-one cut choose *both* pointed profiles, rather than just grow k.

`python replay.py` independently counts the three displayed tournaments by subset DP and verifies their excess-one cuts, then verifies the claimed 343 sporadic seed by its literal adjacency rows. The multiplication-by-3 chain beginning at 63 is from the preceding checkpoint.
