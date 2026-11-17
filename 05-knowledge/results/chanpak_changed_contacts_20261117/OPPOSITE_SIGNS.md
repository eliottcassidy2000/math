# Opposite-sign sevenfold readers with a common block

The [Mersenne seam](MERSENNE_MODULES.md) has a strong companion with the **opposite sign on the Hamiltonian count of the block**:

    H(final tournament) = 7(5P₂(B)+H(B)).

Both signs can be read from one arbitrary nonempty tournament B. Moreover, varying only the size of a transitive seam reader gives a matched sequence of increasingly different coefficients. This does not establish universal target selection.

## The physical positive-sign seam

Let h=H(B), P=P₂(B) with empty slots allowed, and R=P−2h. The positive left scaffold has outgoing masks (8,1,3,6); replace vertex 0 by B and mark vertex 1 as the terminal port. Vertex 0 has a single external exit, vertex 3, so there are at most two B-runs in either a full Hamiltonian path or an ordered cover by a first path ending at 1 and an optional second path. In the latter count only the second path could end inside B. Reading respectively the ordinary and terminal counts a,b_left at block TT₁ gives (5,9), and at TT₂ gives (11,17). Since one run contributes h and two contribute R, these checks determine all orders:

    a = 5h+3R = 3P−h,       b_left = 9h+4R = 4P+h.

For an actual right transitive tournament TT_r, r≥2, orient all arcs right-to-left except for the left terminal→right source seam. Paths that avoid the seam contribute a. Paths that use it choose a subset of the r−1 other right vertices to visit before the first left run, then continue through the right source and the remaining right vertices before the optional second left run. The exact count is a+2^(r−1)b_left. At r=4 it is 7(5P+h). The left frame 2→0→3→1 lifts by any Hamiltonian path on B and closes through the right block, so the final tournament is strong.

## A common family of right readers

The left negative scaffold [reads](MERSENNE_MODULES.md) (a,b_left)=(64P+h,96P−2h). Choose r=3j+1 for the positive reader and r=3j for the negative reader, with j≥1. With u=(8^j−1)/7 these physical quotients (multiply each by 7 for actual tournament counts) are

    x = (1+4u)P+u h,         y = (16+48u)P−u h.

Thus the right reader grows by three vertices and u grows by u↦8u+1. At a prescribed pair (x,y), the absolute inverse does not involve a search over block counts: P=(x+y)/(17+52u) and h=[x−(1+4u)P]/u. At j=1 this is P=(x+y)/69 and h=(64x−5y)/69. Integrality and positivity are necessary but not sufficient: the recovered pair must be the actual (H,P₂) of a tournament.

Literal simultaneous prime quotients at j=1 include (11,127), (641,8191), and (2011,25589), with (h,P) respectively (1,2), (1,128), and (11,400). The last is nontransitive: take masks (24,1,3,6,14) followed by two ordinal sinks. The identity 2^32+1=641·6700417 makes the middle pair a Fermat-factor/Mersenne-prime example; this is an exact arithmetic comparison, not a claim about infinitely many such pairs.

## Why the paired prime examples need an odd block at odd reader levels

The odd-cycle pack expansion gives H(B)=Σ 2^c and P₂(B)=Σ 4^c2^s over the same packs, with c directed odd cycles and s uncovered singletons. Since s≡n−c (mod 2), reducing termwise modulo 3 gives the all-order identity `P₂(B)≡(−1)^n H(B) (mod 3)`. This uses the [Rédei–Berge expansion](https://arxiv.org/abs/2307.05569) specialized to one and two path slots.

At odd j, u≡1 (mod 3), so (x,y)≡(2P+h,P−h). If n is even both are divisible by 3. If n is odd their residues are (−h,h); they avoid 3 exactly when 3∤h. At even j, u≡0, both residues are P. Thus selecting the reader height changes whether block parity locks both quotients out or makes them complementary. This is a joint prime *necessary* condition, not a prime-production theorem.

## Local reach does not give absolute placement, even in the enlarged negative bank

The same smallest negative reader has explicit prime quotients in every odd residue mod 8: 38113, 4091, 6389 and 8191. See the replay for the blocks. Nevertheless **65537 is missed by every negative reader above**, even allowing arbitrary B.

Indeed y≥32+95u, since P≥2h and h≥1. Thus only u=1,9,73,585 can matter at this target. A singleton B reads exactly 32+95u, never 65537. At order n≥2, P≡0 (mod 4): pair each subset with its complement in the sum P=Σ_S H(B[S])H(B−S), with each summand odd and an even number of complement pairs. Writing A=16+48u, we need `uh≡−65537 (mod 4A)`. The multiplier u is invertible because u is odd and gcd(u,1+3u)=1. The least positive residues are 255,1479,1031,96519 respectively. In particular h>1, so n≥3. The split bound P≥(n+1)h then gives h≤floor(65537/(4A−u)), namely 257,36,4,0. Only h=255 at u=1 survives. Since H(B)≤n! and 5!<255, n≥6; the split bound then gives y≥447·255=113985, contradiction.

This is an all-order exclusion from the **negative reader bank**, not from the tournament spectrum and not from the positive bank. It shows why full residue reach cannot on its own settle exact placement.

Replay: `python replay_opposite.py` (same results under `-O`). The arithmetic obstruction is proved independently of the literal witness checks.
