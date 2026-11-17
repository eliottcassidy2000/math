# A smaller core with two runs and independent neutral addresses

**Construction.** The 13- and 15-path tournaments below have the same one-root continuation data. Each appended contact can be chosen independently, the count difference stays two, and the appended tape has at most two path runs. The deletion count at the root is **3**, the smallest possible for such a pair in this continuation system. A separate one-live-port *edge switch* uses two compensating addresses; see [balanced tapes](BALANCED_TAPES.md). Neither yields universal selection.

## 1. Replace the two live ports by a rooted pair

Let r be a retained root, g=H(T−r), and E,O the weighted counts of r-starting paths of even or odd vertex length, where a path is weighted by the Hamiltonian count of its complement (the empty complement counts as 1). For the shared continuation alphabet:

- S: append a vertex losing to everyone.
- N: append a vertex beating just r and losing to every other vertex.

The exact transfers, in coordinates (H,g,E,O), are S=(H,g,E+O,E+O) and N=(H+2E,g,E+O−g,E+O). See [the curried rooted reader](../curried_births_tooth_diagram_20261117/README.md) and [the earlier 41/43 binary controller](../birth_fibers_binary_growth_20261116/README.md). The construction here improves the deletion core, not the inherited transfer.

The following two order-five hosts, with outgoing masks read in order and root 4, have profiles (13,3,6,9) and (15,3,6,9):

    (26,12,9,16,6)        (14,20,24,18,1).

Apply any shared word in S,N retaining the root. The profiles still differ only in H, by exactly 2. Removing the root leaves the original deletion tournament followed by a transitive tape. In a Hamiltonian path each side of the root can visit that tape in only one increasing run. Thus at most **two** runs are available, even when every one of the m tape contacts is freely chosen. These are 2^m physical address choices, not a claim that their absolute centers are all distinct.

A second pair, with the same root convention, has profiles (15,3,9,12) and (17,3,9,12):

    (62,36,40,34,14,16)      (42,60,57,32,9,16).

It comes from making one N move on the (3,5) pair at the common profile (g,E,O)=(3,6,6).

## 2. Why this core works, and where the next width stops

Here is a classification at the two smallest deletion counts. A tournament of count one is transitive. If T−r is transitive, let b_i be the bit for r beating its i-th vertex. For its final bit b_last, O−E=1−b_last. The remaining contributions to E are distinct descending powers of two, plus the known b_last. Equal E,O thus recover the entire bit word up to leading zeros. Those vertices are global sources and have no effect on H. So g=1 cannot give a nonzero difference with the same continuation.

At g=3, the deletion is a directed triangle between transitive prefix and suffix. Indeed the [odd-cycle expansion](https://arxiv.org/abs/2307.05569) forces exactly one triangle, and every vertex in a nontrivial strong component lies on a triangle. Let k=0,1,2,3 be the number of triangle vertices beaten by the root; rotations identify all choices of the same k. Replace that triangle in the address by the two binary digits u,v of k=2u+v, keeping ordinary contact bits of prefix and suffix. Call the resulting binary integer b, allowing leading zeros. Let P be the integer prefix before u,v and let t=1 if the suffix exists and begins with zero, else 0. Then

    E+O = 3(1+b),
    O−E = 3(1−last bit) if there is a suffix, else 3−k.

For the signed parity read, a root-starting path at a prefix singleton has zero contribution to O−E because the weighted odd/even subsets of the triangle balance. At a triangle vertex the signed contribution is −1 when the suffix is empty and zero otherwise. In the suffix only the last vertex has a nonzero signed contribution (−3). The empty old start, i.e. the root alone, contributes 3. For the unsigned read, a prefix singleton contributes 12·2^(later singleton count), each triangle vertex 3·2^(suffix length), and the suffix singleton at i contributes 3·2^(later singleton count). This proves both equations.

For the count itself, let h_0 be the one-root Hamiltonian count with the *compressed* binary word over an entirely transitive deletion. Put L=3 floor(P/2)+(P mod 2). The physical triangle satisfies

    H = 3 h_0 + 2e,
    e = 0                 if uv = 00 or 11,
        L−t+1             if uv = 01,
        −L+t−2            if uv = 10.                         (1)

**Counting proof.** Subtract the paths counted by the deletion, using odd cycles through r: their other vertices are an even directed path from a root outcontact to a root incontact, with weight H of the remaining deletion. For an ordered pair of endpoints in the transitive word with d vertices between, the physical (unweighted) count is 2^max(0,d−1). With both endpoints outside the two encoded slots, the triangle case gives exactly three times that count. If the interval crosses the triangle, optional triangle subsets contribute total 12 with even/odd totals 6 each: an unused triangle weighs 3 and a full directed triangle contributes its 3 paths.

For a single triangle endpoint, fixing it gives exactly one directed path ending or starting there at each subset size 1,2,3 (summing over subsets of that size). With d outside optional vertices the weighted count for that endpoint is 2 at d=0 and 3·2^(d−1) otherwise. Comparing with the two compressed endpoints, an external prefix outcontact contributes v−u at distance zero from the triangle and 3(v−u)·2^(d−1) at larger distances. Their sum is (v−u)L. On the suffix side everything cancels except the adjacent zero contact, which contributes (u−v)t. With both endpoints inside the triangle, the count is one if k=1,2 and zero otherwise; the compressed diagram instead contributes three for uv=10 and zero otherwise. These corrections give (1).

**Sharp width.** At fixed g=3 and fixed E,O, b is fixed, so all count choices come from sliding the two-slot position in one binary word. Write the offsets as e. Distinct 10 positions do not overlap and their prefix integers satisfy P'≥4P+2; their negative offsets are separated by at least 2. For 01 positions, P'≥4P+1, and nonnegative offsets are separated by at least 2 except possibly the first two positions with P=0,P'=1. In that exception the word starts with 0101 up to leading zeros; the first offset is 0, and the second is 1 or 2, not both. If P=0 and P'≥2, the larger offset is at least 3, and if P≥1 it is farther still. A cross-sign consecutive triple would need −1,0,1. The offset −1 requires the first significant bits 100. But +1 requires either P=0 with the first significant 1 followed by 1 (or ending the word), or P=1 with first significant bits 101. They cannot coexist. This proves at most **two consecutive odd H values** in a g=3 continuation profile, even allowing every triangle position and every word length. The displayed pairs attain it.

In the normalized sector O−E=g=3 the least E supporting a unit pair is 6. Indeed E=3b/2; the smaller possibilities b=0 and b=2 have only one eligible count each by (1) and the parity read. Thus the first displayed pair attains the minimal deletion count and minimal normalized slope within it.

## 3. Exact placement tests

Use address z with bit i set when the i-th continuation (chronological order, starting at zero) is N. Let f(0)=1, f(2u+1)=f(u), f(2u)=f(u)+u+(u mod 2) for u>0. In the normalized profile O=E+g, the [binary inverse](../prime_packing_inverse_20261116/README.md) gives H(z)=2E(z+1)+g f(z)+(H(0)−2E−g). The first bank here gives 3[4(z+1)+f(z)]−2 and the same output plus 2. The second gives 3[6(z+1)+f(z)−2] and that output plus 2. In particular:

| Prescribed endpoint | Bank | Address |
|---:|---|---:|
| 343 | 13/15 | 25 |
| 10009 | 13/15 | 761 |
| 69 | 15/17 | 3 |
| 6017 | 15/17 | 313 or 327 |
| 30011 | 15/17 | 1561 or 1639 |
| 1000000007 | 15/17 | 47827936 |

The billion-prime continuation slope and address agree with the earlier single-host witness. The added result is the matching lower count and permanent exact unit over all future words. This is not a new claim of the prime's attainability. The inverse in replay.py checks all addresses of these banks at each chosen target by splitting high/low digits and bounding the remaining f-cofactor. Neither bank is universal; even counting with multiplicity both together have ≤(1/6+1/9)X+O(1) endpoints up to X, by their linear lower bounds in z.

## 4. If the seam must remain a single edge

The construction above uses two physical *core* choices rather than requiring them to differ in just one edge. The separate [one-live-port edge switch](BALANCED_TAPES.md) preserves a literal reversal: make one port inactive and pay the left-right drift by choosing binary addresses (6t,4t). Its exact count grows quadratically, with Θ(√X) distinct endpoints. There is no third run there either, but the address correlation is an additional cost. This distinguishes two ways to escape the old three-run obstruction, and identifies what still has to change to obtain universal placement.

Run python replay.py for the general triangle formula, literal profiles and inverse controls; python verify_balanced.py checks the separate single-edge construction.
