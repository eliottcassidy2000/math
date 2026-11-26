# Chan–Pak: ten million exact witnesses and a low-spoke obstruction

2026-11-26. **FINITE-EXACT** interval certificate and **PROVED** scoped obstructions. [Chan–Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf) remains open.

## 1. The detached interval

Every odd integer in **[34,359,738,369, 34,379,738,367]** is the total Hamiltonian-path count of a tournament on at most **39 vertices**. There are 10,000,000 targets. The [November 18 selector](../chanpak_prime_frontier_selector_20261118/targeted/PROOF.md) certified the first 1,000; the present interval goes 9,999,000 targets farther. It is an unconditional interval. Joining it to a full prefix depends on the separate [reported THM-4885 census](../../../01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md) through 2³⁵−1, which this packet does not re-run.

The program uses F1 (two rooted binary blocks) first, with the smaller positive address capped at 2048, and fills **9,942,098** targets (largest order 38). The remaining **57,902** use F2 (one rooted binary block and an actual bank tournament). The original search used an ordered subset of the inherited [independent S203 bank](../../../session203/coverage-recertification/session203-independent-audit/aud_bank7_inherited.txt); the compact replay keeps exactly the 377 profiles that supplied a first witness, in original order, with their literal outgoing masks. Unused profiles are not needed. A miss at either template stage is not a missing tournament count.

Run **bash replay.sh** from this directory (C++17; writes about 245 MiB temporarily, then removes it). The two programs are deliberately different. The producer uses the binary recurrence and the suffix-window arithmetic below. The verifier requires the entire interval in order, recomputes the binary count by its weighted inversion sum, and checks every rooted bank profile by listing actual paths of all subsets of the saved matrices. It also checks their orientations. The generated witness stream itself has a fixed SHA-256. The profile files retain root 0, with bit j in row i denoting i→j.

## 2. Exact windowed inverse

Let C(v) have a transitive spine in increasing order and a root x which loses to spine vertex i precisely when bit i in v is 1. Use the bit length of v for the spine order. Extra high zero bits are permissible for comparing same-order flips. Write its count F(v). In descending bit order, every 1 before a 0 contributes an inverted-pair term: for distance j the weight is 2^(j−1)+[j=1]. One way to obtain it is the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569): all directed cycles use the apex, and a pair of extremal spine vertices at distance j supports an even number of interior vertices in the cycle. This gives F(0)=1, F(2q+1)=F(q), and F(2q)=F(q)+q+[q odd]. Also 1≤F(v)≤v+1 by induction (when q is odd use F(q)=F(⌊q/2⌋)≤q).

For a right block G with root y, let h=H(G), e=number of spanning paths ending at y, and R count ordered covers by an optional path and a path starting at y. Orient all interblock arcs right→left except x→y. The count is

    h F(v) − e[v even] + R(v+1).

Without the bridge, paths are an ordinal pair, omitting exactly those with y at the end of the right path and x at the start of the left path. The latter root-start reading is [v even]. With the bridge, delete the edge x→y and read the four possible runs: optional right, left ending at x, right starting at y, optional left. A left cover with its first slot ending at x has v+1 choices: one with the slot equal to x, and 2^i choices for an old maximal predecessor i with bit i set. The right slots give R. For the converse of C(b−1) as the right block, this specializes to

    H(b,c) = bc + F(b−1)F(c−1) − [b odd][c odd], c=v+1.

Swapping address roles preserves the scalar height by this identity; the physical reversal analysis below keeps the roles explicit.

Write v=Lq+r, 0≤r<L for a dyadic L, and let Z be the bit-reverse complement of the padded suffix r. The inherited concatenation identity is visible directly from the inverted pairs: the two internal sums give F(q)+F(r)−1, the crossing geometric weights give qZ, and the boundary adjacent surcharge is [r<L/2][q odd]. Thus the same integer N has the exact linear-in-prefix window

    N = K+Bq+h (F(q)+[r<L/2][q odd]),
    K = R(r+1)+h(F(r)−1)−e[r even],   B = LR+hZ.

For an entire interval [X,Y], every possible q therefore lies between the ceiling of (X−K−2h)/(B+h) and the floor of (Y−K−h)/B, also respecting nonnegativity and any minimum v. Enumerate these q and check the exact integer; this reuses each suffix over the entire interval instead of solving for each target separately. The padded reverse complement can be updated by a most-significant-bit-first decrement with borrow. In F1, K≥0 and h≤b=R. With L the least power of two ≥2 satisfying L²b≥Y, the real q-window has width at most (Y−X+h)/B+hY/B², unless it is empty. This is O(1+(Y−X)/(bL)). Summing over suffixes and b≤B₀≤√Y gives O(√(YB₀)+(Y−X+1)log(2B₀)) candidates or empty-window checks before the O(log Y) cost per reading. This is an exact inverse for a finite bank, not a uniform image theorem.

## 3. Two lowest spokes cannot produce a growing fork in F1

A physical flip of the lowest address bit changes the root's arc with the first spine vertex. To include 0↔1, keep that vertex present and pad with a zero bit when necessary. Mark one lowest spoke on **each** side. The only possible central total with these two reversals giving opposite unit shifts (−2 and +2) is **3**, from (b,c)=(1,2) or (2,1).

**Proof.** Write b=2p+1+x, c=2q+1+y, with x,y bits; put α=p+[p odd], β=q+[q odd], f=F(p), g=F(q). The two physical increments at the central bits (x,y) are

    D_b = (1−2x)(2(q+1)−α(g+(1−y)β)),
    D_c = (1−2y)(2(p+1)−β(f+(1−x)α)).

If p=0 and |D_b|=2 then q=0, and the square of heights is (1,3,3,5). It remains to take p,q≥1. At x=y=0 neither increment can be +2, because α≥2 and g+β≥q+1, and symmetrically. At (1,0), the b increment is nonnegative; to be +2 it requires α(g+β)=2q+4. This forces α=2 since α≥4 would already exceed the right side. Then g=2−[q odd], so g odd forces q odd and g=1. A reading F(q)=1 means no inverted pairs, so q=2^j−1. The other exit would require βf=2p+4 with p=1 or 2. This is impossible: β is a power of two, and f is respectively 1 or 3, giving β=6 or 8/3. Interchange the roles for (0,1).

At (1,1) we can exchange the two arms to arrange αg=2q and βf=2p+4. Since α≥p and β≥q, fg≤4+8/p. For p=1 we already have f=1, and for p≥2 the product is below 9. Thus f=1 or g=1. If f=1, p=2^i−1; q is either 2p+4 or 2p+3. Accordingly g is 4+4/(p+1) or 4+2/(p+1). The only positive p yielding an odd integer in those cases gives respectively (p,q,g)=(3,10,5) or (1,5,5), but F(10)=9 and F(5)=3. If g=1, q=2^j−1 and p is 2q or 2q−1. Then f is 4 or 4−2/(q+1). The only odd possibility requires q=p=1 and f=3, while F(1)=1. This proves the exclusion on every order. It does not exclude two spokes at other bits, two flips on the same side, or larger core banks.

## 4. A perfect-number barrier for a fixed left core and its dyadic padding

A single right lowest-spoke edit has central pair c=2q+1 and c+1. For b=2^s d with s≥1 and d odd, the preserved left count is h=F(d−1). The signed increment is b−hβ with β=q+[q odd]. Consequently the complete unit-edge list in this lane is obtained by testing the two even-integer candidates β=(b±2)/h and the at most two prefixes q=β or β−1 (keeping q≥0 and q+[q odd]=β).

Let M be the midpoint of such an edge, and E_p=2^(p−1)(2^p−1). The size of the correction is bounded by

    |M − 2b²/h| ≤ 15b.

Indeed M differs by 1 from b(2q+1)+h(F(q)+β), while |q−b/h|≤1+2/h, hF(q)≤b+2+h, hβ≤b+2 and h≤d≤b. Thus for any fixed d an infinite sequence with M=E_p forces p−1−s bounded and, on a subsequence with that difference constant, 4^(p−1−s)=d²/h. Since d and h are odd, the only equality is 1=d²/h; since h≤d, d=1. This is an arithmetic obstruction on actual shifts, not just on their modular phases.

When d=1 and b≥8, the four possibilities at q=b−2,b−3,b+2,b+1 have odd readings before the flip respectively

    2b² − 3b/2 − 1,   2b² − 15b/4 − 1,
    2b² + 27b/4 + 5,  2b² + 9b/2 + 3,

and the first two flips give +2, the last two −2. These follow from the binary recurrence. In the forced p=s+1 branch none has midpoint 2b²−b. So no fixed d crosses E_p infinitely often via this lane.

**Explicit at exponent 61.** For every p≥61 the same lane fails uniformly for every odd d≤32768. If M=E_p, the error bound and h≤32768 first imply 2^(p−2)≤b<2^p√h: smaller b would give M≤2b²+15b<E_p, and at or above the upper limit b≥15h gives M≥b²/h≥2^(2p). If the two leading terms 2b²/h and 2^(2p−1) are unequal, their difference is a nonzero integer multiple of 2^(2 min(s,p−1)+1)/h, hence at least 2^(2p−3)/d³. For p≥61 this exceeds 15·2^p√h+2^(p−1) using d,h≤2^15, contradicting M=E_p. Equality of the leading terms again gives d=1 and p=s+1, excluded by the four exact rays. This rules out only this prescribed lane of E₆₁±1 crossings, not general F1, other edits, or other tournaments.

Run **python3 low_spokes.py** for 289 literal seam counts, 262,144 finite fork settings, and the four arithmetic rays at exponents s=3 through 80. The finite controls are separate from the proofs.

## What remains

The direct numerical advance is a longer exact band, not a height-independent construction. The negative result identifies a specific place not to look for a growing ternary fork: the two first contacts in the simple F1 seam. If extending the current inverse into a proof, the useful next check is a bank with additional physical ports whose shifts depend on independent suffix readings. Two choices that collapse to the same lowest-spoke coefficient cannot supply this missing degree of freedom.
