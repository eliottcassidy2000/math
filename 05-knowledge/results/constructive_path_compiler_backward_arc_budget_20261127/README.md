# Compile Hamiltonian paths by deleting to a total order

2026-11-27. This generalizes the repository's [two-reversed-arc run assignments](../../../04-computation/experiments/repunits_hyp9357_20261021/two_arc_formulas.py) and complements [one-apex and finite-core interfaces](../tournament_computational_interfaces_fractran_20261127/README.md). The input/output below concern **exact** Hamiltonian path counts, not residues. Chan–Pak remains open.

## 1. Three states per gap, exact counting and unranking

Let F be q exceptional vertices such that U=T−F is transitive. Give U its forward total order. Enumerate permutations p of F. A Hamiltonian path with that exceptional subsequence is uniquely specified by assigning every vertex of U to one of q+1 slots before/between/after the exceptional vertices, then reading each slot in U order.

For an empty internal slot the two exceptional ends must have the forward arc in p order. A nonempty slot only needs the incoming arc to its first vertex and the outgoing arc from its last vertex, where present. It would be **wrong** in a general feedback-vertex set to impose both incidences on every vertex of the slot: the intervening vertices may have arbitrary orientations to F.

Scan U in order. A slot stores a trit: empty, nonempty with an allowed last arc, or nonempty with a forbidden last arc. On first entry require the left incidence; every subsequent entry is permitted and updates the last-arc trit. At the end reject forbidden last arcs and empty slots with a missing direct arc. The scan is a bijection between accepting assignments and Hamiltonian paths with subsequence p. With suffix completion counts, one can unrank: try p's and then slot choices in fixed order, subtracting the size of each preceding branch. This gives an exact uniform sampler as well if the input rank is uniform.

Counting and building the unranking table uses O(q!·3^(q+1)·(q+1)·n) arithmetic operations beyond reading the tournament. Counts fit in O(n log n) bits, so this is FPT in q, including bit complexity. If a promised bound q is given but F is not, find a directed triangle, branch on deleting one of its three vertices, and repeat for at most q levels; triangle-free tournaments are transitive. This elementary search is not a claim about the fastest feedback-vertex algorithm (see [Kumar–Lokshtanov](https://arxiv.org/abs/1510.07676)). There is no bounded-q claim for general input tournaments.

## 2. Exact inverse search at a fixed backward-arc budget

Suppose an order has at most b≥1 backward arcs and let E be the ordered set of their at most 2b endpoints. The ordinary vertices are distributed among the gaps in E; each gap is a transitive chain. Here its vertices have forward incidences in the given order to **every** exceptional vertex, so for fixed p an ordinary vertex may go into a slot exactly when its index lies between those of its two exceptional ends (using exterior sentinels). Indeed the first and last inequalities in a nonempty slot force all vertices in it to lie between the ends. Every gap has one eligibility mask, independent of its length.

For p, let R be the set of internal slots whose exceptional ends have no forward arc. Such slots must be nonempty. For any D⊂R let c_g(D) count the eligible slots for gap g once slots in D have been disallowed. The contribution of p is exactly

    Σ_(D⊂R) (−1)^|D| Π_g c_g(D)^(length_g),

with 0^0=1. Sum over p. This is compiled before the chain lengths are specified. Before deleting D, the base c_g is 1 plus the number of descending jumps of p across g: scan the permutation from the left sentinel to the right sentinel and count signed crossings of the cut. A permutation with a descending pair not backed by an actual backward arc contributes zero (that slot has no allowed ordinary vertex either); the remaining bases are at most b+1. Thus every fixed pattern has a finite sum of products of explicit small-base exponentials. In each gap variable the sequence separately obeys the linear recurrence with characteristic polynomial dividing Π_(j=0)^(b+1)(x−j).

A gap crossed by no backward arc can be set to zero without changing H: all paths traverse it as an ordinal chain, between the portions before and after it. This removes exterior padding as well. For an active gap of length l, induce on that gap and the two endpoints of a backward arc crossing it. The resulting near-transitive tournament has 1+2^l paths. Inserting one vertex into a Hamiltonian path of an induced subtournament is always possible by the ordinary first-fit insertion rule, and deleting that vertex recovers the original path. Hence H≥1+2^l, including if the full tournament is not strong.

Therefore **exact membership through height X at budget b** is obtained by enumerating the finitely many ordered endpoint patterns and 0≤l≤floor(log₂(X−1)) on their active gaps. At most 2b−1 gaps are active. This gives O_b((1+log X)^(2b−1)) candidate tournaments and distinct heights through X, and reconstructs a tournament at each successful target. It explains why any fixed number of edits of a total order cannot prove Chan–Pak. The algorithm only rules out the bounded-backward-arc model, not tournaments of bounded feedback vertex number in general.

## 3. The two-arc logarithmic exponent is sharp, with a direct inverse

Use the inherited crossing formula with positive ordinary gaps p,r,q, and restrict them to 1≤r<p<q. For two reversed arcs whose endpoints interlace,

    H−1 = 2^(p+r+1) + 2^(q+r+1) + 2^(p+q+1)·3^r.

The three 2-adic orders are strictly increasing in that displayed order. Given H, strip its lowest power of two after subtracting 1, then strip the lowest power of two from the remainder. The final nonzero remainder has 3-adic valuation r. Subtract r+1 from the first two 2-adic orders to recover p and q. Validate the inequality and the formula; this gives a deterministic inverse on this subfamily and rejects all other targets. Taking 1≤r<p<q≤c log X with a sufficiently small absolute c puts H≤X and gives Ω((log X)^3) distinct heights. With §2 the two-arc spectrum has exactly order Θ((log X)^3) growth. For comparison b=1 has the familiar 2^l+1 spectrum of order log X.

## 4. The entire third chain needs three or more backward arcs

Take the chain 7·3^t. The root and first step 7,21 are the known unrestricted gaps. Every later term is attainable (begin with the 8-vertex one-apex H=63 witness, then ordinal-join copies of the directed triangle). But no term on the chain can arise with ≤2 backward arcs. Use [THM-4748 and its formulas](../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md). We may assume a minimum order, which has no reverse arc between consecutive positions since swapping them would reduce the budget.

For even t≥2, the target is 7 mod 8. The known two-arc candidates in that residue are 15,23 and 3·2^j−1,5·2^j−1. The target plus one is not divisible by 3, nor by 5 since 3^t=±1 mod 5. One arc cannot give a multiple of 7 because 2 has odd order mod 7.

For odd t≥3 the target is 5 mod 8 and divisible by 9 and 7. For clarity, here is the residue elimination in the exact two-arc formulas:

- Disjoint or touching reversed intervals give a product of two factors 2^j+1, neither divisible by 7.
- With a shared left or right endpoint, put B=s−1, C=L−1−s in the notation of the file above. For C=0 the height is 3·2^B−1, not divisible by 3. For C≥1, it is 1+2^(B+C+1)−2^C+2^B. Residue 5 mod 8 allows only B=2,C≥3 or C=2,B≥3. The first gives 5+7·2^C, not divisible by 7; the second gives 9·2^B−3, not divisible by 9.
- With strict nesting, write d=s−1≥1 and w=x+y≥0. The height is 1+2^d+2^(w+1)(2^(d+1)−1+3^d). Modulo 8 the only possibilities for 5 are d=2 (giving 5+2^(w+5)) or odd d≥3 with w=0 (giving 5·2^d+2·3^d−1). In the first, divisibility by 7 requires w+5≡1 mod 3, whereas divisibility by 9 requires w+5≡2 mod 6. In the second, d≡1,3,5 mod 6 gives residues 1,2,1 mod 7.
- With crossing intervals, take p≤q without loss. If p≥1, residue 5 mod 8 requires r=0,p=1,q≥2, giving 5+3·2^(q+1), not divisible by 3. For p=q=0 the value is 2^(r+2)+2·3^r−1; residue 5 requires r=0 or r odd, neither gives divisibility by 3. For p=0,q≥1 write the formula as 1+2^(r+1)+2^q(2^(r+1)−1+2·3^r). If q≥3, residue 5 forces r=1 and the value is 5+9·2^q. If q=2 it forces r≥2 and the value is 5·2^(r+1) mod 3. If q=1 it is −1 mod 3 for r≥1 and is 9 for r=0. None is a target.

In particular the exact minimum at 63 is three: in order 0,1,…,7 reverse just the arcs from 0 to 3,4,7. For 7·3^t, t≥2, the construction gives upper bound t+1 on the minimum. This leaves open whether **one common budget of three** suffices for the whole attained tail. It is a constrained inverse problem inside the constructive search of §2, not a claim that the rest of the odd spectrum fits a fixed budget.

Run **node replay.js** to verify counting and unranking against direct path permutations, test the 80-vertex sparse case, check the three-arc witness at 63, and replay the arithmetic inverse for the crossing family. The replay uses safe-integer trit indexing; the mathematical algorithm has no such implementation cutoff.
