# Exact prime-frontier witnesses, a terminating selector, and arc usage

This has three separate claims: new **finite** coverage, an **unconditionally terminating** decision and minimum-order procedure, and fast witness/profile algorithms on explicit families. It does not turn Chan–Pak into a solved theorem: the finite decision procedure may output absence, and its worst case is not fast.

## 1. Prime targets and two newly certified bands

The [corrected S203 independent replay](https://github.com/eliottcassidy2000/math/blob/main/session203/coverage-recertification/README.md) certifies the old prefix through 2^32−1. [THM-4885's addendum](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md) separately reports a single bucketed run through 2^35−1; that larger old prefix was not independently replayed. I found no higher prior contiguous certificate in the repo search. A frontier is the boundary of a named result, not a claim that every individual count above it was unconstructed.

[Targeted replay and witness parameters](targeted/PROOF.md) certify the next 1,000 odd targets above **each** cap: through 4,294,969,295 and 34,359,740,367 respectively. Every row materializes an actual strong tournament; the new high band can be checked independently of the global prefix. The first prime and first 7-times-prime just above each reference cap, plus the lone gap in the two-apex search in the high band, are:

| Target H | Role | Construction parameters | Witness order | Unused arcs |
|---:|---|---|---:|---:|
| 4,294,967,299 | 7 × 613,566,757 | F1(2,1502296101) | 34 | 1 |
| 4,294,967,311 | prime | F1(2,1847524150) | 34 | 0 |
| 34,359,738,371 | 7 × 4,908,534,053 | F1(9,2371737543) | 38 | 0 |
| 34,359,738,421 | prime | F1(71,297178177) | 38 | 0 |
| 34,359,739,333 | prime, F2 repair | rooted bank (2,4,9,3,15), address 1012954554 | 36 | 0 |

F1 here uses positive addresses (b,e), so it reads the words b−1,e−1. The [literal outgoing masks](frontier_witnesses.tsv) give the four F1 and one F2 fixtures, which can be printed with all outgoing masks and all arc loads (see commands below). These are upper bounds on first order, not exact first-order certificates. A simple independent lower bound is 16 for the two 2^32 targets and 17 for the 2^35 targets. Indeed partition a Hamilton order into consecutive blocks of sizes m_i. Counting ordered vertex-set partitions, then within-block paths, gives H≤n! ∏ M(m_i)/m_i!, where M is the maximum at that order. Tiny DP gives M(5)=15 and M(6)=45. At n=15 take 5+5+5 (upper bound 2,554,051,500); at n=16 take 5+5+6 (20,432,412,000). A smaller witness would pad to that order, so those sizes are excluded too.

## 2. What termination actually requires

Let N>1 be odd. Sort a tournament into its strong components, then drop all singleton components. Paths pass through the component order exactly once, so their count is a product over the remaining strong parts. [Busch's exact strong minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) gives H(C)≥3^(|C|/3) for each nontrivial strong component (check sizes 3,4,5, then multiply by 5 per additional three vertices). Therefore **every normalized witness has at most floor(3 log_3 N) vertices**. This finite decidability reduction appeared previously in the repo birth-fibre work and is implemented here, not claimed as new.

In selector.py, factor splits combine smaller witnesses by ordinal join. Strong candidates are scanned only while their exact Busch minimum permits the prescribed count, and only orders whose factorial can reach N. A fixed Hamilton path can be made the order 0,…,n−1 by labeling, so only the remaining chord orientations need enumeration. Minimum-order enumeration has a rigorous absence certificate at every smaller size, unlike finding an early witness in a single family. `--all-classes` deduplicates by minimum code over degree-preserving label orderings. The fallback returns a literal witness or global absence and always terminates (roughly 2^O((log N)^2) time with the elementary enumeration). What remains open is whether it always returns a witness on the permitted inputs, and whether it can be replaced by a fast universal constructive decoder.

construct.py first uses the verified frontier receipts when applicable; otherwise it tries the inherited exact F1 branch and then a small default F2 bank (or extra supplied rooted masks); `--exhaustive` enables the finite fallback. The fast branch uses a carry split rather than scanning all right words. Write the free address minus one as v=Lq+r, L a power of two, and g for the one-apex recurrence. Unroll the k low bits:

    g(Lq+r) = g(q)+s_r q+g(r)-1 + [high bit of r is zero][q odd],

where s_r is the k-bit reversal of the complement of r. This holds because appending a zero adds the current prefix plus its parity; only the very first appended bit reads parity of q. Holding r fixed bounds q from 1≤g(q)+epsilon[q odd]≤q+2. Then test the exact equality on the resulting short interval. This inverse is exact within its search domain; a bounded miss has no global implication.

## 3. All arc loads, including what their reversals would do

Let E(S,u) and B(S,v) be the count of paths on S ending at u and starting at v. For distinct u,v define a directional potential even if u→v is not present:

    W(u,v) = sum_(S contains u, excludes v)
                 E(S,u) B(V−S,v).

Each path using u→v splits at that arc; the recurrence therefore gives its actual load if present. Neither subpath can read the arc joining their disjoint vertex sets. Consequently every one-arc reversal count is available after the same two DPs:

    H(flip(u,v)) = H(T) + W(v,u) − W(u,v).

This is O(n²2^n) arithmetic operations **for all flips together**. An unused arc can be reversed without killing a path, either adding an even number of paths or doing nothing. The two directional potentials sum over all ordered nonempty two-covers; the actual loads sum to (n−1)H. Endpoint consistency gives sum_out(u)=H−end(u) and sum_in(u)=H−start(u). Intercomponent arcs have positive load exactly between consecutive strong components; within a component the local loads are multiplied by the other component counts.

There is a companion birth test: for deleting v, H(T)−H(T−v) is the sum of opposite arc loads in T−v over directed triangles v→x→y→v. Indeed along each old path there is one default insertion and one extra for each adjacent transition from a vertex beaten by v to one beating v. At a true first birth, every vertex has at least one such loaded opposite arc. The test restates the attachment law; passing it is not sufficient for globally minimal order.

The class audit inherits the complete induction of the [order-nine birth atlas](../chanpak_birth_atlas_metagraph_20261117/metagraph/PROOF.md), regenerates every first-born class through order 8, and computes both directions of every arc. At orders 1,3,4,5,6,7,8 it profiles respectively 1,1,1,6,33,328,5457 first-born classes. **Only one** has an unused arc: an order-five H=9 class, masks (10,28,25,16,1), whose arc 1→4 has zero load and reverses to H=13. All 1,520 single representatives from the inherited order-nine value table happen to be fully used; this does not say all order-nine classes are. The small generic load routine was independently checked against literal permutations and recomputed flips through order six.

## 4. Polynomial load scan for the large witnesses

join_loads.py does not allocate a 2^38 table. On each one-apex core, a path or marked cover assigns every spine vertex to one of two ordered subsequences. Scan increasing index, remembering the last vertex assigned to each. There are O(n²) states across all layers (one last vertex is always the current index). Each transition contributes a labeled internal arc or exposed endpoint event. Charge it by the number of partial prefixes times the number of valid completions. Joining the two cores uses the two exact cases: paths not using the reversed bridge (subtract the formerly available forbidden endpoint), and paths using it (two independent rooted covers). A small explicit left bank is profiled via permutations instead; the right side still uses the two-slot scan. The full load tables agree with generic subset DP on 256 F1 and 16 F2 controls.

**An all-width unused-arc law for the b=2 family.** Label the small left module a→x and the right spine with decreasing orientation on indices 0,…,n−1. Let y be the right apex, with y→i exactly on the 1 bits of v=e−1 in LSB indexing. Right beats left except the bridge x→y. For n≥2 the only possibly unused arc is (n−1)→x, and it is unused exactly when the second most significant bit is 0.

To see this, a spanning path using the bridge has one of the two forms P,a,x,S or P,x,S,a. Here P is optional, and S starts at y and descends in a right-spine subset. Its subset is empty or has largest index on a 1 bit. For arc i→x, P must end at i, hence every index below i is in S; if i<n−1 put P={i}, so the top bit supports S. For i=n−1, P can only consist of that vertex and the second bit decides. For n=1 this complement is empty, so the arc is used. Every spine arc j→i is used by placing the pair {i,j} alone in P if j is not top, in the spine part of S if it is top. Each y→i can be used by putting i alone in S, and each i→y is used in a no-bridge path with i on one side of y and the complement (including top) on the other side. All edges into a are used by putting their tail at the end of P or S, and a→x and x→y occur in the displayed forms. This isolates a literal structural unused edge in the first 2^32 7p witness, rather than merely reporting a histogram.

## Reproduction

```sh
python construct.py 34359738421
python construct.py 34359739333 --b-limit 1
python construct.py 9 --all-classes
python join_loads.py --check 2:1502296101 2:1847524150 9:2371737543 71:297178177 2,4,9,3,15:1012954554
g++ -O3 -std=c++17 first_classes.cpp -o /tmp/first-classes
/tmp/first-classes 8 --witnesses > /tmp/classes8.tsv
python atlas_audit.py /tmp/classes8.tsv | diff -u atlas_audit.out -
(cd targeted && python replay.py | diff -u replay.out -)
```

The absolute selection step is still missing. See [prime-bearing C7 rail](PRIME_RAIL.md) for a complementary construction that lifts some prime inputs exactly but cannot yet choose its width from an arbitrary prescribed prime.

One known topological constraint on unused arcs is relevant for future pruning: [Busch, Jacobson and Reid's arc-traceability work](https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.20173) implies that both endpoints of an unused arc of a strong tournament must be strong cut vertices. In the first 2^32 7p example above the unused arc's two endpoints are exactly of this kind; deleting the top spine vertex exposes a zero top bit, and deleting the left root destroys the sole inter-block return. The largest witnesses with no unused arcs can still have cut vertices, so this constraint is necessary, not sufficient.
