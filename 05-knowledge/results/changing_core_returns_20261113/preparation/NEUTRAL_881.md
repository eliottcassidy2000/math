# Smallest arc repairs at 881, and a smaller strong core

**PROVED from exact finite certificates:** the specified minimum patch widths and all-branch arithmetic guards. **FINITE-EXACT:** the literal source, patch, component, and return counts. These statements use a fixed named root and distinguish one-arc neutral paths, simultaneous endpoint patches, and subsequent tail deletion. The desired repeatable unit-ready return and Chan–Pak remain OPEN.

## N1. The full one-arc neutral components are too small

Take the two strong unit-ready H881 sources in the [selector receipt](../selector/receipt.json): A has order10 and root 9, from parent(order9,rank168); B has order11 and root 10, from parent(order10,rank174). Their original states are

```
A: (881,74,75,882),       B: (881,218,219,882).
```

The independent `neutral_881_patch.py` checks every one-arc neighbour and traverses the entire H881-neutral component. For A the only source neighbours reverse {4,5} or {7,8}; for B the only source neighbour reverses {5,6}. Every such pair is homogeneous, avoids the named root, and is an adjacent pair in a transitive two-vertex module. The complete components contain respectively four and two labelled graphs, with constant full marked profiles.

The general [neutral-component theorem](../geometry/NEUTRAL_SCC.md) explains why this finite local audit closes every possible route: homogeneous swaps are relabellings, the neighbour property is equivariant, and no nonhomogeneous neutral edge can first appear later. Its SCC factorization also extends the obstruction to arbitrary transitive tails. Thus no sequence of H-neutral single-arc flips at these fixed vertex sets can perform either repair below.

For the selector's particular five-arc role exchanges, the independent checker also enumerates all 32 partial patches. Only the empty and full patch have H881. Within those two five-dimensional cubes, all five arcs must change simultaneously to stay neutral at every macro endpoint. The minimum possible maximum |H-881| along single-arc paths is332 for A and314 for B; monotone paths already achieve these optima:

```
A: 881,1033,549,1103,771,881;
B: 881,1015,687,1195,891,881.
```

Those excursion bounds concern only the prescribed patch cubes. The full-component obstruction above has no such restriction.

## N2. A genuinely different core repairs A with only two arcs

In A, reverse {0,2} and {3,9}, keeping physical root 9. The prepared rows and state are

```
(518,388,456,1011,487,455,387,257,1,502),
(h,s,q,F)=(881,114,141,1066).
```

Both endpoints are strong. Their sorted outdegree lists differ:

```
before: 1,2,2,3,4,5,6,7,7,8;
after:  1,2,3,3,4,4,6,7,7,8.
```

Hence this changes the unrooted strong core; it is not an isomorphic role exchange. The prepared root has indegree 2 and q-s=27. The operation preserves the current full target H881, while losing unit readiness.

Here d=F-q=925 and Delta=F+s-h-q=158. The inherited physical return gives

```
t_k=925*2^k-17,       b_k=2775*2^k+265,       k>=1.
```

Since t_k is odd, gcd(b_k,t_k)=gcd(t_k,316) can fail only at 79. Modulo 79 the obstruction is 2^k=13. The order of 2 is 39, with 2^28=13, so exactly k=28 modulo 39 fails. The native divisibility guard requires k odd. Thus the complete successful set is

```
k odd, except k=67 modulo 78.
```

In particular actual RNR returns (5815,1833,1833,5815). The literal graph, physical new root, and independently counted full state are in the receipt. The final R creates the required root-child interface, and the inherited balanced controller is then repeatable.

Every pair of arcs from A was checked: 990 patches, exactly three H881 endpoints, exactly one with a nonempty native-entry wheel. Combined with N1, two is the minimum simultaneous arc width for this fixed named source to reach the declared R^kNR native-entry guard at unchanged H. This minimum does not impose preservation of a future target.

## N3. A five-arc repair preserves more than the two-arc repair

The incoming [continuation-fibre theorem](../../h17_exact_repairs_20261112/history/CONTINUATION_FIBRES.md#2-four-literal-count-observations-recover-the-rooted-profile) proves that preserving both h and the first birth count h+F-s requires exactly a difference (0,v,w,v). The selector's original A role exchange has

```
(881,74,75,882) -> (881,169,195,977),
difference = (0,95,120,95).
```

It is therefore an actual repair inside that stronger two-observation fibre. Its larger patch is not superseded by the two-arc result:

| A endpoint | Arc width | H | First birth H+F-s | Delta |
|---|---:|---:|---:|---:|
| Original source | 0 | 881 | 1689 | 0 |
| Role exchange | 5 | 881 | 1689 | 70 |
| New changing core | 2 | 881 | 1833 | 158 |

The five-arc preparation preserves the first birth target 1689 and enables the later primitive balanced entry (5207,1689). The two-arc preparation minimizes the H-only repair width and changes that first target. Neither yields the desired ready1691 source.

## N4. B needs four arcs if a tail may be removed, five for a strong endpoint of the same order

All 1,485 two-arc patches of B have only two H881 endpoints; neither meets the native-entry guard for any branch. All 26,235 three-arc patches have 19 H881 endpoints and again no successful branch. Every one of those 19 endpoints has root indegree at most 1. Normal and optimized Python scans agree exactly.

A native-entry repair must have root indegree at least 2. If the root is a source, F=s=q=h, hence Delta=0. If it has indegree 1, the original graph is an N birth at its unique in-neighbour from the graph obtained by deleting the root; the N map also gives Delta=0. In either case b_k=3t_k, excluding the guard. This is a proved filter, not a sampled heuristic.

Among all 341,055 four-arc patches,168,960 pass this necessary indegree filter. An independent C++ subset counter finds 16 H881 endpoints. All 16 full marked states are independently recounted in Python. Four have nonempty exact arithmetic wheels, but each is a strong rooted core followed by a nonempty transitive tail. The only strong endpoint in the filtered list has state(881,97,97,941), Delta=60, and is blocked. The C++ implementation also reproduces the complete filtered three-arc list of the independent Python scan.

Therefore:

- Four is the minimum simultaneous arc width if weak H881 endpoints are allowed.
- Five is the minimum width when the prepared endpoint must remain strong on the original eleven vertices: widths<=4 are excluded, while the selector supplies a five-arc strong repair.

There is a useful four-arc **changing-order** repair. Reverse

```
{0,9}, {2,6}, {2,10}, {3,6}.
```

The resulting H881 graph is a strong nine-vertex core followed by the transitive tail with old names 8,9. Delete those two tail vertices. Each tail deletion preserves all four rooted coordinates by the [tail lemma](../geometry/NEUTRAL_SCC.md#5-tail-padding-preserves-the-full-rooted-state-and-commutes-with-the-repair). Retain the old names 0..7,10 in that order; the same physical root 10 now has index 8. The resulting strong core is

```
rows=(354,228,425,179,231,192,140,1,250),
state=(881,169,201,977),   root indegree 2,   q-s=32.
```

Thus both the original strong eleven-vertex source and the final strong nine-vertex source have full H881. The simultaneous four-arc patch is endpoint-neutral; the following two individual tail deletions are each neutral. There is no claim that the four arcs can be reversed separately while preserving H. The pair {0,9} disappears after deletion; describing only the final surviving adjacency edits would use three flips plus two deletions, but would no longer certify the same intermediate H881 stages. This does not change the fixed-order minimum statements.

Here Delta=64 and d=776, so the all-branch law is especially simple:

```
t_k=776*2^k+137,       b_k=2328*2^k+539.
```

Both are odd; gcd(b_k,t_k)=gcd(t_k,128)=1 for every positive k. Modulo3, t_k=2(2^k+1), so exactly the odd k meet the native guard, and b_k=2 modulo3. The literal k1 exit from the nine-vertex core is

```
(5195,1689,1689,5195),
```

again strong and physically balanced. Its first birth 1689 agrees with A's five-arc role repair, while Delta=64 rather than 70 changes the final b by 12. Equal current and first-next targets still leave a meaningful continuation coordinate.

## N5. Reproduction and exact limits

Run from the repository root:

```
python -B 05-knowledge/results/changing_core_returns_20261113/preparation/neutral_881_patch.py
python -O -B 05-knowledge/results/changing_core_returns_20261113/preparation/neutral_881_patch.py
python -B 05-knowledge/results/changing_core_returns_20261113/preparation/neutral_881_three.py
python -O -B 05-knowledge/results/changing_core_returns_20261113/preparation/neutral_881_three.py
python -B 05-knowledge/results/changing_core_returns_20261113/preparation/neutral_881_four.py
python -O -B 05-knowledge/results/changing_core_returns_20261113/preparation/neutral_881_four.py
```

The last driver compiles its checked-in C++ source to a temporary directory, records all filtered candidates, compares its three-arc universe with the independent exhaustive Python scan, and recounts every four-arc candidate and selected literal return. No executable is checked in. The patch checker independently validates 285 complete rooted profiles through order4 by permutations and cuts. Source labels, patch universes, inherited wheels, and allowed operations are all explicit.

These are exact repairs of declared continuation capabilities at H881. They do not preserve unit readiness, construct881->1691, supply the ready405*2^k+71 family, or prove arbitrary-target coverage. The bounded neutral components explain why genuine changes of core, simultaneous patches, or a change of order are necessary within these menus. The reproducible `next_1691_probe.py` scans 44 recorded H881 endpoints (23 distinct marked profiles) and finds no first-birth target 1691; nearest is1689. This is a statement about those declared endpoint universes, not an all-patch exclusion.
