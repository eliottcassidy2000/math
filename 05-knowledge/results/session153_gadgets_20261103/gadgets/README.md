# Three-vertex attachment obstruction and a two-vertex repair

Status: finite, representative-specific statements, not a spectrum theorem and not a closed induction. We use H for the number of directed Hamiltonian paths. Vertex labels in an n-tournament are 0,...,n−1. Tournament code bits enumerate pairs i<j lexicographically, with bit 1 meaning i→j. An ear of mask m adds vertex x with x→i exactly when bit i of m is 1.

Let A be the order-six tournament formed from code 759 (order five) by ear mask 2. Let P=ear(A,50) and Q=ear(A,24). Direct literal permutation counting yields H(A)=31, H(P)=93, and H(Q)=95. All three are strongly connected. These are exactly the representatives in Session 152's ternary branch test.

**Finite obstruction.** No tournament formed from P by adjoining at most three arbitrary vertices, preserving all old arcs, has H=283. No tournament formed from Q by adjoining at most three arbitrary vertices has H=287 or H=289. In particular these representatives cannot realize their missing 3H+4, respectively 3H+2 or 3H+4, using at most three new vertices. The exclusions do not apply to other H=93 or H=95 representatives, or to larger gadgets.

**Positive two-vertex constructions.** Starting from P, ears 32,32 give a strong tournament with H=281=3·93+2; starting from Q, ears 126,126 give a strong tournament with H=285=3·95. The second mask omits the first added vertex, so the first added vertex points to the second. Each pair has identical old neighborhoods, and is a Session 150 rooted transitive-block gadget or its converse, not a new operation. Literal permutation enumeration independently verifies these positive counts and strongness. The result shows a two-vertex block can overcome a one-ear obstruction for these particular branches, while the preceding exclusions show it does not close the bank.

## Complete search and audit

The direct exhaustive program `independent.cpp` checks all additions of exactly three vertices, using recursive memoization of path completions, without the quadratic identity below. Among any three new vertices, relabel two so x→y; this loses nothing, since every tournament has an arc. Their old-neighborhood masks u and v range independently over 0,...,127, with x→y fixed. The third vertex has an arbitrary 9-bit mask w, encoding all nine arcs to the old seven vertices and x,y. Therefore every three-vertex extension occurs. The code discards a partial extension only when H already exceeds the target; adding a vertex cannot decrease H: canonically insert it into each old path immediately before the first old vertex it beats, or at the end if it beats none. Deleting the inserted vertex recovers the path, so this is injective.

For each target the direct audit prints the surviving (u,v) pair count and checks every w, with results:

| Parent | Target | pairs | 9-bit masks checked | hits |
|---|---:|---:|---:|---:|
| P | 283 | 70 | 35,840 | 0 |
| Q | 287 | 58 | 29,696 | 0 |
| Q | 289 | 58 | 29,696 | 0 |

Adding a new source to any tournament preserves H, since it must be the first vertex in any Hamiltonian path. Thus an extension by zero, one, or two vertices could be padded with sources to exactly three, establishing the stated “at most three” result.

The second implementation `three.cpp` computes the third-ear response using the known exact quadratic dependence on the ear mask, reconstructing its coefficients from the empty mask, singletons and pairs; it independently produces the same zero-hit pair counts. `audit.cpp` compares that quadratic against every 9-bit mask for five surviving two-ear prefixes for each parent. `independent.cpp` directly tests all 95,232 final masks, not a sample; its recursive completion DP is distinct from the forward subset-last DP in `helper.hpp`. The 32-bit integer counts are safe here: each considered tournament has at most 10 vertices and thus at most 10! paths.

## Reproduction

Run from this directory:

```
g++ -O3 -std=c++17 independent.cpp -o /tmp/independent153
/tmp/independent153
g++ -O0 -std=c++17 independent.cpp -o /tmp/independent153_O0
/tmp/independent153_O0
g++ -O3 -std=c++17 three.cpp -o /tmp/three153
/tmp/three153
g++ -O0 -std=c++17 three.cpp -o /tmp/three153_O0
/tmp/three153_O0
g++ -O3 -std=c++17 audit.cpp -o /tmp/audit153
/tmp/audit153
python3 verify_witness.py
python3 -O verify_witness.py
```

Frozen outputs are `independent.out`, `independent_O0.out`, `three_O3.out`, `three_O0.out`, `audit.out`, `verify_witness.out`, and `verify_witness.opt.out`. On the session machine the direct `-O3` check took 9.590 seconds; the quadratic check took 2.817 seconds at `-O3` and 5.880 seconds at `-O0`; the literal witness check took 0.465 seconds. Times are indicative, not performance guarantees. Optimization modes produced identical text in each comparison. The unoptimized direct check took 14.074 seconds.
