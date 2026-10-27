# C2 independent audit: exact Paley11-to-transitive distance fifteen

**Verdict: PASS.** There is a literal legal fifteen-move Hamiltonian-path-
complement word from Paley11 to a transitive tournament. Its score displacement
L decreases by one at every step, from fifteen to zero. Therefore its length is
optimal, even when arbitrary relabellings and full converse are free.

This is a separate C2 supplement. The C1 audit directory and its frozen bytes
were not modified. This audit uses a new adjacency-matrix implementation;
producer programs and earlier audit routines are not imported.

## 1. Exact statement and proof

The fixed carrier is `{0,...,10}`. Initially `u -> v` iff
`v-u (mod 11) in {1,3,4,5,9}`. For a directed Hamiltonian path P, the legal
operation keeps all arcs of P and reverses every other arc.

For sorted outdegrees `d_(0)<=...<=d_(n-1)`, define

`L(T) = (1/2) sum_i |d_(i)-i|`.

The transitive class has L=0, while P11 has every outdegree five and L=15.
The sorted matching inequality identifies `2L` with the minimum L1 assignment
cost to `0,...,n-1`. Full converse preserves it. Relative to the full converse,
a legal complement changes only the two endpoint degrees, by +1 and -1.
Thus a legal complement changes L by at most one. Relabelling also preserves L.
Every word from P11 to any transitive tournament therefore has at least fifteen
paid moves, with or without these free operations.

The independently checked fifteen-step word below attains that bound. Its
final forward order is

**`(6,0,2,9,7,4,8,5,10,3,1)`**.

All 55 pairs point forward in that order. Consequently the exact distance to
the transitive class is **15**. Every prefix is shortest to its literal state;
every suffix is shortest from its starting state to the transitive class.
The operation is an involution along the same retained path, so the reverse
word is legal too.

## 2. Certificate assembly and independent checks

The word is exactly:

1. All six paths in the frozen C1 certificate `C1-six-moves.pinned.out`.
2. The **first six paths only** in `extend-greedy.pinned.out`.
3. All three paths in `backtrack.pinned.out`.

The final two paths of the full greedy attempt are deliberately discarded.
They lead to a verified one-step descent trap, not to the final certificate.

The nine extension paths are:

| depth | directed Hamiltonian path in the preceding literal state | resulting L |
|---|---|---:|
| 7 | 2,0,4,6,5,9,1,7,3,10,8 | 8 |
| 8 | 3,1,7,4,5,8,6,10,2,0,9 | 7 |
| 9 | 2,0,9,5,3,1,6,7,4,8,10 | 6 |
| 10 | 3,2,0,4,8,5,7,9,10,1,6 | 5 |
| 11 | 2,0,4,5,7,9,6,10,1,8,3 | 4 |
| 12 | 3,1,8,2,0,4,5,7,10,9,6 | 3 |
| 13 | 6,2,9,7,10,3,1,8,0,4,5 | 2 |
| 14 | 5,10,3,1,8,0,2,9,7,4,6 | 1 |
| 15 | 6,0,2,9,7,4,8,5,10,3,1 | 0 |

For all fifteen steps, the verifier checks:

- The supplied path is a permutation of all eleven vertices.
- Every consecutive path pair is an arc of the current tournament.
- The output is formed by keeping exactly path arcs and reversing all others.
- Every pair remains oriented exactly once.
- All eleven resulting adjacency rows match the pinned plain-text certificate.
- The sorted-score displacement is exactly `15-depth`.

Final adjacency row bitmasks, where bit v of row u indicates `u -> v`, are

`[1982,0,1978,2,1322,1034,1983,1338,1066,1466,10]`.

`normal.json` contains every full path, adjacency vector, score vector, sorted
score vector, L value and exact prefix/suffix distance.

The depth-six C1 endpoint is independently rechecked to have a connected
10-edge backward graph in forward order `(1,10,4,6,8,7,5,2,3,9,0)`. The new
word therefore simultaneously realizes the established exact forest distance
six and extends that forest endpoint to the transitive class in exactly nine
further moves. Exact forest distance uses the separately proved C1 ceiling
`L<=n-2` for forest frames; the distance-fifteen result does not require that
new ceiling.

## 3. Failed-greedy hostile: verified, not overinterpreted

From the common depth-twelve state with L=3, the discarded greedy branch uses

- `(2,9,6,0,4,5,7,10,3,1,8)`, reaching L=2;
- `(5,7,10,3,1,8,4,6,0,2,9)`, reaching L=1.

Both paths and every output adjacency row are independently checked. The L1
state has exactly **243 directed Hamiltonian paths**, counted by two separate
methods: explicit path DFS and an independent subset-DP recurrence. Exhausting
all 243 legal complements gives successor-L histogram

**`{1:8, 2:235}`**.

There is no L0 successor and hence no one-step transitive successor. Thus
strict L descent is not always available, and L is not universally equal to
the distance to the transitive class. This does **not** prove that the state
is unreachable, nor does it compute its exact transitive distance.

## 4. Reproduction and scope

Run `bash reproduce.sh`. The normal and optimized (`python3 -O`) outputs must
be byte-identical. Every validation gate raises explicitly; none relies on a
Python assertion removed by optimization. `SHA256SUMS` records source and data
hashes. All producer input is pinned plaintext, not executable code.

No complement-state BFS, broad tournament census, mutation of C1, repository
edit, or claim about vertex-changing operations is involved. The proof and
certificate suffice for exact distance fifteen; historical Hamiltonian-path
counts and P11 feedback/cycle-rank censuses are not dependencies.
