# Independent audit: forest score ceiling, cycle rank, and Paley11 distance

**Verdict: PASS.** The all-order theorem is proved below. The complete six-move
certificate independently verifies, so the legal Hamiltonian-path-complement
distance from Paley11 to the forest-framed class is **exactly six**, not merely
at least six. Relabellings and the full converse may be free.

Audit baseline: repository commit `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`.
The auditor did not edit the shared repository and did not import producer code.
Only the six-move plain-text witness was copied as input. All arithmetic and
adjacency operations below use independently written code.

## 1. Definitions and convention

For an n-vertex tournament, let d(v) be its outdegree and let

`L(T) = (1/2) sum_{i=0}^{n-1} |d_(i) - i|`,

where the scores are sorted increasingly. Equivalently, `2L` is the minimum
L1 assignment cost from the labelled degree vector to the targets `0,...,n-1`.
This is the elementary sorted matching inequality, obtained by uncrossing
assignments. The signed deviations sum to zero, so L is an integer.

Use the transitive baseline `j -> i` for `i < j`; its score at vertex i is i.
A backward-pair set D reverses those baseline edges. This is the reverse of
writing a forward order left-to-right, but it is the same undirected backward
pair graph after reversing the order. Let m be its edge count, c its number of
connected components **including isolated vertices**, and `b=m-n+c` its cycle
rank. Define intrinsic `beta(T)=min_order b(D_order(T))`.

## 2. All-order proof, n >= 2

At vertex i let p_i count reversed pairs from i to larger indices, and q_i
reversed pairs from smaller indices to i. Then

`d_i=i+p_i-q_i`, and `sum_i(p_i+q_i)=2m`.

Consequently `2L <= sum_i |p_i-q_i| <= 2m`.

**Connected strictness.** If D is connected, then `L <= m-1`. Suppose instead
that `L=m`. Equality throughout the preceding display requires that each
vertex has only one of p_i,q_i nonzero. Connectivity and n>=2 eliminate
isolates, so every deviation `Delta_i=p_i-q_i` is nonzero. Vertex 0 has positive
deviation and vertex n-1 has negative deviation. Thus some consecutive indices
i,i+1 have `Delta_i>0>Delta_(i+1)`. In particular `d_i>=i+1` and
`d_(i+1)<=i`. Swapping just their assigned targets i,i+1 lowers total L1 cost
by exactly two, contradicting minimal cost `2L=2m`. Integrality gives the
claimed `L<=m-1`. **Acyclicity was not used in this lemma.**

Now if c=1, connected strictness yields
`L<=m-1=n-2+b`. If c>=2, the elementary bound yields
`L<=m=n-c+b<=n-2+b`. Therefore, for every backward graph,

**`L(T) <= n-2+b(D)`; hence `beta(T) >= max(0,L(T)-n+2)`.**

The minimization is over b itself; this argument never assumes that an
f-minimizing order also minimizes b. Empty D has c=n, b=L=0 and is covered.
At n=1 the displayed n-2 bound fails, so the n>=2 hypothesis is essential.

A forest has b=0, yielding the proposed **`L(T)<=n-2`**. In particular,
`L(P11)=15` while every order-eleven forest frame has `L<=9`.

## 3. The lawful operation and its consequences

For a directed Hamiltonian path P, the legal complement reverses exactly arcs
off P. Relative to the full converse it restores path arcs. Only the two
endpoint degrees change: `d'_start=n-1-d_start+1`,
`d'_end=n-1-d_end-1`, and `d'_v=n-1-d_v` internally. Converse fixes L, while
the two endpoint changes alter the assignment minimum by at most two. Thus
`|L(tau_P T)-L(T)|<=1`. Relabelling also fixes L.

After k legal complements from P11, `L>=15-k`, and therefore

**`beta(T_k)>=max(0,6-k)`**, and **forest distance >=6**.

These statements hold on the fixed eleven-vertex carrier, including free
relabelling/converse. They do not concern vertex additions, deletions,
substitutions, or arbitrary arc edits.

## 4. Six-move upper certificate: independently verified

Start with vertices 0,...,10 and arc `u->v` exactly when
`v-u (mod 11) in {1,3,4,5,9}`. Every following path was independently checked
in its current literal tournament, followed by the exact off-path complement.
All resulting adjacency rows match `producer-greedy.pinned.out`.

| depth | legal path in the preceding tournament | L | f (independent subset DP) | frame m,c | intrinsic beta |
|---|---|---:|---:|---|---:|
| 1 | 0,1,2,3,4,5,6,7,8,9,10 | 14 | 14 | 14,2 | 5 |
| 2 | 1,3,9,0,10,5,2,8,4,6,7 | 13 | 14 | 14,1 | 4 |
| 3 | 0,3,2,1,6,5,9,8,4,7,10 | 12 | 12 | 12,2 | 3 |
| 4 | 1,6,2,4,7,8,3,5,9,10,0 | 11 | 12 | 12,1 | 2 |
| 5 | 0,4,3,5,8,6,9,10,2,7,1 | 10 | 11 | 11,1 | 1 |
| 6 | 1,8,7,3,10,6,5,0,4,2,9 | 9 | 10 | 10,1 | 0 |

The final **forward** order is `(1,10,4,6,8,7,5,2,3,9,0)`. Its backward graph
has ten edges and is connected and acyclic. Therefore the lower bound is
attained and the distance is **EXACTLY SIX**. Each prefix is shortest from
P11 to that state, and each suffix is shortest from that state to the forest
class, since their L differences attain the one-per-move bound.

Each listed intrinsic beta is exact: the displayed frame is an upper witness
and the new universal bound supplies the matching lower bound. In particular
the minimum beta over all one-step P11 images is exactly five, without a
one-step census. The independent subset-DP recurrence for f is

`F(S)=min_{v in S} [F(S\{v}) + |out(v) intersect (S\{v})|]`,

where v is placed last in a forward order.

## 5. Sharpness and hostile controls

* The n11 backward Hamiltonian tree with undirected vertex sequence
  `(4,7,0,8,1,9,2,10,3,5,6)` has sorted scores
  `(2,3,4,5,5,5,5,5,6,7,8)` and L=9. Thus the forest L ceiling is sharp at 11.
* A blanket even-order improvement to n-3 is false: at n8 reverse the forest
  edges `{0-7,1-6,1-7,2-5,2-7,3-6}`; L=6=n-2.
* `L<=m-1` for **every** m-edge forest is false: n3, single reversed edge
  `{0-2}`, has L=m=1. Connectedness across the full carrier is load-bearing.
* Replacing forest by merely `m<=n-1` is false: at n7 reverse
  `{0-4,0-5,0-6,1-5,1-6,2-6}`. All seven scores are three, so L=6>5. The graph
  has a four-cycle and an isolate.
* The final six-move target has f=10=n-1. Thus improving the **feedback-edge**
  forest ceiling to n-2 is false even though the **score** ceiling improves.
* For n=2r+1, r>=1, reverse precisely `(i,j)` with
  `0<=i<=r-1` and `i+r+1<=j<=2r`. All scores become r, the noncentral vertices
  form one connected component and vertex r is isolated. Thus
  `m=L=r(r+1)/2`, c=2, and
  `beta=L-n+2=(r-1)(r-2)/2` exactly. This is an all-order sharpness family for
  the cycle-rank lower bound. Its n11 member has beta6; it is not P11.

The inherited claim `beta(P11)=10` is in THM-4810's Paley census at the pinned
repository commit, not independently recomputed by this audit. Combined with
our independently verified first step beta5, it shows beta itself is not
1-Lipschitz. **The exact distance-six proof does not depend on beta(P11)=10,
f(P11)=20, or any old P11 census.**

## 6. Finite exact checks and reproduction

* C++: every backward forest on a fixed labelled order n=2..9, without orbit
  quotienting. This is 10,628,684 forest instances, including 10,026,505 at n9.
  Maximal L values are `0,1,1,3,3,5,6,7` for n=2..9. Every cap check passes.
* AddressSanitizer/UndefinedBehaviorSanitizer replay matches the C++ output;
  leak detection was disabled because the execution environment uses ptrace
  (the first LeakSanitizer exit reported that environment limitation).
* Independent Python raw-edge-mask / DFS check: every forest n=2..6, matching
  the full C++ L histogram and edge-count histogram. It brute-forces all target
  assignments for every forest n<=5, checking sorted matching directly.
* C++ forest totals n<=9 are independently matched by the component recurrence
  `F_n=sum_k binom(n-1,k-1) t_k F_(n-k)`, where `t_1=1` and
  `t_k=k^(k-2)` for k>=2; tree counts are independently checked against Cayley.
* Every arbitrary backward graph n=2..6 (33,866 graphs) satisfies both the
  universal beta inequality and connected strictness where applicable.
* Every legal Hamiltonian-path complement of every tournament n=2..5
  (7,886 directed moves) satisfies the endpoint law, involution and L bound.
* Sharp regular staircase examples are checked through r15 / n31.
* Six literal P11 moves and their supplied frame orders are checked directly;
  f is independently computed by 2^11-state DP. No complement-state BFS.
* Python normal and `-O` runs produce byte-identical JSON; all gates use explicit
  exceptions rather than assertions that optimization would erase.

Run `bash reproduce.sh` from this directory. See `SHA256SUMS` for hashes.
All claims of exact distance are theorem plus explicit certificate, not an
inference from the bounded forest enumeration.
