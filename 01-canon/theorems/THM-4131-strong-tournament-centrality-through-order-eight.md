---
id: THM-4131
title: "Strong-tournament Johnson centrality through order eight"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every strong
  tournament of orders four through eight has only central maximizers of the
  THM-4128 rational Johnson-slice support floor and only central maximizers
  after THM-4123 exact layer-coset rounding. The central-versus-outer coset
  margin is uniformly positive in every audited order. Actual ear-response
  maxima need not be central: four order-six and 1,702 order-eight strong
  isomorphism classes are explicit finite boundaries. THM-4135 later extends
  support-floor centrality through order nine, while THM-4133 refutes its
  all-order extension at order twelve. No all-order or central actual-
  maximizer theorem is claimed.
source: codex-frontier-synthesis-creative-20260825m
depends_on:
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4127-johnson-slice-hoeffding-variance-and-central-support-dominance
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
related:
  - THM-4118-ear-response-lattice-and-stateful-unit-component-intervals
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4135-strong-tournament-centrality-complete-order-nine
  - HYP-2879-strong-ear-atom-calculus
  - HYP-9029-strong-interval-tiling-law
script: 04-computation/tournament_strong_centrality_through_order_eight_thm4131.py
output: 05-knowledge/results/tournament_strong_centrality_through_order_eight_thm4131.out
independent_audit_script: 04-computation/tournament_strong_centrality_through_order_eight_thm4131_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_strong_centrality_through_order_eight_thm4131_independent_audit.out
script_sha256: 6b195b0379d1ae3e5d215aa1c495f7180daeecae189df86269d07ef855867881
output_sha256: 1e6f3258fd740a329451f214302d56eb6754795ee211dcac39df4a4eeb96157e
semantic_sha256: 8e5ef116e577c3d3ab5dd4bea50953581a457ffbaae4f80504060c8793fed578
independent_audit_script_sha256: 693484165663f5e05d35fb7b6dfd9ae0cc01a49c72f2db74db74017d93ab5b94
independent_audit_output_sha256: 6efb87263c698375a1ce5da57251434ad3f71fe01060db96fa54d19b8a340542
independent_semantic_sha256: 8e5ef116e577c3d3ab5dd4bea50953581a457ffbaae4f80504060c8793fed578
hash_basis: raw LF bytes
primary_audit: >
  PASS. Homebrew nauty 2.9.3 gentourng -q -c supplies one strong
  isomorphism-class representative in each order. A standalone
  Start/End/exposed-gap subset DP reconstructs every ear response, support
  floor, layer lattice, and actual layer maximum for all 6,403 classes.
  The five raw generator streams and five full profile multisets are
  fingerprinted. Normal, optimized, hash-seeded, and frozen replays agree.
independent_audit: >
  ACCEPT. A pure-standard-library deletion-closed canonical augmentation
  recovers all tournament class counts 4,12,56,456,6880 and independently
  filters strong counts 1,6,35,353,6008. Exposure uses a contracted-adjacent-
  block Hamilton DP, adding each good arc block to its reversed one-defect
  block, rather than the primary exposed-gap carrier. All order rows,
  controls, and profile fingerprints agree with the primary semantic ledger.
  A third no-edit generator referee independently reproduced the census and
  accepted the theorem. Normal, optimized, hash-seeded, and frozen replays
  agree.
---

# THM-4131 -- strong-tournament Johnson centrality through order eight

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4128 leaves strong centrality open because its exact criterion is a
constraint on the directed exposure packet, not a consequence of strongness.
THM-4123 leaves a second gap: layer-specific response lattices can reorder the
rational floors. The complete strong isomorphism-class census through order
eight closes both finite questions. It deliberately retains THM-4123's
warning that a support floor need not locate an actual maximizing ear.

## 1. Statement and coordinates

Let `T` be a strong tournament of order `4<=n<=8`, and let

```text
F_T(S)=H(T+x_S),                         S subseteq V(T),       (1)
```

with the ear convention of THM-4123. For `1<=m<=n-1`, put

```text
J_m=E_m F_T+Var_m(F_T)/(E_m F_T-H),                         (2)
L_m=ceil_(a_(T,m)+d_(T,m)Z)(J_m),                           (3)
t=n-2m.                                                      (4)
```

Here `J_m` is THM-4127's rational slice-support floor, while `(a,d)` is
THM-4123's exact Johnson-layer coset. Call `t=0` central for even `n`, and
`t=+-1` central for odd `n`.

> **Finite centrality theorem.** For every strong tournament of orders four
> through eight:
>
> 1. every layer maximizing `J_m` is central;
> 2. every layer maximizing `L_m` is central; and
> 3. the best central `L_m` is strictly larger than the best noncentral
>    `L_m`.

The third assertion is stronger than existence of a central optimizer, but
only in this finite universe.

For the rational part, THM-4128 gives the vertex and normalized centrality
load

```text
theta=(n-3)C_hd/(2D_4),                                      (5)

rho(T)= (n-3)|C_hd|/(2D_4)                 n even,
        (n-3)|C_hd|/(4D_4)                 n odd.             (6)
```

Central rational optimality is equivalent to `rho<=1`. The census proves the
strict inequalities below; the exact-coset assertion then needs the separate
layer-by-layer computation `(3)`.

## 2. Complete order rows

The all-class counts, strong-class counts, rational/coset failures, actual-
maximum failures, zero-tilt counts, worst normalized tilts, and minimum
central-versus-outer coset margins are:

| `n` | all / strong classes | rational / coset failures | actual failures | `theta=0` | `max rho` | `min(best central L-best outer L)` |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | `4 / 1` | `0 / 0` | `0` | `1` | `0` | `2` |
| 5 | `12 / 6` | `0 / 0` | `0` | `6` | `0` | `6` |
| 6 | `56 / 35` | `0 / 0` | `4` | `11` | `198/901` | `6` |
| 7 | `456 / 353` | `0 / 0` | `0` | `79` | `711/3959` | `22` |
| 8 | `6880 / 6008` | `0 / 0` | `1702` | `162` | `14809/29733` | `16` |

Thus `6,403` strong isomorphism classes are audited. The rational optimizer
histograms, written as `optimizer-tuple:class-count`, are

```text
n=4: (0):1
n=5: (-1,1):6
n=6: (0):35
n=7: (-1):137, (-1,1):79, (1):137
n=8: (0):6008.                                               (7)
```

The exact-coset histograms agree at even orders and are

```text
n=5: (-1,1):6,
n=7: (-1):100, (-1,1):153, (1):100.                         (8)
```

At order seven, `74` classes have their rational optimizer tuple changed by
coset rounding, but every new optimizer remains central. This is precisely
why rational centrality alone does not prove the second assertion.

The worst-tilt classes occur in reversal pairs. Representatives have packets

```text
n=6: (H,W,D_4,C_hd)=(27,117,2703,-396), theta=-198/901;
n=7: (H,W,D_4,C_hd)=(45,321,23754,-4266), theta=-1422/3959;
n=8: (H,W,D_4,C_hd)=(111,1077,297330,-59236),
     theta=-14809/29733.                                    (9)
```

Even the order-eight extremum uses less than half of the allowed rational
centrality budget.

## 3. Actual maxima are a different coordinate

Let `A_m=max_(|S|=m)F_T(S)`. The exact actual-optimizer histograms are

```text
n=4: (0):1
n=5: (-1,1):6
n=6: (-2):2, (-2,0):4, (0):23, (0,2):4, (2):2
n=7: (-1):128, (-1,1):97, (1):128
n=8: (-4):1, (-4,-2):1, (-2):849, (-2,0):68, (0):4170,
     (0,2):68, (2):849, (2,4):1, (4):1.                  (10)
```

Accordingly, four order-six and `1,702` order-eight classes have no central
actual maximizer. Strong labelled code `20` is the inherited minimal hostile:

```text
argmax J_m=argmax L_m={t=0},
J(0)=10399/109>J(2)=9671/109,
max_(t=0)F_T=131<133=max_(t=2)F_T.                        (11)
```

This does not contradict the theorem: `(2)` and `(3)` are certified lower
support floors, not exact maxima.

Two nonstrong controls show that strongness and the coset sidecar are both
load-bearing. Labelled code `2` has

```text
rho=theta=141/125>1,             argmax J_m={t=2},        (12)
```

while labelled code `140` has

```text
argmax J_m={t=2},                argmax L_m={t=4};         (13)
```

its `t=4` layer has lattice `102`, exact floor `135`, and actual maximum
`135`. These are the THM-4128 hostile mechanisms, not merely random negative
tests.

## 4. Primary exhaustive universe

The primary audit invokes

```text
gentourng -q -c n,                         4<=n<=8,        (14)
```

using Homebrew nauty `2.9.3`. The environment variable
`THM4131_GENTOURNG` may override the executable path. The audited binary has
SHA-256

```text
89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110. (15)
```

The five raw strong-stream fingerprints are

```text
n=4  76650554f5ab120115e47364bbe1822257753e89b72c0731b69edd77c0cd9404
n=5  61a69b4844a4f1ec611b88250fb68654626be205c3c9245e30ad341eca745972
n=6  3d9ab51665d390c367ced4940d6b80a1233ac8d1a8e158d07134f9abd7bb9ab2
n=7  9b96ef048acfddb3b5b1ea29a3964f0987052cdf140f9dc1e150cbcd255c21bf
n=8  6900758e8f64444dd2a75450f35e05a1f1f2e00bdf6a77f1395fecd26d689e5a. (16)
```

For each representative, a subset DP computes the number of Hamilton paths
on every vertex set with each possible start and end. Splitting an old word
at every exposed gap reconstructs every value `F_T(S)`. The audit then
computes `(2)`, the THM-4128 nearest-grid optimizer, every layer gcd and
coset ceiling `(3)`, and every actual layer maximum. It also independently
checks strong reachability despite the generator flag.

## 5. Pure-standard-library independent universe and exposure

The independent audit begins with the one-vertex tournament. From every
canonical order-`n-1` representative it attaches a new vertex in all
`2^(n-1)` possible orientations. Its canonical certificate first applies
invariant directed color refinement, then takes the least upper-triangle word
over every permutation inside the stable color cells. Isomorphisms preserve
those cells, so certificate equality is equivalent to isomorphism. Canonical
deduplication returns exactly

```text
4,12,56,456,6880                                             (17)
```

classes in orders four through eight. An internal forward/reverse
reachability filter returns exactly

```text
1,6,35,353,6008                                              (18)
```

strong classes. This path uses only the Python standard library and is the
replay fallback when nauty is unavailable.

It also avoids the primary exposed-gap carrier. Let `E(A,p)` count Hamilton
paths on `A` ending at `p`, and `S(B,q)` count those on `B` starting at `q`.
For distinct `u,v`, put `R=V minus {u,v}` and

```text
B(u,v)=sum_(L subseteq R)
 [1_(L=empty)+sum_(p in L,p->u)E(L,p)]
 [1_(R\L=empty)+sum_(q in R\L,v->q)S(R\L,q)].             (19)
```

The empty-side indicators in `(19)` replace the corresponding sums rather
than add to nonempty sides. If `u->v`, then the directed exposure capacity is

```text
c_uv=B(u,v)+B(v,u).                                        (20)
```

The first summand contracts the good adjacent block `(u,v)`; the second
contracts the reversed one-defect block `(v,u)`. Every remaining adjacency
is good, so this partitions exactly the exposed zero/one-defect words. Hence

```text
F_T(S)=H+sum_(u in S,v notin S,u->v)c_uv.                  (21)
```

The independent implementation checks `(21)` against literal child
Hamilton-path DPs at codes `2`, `140`, and `20`. Its five invariant profile
fingerprints agree byte-for-byte with the primary values:

```text
n=4 d04c46f880d2b04a05ae4be53e5d68d825a117b3ac92eb23d439242cd2cd1101
n=5 13fb1174d8dabac1d4ef351fc69794ef4184a2935ddf0b672a55aac0ee55f1c6
n=6 250726d6a4b577c5514f3e6578a38da1a72b3a3755ebbf3b5ede6ee811bea28e
n=7 2aeca97db51de7dcd88ece3b3aba96289f2a1e0208119c6deae5855dd5814524
n=8 a07ed552eb832d7e84bf2e822311c8b722ca48927a797f1e6a141757f03dad31. (22)
```

A third no-edit generator referee used a separate canonical convention. Its
combined strong-universe fingerprint was
`3faa5771d55b15faf2405c22d6ef6f7792952ebd67ba239ee8e54116e6af8e4d`,
and it independently reproduced `(17)--(18)`, every failure row, the five
worst tilts, zero-tilt counts, and strict coset margins.

## 6. Replay and scope

Run

```text
python3 -B 04-computation/tournament_strong_centrality_through_order_eight_thm4131.py
python3 -B -O 04-computation/tournament_strong_centrality_through_order_eight_thm4131.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_strong_centrality_through_order_eight_thm4131.py
python3 -B 04-computation/tournament_strong_centrality_through_order_eight_thm4131_independent_audit.py
python3 -B -O 04-computation/tournament_strong_centrality_through_order_eight_thm4131_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_strong_centrality_through_order_eight_thm4131_independent_audit.py
```

All six streams match their corresponding frozen outputs. The shared semantic
digest is
`8e5ef116e577c3d3ab5dd4bea50953581a457ffbaae4f80504060c8793fed578`.

The source is the full strong-tournament isomorphism-class universe through
order eight; the preserved coordinates are the complete response vector,
the THM-4128 exposure packet, and every THM-4123 layer coset. Passing to a
support floor still loses the actual maximizer, as `(10)--(11)` demonstrate.
This theorem itself proves no strong centrality statement in order nine or
higher. THM-4135 later closes order nine and THM-4133 refutes the all-order
extension at order twelve; orders ten and eleven remain unclassified. No
central actual-maximizer theorem, slice interval, response connectivity, or
global `H`-spectrum result follows. **QED.**
