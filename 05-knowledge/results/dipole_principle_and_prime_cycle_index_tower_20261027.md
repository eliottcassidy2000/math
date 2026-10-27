# The dipole principle, cut or colour, and the prime cycle-index tower

2026-10-27, mac-mini-2026-10-27-dipole (worktree `codex/session-pyth-cuboids`), owner directive: "pursue the new
reflection [which] argues that every distance lower bound in the repo comes from one identity: a move negates the
tournament's matrix and adds the path. Where the bound fails to be exact, a forced symmetry explains why. also the
Order-11 prime cycle-index types and the order-12 ledger". Canonical records: THM-4815 (dipole principle), THM-4816
(prime cycle-index tower; order-11, 12, 13 types), THM-4817 (order-12, 13, 14 ledgers), HYP-9423; MISTAKE-682.
Not independently audited at filing.

| Part | Content | Status |
|---|---|---|
| 1. One identity (THM-4815 A) | `K -> -K + 2A_Q` read through converse-invariant functionals gives every native-move lower bound in the repository (survey of all THM/HYP/result packets); so do its square (two-move XOR law), its reduction mod 2 (forward-arc parity) and its endpoint form (non-adjacency, module shattering, source/sink trap) | PROVED; survey |
| 2. Exactness (THM-4815 B1-B6) | for score and cut functionals only the endpoints of the path matter; a state descends iff the Hamiltonian endpoint relation `Pi(T)` meets the descending pairs; non-strong `T` has `Pi = C_1 x C_k`; the `delta0` local minima are exactly the non-strong tournaments without source or sink; the non-strong `L` local minima are exactly those with a source or sink; a colour obstruction forces an anti-automorphism | PROVED |
| 3. Max-to-min paths (THM-4815 B7, B8) | every tournament has a Hamiltonian path from a maximum-score to a minimum-score vertex; hence at odd orders the distance to the regular class is exactly `R = (1/2) sum |s_v - (n-1)/2|` | PROVED |
| 4. P_7 (THM-4815 C) | forward-arc parity alternates under every move at `n ≡ 0, 3 (mod 4)`; `P_7` is 6 = L moves from the transitive class (explicit word); its 7 is the labelled distance to one fixed transitive labelling, forced by parity | PROVED; FINITE-EXACT |
| 5. Census (THM-4815 D) | class graph orders 4-8 (9 in progress at filing), endpoint relation orders 4-10 (9,733,056 classes), random tests orders 11-16 | FINITE-EXACT |
| 6. Prime tower (THM-4816) | rigid primes from the ledger, non-rigid primes by invariant enumeration over the cycle types the cycle lemma allows; order-11 types (12 rows, 12,155 non-rigid), order-12 (14 rows, 285,705), order-13 (22 rows, 11,256,145); reproduces klein's order-3..10 census row for row with no gentourng input; labelled-count identity at every order 3-13 | PROVED (method); FINITE-EXACT |
| 7. Ledgers (THM-4817) | `P(12) = 142,980,435,791`, `SCP(12) = 1,170,949`, `c(12) = 76,601,960,166`, `P_lab(12)`, `l(12)`; the order-13 and order-14 ledgers; strong and labelled-strong counts agree with OEIS A051337 / A054946; order-14 defect `3,574,801,896/7` pre-registered | FINITE-EXACT |

**Credits.** The identity is THM-4811 (A) (this machine's previous session) and Session75 section 6 (`K -> -K + 2L`).
Theorem L, the parity theorem and the labelled distances are THM-4810 / S63. The source/sink trap and the `P_11`
certificate are Session74. The cycle lemma, anti-automorphism lemma, self-converse identity and order-11 ledger are
klein (THM-4814). The species identity is session56 / THM-4803. Lu (Graphs Combin. 34 (2018) 613-618) studies
Hamiltonian paths between min-max pairs (abstract only consulted). The survey of distance bounds was done by a
read-only agent; its flag that `P_7`'s 7 is a labelling artefact matched this session's independent BFS.

## 1. One identity

A native move along a Hamiltonian path `Q` from `a` to `b` sends `K` to `-K + 2A_Q`. Equivalently it takes the
converse and then reverses one Hamiltonian path (`tau_Q T = rho_(Q^rev)(T^op)`). So any functional invariant under the
converse moves under a native move exactly as it moves under reversing one Hamiltonian path. That is the whole
mechanism behind the potentials:

- *Row sums* (`d -> -d + e_a - e_b`): `L`, `delta0`, `R`, `c3` (exact change `s(a) - s(b) - 1`).
- *Cut charges* (`f(X) -> f(V\X) + [a in X] - [b in X]`): `lambda`.
- *Hamming distance* to a converse-closed bank (at most `n - 1` per move): feedback number, exception budgets,
  forest frames.
- *Unitarily invariant norms* (at most `2 ||A_Q||` per move): spectral and trace-norm budgets.

The identity also gives the non-potential bounds:

- *Summed over the upper triangle mod 2*, `F' = C(n,2) + n - 1 - F - 2 back(Q)`: THM-4810's parity theorem, and its
  missing half, that at `n ≡ 0, 3 (mod 4)` every move flips the parity, so the labelled move graph is bipartite.
- *Applied twice*, `K'' = K - 2A_P + 2A_Q`: the XOR support law of THM-4810 Lemma 1 / C2 and THM-4732 Lemma 1.
- *Restricted to endpoints*:
  - non-strong tournaments have all their Hamiltonian paths from the top component to the bottom one (THM-4729 L3's
    non-adjacency, Session74's source/sink trap);
  - a path enters and leaves every module, so the move breaks it (THM-4804).

The survey agent read every THM 4720-4814, the HYP 93xx-94xx and the 2026-10-18..27 result packets. Every other
distance bound it found is for a different move set (arc flips, circulant flips, triangle reversals) or is a static
count feeding one of the functionals.

## 2. Exactness: cut or colour

For score and cut functionals the value after the move is a function of `(T, a, b)` alone. So the bound is attained
from `T` iff there is a chain of states, each with a Hamiltonian path joining a *descending pair*. A state without one
is a local minimum, and every failure of exactness lies above one. Two different things can stop a descent:

- **cut**: no Hamiltonian path joins a descending pair. For a non-strong tournament the endpoint relation is exactly
  `C_1 x C_k` (top component times bottom component). For a strong one with a cut vertex it is restricted by forced
  arcs.
- **colour**: descending paths exist, but every one is blue (black moves only). Then the positional reversal of each is
  an anti-automorphism. That forced symmetry is exactly what the reflection saw (`R5`, the family E of THM-4811, the
  pure-blue alphabet of THM-4730).

The reflection's claim that "a forced symmetry explains why" is therefore true for colour failures and false for cut
failures. The census shows that cut failures are the common ones. The smallest rigid example is the strong
6-tournament with scores `1,1,2,3,4,4` (gentourng `111011111111101`, `|Aut| = 1`, not self-converse). Its arcs
`5 -> 3 -> 4` are forced, so no Hamiltonian path joins a minimum-score vertex to a maximum-score vertex. So `L = 1`
while the distance to the transitive class is 3.

Two classifications are exact at every order:

- (`delta0`, singleton target) the local minima are exactly the non-strong tournaments without a source or sink;
- (`L`, transitive target) a non-strong non-transitive tournament descends iff it has neither a source nor a sink.

The strong `L`-local minima number 0, 1, 3, 11, 49, 242, 1,302 at orders 4-10. All have a cut vertex, and HYP-9423 (1)
conjectures that 2-connected tournaments always descend. Its natural special case, a Hamiltonian path from a
minimum-score to a maximum-score vertex, fails for 1, 26, 590 two-connected classes at orders 7-9, though they still
descend through other pairs.

## 3. Max-to-min Hamiltonian paths and the regular class

**Theorem (THM-4815 B7).** Every tournament has a Hamiltonian path from a vertex of maximum score to a vertex of
minimum score. Let `x` be a max vertex, `y` a min vertex, and `T' = T - {x, y}`. If `y -> x`: when `T'` is strong,
injecting the out-neighbours of `x` into the non-in-neighbours of `y` through a Hamiltonian cycle of `T'` would need
`Delta <= delta - 1`, which is impossible; when `T'` is not strong, extremality puts an out-neighbour of `x` in its top
component and an in-neighbour of `y` in its bottom one. If instead every max vertex beats every min vertex, the same two
arguments go through, because each exceptional configuration produces a min vertex beating a max vertex.

The arc `x -> y` can fail to lie on a max-to-min Hamiltonian path (526 arcs in 514 classes at order 9, always with tied
extremes); the arc `y -> x` never fails (67,962 arcs at order 9). This matches Lu's min-max pair theorem with exceptions.

**Corollary (THM-4815 B8).** At odd orders the native-move distance to the regular class is exactly
`R(T) = (1/2) sum_v |s_v - (n-1)/2|`: repeatedly move along a max-to-min Hamiltonian path. So the regular target, unlike
the transitive one, never fails. Its descending pairs run in the natural direction of Hamiltonian paths (high score to
low), while the transitive target's run against it.

## 4. P_7 and parity

`P_7` (Paley labelling, `F = 14` forward arcs) is 7 moves from `TT_id` (`F = 21`, THM-4810 (D)). Every move flips
the parity at `n = 7`, so this labelled distance is odd. The distance to the transitive *class* is `L = 6`: labelled_bfs.c
(BFS over all `2^21` labelled tournaments) finds the word of paths
`0123456, 1203645, 0352641, 2316405, 4053162, 2163405`. It ends at a transitive tournament with `F = 12`, and every
step is re-verified. The "7 > 6" in the reflection part II table, the dipole note's section 1 table, Session74's
supplement and HYP-9420 (4) is the labelled number read as a class distance (MISTAKE-682). `P_11` is likewise 15 moves
from the class (Session74) and at least 16 from `TT_id`.

## 5. Census

| order | classes | class graph (any / black) | `L` local minima (non-strong + strong) | `delta0` local minima | `R` local minima (odd) |
|---|---|---|---|---|---|
| 4 | 4 | exact | 2 + 0 | 0 | – |
| 5 | 12 | exact | 5 + 1 | 0 | 0 |
| 6 | 56 | exact | 19 + 3 | 1 | – |
| 7 | 456 | exact | 99 + 11 | 3 | 0 |
| 8 | 6,880 | exact | 855 + 49 | 16 | – |
| 9 | 191,536 | in progress at filing | 13,303 + 242 | 99 | 0 |
| 10 | 9,733,056 | – | 376,191 + 1,302 | 915 | – |

Failures (distance above the bound) at orders 4-8: `L` 2, 6, 33, 182, 1,431 (symmetric 2, 4, 13, 52, 173); `delta0`
0, 0, 1, 3, 16 (all symmetric so far, each containing `C_3`- or `S_4`-blocks); `lambda` none; `R` none at odd orders.
The black class graph at orders 4-8 is one component plus the isolated pure-blue classes. The black-only failures
beyond the unrestricted ones are colour obstructions (self-converse: two at order 6, one at order 7 for the regular
target). Random endpoint tests at orders 11-16 found no counterexample to B7 or to HYP-9423 (1)
(endpoint_random_n11_16.out).

## 6. The prime cycle-index tower (THM-4816)

klein's THM-4814 (A) took the order-10 types from a gentourng census of all 9,733,056 classes and estimated about
12,200 non-rigid primes at order 11. HYP-9422 (4) called the order-11 types a multi-hour census of 903,753,248 classes.
The tower avoids every census:

- the rigid rows come from the ledger (`P(n)` and `SCP(n)` follow from the lower types, `A(n)` and `SC(n)`);
- a rigid self-converse prime has anti-automorphism type `1^1 2^((n-1)/2)` or `2^(n/2)`;
- the non-rigid primes have an automorphism of type `1^f p^c` with `c >= 2` or `f = 0` (cycle lemma), so they are
  enumerated as invariant tournaments, one bit per pair orbit, with `T[F]` up to isomorphism.

The orders 3-13 take 2, 0, 4, 32, 136, 512, 8,192, 66,048, 657,440, 15,736,832 and 373,358,656 invariant tournaments,
which is 1 s at order 11, 35 s at 12 and 26 min at 13. Checks:

- at orders 3-10 the tower's tables equal klein's census files line for line;
- at every order 3-13 the labelled prime count from `2^C(n,2)` equals `sum n!/|Aut|` over the table.

Results: order 11, 12 rows, `|Aut|` in {1, 3, 5, 11, 55}, 12,155 non-rigid; order 12, 14 rows, `|Aut|` in
{1, 3, 5, 9} with the first `Z_3^2` groups and the first `Z_9` generated by a non-regular (`3^1 9^1`) element (the circulant `Z_9` already occurs at order 9), 285,705 non-rigid; order 13, 22 rows, `|Aut|` in
{1, 3, 5, 9, 13, 39}, 11,256,145 non-rigid. Tables in THM-4816 and
`04-computation/experiments/prime_cycle_index_tower_20261027/types_tables_n11_n13.md`.

## 7. The ledgers (THM-4817)

| | n = 12 | n = 13 | n = 14 |
|---|---|---|---|
| prime `P(n)` | 142,980,435,791 | 46,453,867,277,461 | 27,687,212,017,612,715 |
| self-converse primes | 1,170,949 | 92,996,193 | 767,626,075 |
| intransitivity classes `c(n)` | 76,601,960,166 | 24,193,986,713,786 | 14,176,434,888,156,918 |
| labelled primes | 68,487,766,269,146,864,640 | 289,269,151,048,145,181,957,120 | 2,413,723,787,266,849,272,549,734,400 |

Full rows (non-strong, strong non-prime, self-converse parts, `c_s`, `s_lab`, `l`) are in THM-4817 and
`tower_n14.out`. A280772 (OEIS, to order 8) extends by 136,655, 7,906,161, 796,876,429 and the three values above. The
self-converse-prime, labelled-prime and intransitivity sequences are not in the OEIS. The order-14 Burnside defect
`3,574,801,896/7` is pre-registered (HYP-9423 (6)). Its denominator says order-14 primes with an automorphism of
order 7 exist.

## 8. Reproduction

```
cd 04-computation/experiments/prime_cycle_index_tower_20261027
cc -O3 -I/opt/homebrew/opt/nauty/include nonrigid_primes.c /opt/homebrew/opt/nauty/lib/libnauty.a -o nonrigid_primes
for n in 3 4 5 6 7 8 9 10 11 12 13; do ./nonrigid_primes $n out/nr$n.txt - ; done      # 13: 26 min
python3 burnside_counts.py 15 ; python3 tower.py 14 out > tower_n14.out ; python3 format_types.py 11 12 13
cd ../dipole_exactness_20261027
g++ -O3 -std=c++17 -I/opt/homebrew/opt/nauty/include -o dipole_exactness dipole_exactness.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
for n in 4 5 6 7 8; do gentourng -q $n | ./dipole_exactness $n details > exact_n${n}_details.out; done
cc -O3 endpoint_checks.c -o endpoint_checks ; for n in 4 5 6 7 8 9 10; do gentourng -q $n | ./endpoint_checks $n; done
cc -O2 labelled_bfs.c -o labelled_bfs ; ./labelled_bfs 7 > labelled_bfs_p7.out
cc -O3 endpoint_random.c -o endpoint_random ; ./endpoint_random 13 100000 2 132
```

## 9. Next

1. HYP-9423 (1), 2-connected `L`-descent: Hamiltonian paths against the score direction.
2. The order-9 class-graph census (black distances to the regular class), running at filing.
3. An order-14 non-rigid census needs about `1.7 * 10^10` invariant tournaments and a registry of about `7.7 * 10^8`
   classes. It wants canonical augmentation instead of a registry, and it would test HYP-9423 (6).
4. The even-order regular target (HYP-9423 (3)).
5. An audit of THM-4815 B7/B8 and of the tower's code (two independent checks are built in: klein's census and the
   labelled identity).
