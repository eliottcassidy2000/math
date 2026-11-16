# The out-degree-1 law for digit-0 vertices, proved for every tournament

Session klein-2026-11-16-outdeg (2026-11-16). Owner's request: *"Prove the out-degree-1 law for digit-0 vertices of
strong tournaments."* THM-5263 had it FINITE-EXACT (orders ≤ 8). Here it is proved at every order, together with the
two-path lemma it needs, the complete characterisation of digit-0 vertices, and the uniqueness of the digit-0 vertex in a
strong tournament. Canon: THM-5264 (new); THM-5263's status upgraded.

## 0. Setting (THM-4865, THM-5263)

`S` is a tournament on `V`, `|V| = n ≥ 2`, `u ∈ V`. `h(X) = H(X)` is the number of Hamiltonian paths of `X` (`h(∅) = 1`),
`s_u(X)` the number starting at `u`, `e_u(X)` the number ending at `u`. The root response is
`B_u(S) = Σ_{A ⊆ V∖u} h(S[A]) · s_u(S[V∖A])`: the number of ordered pairs `(P_1, P_2)` with `P_1` a Hamiltonian path of
`S[A]` (empty when `A = ∅`) and `P_2` a Hamiltonian path of `S[V∖A]` starting at `u`. The digit is `σ_u = B_u − H(S) ≥ 0`.
THM-5263 (I1): `σ_u` is the number of such pairs with `P_1` nonempty and `u → last(P_1)` (a *backward junction*).
THM-5263 (I2): if `d⁺(u) = 1` with out-neighbour `t_0`, then `σ_u = e_{t_0}(S − u)`.

Classical facts used throughout: every tournament has a Hamiltonian path (Rédei); a strong tournament of order ≥ 3 has a
Hamiltonian cycle (Camion), hence a Hamiltonian path starting at any prescribed vertex and one ending at any prescribed
vertex (trivially also for a single vertex). If `C_1 → … → C_m` are the strong components of a tournament in condensation
order (every arc between distinct components goes from the lower index to the higher), concatenating Hamiltonian paths of
`C_1, …, C_m` gives a Hamiltonian path; every Hamiltonian path ends in `C_m`, and every vertex of `C_m` ends one.

## 1. The two-path lemma

**Lemma 1.** Let `C` be a strong tournament of order ≥ 3 and `x ≠ x′` two of its vertices. Then `V(C)` is the disjoint
union of the vertex sets of two directed paths, one of which starts at one of `x, x′` while the other ends at the other
one. Either path may be a single vertex.

*Proof.* Let `E_1 → … → E_q` be the strong components of `C − x′` (`q ≥ 1`).
*Case 1, `x ∈ E_1`.* `P_a` = a Hamiltonian path of `E_1` starting at `x`, followed by Hamiltonian paths of `E_2, …, E_q`;
`P_b = (x′)`. `P_a` starts at `x`, `P_b` ends at `x′`.
*Case 2, `x ∈ E_b`, `b ≥ 2`.* No vertex of `E_2 ∪ … ∪ E_q` beats a vertex of `E_1`, and `C` is strong, so `x′ → e` for some
`e ∈ E_1`. `P_a = (x′`, then a Hamiltonian path of `E_1` starting at `e`, then Hamiltonian paths of
`E_2, …, E_{b−1}, E_{b+1}, …, E_q)`; `P_b` = a Hamiltonian path of `E_b` ending at `x`. Consecutive blocks of `P_a` are
joined by forward arcs (`E_i → E_j` for `i < j`, in particular `E_{b−1} → E_{b+1}`). `P_a` starts at `x′`, `P_b` ends at
`x`. ∎

The ordered form — "a path from `x` and a path ending at `x′`, for every ordered pair" — is false: `twopath.py` finds 726
ordered pairs `(x, x′)` among the 352,462 in strong tournaments of orders 3–8 with no such partition, and no failure of
the unordered form ([receipts/twopath.out](receipts/twopath.out)). The unordered form is all the law needs, because the
two vertices it is applied to play symmetric roles.

## 2. The law

**Theorem 2 (out-degree-1 law).** If a vertex `u` of a tournament `S` has out-degree at least 2, then `σ_u(S) ≥ 1`.
Equivalently: a digit-0 vertex has out-degree at most 1 — in a strong tournament, exactly 1.

*Proof.* Let `C_1 → … → C_m` be the strong components of `S − u`, `j` the largest index with `N⁺(u) ∩ C_j ≠ ∅`, and
`x′ ∈ N⁺(u) ∩ C_j`. By (I1) it suffices to exhibit paths `P_1, P_2` partitioning `V`, `P_2` starting at `u`, `P_1`
nonempty, with `u → last(P_1)`.
*Case A: `u` has an out-neighbour `x ∈ C_i` with `i < j`.* `P_2 = (u, x, Ham(C_{j+1}), …, Ham(C_m))` (`x` beats all of
`C_{j+1}` because `i < j + 1`; the tail is empty when `j = m`) and
`P_1 = (Ham(C_1), …, Ham(C_{i−1}), Ham(C_i − x), Ham(C_{i+1}), …, Ham(C_{j−1}), Ham(C_j ending at x′))`, the block
`C_i − x` omitted when empty. Every block is a tournament, so has a Hamiltonian path; `C_j`, strong or a single vertex,
has one ending at `x′`; consecutive blocks are joined by forward arcs. `last(P_1) = x′` and `u → x′`.
*Case B: all out-neighbours of `u` lie in `C_j`.* There are at least two, so `C_j` is strong of order ≥ 3. Lemma 1 for
`C_j` and two out-neighbours `x, x′` of `u` gives disjoint paths `Q_a, Q_b` covering `C_j`, `Q_a` starting at one of them
and `Q_b` ending at the other, call it `y`. `P_2 = (u, Q_a, Ham(C_{j+1}), …, Ham(C_m))`,
`P_1 = (Ham(C_1), …, Ham(C_{j−1}), Q_b)`; `last(P_1) = y` and `u → y`. ∎

**Theorem 3 (digit-0 vertices).** `σ_u(S) = 0` if and only if `u` is a sink of `S`, or `d⁺(u) = 1` and the out-neighbour
`t_0` of `u` lies outside the terminal strong component of `S − u` (equivalently `e_{t_0}(S − u) = 0`).

*Proof.* A sink admits no backward junction. For `d⁺(u) = 1`, (I2) gives `σ_u = e_{t_0}(S − u)`, which vanishes exactly
when `t_0` is not in the terminal component. For `d⁺(u) ≥ 2`, Theorem 2. ∎

**Corollary 4 (strong tournaments).** Let `S` be strong, `z` a vertex with `d⁺(z) = 1` and out-neighbour `t_0`, and
`C_1 → … → C_m` the components of `S − z`. Then `t_0 ∈ C_1` (every vertex other than `t_0` beats `z`, and the only arc out
of `z` goes to `t_0`, so strongness forces every vertex of `S − z` to be reachable from `t_0` inside `S − z`), and `z` is a
digit-0 vertex iff `t_0 ∉ C_m`, i.e. iff `m ≥ 2`. For such `z`, `H(S) = h(S − z) + B_{t_0}(S − z)` (THM-5263 (I2) with
`e_{t_0} = 0`): **every strong tournament with a digit-0 vertex is a dead-end append, at every order.** Conversely, for a
tournament `B′` and a vertex `t_0` of `B′` outside its terminal component, appending `z` with `z → t_0` and `B′ − t_0 → z`
gives a tournament with digit-0 vertex `z`, strong iff `t_0` lies in the initial component of `B′`.

**Theorem 5 (uniqueness).** A strong tournament of order ≥ 4 has at most one digit-0 vertex. (`C_3` has three.)

*Proof.* Let `z_1 ≠ z_2` be digit-0 vertices; both have out-degree 1 (Theorem 2). Say `z_1 → z_2`, so `N⁺(z_1) = {z_2}`;
let `t_2` be the out-neighbour of `z_2` (`t_2 ≠ z_1`). Digit 0 at `z_1` means `z_2 ∉ C_m(S − z_1)`. Any `w ∈ C_m(S − z_1)`
other than `t_2` beats `z_2` (the only vertex `z_2` beats is `t_2`), which would be an arc leaving the terminal component;
so `C_m(S − z_1) = {t_2}`, `t_2` is a sink of `S − z_1`, and strongness of `S` forces `t_2 → z_1`. Now `z_1, z_2, t_2` form a
3-cycle in which each vertex has out-degree exactly 1, so every other vertex beats all three: `S = (S − {z_1, z_2, t_2}) → C_3`,
not strong when the order is ≥ 4. ∎

Consequences. (a) Any strong tournament with minimum out-degree ≥ 2 — regular tournaments, every vertex-transitive
tournament of order ≥ 5, the Paley tournaments — has no digit-0 vertex: it is *sterile* ([DESCENT.md](DESCENT.md)).
(b) The one-exit law used in THM-5263 (B) holds at every order, so the dead-end-append representation of strong digit-0
counts and the dichotomy of [CONVERSE.md](CONVERSE.md) are unconditional. (c) On strong tournaments of order ≥ 4 the
map "delete the digit-0 vertex" has no choice in it: the descent of DESCENT.md is a genuine one-outward-edge system.

## 3. Replay

- `converse_descent.py` (1): for every tournament of orders 3–8 and every vertex (58,650 vertices), "`σ_u = 0`" agrees
  with "sink, or sole exit whose port ends no Hamiltonian path of `S − u`": 0 mismatches
  ([receipts/converse_descent.out](receipts/converse_descent.out)).
- `twopath.py`: Lemma 1 has no failure on the 352,462 ordered pairs of strong tournaments of orders 3–8; the ordered form
  fails 726 times.
- `closure_check.py` (P-12) and `sterile9.c`: no strong tournament of orders 4–9 has two digit-0 vertices (the largest
  number of digit-0 vertices in a strong class is 1 at every order 4–9; `C_3` has 3)
  ([receipts/closure_check.out](receipts/closure_check.out), [receipts/sterile9.out](receipts/sterile9.out)).
