# The free gas and the B₂ atom: the OCF as a lattice gas, and triangular numbers ≅ the four 4-tournaments

**Source:** mac-mini-2026-06-15-S2 (disk-constrained — hand-verified at n=4, scripts
deferred). Dispatch: the mean-field "free board" vs interacting tournament; and the
triangular-number factorization T(x)=f·g ≅ the 4 iso classes of T₄. Canon: T821,
THM-510, HYP-2529..2531, OPEN-Q-100. Builds on THM-002/T006 (OCF = hard-core lattice
gas), THM-029 (H=7 forbidden), THM-509 (det/permanent wall), THM-468 (diamond det=9).

## Part 1 — The OCF is a lattice gas; the free board is the ideal gas

`H(T) = I(Ω, 2)` (THM-002): the independence polynomial of the odd-cycle conflict
graph Ω (vertices = directed odd cycles, edges = shared vertex) at fugacity z = 2 —
a **hard-core lattice gas** (T006). Read statistical-mechanically:

- **The free board = the ideal gas.** If Ω had NO edges (all α₁ odd cycles mutually
  disjoint-able), `H_free = Σ_k C(α₁,k) 2^k = (1+2)^{α₁} = 3^{α₁}` — the **mean-field
  reference**, every cycle an independent particle.
- **The tournament is the interacting gas.** Real Ω has edges (cycles forced to share
  vertices), so `H = Σ_k 2^k α_k ≤ 3^{α₁}`; the deficit `3^{α₁} − H` is the
  **interaction energy** (the cost of forced overlaps).
- **Forbidden values = interaction obstructions** — values no *real* Ω realizes even
  though the *free* board's α-vector formally gives them.

**"Why is 7 forbidden?"** `7 = 1 + 2·3` forces `α = (3,0,…)`: three odd cycles, no two
disjoint ⟹ `Ω = K₃`. The free board writes `α=(3,0)` happily; but (THM-029) no
tournament's conflict graph is exactly K₃ — three pairwise-intersecting odd cycles
always force a disjoint pair (`α₂>0`) or a fourth cycle (`α₁>3`). **The
three-pairwise-intersecting-cycles obstruction is an interaction the free board lacks**
— exactly the dispatch. 7 and 21 are the lattice gas's **excluded volumes**.

- **Mayer/virial dictionary.** The overlap/Witt defects (THM-502/505): `p33` =
  intersecting triangle pairs = the **edges of Ω** = the first **Mayer bond**; `TQ`,
  `Q44`, … = higher cluster integrals. The non-spectral carriers ARE the cluster
  expansion of `log H` about the ideal gas.
- **The wall (THM-509).** Ideal gas = product form = spectral / determinant / P side;
  interacting independence polynomial = packing count = permanent / #P side. Forbidden
  values live where the #P-side image misses a free-gas value. The baby-Hodge holes
  ARE the lattice gas's interaction-forbidden free-α-vectors.

## Part 2 — Triangular numbers and the B₂ atom: T(x)=f·g ≅ the four 4-tournaments

`T(x)=x(x+1)/2`. Two atomic ops: **a(x)=x+1** (successor/additive), **b(x)=x/2**
(halving/dyadic); `T(x)=x·b(a(x))`. Split by parity to keep factors integral:

| | f(x) | g(x) | f as op-set | g as op-set |
|---|---|---|---|---|
| x even | x+1 = a | x/2 = b | **{a}** | **{b}** |
| x odd | x = id | (x+1)/2 = b∘a | **∅** | **{a,b}** |

The four faces are the **power set B₂ = {∅,{a},{b},{a,b}}**. The dispatch's pattern
(*empty / both-combined / one-of-each*) = the **four iso classes of T₄**, graded by
**c₃ = |subset|** (all hand-verified):

| subset | c₃ | 4-tournament | score | transpose | det(S) |
|---|---|---|---|---|---|
| ∅ | 0 | transitive | (0,1,2,3) | self | 1 |
| {a} | 1 | out-vortex diamond | (1,1,1,3) | ↔ {b} | 9 |
| {b} | 1 | in-vortex diamond | (0,2,2,2) | ↔ {a} | 9 |
| {a,b} | 2 | strong | (1,1,2,2) | self | 1 |

1. **c₃ = |subset|.** Counts `0,1,1,2` = subset sizes of a 2-set; distribution
   `{0:1,1:2,2:1} = (1,2,1)` = **Pascal row 2** (`C(2,k)`, the d=1 row of the
   Pascal-slope family T819), summing to `4 = 2² = |B₂| = A000568(4)`.
2. **Transpose = the a↔b swap.** Arc reversal fixes transitive (∅) and strong ({a,b}),
   **swaps the two diamonds**. The B₂ involution fixing ∅,{a,b} and swapping {a},{b}
   IS "swap a↔b." So **transposing a tournament = swapping the additive (+1) and dyadic
   (/2) operations** in the triangular factorization; self-transpose classes = the
   symmetric subsets, diamonds = the asymmetric pair.
3. **Parity of x ↔ transpose-type.** Odd x uses {∅,{a,b}} = the **self-transpose** pair
   (transitive, strong); even x uses {{a},{b}} = the **transpose-swapped** diamonds.
   The Z₂ governing T(x)'s integer split = the transpose symmetry of the matching pair.

Bonus (THM-468, Belkouche): `det(S)=9` (`|Pf|=3`) for the diamonds, `det(S)=1`
(`|Pf|=1`) for transitive & strong. The Pfaffian flags the **single-operation
(off-diagonal) pair** — the **vortices**, the obstructions to local-orderness (the d=1
floor code). The diamond is where the **interaction first appears** (a vertex glued to
a 3-cycle): n=4 is the smallest stage where free board and interacting tournament part.

## The synthesis

a = successor (additive/Peano), b = halving (dyadic) — the **additive × multiplicative
atom**, the same pair behind Goldbach-vs-Euler duality and the doubled-prime ops. The
four 4-tournaments are the 2×2 truth table "apply +1?" × "apply /2?", transpose = the
swap of axes. T(x) (arc-count of K_{x+1}, the staircase size) is assembled by
distributing these two ops, and the *way* they distribute (split vs concentrated) is
read by the parity of x = the transpose-symmetry of the corresponding 4-tournament.

Free gas (Part 1) and the B₂ atom (Part 2) meet: the **diamonds are the first
interaction** (c₃=1, det=9 vortices), the **transitive is the free/acyclic ground
state** (c₃=0), the **strong is fully interacting** (c₃=2). The n=4 Boolean atom is the
lattice gas's first nontrivial cell — where the mean-field free board and the
interacting Ω first diverge.

## Open

- **OPEN-Q-100:** make the Mayer cluster expansion of `H=I(Ω,2)` about `3^{α₁}`
  explicit — Witt/overlap defects (p33, TQ, …) as virial coefficients, forbidden {7,21}
  as excluded volumes; does it give a cleaner proof that {7,21} are the only permanent
  gaps (HYP-2271) than the numerical-semigroup argument?
- Does B₂ ≅ T₄ seed a general-n structure, or is the *count* match special to n=4
  (`A000568(4)=4=2²`, Pascal row 2)? The *gradings* match (c₃=|subset|, transpose=swap,
  parity↔type) is the real, robust content.
