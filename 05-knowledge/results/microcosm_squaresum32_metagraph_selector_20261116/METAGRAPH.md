# What a proof of Chan–Pak could use if it could only use the metagraph: the count as a potential on flips, the seam ladder, and the image selector

Owner: *"Consider what proving chan-pak would work with if we could only use the metagraph, try to find a convoluted
connection. Aim for a selector that fills every count except 7 and 21."*

## 1. The metagraph, as the repo has it

The tournament metagraph `G_n` (THM-1965, THM-282, THM-4730) has the isomorphism classes of order `n` as vertices and
flips as edges; its nodes are coloured by the symmetry of their Hamiltonian paths (pure-blue nodes = all paths
palindromic = `K_1, C_3, T_5`; blue lives on the self-converse spine, the ribs and the sea are black); it is not
vertex-transitive, and the transitive tournament is its distinguished corner. Two of this programme's results are
metagraph statements in that sense:

- **The count is a potential whose edge differences are pair-cover numbers.** Along a one-arc flip `a → b ↦ b → a`,
  `H` changes by `M(b, a) − M(a, b)` (THM-5267 (A)); along the one-vertex extension edges between `G_n` and `G_{n+1}`,
  `H(T + v_X) = h + ⟨B − e, x⟩ − xᵀMx` (THM-5266 (A)). Every edge difference is even, the potential is odd everywhere,
  and `H = 1` exactly at the corner.
- **Deletion monotonicity.** Along the extension edges `H` never decreases (THM-5267 (B)); the heaviest-parent tree is
  the gradient flow of the potential, and its trunk is the maximal-count class, not the corner (THM-5267 (E)).

## 2. The seam ladder: the convoluted connection

Call `c` a *rung* if some flip edge of the metagraph joins a class of count `c − 1` to a class of count `c + 1` — a
one-arc seam at centre `c` (THM-5267). Since every tournament of order `n + 1` is a one-contact extension of one of
order `n`, and a contact flip is an arc flip one order up, the rungs are exactly the one-contact seam centres of
THM-5267 §1.

**Proposition (PROVED, trivial).** If every even `c ≥ 24` is a rung, then every odd number `≥ 23` is a count. *Proof.*
23 is a count (order 6); if `c − 1` is a count and `c` is a rung, `c + 1` is a count. ∎

So Chan–Pak would follow from a statement purely about edges of the flip metagraph: *the count potential has a step of
exactly 2 across every even centre from 24 on*. The census (THM-5267 (C)): from bases of order ≤ 9 every even
`c ∈ [264, 10000]` is a rung, and below 264 the only non-rungs are `6, 8, 20, 22` (forced: a neighbour is a hole) and
`36, 40, 54, 56, 60, 64, 84, 94, 108, 148, 188, 196, 262`, which the ladder jumps by other realisers (each of their
neighbours is a count anyway). The two solo separators of [MICROCOSM.md](MICROCOSM.md) are the two forced jumps of the
ladder: `5 → 9` and `19 → 23`, the only places where a monotone walk of the potential must step by 4 instead of 2.
Conjecture (HYP-9864): every even `c ≥ 264` is a rung; it implies Chan–Pak and is FINITE-EXACT to 10,000. It is a
density statement about the antisymmetric part of the pair matrix, not a structural one, exactly as THM-5265 found for
the near-join closure; the metagraph gives the ladder its shape, not its proof.

## 3. The image selector at one level ([receipts/cover.out](receipts/cover.out), `cover.py`)

The new counts of order 9 (1,200 odd numbers first born there) are, by completeness of the extension grammar, values of
the quadratic forms of order-8 classes. Which classes select them?

- The 5,457 first creators of order 8 cover all 1,200 (M-6 confirmed; the weak form of HYP-9863).
- A greedy cover needs only **31 classes** of order 8 (M-5's "fewer than 60" confirmed); the best single class, of count
  393, covers 226 of the 1,200 (18.8 %, M-5's "25 %" refuted), the next of counts 597, 285, 261, 627 add 180, 164, 119,
  102.
- The six classes of maximal count 661 are poor selectors: their images hold 106–172 new counts and rank around 4,843 of
  6,880 (M-5's last clause refuted). The heaviest lineage is the trunk of the gradient flow but not the best bank for
  selection; mid-count classes with many distinct extension values are.

So a metagraph-only selector at one level is a bank of a few dozen classes of order `n` whose cube images cover the new
counts of order `n + 1`; across levels the bank grows with the order and the proof of coverage would again be a density
statement. [SELECTOR.md](SELECTOR.md) assembles the best selector this programme can state and measures its reach.
