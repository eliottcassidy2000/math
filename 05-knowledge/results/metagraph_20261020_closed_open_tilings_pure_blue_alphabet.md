# Merged tournament metagraph, 2026-10-20 (mac-mini-2026-10-20-metagraph): closed/open tilings, the forced arrangement, and the complete strong pure-blue alphabet

Status: pre-audit. Canon written: THM-4729, THM-4730, HYP-9345. Scripts and outputs in
`04-computation/experiments/metagraph_20261020/`.

## 1. What the owner asked and what was already canon

The owner described the merged metagraph (tilings over a fixed Hamiltonian path on the
triangular grid, complement lines coloured blue/black by mirror symmetry, converse classes
merged, points of symmetry = nodes with only blue lines, self-lines) and asked for
comprehensive patterns, count formulas, new lines, forced arrangements and terminology.
Before writing, the canon inventory (Explore agent) showed most elementary counting laws
already proved: fibre `H/|Aut|` odd (THM-643 T1), pure-black = non-SC (T2), blue parities
(T4, THM-796 (4.8)), black multiplicities even (THM-796 (4.7)), blue/black totals (T5), blue
self-lines only at even `n` (THM-023/648), the `n = 10` blue self-loop count `24` refuting
`2^{n/2-2}` (klein-S173), the pure-blue (= point of symmetry) alphabet law with centre law
and census to `n = 13` (THM-2444/2453/2454), rigid strong = C3 through 62 (THM-2477). The
repo's name for a point of symmetry is **pure-blue**; "POS" elsewhere in the repo means
unrelated things. The genuinely open items were: completeness of the strong pure-blue
alphabet `{1, C3, T5}` beyond `n = 13` (first odd case `n = 15`, all even sizes via
centres), the mechanism of the blue self-loops, black self-lines `= SC/2` at odd `n`.

## 2. THM-4730: the strong pure-blue alphabet is complete (PROVED, machine-checked)

A strong tournament all of whose Hamiltonian paths are palindromic (positional reversal
an anti-automorphism) is `K1`, `C3` or `T5`.

- Rotation lemma: for a Hamiltonian cycle `(v_0..v_{n-1})` of a pure-blue tournament the
  reflections `v_x -> v_{2i-1-x}` are anti-automorphisms, so `v_x -> v_{x+2}` is an
  automorphism; for odd `n` the shift itself is, so the tournament is circulant with
  `1 in S` with respect to every Hamiltonian cycle.
- Odd order `>= 7`: every circulant with `1 in S` has an explicit non-palindromic
  Hamiltonian path: the three-interval path `(0..a-2, a+b-2..n-1, a-1..a+b-3)` whenever two
  distinct `a, b in S\{1}` have `a + b <= n+1` (always, except for
  `S = {1} u [(n+1)/2, n-2]`), and the descending-pairs path `(0,1,n-2,n-1,n-4,n-3,...,3,4,2)`
  for that one set. Checked for all 2,047 sets at `n = 7..23`.
- Even order: `n = 0 (mod 4)` gives an involutory automorphism; `n = 2 (mod 4)` makes the
  tournament a `Z_{n/2}`-bicirculant `(S', D)` with a mixing Hamiltonian cycle
  `(0, 2s, ..., -2s, d-2s, d, ..., d-4s)` whose two-step rotation moves an even vertex to
  an odd one, so `Aut` is vertex-transitive -- impossible at even order. Checked on all
  1,198,180 bicirculants at `n = 6, 10, ..., 30`; brute force at `n = 6, 10`.
- Consequences: THM-2454's law `pure-blue(n) = c(n/2)` (even), `c((n-1)/2) + c((n-3)/2)
  + c((n-5)/2)` (odd), `c = A000930`, is PROVED for all `n <= 125` (THM-2477 supplies the
  rigid halves) and equivalent to THM-2477's lemma ML beyond; THM-643 C3 (pure-blue fibres
  `1` or `3`) proved; `T_p` is not pure-blue for any prime `p >= 7`.
- Independent finite checks: no strong pure-blue tournament of order 6, 7, 8 (all strong
  classes); circulants of odd order `5..17`: only `T5` (affine anti-automorphisms, exact
  at prime order by Burnside's prime-degree theorem; composite orders 9, 15 redone with
  all involutory anti-automorphisms); `Z_3 x Z_3` Cayley tournaments not pure-blue.

## 3. THM-4729: closed/open tilings and the forced arrangement (PROVED + FINITE-EXACT n <= 8)

- A tiling is **closed** when the tile `(n,1)` carries `1 -> n` (base path + arc =
  Hamiltonian cycle), **open** otherwise. The complement flips `(n,1)`, so every line joins
  a closed to an open tiling; `2^{m-1}` closed tilings, `2^{(m+f)/2-1}` closed symmetric.
- **Closure number** `cl(C) = n HC(C)/|Aut C|` (rooted Hamiltonian cycles modulo the free
  `Aut` action); positive iff strong (Camion).
- **Forced arrangement**: non-strong nodes are fully open, hence pairwise non-adjacent and
  loopless (all their lines end at strong nodes); fully closed nodes (`H = n HC`) likewise;
  every self-line sits at a strong node with both closed and open tilings; for `n >= 6`
  every point of symmetry is non-strong, so all its blue lines end at strong mixed nodes.
- **Strong/non-strong line count** = number of non-strong tilings
  `NS(n) = sum_k a(k) 2^{C(n-k-1,2)}` (first cut), `a(n) = 2^{C(n-1,2)} - NS(n)` strong
  tilings: `a = 1,0,1,5,50,903,30773,2032504`, `NS = 0,1,1,3,14,121,1995,64648` (`n = 1..8`;
  not in OEIS); census agrees at `n = 4..8`.
- **Convention independence**: complementing the non-path arcs, reversing the path arcs,
  and complement-then-mirror induce the same line system on the merged metagraph and
  different ones on the unmerged class graph: the merge is the convention-free quotient.
- Census `n = 4..8` (3,528 nodes, 1,048,576 lines at `n = 8`): zero violations;
  metagraph and blue graph connected; fully closed classes only `C3` (`n = 3`) and
  `Q5 = a -> C3 -> b -> a` (`n = 5`, `(H,|Aut|) = (15,3)`); `n = 4` metagraph forced by
  the laws; `n = 5` table with closure column (blue system not forced by degrees alone).
- Blue self-line classes `n = 4, 6, 8, 10` (`1, 2, 4, 24`, reproducing klein-S173): one
  self-line per class, one closed and one open tiling, all strong and near-regular,
  `|Aut| = 1` except one class each at `n = 6, 8, 10` with `|Aut| = 3, 3, 5`.
- Negative results: the staircase transpose and rotation, and half flips of the grid,
  preserve neither class nor `H` nor scores (`staircase_symmetries.out`, `half_flips.out`).

## 4. HYP-9345 (open)

(a) closure conjecture `H > n HC` for strong `n >= 6` (exact through n = 9, section 6);
(b) connectivity for all `n`; (c) blue self-line classes: one per class, near-regular,
count mechanism (`n = 12` discriminates); (d) black self-lines `= SC/2` at odd `n`
(`n = 9` untested); (e) THM-2477's ML (rigid strong composite orders `>= 63`; a rigid
strong tournament is circulant w.r.t. every Hamiltonian cycle with `S n U = M` an
index-two multiplier group); (f) the sequence `a(n)`; (g) the determination problem;
(h) rotation lines.

## 5. Files

`merged_metagraph.py` (`_n7.out`, `_n8.out`), `closed_open_laws.py` (`_n7.out`, `_n8.out`),
`blue_self_lines.py`, `blue_self_lines_sym.py`, `blue_self_lines_classes.py` (`.out`),
`atom_alphabet_proofcheck.py` (`.out`), `pos_exact.py` (`.out`), `circulant_pos_fast.py`
(`.out`), `palindromic_paths.py` (`.out`), `staircase_symmetries.py`, `half_flips.py`,
`fully_closed_n9.py` (`.out`).

## 6. Addendum: the n = 9 closure census

`fully_closed_n9.py` (gentourng -c 9, bitmask DP for H and HC on every class): among the 178,133 strong tournaments of order 9 none has H = 9 HC. The closure conjecture of HYP-9345 (a) now stands FINITE-EXACT for n <= 9 (fully closed strong classes: C3 at n = 3, Q5 at n = 5, none at n = 4, 6, 7, 8, 9).
