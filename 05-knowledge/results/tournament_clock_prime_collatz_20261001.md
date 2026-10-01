# Tournament Clock Prime Collatz: residue clocks with one-way, doubled, missing and looped pairs

**Thread session `tcpc-thread-20261001` (temporary identity, no `.machine-id`), 2026-10-01.**
Branch `claude/project-thread-au0gxe`.

- Script: `04-computation/experiments/tournament_clock_prime_collatz_20261001.py` (+ `.out`,
  56 checks, ALL CHECKS PASSED, about a minute, standard library only).
- Census script: `04-computation/experiments/tournament_clock_teichmuller_census_20261001.py` (+ `.out`,
  ALL CHECKS PASSED, about 6 minutes; §6.2 and §6.6).
- Prior work used: DB1 of [`arithmetic_braids_20260917_divisors.md`](arithmetic_braids_20260917_divisors.md)
  (`F = S + U` iff `N = p, p^3, p^2qr`); the `F_2^3`/Fano sidecar of
  [`arithmetic_braids2_20260917_squarefree_symmetry.md`](arithmetic_braids2_20260917_squarefree_symmetry.md);
  the mod-9 inverse chain of [`collatz_mod6_20260917_three_adic_g_map.md`](collatz_mod6_20260917_three_adic_g_map.md);
  THM-4520 (`collatz-level-operator-is-a-gauss-twisted-circulant…`); THM-4523
  (`collatz-functional-graph-is-rigid…`, Theorem R_a); guardrails 4, 27, 41, 42 and MISTAKE-232.
- No canon ID is reserved. This is a research note.

**Status.**
- PROVED (elementary; several statements are classical and are marked so): Theorems 1–3, 5, 7(a, b),
  9, 10, 12–14; Propositions 4, 6, 8, 11, 15, 16.
- FINITE-EXACT: Theorem 7(c) (the converse prime-atom rule), the "only if" half of the clique-module
  rule, all tables and censuses.
- HEURISTIC: the Poisson(1/6) model for Fermat case I modulo `p^2`; the `A/2` density of Paley
  doubling primes (a GRH proof is expected in the near-primitive-root framework but is not cited here).
- Literature: CITED statements in §6 were checked against the sources listed in §6.6 on 2026-10-01
  (several only through abstracts or secondary sources, as noted). Theorem 14 is an explicit evaluation
  inside Tao's published framework and is not claimed as new.
- **No open problem is solved here.** TCPC computes the local (congruence) layer exactly. Every open
  problem it touches sits where that layer stops (§7).
- Independent audit (2026-10-01): SOUND WITH CORRECTIONS. Every correction is applied (§11).

The owner's prompt (2026-10-01) asked for four things:
1. a principle for triplets: one member primarily differentiated, the other two subtly differentiated
   "in an even/odd halving/doubling manner", like the legs of a Pythagorean triple;
2. the five complementary pairs of `p^2qr`, the 3-4-3 split of its ten proper divisors, and
   "`p^2qr` is `p` and `p^3` combined … we can employ its fractal properties";
3. squares mod 9 (`0,1,4,0,7,7,0,4,1,0`) versus cubes mod 9 (`0,1,8,0`), and why mod 9 matters for
   Collatz;
4. "Tournament Clock Prime Collatz": tournaments with missing, doubled or looped edges on partitions of
   a modulus, laws for how they combine, many families of numbers, grounded in open problems.

| Question | Answer | Label |
|---|---|---|
| What is a clock? | `Cay(Z/n, D)`: each pair is one-way, doubled or missing; loops iff `0 ∈ D` | definition (§0) |
| How do clocks combine? | CRT is the AND of pair types: doubled/loop is the identity, missing absorbs, opposite arcs annihilate | PROVED (Thm 1) |
| Which power clocks are tournaments? | Only the Paley ones: `n = p ≡ 3 (mod 4)`, `gcd(k, p−1) = 2` | PROVED (Thm 2) |
| Do orientations add or multiply? | Set products AND them; character products multiply them: two Paley tournaments make an undirected graph | PROVED (Thm 3) |
| Where is the fractal? | At prime powers: squares mod `3^m` = `C3[K̄3[C3[…]]]`, cubes mod `3^m` = `C9[K̄3[C9[…]]]` | PROVED (Thm 5) |
| The periods 9 and 3 | On units mod `p^m` with `m > v_p(k)`, `x^k` forgets `v_p(k)` input digits, and its image needs `v_p(k)` extra digits | PROVED (Prop 6) |
| "Prime" | Prime moduli give prime (indecomposable) digraphs; CRT dissolves towers unless a factor has a clique module | PROVED (7a, 7b) + FINITE-EXACT (7c) |
| The triplet principle | Exactly the transitive triangle TT3: the closure arc is fixed by its unique anti-automorphism, the two path arcs are swapped | PROVED (Prop 8) |
| Pythagorean legs | Exchangeable mod `p` iff `p ≢ 5 (mod 8)`; the even leg is a quadratic residue modulo every prime of the hypotenuse, the odd leg only modulo those `≡ 1 (mod 8)` | PROVED (Thm 10) |
| 3-4-3 at `p^2qr` | Two 3-4-3 splits differing by `p ↔ qr`; the five pairs form a 3-vertex multigraph with a doubled edge and a loop; `r = 3` is forced | PROVED (Prop 11) |
| Squares vs cubes mod 9 | `(Z/9)^× = μ2 × U1` = cubes × squares; a unit generates iff it is neither a square nor a cube, and 2 is neither | PROVED, classical (Thm 12) |
| Collatz on the mod-9 clocks | `T(n) ≡ (−1)^v 4^(n+v) (mod 9)`; orbit law `(8,16,11,4,2,22)/63` | PROVED (Thms 13, 14) |
| `an + 1` maps | `(an+1)/2^v ≡ ω(2)^(−v)(1+a)^(n+q_a(2)v) (mod a^2)`: halving turns the Artin factor by `ω(2)^(−1)` and the Wieferich factor by the Fermat quotient | PROVED (Thm 13(b)) |
| Can clocks decide Collatz? | No: rational cycles carry every clock reading modulo `2^a 3^b` | PROVED, classical (Prop 16) |
| Open problems | Artin and Wieferich (the condition of THM-4523's Theorem R_a), Fermat case I mod `p^2`, `G(3)`, three cubes, doubling primes | §6 |

## 0. Definitions, naming and loss ledger

**Residue clock.** For `n ≥ 2` and `D ⊆ Z/n`, the clock `Cl(n, D)` is the Cayley digraph `Cay(Z/n, D)`:
`x → y` iff `y − x ∈ D`. A pair `{x, y}` with `x ≠ y` has the type `τ(d) = ([d ∈ D], [−d ∈ D])` of
`d = y − x`:

- `10` or `01`: one-way;
- `11`: doubled;
- `00`: missing.

There is a loop at every vertex iff `0 ∈ D`. Since `τ(−d)` is `τ(d)` with its bits swapped, all four
types are honest data. Nothing is forced into a tournament (guardrail 41).

Families used below:
- **power clocks** `D_k(n) = {x^k mod n}` (they contain 0, so they are looped);
- **character clocks** (§1);
- **doubling clocks** `Cay(Z/p, ⟨2⟩)` (§6.4);
- the **Teichmüller clock** `D_p(p^2)` (§6.2).

**Per-vertex counts.**
- one-way out-arcs: `|D| − |D ∩ −D|`;
- doubled neighbours: `|D ∩ −D| − [0 ∈ D]`;
- missing neighbours: `n − 1 − 2|D| + |D ∩ −D| + [0 ∈ D]` (for power clocks, which contain 0, this is
  `n − 2|D| + |D ∩ −D|`).

Examples (Check A):
- squares mod 9, `D = {0,1,4,7}`: 3 one-way, 0 doubled, 2 missing (`d = ±3`);
- cubes mod 9, `D = {0,1,8}`: the 9-cycle with loops (2 doubled, 6 missing);
- squares mod 7: the Paley tournament `P7`;
- squares mod 5: the 5-cycle.

**Naming.** The repository already uses "clock" for other objects: LRC owner clocks (guardrails 14–27),
the necklace clock `2^K − q^X`, THM-3356's "prime-clock resultants", and an earlier "tournament clock",
the lonely-runner phase comparator of oracle S24 (HYP-1951 `runner-tournament-clock-circular-menu`,
THM-373 `runner-phase-clock-wall-decomposition`, and the S502 wall overlay
`lrc_tournament_clock_overlay_s502.py`). This note's object is a *residue clock*, and no
identification with any of those is claimed (MISTAKE-233: clocks are not modular cusps).

**Loss ledger (guardrail 4).**
- Source: `Z` (or `Z_p`) with its order.
- Target: the pair-type function on `Z/n`.
- Map: reduction mod `n`, then `d ↦ τ(d)`.
- Preserved: which differences are `k`-th powers mod `n`, hence every congruence obstruction.
- Lost: size, sign and order (the archimedean place); carries beyond the modulus; multiplicities
  (how many `x` have `x^k = d`).
- Restoration sidecar: the tower over `p^m` (§2) restores `p`-adic digits. Nothing in TCPC restores the
  archimedean data.
- Hostile: the rational Collatz cycles of Proposition 16 and the negative integer cycles. They carry
  every reading modulo `2^a 3^b` of a hypothetical positive cycle with the same halving word.

## 1. Combination laws

**Theorem 1 (AND law; PROVED).** Let `n = ab` with `gcd(a, b) = 1`. Under CRT,
`D_k(n) = D_k(a) × D_k(b)`, and the type of every difference is the componentwise AND:

```text
τ_n(d) = τ_a(d mod a) ∧ τ_b(d mod b).
```

In the monoid `({11, 10, 01, 00}, ∧)`:
- doubled/loop `11` is the identity;
- missing `00` is absorbing;
- `10 ∧ 01 = 00`: opposite arcs annihilate.

Consequently `|D_k(n)|` and `|D_k(n) ∩ −D_k(n)|` are multiplicative, and so are the per-vertex counts'
ingredients.

*Proof.* `x ↦ x^k` commutes with CRT, and `d ∈ D_k(n)` iff both components are `k`-th powers. ∎

Checked for the moduli 63, 72, 35, 80, 108, 143 with `k = 2, 3, 4`; multiplicativity for `n < 400`,
`k = 2..6`. The same AND law holds for any product set `D_a × D_b`: the clock of the product is the
categorical product of the factor clocks.

*Example.* In the square clock mod `63 = 9·7`, `d = 1` is one-way (`10 ∧ 10`). But `d = 55`, which is
`1 mod 9` and `−1 mod 7`, is missing (`10 ∧ 01`): the two orientations disagree and cancel.

**Theorem 2 (tournament rigidity; PROVED).** For `n ≥ 2` and `k ≥ 2`, the power clock is a tournament
(every nonzero difference one-way) iff `n = p` is prime, `p ≡ 3 (mod 4)` and `gcd(k, p−1) = 2`. These are
exactly the Paley tournaments.

*Proof.*
- Two coprime factors `a, b ≥ 2`. If `−1 ∈ D_k(a)`, the difference `(1 mod a, 0 mod b)` has type
  `11 ∧ 11`, doubled. The case `−1 ∈ D_k(b)` is symmetric. If neither, `(1 mod a, −1 mod b)` has type
  `10 ∧ 01`, missing.
- A prime power `p^e`, `e ≥ 2`: `d = p` is missing, because a `k`-th power (`k ≥ 2`) of a non-unit has
  valuation at least 2 or vanishes mod `p^e`.
- A prime `p`: `D_k(p) ∖ {0}` is the subgroup `H` of index `g = gcd(k, p−1)`. A tournament needs
  `H ⊔ −H = F_p^×`, i.e. `g = 2` and `−1 ∉ H`, i.e. `p ≡ 3 (mod 4)`. At `p = 2`, `1 = −1` is doubled. ∎

Checked for `n < 3000`, `k = 2..8`. This is the clock form of guardrail 42 (MISTAKE-232): composite
moduli never carry Paley tournaments.

**Theorem 3 (orientation grading of character clocks; PROVED).** For distinct odd primes `p, q` let
`J(pq) = {d : (d/p)(d/q) = +1}` (the Jacobi clock).
- If `(−1/p)(−1/q) = +1`, `J(pq)` is symmetric: an undirected graph.
- Otherwise every unit difference is one-way (a tournament on unit differences), and non-unit
  differences are missing.

So two Paley tournaments (`p ≡ q ≡ 3 mod 4`) multiply to an undirected graph, and Paley times symmetric
is oriented. Under character products the orientation `ε(n) = (−1/n)` is multiplicative (a `Z/2`
grading). Under the AND law of Theorem 1, symmetry is an AND instead.

*Proof.* `((−d)/pq) = (−1/pq)(d/pq)`. ∎ Checked for all pairs `p, q < 60`.

**Proposition 4 (distance law; PROVED).** If `gcd(a, b) = 1`, `0 ∈ D_a` and `0 ∈ D_b`, then
`Cl(ab, D_a × D_b)` is the strong product of the loopless factor clocks. The `s`-fold sumset of `D_a × D_b` is `sD_a × sD_b`, so

```text
diam Cl(ab) = max(diam Cl(a), diam Cl(b)).
```

For a power clock the diameter is the local Waring number `g_k(n)`: the least `s` such that every residue
mod `n` is a sum of `s` `k`-th powers (0 allowed). So `g_k(n)` is the maximum over the prime powers of `n`.
Checked for `n < 600` (`k = 2, 3`) and `n < 400` (`k = 4`). The maxima are 4 (first at `n = 8`), 4 (first
at 9) and 15 (first at 16).

## 2. The fractal law at prime powers

**Theorem 5 (tower law; PROVED).** Let `p` be odd, `n = p^m`, `a = v_p(k)`, and write a nonzero difference
as `d = p^v w` with `p ∤ w` and `v < m`. Then `d ∈ D_k(p^m)` iff `k | v` and `w` is a `k`-th power unit
mod `p^(m−v)`. The latter is decided by `w mod p^min(m−v, 1+a)`. Hence:
- the type of `d` depends only on `v` and on the block of `min(m−v, 1+a)` base-`p` digits of `d` starting
  at position `v`;
- `d` is missing whenever `k ∤ v`.

So the clock mod `p^m` is a lexicographic tower. Read the digits from the least significant. Level `v` is
active iff `k | v`, and it is then decided by a block of `1 + a` digits; inactive levels are empty. The
partition into cosets mod `p^j` is a module partition iff no active block straddles position `j`.

- Squares mod `3^m` (`k = 2`, `a = 0`): `C3[K̄3[C3[K̄3[…]]]]` with loops. Every `3^j`-coset partition is a
  module partition (checked `m ≤ 5`).
- Cubes mod `3^m` (`k = 3`, `a = 1`): two-digit blocks, `C9[K̄3[C9[K̄3[…]]]]`, where `C9` is the undirected
  9-cycle (cubes mod 9 are `{0, ±1}`); a final one-digit active block is complete. The `3^j`-coset
  partition is a module partition iff `j ≢ 1 (mod 3)` (checked `m ≤ 6`).
- The general type statement is checked for `p = 3, 5, 7, 11`, `k = 2..6`, `p^m ≤ 2500`.

*Proof.*
- `v_p(x^k) = k·v_p(x)`; so for `v < m`, `d ≡ x^k` forces `x = p^s y` with `ks = v` and `w ≡ y^k`
  mod `p^(m−v)`, and conversely.
- For odd `p`, `(Z/p^j)^× = μ_(p−1) × (1 + pZ/p^j)`. If `k = p^a k'` with `p ∤ k'`, the `k`-th powers are
  `μ_(p−1)^k` (decided mod `p`) times `1 + p^(1+a)Z/p^j` (decided mod `p^(1+a)`).
- Modules: for `z` outside a coset `c + p^j Z`, the differences `m' − z` (`m'` in the coset) share `v` and
  all digits below position `j`. They have a common type iff the deciding block lies below `j`; if it
  straddles `j`, varying the digit at position `j` changes membership for suitable `c − z`. ∎

**Proposition 6 (digit conservation; PROVED).** For odd `p`, on units mod `p^m`, with `a = v_p(k)`:
- if `m ≥ a + 1`, `x ↦ x^k` depends only on `x mod p^(m−a)`: it forgets `a` input digits, and for
  `m ≥ a + 2` it really depends on the digit at position `m − a − 1`;
- if `m ≤ a + 1`, it depends only on `x mod p` (on nothing when `(p − 1) | k`);
- its image is decided mod `p^(1+a)`: it needs `a` extra digits.

For `p = 3`:
- `k = 2` (`a = 0`): squares mod 9 have period 9 and are read mod 3. The unit squares mod `3^m` are
  `{1 mod 3}`.
- `k = 3` (`a = 1`): cubes mod 9 have period 3 and are read mod 9. The unit cubes mod `3^m` are
  `{±1 mod 9}`.

These are the owner's two cycles, `0,1,4,0,7,7,0,4,1,0` (period 9) and `0,1,8,0` (period 3). Checked
`m ≤ 9` for `p = 3`, and the exact input dependence for `p = 3, 5, 7` and several `k` with `p^m ≤ 20000`.

**Theorem 7 (prime atoms).**

(a) PROVED. For prime `p`, every circulant `Cay(Z/p, D)` is a prime digraph (only trivial modules) unless
`D ∖ {0}` is empty or everything. In particular `Cl(p, D_k(p))` is prime unless `gcd(k, p−1) = 1`, in
which case it is complete.

*Proof.* By Gallai's modular decomposition theorem for digraphs, the maximal proper strong modules
partition the vertex set. The partition is invariant under automorphisms, and the quotient is complete,
linear or prime. The translations act transitively on `p` points, hence primitively. So the partition is
into singletons and the quotient is the clock itself. A vertex-transitive digraph on at least 3 vertices
is not linear. "Complete" means all pairs have one type, i.e. `D ∖ {0}` is empty or everything
(all-one-way is impossible). ∎ Checked for power clocks with `3 ≤ p < 200`, `k = 2..8`, and for 25 random
`D` at each prime `3 ≤ p < 60`. (At `p = 2`, `D ∖ {0}` is always empty or everything.)

(b) PROVED. Let `n = p^e m` with `gcd(p, m) = 1` and `m ≥ 2`. If `Cl(p^e, D_k(p^e))` has a clique module
`C` (`|C| ≥ 2`, every pair inside doubled), then `C × {b}` is a non-trivial module of `Cl(n, D_k(n))`, so
the composite clock is not prime.

*Proof.* Take `z = (z_1, z_2)` outside `C × {b}`. If `z_1 ∉ C`, the first factor's type is constant over
`C`. If `z_1 ∈ C`, then `z_2 ≠ b`, and the first factor's type is `11` (loop or doubled) against every
element of `C`. Either way the AND is constant. ∎

(c) FINITE-EXACT. Conversely, in all 1100 cases with `n ≤ 300`, at least two distinct primes and
`k = 2..6`, the clock is prime iff no prime-power factor has a clique module. The independent audit
repeated this with its own code for 8624 cases (`n ≤ 1500`, `k = 2..8`; not part of the committed
scripts) and found no exception. The other half of this equivalence is 7(b). A clique module may be the
whole factor: if `Cl(p^e, D_k(p^e))` is complete (for example `e = 1` and `gcd(k, p−1) = 1`), all of
`Z/p^e` is a clique module.

**Clique-module rule for odd prime powers.** For odd `p`, the clocks `Cl(p^e, D_k(p^e))` up to 300 with a
clique module are exactly those with `k | e − 1` and `gcd(k, p−1) = 1` (the audit confirmed this up to
1500; at `p = 2` the pattern is different).
- "If" is PROVED: the top coset `p^(e−1)Z/p^e` is then a clique module. Its differences `p^(e−1)w` are
  all doubled, and no active block straddles position `e − 1`, because consecutive active levels are at
  least `k ≥ p^a ≥ 1 + a` apart.
- "Only if" is FINITE-EXACT.

*Example.* Squares mod 63, 45, 225 and cubes mod 175, 189, 225 are prime, although each has a tower
factor (squares mod 9 and 25, cubes mod 25 and 27). CRT dissolves towers. Cubes mod 9 is not a tower: it
is the 9-cycle, which is prime (Theorem 5 with `m = 2` has no coset module), so cubes mod 63 involve no
tower at all.

So in TCPC:
- prime moduli give prime (indecomposable) digraphs;
- prime powers give towers whenever Theorem 5 allows a coset module (the fractal); the cube clock mod 9
  is too short for one and is prime;
- composites built by CRT are prime again unless a clique module survives.

## 3. Triplets

**Proposition 8 (the triplet principle is the transitive triangle; PROVED).** A 3-vertex tournament is
either the transitive triangle `TT3` or the cyclic triangle `C3`.
- `TT3` has no non-trivial automorphism and exactly one anti-automorphism (an isomorphism onto its
  converse). It swaps source and sink and fixes the middle vertex. On arcs, it fixes the closure arc
  source → sink and swaps the two path arcs.
- `C3` has the rotation group of order 3, so all three arcs are equivalent.

In a residue clock, the reflection `x ↦ C − x` is an anti-automorphism. On a transitive triangle
`0 → A → C` with closure arc `0 → C`, it fixes the closure arc and exchanges the two orders `A, C − A`
of the path. So "one member primarily differentiated, the other two differentiated only by orientation"
is exactly the `TT3` symmetry type. A `C3` has no distinguished member.

*Proof.* Direct. ∎

**Theorem 9 (Pythagorean triples are transitive triangles; loop placement; PROVED, classical content).**
A Pythagorean triple `a^2 + b^2 = c^2` is the triangle `0 → a^2 → c^2`, `0 → c^2` in the square clock
modulo any `n`. The hypotenuse is the closure arc and the legs are the two path arcs.
- At a prime `p ≡ 3 (mod 4)` with `p ∤ abc` (where the square clock is the Paley tournament), it is a
  `TT3`, never a `C3`. A `C3` would need `x + y + z ≡ 0` with `x, y, z` nonzero squares, i.e. the closure
  arc reversed.
- At `p ≡ 1 (mod 4)` the square clock is symmetric and the triangle is undirected.

(i) For odd `p`, the square clock mod `p` has a unit transitive triangle (`x + y = z`, all nonzero
squares) iff `p ≥ 7`. Hence a prime `p` divides `abc` for every primitive triple iff `p ∈ {2, 3, 5}`
(`p = 2` because `b = 2mn`). With `4 | b` (one of `m, n` is even), this gives the classical `60 | abc`.
*Proof.* For odd `p` the conic `X^2 + Y^2 = Z^2` has `p + 1` points over `F_p`, of which
`4 + 2[p ≡ 1 mod 4]` have `XYZ = 0`. So unit points exist iff `p ≥ 7`; and `(3, 4, 5)` avoids every
`p ≥ 7`.

(ii) For odd `p`, the hypotenuse can carry the loop (`p | c`) iff `−1` is a square mod `p` iff the square
clock mod `p` is symmetric (`p ≡ 1 mod 4`). At `p = 2` the clock is symmetric, but `c` is always odd.
- At `p = 3` the clock is the Paley `C3`, and the loop always sits on a leg.
- At `p = 5` the clock is symmetric, and the loop sits on each member equally often: the residue count of
  `(m, n) mod 5` is `8 : 8 : 8`. Observed 10607 (odd leg), 10617 (even leg), 10595 (hypotenuse) among the
  31819 primitive triples with `c ≤ 2·10^5`.
- At 2 the loop is always on the even leg (`4 | b`).

(iii) Every prime of `c` is `≡ 1 (mod 4)`.

Checked for all 31819 primitive triples with `c ≤ 2·10^5`, and for all `p < 400`.

**Theorem 10 (leg law: the even/odd differentiation is the quadratic character of 2; PROVED).** Write
primitive triples as `a = m^2 − n^2` (odd leg), `b = 2mn` (even leg), `c = m^2 + n^2`.

(i) For an odd prime `p`, let `N_p(x, y) = #{(m, n) mod p : (m^2 − n^2, 2mn) ≡ (x, y)}`. Then
`N_p(x, y) = N_p(y, x)` for all `x, y` iff `p ≢ 5 (mod 8)`. So modulo `p` the two legs are statistically
exchangeable exactly when `p ≢ 5 (mod 8)`.

(ii) For a prime `q ≡ 1 (mod 4)`:
- `q | a` implies `(b/q) = (2/q)`;
- `q | b` implies `(a/q) = +1`;
- `q | c` implies `(a/q) = (2/q)` and `(b/q) = +1`.

(iii) In Jacobi symbols, `(b/c) = +1` and `(a/c) = (2/c) = (2/a)` for every primitive triple.

So the even leg is a quadratic residue modulo every prime of the hypotenuse, while the odd leg is one
exactly at the primes `≡ 1 (mod 8)`. The asymmetry is carried by the factor 2 in `b = 2mn`. It is
invisible where 2 is a residue (`p ≡ ±1 mod 8`) or where `−1` is not (`p ≡ 3 mod 4`), and visible exactly
at `p ≡ 5 (mod 8)`.

*Proof.*
- (i) `N_p(x, y)` is the number of square roots of `x + iy` in `F_p[i] = F_p[X]/(X^2 + 1)`, and
  `(x, y) ↦ (y, x)` is `z ↦ i·z̄`. Conjugation preserves root counts; multiplication by `i` does iff `i` is
  a square in `F_p[i]`. For `p ≡ 3 (mod 4)`, `F_p[i] = F_(p^2)`, and `i` is a square because
  `8 | p^2 − 1`. For `p ≡ 1 (mod 4)`, `F_p[i] ≅ F_p × F_p` with `i ↦ (ρ, −ρ)`, `ρ^2 = −1`, and `i` is a
  square iff `(ρ/p) = (−1)^((p−1)/4) = +1` iff `p ≡ 1 (mod 8)`. For `p ≡ 5 (mod 8)`,
  `N_p(1, 0) = 4 ≠ 0 = N_p(0, 1)`.
- (ii) `q | a` gives `m ≡ ±n`, so `b ≡ ±2n^2`. `q | b` gives `q | m` or `q | n`, so `a ≡ m^2` or `−n^2`.
  `q | c` gives `m ≡ ρn`, so `a ≡ −2n^2` and `b ≡ 2ρn^2`; and `(ρ/q) = (2/q)` for `q ≡ 1 (mod 4)`.
- (iii) Multiply (ii) over the primes of `c`, all `≡ 1 (mod 4)`. Finally `16 | b^2` gives
  `c ≡ ±a (mod 8)`. ∎

*Examples.* `(3, 4, 5)`: `(4/5) = +1` and `(3/5) = −1 = (2/5)`. `(5, 12, 13)`: `5 | a`, and
`12 ≡ 2 (mod 5)` is a non-residue.

Checked: (i) for odd `p < 200`; (ii) for every prime `q ≡ 1 (mod 4)` dividing `abc` (137630 incidences)
and (iii) including `(2/c) = (2/a)`, for all primitive triples with `c ≤ 2·10^5`.

This is the precise sense in which "each leg is subtly differentiated from the other in an even/odd
halving/doubling manner". The reflection of Proposition 8 exchanges the legs as arcs. The arithmetic tells
them apart only through the 2 of `b = 2mn`: 2-adically (`4 | b`), and through `(2/q)` at
`q ≡ 5 (mod 8)`. Novelty is not checked; the statements are elementary.

## 4. F = S + U and p^2qr

Definitions follow DB1: `F`, `S`, `U` are the proper nontrivial divisors, the squarefree ones and the
prime ones. Put `Q = F ∖ S` (square-containing) and `M = S ∖ U`.

**Proposition 11 (PROVED).**
(i) `F = S + U` iff `|Q| = |U|`, since `|F| = |S| + |Q|`. By DB1 this holds iff `N ∈ {p, p^3, p^2qr}`.
Checked for `N < 2·10^4`.

(ii) At `N = p^2qr` the classes are `U = {p, q, r}`, `M = {pq, pr, qr, pqr}` and
`Q = {p^2, p^2q, p^2r}`: the owner's 3-4-3.

(iii) The five complementary pairs `d ↔ N/d` are
`(p, pqr), (q, p^2r), (r, p^2q), (p^2, qr), (pq, pr)`. On the three classes they form a multigraph:
- `U = Q` doubled;
- `U − M` and `M − Q` single;
- a loop at `M`, and no loop at `U` or `Q`.

This is an undirected three-vertex multigraph with a double edge and a loop. It uses the owner's
vocabulary of doubled and looped edges, but it is not a clock: complementation pairs are unordered, so
nothing is oriented, and "doubled" here means multiplicity 2, not two opposite arcs as in §0. At `N = 60`
the pairs are `(2,30), (3,20), (5,12), (4,15), (6,10)`.

(iv) The `p`-valuation layers `L_0 = {q, r, qr}`, `L_1 = {p, pq, pr, pqr}`, `L_2 = {p^2, p^2q, p^2r}` are a
second 3-4-3 split. Complementation swaps `L_0 ↔ L_2` and maps `L_1` to itself. The two splits agree on
`Q = L_2` and differ by the transposition `p ↔ qr` between their first two blocks.

(v) `Q = lcm(U, p^2)`: raising the repeated prime to its full power carries `U` bijectively onto `Q`.

(vi) In squareclass coordinates `F_2^3` (exponents of `p, q, r` mod 2), complementation is translation by
`[N] = [qr]`. Its orbits are `{0, [qr]}` and the three Fano lines through `[qr]` with `[qr]` removed. This
is the Fano sidecar of arithmetic_braids2.

(vii) For the profile `p^2 q_2 ⋯ q_r`, the classes have sizes `(r, 2^r − 1 − r, 2^(r−1) − 1)` and the
layers `(2^(r−1) − 1, 2^(r−1), 2^(r−1) − 1)`. They have the same shape iff `r = 3`, which is also the only
`r` with `F = S + U`. So 3-4-3 singles out `p^2qr` within this family. Checked `r ≤ 6`.

(viii) At `N = p^3`, `U = {p}` and `Q = {p^2}`, and complementation swaps them. At `N = p` everything is
empty.

**Reading (DICTIONARY, not a theorem).**
- The layer split of (iv) is a triplet in the sense of Proposition 8. Complementation is an involution
  that fixes the middle layer `L_1` (size 4, "primarily differentiated") and exchanges the outer layers
  `L_0` and `L_2` (size 3 each). The outer layers differ by doubling the `p`-exponent, `0 ↔ 2`.
- The class split is the same triplet only up to the transposition `p ↔ qr`.
- We read "`p^2qr` is `p` and `p^3` combined" against DB1's three solution shapes `{p, p^3, p^2qr}`. The
  `U ⇄ Q` exchange of `p^3` occurs twice inside `p^2qr` (`q ↔ p^2r`, `r ↔ p^2q`), while the remaining
  pairs (`p ↔ pqr`, `p^2 ↔ qr`, `pq ↔ pr`) involve `M`. Whether this is the intended reading is for the
  owner to say.
- The divisor poset of `p^2qr` is a finite box (DB1's mechanism), so by itself it carries no infinite
  self-similarity. The fractal behaviour available to TCPC is the prime-power tower of Theorem 5, and
  Theorem 7 says that a CRT product with another prime dissolves it.

## 5. Mod 9: the owner's two cycles are the Teichmüller split, and they drive Collatz

**Theorem 12 (PROVED, classical).** `(Z/9)^× ≅ C6 = μ2 × U1`, where `μ2 = {1, 8}` is the set of unit cubes
and `U1 = {1, 4, 7} = 1 + 3Z/9` is the set of unit squares. These are the two maximal subgroups. So `g`
generates `(Z/9)^×` iff `g` is neither a square nor a cube mod 9. Since `2 = (−1)·7` is such a `g`, 2
generates `(Z/3^k)^×` for every `k`.

For an odd prime `a`: 2 is a primitive root mod `a^2` iff 2 lies in no `q`-th power clock mod `a^2` for
the primes `q | a(a−1)`.
- The conditions with `q | a − 1` say that 2 is primitive mod `a` (Artin type).
- The condition `q = a` says that `2^(a−1) ≢ 1 (mod a^2)`, i.e. `a` is not a base-2 Wieferich prime,
  because the `a`-th power units mod `a^2` are the Teichmüller lifts.
- At `a = 3` these two conditions are exactly "2 is not a square mod 9" and "2 is not a cube mod 9".

Checked for `a < 400`. Among `a < 40`, `a = 7, 17, 23` fail through the square clock (`q = 2`) and
`a = 31` through the square and cube clocks (`q = 2, 3`, both Artin type). No `a < 400` fails through
the Wieferich clock: the first base-2 Wieferich prime is 1093.

So **the squares mod 9 are the Artin clock and the cubes mod 9 the Wieferich clock at the prime 3.**

**Theorem 13 (exact Syracuse clock formula; PROVED).** For odd `n`, let `T(n) = (3n + 1)/2^v` with
`v = v_2(3n + 1)`. Then

```text
T(n) ≡ (−1)^v · 4^(n+v)   (mod 9).
```

So Syracuse acts on the cube clock `{1, 8}` by the parity of `v`, and on the square clock
`{1, 4, 7} ≅ Z/3` by the rotation `n + v mod 3`. Consequences:
- `3n + 1 ≡ 4^n` is a unit square (digital root 1, 4 or 7);
- `T(n)^3 ≡ (−1)^v`;
- `T(n)` is a square mod 9 iff `v` is even;
- for `n` uniform among odd integers (`n mod 3` uniform, `v` geometric), the image law is `1/9` on each
  square and `2/9` on each non-square.

*Proof.* `4^n = (1 + 3)^n ≡ 1 + 3n (mod 9)`, and `2^(−1) ≡ 5 ≡ −4 (mod 9)`. ∎ Checked for all odd
`n < 2·10^6`, and the image law to `10^(−3)`.

**Theorem 13(b) (the `an + 1` maps mod `a^2`; PROVED).** For an odd prime `a` and odd `n`, with
`v = v_2(an + 1)`:

```text
(an + 1)/2^v ≡ ω(2)^(−v) · (1 + a)^(n + q_a(2)·v)   (mod a^2),
```

where `ω(2) = 2^a mod a^2` is the Teichmüller lift of 2 and `q_a(2) = (2^(a−1) − 1)/a` is the Fermat
quotient (only `q_a(2) mod a` matters). In the split `(Z/a^2)^× = μ_(a−1) × (1 + aZ/a^2)`:
- each halving turns the first factor by `ω(2)^(−1)`; this rotation generates iff 2 is a primitive root
  mod `a` (Artin's condition);
- each halving turns the second factor by `q_a(2)`; this rotation is frozen iff `a` is a base-2
  Wieferich prime.

Both rotations generate iff 2 is a primitive root mod `a^2`, which is the condition of THM-4523's
Theorem R_a (an iff). At `a = 3`, `ω(2) = −1` and `q_3(2) = 1`, and the formula is Theorem 13.

*Proof.* `an + 1 ≡ (1 + a)^n (mod a^2)` by the binomial theorem. Also `2 = ω(2)·2^(1−a)` and
`2^(1−a) = (1 + a·q_a(2))^(−1) ≡ (1 + a)^(−q_a(2)) (mod a^2)`. ∎ Checked for
`a = 3, 5, 7, 11, 13, 17, 1093, 3511` and all odd `n < 2·10^5`. The values of `q_a(2) mod a` are
1, 3, 2, 5, 3, 13, 0, 0: the rotation is frozen exactly at the two Wieferich primes.

**Theorem 14 (orbit law mod 9; PROVED).** For each `j ≥ 2`, the natural density of odd `n` with
`T^j(n) ≡ x (mod 9)` is `c_x/63`, where

```text
(c_1, c_2, c_4, c_5, c_7, c_8) = (8, 16, 11, 4, 2, 22).
```

*Proof.* The halving counts `v_1, …, v_j` of the first `j` steps are independent with
`P(v = t) = 2^(−t)` in natural density. Indeed, the odd `n` with a prescribed halving word form one
residue class mod `2^(v_1+⋯+v_j+1)` (induction on `j`; standard). Over all words of length `j` these
classes partition the odd integers, and their densities sum to 1. So for any union of them, the lower
densities of the union and of its complement sum to 1, and the natural density of the union exists and
is the sum of the class densities. For `i ≥ 2`,
`n_i ≡ (−1)^(v_(i−1)) (mod 3)`, so by Theorem 13, `n_(i+1) mod 9 = (−1)^(v_i) 4^(r + v_i)` with
`r = 1 + [v_(i−1) odd]`. Use `P(v odd) = 2/3` and `P(v ≡ t mod 6) = 2^(6−t)/63` (`t = 1..6`); the script
evaluates the sum in exact rationals. ∎ Observed for the second iterates of all odd `n < 2·10^6`, to
`10^(−3)`.

**Not new.** In Tao's framework (*Almost all orbits of the Collatz map attain almost bounded values*,
Forum Math. Pi 10 (2022) e12, eq. (1.22)), `T^j(n) mod 9` for `j ≥ 2` has the law of
`2^(−a_1) + 3·2^(−a_1−a_2) mod 9` with `a_i` independent geometric variables. Theorem 14 is the explicit
mod-9 table of that law (CITED framework; the literature check found the mod-3 law `(0, 1/3, 2/3)` printed
there, not the mod-9 table). Tao's Remark 1.15 stresses exactly this irregularity at coarse 3-adic
scales.

The image law and the orbit law agree on the square/non-square split (`1/3` against `2/3`). Inside each
class the orbit law is far from uniform: 7 mod 9 is rare (`2/63`) and 8 mod 9 is common (`22/63`). For
contrast, the greedy inverse map's stationary law is `2/9` on the squares and `1/9` on the non-squares
([three_adic_g_map](collatz_mod6_20260917_three_adic_g_map.md), Theorem 2.2).

**Proposition 15 (sliding window; PROVED).** For `i ≥ k`, `n_(i+1) mod 3^k` is a function of

```text
(v_(i−k+1) mod 2, v_(i−k+2) mod 6, …, v_i mod 2·3^(k−1)).
```

So the 3-adic clock reading of an orbit forgets the start after `k` steps, and its law at later times is
the push-forward of `k` independent geometric halving counts.

*Proof.* Induction on `k`: `n_(i+1) ≡ (3n_i + 1)2^(−v_i) (mod 3^k)`; `3n_i mod 3^k` depends on
`n_i mod 3^(k−1)`; and `2^(−v_i) mod 3^k` depends on `v_i` modulo `ord(2) = 2·3^(k−1)`. ∎ Checked for
`k ≤ 5` on all odd starts below `10^5` (60 steps each).

This push-forward is Tao's Syracuse random variable `Syrac(Z/3^k Z)`. His Proposition 1.14 (fine-scale
mixing) bounds its oscillation within cosets of `3^m` by `O_A(m^(−A))` (CITED): the deep 3-adic digits
are nearly uniform, while the coarse digits are not (Theorem 14). That mixing is the key input to his
theorem that almost all orbits attain almost bounded values.

**Duality (PROVED; checked `2 ≤ m ≤ 12`).**
- `⟨2⟩` mod `3^m` is all of `(Z/3^m)^×` and contains `−1`. So `Cay(Z/3^m, ⟨2⟩)` is symmetric: the complete
  3-partite graph on the residues mod 3.
- For `m ≥ 2`, `⟨3⟩` mod `2^(m+1)` is `{1, 3 mod 8}` and misses `−1`. So `Cay(Z/2^(m+1), ⟨3⟩)` is
  oriented: a bipartite tournament between even and odd residues. (Mod 4, `⟨3⟩ = {1, 3}` contains `−1`.)

Halving acts 3-adically as a doubled clock; tripling acts 2-adically as an oriented one.

**Where this enters the repository's Collatz theorems.**
- THM-4520 (`…gauss-twisted-circulant…`) uses that `{2^k}` runs through all units mod `3^n`. In clock
  terms: 2 avoids the square clock and the cube clock mod 9.
- THM-4523, Theorem R_a (`…functional-graph-is-rigid…`): the backward trees of `n/2, an + 1` separate the
  vertices prime to `a` iff 2 is a primitive root mod `a^2`. In clock terms (Theorem 12): 2 avoids the
  Artin clocks and the Wieferich clock mod `a^2`.

**Proposition 16 (the clock barrier; PROVED, classical).** For every finite sequence `v_1, …, v_k` of
positive integers with `L = Σ v_i`, the rational number

```text
n_1 = Σ_(i=1..k) 3^(k−i) 2^(v_1+⋯+v_(i−1)) / (2^L − 3^k)
```

has denominator prime to 6 and a periodic Syracuse orbit with exactly these halving counts (Böhm and
Sontacchi 1978; Lagarias, *The set of rational cycles for the 3x+1 problem*, Acta Arith. 56 (1990)
33–53; CITED via the literature check). Theorem 13, Proposition 15 and the one-step formula in the proof
of Theorem 14 hold verbatim on the odd elements of `Z_(6)`, reading `n mod 3` 3-adically (Theorem 14
itself is a density statement about integers). A cycle is determined by its halving word, and every word
has its rational cycle, whose elements have readings modulo every `2^a 3^b` just as integers do. So an
argument that uses only these readings would also exclude the rational cycles, which exist: it cannot
exclude a cycle. Integer instances are the negative cycles `{−1}`, `{−5, −7}` and the one through
`−17`. What the clocks at 2 and 3 cannot see is integrality: `2^L − 3^k` must divide the numerator. That
is a condition at the primes of `2^L − 3^k`, and the known exclusions of cycles control it through size,
i.e. how close `2^L` is to `3^k` (INFERENCE about the method; no theorem of this note).

*Proof.* Composing `k` steps gives `n_(k+1) = (3^k n_1 + Σ_i 3^(k−i) 2^(v_1+⋯+v_(i−1)))/2^L`. The numerator
sum is odd and `2^L − 3^k` is prime to 6. Applying the same formula to the cyclic shifts gives the other
cycle points, each with odd numerator, so `v_2(3n_i + 1) = v_i`. ∎ (`(1, 2) ↦ −5`, `(2) ↦ 1`, `(1) ↦ −1`.)

## 6. Open problems seen through clocks

### 6.1 Artin and Wieferich: the two halves of THM-4523's condition

By Theorem 12 and Theorem 13(b), "2 is a primitive root mod `a^2`" splits into Artin's condition (2
primitive mod `a`: the `μ_(a−1)` rotation generates) and the non-Wieferich condition (the `1 + aZ`
rotation `q_a(2)` is not frozen). Whether THM-4523's separation holds for infinitely many maps
`n/2, an + 1` is therefore entangled with two famous problems.
- Artin's conjecture for base 2 is OPEN unconditionally. Under GRH the density is `A = 0.3739558…`
  (Hooley, J. reine angew. Math. 225 (1967); CITED, CONDITIONAL). Heath-Brown (Quart. J. Math. 37 (1986))
  shows that at most two primes fail Artin's conjecture, but this does not decide the base 2 (CITED).
- Infinitely many non-Wieferich primes is OPEN. The abc conjecture implies `≫ log x` of them up to `x`
  (Silverman, J. Number Theory 30 (1988); CITED). Only 1093 and 3511 are known Wieferich primes, and the
  search is complete to `2^64` (PrimeGrid, 2022; CITED via OEIS A001220).

Even under GRH alone the intersection is not controlled, because nothing bounds the Wieferich primes
(INFERENCE). At `a = 3` both conditions are the owner's mod-9 cycles.

### 6.2 Fermat case I modulo `p^2`: triangles in the Teichmüller clock

The `p`-th power clock mod `p^2` is `{0} ∪ μ_(p−1)` (the Teichmüller lifts); at `p = 3` it is the owner's
cube clock `{0, 1, 8}`. A unit triangle `x^p + y^p ≡ z^p (mod p^2)`, `p ∤ xyz`, is a solution of
`ω(1) + ω(t) = ω(1 + t)`, i.e. a root `t ≢ 0, −1` of `((1 + t)^p − 1 − t^p)/p mod p` (checked against
direct powering for `p < 400`). If there is no such root, case I of Fermat's equation for `p` follows
from this congruence alone (reduce a hypothetical solution mod `p^2`).

**The clock is the whole `p`-adic obstruction (PROVED).** The unit `p`-th powers of `Z_p` are
`μ_(p−1) × (1 + p^2 Z_p)`, i.e. the units congruent to a Teichmüller lift mod `p^2`. So a unit triangle in
the clock mod `p^2` lifts to a solution of `x^p + y^p = z^p` in `Z_p` with `p ∤ xyz`, and conversely. Where
the clock has a triangle, no `p`-adic argument can prove case I.

- **Orbit law (PROVED; checked `5 ≤ p < 30000` in the main script and `5 ≤ p < 150000` in the census
  script, which test the invariance of the root set, the Eisenstein pair and "`t = 1` iff Wieferich"
  directly).** The roots are permuted by the `S3` generated by `t ↦ 1/t` and `t ↦ −1 − t`. So
  `#roots = 2[p ≡ 1 mod 3] + 3[t = 1 is a root] + 6j`.
  - The 2-orbit is the Eisenstein triangle `1 + ω = −ω^2`. It comes from the classical Cauchy–Liouville
    factorization `(X+1)^p − X^p − 1 = pX(X+1)(X^2+X+1)^e E_p(X)`, with `e = 2` for `p ≡ 1 (mod 6)` and
    `e = 1` for `p ≡ 5 (mod 6)` (CITED; treated in Ribenboim, *13 Lectures on Fermat's Last Theorem*
    (1979), which cites Klösgen 1970; we saw the contents, not the text).
  - The 3-orbit `{1, −2, −1/2}` is the isosceles triangle `1 + 1 = 2`. It lies in the clock iff
    `2^p ≡ 2 (mod p^2)`, i.e. iff `p` is a base-2 Wieferich prime: the owner's "doubling" appears here
    literally. Observed at exactly `p = 1093, 3511`.
- **Census (FINITE-EXACT).** For `5 ≤ p < 30000` and `p ≡ 2 (mod 3)`, the root counts are 0 (1371 primes),
  6 (243), 12 (18) and 18 (1). So 1371/1633 = 0.8396 of these primes have a triangle-free Teichmüller
  clock; the first are 5, 11, 17, 23, 29, 41, 47, 53, 71, 89, and the first with a triangle is 59. The
  literature check independently found 118 of the 725 primes `p ≡ 2 (mod 3)` below 12000 with a solution;
  we reproduce exactly that count. For `p ≡ 1 (mod 3)` the counts are 2 (1401), 5 (1, `p = 1093`),
  8 (192), 14 (15) and 17 (1, `p = 3511`). The census script extends the range to `p < 150000`.
- **Model (HEURISTIC).** There are about `p/6` generic `S3`-orbits of candidates. One representative
  decides its orbit, so if it is a root with probability about `1/p`, the number of root 6-orbits is about
  Poisson(1/6), and `P(no triangle) ≈ e^(−1/6) = 0.8465`. The census results for both residue classes are
  recorded in §6.6 against this model.
- **Open.** Whether infinitely many primes have a triangle-free Teichmüller clock is, as far as we know,
  open (UNCITED). It would give case I for infinitely many `p` from a single congruence mod `p^2`. That
  is a new route, not a new result: case I for infinitely many `p` is known (Adleman and Heath-Brown,
  with Fouvry, 1985; UNCITED-RECOLLECTION), and FLT itself is a theorem (Wiles 1995).

### 6.3 Waring and three cubes: diameters and antipodes

The cube clock mod 9 is the 9-cycle with loops. Its diameter is 4, attained at the antipodes 4 and 5, so
`±4 mod 9` are not sums of three cubes (checked). By Proposition 4, local Waring numbers are maxima over
prime powers: 4 for squares (at 8), 4 for cubes (at 9), and 15 for fourth powers (at 16; see the last
bullet for why the integer answer is 16).
- `G(3)`: the mod-9 clock gives `G(3) ≥ 4`. Linnik (Mat. Sb. 12(54) (1943)) gives `G(3) ≤ 7`, still the
  best bound; Siksek (Algebra & Number Theory 10 (2016)) made it explicit: every integer above 454 is a
  sum of at most seven positive cubes (CITED). Whether the clock diameter is the truth, `G(3) = 4`, is
  OPEN; Deshouillers, Hennecart and Landreau (Math. Comp. 69 (2000)) conjecture that 7373170279850 is the
  largest integer that is not a sum of four positive cubes (CITED).
- Three integer cubes: `k ≡ ±4 (mod 9)` has no representation. Heath-Brown (Math. Comp. 59 (1992))
  conjectures infinitely many representations for every other `k`. Booker (2019) found 33 and 795;
  Booker and Sutherland (PNAS 118 (2021)) found 42, 165, 579, 906 and a new representation of 3. The
  open cases `k ≤ 1000` are 114, 390, 627, 633, 732, 921 and 975 (Booker–Sutherland's list; no later
  solution was found in the literature check). All CITED.
- Four integer cubes: every `n ≢ ±4 (mod 9)` is a sum of four cubes (Dem'janenko 1966; CITED). Whether
  every integer is, i.e. whether the antipodes of the 9-cycle are reached with one more cube, is OPEN.
- `G(2) = 4` (Lagrange) equals the clock diameter mod 8. `G(4) = 16` (Davenport 1939;
  UNCITED-RECOLLECTION) is one more than the clock diameter 15 mod 16. The extra power is
  Hardy–Littlewood's `Γ(4) = 16`, which asks for non-singular local solutions (some term odd); in the
  integers it is forced by `31·16^m` (UNCITED-RECOLLECTION). Here the clock diameter is a lower bound that
  is not attained.

### 6.4 Doubling clocks: the boundary between provable and open

For an odd prime `p`, `Cay(Z/p, ⟨2⟩)` is
- complete iff 2 is a primitive root;
- the Paley tournament iff `p ≡ 7 (mod 8)` and `ord_p(2) = (p − 1)/2` (first: 7, 23, 47, 71, 79, 103, 167,
  191);
- symmetric iff `ord_p(2)` is even (`−1 ∈ ⟨2⟩`);
- otherwise oriented with missing pairs.

Census for odd `p < 2·10^6`: complete 0.37429, Paley 0.18694, symmetric non-complete 0.33377, oriented
with missing pairs 0.10500. Comparisons:
- complete against `A = 0.373956` (Artin; CONDITIONAL on GRH, Hooley 1967);
- Paley against `A/2 = 0.186978` (HEURISTIC: `p ≡ 7 mod 8` has density 1/4, and no odd prime `ℓ | p − 1`
  may have 2 as an `ℓ`-th power, which has density `Π_(ℓ odd)(1 − 1/(ℓ(ℓ−1))) = 2A`);
- symmetric in total, 0.70806, against `17/24 = 0.70833`. That density is an unconditional theorem of
  Hasse (Math. Ann. 166 (1966); CITED): an odd `p` divides some `2^n + 1` iff `ord_p(2)` is even.

For `p ≡ 7 (mod 8)`, index 2 is the same as "`−2` is a primitive root mod `p`", so the Paley doubling
primes are the primes `p ≡ 7 (mod 8)` with primitive root `−2`. The other index-2 primes (`p ≡ 1 mod 8`:
17, 41, 97, …) give the Paley graph, which is symmetric. Under GRH the density of all index-2 primes for
base 2 is `(3/4)A` in the near-primitive-root framework (Wagstaff 1982; Moree's survey, Integers 12A
(2012)); the value `(3/4)A` is the literature check's derivation, matching 0.2807 below `3·10^6`.

So the clock type itself stratifies known from open. Orientation is not read off `p − 1` alone (`p = 7`
and `p = 11` both have `v_2(p − 1) = 1`, but `ord_7(2) = 3` is odd and `ord_11(2) = 10` is even). It is
decided by whether 2 is a `2^(v_2(p−1))`-th power residue mod `p`: Frobenius data in the single Kummer
tower `Q(ζ_(2^k), 2^(1/2^k))` at the prime 2. The primes with `v_2(p − 1) ≥ k` have density `2^(1−k)`, so
the tower's tail is trivially small, and Hasse's density is unconditional. Completeness and the Paley
type involve every prime factor of `p − 1`; the density of complete clocks needs GRH, and the Paley
density is heuristic here. No unconditional infinitude result is known for either type in base 2. The
classification itself is checked against the actual pair types for `p < 3000`. The Wieferich prime 3511
is a Paley doubling prime (`ord = 1755`) and 1093 is not (`ord = 364`): NUMEROLOGY, recorded as data.

### 6.5 Inside a prime atom

Theorem 7 makes prime-level clocks the atoms of TCPC. Their own global combinatorics is not reached by
any combination law. The clique number of the Paley graph (squares mod `p ≡ 1 mod 4`) is at most
`(sqrt(2p − 1) + 1)/2` (Hanson–Petridis, Proc. LMS 122 (2021); CITED), still the best bound for primes,
and it is conjectured to be polylogarithmic: OPEN.

### 6.6 Literature check and the extended census

A literature subagent checked the cited statements on 2026-10-01 against journal pages, arXiv abstracts
and secondary sources. Primary texts were not all read: Ribenboim and Klösgen (§6.2) only through tables
of contents, and Lagarias 1990 through Brent's exposition (arXiv:math/0204170). Near misses recorded:
- A working assumption at the start of this session, that unit solutions of
  `x^p + y^p ≡ z^p (mod p^2)` exist for every `p ≥ 7`, is false. Our census and the literature check
  agree: they exist for every `p ≡ 1 (mod 6)` and only for a minority of `p ≡ 5 (mod 6)` (§6.2).
- The density `17/24` is Hasse's 1966 theorem. Lagarias 1985 is the Lucas-number analogue (`2/3`).
- Wirsching's monograph (LNM 1681) is not cited here: the check could not confirm that it uses the
  primitivity of 2 modulo powers of 3.
- Theorem 14 follows directly from Tao's framework and is not new.

**Extended census (FINITE-EXACT counts, HEURISTIC comparison).** The census script covers all primes
`5 ≤ p < 150000`. The orbit law holds with no exception, and `t = 1` is a root exactly at 1093 and 3511.
Here `j` is the number of generic root 6-orbits.

| Class | Primes | `j = 0, 1, 2, 3, 4` | Share `j = 0` | Mean `j` |
|---|---|---|---|---|
| `p ≡ 2 (mod 3)` | 6940 | 5842, 999, 95, 4, 0 | 0.8418 | 0.1731 |
| `p ≡ 1 (mod 3)` | 6906 | 5924, 903, 78, 0, 1 | 0.8578 | 0.1539 |
| Poisson(1/6) | | | 0.8465 | 0.1667 |

- For `p ≡ 2 (mod 3)` the share `j = 0` is exactly the share of triangle-free Teichmüller clocks, and it
  fits the model.
- Pooled over both classes the mean of `j` is 0.1635, 0.9 standard errors below 1/6.
- The classes differ by 2.8 standard errors: generic orbits are rarer when `p ≡ 1 (mod 3)`, where the
  Eisenstein triangle is always present. The deficit is not stable: in six ranges of width 25000 the
  `p ≡ 1 (mod 3)` mean is 0.1415, 0.1468, 0.1502, 0.1732, 0.1437, 0.1720, below 1/6 in four ranges and
  above it in two. We have no explanation, and with two classes compared it may be noise. It is recorded
  as an open observation, not a claim.
- `j ≥ 3` occurs at `p = 20123, 32909, 57149, 137201` (`j = 3`) and `p = 36847` (`j = 4`, so 26 roots).

## 7. Where the ungraspable structure lives

TCPC is exact on the combination layer: CRT (Theorem 1), towers (Theorem 5), atoms (Theorem 7) and
distances (Proposition 4). Every open problem met in this session sits at one of three places that these
laws cannot reach.
1. **Inside a prime atom:** the global combinatorics of one prime clock, such as the Paley clique number.
2. **Across infinitely many primes:** Artin (2 primitive mod `p`), Wieferich (2 a `p`-th power mod
   `p^2`), triangle-free Teichmüller clocks (case I mod `p^2`), Paley doubling primes. Each is a clock-type
   condition that is easy to evaluate at one prime and hard to control for infinitely many. Where the
   condition involves only the prime 2 (orientation of the doubling clock: one Kummer tower with a
   trivially small tail), the density is a theorem (`17/24`); where it involves every prime factor of
   `p − 1` (completeness, Paley), it is GRH-conditional or heuristic.
3. **From residues to integers:** Waring's `G(3)` against the cube clock's diameter 4; three cubes
   against the 9-cycle's antipodes; Collatz cycles against rational cycles (Proposition 16).

This is the honest answer to "how is TCPC useful for ungraspable structure". It places each problem in one
of these three layers, and it supplies the exact local data (clock types, densities, towers, orbit laws)
that any proof must respect. It also shows that the local data alone are realized by objects that the open
problems must exclude: rational cycles, `p`-adic solutions, residue-class sums.

## 8. Non-consequences

- Nothing here proves Collatz, Artin, Wieferich, the infinitude of triangle-free Teichmüller clocks,
  `G(3) = 4`, or any three-cubes case.
- The AND law, the towers and the mod-9 group structure are elementary and largely classical (CRT,
  Hensel, Teichmüller). The contribution is the uniform pair-type bookkeeping and the specific laws
  (Theorems 3, 7, 10, 13, 14). Their novelty is not established beyond §6.6.
- The 3-4-3 readings and "`p` and `p^3` combined" are a DICTIONARY. No map from divisor classes to clock
  vertices is claimed.
- The two 3-4-3 splits are equinumerous only at `r = 3` (Proposition 11(vii)). No general principle "any
  triplet has a TT3 structure" is claimed: Proposition 8 says which triplets do.

## 9. Directions

1. Prove Theorem 7(c), the converse prime-atom rule, probably through Gallai's theorem for categorical
   products of 2-structures.
2. Tower laws at `p = 2` (blocks of `2 + v_2(k)` digits) and for multiplicative subgroups other than
   `k`-th powers.
3. Exact orbit laws mod 27 and mod 81 (Proposition 15), compared with fine-scale mixing results for the
   3-adic Syracuse variable.
4. For `an + 1`: the `a`-adic analogue of Theorem 13 in Teichmüller coordinates mod `a^2`, its orbit law,
   and a comparison with THM-4523's separation exponents.
5. Teichmüller census beyond `1.5·10^5`: test Poisson(1/6) and whether the `p ≡ 1 (mod 3)` deficit of
   §6.6 persists.
6. Other families: prime clocks (residues of primes: local Goldbach diameter 2 against the open global
   problem), triangular numbers (Gauss's three-triangle theorem against the local diameter),
   Fibonacci and Lucas residues.

## 10. Reproduction

```text
python3 04-computation/experiments/tournament_clock_prime_collatz_20261001.py
python3 04-computation/experiments/tournament_clock_teichmuller_census_20261001.py
```

Standard library only. The main script takes about a minute and the census about 6 minutes. Each
output is saved beside its script and ends with `ALL CHECKS PASSED`. Sections A–H of the main script
match §§1–6 of this note; the census script is §6.2 and §6.6.

## 11. Audit

An independent audit subagent (2026-10-01) re-derived the statements with its own standard-library code
at commit `b7d759a4` and then re-read the final text. **Verdict: SOUND WITH CORRECTIONS.** No statement
labelled PROVED, FINITE-EXACT or HEURISTIC was false. Every correction below is applied in this version.

Errors, fixed:
1. §6.4 and §7 said that the orientation of the doubling clock depends only on the 2-part of `p − 1`.
   That is false (`p = 7` against `p = 11`). Orientation depends on whether 2 is a `2^(v_2(p−1))`-th power
   residue, i.e. on the single Kummer tower at the prime 2.
2. §0 gave the missing-neighbour count `n − 2|D| + |D ∩ −D|`, which needs `0 ∈ D` (the loopless `P7`
   would get 1 instead of 0). It is replaced by the general count, which the script now tests on random
   `D`.
3. The Theorem 7 example called the cube clock mod 9 a tower. It is the prime 9-cycle. Cubes mod 63 are
   replaced by cubes mod `189 = 27·7`, and the script now checks every modulus in the example.

Smaller corrections, fixed:
- Theorem 9: odd `p` in (i) and (ii), and `4 | b` for `60 | abc`.
- §6.3: `G(4) = 16` is one more than the clock diameter (Hardy–Littlewood's `Γ(4)`).
- Hypotheses: `m ≥ 2` in the duality, `gcd(a, b) = 1` in Proposition 4, the range of Proposition 6.
- Proposition 11(iii) is an undirected multigraph, not a clock.
- Proposition 16, the table and the loss ledger now say "modulo `2^a 3^b`". Only Theorem 13,
  Proposition 15 and the one-step formula transfer to `Z_(6)`.
- Theorem 7(c) is restructured: the two halves of the clique-module rule are separated, and a clique
  module may be the whole factor.
- Theorem 14: the step showing that the natural density exists is added.
- THM-4523's Theorem R_a is an iff, so 2-primitivity mod `a^2` is its condition, not its hypothesis.
- §6.2: case I for infinitely many `p` is already known.
- The `A/2` density is relabelled HEURISTIC.

Script checks strengthened (all pass):
- The hypotenuse law uses the actual triples.
- The doubling-clock classification is read off the actual pair types for `p < 3000`.
- Theorem 10 tests `(2/c) = (2/a)` and every prime `q ≡ 1 (mod 4)` dividing `abc`.
- Theorem 7(a) and the per-vertex counts are tested on random `D`.
- Proposition 6's exact input dependence is tested.
- Both Teichmüller scripts test the invariance of the root set, the Eisenstein pair and "`t = 1` iff
  Wieferich" directly, instead of inferring the Wieferich term from `t = 1`.

Reproduced independently by the audit (its own code, not committed):
- the main script's output, byte for byte;
- Theorems 1–3 on wider ranges;
- Theorem 5 on 220864 differences and 1273 module cases (`p ≤ 11`, `k ≤ 81`);
- Theorem 7(b, c) on 8624 composite cases (`n ≤ 1500`, `k = 2..8`);
- Theorem 9 through the Berggren tree (31819 triples, the same 10607/10617/10595 split);
- Theorem 10 on all 137630 incidences;
- Theorem 12 for `a < 3000`;
- Theorem 13 for odd `|n| < 2·10^6`, negative `n` included, and 13(b) for all odd primes `a < 200` and
  random 100-bit `n`;
- Theorem 14 from Tao's formula. For `j = 1` the law is `(7, 14, 7, 14, 7, 14)/63`, so `j ≥ 2` is sharp;
- the Teichmüller census to 150000 with a different algorithm, matching every count in §6.2 and §6.6;
- the index-2 share 0.28067 below `3·10^6`, against `(3/4)A = 0.28047`.

Not checked by the audit: novelty, and the citations labelled UNCITED-RECOLLECTION (Davenport, `Γ(4)`,
Adleman–Heath-Brown and Fouvry).
