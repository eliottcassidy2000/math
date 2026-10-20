---
id: THM-4730
title: "The strong pure-blue alphabet is complete: a strongly connected tournament all of whose Hamiltonian paths are palindromic (a pure-blue node = a point of symmetry of the merged metagraph) is K1, C3 or T5. (A) Rotation lemma: for a pure-blue tournament and any Hamiltonian cycle (v_0 ... v_{n-1}), every reflection v_x -> v_{2i-1-x} is an anti-automorphism and the rotation by two steps v_x -> v_{x+2} is an automorphism; for odd n the rotation itself is an automorphism, so a strong pure-blue tournament of odd order is a circulant on Z_n with 1 in S (with respect to EVERY Hamiltonian cycle). (B) No circulant tournament of odd order n >= 7 is pure-blue: with 1 in S, either two distinct a, b in S\{1} have a + b <= n + 1 and the three-interval path (0..a-2, a+b-2..n-1, a-1..a+b-3) is a non-palindromic Hamiltonian path, or S = {1} u [(n+1)/2, n-2] and the descending-pairs path (0,1, n-2,n-1, n-4,n-3, ..., 3,4, 2) is one. (C) No strong pure-blue tournament has even order: for n = 0 mod 4 the two-step rotation has even order and a tournament has no involutory automorphism; for n = 2 mod 4 the tournament is a Z_{n/2}-bicirculant (S', D) whose mixing Hamiltonian cycle (0, 2s, ..., -2s, d-2s, d, ..., d-4s) (s = +-1 in S', d in D with d - 4s not in D, which exists because D is not a union of even cosets) has a two-step rotation carrying an even vertex to an odd one, making Aut vertex-transitive, impossible for even n. (D) Hence THM-2454's alphabet {1, C3, T5} is complete for all n, THM-643's conjecture C3 holds (pure-blue fibres are 1 or 3), the Paley tournament T_p is not pure-blue for any prime p >= 7, and pure-blue(n) = c(n/2) (n even), c((n-1)/2) + c((n-3)/2) + c((n-5)/2) (n odd), c = A000930, is PROVED for every n <= 125 (THM-2477's rigid classification through 62 supplies the halves) and for all n conditional on THM-2477's lemma ML"
status: >
  AUDIT OVERRIDE (2026-10-20, weighted-core session): Theorem B's implication
  from nonconstant partner sums to a non-palindromic path is invalid for
  general circulants; it tests only affine reflections. The all-order
  odd-order proof and dependent alphabet/count consequences are OPEN here
  pending a direct anti-automorphism argument. The specific path claims
  remain unrefuted. Theorem A and the even-order argument C are unaffected.
  The following records the original pre-audit claims, not their promotion:
  PROVED (pre-audit): Theorem A (rotation lemma), Theorem B (odd circulants of order >= 7 are not pure-blue;
  the two explicit Hamiltonian paths), Theorem C (no strong pure-blue tournament of even order), Corollary D
  (alphabet completeness; pure-blue count law for n <= 125; fibre law; T_p not pure-blue for p >= 7).
  MACHINE-CHECKED (atom_alphabet_proofcheck.out): the proof's explicit paths are Hamiltonian and non-palindromic for
  every connection set with 1 in S at every odd n = 7..23 (2,044 sets: 2,035 three-interval, 9 descending-pairs);
  the mixing cycle C' is a Hamiltonian cycle with an even-even arc in every bicirculant (S', D) at n = 6..30
  (1,198,372 pairs); brute force confirms no pure-blue bicirculant at n = 6, 10.
  FINITE-EXACT (independent of the proof): no strong pure-blue tournament of order 6, 7, 8 (gentourng -c census,
  pos_exact.out: 35, 353, 6008 strong classes); among circulant tournaments of odd order 5..17 only T5 is
  pure-blue (circulant_pos_fast.out, affine anti-automorphisms, exact for prime orders by Burnside's prime-degree
  theorem; pos_exact.out redoes the composite orders 9 and 15 with all involutory anti-automorphisms, 27 and 135
  of them for the two Z_3 x Z_3-type sets); the Z_3 x Z_3 Cayley tournaments at n = 9 are not pure-blue
  (palindromic_paths.out).
  CITED: THM-643 (pure-blue, blue parities), THM-644(a) (grid-symmetric iff positional reversal is an
  anti-automorphism), LEM-003 / THM-2477 LEM-003+ (Aut acts freely on Hamiltonian paths), THM-2444 (rigidity lemma;
  T5 = (15,5,3)), THM-2454 (center law and the alphabet law (2), census n <= 13), THM-2477 (rigid strong = C3 for
  2 <= m <= 62; lemma ML open beyond), Camion 1959 (strong iff Hamiltonian cycle), Redei 1934, odd order of
  tournament automorphism groups, Burnside's theorem on transitive groups of prime degree.
  NOT PROVED here: lemma ML of THM-2477 (rigid strong tournaments of order >= 63), hence the count law beyond
  n = 125 is conditional; any statement about pure-blue classes beyond the alphabet (the mixed/pure-black
  structure is THM-4729's subject).
  Terminology: pure-blue node (THM-643) = point of symmetry (the owner's name) = a class all of whose tilings
  are grid-symmetric = a tournament all of whose Hamiltonian paths are palindromic. T5 = the regular
  tournament of order 5 = the circulant {1, 2} on Z_5.
source: mac-mini-2026-10-20-metagraph (owner's prompt: comprehensive patterns of the merged metagraph, points of symmetry, forced arrangements, new terminology)
depends_on:
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-644-gridsym-fiber-law-and-anti-redei.md
  - 01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md
  - 01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md
related:
  - 01-canon/theorems/THM-2444-pure-blue-count-refutation-and-rigid-class-law.md
  - 01-canon/theorems/THM-2453-palindromic-narayana-law-for-the-rigid-stratum.md
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
scripts:
  - 04-computation/experiments/metagraph_20261020/atom_alphabet_proofcheck.py
  - 04-computation/experiments/metagraph_20261020/pos_exact.py
  - 04-computation/experiments/metagraph_20261020/circulant_pos_fast.py
  - 04-computation/experiments/metagraph_20261020/palindromic_paths.py
outputs:
  - 04-computation/experiments/metagraph_20261020/atom_alphabet_proofcheck.out
  - 04-computation/experiments/metagraph_20261020/pos_exact.out
  - 04-computation/experiments/metagraph_20261020/circulant_pos_fast.out
  - 04-computation/experiments/metagraph_20261020/palindromic_paths.out
---

# THM-4730 -- strong points of symmetry are K1, C3 and T5

**AUDIT OVERRIDE, 2026-10-20: odd-order proof gap.** Nonconstant reversal-partner
sums exclude affine reflection, but do not exclude a nonaffine positional
anti-automorphism. In `Circ(Z_9,{1,3,4,7})`, the Hamiltonian path
`(0,1,2,3,7,5,6,4,8)` is palindromic although its partner sums are `{5,8}`.
The original companion's odd-order `palindromic` test used the weaker
criterion. [Exact witness and audit](../../05-knowledge/results/weighted_core_transport_20261020/inheritance.md)
preserve this correction. The special three-interval and descending-pairs
claims are not refuted; they need an actual arc argument. Theorem B and
the all-order alphabet consequences below are therefore provisional.
The rotation lemma A, even-order argument C, and properly stated finite
censuses are unaffected. Original pre-audit prose is retained for lineage.

**PROVED (pre-audit) + MACHINE-CHECKED + FINITE-EXACT** as itemized in the status.

## 0. Setting and the question

Tiling model (01-canon/definitions.md): base path `n -> n-1 -> ... -> 1`, tiles `(a,b)`
with `a >= b+2`, `m = C(n-1,2)` tiles; the grid involution `sigma` (THM-643) is the
relabelling `i -> n+1-i`. By THM-644(a) a tiling `t` is grid-symmetric iff `i -> n+1-i`
is an anti-automorphism of `T_t`, i.e. iff the base path read backwards is carried onto
itself by an anti-automorphism. Call a Hamiltonian path `(v_1, ..., v_n)` of a
tournament `T` **palindromic** when the positional reversal `v_i -> v_{n+1-i}` is an
anti-automorphism of `T`. The tilings of the class of `T` are the Hamiltonian paths of `T`
up to `Aut(T)` (LEM-003), and a tiling is blue (grid-symmetric) iff its path is
palindromic. So:

```text
C is pure-blue (THM-643)  <=>  every Hamiltonian path of C is palindromic
                          <=>  every line at the node of C is blue (the owner's "point of symmetry").
```

THM-2454 proved the center law: a pure-blue class is a palindromic stack
`(S_1, ..., S_k)` of strong components with `S_{k+1-i} = S_i^rev`, every half part rigid
(`H = |Aut|`, one tiling), and the center (if `k` is odd) a STRONG pure-blue tournament;
through `n = 13` the strong pure-blue atoms were `{1, C3, T5}` and completeness was open
(`n = 15` the first open odd size, every even size via centers). This theorem closes it
for all `n`.

## 1. Theorem A (rotation lemma)

Let `T` be pure-blue and `C = (v_0, v_1, ..., v_{n-1})` a Hamiltonian cycle of `T`
(indices mod `n`). For `i = 0, ..., n-1` the rotation `(v_i, v_{i+1}, ..., v_{i-1})` is a
Hamiltonian path, hence palindromic: its positional reversal is

```text
sigma_i : v_x -> v_{2i-1-x}          (an anti-automorphism of T).            (A1)
```

Composing two anti-automorphisms gives an automorphism:

```text
sigma_i o sigma_j : v_x -> v_{x + 2(i-j)}  in Aut(T)  for all i, j.          (A2)
```

So the two-step rotation `rho^2 : v_x -> v_{x+2}` is an automorphism. If `n` is odd, `2`
is invertible mod `n`, so `rho = (rho^2)^{(n+1)/2}` is an automorphism: `T` is invariant
under the cyclic shift along `C`, i.e. `T` is the circulant tournament on `Z_n` with
connection set `S = {v_y - v_0 : v_0 -> v_y}` in the coordinates `v_x = x`, and `1 in S`
because `v_x -> v_{x+1}` are the arcs of `C`. This holds for every Hamiltonian cycle of
`T`: a strong pure-blue tournament of odd order is circulant with respect to each of its
Hamiltonian cycles. (Camion: a strong tournament has a Hamiltonian cycle.) `#`

## 2. Theorem B (odd order): no circulant tournament of odd order n >= 7 is pure-blue

Let `T = Circ(Z_n, S)`, `n >= 7` odd, and (relabelling by `x -> -x` if necessary) `1 in S`.
For each `d in [2, (n-1)/2]` exactly one of `d, n-d` lies in `S`; call `d in S` a
*small* element, `n - d in S` a *large* one. Write `arc(u,v)` for `v - u in S`.

**Three-interval paths.** For `x, y >= 1` with `x + y <= n - 1` put `X = [0, x-1]`,
`Y = [x, x+y-1]`, `Z = [x+y, n-1]` (consecutive intervals of `Z_n`) and let

```text
P(x,y) = (0, 1, ..., x-1,  x+y, x+y+1, ..., n-1,  x, x+1, ..., x+y-1)   (X, then Z, then Y).
```

Inside each interval the steps are `+1 in S`. The two jumps are `(x-1) -> (x+y)` of
difference `y + 1` and `(n-1) -> x` of difference `x + 1`. So `P(x,y)` is a Hamiltonian
path of `T` whenever `x + 1, y + 1 in S`. Its palindromic test: `v_q = q - 1` for
`q <= x` and `v_{n+1-q} = x + y - q` for `q <= y`, so `v_q + v_{n+1-q} = x + y - 1` for
`q <= min(x,y)`; the middle block `Z` is read increasingly, and if `x = y` its positions
are symmetric and `v_{x+i} + v_{n+1-x-i} = (x+y+i-1) + (n-i) = x + y - 1 (mod n)`: the path
IS palindromic when `x = y`. If `x < y` then position `x+1` carries `x + y` (the first
vertex of `Z`) and position `n - x` lies in `Y`, carrying `v_{n-x} = x + y - (x+1) = y - 1`;
their sum `x + 2y - 1` differs from `x + y - 1` by `y`, which is not `0 (mod n)` since
`1 <= y <= n-2`. Symmetrically for `x > y`. Hence

```text
a, b in S \ {1},  a != b,  a + b <= n + 1   ==>  P(a-1, b-1) is a NON-palindromic Hamiltonian path.   (B1)
```

**Which S admit (B1).** If `S\{1}` has two small elements `a < b <= (n-1)/2`, then
`a + b <= n - 1` and (B1) applies. If it has exactly one small element `a`: for `a >= 3`
take `d = a - 1 in [2, (n-1)/2]`, `d != a`, so `n - d in S` is large and
`a + (n - d) = n + 1`; for `a = 2` take `d = 3 <= (n-1)/2` (this is where `n >= 7` is
used), `n - 3 in S`, `2 + (n-3) = n - 1`. Either way (B1) applies. The only remaining
case is **no small element**:

```text
S = {1} u {n - d : 2 <= d <= (n-1)/2} = {1} u [(n+1)/2, n-2]       (x -> x+1 and x -> x-d, 2 <= d <= (n-1)/2).
```

**Descending-pairs path.** For this `S` and odd `n >= 7` let

```text
Q = (0, 1,  n-2, n-1,  n-4, n-3,  ...,  3, 4,  2).
```

The runs `(n-2j, n-2j+1)`, `j = 1, ..., (n-3)/2`, use the step `+1`; the jumps
`1 -> n-2` and `(n-2j+1) -> (n-2j-2)` have difference `-3 in S` (as `3 <= (n-1)/2`), and
the last jump `4 -> 2` has difference `-2 in S`. The runs cover `{3, ..., n-1}` once each
(their left ends are the odd numbers `n-2, n-4, ..., 3`), so `Q` is a Hamiltonian path.
It is not palindromic: `v_1 + v_n = 0 + 2 = 2` but `v_2 + v_{n-1} = 1 + 4 = 5`, and
`5 != 2 (mod n)` for `n >= 7`. `#`

Together: every circulant tournament of odd order `n >= 7` has a non-palindromic
Hamiltonian path, so by Theorem A no strong pure-blue tournament has odd order `>= 7`.
The strong tournaments of odd order `<= 5` which are circulant are `K1`, `C3` and the
two connection sets `{1,2}`, `{1,3}` on `Z_5`, both `T5`; `T5` is pure-blue (all `15`
Hamiltonian paths palindromic: THM-2444's `(15, 5, 3)`, re-verified in pos_exact.out
and circulant_pos_fast.out).

## 3. Theorem C (even order): no strong pure-blue tournament has even order

Let `T` be strong, pure-blue, of even order `n >= 4`, with Hamiltonian cycle `C`
labelled `v_x = x in Z_n`. By (A2) the even translations `x -> x + 2j` are automorphisms
and by (A1) the reflections `x -> c - x`, `c` odd, are anti-automorphisms.

**n = 0 (mod 4).** `x -> x + n/2` is an automorphism of order `2`. A permutation of order
`2` has a 2-cycle `(u v)`, and an automorphism swapping `u` and `v` would reverse the arc
between them: impossible. `#`

**n = 2 (mod 4), k = n/2 odd.** Let `E` and `O` be the even and odd residues; the even
translations act transitively on each. The arcs inside `E` are `e -> e'` iff
`(e' - e)/2 in S'` for a set `S' subset Z_k` with exactly one of `y, -y` in `S'` for each
`y != 0` (a circulant tournament `R = Circ(Z_k, S')`); applying the anti-automorphism
`x -> 1 - x` shows the arcs inside `O` are governed by the SAME `S'`. The arcs between
`E` and `O`: `e -> o` iff `o - e in D` for a set `D` of odd residues, and then
`o -> e` iff `o - e not in D`. The cycle `C` gives `1 in D` (arcs `e -> e+1`) and
`-1 not in D` (arcs `o -> o+1`). Let `s in {1, -1}` be the element of `S'` among `+-1`.

*Claim: there is `d in D` with `d - 4s not in D`.* Otherwise `D - 4s = D`, so `D` is a
union of cosets of the subgroup generated by `4s = +-4`, which is the subgroup `2Z_n` of
all even residues because `gcd(4, n) = 2`; then `1 in D` forces `-1 in D`. Contradiction.

With such `d` consider

```text
C' = (0, 2s, 4s, ..., (k-1)2s = -2s,   d - 2s, d, d + 2s, ..., d - 4s).
```

The first half lists `E` (steps `2s`, allowed since `s in S'`), the second half lists `O`
(steps `2s`, same `S'`); the junction `-2s -> d - 2s` is an arc because
`(d - 2s) - (-2s) = d in D`; the closing step `(d - 4s) -> 0` is an arc because
`(d - 4s) - 0 not in D`. So `C'` is a Hamiltonian cycle of `T` containing the arc
`0 -> 2s` between two even vertices. By (A2) applied to `C'`, the two-step rotation
along `C'` is an automorphism; it carries the vertex two steps before `d - 2s`, namely
`-4s in E`, to `d - 2s in O`. Hence `Aut(T)` is transitive on `E u O = V(T)`, so all
scores are equal to `(n-1)/2`, which is not an integer for even `n`. Contradiction. `#`

(For `n = 6` and `n = 10` the whole family `(S', D)` was also brute-forced: no member has
all Hamiltonian paths palindromic.)

## 4. Corollary D

1. **Alphabet completeness.** The strong pure-blue tournaments are exactly `K1`, `C3`
   and `T5`. Hence (THM-2454, center law) the pure-blue classes of order `n` are the
   palindromic stacks whose half parts are rigid strong tournaments and whose center is
   `none`, `K1`, `C3` or `T5`.
2. **The count law.** With THM-2477 (a rigid strong tournament of order `2 <= m <= 62` is
   `C3`) the half parts of any pure-blue class of order `n <= 125` are `K1` or `C3`, so
   `pure-blue(n) = c(n/2)` for even `n` and `c((n-1)/2) + c((n-3)/2) + c((n-5)/2)` for
   odd `n`, `c = A000930` (Narayana's cows), holds for every `n <= 125`; for all `n` it is
   equivalent to THM-2477's lemma ML.
3. **Fibres (THM-643 C3 PROVED).** By THM-2454's factorization `tc = prod tc_i^2 * tc(center)`
   the fibre of a pure-blue class is `tc(center) in {1, 3}`, and
   `(H, |Aut|) = (3^a, 3^a)` or `(15 * 3^a, 5 * 3^a)`, `a = #C3 parts`.
4. **Rotational family.** The Paley tournament `T_p` (and every circulant of odd order
   `>= 7`) is not pure-blue; THM-2454's rotational-purity question is settled at the
   blue level, extending THM-2477's `tc(T_p) > 1`.
5. **Both faces of a strong point of symmetry.** For `C3` and `T5` every Hamiltonian
   cycle is a rotation class of an arithmetic progression (`T5` has exactly `2`
   Hamiltonian cycles, `C3` one), consistent with Theorem A's demand that the shift
   along each Hamiltonian cycle be an automorphism (`Aut T5 = Z_5`, `Aut C3 = Z_3`).

## 5. Terminology fixed by this theorem

- **palindromic path**: a Hamiltonian path whose positional reversal is an
  anti-automorphism (= a blue tiling of its class).
- **point of symmetry** (owner) = **pure-blue node** (THM-643): all paths palindromic.
- **strong point of symmetry** = atom of the alphabet: `K1, C3, T5` only.
- **mixing cycle**: a Hamiltonian cycle with two consecutive vertices of equal parity
  with respect to a reference cycle (the object that kills even orders).

## 6. Files

- `04-computation/experiments/metagraph_20261020/atom_alphabet_proofcheck.py` / `.out`:
  the proof's paths and cycles checked at odd `n = 7..23` and `n = 6, 10, ..., 30`.
- `pos_exact.py` / `.out`: strong classes of order 6, 7, 8 (none pure-blue); composite
  circulants `n = 9, 15` with every involutory anti-automorphism.
- `circulant_pos_fast.py` / `.out`: circulants of odd order `5..17`, only `T5` pure-blue.
- `palindromic_paths.py` / `.out`: circulants `n <= 13`, `Z_3 x Z_3` Cayley tournaments,
  all strong 7-tournaments.
