# The three recursion modes are parity-stratified, and LRC(14) is Eisenstein(even) ∘ Legendre(odd) on the Möbius skeleton

*kind-pasteur-2026-06-22-S31r. The owner's clarification: the Möbius mode applies to ALL sizes, the
Legendre mode only to ODD sizes, the Eisenstein mode only to EVEN sizes — and the subtournament sizes
behind each letter differ per mode. Composing them in this parity-stratified way builds exactly the
structure LRC(14) lives in. Made explicit and comprehensive.*

## The domain: the half-tiling = the complement quotient
All three recursions live on the **half-tiling** `h(n) = ⌊(n−1)²/4⌋` — the fundamental domain of the
complement involution `σ(x,y)=(n+1−y,n+1−x) = φ(T^op)` (THM-549). Because every complement-invariant
invariant (`c_3`, `HP`, `H`, the OCF, and the LRC cover/floor — all verified complement-invariant)
is computable on this half, the three modes are recursions for *the objects we actually care about*,
not just for tile counts.

## The three modes, explicitly (sizes verified, parity verified)
`lrc_three_modes_parity_composition_kps.py` confirms each claim.

| mode | sign word | char poly / order | holds for | half-tiling shape | reduction sizes (the letters) |
|---|---|---|---|---|---|
| **Möbius** | `A+B+C−D−E−F+G` `+++---+` | incl–excl lattice of {A,B,C} | **ALL n** | — (the skeleton) | `{n−1, n−2, n−3}` = Modes A,B,C + overlaps |
| **Legendre** | `A+B−C+D−E−F+G` `++-+--+` | `(x−1)³(x+1)`, order 4 | **ALL n, but the ONLY form for ODD** | **SQUARE `k²`** (3 corners) | `{n−1×2, n−2×2, n−3×2, n−4}` — 3-set Venn (see correction) |
| **Eisenstein** | `A+B−C` `++-` | `(x−1)²`, order 2 | **EVEN n ONLY** | **PRONIC `k(k−1)`** (no corner) | `{n−1, n−2}` |

VERIFIED: the order-2 Eisenstein recurrence `h(n)=2h(n−1)−h(n−2)` **fails at every odd `n`**
(`5→3, 7→8, 9→15, …`) and holds at every even `n` — exactly "Eisenstein only for evens." The
order-4 Legendre recurrence holds for all `n` but is the *minimal* form an odd half-tiling admits
(odd = square = 3-corner; even = pronic = degenerate, gets the extra order-2). So the **sign words and
the subtournament sizes are different objects per mode** — `{n−1,n−2}` vs `{n−1,n−3,n−4}` vs
`{n−1,n−2,n−3}`. The Möbius `+++---+` is the inclusion–exclusion over the three reduction modes
themselves (Mode A `n→n−1`, Mode B `n→n−2`, Mode C `n→n−3`; `7=2³−1` subsets, 3+3+1, HYP-2689).

## CORRECTION (owner, S31s): the Legendre recursion is a 3-set inclusion–exclusion (Venn)
The Legendre letters are NOT a sparse `{n−1,n−3,n−4}` — they are a full **3-set Venn** over three
sub-tilings, with sizes `A,B = n−1`, `C,D = n−2`, `E,F = n−3`, `G = n−4` (`lrc_legendre_venn_correction_kps.py`):
- **3 corners** (the generating sets): `A, B` (depth `n−1`, the `ω,ω²` conjugate pair) and `D` (depth
  `n−2`, the lone generator — the `S₃`/Eisenstein cube-root split, HYP-2689).
- **3 edges** (pairwise *unions* `|X∪Y|`): `A+B−C`, `A+D−E`, `B+D−F` ⟹ the pairwise *intersections* are
  `A∩B = C` (`n−2`), `A∩D = E` (`n−3`), `B∩D = F` (`n−3`).
- **center** (full union `|A∪B∪D|`): `A+B+D−C−E−F+G = h(n)`.

VERIFIED: the 7 Venn regions **partition** the odd square (`n=7`: `9 = 1+1+1 (corners) + 3+1+1 (edges)
+ 1 (center)`), and `A+B+D−C−E−F+G = h(n)` for all `n`. So **all four depths `n−1..n−4` are present**;
the net coefficients `(2,0,−2,1)` arise because the `n−2` corner `D` (+) and the `n−2` overlap `C` (−)
**cancel** — same size, geometrically distinct (corner vs edge-overlap, THM-549). My earlier
"Legendre skips `n−2`" was the net-coefficient view, not the geometry.

Two consequences:
1. **The even (Eisenstein) mode IS the `A∩B` edge of the odd Venn.** `A+B−C = |A∪B|`, the conjugate-pair
   union, with the lone corner `D` and the triple `G` **folded away** by the complement symmetry (the
   `(x+1)` factor) — pronic, not square. So Eisenstein (even) = the *degenerate* Legendre, the conjugate
   pair without the lone generator.
2. **This Venn IS the LRC coverage recursion (THM-548).** The 3 sub-tilings ↔ the **3 far runners**
   `{u,v,w}`; corners/edges/center ↔ **one-far / two-far / three-far** Newton packets `Δ_S` (`S⊆{u,v,w}`),
   `Σ_S Δ_S = p0(B∪{u,v,w})−p0(B)`. So the half-tiling Legendre Venn and the LRC three-far cover are the
   *same* 3-set inclusion–exclusion — the "relevant structure": **Möbius is the inclusion–exclusion
   principle; Legendre (odd, full square) and Eisenstein (even, degenerate edge) are its two parity
   realizations, and the LRC cover is that same Venn with the far runners as the three sets.**

## Why these three characters, by parity
- **Möbius (always)** is the combinatorial skeleton — the subset lattice of the three reduction
  depths. It is the `μ`/inclusion–exclusion that gives `φ = μ*id` and the coprime-density floor
  (`1/ζ(2)`, kps-S31q). It cannot "not apply": any depth-≤3 local recursion obeys it.
- **Legendre (odd)** is the quadratic character. An odd half-tiling is a **perfect square `k²`** with a
  genuine 3-corner (S₃) decomposition; the doubling `z→2z` of order 3 splits `F_p^*` into
  `QR={1,2,4}` vs `NQR={3,5,6}` (the `χ_7` of `2S/U=−43−7χ_7(a)`). Squares and quadratic residues are
  the same parity fact: **odd ⇒ square ⇒ Legendre**.
- **Eisenstein (even)** is the complement-fold. An even half-tiling is a **pronic rectangle `k(k−1)`**
  with no 3-corner; its `(x−1)²` order-2 recurrence is the binary/`Z[i]` Cayley–Dickson step, the
  literal "÷2." **even ⇒ pronic ⇒ Eisenstein-fold.**

## The composition IS the structure we care about: 14 = 2·7
For `LRC(14)`:
```
   14 (EVEN)  --Eisenstein-->  half-tiling = PRONIC 7·6 = 42  ==>  exposes k = 7 = 14/2  (the apex prime)
    7 (ODD)   --Legendre--->   half-tiling = SQUARE 3² = 9    ==>  the χ_7 character, QR{1,2,4}/NQR{3,5,6}
```
The even Eisenstein mode **factors out the 2** — the complement `Z/2` is *the same involution* as the
LRC sector-halving `x→−x` that takes `14→7` (THM-280/549, kps-S31i) — and the resulting pronic
`k(k−1)` has `k = n/2 = 7`, **handing the apex prime to the odd Legendre mode**. So the arithmetic
factorization `14 = 2·7` is *literally* the operator composition **Eisenstein(even) ∘ Legendre(odd)**,
sitting on the Möbius skeleton. Parity alternates down the whole tower `14e, 13o, 12e, 11o, …`, each
size taking its own mode; the LRC only needs the top fold (even, ÷2) and the apex (odd, χ).

## What it gives for proof and disproof
This is the **structural source of the repo's q-uniform proof** (LRC(2q) for all q): every even case
`2q` is Eisenstein-folded to its odd apex `q`, where the Legendre `χ_q` character sits on the Möbius
totient floor `→ 3/π² = 1/(2ζ(2))`.

> **LRC(2q) = Eisenstein(even: 2q→q, ÷2) ∘ Legendre(odd: q, χ_q) ∘ Möbius(always: φ-floor > 0).**

- **Proof side:** the three modes are the three *necessary* steps, and each is benign — the even fold
  is exact (complement-invariance), the Möbius floor is a positive totient sum, and the Legendre `χ_q`
  only *biases* the floor (`osc/floor < 1` q-uniformly, kps-S31q) — never flips it. The composition
  closes for every q; the apex prime is the only place the bias is largest (smallest q = LRC(6)
  tightest).
- **Disproof side:** a counterexample at `2q` must survive the even fold to the apex `q` and then beat
  the Legendre-biased Möbius floor — but the bias cannot overcome the principal, so it is pushed onto
  the cap margin (OPEN-Q-108) at the apex, exactly where `q` prime is hardest. **The parity
  composition explains why `14 = 2·7` is the first open case: it is the smallest even fold onto the
  unique permanent-gap apex prime `7` (the `H=7`/`K_3`-forbidden value).**

## Net (the comprehensive picture)
The owner's three modes are not three competing recursions — they are a **parity-stratified
factorization** of one structure: Möbius (the always-on inclusion–exclusion skeleton / coprime
floor), Legendre (the odd-size square / quadratic-residue character), Eisenstein (the even-size
pronic / complement-fold). Their subtournament sizes differ exactly because square vs pronic vs
lattice are different geometries. Composing them is the arithmetic factorization of `n`:
`14 = 2 · 7 = Eisenstein(even) ∘ Legendre(odd)`, on the Möbius skeleton — which is precisely the
LRC(14) reduction to 7 sectors and the apex-7 `χ_7` floor. The difficulty of `n` is the depth of its
even folds onto a defective odd apex prime.

→ THM-549 (half-tiling parity recursions), HYP-2689 (3 modes incl–excl), THM-291 (Mode B), THM-280
(reflection = complement), HYP-2882 (φ=μ*id totient bridge), HYP-2856 (3/π² floor),
`three-modes-three-characters-mobius-legendre-eisenstein-and-the-lonely-runner-floor.md`,
`the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner.md`, [[lrc14-thread]].
