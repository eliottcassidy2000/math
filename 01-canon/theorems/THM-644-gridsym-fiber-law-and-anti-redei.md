# THM-644 — The grid-symmetry fiber law and the Anti-Rédei parity

**Source:** opus-2026-07-07-S139 (fiber-allocation census + targeted verification).
**Status:** part (a) **PROVED** (two lines); part (b) **VERIFIED-exact** on all 72 classes
n = 4,5,6 with proof sketch via LEM-003 (formalizable); part (c) **CONJECTURE** (the new
parity theorem), verified all classes n ≤ 6 and via the all-odd g-spectra at n = 7;
part (d) **PROVED** (one line from Rédei + odd automorphism groups).
**Companions:** mac-mini-S46 (THM-643 reserve, line-side parity/node-type formulas) and
klein-S161 (line structure) — this file is the FIBER side ("which tilings map to which
class"); strict blue/black definitions per CLAUDE.md.
**Scripts:** `metagraph_fiber_allocation_opus_S139.py`, `metagraph_anti_redei_test_opus_S139.py`
(+ .out in 05-knowledge/results).

Setting: base path `n → n−1 → ⋯ → 1`; tiles `(x,y)`, `x ≥ y+2`, `m = C(n−1,2)`; grid
involution `σ(x,y) = (n−y+1, n−x+1)` with `f = ⌊(n−1)/2⌋` fixed tiles; `flip` = complement
all tiles (a LINE = a `{t, flip(t)}` pair; BLUE iff gridsym). For a merged class `C`:
`N(C)` = #tilings, `g(C)` = #gridsym tilings, `b = N − g`.

## (a) Grid symmetry IS the reversal anti-automorphism (PROVED)

> A tiling `t` is grid-symmetric **iff** `ρ(i) = n+1−i` is an **anti-automorphism** of `T_t`
> (`T(ρu, ρv) = T(v,u)`).

*Proof.* The base path is ρ-anti-symmetric (`ρ` maps the base arc `k → k−1` to the pair
`(n+1−k, n+2−k)`, whose anti-image is the base arc `n+2−k → n+1−k`). For a tile `(x,y)`:
`T(ρx, ρy) = T(y,x)` unpacks, via `{ρx, ρy} = σ(x,y)`, to `bit(σ(x,y)) = bit(x,y)`. ∎

Corollary (upgrades CLAUDE.md's verified fact to a theorem): `g(C) > 0 ⟹ C` admits an
anti-automorphism ⟹ `C` is transpose-self. **Pure-black ⊇ non-transpose-self**; the census
converse (no transpose-self class is pure black, n ≤ 7) is part (c)'s consequence.

## (b) The fiber law (VERIFIED 72/72; the anti-symmetric LEM-003)

> `g(C) · |Aut(C)| = H_anti(C) := #{ Hamiltonian paths P = (p₁…p_n) of T :
> (p_j ↦ p_{n+1−j}) ∈ AntiAut(T) }.`

Verified exactly on every class at n = 4, 5, 6 (with LEM-003 `N·|Aut| = H` re-checked in the
same run). Proof shape: LEM-003's free action of `Aut` on (tournament, Hamiltonian-path)
pairs restricts to the ρ-equivariant locus — a gridsym tiling in `C` is exactly a pair
`(T, P)` with `P`'s positional reversal an anti-automorphism, and `Aut` acts freely on these.

## (c) THE ANTI-RÉDEI CONJECTURE (new parity theorem candidate)

> For every tournament `T` on `n` vertices admitting an anti-automorphism,
> **`H_anti(T)` is odd.** (In particular ≥ 1: every self-converse class contains an
> anti-symmetric Hamiltonian path — equivalently, no transpose-self class is pure black.)

Verified: all transpose-self classes, n = 4, 5, 6 (H_anti odd, every class), and at n = 7
the full g-spectrum `{1:10, 3:12, 5:16, 7:32, 9:18}` is all-odd (via (b), `H_anti/|Aut|`
odd, and `|Aut|` odd). This is the σ-equivariant refinement of Rédei's theorem — precisely
the shape THM-587's frame (Rédei is σ-odd, lonely is σ-even) called for.

### (c′) The reduction structure (opus-S140 — proved pieces + two verified sub-lemmas)

**PROVED:**
1. *(Partition.)* For any Hamiltonian path `P`, the positional reversal `α_P` is an
   involution with `n mod 2` fixed points; a path contributes to `H_anti` iff `α_P` is an
   anti-involution, so `H_anti = Σ_β h_β` over **anti-involutions** β,
   `h_β = #{P : α_P = β}`. (An anti-involution has ≤ 1 fixed point: two fixed points u,v
   would force `T(u,v) = T(v,u)`.)
2. *(Lemma A — oddly many β.)* If `AntiAut(T) ≠ ∅` it is a coset `α·Aut(T)` of odd size;
   inversion `x ↦ x⁻¹` is an involution of this coset whose fixed points are exactly the
   anti-involutions; an involution on an odd-cardinality set has an odd number of fixed
   points. **#anti-involutions is odd.** ∎
3. *(Center/mirror structure.)* In a β-symmetric Hamiltonian path (n = 2h): an internal
   pair-arc `(a, βa)` can occur **only as the center arc** (its positional mirror is itself:
   any other occurrence would repeat a vertex), and every cross arc is used **jointly with
   its β-tie partner** at mirrored positions. Hence flipping an internal pair bit changes
   `h_β` by exactly `W_{βa} − W_a`, where `W_x` = #(one-representative-per-other-pair
   T-paths whose last vertex beats `x`).
4. *(Base case.)* The transitive tournament admits ρ and has `h_ρ = 1`.
5. *(Fold algebra.)* β-symmetric tournaments with fixed β form a GF(2) cube (internal bits +
   tied cross bits `s = T(a,b) ⟷ T(βb,βa)`, `t = T(a,βb) ⟷ T(b,βa)`); the bit
   `q(A,B) = s ⊕ t` is representative-independent (the even-graph-layer invariant), and
   over GF(2): `W_a + W_{βa} ≡ Σ_Q q(pair(last(Q)), A)`.

**Hence (proved reduction):** Anti-Rédei follows from Lemma A + *(Main Lemma)* `h_β` odd for
each anti-involution β; and the Main Lemma follows, via cube-connectivity + the base case,
from two GF(2) sub-lemmas:
- **(M) Mirror parity:** `W_a ≡ W_{βa} (mod 2)` — equivalently the q-weighted rep-path count
  `Σ_Q q(pair(last Q), A)` is even. *Verified 1832/1832 random instances (n = 4,6,8).*
- **(S/T) Tie-flip parity:** flipping a tied cross pair changes `h_β` evenly — equivalently
  `h_β` of the two β-anti-symmetric block-contraction digraphs agree mod 2 (contraction of
  `[a,b]` with mirror `[βb,βa]` stays anti-symmetric — proved). *Verified 1200/1200.*

All flip types tested parity-safe (3032/3032); the Main Lemma itself verified on every
self-converse class n ≤ 6 (per-β!) with spectra `{1,3}, {1,3}, {1,3,5,7,9}` and on 400
random symmetric tournaments at n = 5, 7. **Status: Anti-Rédei = two verified GF(2)
path-parity sub-lemmas on the fold.** Scripts: `anti_redei_decomposition_opus_S140.py`,
`anti_redei_flip_parity_opus_S140.py`.

## (d) The parity law of the fiber allocation (PROVED / structured)

> **`N(C)` is odd for every class** (*proof:* `N = H/|Aut|`, `H` odd by **Rédei**, `|Aut|`
> odd since tournament automorphism groups have odd order). Assuming (c):
> **`g(C)` is odd on transpose-self classes and 0 otherwise**; hence
> `b(C) = N − g` is **even exactly on the transpose-self (blue-touched) classes** and odd
> (= N) on pure-black ones.

This is the precise form of the owner's "blues contribute odd amounts and blacks even":
per class, the blue (gridsym) count is odd wherever nonzero; the black count is even
exactly there.

## Line-level closed forms and small structure (verified n = 4..7)

- `#gridsym = 2^{(m+f)/2}`; **blue lines** `= 2^{(m+f)/2 − 1}`; **black lines**
  `= 2^{m−1} − 2^{(m+f)/2 − 1}` (exact at n = 4..7).
- Node types (pure-blue / mixed / pure-black): `1/1/2`, `3/5/4`, `2/10/44`, `4/84/368`
  (n = 4..7); pure-black = non-transpose-self count exactly; the pure-blue vs mixed split is
  `H_anti = H` vs `H_anti < H`.
- **Blue class-self lines vanish at odd n** (0 at n = 5, 7; 1, 2 at n = 4, 6) — conjecture:
  parity alternation (the ½-shift/zero-halo flavor).
- The flip-partner class map is genuinely non-constant on most classes (8/12 at n=5,
  50/56 at n=6, 447/456 at n=7): the line-multigraph on classes is rich; its full
  distribution is the remaining descriptive layer (fiber data files in results/).

## The restricting picture (assembled)

`N(C) = H/|Aut|` (LEM-003) says **how many** tilings hit each node; part (a)+(b) say
**which ones are blue** (`H_anti/|Aut|`, the anti-symmetric layer); Rédei + odd-Aut force
`N` odd; Anti-Rédei forces `g` odd-or-zero; the line pairing `{t, flip(t)}` with these
parities constrains the allocation class-by-class. Two parity theorems — the classical
Rédei and its σ-equivariant refinement — govern the entire blue/black fiber structure.
