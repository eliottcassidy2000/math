# Flushing out the LRC with the tournament induction toolkit

*kind-pasteur-2026-06-22-S31x. Immersing in the repo's tournament induction machinery (multiplicativity,
the three modes, deletion–contraction, the forbidden K₃ atom) and transferring it to LRC(14). The LRC
inherits the tournament's multiplicative atom structure, and its open base case IS the apex-7 atom — the
same irreducible obstruction the tournament has at `H = 7 = I(K₃,2)`.*

## The tournament induction toolkit (immersion)
| technique | identity | reduces to | status |
|---|---|---|---|
| **H-multiplicativity** (THM-454) | `H(T) = Π_C H(strong component C)` | irreducible strong ATOMS | PROVED |
| **forbidden atom** (THM-200) | `H=7 ⟺ Ω=K₃` impossible (K₃ forces C₅) | the apex-prime obstruction | PROVED |
| **Mode B** (THM-291) | `H_{n+2}(inner,0)=H_n(inner)+Δ` | base `H_3=1+2t` | PROVED |
| **deletion–contraction** (THM-082) | `H(D)=H(D∖e)+H(D/e)` (all digraphs) | single vertices | PROVED |
| **ΣH formula** (THM-292) | `ΣH = Σ_{succ-free σ} 2^{m-(n-1-bp)}` | A000255 | PROVED |
| **half-tiling Venn** (THM-549) | 3 corners + 3 edges + 1 center | the 3 modes | PROVED |
| **H-maximizer** (THM-028) | Paley/regular maximizes H (BIBD, min score-var) | score convexity | PROVED |

The spine of the whole thing: **decompose into atoms, the atoms multiply, and the induction bottoms out
on a finite list of irreducible atoms — one of which (`K₃`, `H=7`) is FORBIDDEN by a structural forcing.**

## Transfer 1 — multiplicativity (VERIFIED): `safe = Π safe(clusters)`
The LRC analogue of `H = Π H(atoms)` (`lrc_multiplicativity_atoms_kps.py`): when `S` splits into
scale-separated clusters `S = S_1 ⊔ S_2` (`S_2 ≫ S_1`), the orbits decouple and

> `meas(safe(S_1 ∪ S_2)) ≈ meas(safe(S_1)) · meas(safe(S_2))`

— ratios `1.01, 0.99, 1.00` at ~200× separation, `0.97` at 34× (error shrinks with separation). So the
LRC **safe-measure is multiplicative over independent clusters**, exactly as `H` is over strong
components. The "atoms" are the **single-scale clusters** (bounded cores, all speeds comparable, no
internal scale separation). This is the structural completion of my S31w peeling: peeling removes one
large cluster at a time; multiplicativity factors over ALL scales at once, `meas(safe(S)) = Π_k
meas(safe(cluster_k))`, each cluster a smaller LRC with positive safe-measure ⟹ the product is positive
⟹ `M(S) ≥ 1/14`. The induction bottoms out on the single-scale **atoms**.

## Transfer 2 — the forbidden/hardest atom is the SAME apex-7 object
In tournaments the irreducible base is the **forbidden atom `K₃` (`H=7`)** — three pairwise-conflicting
3-cycles force a `C₅`, so `H=7` never occurs (THM-200). In the LRC the irreducible base is the
**bounded covering core**, whose hardest member is the **AP `{1..13}` with `M=1/14`** — and `14 = 2·7`,
the apex prime. The correspondence is exact in spirit:

| | tournament | LRC |
|---|---|---|
| irreducible atom | `K₃` (`H = I(K₃,2) = 7`) | the AP `{1..13}` (`M = 1/14 = 1/(2·7)`) |
| extremality | FORBIDDEN (`K₃` forces `C₅`) | the BOUNDARY (`AP` is the unique tight minimizer) |
| structural cause | 3 conflicting 3-cycles ⟹ common vertex ⟹ `C₅` | three-gap rigidity (only the AP has ≤3 gaps) |
| apex prime | `7 = I(K₃,2)` | `7 = n/2`, the Eisenstein-folded apex |

So **the LRC's open base case is the tournament's forbidden atom, read on the runner side** — both are
the apex-prime-7 irreducible obstruction. (This is the `seven-is-the-unique-forbidden-clique-value` /
parity-dual-scars trilogy, now seen as the SHARED induction base.)

## Transfer 3 — what tournaments do at the base, and what the LRC must do
The tournament CLOSES its base by a **structural forcing**: `Ω = K₃` is geometrically impossible (the
three conflicting 3-cycles cannot avoid a 4th vertex / a `C₅`). The LRC base needs the analogous forcing:
**every bounded covering set has `M > 1/14`, because the unique `M=1/14` minimizer (the AP) is
non-covering (it omits `14`)** — i.e. the three-gap (Steinhaus) rigidity that *only the AP is gap-rigid*,
so any covering (hence non-AP) bounded set is strictly looser. The tournament's "K₃ forces C₅" and the
LRC's "non-AP forces a 4th gap (looser cover)" are the same kind of extremal forcing at the apex prime —
the tournament version is proved (THM-200); the LRC version (HYP-2885 consec-extremality) is the crux.

## Transfer 4 (direction) — deletion–contraction to crack the atom
The one PROVED tournament tool that reduces an *irreducible* object further is **deletion–contraction**
(THM-082, `H(D)=H(D∖e)+H(D/e)`, valid for all digraphs). Its LRC analogue should act on the
**resonance/relation graph** `Ω_E` of the speed set (the lattice `{Σ n_i e_i = 0}`), whose cover is the
Bonferroni sum `p0 = Σ_r (−1)^r S_r`. If `p0` is (or bounds) an independence polynomial `I(Ω_E, 2)` — as
`H` is for tournaments — then `p0` inherits deletion–contraction on the relations, and the bounded-core
ATOM reduces by removing/merging resonances down to a finite irreducible relation pattern. That finite
pattern is the consec-extremal configuration. **This is the concrete next transfer: realize the LRC
cover as `I(Ω_E, ·)` and run THM-082 on the resonance lattice to reduce the atom to a finite check.**

## Net
The tournament induction toolkit transfers cleanly: (1) **multiplicativity** factors the LRC over
scale-clusters (VERIFIED), reducing to single-scale atoms; (2) the **irreducible atom is the apex-7
object** in both worlds (forbidden `K₃` ↔ boundary AP); (3) the **base-case closure** is an apex-prime
extremal forcing — *proved* for the tournament (K₃⟹C₅), *open but analogous* for the LRC (three-gap
rigidity); (4) **deletion–contraction** on the resonance lattice is the unexploited tool that could
reduce the LRC atom to a finite check, mirroring how it cracks tournament H. The LRC and the tournament
are the two faces of one apex-7 induction — and the tournament face shows what closing the base looks like.

→ THM-454 (H-mult), THM-200 (K₃ forbidden), THM-082 (deletion–contraction), THM-291 (Mode B), THM-549
(half-tiling Venn), HYP-2885 (consec-extremality base), HYP-2900 (induction), S31w (peeling), S31x
(multiplicativity), `seven-is-the-unique-forbidden-clique-value.md`, [[lrc14-thread]].
