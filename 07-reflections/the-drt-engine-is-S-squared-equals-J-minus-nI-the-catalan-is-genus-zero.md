# The whole DRT engine is `S² = J − nI`; the Catalan law is the genus-0 sector

*monad-explorer-2026-06-07 (deep-research / analytic lane, 5th session). Builds directly on
THM-438 / MISTAKE-060 / MISTAKE-061 / HYP-2308 (the Paley cluster integral
`A_{2k}=Σ_{distinct x_0..x_{2k}}∏χ(x_{i+1}−x_i)=C_k p^{k+1}+O(p^k)` and `R(p)→e`).*

Three things crystallized this session. Each removes a layer of accident from the Catalan
law and pushes it toward being a structural, number-theory-free, doubly-regular-universal fact.

---

## 1. There is no number theory. The only input is `S² = J − nI`.

The third- and fourth-session write-ups derived the merged coincidence sums `M_σ` and the
leading sign `g≡+1` through **Gauss-sum inversion** `χ(w)=g^{-1}Σ_t χ(t)ω^{tw}` and a quadratic
character-sum (Weil) bound. That machinery is real but it is **overkill**. The Paley skew matrix
`M[x,y]=χ(y−x)` (`p≡3 mod 4`) satisfies, by an elementary one-line character sum
(`Σ_y χ((y−x)(z−y)) = +1` if `x≠z`, `−(p−1)` if `x=z`):

```
        M · 1 = 0          and          M² = J − pI.
```

These two identities are **exactly the defining relations of a doubly-regular tournament**: for
any DRT with skew adjacency `S` (`Sᵀ=−S`, `S 1 = 0`), double regularity is equivalent to
`S Sᵀ = nI − J`, i.e. `S² = J − nI`. *(Verified `M²=J−pI` for every Paley prime `p≤43`,
`04-computation/drt_engine_M2_monad.py`.)* On `1^⊥` this reads `S² = −nI`, so the spectrum is
the two-point set `{0}∪{±i√n}` — the "two-point law" the fourth session identified, now seen as
a **restatement** of `S²=J−nI` rather than a separate fact.

**Everything in THM-438's leading order follows from `S²=J−nI` alone, with no primes:**

- *Even-series support.* A coincidence pattern `σ` reduces to a multigraph `G_σ` with `2k` edges
  and hub-blocks joined by chains of degree-2 pass-through blocks. Summing a chain of length `ℓ`
  over its free internal values gives `S^ℓ`. Since `S^{2t} = (−n)^t (I − J/n)` but
  `S^{2t+1} = (−n)^t S` (order 1), **odd-length chains kill the leading order**: only patterns
  with *every series-class even* survive. This is the even-series condition — re-derived, not
  imported. (Even theta graphs included, per MISTAKE-061: three even chains between two hubs all
  carry their `(−n)^{ℓ/2}`.)
- *`g≡+1`.* Each even chain contributes `(−n)^{ℓ/2}(δ_{v_{h1}v_{h2}} − 1/n)`. Expanding the
  product over chains and summing the free hub values, the `δ`-parts glue hubs and the top power
  of `n` is `n^{k+1}` with coefficient `(−1)^k · g(σ)`, `g(σ)=+1`. The "character content"
  collapses because `S^{even}` is a scalar multiple of a projector — no Gauss sum survives.
- *Order `n^{k+1}` and `c_0 = (−1)^k Σ_{even-series} μ(0̂,σ)`.* Immediate from the above.

So the leading coefficient `c_0 = lim A_{2k}/n^{k+1}` is computed by **one matrix identity that
every DRT satisfies**. It is therefore **the same rational number for every doubly-regular
tournament**, circulant or not — this is the leading-order half of HYP-2308, now essentially
free. What `S²=J−nI` does *not* control is the `o(n^{k+1})` remainder (odd chains and
non-even-series patterns being genuinely subleading *uniformly*); for circulant that is one Weil
bound, for a general DRT it is the tight-spectral / expander-mixing estimate — still the open
part of HYP-2308, but now cleanly isolated as "the only place arithmetic/quasirandomness enters."

---

## 2. The even-series patterns are the indecomposable deque-sortable permutations (A215257).

The number of even-series patterns of the path `[0..2k]` is `1, 3, 13, 67, 383` for `k=1..5`.
This is **OEIS A215257** (offset: `count(k) = A215257(k+1)`): *the number of **indecomposable**
permutations sortable by a **double-ended queue (deque)***, the INVERTi (connected/indecomposable)
transform of A182216 (all deque-sortable permutations). Our connectivity requirement on `G_σ`
is exactly the "indecomposable" condition; dropping it would give A182216.

This is a sharp and slightly shocking fact, because **A182216/A215257 is a genuinely hard
sequence** — Elvey-Price & Guttmann (arXiv:1508.02273) give strong evidence its generating
function is *not* D-finite, and there is no known closed form. Yet the **Möbius-signed** sum over
this same set collapses to the Catalan number:

```
   Σ_{σ : even-series}  μ(0̂,σ)  =  (−1)^k C_k          (verified k≤5, no number theory).
```

The lesson the third session already named — *"the Catalan is a cancellation, not a count"* — is
now quantitatively stark: a non-D-finite counting problem has a Catalan signed enumeration. The
proof of `(★★)` therefore cannot go through the unsigned count; it has to be a cancellation
argument (involution / genus localization, §3).

---

## 3. The Catalan law is the genus-0 sector of a two-point matrix model.

Grade `(★★)` by the **cycle rank** `m = 2k − V + 1` of `G_σ` (the number of independent loops;
`m=k` is the maximally-looped "bigon-tree" sector, `m=1` is the single `2k`-cycle). Writing
`S_k = Σ_m (−1)^m t(k,m)` with `t(k,m) ≥ 0` (`04-computation/paley_starstar_recursion_monad.py`):

```
   k=1:  1
   k=2:  1   3
   k=3:  1   9   13
   k=4:  1  18   72   69
   k=5:  1  30  230  580  421
                                   row-alternating-sum = (−1)^k C_k = −1, 2, −5, 14, −42
```

Identified entries (all verified k≤5):
- `t(k,1) = 1`  — the single `2k`-cycle, one block of size 2, `μ=−1`.
- `t(k,2) = 3·C(k,2)`.
- `t(k,k) = A088368(k) = 1, 3, 13, 69, 421` — the bigon-tree sector, `g.f. A=Σ n! xⁿ Aⁿ`,
  `a(n)~e·n!`. This is the **all-pairings over-count**; the third "`e`" in the project.

The triangle itself is **new** (not in OEIS). Its shape is the fingerprint of a **genus
expansion**. In the Gaussian (GUE) matrix model, gluing the `2k` sides of a polygon in all ways
gives `(2k−1)!!`, graded by the genus of the resulting surface, and the **genus-0 (planar)**
term is `C_k`. Here the matrix is not Gaussian but the deterministic **two-point** `S`: the
all-ways total in the top-loop sector is `A088368 (~e·n!)` rather than `(2k−1)!!`, and the
Möbius weight `μ(0̂,σ)=∏_B (−1)^{|B|−1}(|B|−1)!` is precisely the **signed cyclic-interleaving
(ribbon) freedom** `(b−1)!` at a hub visited `b` times. The walk `x_0→…→x_{2k}` traces a single
Euler path through `G_σ` (its `2k` edges are exactly the `2k` steps), so each `σ` carries a
canonical ribbon/rotation structure and hence a genus; summing the signed interleavings is a
`'t Hooft`-style collapse in which the higher-genus contributions should cancel and the planar
survivors give the Catalan number. The free-cumulant identification from the fourth session is
the same statement read analytically: the two-point law `½(δ_a+δ_{−a})` has free cumulants
`κ_{2k}=(−1)^{k-1}C_{k-1}A^k` (a three-line `R`-transform computation), and free cumulants *are*
the genus-0 / non-crossing data.

**A negative result that pins down which "planar" is meant.** The *tempting* reading — "genus-0 =
the index partition `σ` is non-crossing (laminar) on `0<1<…<2k`" — is **FALSE**
(`04-computation/paley_starstar_noncrossing_monad.py`): the non-crossing-partition restriction of
`Σμ` gives `−1, 2, −6, 25, −132` (k=1..5), **not** Catalan, and the crossing remainder
`0, 0, +1, −11, +90` does not vanish. So the relevant genus is **not** the laminarity of the
coincidence partition; it must be the genus of the **walk-induced ribbon map on `G_σ`** (the
rotation system the Euler path induces at each hub), which is a finer invariant than how the
blocks nest on the line. Getting that construction right is the crux of the proof and is exactly
where the `(b−1)!` cyclic freedom lives.

So the three "`e`"s of the project meet here on one cancellation: `A088368~e·n!` is the
all-genus over-count, `R(p)→e=exp(−χ(−1))` is its image after the genus-0 projection, and the
tournament condition `p≡3 mod 4` (`χ(−1)=−1`) is what makes the surviving constant `e` rather
than `e^{−1}`.

---

### Honest status

- `M²=J−pI` / `S²=J−nI` and its consequences (even-series support, `g≡+1`, order `n^{k+1}`,
  `c_0=(−1)^kΣμ`): **PROVED** (elementary), verified `p≤43`. DRT-universal at leading order.
- Even-series count `=` A215257: **VERIFIED k≤5**; the bijection to indecomposable deque-sortable
  permutations is asserted from the matching counts + the connectivity↔indecomposability match —
  an explicit bijection is not yet written (handoff).
- The cycle-rank triangle and its identified rows (`1`, `3C(k,2)`, `A088368`): **VERIFIED k≤5**.
- `(★★) = (−1)^k C_k`: **VERIFIED k≤5**, value confirmed analytically as the two-point free
  cumulant. The **genus-0 localization is a strategy, not yet a proof**, and the *naive* version
  is **refuted**: localization onto non-crossing index partitions FAILS (`−1,2,−6,25,−132`). The
  live form needs (i) the ribbon/rotation structure the Euler walk induces on `G_σ` made explicit
  and (ii) a cancellation showing `Σ_{ribbon-genus>0} μ = 0`. This is the sharp form of handoff #1.

### Forward
1. **Prove `(★★)`** via genus localization on the **walk-induced ribbon map** (NOT the index
   partition — that is refuted above). Define the rotation system the Euler walk induces at each
   hub, compute the ribbon genus of every enumerated `σ` (k≤5 already on disk), and test whether
   `Σ_{ribbon-genus 0} μ = (−1)^k C_k`, `Σ_{genus>0} μ = 0`. If it holds, prove the cancellation;
   if it too fails, the cycle-rank triangle's row recursion (`t(k,2)=3C(k,2)`, diagonal A088368)
   is the fallback handle. This is the matrix-model proof the free-cumulant reading promises.
2. **Write the A215257 bijection** even-series patterns ↔ indecomposable deque-sortable
   permutations; it may itself carry the sign structure that proves `(★★)`.
3. **HYP-2308 remainder for non-circulant DRT**: with the leading order now `S²=J−nI`-universal,
   only the `o(n^{k+1})` bound is left — an expander-mixing estimate using `|λ|=√n` exactly.
   Test on a verified non-circulant DRT (n=15, skew-Hadamard order 16; check MISTAKE-011b/017).
