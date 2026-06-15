---
id: THM-505
title: The OCF non-spectral defect — H = (spectral skeleton) + (integer-linear combination of Witt/census defects); at n=7, H = [1+2c3+2c5+4C(c3,2)-4W6] + 4c6 + 2c7
status: PROVED for n=7, n=8 AND n=9 (closed forms derived by substitution from THM-500 + THM-502; verified 60000/60000 n=7, 12000/12000 n=8, 45000/45000 n=9). non-spectral dimension = n−5 ONLY for n<=8 (n=8 chain at 60000: (c6,c7,c8) determines H, dim=3); at n=9 the dimension exceeds n−5. INTRINSIC dimension (basis-independent rank, OCF packing basis, monad-explorer-S3) = 3,5,7 at n=8,9,10 — the n=9 trace-basis "6" was an OVER-COUNT (c8 and Q44 enter H only via their sum D35). GROWTH LAW: dim_nonspec(H)(n) = #{partitions of s≤n into odd parts ≥3} − 3 = (1,2,3,5,7,9,12,…) for n=6,7,8,9,10,11,12 — a restricted-partition generating function Π_{k odd≥3}1/(1−x^k); equals n−5 only for n≤8. VERIFIED n≤11 (rank 3,5,7,9 at n=8..11, all carriers independent, H in span); dim ≤ #{λ}−3 PROVED, equality CONJECTURE (no N_λ spectrally pinned). CLOSED FORM (monad-S4): the cumulative restricted-partition count is the named sequence A000009 — `#{λ:odd≥3,Σλ≤n} = q(n) = `#partitions of n into distinct(=odd) parts, so **dim(packing-vector) = A000009(n)−3 ~ exp(π√(n/3))**. TWO-DIMENSIONS CORRECTION (monad-S4): A000009(n)−3 is the non-spectral dim of the PACKING-COUNT VECTOR (N_λ), NOT of H. Since `H = 1+Σ_{j≤⌊n/3⌋} 2^j α_j` factors through the LEVEL-SUMS α_j (H never sees the split of α_j into length-types), **dim_func(H)(n) ≤ ⌊n/3⌋ (LINEAR)**, = ⌊n/3⌋ for n≥7; PROVED ≤, verified =2 at n=8 (vs carrier-rank 3). The fugacity-2 evaluation compresses exp(√n)→n/3. LINEARITY: H is universal-integer-linear in the FULL carrier set but NOT a bounded-degree polynomial in the simple cycles alone past n=7. (HYP-2513, dimension claim CORRECTED → partition law for the packing vector; H itself is ⌊n/3⌋-dim)
source: monad-explorer-2026-06-15
depends_on:
  - THM-499   # H = 1 + 2(c3+c5) + 4*alpha_2 at n<=6; alpha_2 first non-spectral OCF ingredient
  - THM-500   # H = 1 + 2*alpha_1 + 4*alpha_2 at n=7; alpha_1 = c3+c5+c7, c7 non-spectral
  - THM-502   # closed-walk census ladder: W_k = (1/k)sum mu(d) trA^{k/d} = simple cycles + overlap configs
  - THM-118   # c_k = tr(A^k)/k for k <= 5 (spectral)
related:
  - HYP-2513  # the general-n non-spectral defect principle
  - OPEN-Q-093
---

# THM-505 — the OCF non-spectral defect

THM-499/500 located the two spectral boundaries (where the eigenvalue spectrum of `A`
loses `H` at n=6, then `alpha_1` at n=7) and THM-502 explained the *mechanism*: the
spectrum fixes the Witt sums `W_k = (1/k) sum_{d|k} mu(d) tr A^{k/d}` (manifestly a
Z-combination of traces) but not the **split** of `W_k` into simple cycles `c_k` and
overlap configurations. This theorem turns that mechanism into an **exact closed form
for `H`**: it separates the OCF Hamiltonian-path count into a *spectral skeleton* plus
an explicit integer-linear combination of the **non-spectral carriers**, and pins the
coefficients via the independence-polynomial fugacity `x = 2`.

## Statement (n=7, PROVED)

For every tournament `T` on 7 vertices,

> **H(T) = [ 1 + 2 c3 + 2 c5 + 4·C(c3,2) − 4 W6 ]  +  4 c6 + 2 c7.**
>
> - `c3 = tr(A^3)/3`, `c5 = tr(A^5)/5`, `W6 = (tr A^6 − tr A^3)/6` are **spectral**
>   (Z-combinations of trace power-sums = symmetric functions of the A-spectrum), and
>   `C(c3,2) = c3(c3−1)/2` is a polynomial in `c3`. The bracket is the **spectral
>   skeleton** `S(T)`, constant on every cospectral class.
> - `c6` (number of directed 6-cycles) and `c7` (number of directed Hamiltonian
>   7-cycles) are the **two non-spectral carriers**.

**Corollary (the non-spectrality of H, made exact).** Within any cospectral class at
n=7 the skeleton `S` is constant, so

> **ΔH = 4 Δc6 + 2 Δc7.**

The entire non-spectral content of `H` at n=7 is the 2-dimensional vector `(c6, c7)`
read against the fixed weights `(4, 2)`.

## Proof (n=7)

1. **OCF (THM-500).** `H = 1 + 2·alpha_1 + 4·alpha_2`, where `alpha_1 = c3+c5+c7` is
   the number of directed odd cycles and `alpha_2` is the number of vertex-disjoint
   *pairs* of odd cycles. (`alpha_k = 0` for `k>=3` at n=7: three disjoint odd cycles
   need `>= 9` vertices.)
2. **alpha_2 at n=7 is triangle-pairs only.** A disjoint odd pair needs total length
   `<= 7`; the only option is `(3,3)` (since `3+5 = 8 > 7`). Hence
   `alpha_2 = #disjoint triangle pairs = C(c3,2) − p33`, where `p33` is the number of
   *intersecting* (vertex-sharing) triangle pairs.
3. **Census defect (THM-502).** `tr A^6 = 6 c6 + 3 c3 + 6 p33` and `tr A^3 = 3 c3`,
   so `W6 = (tr A^6 − tr A^3)/6 = c6 + p33`, i.e. **`p33 = W6 − c6`**. (This is the
   k=6 Witt defect `delta_6 = W6 − c6 = p33`.)
4. **Substitute.** `alpha_2 = C(c3,2) − (W6 − c6) = C(c3,2) − W6 + c6`. Then
   `H = 1 + 2(c3+c5+c7) + 4(C(c3,2) − W6 + c6) = [1 + 2c3 + 2c5 + 4C(c3,2) − 4W6] + 4c6 + 2c7.` ∎

Each step is an established theorem or a definition; the closed form is their
composition. **Verification:** 60000/60000 random n=7 tournaments satisfy all of
`H=1+2α1+4α2`, `α2=C(c3,2)−p33`, `p33=W6−c6`, and the closed form; the skeleton is
constant (|skel|=1) on every one of the 168 sampled cospectral classes; the within-class
law `ΔH = 4Δc6+2Δc7` holds on all 47 split classes, and `c6` co-varies with `c7` in
46/47 (the one exception is a class whose H-non-spectrality is carried by `c7` alone).
(`04-computation/ocf_nonspectral_defect_monad.py`, `05-knowledge/results/ocf_nonspectral_defect_n7_monad.out`.)

## Extension (n=8, PROVED, same substitution)

At n=8, `alpha_3 = 0` (three disjoint odd cycles need `>= 9` vertices), so the OCF still
truncates: `H = 1 + 2·alpha_1 + 4·alpha_2`. The disjoint odd pairs are `(3,3)` and now
also `(3,5)` (`3+5 = 8`): `alpha_2 = D33 + D35`, `D33 = C(c3,2) − p33`,
`D35 = c3·c5 − TF`. Using the census defects `p33 = W6 − c6` and `TF = W8 − c8 − Q44`
(from `tr A^8 = 8 c8 + 4 c4 + 8 Q44 + 8 TF`, `W8 = (tr A^8 − tr A^4)/8 = c8+Q44+TF`):

> **H = [ 1 + 2c3 + 2c5 + 4·C(c3,2) + 4·c3·c5 − 4 W6 − 4 W8 ]  +  2 c7 + 4 c6 + 4 c8 + 4 Q44.**

Equivalently, the **minimal-defect form** (fold `4c8+4Q44 = 4W8 − 4TF`, drop the
spectral `4W8` into the skeleton) uses only three carriers `c6, c7, TF`:

> **H = [ 1 + 2c3 + 2c5 + 4·C(c3,2) + 4·c3·c5 − 4 W6 ]  +  4 c6 + 2 c7 − 4 TF.**

Both forms verified 12000/12000 random n=8 tournaments. The coefficients of `c6` and
`c7` are unchanged (4 and 2) from n=7: `c6` always enters at the disjoint-pair level
(`2^2 = 4`), `c7` at the single-cycle level (`2^1 = 2`).

## Extension (n=9, PROVED, same substitution) — the TRIPLE level switches on

At n=9 the OCF gains its **first triple term**: `alpha_3 != 0` because three disjoint
triangles need exactly `3+3+3 = 9` vertices. (`alpha_4 = 0`: four disjoint odd cycles
need `>= 12`.) So `H = 1 + 2 alpha_1 + 4 alpha_2 + 8 alpha_3`. The disjoint odd pairs are
still `(3,3)` and `(3,5)` (as at n=8; `(5,5)`,`(3,7)` need `>= 10`), and the new
`alpha_3 = T333` = the number of **vertex-disjoint triangle TRIPLES**. With the same
defect substitutions (`p33 = W6−c6`, `TF = W8−c8−Q44`):

> **H = [ 1 + 2c3 + 2c5 + 4·C(c3,2) + 4·c3·c5 − 4 W6 − 4 W8 ]  +  2 c7 + 2 c9 + 4 c6 + 4 c8 + 4 Q44 + 8·T333.**

Verified **45000/45000** random n=9 tournaments (`ocf_nonspectral_n9_monad.py form9`,
`05-knowledge/results/ocf_nonspectral_n9_form_monad.out`). The new carrier `T333` enters
with coefficient **`8 = 2^3`** — the fugacity-2 weight of the **triple** independent-set
level — confirming the `2^level` rule at the next rung. The odd cycle `c9` enters at
level 1 (weight `2`), `c8` (even) at level 2 (weight `4`), exactly as `c7`, `c6` did.

## The fugacity-polynomial form (the clean generalization)

The `x = 2` evaluation is a shadow of a **polynomial identity in the fugacity `x`**. The
full independence polynomial `I(Ω, x) = Σ_k α_k x^k` of the odd-cycle conflict graph
decomposes (n=9):

> **I(Ω, x) = SKEL(x)  +  (c7 + c9)·x  +  (c6 + c8 + Q44)·x²  +  T333·x³**,
>
> `SKEL(x) = 1 + (c3 + c5)·x + (C(c3,2) + c3·c5 − W6 − W8)·x²`  (spectral coefficients).

The non-spectral carriers are **sorted by their independent-set level = power of `x`**:
level 1 (`x¹`) = the new odd cycles `c7, c9`; level 2 (`x²`) = `c6, c8, Q44` (entering via
the disjoint-odd-**pair** overlap defects); level 3 (`x³`) = `T333` (the disjoint triangle
**triple**). `H = I(Ω, 2)` sets the weights `(2, 4, 8) = (2¹, 2², 2³)`, recovering the
closed form. Other special values: `I(Ω, 1) = SKEL(1) + (c7+c9) + (c6+c8+Q44) + T333` is
the **total number of odd-cycle packings** (independent sets in Ω), the same carriers at
weight 1; `I'(Ω, 0) = α_1 = c3+c5+c7+c9`. The `x = 2` of the OCF is the **first nontrivial
oblong fugacity** (`x = n(n−1)` at `n=2`, roots `{2,−1}`; cf.
`fugacity-axis-and-vanishing-theorem`), which is why the carrier weights are clean powers
of 2.

## The non-spectral dimension: `n−5` for n<=8, but it BREAKS ABOVE `n−5` at n=9

How many independent non-spectral degrees of freedom does `H` carry? Probe by grouping
labelled tournaments by trace-vector (`tr A^3,...,tr A^n` = the full characteristic
polynomial, since `tr A^1 = tr A^2 = 0`) and asking which cycle counts must be added to
determine `H` within a cospectral class.

- **n=7:** `(c6, c7)` determines `H` (the closed form). `dim = 2`.
- **n=8 (chain, 60000 samples, 718 split cospectral classes):** the nested chain
  `sig → +c6 → +c7 → +c8` reduces the H-split count `718 → 333 → 163 → 0`. So
  `(c6, c7, c8)` determines `H`, and the overlap config `Q44` is **spectrally dependent**
  on `(c6, c7, c8)` (0 free buckets). `dim = 3`. **`n − 5` holds.**
- **n=9 (chain, 130000 samples, 29354 cospectral classes, 14804 split):** the chain
  `sig → +(c6,c7,c8) → +c9 → +Q44 → +T333` gives split counts
  `14804 → 482 → 24 → 1 → 0`. **`(c6, c7, c8, c9)` does NOT determine `H`** (24 cospectral
  witness-buckets with *identical* `(c6,c7,c8,c9)` but **different `H`**), and **neither
  does `(c6,c7,c8,c9,Q44)`** (1 residual split). Only adding BOTH the overlap config `Q44`
  AND the triple packing `T333` closes it. Every witness satisfies the closed form
  `ΔH = 4·ΔQ44 + 8·ΔT333` exactly (e.g. `(c6,c7,c8,c9)=(80,85,62,23)` carries
  `(H,Q44,T333) ∈ {(611,513,2),(615,514,2)}`; `(127,150,113,40)` carries
  `(1001,795,1),(1005,794,2)`). Cross-checks: `(sig,c6..c9,Q44)` leaves `T333` free in 1
  bucket; `(sig,c6..c9,T333)` still splits `H` in 5 (so `Q44` is the *finer* of the two
  overlap carriers, but both are needed).

> **CORRECTED dimension picture.** `dim_nonspec(H)` = 0, 1, 2, 3 for n = 5,6,7,8 (`= n−5`),
> then **= 6 at n=9** (`> n−5 = 4`). The minimal carrier set is the FULL
> `{c6, c7, c8, c9, Q44, T333}` (dimension capped at 6 by the closed form; lower bound 6 by
> cospectral witnesses showing `c9`, `Q44`, `T333` each detach). The simple-cycle counts
> `(c6, ..., c_n)` do **NOT** determine `H` past n=8: BOTH the overlap config `Q44` and the
> triple packing `T333` are genuinely **independent** non-spectral carriers. The dimension
> JUMPS `3 → 6` at n=9 — three carriers detach at once: the new odd cycle `c9`, the
> previously-dependent `Q44` (pinned at n=8 for lack of room, free at n=9), and the brand-new
> triple `T333`. The break occurs **exactly** when the triple independent-set level `α_3`
> (and the `(3,5)`-pair structure inside `α_2`) gains full room. **So the earlier claim "the
> overlaps are spectral shadows of the simple cycles" is true only for n ≤ 8** (corrects
> HYP-2513 / reflection `the-zeta-function-and-the-ocf-read-complementary-halves`).

## INTRINSIC dimension and the GROWTH LAW (monad-explorer-2026-06-15-S3) — the n=9 "6" is a trace-basis over-count; the true law is a partition function

The chain above works in the **trace basis** `{c6,c7,c8,c9,Q44,T333}` and counts 6 at n=9.
But that basis is **over-complete**: `c8` and `Q44` enter `H` *only through their sum*. The
closed form has `+4c8 + 4Q44`, and both come from `4·D35` via `D35 = c3c5 − W8 + c8 + Q44`
(`W8` spectral). So `H` depends only on the single quantity `c8+Q44 = D35` (mod spectrum);
adding `c8` then `Q44` as two chain steps double-counts one degree of freedom.

Working instead in the **OCF packing basis** — expand `H = Σ_λ 2^{|λ|} N_λ` by the
length-multiset `λ` (parts odd ≥3) of each disjoint-cycle packing, `N_λ` = #packings with
that multiset — gives a basis-independent count via the rank of the within-class
carrier-delta matrix (`ocf_nonspectral_n10_monad.py`, rank over ℚ, robust to small cospectral
classes since deltas pool across classes):

| n | OCF non-spectral carriers | intrinsic dim (RANK) | H in carrier span? |
|---|---|---|---|
| 8 | `{c7, D33, D35}` | **3** | yes |
| 9 | `{c7, c9, D33, D35, T333}` | **5** | yes |
| 10 | `{c7, c9, D33, D35, D37, D55, T333}` | **7** | yes |
| 11 | `{c7, c9, c11, D33, D35, D37, D55, T333, T335}` | **9** | yes |

Verified by cospectral sampling: rank `3,5,7,9` at n=8,9,10,11; every carrier drop-one
independent; `H` always in the carrier span (so dim is capped by the OCF closed form); OCF
identity holds `6000/6000` (n=10), `5000/5000` (n=11, 704 split cospectral classes). At n=9 the explicit reconciliation: rank`{c6,c7,c8,c9,Q44,T333}`=6 but
rank`{c6,c7,(c8+Q44),c9,T333}`=5 and already contains `H`.

> **GROWTH LAW (the corrected dimension claim).** With `N_∅=1`, `N_{3}=c3=trA³/3`,
> `N_{5}=c5=trA⁵/5` the only spectral/trivial packings (THM-118), and every other `N_λ`
> non-spectral and (verified n≤10) independent:
> ```
>   dim_nonspec(H)(n) = #{ partitions of s into odd parts ≥3, 0 ≤ s ≤ n } − 3
>                     = ( Σ_{s≤n} [x^s] Π_{k odd ≥3} 1/(1−x^k) ) − 3
> ```
> Values `n=6..14`: **1, 2, 3, 5, 7, 9, 12, 15, 19** (the n=8,9,10,11 entries 3,5,7,9 VERIFIED). The increment
> `dim(n)−dim(n−1) = p_{odd≥3}(n)` (partitions of `n` into odd parts ≥3). The sequence
> equals `n−5` *only* for `n ≤ 8` (each integer has ≤1 odd-≥3 partition there); it first
> exceeds `n−5` at `n=9`, where `9={9}={3,3,3}` is the first integer with **two** such
> partitions. So the "break above `n−5`" is exactly the point where the restricted partition
> function first exceeds 1. Asymptotically `dim ~ exp(c√n)` (Hardy–Ramanujan), sub-exponential
> but super-polynomial. The trace-basis chain count (`6` at n=9) over-counts because `c_even`
> and overlap configs `Q44,TF` are *change-of-basis coefficients* (inclusion–exclusion
> converting `N_λ` to spectral + lower carriers), not independent carriers.
>
> **Status:** decomposition `H=Σ_λ 2^{|λ|}N_λ` PROVED (regroup the OCF closed form);
> `dim ≤ #{λ}−3` PROVED (upper bound — `H` can't depend on more than its own carriers);
> EQUALITY (no `N_λ` spectrally pinned) VERIFIED n≤11, CONJECTURE general. See reflection
> `the-non-spectral-dimension-of-H-is-a-partition-function`.

## TWO DIMENSIONS: the growth law counts the PACKING VECTOR, not H (monad-explorer-S4)

The growth-law section measures the rank of the **individual packing carriers**
`{c7, D33, D35, …}` and gets `q(n)−3`. But `H` does **not** depend on the individual
`N_λ` — only on the **level-sums** `α_j = Σ_{|λ|=j} N_λ`. The OCF is literally
`H = I(Ω,2) = 1 + Σ_{j≥1} 2^j α_j`, and `α_j = 0` for `j > ⌊n/3⌋` (j disjoint odd cycles
need `≥3j` vertices). So `H` is a function of **at most `⌊n/3⌋`** quantities, and within a
cospectral class `ΔH = Σ_j 2^j Δα_j` is an *identity*. Hence two genuinely different
non-spectral dimensions:

| dimension | object | value | growth |
|---|---|---|---|
| `dim(fine)` | packing-count vector `(N_λ)` | `A000009(n) − 3` | `~exp(π√(n/3))` |
| `dim_func(H)` | what `H` actually depends on (level-sums `α_j`) | `⌊n/3⌋` (n≥7) | **linear** |

- **`dim_func(H) ≤ ⌊n/3⌋` is PROVED** (the OCF + the 3j-vertex floor; no computation
  needed). It is `< A000009(n)−3` for all `n ≥ 8`. So the S3 label "dim_nonspec(H) =
  q(n)−3" over-attributes the packing vector's complexity to `H`.
- **Closed form for `dim(fine)`** (resolves the S3 OEIS frontier): the cumulative
  restricted partition is A000009 by a one-line GF identity —
  `Σ_{s≤n}[x^s]Π_{k odd≥3}1/(1−x^k) = [x^n] Π_{k odd≥1}1/(1−x^k) = q(n)` (the `1/(1−x)`
  cumulative factor is the missing odd part `k=1`), and `q(n)` = partitions of `n` into
  odd `=` distinct parts `=` A000009`(n)`. Asymptotic `~ exp(π√(n/3))/(4·3^{1/4}n^{3/4})`.
  Bijective: carrier `λ` (`Σλ=s≤n`) ↔ odd-part partition `λ∪{1^{n−s}}` of `n`, the 1's =
  uncovered vertices; the `−3` removes `{1^n},{3,1^{n−3}},{5,1^{n−5}}`.
- **Verified** (`04-computation/ocf_two_dimensions_monad.py`, n=8, 159 286 members /
  1522 cospectral classes): carrier basis `{c7,D33,D35}` rank **3** (reproduces S3);
  level-sum basis `{c7, D33+D35}` rank **2**, `H` in span; `{D33,D35}` independent (rank
  2) yet `H` reads only their sum. So `dim(fine)=3` but `dim_func(H)=2`.

This is the same over-counting S3 caught at `6→5` (trace→packing basis), pushed one level
deeper: even the packing basis over-counts *for H*, because `H` sees only the level
marginals. The terminal honest answer for `H` is the level grading `⌊n/3⌋`; `A000009−3` is
the honest answer for the finer packing vector. The full fugacity polynomial `I(Ω,x)` (all
`x`) also carries only `⌊n/3⌋` non-spectral DOF — its coefficients *are* the `α_j`. To see
the `A000009−3` dimensions one needs the length-refined multivariate packing GF, not the
scalar `I(Ω,x)`. See reflection `H-reads-only-the-level-grading`.

## The linearity dichotomy (resolves the "is H linear in the carriers" question)

There are **two different** notions of "carrier," and they diverge past n=7:

1. **Universal-LINEAR carriers.** `H = (spectral skeleton) + Σ w_c · c` with *universal
   integer* weights `w_c = ±2^{level}`. By the OCF substitution this *always* exists, but
   the carrier set **must include the overlap/triple configs** `Q44, TF, T333, ...` — they
   are irreducible. Verified universal weights (n=8, exact within every cospectral class):
   `(c6,c7,c8,Q44) → (4,2,4,4)` and minimal `(c6,c7,TF) → (4,2,−4)`.
2. **FUNCTIONAL carriers.** The minimal cycle-count set `S` with `H = f(spectrum, S)` for
   *some* function `f` (not necessarily polynomial). At n≤8 `S = {simple cycles}` of size
   `n−5`; at n=9 it must include `Q44`.

These coincide at n≤7. At n=8 they have the same *size* (3) but different *sets*:
`{c6,c7,c8}` (functional) vs `{c6,c7,TF}` (linear). The reason: **`H` is a FUNCTION of
the simple cycles `(c6,c7,c8)` but NOT a bounded-degree polynomial in them.** Within
cospectral classes, the within-class fit `ΔH = lin(Δc6,Δc7,Δc8)` is **not exact** (best
linear residual ≈ 4), and neither is the degree-2 fit; only `120/645` classes admit even a
*per-class* linear `Q44~(c6,c7,c8)` fit. The map `(c6,c7,c8) ↦ Q44` is a genuine
non-polynomial function (`ocf_nonspectral_n8_linearity_monad.out`). At n=9 the simple
cycles fail even *functionally*, so `Q44` enters as its own universal-linear carrier.

> **Three-stage degradation of "what the simple cycles tell you about `H`":**
> n≤7 — simple cycles determine `H` *linearly*; n=8 — simple cycles determine `H`
> *functionally but non-polynomially*; n=9 — simple cycles do **not** determine `H` (the
> cycle *correlations* `Q44, T333` carry independent non-spectral information). Each step
> up in `n` is a step deeper into the correlation structure the spectrum cannot see.

## The fugacity-2 coefficient rule and the zeta/Euler-product picture

`H = I(Ω, 2) = sum_k 2^k · alpha_k`, where `alpha_k` = number of vertex-disjoint
odd-cycle k-collections (independent k-sets in the conflict graph `Ω`). Expanding each
`alpha_k` by inclusion–exclusion turns products of cycle counts into *overlap
corrections*, and those corrections are exactly the **Witt/census defects**
`delta_j = W_j − (spectral cycle counts)` — the non-simple primitive closed orbits.
The fugacity `x = 2` supplies the weights: **a defect entering at independent-set
level `j` carries coefficient `2^j`.** Hence the universal shape

> `H = (spectral skeleton in the W_k) + Σ (2^{level}) · (non-spectral orbit-split count)`.

This is the quantitative form of the "spectrum is mean-field, OCF is correlation"
principle (reflection `the-spectral-resolution-ladder-of-the-ocf`). The eigenvalue
spectrum is the **Bowen–Lanford / Artin–Mazur zeta function**

> `ζ_T(u) = exp(Σ_k tr(A^k) u^k / k) = 1/det(I − uA) = Π_k (1 − u^k)^{−W_k}`,

an Euler product over **primitive closed orbits**, with `W_k` = the count of
primitive (aperiodic) closed k-walks up to rotation (clean proof:
`tr A^k = Σ_{d|k} d·P_d` by orbit-size counting ⟹ `P_k = W_k` by Möbius). The zeta
function records only the *total* primitive-orbit count `W_k`; it is blind to the
simple-vs-overlap split inside each `W_k`. `H` reads precisely that split. The two are
**complementary readings of the same closed-orbit data** — which is why the
determinant/spectral lens (THM-468/472) is orthogonal to `H` (THM-499): they sit on
opposite sides of the simple/overlap partition of `W_k`.

## General-n principle (CONJECTURE — HYP-2513)

For every `n`, `H(T) = (spectral skeleton) + Σ m_k·(non-spectral carrier_k)`, an
**explicit universal-integer-linear functional** of the carriers, where the carriers are
the simple-cycle counts `c_k` (`k>=6`) AND the overlap/independent-set configs
(`Q44`, `TF`, `T333`, triple-overlaps, ...) that first acquire room, each weighted by
`±2^{(independent-set level at which it enters Ω)}`. The full independence polynomial
`I(Ω, x) = SKEL(x) + Σ (carrier)·x^{level}` is the fugacity generalization, with `x = 2`
the oblong-fugacity evaluation. PROVED (by substitution from the OCF + THM-502 census
defects) and verified for n=7, n=8, n=9.

**CORRECTED open content** (the n=9 break sharpens the conjecture):
- (i) The exact carrier list and `±2^{level}` weights at general `n` (driven by which
  disjoint odd-collections fit in `n` vertices and which overlap configs appear). Weights
  observed: pure `±2^{level}` through n=9.
- (ii) **RESOLVED → the GROWTH LAW (monad-explorer-S3, see section above).** The intrinsic
  (basis-independent) non-spectral dimension, in the OCF packing basis `N_λ`, is
  `dim_nonspec(H)(n) = #{partitions of s≤n into odd parts ≥3} − 3` = `1,2,3,5,7,9,12,…` for
  `n=6,7,8,9,10,…`. Equals `n−5` only for `n≤8`. The n=9 trace-basis count `6` was an
  over-count (`c8,Q44` enter `H` only via their sum `D35`); intrinsic dim is `5`. VERIFIED
  n≤10 (rank `3,5,7`, all carriers independent, `H` in span); the only general-`n` gap is
  proving no `N_λ` is spectrally pinned (then the law is a theorem).
- (iii) `H` is universal-linear in the *full* carrier set but is **not** a bounded-degree
  polynomial in the simple cycles alone past n=7 (the linearity dichotomy above).
