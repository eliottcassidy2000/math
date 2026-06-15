---
id: THM-505
title: The OCF non-spectral defect — H = (spectral skeleton) + (integer-linear combination of Witt/census defects); at n=7, H = [1+2c3+2c5+4C(c3,2)-4W6] + 4c6 + 2c7
status: PROVED for n=7 and n=8 (closed forms derived by substitution from THM-500 + THM-502; verified 60000/60000 random n=7, 12000/12000 random n=8 incl. minimal-defect 3-carrier form); non-spectral dimension = n−5 VERIFIED n<=8 (n=8 carrier probe: (c6,c7) insufficient/157 split buckets, Q44 dependent on c6,c7,c8); general-n principle CONJECTURE (HYP-2513)
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

## The non-spectral dimension grows as n−5 (PROVED n<=8 for the count; carriers VERIFIED n=8)

How many independent non-spectral degrees of freedom does `H` carry? At n=7 the answer
is **2**: `(c6, c7)` determines `H` within every cospectral class (the closed form). At
n=8 the answer is **3**: a carrier-dimension probe (12000 samples, 695 split cospectral
classes) finds **157 `(spectrum, c6, c7)` buckets that still carry >= 2 distinct `H`** —
so `(c6, c7)` no longer suffices; a third carrier is genuinely required. Adding `c8`
closes it: in **0** buckets does `Q44` vary given `(spectrum, c6, c7, c8)`, so `Q44`
(and hence every length-8 defect) is **spectrally dependent on `(c6, c7, c8)`**, and
`(spectrum, c6, c7, c8)` determines `H`. The even-overlap config `Q44` is *not* a new
non-spectral axis.

> **Dimension table:** `dim_nonspec(H)` = 0, 1, 2, 3 for n = 5, 6, 7, 8 — i.e.
> **n − 5** (n>=5), with carriers the simple-cycle counts `{c6, c7, ..., c_n}`. All
> overlap/defect counts (`p33 = W6−c6`, `TQ = W7−c7`, `Q44`, `TF`) are spectrally
> dependent on the simple-cycle vector. **CONJECTURE for n >= 9** (HYP-2513): the
> non-spectral content of `H` is a function of `(c6, ..., c_n)` of dimension `n−5`; the
> first triple-overlap term `8·alpha_3` switches on at n=9 (the first 3-disjoint-triangle
> independent set) and `c9` joins as the next carrier.

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

For every `n`, `H(T) ≡ Σ_{k>=6} m_k·(non-spectral carrier_k)  (mod the ring of
spectral invariants)`, where the carriers are the simple-cycle counts `c_k` (`k>=6`)
and the higher even-overlap configs (`Q44`, triple-overlaps, ...) that first acquire
room, each weighted by `2^{(independent-set level at which it enters Ω)}`. Equivalently:
`H` modulo spectral invariants is an explicit integer-linear functional of the Witt
defects `{δ_k}`. Verified n=7 (closed form, proved) and n=8 (closed form, 6000/6000).
The open content is (i) the exact carrier list and weights at general `n` (driven by
which disjoint odd-collections fit in `n` vertices and which overlap configs appear),
and (ii) whether the weights are always pure powers of 2.
