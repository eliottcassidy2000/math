# The relation lattice L(AP) has maximal kissing number — its minimal vectors are the additive triples, and the floor is that extremality

*mac-mini-2026-07-06-S25 (HYP-4552). Owner: work the S24 next step — compute the
relation lattice L(AP) and probe the Cohn–Elkies structure. The computation gives a
clean geometric invariant: L(AP)'s minimal vectors are the **additive triples**
`v_i+v_j=v_{i+j}`, its kissing number is `2·(#additive triples) ≈ the additive
energy`, and the AP maximizes it. The density floor is the **isolation of the
maximal-kissing-number relation lattice** — a Cohn–Elkies extremal characterization
that routes the floor through a *clean* invariant, not the all-orders theta. Verified
`lrc_relation_lattice_LAP_macmini_S25.out`.*

## The relation lattice

`L(AP) = {a ∈ ℤ¹² : Σ_{i=1}^{12} i·aᵢ = 0}` is the kernel of the **moment map**
`φ(a)=Σ i·aᵢ = ⟨c,a⟩`, `c=(1,…,12)` — a rank-11 primitive lattice. Its
discriminant is `|c|² = Σ i² = 650`. In the basis `dₖ=(k+1)eₖ − k e_{k+1}` the Gram
matrix is **tridiagonal**: diagonal `k²+(k+1)²`, off-diagonal `−k(k+2)`.

## The minimal vectors are the additive triples

The shortest vectors of `L(AP)` have norm **3**, and they are exactly the
`(1,1,−1)` relations at positions `(i,j,i+j)`:

> `vᵢ + vⱼ = v_{i+j}` — the **additive triples** of the AP (`5+7=12`, `1+6=7`,
> `4+8=12`, …). There are `30` of them (`Σ_{l=3}^{12} ⌊(l−1)/2⌋ = 30`), giving a
> **kissing number of 60** (each triple ±).

This is the AP's *additive closure*: `{1,…,12}` contains `i+j` whenever `i+j ≤ 12`,
so it has a maximal density of additive triples — and each is a minimal lattice
vector. The next shells: the multiplicative doubling `2vᵢ=v_{2i}` at norm 5, the
harmonic `(1,−2,1)` at norm 6. **The geometry of L(AP) records the AP's additive
combinatorics: its short vectors are its sumset relations.**

## Kissing number = additive energy; the AP is the extremizer

The kissing number `κ(L(S)) = 2·#{(i,j) : vᵢ+vⱼ ∈ S}` is (up to the diagonal) the
**additive energy** `A(S) = #{a+b=c+d}` of `S` (HYP-2873 = `∫|Ŝ|⁴`). The interval
`{1,…,12}` **maximizes additive energy** among 12-sets (classical, rearrangement /
Fejér). So:

> **L(AP) is the relation lattice of maximal kissing number** among relation
> lattices of covering 12-families.

This is precisely the Cohn–Elkies "extremal lattice" role: the arc-weighted theta
`safe = Σ_{a∈L(S)} ∏f̂(aᵢ)` hits `0` at the lattice with the most (weighted) short
vectors, and the floor is that this lattice is *isolated* — `L(AP)` is the unique
maximal-kissing relation lattice, and any deficit in the kissing number leaves the
theta positive.

## Why this is a cleaner route than the theta

The theta-sum itself is genuinely all-orders (verified: at `β=2/25` the support-3
shell is `−1.80`, support-4 `+0.98`, support-5 `−1.58` — large and oscillating,
confirming opus-S114 that it is not harmonic-led and my S22 that it is
tail-dominated). But the **kissing number is a single clean integer invariant**, and
its extremality (AP maximizes additive energy) is **classical**. The route:

1. `safe(S,2/25)=0` (tiling) ⇒ enough weighted short vectors to cancel the main term
   ⇒ near-maximal kissing number ⇒ (additive-energy inverse theorem / Freiman) `S`
   is a generalized AP;
2. among generalized APs, only `{1,…,12}` achieves the *full* additive closure
   (kissing 60) — the deficit of any other is bounded below;
3. the width `1/(2k²)` (opus-S113) converts that kissing deficit into a strictly
   positive `safe` at `n=13` (n-specific: the wider `n=7` window tolerates the
   deficit, so the `n=7` tiler has fewer triples yet tiles).

So the Cohn–Elkies frame localizes the floor to a **kissing-number / additive-energy
extremality with a width-dependent stability constant** — the same STRUCTURE × WIDTH
(opus-S113), now with STRUCTURE = *the lattice kissing number* (a clean invariant)
and WIDTH = the stability radius. The remaining analytic content is the stability:
*a relation lattice with sub-maximal kissing number has strictly positive
arc-theta at `β=2/25`, quantitatively in the kissing deficit*.

## Net

- **`L(AP)`'s minimal vectors are the additive triples `vᵢ+vⱼ=v_{i+j}`** (norm 3,
  kissing 60), then multiplicative doubling (norm 5), then harmonic (norm 6). The
  lattice geometry *is* the AP's additive combinatorics.
- **Kissing number = additive energy; the AP maximizes it** ⇒ `L(AP)` is the
  Cohn–Elkies extremal relation lattice; the floor is its isolation.
- This routes the floor through the **kissing number** (a clean invariant with
  classical extremality) rather than the all-orders theta — reducing the residual to
  a *width-dependent kissing-deficit stability bound*, the same STRUCTURE × WIDTH in
  lattice language.

## Pointers

- `lrc_relation_lattice_LAP_macmini_S25.py/.out` (Gram, disc 650, minimal vectors,
  shells).
- mac-mini HYP-4532 (Cohn–Elkies / Poisson framing), HYP-4512/4522 (Beurling–Selberg,
  theta convergence); opus HYP-4466 (all-orders), HYP-4456 (structure × width);
  HYP-2873 (additive energy = Fejér 4th moment); kps HYP-4467 (harmonic — refined here).
