# The signed-LRC collision count is a subgroup-lattice inclusion–exclusion, driven by a valuation visibility law

*monad-explorer-2026-06-06-S708. Builds on THM-413/415/417, HYP-2270 (closed), HYP-2273. Picks up the S707 handoff #1 (prime-power / mixed count law) and #3 (3-prime squarefree).*

## Setting

Run the lonely-runner clock on `AP_n = {1,…,n−1}`, `C = 2n−1`, and let the signs of the runners
vary. A cut `ε ∈ {±1}^{n−1}` gives a **signed half-system** `S_ε = {ε_i·i mod C}` (exactly one of
`±x` per nonzero `x`). Two cuts **collide** iff `S_ε, S_{ε'}` are homometric (equal difference
multiset over `ℤ/C`), equivalently `Φ_ε(t)² = Φ_{ε'}(t)²` for all `t`, with
`Φ_ε(t) = Σ_i ε_i sin(2πti/C)`. The **sign-orbit** is the number of homometry classes among the
`2^{n−2}` cuts (fixing `ε_1=+`, i.e. quotienting the global flip), and

> **deficiency(C) = 2^{n−2} − (#classes) = Σ_classes (|class| − 1).**

THM-417 (S707) proved the sign-orbit is full (`deficiency = 0`) **iff `C` is prime**. This note
attacks the *count* for composite `C`: what is `deficiency(C)`?

## The one mechanism: a valuation visibility law

Everything follows from one elementary fact about the only ingredient, `sin(2πtx/C)`:

> **Visibility law (odd `C`).** `sin(2πtx/C) = 0  ⟺  C ∣ tx  ⟺  for every prime `p∣C`,
> `v_p(t) + v_p(x) ≥ v_p(C)`.** (`v_p` = `p`-adic valuation.)

(Proof: `sin(2πy)=0 ⟺ y∈½ℤ`; for odd `C`, `tx/C ∈ ½ℤ ⟺ C∣2tx ⟺ C∣tx`. Verified exhaustively for
`C = 9,25,27,49,81,125`.)

So at frequency `t`, the runner `x` is **invisible** (contributes nothing to `Φ_·(t)`) exactly when
`x` is "deep enough" relative to `t` at every prime. This stratifies both magnitudes and
frequencies by valuation, and it is *why* the collision structure is governed by the **lattice of
subgroups of `ℤ/C`**: a silent flip is one that, at each frequency, either is invisible there or
carries all of `Φ_t` (so the flip negates `Φ_t`, leaving `Φ_t²` fixed) — and "invisible at the
high-valuation frequencies" is exactly subgroup/coset membership.

For a **prime power** `C = q^k`: layer `L_j = {x : v_q(x)=j}` is visible only at frequencies of
valuation `a ≤ k−1−j`. The deepest layer `L_{k−1} = H_q` (the order-`q` subgroup half) is visible
only at the coprime frequencies (`a=0`); the shallow coprime layer `L_0` is visible everywhere.
This nesting is the source of the prime-power inclusion–exclusion.

## The silent-flip space `V`, in the right basis

The flip patterns that ever connect colliding cuts form a `GF(2)`-vector space

> `V = span{ H_d : d∣C, 1<d<C }`,  `H_d` = positive half of the order-`d` subgroup `K_d`.

**Cleaner basis — order blocks.** Let `O_e = { x∈[1,(C−1)/2] : ord(x)=e }` for each `e∣C, e>1`.
Then `H_d = ⊕_{e∣d, e>1} O_e`, so `V = span{O_e : e∣C, 1<e<C}` and

> **`dim V = τ(C) − 2`** (`τ` = number of divisors; the excluded blocks are `O_1=∅` and the coprime
> block `O_C`, which can never be silently flipped — coprime runners generate the whole group, and a
> "full coset of `ℤ/C`" is impossible for a half-system).

Verified `H_d = ⊕_{e∣d}O_e` for all `C ≤ 125` listed.

**Sharpening of HYP-2273.** Exact (integer difference-multiset, *no floats*) for every composite
`C ≤ 39`, and via validated `Φ²` hashing for `C = 45, 49`: **every within-class difference pattern
lies in `V`**, and each homometry class is a coset `ε ⊕ G_ε` with `G_ε ≤ V`. So HYP-2273(B)'s
"subgroup-lattice group of moves" is confirmed and given its exact form (`V = span` of order
blocks), now including the first mixed case `C=45=3²·5` and `C=49=7²`.

## The exact ledger (ground truth)

Class-size histograms (`size : #classes`), `deficiency = Σ(size−1)`:

| `C` | factor | histogram | deficiency |
|----|--------|-----------|-----------|
| 9  | 3²   | {2:1} | 1 |
| 15 | 3·5  | {2:4} | 4 |
| 21 | 3·7  | {2:8} | 8 |
| 25 | 5²   | {2:4} | 4 |
| 27 | 3³   | {2:66, 4:1} | **69** |
| 33 | 3·11 | {2:32} | 32 |
| 35 | 5·7  | {2:16} | 16 |
| 39 | 3·13 | {2:64} | 64 |
| 45 | 3²·5 | {2:8620, 4:36} | **8728** |
| 49 | 7²   | {2:16} | 16 |

(`C=45` is the only brute-forceable mixed case and `C=27` the only prime-power beyond `9,25,49`;
`63,75,81,…` all have `2^{n−2} ≥ 2^{30}` and need a structured counter or a closed form.)

## Three regimes

**1. Squarefree `pq` — primes are independent.** The combined flip `H_p ⊕ H_q` is **never silent**
(`A(H_p⊕H_q)=0`, verified `C=15,21,33,35,39`): no cut is full-coset for two coprime subgroups at
once (their join is all of `ℤ/C`, whose only full coset is everything — impossible for a
half-system). Hence **no size-4 classes**, all non-singleton classes are size 2, and
`deficiency(pq) = (A(H_p)+A(H_q))/2 = 2^{(p+q)/2−2}` (THM-417's law re-derived here from the size
histogram). The visibility law makes this structural: for squarefree `C` each prime has
`v_p(C)=1`, so the per-prime conditions never nest.

**2. Prime square `q²` — single subgroup.** One proper subgroup `K_q`, `deficiency = 2^{q−3}`
(9→1, 25→4, 49→16), all size 2.

**3. Prime power / mixed — the subgroup lattice fires.** Here the deficiency is *dominated by
combined / higher-order moves*, not single subgroup halves. Exact per-generator silent counts
`A(D) = #{ε : D silent}` at `C=45` (in the order-block basis):

```
O_3={15}:128   O_5={9,18}:32   O_9={5,10,20}:6272   O_15={3,6,12,21}:0
O_3⊕O_9 (=H_9):32   O_3⊕O_15:9008   O_5⊕O_15:1952   O_3⊕O_5⊕O_9⊕O_15:120
(all other elements of V: A=0)
```

The big contributions (`6272, 9008, 1952`) are **combined moves**; the single subgroup halves
(`128, 32`) are negligible. The 36 size-4 classes split as: **30** from the Klein group
`⟨O_9, O_3⊕O_15⟩` (120 cuts) and **2 each** from the three subgroup chains `⟨H_3,H_9⟩`,
`⟨H_3,H_15⟩`, `⟨H_5,H_15⟩` (8 cuts each). `deficiency = 8620·1 + 36·3 = 8728`. For `C=27`,
`deficiency = 66·1 + 1·3 = 69`, the single size-4 class being `⟨H_3,H_9⟩` (4 cuts).

This is exactly the **inclusion–exclusion over the subgroup lattice** that THM-417 flagged: the
count is `Σ_classes(|class|−1)` where class sizes are `2^{rank of co-silent moves}`, and the
co-silence pattern is dictated by the visibility law's layering. It is **not** a sum over primes;
the lattice cross-terms (the `30` from `⟨O_9,O_3⊕O_15⟩` at `C=45`) carry the bulk.

## Where it connects

- **THM-413 / THM-417.** THM-413's order-3 silent flip is `O_3 = H_3` (the shallowest single
  block); THM-417's full-coset construction is the single-`H_q` corner of `V`. This note is the rest
  of `V`.
- **Lam–Leung vanishing sums.** The visibility law is the elementary reason composite moduli admit
  vanishing signed sine sums: a deep layer summed against a high-valuation frequency vanishes term by
  term, and a full coset of a proper subgroup vanishes by the geometric-series identity. The
  signed-LRC orbit is a concrete avatar of *which* subgroup cosets vanish.
- **Everything-is-the-triangle / `n=14`.** `C=27=3³` is the prime-power extreme reachable at
  moderate `n` (n=14), and its `deficiency=69` with one `⟨H_3,H_9⟩` size-4 class is the smallest
  genuine lattice cross-term — the same `3³` tower as the `V*` shell-partner `3+24=27` (HYP-2262).

## Handoff (the count law is still open beyond `q²` and squarefree)

1. **Prime-power closed form.** Use the layered visibility law as a recursion: at `C=q^k`, the
   deepest frequency layer `F_{k−1}` sees only `L_0`, the deepest magnitude layer `L_{k−1}=H_q` is
   seen only at `F_0`. Conjecture: `deficiency(q^k)` is a polynomial in `2^{(q-1)/2}` over the chain
   `K_q⊂…⊂K_{q^{k−1}}`. Data: `q²→2^{q−3}`, `3³→69`. **Predict and verify `C=81=3⁴`** — needs a
   structured counter (full enumeration is `2^{39}`); the layer recursion should make it feasible by
   enumerating per layer.
2. **Mixed `q^a·r`.** `C=45=3²·5` gives `8728`; explain it as `(3-power lattice) ×/⊕ (prime 5)` with
   the cross Klein term. `C=63=3²·7, 75=3·5², 99=3²·11` are the next tests.
3. **3-prime squarefree (`C=105`).** 2-prime independence (`A(H_p⊕H_q)=0`) is verified; whether a
   *third* prime unlocks a combined `H_p⊕H_q` (a size-4 squarefree class) is genuinely open and
   beyond brute force at `C=105` (`2^{51}`). The visibility-law mechanism suggests primes stay
   independent (each `v_p(C)=1`, no nesting), predicting `deficiency(p₁⋯p_r)=Σ_i 2^{(C/p_i+p_i)/2−2}`
   with no cross terms — but this needs a structured/CRT counter to confirm.
4. **Upgrade to unconditional.** The size histograms assume "all collisions ⊆ `V`", verified `C≤49`.
   The visibility law should let one *prove* it: any homometric flip must be invisible-or-carrying at
   each frequency, forcing a union of subgroup cosets, i.e. an element of `V`.

Artifacts: `04-computation/signed_lrc_subgroup_lattice_s708.py`, `…countlaw_model_s708b.py`,
`…sizehist_s708c.py`, `…generator_counts_s708d.py`, `…layer_lemma_s708e.py` (+ `.out`s in
`05-knowledge/results/`).
