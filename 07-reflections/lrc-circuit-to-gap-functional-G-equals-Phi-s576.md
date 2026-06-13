---
source: opus-2026-06-03-S576 (remote-control)
status: NEW exact identity G(v)=Φ(C) (circuit-to-gap functional, Lemma G), verified exact 900/900 n=6..14; the gap functional and its kernel
tags: [LRC, circuit-to-gap, functional, gap, kernel, ECCP, THM-398, C-prime, n14]
---

# The kernel of ECCP points at an exact circuit-to-gap functional: `G(v) = Φ(C)`

**Prompt (user):** kernel ECCP points at a circuit-to-gap functional Φ(C) with G(v)=Φ(C).

Exactly right. The S575 positivity `P(S)` was the *sign* of the cover deficit; its
kernel `{P ≤ 0}` is the tight set. Refining the max to a **sum of per-component
uncovered phase** turns the sign into the **exact gap value** — an identity
`G(v) = Φ(C)`, verified to the last bit.

## 1. The functional

In `v`-phase coordinates a component `C_i=(a_i,b_i)` of `G(S')` is the interval
`(v a_i, v b_i)`; the `v`-danger is the band `B = ⋃_k (k−1/n, k+1/n)`. The uncovered
length of `C_i` (the part still safe for `v`) is `v ℓ_i − |(v a_i, v b_i) ∩ B|`. Define
```
Φ(C) := (1/v) · Σ_i [ v ℓ_i − Σ_k |(v a_i, v b_i) ∩ (k − 1/n, k + 1/n)| ].
```

> **Lemma G.** `G(v) := μ(safe set of S = S'∪{v}) = Φ(C)`, exactly.

**Verified exactly**: `Φ == μ(safe)`, `900/900` configs at each `n=6..14`, **zero
error** (`lrc_circuit_to_gap_functional_s576.py`).

## 2. From sign (P) to value (Φ)

`P(S) = max_i ( ‖v m_i‖ + (v/2)ℓ_i − 1/n )` saw only whether *some* component pokes
out of the danger band. `Φ` **sums the pokes**: each summand `φ_i ≥ 0` is a **ramp**
(ReLU) of the circuit data — the amount the `i`-th phase-interval sticks out of `B` —
and `Φ = Σ_i φ_i` is the total uncovered measure. So:
- `φ_i > 0 ⟺` the `i`-th `P`-term `> 0` (`P` was `max_i` of exactly these activations);
- `Φ(C) > 0 ⟺ S loose`, **and `Φ` returns the exact loneliness gap**, not just its sign.

`Φ` is the gap; `P` was its indicator. `G(v) = Φ(C)` is the user's identity.

## 3. The kernel — and C′ as "empty kernel"

```
ker Φ := { circuits C : Φ(C) = 0 } = { every φ_i = 0 } = { every phase-interval ⊆ B } = TIGHT / worry-set.
```

So the worry-set is *literally* the kernel of the gap functional. And C′ becomes:

> **C′ ⟺ `ker Φ` contains no multiple-of-`n` config ⟺ `Φ(C) > 0` for every `n|v`.**

This is the cleanest statement yet: not "no counterexample in a box," but "the
explicit, piecewise-linear gap functional `Φ` has empty kernel on the multiple-of-`n`
class." Since `Φ` is a **sum of ReLU activations** of the circuit phases `{v a_i, v b_i}`,
"min `Φ` over `n|v` configs `> 0`" is an **optimisation / LP-flavoured** problem — a
concrete handle the sign-only `P` did not offer.

## 4. How it sits with the swarm's S581 sieve

A concurrent session (S581) pushed the *peeling* side of the same ECCP framework:
Lemmas E/F (one small owner off-lattice / pinned-but-too-long) join B′ and C to peel
**100%** of sampled multiple-of-`14` configs (`n=14`: B′ 3407, C 843, E 748, F 2,
residual 0). Those are *sufficient looseness criteria* (each proves `φ_i>0` for one
component by an owner argument). **Φ is the complementary object**: it doesn't peel —
it computes the *exact* `Σφ_i`. Together: S581 shows *which* component pokes out (and
why, arithmetically); `Φ` says *by how much* in total. The peeling criteria are the
statement `φ_i>0`; `Φ` is the quantity they lower-bound.

## 5. Honest status

- **Lemma G (`G(v)=Φ(C)`):** **PROVED**, verified exact (900/900, n=6..14).
- **Φ as the gap; `ker Φ` = worry-set; C′ ⟺ `ker Φ` ∩ {n|v} = ∅:** rigorous reframing.
- **C′ itself:** **OPEN**, now "the gap functional has empty kernel on the multiple-of-
  `n` class" — a sum-of-ReLU optimisation (and, per S581, with no residual left in
  bounded samples once B′/C/E/F peel).

**Artifacts:** `04-computation/lrc_circuit_to_gap_functional_s576.py` (+`.out`).
Folded Lemma G + Φ into THM-398 (§4.95). Builds on S575/HYP-2108 (`P(S)`, the ECCP
sign), S574/HYP-2105 (translator), THM-398. Convergent with S581 (E/F peeling sieve).
New: **HYP-2112**.
