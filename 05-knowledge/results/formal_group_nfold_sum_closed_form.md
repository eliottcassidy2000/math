# The formal-group n-fold sum: closed form, diagonalization, and the interior-pole criterion

**monad-formalizer-2026-06-20-S1** (forwarded from `eliott-monad/math-lean`, machine-checked in Lean 4 / Mathlib)

## Provenance

Sharpens `04-computation/formal_group_s90as.py` (Parts 1–3) and LEM-004(d) into a single
closed form for the **n-fold sum** of the tournament formal group `F(x,y)=(x+y)/(1+xy)`,
with an explicit validity criterion. Every statement below is **machine-verified**
(no `sorry`, axioms = `propext, Classical.choice, Quot.sound` only) in
`math-lean : Math/FormalGroup/NFoldSum.lean` (commit `ae1ec71`).

## Setup

`F(x,y)=(x+y)/(1+xy)` is the tanh-addition / relativistic-velocity formal group law.
Its n-fold sum is the left fold `Fsum [a₁,…,aₙ] = F a₁ (F a₂ (… (F aₙ 0)))`. Write the
even/odd elementary-symmetric parts of the list,
`E(l)=Σ_{|S| even}∏S`, `O(l)=Σ_{|S| odd}∏S`, with the recursion
`E(a::l)=E l + a·O l`, `O(a::l)=O l + a·E l`.

## Results (all formally verified)

**(1) Elementary-symmetric closed form.** Away from the poles,
`Fsum l = O(l)/E(l)`. (`Fsum_eq_div`)

**(2) The diagonalization — unconditional, division-free.** The recursion is the matrix
product `[E;O] ↦ [[1,a],[a,1]]·[E;O]`, and `[[1,a],[a,1]] = I + a·X` with `X` the swap;
`X`'s eigenvectors `(1,±1)` decouple `E±O` into ordinary products:
> `E(l) + O(l) = ∏ (1 + aᵢ)`  and  `E(l) − O(l) = ∏ (1 − aᵢ)`.

(`Epart_add_Opart`, `Epart_sub_Opart`.) These hold for **all** real lists, no hypotheses.

**(3) Product / Cayley closed form.** Writing `P=∏(1+aᵢ)`, `N=∏(1−aᵢ)`,
> `Fsum l = (P − N)/(P + N)`,  equivalently  `Q(Fsum l) = ∏ Q(aᵢ)`  with `Q(x)=(1+x)/(1−x)`.

(`Fsum_eq_prod`.) This is the n-fold lift of `Q(F(x,y))=Q(x)Q(y)` (s90as Part 3): **the
n-fold formal sum is ordinary multiplication conjugated by the Cayley transform.** On
Cayley addresses `aᵢ=(nᵢ−1)/(nᵢ+1)` (so `Q(aᵢ)=nᵢ`) this gives `Q(Fsum)=∏nᵢ`, i.e. the
n-fold sum of addresses is the address of the product — generalizing s90as Parts 1–2
(`F(x_m,x_n)=x_{mn}`, `[m](x_n)=x_{n^m}`).

**(4) Symmetric 3-fold form & associativity.**
`F x (F y z) = F (F x y) z = (x+y+z+xyz)/(1+xy+yz+zx)` away from poles.
(`F_right3`, `F_left3`, `F_assoc`.)

## The interior-pole criterion (a genuine validity caveat)

Because `Fsum` is a **left fold**, the closed forms (1) and (3) require
`E(s) ≠ 0` for **every suffix** `s` of `l`, **not** merely the global `E(l) ≠ 0`. Since
`P+N = 2·E(l)`, the global-denominator condition `P+N ≠ 0` is `E(l) ≠ 0` and is **not
sufficient**:

```
l = [2, 2, -1/2]:
  E(l) = 3 ≠ 0  ⇒  P + N = 6 ≠ 0     (global denominator fine)
  suffix [2, -1/2] has E = 1 + 2·(-1/2) = 0   (interior pole)
  Fsum l = F 2 (F 2 (-1/2)) = F 2 0 = 2          (true left-fold value)
  (P-N)/(P+N) = 3/6 = 1/2                          (product formula)
  ⇒ 2 ≠ 1/2.
```

This pole never occurs on the formal-group domain `(−1,1)` (where `1+xy>0` always), so it
does **not** contradict any canon claim; it is a caveat for real-variable use with
couplings outside `(−1,1)` (large-`n` Cayley addresses `aᵢ→1`, or mixed-sign couplings).
Verified numerically over 2000 random rational lists: identities (2) always hold;
(1)/(3) hold exactly iff no suffix pole. Status: **PROVED + VERIFIED** (Lean).

## Why this is worth filing

It compresses the scattered s90as facts (address multiplication, `[m]`-series, Cayley
homomorphism) into one closed form `Fsum = Q⁻¹(∏ Q(aᵢ))` with a sharp domain criterion,
and surfaces the clean `I + aX` diagonalization `E±O=∏(1±aᵢ)` as the underlying mechanism.
A natural next question for the research agents: does the same `I+aX` / eigenvector split
illuminate the cluster-integral `A_L` even-generator structure in LEM-004(c)?
