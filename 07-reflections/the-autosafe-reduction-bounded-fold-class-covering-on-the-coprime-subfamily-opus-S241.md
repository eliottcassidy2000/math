---
source: opus-2026-07-11-S241
status: A PROVED reduction lemma (owner-directed) + honest scope. The AUTO-SAFE lemma turns clearing at a
  composite modulus into a BOUNDED fold-class covering on the small coprime-to-q sub-family — the structured
  (mult-of-a-factor) speeds drop out. Verified 0 violations / 372830, 0-mismatch reduction. The odd-composite
  window handles 100% of random divisor-complete families, but adversarially no bounded window is a shortcut
  (a family blocks any fixed modulus by carrying a multiple of it) — reconfirming S238.
tags:
  - lrc14
  - divisor-complete
  - auto-safe
  - coprime-subfamily
  - fold-class-covering
  - proved-lemma
---

# The auto-safe reduction: bounded fold-class covering on the coprime sub-family

**opus-2026-07-11-S241.** Owner: attack the residual 8.5% (divisor-complete families) as a *bounded
fold-class covering on the small coprime sub-family*. This framing yields a clean **proved** reduction.

## The auto-safe lemma (PROVED, elementary; verified 0 / 372830)

> **Lemma.** Let `q` be composite with all prime factors ≤ 13, `q ∈ [15,28]` (danger band `{0,±1}`), and `v`
> a family with **no multiple of `q`**. Then for every **unit** multiplier `p` (`gcd(p,q)=1`) and every speed
> `v_i` with `gcd(v_i,q) > 1`:  `v_i·p mod q ∉ {0, 1, q−1}` (**auto-safe**).
>
> *Proof.* `gcd(v_i p, q) = gcd(v_i, q) = g > 1` (`p` a unit), so `v_i p mod q` shares the factor `g` with `q`;
> since `gcd(1,q) = gcd(q−1,q) = 1` it is `≠ ±1`, and it is `≠ 0` unless `q | v_i` (excluded). ∎
>
> **Consequence.** `bandCount(v,q,p) = #{coprime-to-q speeds with v_i p ≡ ±1}`, so
> **`v` clears at `q` (via a unit `p`) ⟺ the coprime-to-`q` sub-family misses some unit ±-fold-class mod `q`.**

The structured speeds — every multiple of a prime factor of `q` — **drop out**. Clearing becomes a **bounded**
covering of the `φ(q)/2` unit fold-classes by the **smaller** coprime sub-family. For `q = 15` there are only
**4** unit fold-classes `{±1,±2,±4,±7}`, and the coprime-to-15 sub-family is ~6 of 13 speeds. Verified
**0-mismatch** at `q ∈ {15,21,25,27}`. This is exactly the owner's "bounded fold-class covering on the small
coprime sub-family," and it is fully **provable** and **diameter-free**.

## Coverage: 100% of random families via odd composites

The odd-composite window `{15,21,25,27,33,35,39,45,49,55,63,65}` (prime factors ≤ 13, up to 65) clears
**100%** of 3000 random divisor-complete spread families — the entire residual, with **no primes** and **no
unbounded anti-concentration**. On random families, the auto-safe reduction closes it.

## Honest scope: no bounded window is a shortcut (reconfirms S238)

Adversarially the odd-composite window **fails**. A divisor-complete family can carry **multiples of the
window composites** — the multiple sits at residue 0 (the danger center) and blocks clearing there. Found:

`v = [3, 9, 35, 88, 98, 110, 189, 195, 225, 238, 264, 270, 273]`

is divisor-complete with a multiple of **15, 25, and 27**, so it clears at **none** of the 12 odd composites
≤ 65 — yet clears at `q = 16` (even) and primes `19, 23, 29, 31`. So for **any** fixed bounded window an
adversary blocks every modulus by a multiple, and the clearing modulus must **adapt** to the family (a modulus
it happens to lack a multiple of). This is the same obstruction as S238: no fixed bounded set of moduli is a
shortcut.

## Net

The owner's framing delivered a genuine **proved reduction** — the auto-safe lemma converts composite-clearing
into a bounded fold-class covering on the small coprime sub-family, an elementary, diameter-free tool that
handles 100% of random families and shrinks the problem (structured speeds drop out, few fold-classes remain).
It does **not** close the residual: the full disjunction — *the family clears at some modulus it lacks a
multiple of, whose coprime sub-family misses a fold-class* — is still the covering-system / anti-concentration
wall, and no bounded window shortcuts it (an adversary blocks any fixed window via multiples). What is new and
banked is the **structural reduction on the composite part**: clearing there is not a mystery but a bounded
covering by the coprime sub-family, with the divisor-structured speeds provably inert.

→ opus-S238 (no small-modulus shortcut — reconfirmed), opus-S232 (fold-class / summand-shell structure at
prime moduli — the composite analogue here), opus-S235 (band-edge margin — the strict `M > 1/14` this feeds),
THM-366 (covering ⟹ divisor-complete), opus-S239 (the shared crux). Files:
`lrc14_autosafe_coprime_covering_opus_S241.py` (+`.out`).
