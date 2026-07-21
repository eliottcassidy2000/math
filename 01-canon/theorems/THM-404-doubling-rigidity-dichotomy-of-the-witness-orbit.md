---
id: THM-404
title: Doubling-rigidity dichotomy — the AP witness orbit is ⟨×2⟩-connected iff n is odd
status: PROVED
source: opus-2026-06-03-S594
depends_on:
  - THM-403   # the witness orbit is (ℤ/n)*
related:
  - HYP-2126  # the dynamical-rigidity dichotomy (S585)
  - THM-398   # C′ (the even-n residual)
---

# THM-404 — doubling rigidity: the witness orbit is ⟨×2⟩-connected iff n is odd

## Statement

Let `W = { j/n : j ∈ (ℤ/n)^* }` be the AP's floor-witness orbit (THM-403). The
**doubling map** `t ↦ 2t` acts on the clock as `j ↦ 2j (mod n)`. Then:

> `2 · W ⊆ W` **iff `2 ∈ (ℤ/n)^*` iff `n` is ODD.**
>
> - **`n` odd:** `×2` is a bijection of `(ℤ/n)^*`, permuting `W`; the orbits are the
>   cosets of the cyclic subgroup `⟨2⟩ ≤ (ℤ/n)^*`, each of size `ord_n(2)`. The witness
>   orbit is **dynamically connected** (foliated by doubling cycles).
> - **`n` even:** `×2` maps every unit `j` to `2j`, which is even, hence
>   `gcd(2j, n) ≥ 2`, so `2j ∉ (ℤ/n)^*`. Thus `2·W ∩ W = ∅` — the doubling action
>   **fragments** the witness orbit (no nontrivial doubling cycle lies inside `W`).

## Proof

`2 ∈ (ℤ/n)^* ⟺ gcd(2,n)=1 ⟺ n` odd.

*Odd `n`:* multiplication by the unit `2` is a bijection of `(ℤ/n)^*`; it permutes the
witnesses `j/n`. Its orbits are the right cosets of `⟨2⟩`, each of cardinality the
multiplicative order `ord_n(2)`.

*Even `n`:* for any unit `j` (so `j` odd, `gcd(j,n)=1`), `2j` is even; since `n` is even,
`gcd(2j,n) ≥ 2`, so `2j` is a non-unit, i.e. `2j/n ∉ W`. Hence `×2` carries `W` entirely
out of `W`. ∎

## LRC corollary

The polynomial / sieve method propagates loneliness along the doubling `t ↦ 2t`
(Frobenius at 2). It requires the witness orbit `W` to be **doubling-connected**:

> **Odd `n`:** `W` is `⟨×2⟩`-connected, so the propagation rides the doubling cycles —
> consistent with the odd cases being tractable.
> **Even `n`:** `W` fragments (`2·W ∩ W = ∅`); the doubling propagation **stalls** — the
> dynamical face of the even-`n` residual (C′ / the `2q` apex). For `n=14`,
> `W = {1,3,5,9,11,13}/14` and `2·W = {2,6,10,18≡4,22≡8,26≡12}/14`, all non-units —
> total fragmentation, the prime-2 nonunit obstruction.

So the two rigidities of the same orbit (THM-403 static, THM-404 dynamical) **diverge
exactly at the even-`n` frontier**: static (cyclotomic, `(ℤ/n)^*`) holds for all `n`;
dynamical (`⟨×2⟩`-connectivity) holds iff `n` is odd.

## Verification

`lrc_doubling_rigidity_dichotomy_s585.py`: the unit witness orbit is `⟨×2⟩`-connected for
`n=5,7,9,11,13,15,17` (odd) and fragments into singletons for `n=6,8,10,12,14,16,18`
(even) — exactly the theorem.

**Artifacts:** see S585 (`lrc_doubling_rigidity_dichotomy_s585.out`). Builds on THM-403,
THM-398. (Formalises HYP-2126.)

## 2026-07-21 correction and unit-dilation addendum

The phrase "Frobenius-at-2 ramification" in the original LRC corollary was
field-theoretically inaccurate. In fact

```text
Q(zeta_14)=Q(zeta_7),
```

whose conductor is `7`, so the rational prime `2` is unramified. What fails
is the exact-order-14 group scheme/group algebra in characteristic `2`:
`mu_2` is non-etale, `-1=1`, primitive order-14 points collapse to order `7`,
and `F_2[C_14]` is nonreduced. THM-2041 gives the precise boundary: the
semisimple projector theorem works exactly when the characteristic is
coprime to the period.

The doubling failure is also specific to the multiplier `2`, not to every
Frobenius prime. For any integer `a`,

```text
aW=W  iff  gcd(a,n)=1,
```

and for a unit `a` its orbits are the cosets of `<a>` in `(Z/nZ)^*`. At
`n=14`, multiplication by `3` is one six-cycle:

```text
1 -> 3 -> 9 -> 13 -> 11 -> 5 -> 1.                       (1)
```

Thus characteristic `3` is good for the period-14 exact-period projector,
and THM-2041 preserves its whole primitive layer under cubing. This does not
yet repair the polynomial/sieve proof: no current theorem says its blocker or
strict-safe inequality is covariant under the cubing step. The precise new
target is a ternary propagation lemma carrying an endpoint-labelled safe or
dual certificate around cycle (1). If proved, it would remove the witness-
orbit fragmentation that is fatal to the doubling version.
