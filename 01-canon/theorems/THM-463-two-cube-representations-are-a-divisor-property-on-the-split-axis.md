---
id: THM-463
name: two-cube-representations-are-a-divisor-property-on-the-split-axis
status: PROVED (elementary, complete) + VERIFIED exact-integer (bijection all n <= 10^6; split lemma all primitive pairs <= 2000; full audit of all 5464 doubly-primitive taxicab numbers <= 10^12)
date: 2026-06-10
session: kind-pasteur-2026-06-10-S1 (Thread A lab)
depends_on:
  - THM-434   # Moser-ladder rosette count 12 + r_E(t) = 12 + 6B(t); the split-rich record-rung axis
resolves:
  - "THM-434 §Record-rungs honest flag, TAXICAB SIDE: the taxicab-Moser 1729 resonance IS a proven structural bridge (the tournament side was resolved as coincidence by HYP-2306 / the-1729-resonance-is-isolated reflection)"
  - "HYP-2367 (the taxicab-Moser bridge is structural): PROVED"
---

# THM-463: two-cube representations are a divisor property on the split axis

> **Namespace + provenance note.** Originally claimed and written as THM-461
> (kind-pasteur-2026-06-10-S1 reservation commit e465dea6); the reservation push was
> blocked by a network outage, and `monad-explorer-2026-06-10`'s distinct THM-461
> (unit-distance deletion ladder + floor reduction, commit c335323d) reached origin
> first. Renumbered THM-463 per first-come convention (cf. the THM-433/434 collision
> note in THM-434). All `_kpc1` scripts/outputs renumbered consistently. THM-462
> (cubic spectrum, same session) is unaffected.

> **One-sentence punchline.** 1729 sits at the bottom of both lists because both
> lists reward complete splitting in `ℤ[ζ₆]` — the taxicab–Moser resonance is
> structural; the tournament resonance (HYP-2306) was the coincidence.

## Statement

Fix an integer `n ≥ 1`.

### (A) The divisor criterion (representation ⟷ good divisor bijection)

Unordered pairs `{x, y} ⊂ ℤ` with `x³ + y³ = n` are in bijection with the
**good divisors** of `n`: positive divisors `d | n` such that, writing

```
   s := (d³ − n)/(3d)        (this must be an integer:  3d | d³ − n),
   Δ := d² − 4s = (4n − d³)/(3d),
```

`Δ` is a **perfect square** `e²`. Under the bijection `d = x + y`, `s = xy`,
`e = |x − y|`, and `{x, y}` is the **unique** root pair of `z² − dz + s`.

Sharpenings, all proved below:
- the classical parity side-condition `e ≡ d (mod 2)` is **automatic**
  (`e² = d² − 4s ≡ d² (mod 4)` forces matching parity), so the criterion is
  just "`s` integral and `Δ` a square";
- only divisors with `d³ ≤ 4n` can be good (`Δ ≥ 0`), so two-cube
  representability of `n` is decided by **at most `τ(n)` divisor checks** —
  the cubic-curve membership question collapses to divisor arithmetic;
- edge cases are uniform: `x = y ⟺ e = 0 ⟺ d³ = 4n`; `n = m³` gives the good
  divisor `d = m` with `{x,y} = {m, 0}`; a representation with a negative
  member is exactly a good divisor with `e > d` (and `d > 0` always).

### (B) The split lemma (the Eisenstein cofactor)

`x³ + y³ = (x + y)·q` with `q = x² − xy + y² = N(x + yω)`, `ω = ζ₃` a
primitive cube root of unity — **the Eisenstein norm form**, the same form
(up to the `GL₂(ℤ)` change `y ↦ −y`, which turns it into THM-434's
`x² + xy + y² = N(x + yζ₆)`) whose representation count `r_E` counts the
transverse unit vectors of THM-434's Moser-ladder rosettes.

If `gcd(x, y) = 1` (a **primitive** representation) then:
1. `q` is odd, and every prime `p ≠ 3` dividing `q` satisfies
   **`p ≡ 1 (mod 3)`** (split in `ℤ[ω]`);
2. `v₃(q) ≤ 1`, with `v₃(q) = 1 ⟺ 3 | x + y`.

Hence the cofactor `m = n/d = q` of every primitive representation is

```
   m = 3^ε · ∏ pᵢ^{aᵢ},   ε ∈ {0, 1},   all pᵢ ≡ 1 (mod 3),
```

so in THM-434's notation `B(m) = ∏(aᵢ + 1) ≥ 1` and `r_E(m) = 6B(m) ≥ 6 > 0`:
**the cofactor of every primitive two-cube representation lies exactly on
THM-434's split-rich axis** — it is itself a transverse-vector-bearing rung of
the Moser ladder. (Non-primitive scaling: `gcd = g` gives cofactor
`g²·(split-axis number)`.)

### (C) Taxicab corollary (where the inert content must hide)

For every primitive representation with divisor `d`: every **inert** prime
`p ≡ 2 (mod 3)` dividing `n` satisfies `v_p(d) = v_p(n)` (the *full* inert
content of `n` sits inside `d`), and `v₃(d) ≥ v₃(n) − 1`. Consequently if `n`
has `k ≥ 2` primitive representations with divisors `d₁, …, d_k` (a
**doubly-primitive taxicab number** when `k = 2`), then

```
   p^{v_p(n)} | gcd(d₁, …, d_k)      for every inert prime p | n.
```

A doubly-primitive taxicab number with `gcd(d₁, d₂) = 1` is therefore free of
inert primes; `1729` (with `d₁, d₂ = 13, 19`) is the first example.

## Why it matters (the bridge, and what it does NOT say)

THM-434 proved `#units(L_t) = 12 + r_E(t) = 12 + 6·Σ_{d|t} χ₋₃(d)` and flagged,
honestly: *"the shared content [of taxicab-1729 and Moser-1729] is only the
factorization 1729 = 7·13·19, not (yet) a proven structural bridge."* The
companion reflection (`the-1729-resonance-is-isolated…`) then **severed the
tournament side**: `r(p) = H(T_p)/|Aut|` has no splitting law (`p = 19, 23`
ratios carry inert primes), so tournament-1729 is a smoothness accident.

This theorem resolves the **taxicab side in the opposite direction**: being a
sum of two coprime cubes is *literally* a statement in the arithmetic of
`ℤ[ζ₆]` — the splitting (A) is the norm factorization
`n = (x+y)·N(x+yω)`, and the cofactor is constrained (B) to the exact
divisor-character axis that makes THM-434's rosettes dense. The two lists meet
at their bottom element for one proven reason — **complete splitting in
`ℚ(√−3)` is the cheapest way to be divisor-rich**:

- **Taxicab list:** `Ta(2) = 1729` (verified: no `n < 1729` has two positive
  representations). Both representations are primitive; both divisors `13, 19`
  are split primes and both cofactors `133 = 7·19`, `91 = 7·13` are split
  products. Every one of 1729's eight divisors is `≡ 1 (mod 3)`:
  `B(1729) = τ(1729) = 8`, the maximum possible for a 3-prime-factor integer.
- **Moser record-rung list:** `1729 = 7·13·19` (the product of the **three
  smallest split primes**) is the first `t` with `B(t) = 8` — over all `t` the
  `B`-records fall at `t = 1, 7, 49, 91, 637, 1729` (`B = 1,2,3,4,6,8`; `B = 7`
  would need `p⁶ ≥ 7⁶ = 117649`), and restricted to THM-434's non-degenerate
  ladder rungs (excluding `4t−1 = 3·square`, which removes `7` and `91`) they
  fall at `t = 3, 13, 49, 133, 637, 1729` — exactly THM-434's published record
  rungs. So `1729` is the first **60-unit rosette** (`12 + 6·8`).

The same eight `χ₋₃`-positive divisors that make `B(1729) = 8` (48 transverse
unit vectors) are the divisor pool from which the good divisors `{13, 19}` are
drawn. One arithmetic, two records.

**What it does NOT say (honest precision, from the 10¹² audit):** the theorem
does *not* force a multi-representation number to be completely split. Inert
primes can and do divide doubly-primitive taxicab numbers — `1952` of the
`5464` below `10¹²`, the first being `4342914 = 2·3²·31·43·181` — but Corollary
(C) confines them: in **all 10931 primitive representations audited**,
`v_p(dᵢ) = v_p(n)` for every inert `p`, i.e. the inert content hides in
`gcd(d₁, d₂)` and never touches a cofactor. The two lists share their minimum
for a structural reason; they are not the same list (the next doubly-primitive
taxicab number `20683 = 13·37·43` has `B = 8` but is not a `B`-record).

## Proof (complete, elementary)

### (A) The bijection

**Forward.** Let `x³ + y³ = n ≥ 1` and set `d := x + y`, `s := xy`,
`q := x² − xy + y²`. Since `4q = (2x − y)² + 3y² ≥ 0` with equality only at
`x = y = 0` (impossible: `n ≥ 1`), we have `q ≥ 1`. From
`n = (x + y)(x² − xy + y²) = d·q` and `q ≥ 1, n ≥ 1`: `d = n/q > 0`, so `d` is
a positive divisor of `n` with integer cofactor `q`. Newton's identity
`x³ + y³ = (x + y)³ − 3xy(x + y)` gives `n = d³ − 3sd`, so
`s = (d³ − n)/(3d) ∈ ℤ`, i.e. `3d | d³ − n`. Finally
`Δ := d² − 4s = (x + y)² − 4xy = (x − y)² = e²` with `e = |x − y| ≥ 0`, and
`d² − 4(d³ − n)/(3d) = (3d³ − 4d³ + 4n)/(3d) = (4n − d³)/(3d)`, so
`Δ = (4n − d³)/(3d) ≥ 0` forces `d³ ≤ 4n`. So `d` is good.

**Backward.** Let `d` be a good divisor: `s := (d³ − n)/(3d) ∈ ℤ` and
`Δ := d² − 4s = e²`, `e ≥ 0`. *Parity is automatic:* `e² = d² − 4s ≡ d²
(mod 4)`; since a square is `0 (mod 4)` for an even root and `1 (mod 4)` for an
odd root, `e ≡ d (mod 2)`. Hence `x := (d + e)/2` and `y := (d − e)/2` are
integers with `x + y = d` and `xy = (d² − e²)/4 = (d² − Δ)/4 = s`, and

```
   x³ + y³ = (x + y)³ − 3xy(x + y) = d³ − 3sd = d³ − (d³ − n) = n.
```

**Mutually inverse / uniqueness.** A good `d` determines `s` (from `d` and `n`
alone) and hence the unordered root pair of `z² − dz + s` — exactly one
representation. Conversely, two representations `{x, y} ≠ {x′, y′}` with the
same `d` would share `s = (d³ − n)/(3d)`, making both the root pair of the same
quadratic, hence equal — contradiction. So representation `↦ d` is injective,
and the backward construction shows every good divisor is hit. Bijection. ∎

### (B) The split lemma

Let `gcd(x, y) = 1`, `q = x² − xy + y²`.

**Norm form.** With `ω = ζ₃` (so `ω + ω̄ = −1`, `ωω̄ = 1`):
`N(x + yω) = x² + xy(ω + ω̄) + y²ωω̄ = x² − xy + y² = q`. The substitution
`y ↦ −y` carries this to THM-434's `x² + xy + y² = N(x + yζ₆)`: the same set of
represented integers, the same `r_E`.

**(0) `q` is odd.** Mod 2 the admissible cases `(x, y) ≡ (1,0), (0,1), (1,1)`
all give `q ≡ 1`; `(0,0)` is excluded by primitivity.

**(1) Split primes.** Let `p | q` be prime, `p ≠ 3` (so `p` odd by (0)).
First `p ∤ y`: otherwise `p | q + xy − y² = x²`, so `p | x`, contradicting
`gcd(x, y) = 1`. From `4q = (2x − y)² + 3y²`: `(2x − y)² ≡ −3y² (mod p)`, so
`t := (2x − y)·y⁻¹ (mod p)` satisfies `t² ≡ −3 (mod p)`. Now set
`z := (t − 1)/2 (mod p)` (`p` odd). Then

```
   z² + z + 1 = ((t − 1)² + 2(t − 1) + 4)/4 = (t² + 3)/4 ≡ 0 (mod p),
```

so `z³ ≡ 1 (mod p)` and `z ≢ 1` (else `3 ≡ 0`, i.e. `p = 3`, excluded). Thus
`z` has multiplicative order exactly `3` in `F_p^×`, and by Lagrange
`3 | p − 1`, i.e. `p ≡ 1 (mod 3)`. (Fully elementary — no quadratic
reciprocity.) ∎

**(2) The factor 3.** `q = (x + y)² − 3xy`. If `3 | q` then `3 | (x + y)²`, so
`3 | x + y`; conversely `3 | x + y ⟹ 3 | q`. In that case write `x + y = 3k`:
`q/3 = 3k² − xy`, and `3 | q/3` would force `3 | xy`; but `3 | x` together with
`3 | x + y` gives `3 | y` (and symmetrically), contradicting primitivity. So
`v₃(q) = 1` exactly when `3 | x + y`, else `v₃(q) = 0`. ∎

### (C) The corollary

For a primitive representation with divisor `d` and cofactor `m`: by (B)(1) an
inert prime `p ≡ 2 (mod 3)` has `v_p(m) = 0`, so `v_p(d) = v_p(n)`; by (B)(2)
`v₃(m) ≤ 1`, so `v₃(d) ≥ v₃(n) − 1`. With `k ≥ 2` primitive representations,
`p^{v_p(n)}` divides every `dᵢ`, hence their gcd. ∎

## Verification (exact integers; two independent scripts)

`04-computation/taxicab_moser_bridge_kpc1.py` →
`05-knowledge/results/taxicab_moser_bridge_kpc1.out` (10.3 s):

- **(a) Bijection, all `n ∈ [1, 10⁶]`:** representation side enumerated over
  the box `|x|, |y| ≤ 1000` — **provably complete**: a representation of
  `n ≤ 10⁶` with max-coordinate `x` has `n = d(3x² − 3dx + d²)` non-decreasing
  in `d` (`∂/∂d = 3(x − d)² ≥ 0`), so `n ≥ 3x² − 3x + 1` and `x ≤ 577` (the
  run's observed max coordinate is exactly 577). Result: `12442` unordered
  integer representations vs `12442` good divisors, **0 mismatches at every
  `n`**; parity condition: **0 violations** (automatic, as proved); `11899`
  two-cube integers `≤ 10⁶`; max representation count 3 (first at `n = 728`).
- **(b) Split lemma, all `1216588` primitive pairs `1 ≤ y ≤ x ≤ 2000`:**
  0 violations (oddness, split-primes-only, `v₃` law).
- **(c) Doubly-primitive taxicab numbers `≤ 10¹²`** (meet-in-the-middle,
  `x ≤ 10⁴`, 44 164 746 pairs; numpy int64 used only for enumeration, max value
  `2·10¹² ≪ 2⁶³`, every collision re-derived in pure ints): `40426` integers
  with ≥ 2 positive representations, of which **`5464` doubly-primitive**;
  first 25 listed with factorizations (first five: 1729, 20683, 40033, 149389,
  195841 — all completely split). Inert primes occur in `1952` of `5464`
  (first: `4342914 = 2·3²·31·43·181`), always inside `gcd(d₁, d₂)`.
- **(d) The 1729 slice:** table below; `r_E(1729) = 48` by direct norm-form
  count, `B(1729) = 8` by divisor character, `12 + 48 = 60` units (THM-434);
  `B`-records `(1,1)(7,2)(49,3)(91,4)(637,6)(1729,8)`, no `t ≤ 1729` with
  `B = 7`; `Ta(2) = 1729` re-verified (no smaller `n` with two positive reps).

`04-computation/taxicab_inert_audit_kpc1.py` →
`05-knowledge/results/taxicab_inert_audit_kpc1.out` (7.3 s) — removes the one
circularity (the bridge script tested inert primes *through* Corollary C):
**directly factors all `5464` numbers and all cofactors**, auditing every one
of the `10931` primitive representations:

- AUDIT-1 (cofactor `= 3^{0|1}·`split primes, `ε = 1 ⟺ 3 | d`): **0 violations**;
- AUDIT-2 (`v_p(d) = v_p(n)` for every inert `p | n`): **0 violations**;
- AUDIT-3 (`v₃(n) − v₃(d) ≤ 1`): **0 violations**;
- AUDIT-4: direct inert count `1952` = gcd-derived count `1952`; inert prime
  histogram starts `2: 649`, `5: 561`, `11: 224`, `17: 160`, `23: 117`.

### The 1729 slice table (full divisor scan)

`1729 = 7·13·19`, all three `≡ 1 (mod 3)` — completely split in `ℚ(√−3)`.
Every divisor is `≡ 1 (mod 3)`, so the mod-3 gate `3d | d³ − n` is open at
**all eight** divisors (`k = n/d ≡ 1 ≡ d² (mod 3)`); only the square condition
selects:

| `d`  | `s = (d³−1729)/(3d)` | `Δ = (4n−d³)/(3d)` | square?        | rep `{x, y}` | cofactor `m = n/d` |
|------|------|---------|----------------|-----------|-----------|
| 1    | −576 | 2305    | no (48² = 2304) | —         | — |
| 7    | −66  | 313     | no             | —         | — |
| **13** | 12 | **121 = 11²** | **yes**   | `{1, 12}` | `133 = 7·19` (split·split) |
| **19** | 90 | **1 = 1²**    | **yes**   | `{9, 10}` | `91 = 7·13` (split·split) |
| 91   | 2754 | < 0     | no             | —         | — |
| 133  | 5892 | < 0     | no             | —         | — |
| 247  | 20334 | < 0    | no             | —         | — |
| 1729 | 996480 | < 0   | no             | —         | — |

Good divisors `= {13, 19}` exactly; both representations primitive
(`gcd(1,12) = gcd(9,10) = 1`); `gcd(13, 19) = 1` ⟹ no inert prime divides 1729
(Cor. C) — and indeed there is none. Moser side: `B(1729) = 8` ideals of norm
1729 ⟹ `r_E = 48` transverse unit vectors ⟹ the **60-unit rosette**, the first
on the ladder. And `9³ + 10³ = 12³ + 1` — the first taxicab number is the
Ramanujan near-miss sitting one above Klein's `1728 = j(i)` (THM-436 lane; that
adjacency remains a fact about the *number*, per HYP-2306).

## Honest scope

- (A), (B), (C) are **PROVED** in full, elementarily (no reciprocity, no
  computation needed); the scripts are independent exact checks over the stated
  ranges. The bijection counts unordered *integer* representations (negative
  and zero coordinates included); the positive-representation sublist is the
  classical taxicab list.
- The bridge claim is exactly this: the taxicab and Moser-rosette lists are
  governed by the *same* `ℤ[ζ₆]` divisor arithmetic, and their shared bottom
  element 1729 is structural (product of the three smallest split primes ⟹
  first `B = 8` rosette AND a maximally divisor-rich two-rep configuration).
  It is **not** claimed that taxicab numbers are always completely split
  (refuted in-range: 1952/5464 carry inert primes, confined to `gcd(d₁,d₂)`),
  nor that the two lists agree beyond 1729 (they do not: 20683 is not a
  `B`-record).
- `r_E`/`B` here are THM-434's objects; the form appears as `x² − xy + y²`
  rather than `x² + xy + y²` — `GL₂(ℤ)`-equivalent (`y ↦ −y`), same values,
  same counts.
- The contrast with the tournament lane is inherited from HYP-2306 (refuted
  there, p = 19/23 counterexamples), not re-litigated here.
