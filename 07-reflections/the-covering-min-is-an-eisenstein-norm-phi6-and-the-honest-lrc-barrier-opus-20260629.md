# The covering-min denominator is an EISENSTEIN-INTEGER NORM: n²−n+1 = Φ₆(n) = N(n−ζ₆), n has order 6 (n³≡−1), the killer binds via the PROVEN identity (n−1)n²≡−n mod Φ₆(n) — the structure is all provable, but the global lower bound M(covering) ≥ n/Φ₆(n) IS LRC(n) itself (it exceeds 1/n), so that one wall is the conjecture

*opus-2026-06-29. Owner: attempt the global-minimality proof of the covering-min `n/(n²−n+1)`, have fun.
The fun is real — the covering min lives in the Eisenstein integers `Z[ζ₆]`, Heegner `−3`, hexagonal. The
honesty is also real — the lower bound is `>1/n`, so proving it proves LRC(n), and I can't (nor should
claim to) leap that wall. What IS provable, I prove.*

## The honest barrier (stated first)
`M_min(covering) = n/(n²−n+1) > 1/n`. So **`M(covering) ≥ n/(n²−n+1)` for all covering sets would IMPLY
LRC(n)** (which at `n=14` is open). The global lower bound is therefore the conjecture itself — not a
lemma below it. I do NOT claim it. Everything below is the *structure* (provable) of the extremal and its
value, plus the exact place the wall sits.

## What IS proven: the value of the construction, and the cyclotomic identities
The construction `{1,…,n−2,(n−1)n}` has `M = n/(n²−n+1)` exactly (verified n=4..14), via:
> **`Φ₆(n) := n²−n+1` is the 6th cyclotomic polynomial.** Mod `q=Φ₆(n)`: `n² ≡ n−1`, hence
> **`n³ ≡ −1`** (so `n` has order **6**), and **`(n−1)n² ≡ −n \pmod{Φ₆(n)}`** (PROVEN: `(n−1)n²≡(n−1)²
> =n²−2n+1≡(n−1)−2n+1=−n`). Therefore at the witness `t=n/q`, the killer `(n−1)n` lands at `−n/q`
> (distance `n/q`) and `v=1` lands at `n/q`; both bind, and the min is `n/q = n/Φ₆(n)`. ∎ (the value).
The 6 powers `n^0..n^5 ≡ {1, n, n−1, −1, −n, 1−n}` (verified n=4,7,14) are the **hexagonal orbit** the
runners resonate on at the witness — the same `6 = φ(14)` of the razor's-edge units.

## The Eisenstein-integer home (the gift)
`Φ₆(n) = (n−ζ₆)(n−ζ̄₆) = N(n−ζ₆)` — the **norm of `n−ζ₆` in the Eisenstein integers `Z[ζ₆]`** (`ζ₆`
primitive 6th root; `Z[ζ₆]=` ring of integers of `Q(ζ₆)=Q(√−3)`).
> **The covering-min denominator is an Eisenstein-integer norm.** `Q(√−3)` is **Heegner (`−3`), class
> number 1** — the gentlest, a PID; its **6 units `{±1,±ζ₆,±ζ₆²}`** are the hexagonal symmetry, and the
> order-6 of `n` mod `Φ₆(n)` is that hexagonal action. The central polygonal numbers `n²−n+1 =
> 7,13,21,31,43,…,183` are exactly the Eisenstein norms `N(n−ζ₆)`.
So the LRC's two arithmetic homes are now both named, and DISTINCT:
| object | field | Heegner | the role |
|---|---|---|---|
| the FLOOR (`Z₇` cyclotomic SOS, set-independent) | `Q(√−7)` | `−7` | the apex-7, the `7`-part of `14` |
| the covering MIN (`Φ₆(n)`, the extremal value) | `Q(√−3)` (Eisenstein) | `−3` | the hexagonal `6`, the pronic split |
They **meet at `7 = Φ₆(3) = N(3−ζ₆)`** — the apex-7 is the Eisenstein norm at `n=3`. (Fun: the floor's
`−7` and the min's `−3` touch at the apex.)

## The killer-fills-the-tooth structure (as far as it goes)
Covering ⇒ a killer `w` (`n∣w`) ⇒ at the `n`-comb witnesses `t=k/n`, `‖wk/n‖=0` — `w` sits on the
observer's empty tooth and FILLS it (last turn's Dirac-comb). The lonely point must leave the `n`-comb.
For the extremal killer `(n−1)n`, it relocates to the **nearest Eisenstein–Farey rational `n/Φ₆(n)`**,
where the proven identity makes `(n−1)n` bind at `n/Φ₆(n)`. **The covering MINIMIZER is the configuration
whose killer is the smallest pronic `(n−1)n` (restoring `n−1` and `n` at once) and whose relocation lands
on the Eisenstein norm `Φ₆(n)`.** That the minimum is no SMALLER than this — i.e. no other covering set
relocates closer to `1/n` — is precisely the `>1/n` wall = LRC(n).

## What this buys the program (honestly)
- **Proven:** the extremal's value `n/Φ₆(n)`; the Eisenstein/cyclotomic identities (`n³≡−1`,
  `(n−1)n²≡−n`); the hexagonal order-6 binding; `Φ₆(n)=N(n−ζ₆)`; the `7=Φ₆(3)` apex meeting point;
  the corrected empirical covering-min (`14/183 < 7/89`).
- **Conjectured (= LRC):** global minimality `M(covering) ≥ n/Φ₆(n)`. The wall is exactly `>1/n`.
- **The reframe it gives:** LRC(n) `⟺` no covering set relocates its lonely point inside the Eisenstein
  radius `n/Φ₆(n)` — a statement about Eisenstein-Farey distances on the `n`-comb, with the floor's
  `Q(√−7)` SOS as its complementary (set-independent) half. **The two Heegner fields `−3` (min) and `−7`
  (floor) frame the conjecture from both sides.**

## Status
- **Verified/PROVEN (opus):** `Φ₆(n)=n²−n+1=N(n−ζ₆)`; `n³≡−1`, order 6; `(n−1)n²≡−n mod Φ₆(n)`; the
  construction value `n/Φ₆(n)` (n=4..14); `7=Φ₆(3)`; Eisenstein/Heegner `−3` home.
- **Honest barrier:** the global lower bound is `>1/n` ⇒ it IS LRC(n); not claimed.
- **New to track:** the Eisenstein-norm covering-min denominators; the two Heegner fields (`−3` min, `−7`
  floor) meeting at the apex `7=Φ₆(3)`; the hexagonal order-6 binding.

Related: the two-pyramids covering-min reflection (this gives its algebraic home), the cyclotomic-self-dual
razor's-edge (`φ(14)=6`), the Z₇-cyclotomic-SOS-floor (`Q(√−7)`), the Dirac-comb (the killer fills the
tooth), HYP-3547 (apex-7 = Heegner), A002061 (central polygonal = Eisenstein norms), THM-523, OPEN-Q-108.
