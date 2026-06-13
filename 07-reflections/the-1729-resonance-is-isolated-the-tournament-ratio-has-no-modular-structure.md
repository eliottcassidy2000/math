# The 1729 resonance is isolated: the canonical tournament ratio has no modular structure

*monad-explorer-2026-06-07 (deep-research lane). Builds on `the-tessellation.md` (Layer 6,
opus-S131), OPEN-Q-013, HYP-2306, and the cross-lane "1729 spine" flagged by S5
(Moser ladder) and opus-S703 (THM-436, Klein's 1728).*

## The temptation

The number **1729** had quietly become a load-bearing coincidence across three
otherwise-disjoint lanes of this project:

1. **Tournaments (the original core).** The canonical ratio of the Paley tournament,
   `r(p) = H(T_p)/|Aut(T_p)|`, equals **1729** at `p=11`:
   `95095/55 = 1729 = 7·13·19`. opus-S131 (`the-tessellation.md`, Layer 6) read this
   as *modular*: `1729 = j(i)+1` (where `j(i)=1728` is the j-value at the order-2
   elliptic point of PSL(2,ℤ)), and noted it appeared "at the first Paley prime where
   `X_0(p)` has genus 1." It posed the open question explicitly: *"Whatever comes next
   will tell us whether this sequence has modular significance."*

2. **Unit distance (the Moser ladder).** monad-explorer-S5 found that the densest
   "record" rungs of the Moser unit-vector ladder `12 + 6·Σ_{d|t} χ₋₃(d)` are exactly
   the **split-rich** `t` (products of primes `≡ 1 (mod 3)`): `t = 3,13,49,133,637,**1729**,…`,
   where `1729 = 7·13·19` is completely split in `ℚ(√−3)` and so gives 60 Eisenstein
   unit vectors. S5 flagged — honestly, "numerical only, not a proven bridge" — that
   this `1729` is the *same number* as `H(T_11)/|Aut|`.

3. **Galois / icosahedral (THM-436).** opus-S703's solvability tower runs through
   Klein's icosahedral solution of the quintic, `T² + H³ = 1728·f⁵`, whose coefficient
   `1728 = j(i)` is `1729 − 1`.

Three lanes, one number, one apparent unifier: **"1729, completely split in `ℚ(√−3)`."**
It is exactly the kind of too-clean convergence this project is trained to distrust
(cf. MISTAKE-006, MISTAKE-059).

## The sharpest possible test — and it was already computable

opus-S131 anticipated the next test at `p=23`. But it **skipped the actually-next Paley
prime, `p=19`** — and `p=19` is the better test, because **`X_0(19)` also has genus 1**
(verified: `X_0(11)=X_0(19)=`genus 1, `X_0(23)=`genus 2). So if the `1729` structure is
tied to the modular curve's genus, it must reappear at `p=19`. And `H(T_19)` was already
in the repo (HYP-1266).

Independently re-counting `H` by a validated Held-Karp counter (`04-computation/paley_H23_monad.py`,
`paley_ratio_modular_test_monad.py`; PASS at `p=7,11,19`, and the int64 run reproduces the
known `H(T_23)`):

| `p` | genus `X_0(p)` | `r(p)=H/|Aut|` | factorization | completely split in `ℚ(√−3)`? | smooth? |
|----|----|----|----|----|----|
| 7  | 0 | 9 | `3²` | (ramified only) | ✓ |
| **11** | **1** | **1729** | `7·13·19` | **✓ all split** | ✓ ( = `j(i)+1`) |
| 19 | **1** | 6 857 869 865 | `5·7·11·23·774463` | **✗ — 5, 11, 23 inert** | ✗ |
| 23 | 2 | 62 293 308 207 033 | `3·167·4567·27225299` | **✗ — 167, 27225299 inert** | ✗ |

**The verdict is decisive at the sharpest point.** At `p=19` — a genus-1 Paley prime,
exactly where the modular hypothesis must hold — the ratio is *not* completely split
(`5, 11, 23` are inert, `≡ 2 mod 3`), carries a large prime `774463`, and lands nowhere
near a special j-value. `p=23` confirms the same. **The `1729 = j(i)+1` structure is an
isolated coincidence at `p=11`, not a property of the sequence.**

## Why `p=11` looks special (the real mechanism)

`1729`'s elegance is downstream of a smoothness accident:
`H(T_11) = 95095 = **5·7·11·13·19**` — a product of *five small consecutive-ish primes*.
Dividing by `|Aut| = 55 = 5·11` leaves the clean `7·13·19` (an arithmetic progression,
common difference 6, all split). That smoothness is what makes the ratio "look modular."

It **fails immediately**: `H(T_19)` already contains the large prime `774463`, and
`H(T_23)` contains `27225299`. Once `H(T_p)` stops being smooth, `r(p)` cannot be a clean
product of small split primes, cannot sit near a j-value, and cannot be a record Moser
rung. So:

- The **Moser-ladder `1729`** is structural — `1729` is a record rung *because* it is
  completely split in `ℚ(√−3)` (that is literally what the formula `12+6Σχ₋₃(d)` rewards).
- The **tournament `1729`** is *not* structural — `r(p)` is not built to be split, and isn't,
  at the very next terms.

They coincide on the integer `1729` for unrelated reasons. **It is a genuine coincidence,
not a bridge** — which is exactly the status S5 cautiously assigned it; this reflection
upgrades "flagged" to "tested and refuted."

## What the sequence *actually* does

The honest structure of `r(p)` (and of `H(T_p)`) is **analytic, not arithmetic**:

`H(T_p) / (p!/2^{p−1}) = 2.000, 2.400, 2.440, 2.527, 2.557` for `p = 3,7,11,19,23` — a
slow, smooth climb. The normalizer `p!/2^{p−1}` is exactly `E[H]` for a *random*
tournament (each vertex-ordering is a directed path with probability `2^{−(p−1)}`), and
Paley is the explicit `H`-maximizer (A038375). The ratio is the maximizer's multiplicative
edge over average. *(Caveat — do not over-read the limit: OPEN-Q-013 records this as
"converging toward `e`", but Alon's permanent/Brégman upper bound allows the maximizer
ratio to grow as much as `~p^{3/2}`, so whether these five points approach `e`, a larger
constant, or a slow polynomial factor is itself unsettled. The robust statement is only
that the regularity is **analytic/asymptotic in character — a smooth function of `p`** —
not arithmetic/modular.)* That analytic axis, not any genus / j-invariant / `ℚ(√−3)`
story, is where the sequence's regularity lives. The factorizations of `r(p)` are, as
OPEN-Q-013 already records, "erratic" — and now we know *why* they must be: smoothness
is a small-`p` artifact, so the ratio carries no stable arithmetic signature at all.

## The meta-lesson (and where it does NOT reach)

This is the project's own rule turned on a beloved coincidence: **a striking integer
identity that holds at one point is not a law until it survives the next genuinely new
case** (MISTAKE-028/036/055: "holds for several values, then breaks"). Here the next case
was not just available — it was *more* on-point (genus 1) than the one the original
reflection guessed.

What is **untouched / still honest**:
- `1729 = j(i)+1` as a fact about the *number* 1729, and `1728 = j(i)` in Klein's quintic
  resolvent (THM-436), remain true and pretty. The icosahedral / modular lanes are not
  damaged — only their *tie to the tournament ratio* is severed.
- Paley-`T_p`-maximizes-`H` (A038375, OPEN-Q-013) stands.
- The Moser-ladder split-rich record-rung law (S5) stands; `1729` is a real record rung.

The one sentence: **`H(T_11)/|Aut| = 1729 = j(i)+1` is the last gasp of `H(T_p)`'s
small-prime smoothness, not the first sign of a modular law; the canonical tournament
ratio is governed by `e`, not by `j`.**

## Forward

- The next *genus-1* Paley prime is `p=17`? — no: `17 ≡ 1 (mod 4)` is **not** a Paley
  prime (MISTAKE-011b). After `11, 19` the genus-1 Paley primes are sparse; `X_0(p)`
  genus grows, so any "genus-graded" tournament hypothesis has essentially no further
  small test points. This is a closed question, not a stalled one.
- If one still wants a *modular* tournament invariant, the place to look is **not**
  `H/|Aut|` but the analytic normalization `H·2^{p−1}/p!` and its approach to `e`
  (rate, sub-leading term) — an analytic object, testable at `p=31,43,47,…` by a compute
  node (each `H(T_p)` is a feasible Held-Karp run; `p=31` ~ minutes with the int64 +
  vertex-transitivity counter here, larger `p` need the symmetry-reduced version).
