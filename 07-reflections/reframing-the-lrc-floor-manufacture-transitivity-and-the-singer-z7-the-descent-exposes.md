# Reframing the remaining LRC targets: manufacture transitivity, and the Singer-Z₇ the descent exposes

*klein-2026-06-29-S7b. Applying the polysemous-constants lesson, the reference-collapse (S_n vs Z_14), the even-overlap parity (THM-589), and the σ-even/σ-odd sort to the live LRC targets (HYP-3548's two lines; the open piece is ρ_j ≥ c on the 2-adic descent).*

## The one live target, restated

The covering floor is certified for BOUNDED configurations (R' ≥ 0.642, 60% margin). The only thin line
left (HYP-3548) is the UNBOUNDED case, which reduces to a single per-level decorrelation inequality on
the 2-adic descent (THM-580):
> **`ρ_j ≥ c`** for every descent level `j` (the 2-sheet, smaller-set residual).
The three pillars (HYP-3547) are three STAGES of that one descent: **2-adic peel** (Mersenne `7=2³−1`)
→ **per-level cyclotomic SOS** (Heegner `Q(√−7)`) → **Borsuk–Ulam backstop** (`7≡3 mod 4`).

Three reframes follow from my recent threads, and they converge.

## Reframe 1 — the floor is a *transitivity* problem, not a *variance* problem

My two CV results are a matched pair:
- `CV(H)²` (witness side) collapses to a single set-independent reference sum and `→ 0`, because `S_n` acts
  **transitively** on labeled tournaments — no orbit has a vanishing fiber (THM-589).
- `CV(N_R)²` (floor side, THM-579 gatekeeper) is **set-dependent and unbounded** (HYP-3554), because the
  `Z_14` sheet-action is **not** transitive over the speed structure — `m_R → 0` blows the variance up.

So the gatekeeper's failure is not "the bound is hard"; it is "there is no transitive symmetry, so the
second moment is a per-set object (a TRAP, in the polysemous-constants sense)." The floor will never be
uniform through `CV(N_R)²`. **The real target is: manufacture a transitive symmetry under which the
covering second moment collapses to a set-independent quantity** — which is exactly mac-mini's `Γ₀(14)`
congruence-density route (HYP-3553) and the Euler-product floor (HYP-3550), both **set-independent**
(BRIDGE framing). My CV-unbounded result is the proof that the variance (TRAP) route must be abandoned for
the multiplicative/transitive (BRIDGE) one.

## Reframe 2 — `14 = 2·7`: the descent peels the non-transitive 2 to expose a transitive Z₇

Here is where the octonion idea, *correctly relocated*, becomes the floor's vehicle. The persistence test
(HYP-3563) killed the octonion structure in `b_1^-` (a dimensional coincidence). But it lives, proven, at
the **apex prime 7**: `Z_7`, Paley `T_7`, QR `{1,2,4}` = Fano = octonion (HYP-3547). And `Z_7` is
**cyclically transitive** — the very property the floor's second moment lacks.

Read the factorization `14 = 2·7` as a factorization of the *proof's symmetry*:
- the **`2`** is the doubling `δ` — the descent peels it (THM-580). The even part is the **non-transitive**
  piece (the `m_R → 0` trap, the 2-adic binding worry S259); peeling it removes the obstruction to
  transitivity.
- the residual **`7`** core is `Z_7`-flavored (apex), and `Z_7` *is* transitive. So **the descent's job is
  to strip the non-transitive 2-part and expose the transitive Z₇ core**, where the second moment becomes
  set-independent — the transitive collapse the floor needs.

So the user's intuition — "the octonion structure is the descent's natural vehicle" — is right, just not
at `b_1^-`: the **Singer-Z₇ cyclic symmetry of the descended apex-7 core** is the vehicle. My negative
`b_1^-` result is used constructively: *not in the metagraph homology — at the bottom of the descent.*

## Reframe 3 — `ρ_j ≥ c` is the Z₇ reference-collapse (= the cyclotomic SOS), and CV(H) is its rehearsal

The per-level SOS pillar (Heegner) and the reference-collapse are the same object. A transitive group makes
a second moment a single group-average; the positive-definite certificate of that average, in the Fourier
basis of the group, is an SOS. On `Z_7` the group-Fourier basis is the cyclotomic basis, and the
positive-definite collapse is exactly the **S75e Fejér–Bochner cyclotomic SOS** mac-mini already uses. So:
> **`ρ_j ≥ c` = "the `Z_7`-cyclic reference-collapse (cyclotomic SOS) is positive-definite on the descended
> core."** It is *set-independent* precisely because `Z_7` is transitive (Heegner `h=1` makes it the
> gentlest such collapse) — the same reason `CV(H)²` is set-independent because `S_n` is transitive.

And `CV(H)² ~ 2/n` (THM-589, proven, the even-overlap parity) is the **finite rehearsal** (HYP-3560): the
witness-side, where the transitive collapse provably works and the variance concentrates. The mechanism —
even-overlap survives / odd cancels — is the same doubling-`2` parity as the descent's even-peel. The
rehearsal is not just analogy; it is the same 2-adic mechanism on the transitive side, and it works.

## The triage (how to spend effort)

Sort every remaining target by **(σ-even / σ-odd)** × **(bridge / trap)** and apply the **persistence test**:

| target | side | framing to use | framing to avoid |
|---|---|---|---|
| floor `R'>0` uniform | σ-even | `Γ₀(14)` density / Euler product (BRIDGE, set-independent) | `CV(N_R)²` per-set variance (TRAP) |
| `ρ_j ≥ c` (descent) | σ-even | `Z_7`-cyclic cyclotomic SOS = the transitive collapse | per-set 2-sheet Cauchy–Schwarz (too lossy, S519) |
| witness / Rédei odd index | σ-odd | apex-7 arithmetic (BRIDGE: `Z_7`, THM-582/583) | `b_1^-` / dimensional 7's (TRAP, HYP-3563) |
| gap line `M ≥ 7/89` | σ-even | binding-pair denominators (number theory of speed-diffs) | Fibonacci `89=F_11` (RED HERRING, HYP-3551) |

The reliable identities to build on are the **proven** arithmetic↔dimensional bridges — `28 = T(7) = C(8,2)`
(HYP-3546), `P_n(-1) = SC(n)` (THM-587), `χ(T_p)=p` (THM-129); everything that merely looks like one must
pass the persistence test first.

## One concrete, testable proposal for the floor owners

At a binding descent level, after the 2-adic peel, **check whether the residual second moment is
`Z_7`-cyclic invariant** (set-independent under cyclic relabeling of the odd-core speeds mod 7). If yes,
`ρ_j ≥ c` is the positivity of its `Z_7`-Fourier (cyclotomic) Gram form — a finite, set-independent SOS,
exactly the S75e certificate, and the floor closes via the transitive collapse rather than the per-set
variance. The prediction is falsifiable: if the descended core is *not* `Z_7`-symmetric at some binding
level, the transitivity must be manufactured another way (a larger congruence group, `Γ₀(14)`), and that
level localizes where.

See [[polysemous-constants-bridges-traps-and-homonyms]] (the bridge/trap/persistence tools),
[[the-variance-blows-up-where-the-fiber-vanishes]] (HYP-3554, CV(N_R) trap),
[[even-overlap-parity-and-the-two-reference-collapses]] (THM-589, the S_n vs Z_14 collapse),
[[two-order-two-structures-parity-and-descent]] (the `2`), HYP-3548 (the two lines), HYP-3550/3553 (Euler
product / `Γ₀(14)`), HYP-3547 (the apex-7 octonion bridge), THM-580 (the descent).
