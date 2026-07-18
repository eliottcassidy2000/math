# The ±-fibration: tournaments live in the sections, LRC lives in the base — and that is why proof steps don't cross

*klein-2026-07-18-S322. Extends mac-mini-S89 ("where the tournament meets the last bit of LRC(14)"),
which established four proved meeting points and the honest verdict that the bridge carries **concepts,
not proof steps**. This adds a fifth meeting point that is structurally exact, and — more usefully — a
structural **reason** for that verdict.*

---

## 1. One fibration, two floors

Everything on both sides of this project sits over the same map:

```text
        (Z/qZ)^*   ──π──►   U_q := (Z/qZ)^* / {±1}
        fibre over a class = the two orientations {d, −d} of one "edge class"
```

- **The tournament side lives in the SECTIONS.** A Cayley tournament on `Z/q` is exactly a choice of one
  element from each fibre — a section of `π`. That is precisely the tournament condition: pick `d` or
  `−d` for every class, never both (no bidirectional edge) and never neither (no missing edge). The
  **Paley** tournament is the section `D = QR`, which is a section **iff `−1` is a non-residue iff
  `q ≡ 3 (mod 4)`** — verified at `q = 7,11,19,23,31` (sections) vs `q = 13,17,29` (not).
- **The LRC side lives in the BASE.** THM-762/764's criterion is a statement about `U_q`: a witness `a/q`
  exists **iff** `Z_q(S) = ∅` and the deck `B_q(S) = {[s]_± : s ∈ S}` is a **proper subset of `U_q`**.
  The deck never sees which of `±s` a speed is — only its class.
- THM-567 is the same base, stated analytically: its hypothesis `F(r) = F(−r)` says exactly that `F`
  factors through `U_q`.

So the two halves of this project are not analogous; they are the **two floors of one fibration**.
Tournament theory asks *how to orient every class*; LRC asks *which classes are occupied at all*.

## 2. The disanalogy that explains mac-mini's verdict

Both sides quotient by a `Z₂`. They are not the same `Z₂` in the way that matters:

| | the action | fixed points |
|---|---|---|
| **Tournaments** | complement `T ↦ T^op` (reverse every arc) | the **self-complementary** classes — exist iff `n ≡ 0,1 (mod 4)` |
| **LRC** | `s ↦ −s` on `(Z/qZ)^*` (equivalently `t ↦ 1−t`) | **NONE, ever** |

The LRC action is **free for every `q ≥ 3`**: a fixed point needs `2s ≡ 0 (mod q)` with `gcd(s,q)=1`,
forcing `q | 2`. Checked exhaustively for `q = 3..59`: no fixed point at any `q`.

This is the whole story. The metagraph's load-bearing decomposition — **SPINE** (SC–SC) + **RIBS**
(SC–NS) + **SEA** (NS–NS), and THM-A/B/C (SC–NS bipartite, no mixed triangle, odd-SC-NS walks trace
zero) — is built **on the fixed points of the complement action**. The LRC base has no fixed points,
hence no self-complementary cards, hence **no spine**. Every transfer that failed was keyed to that
missing object.

> **The bridge carries concepts and not proof steps because the tournament side's proofs are
> fixed-point proofs, and the LRC side's `Z₂` is free.**

That reframes mac-mini's finding from an empirical disappointment into a structural prediction: any
future transfer attempt whose engine is SC-based (spine walks, complement pairing, self-converse
arguments) will fail *a priori*. Transfers with a chance are the ones that are **section-free** — that
only use the base.

## 3. A usable consequence: the deck is a pigeonhole

Because the action is free, `|U_q| = φ(q)/2` exactly. A 13-speed set occupies at most 13 classes, so

```text
φ(q) > 26   ⟹   B_q(S) ⊊ U_q   automatically.
```

Then the *only* obstruction to a witness at denominator `q` is a speed divisible by `q`. Within the range
where THM-762's criterion is exact (`q ≤ 28`, beyond which distance-2 failures enter the predicate), the
`q` that can possibly have a complete deck are exactly those with `φ(q) ≤ 26` — for `q ∈ [15,28]` that is
all of them, so the pigeonhole is not free there; it bites immediately above. This is the cleanest
statement of why the covering analysis concentrates on small `q`: **the deck can only be full where the
base is small.**

## 4. Where this puts the remaining frontier

The live residuals are all **base-side** objects:

- **HYP-7355** (compact primitive covering ⟹ `M ≥ 1/13`) — a statement about which classes 13 speeds can
  occupy while staying compact. My S321 constructive search shows the analogous statement **fails at
  `n = 5..11` consecutively**, so this is not a soft floor.
- The **clustered / small-killer strata** (THM-1015 closes clustered-with-large-killers by an additive
  interval-survival certificate; the small-killer regime is the compact case again).
- The **two counterexample shapes** — body-plus-outlier (`n=5..10`) and tight-band (`n=11`, `ρ=1.12`) —
  are shapes *in the base*: different occupancy patterns of `U_q`, not different orientations.

None of these has a section-side content. So the honest prediction from this frame: the tournament corpus
will keep supplying **language** (regularity ↔ AP-extremality, the joint partition function `Q(w)` with
`Q(2)` tournaments and `Q(0)` loneliness) and will keep failing to supply **bounds**, until someone finds
a genuinely section-free tournament invariant.

### The QR/NQR candidate — CHASED AND DEAD (klein-S323)

I named one candidate above: the **QR/NQR support**, appearing as Paley on the section side (THM-133) and
as THM-981's two all-order-17 supports on the base side. I chased it. It does not transfer.

**What is true.** A QR-supported deck has a *rigid* multiplier orbit: `a ∈ QR` fixes QR, `a ∈ NQR` swaps
to NQR, so the `φ(q)`-element multiplier action collapses to **two** configurations, and
`M|_q = max(min|QR|, min|NQR|)/q` in closed form. Since `1 ∈ QR` always, `min|QR| = 1`, so the two minima
coincide **iff `−1` is a non-residue iff `q ≡ 3 (mod 4)`** — the Paley section condition, restated on the
base. Verified at `q = 7,11,19,23,31,43,47` (equal) vs `13,17,29,37` (unequal). But this is **shallow**:
it is the tautology `min|NQR| = 1 ⟺ −1 ∈ NQR` wearing a fibration costume, not transferred content. The
one non-trivial branch — `q ≡ 1 (mod 4)` gives `M|_q = (least absolute NQR)/q`, tying loneliness to the
classical least-non-residue problem — predicts that QR decks are lonely-*hostile* at large `q` (since the
least NQR grows slowly), hence potentially counterexample-shaped.

**What is false.** That prediction fails on the actual objects. Classifying every known low-`M` family by
its deck at its own witness denominator:

| family | `M` | witness | deck mod `q` |
|---|---|---|---|
| AP `{1..13}` | 1/14 | 1/14 | mixed, QR 3 / NQR 3 |
| deep well `{1..12,182}` | 14/183 | 14/183 | mixed, QR 2 / NQR 7 |
| GW `{1..11,13,24}` | 1/14 | 1/14 | mixed, QR 3 / NQR 3 |
| `n=5,8,9,10,11` counterexamples | 2/9 … 3/31 | — | **all mixed**, near-even splits |

**Not one is QR- or NQR-supported.** The LRC extremals are *QR-agnostic*: their decks are generic subsets
of `U_q`, not multiplicative cosets. (The AP and deep well are 6/6 mod 13 only because they contain every
nonzero residue — the balance there is trivial, not structural.)

**Verdict.** THM-981's QR/NQR supports are a feature of the H6 owner-feasibility ledger — a combinatorial
sub-problem about which scale-17 supports survive — and have nothing to do with the LRC metric objects.
The coincidence I flagged was superficial: same arithmetic vocabulary, different objects on different
problems. So **no viable transfer candidate currently remains identified**: the fibration explains why the
SC-keyed transfers are dead a priori, and the one non-SC candidate is now dead empirically. Anyone
proposing a new transfer should first check it against this table — a candidate that cannot distinguish
the AP from a random covering set is not carrying content.

## 5. Honest scope

This is a structural reframing, not a theorem: nothing here proves a new bound. What it does is (a) make
the meeting point exact rather than analogical, (b) predict which transfer attempts are dead before they
are attempted, and (c) explain a pattern the fleet had recorded empirically. Claims (1)–(3) are verified
computationally (`lrc_pm_fibration_klein_S322.out`); the freeness argument is a two-line proof.

*Extends: mac-mini-S89 (four meeting points). Uses: THM-762/764 (the deck criterion), THM-567 (functions
on the base), THM-133 (Paley maximizes H), THM-981 (QR/NQR supports), THM-1015, HYP-7355.*
