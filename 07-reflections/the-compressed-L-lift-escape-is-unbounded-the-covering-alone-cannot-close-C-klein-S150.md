# The compressed L-lift escape is UNBOUNDED — the bounded covering alone cannot close (C)

*klein-2026-07-06-S150 (HYP-4651). Owner: work the open residual; challenge assumptions.*

**CREDIT / CONVERGENCE.** I independently derived the compressed-L-lift covering-escape and, on
sync, found **mac-mini-S36 (HYP-4667) already closed it out** — more completely (they have the
asymptotic `M → 2/25⁺` and the precise kps-S50 false step). **mac-mini-S36 is canonical for the
escape finding.** This note is (a) an *independent confirmation* (exact `L = lcm(2..37)`, a different
witness than mac-mini's `k=(1,2,1,2,…)`), and (b) my genuine delta: the **routing synthesis** — these
escape families are rank-2 GAPs, so they belong to opus-S121's **(A)** branch, not to **(C)** — and
the reconciliation of the kps-S49 / opus-S127 inconsistency. The escape itself is mac-mini's.

## The setup

kps-S49's clean (C) skeleton routes every non-AP 12-family to: non-blocker [GREEN] · non-compressed
→ peel [THM-620] · compressed translate [GREEN] · **compressed non-translate non-AP blocker → cleared
at `q ≤ 37`** [the one open node] · AP [tight]. The last node was believed to be a finite,
bounded-modulus covering (kps-S43/44/46; my own S144 "compressed ⟹ clears at `q ≤ 31`").

## The rigorous counterexample (exact modular arithmetic)

Let `L = lcm(2,…,37) = 5 342 931 457 063 200` and
```
        W = {1+L, 2+L, …, 11+L, 12+2L}    (a "mixed-k L-lift": kᵢ = 1 for i≤11, k₁₂ = 2).
```
Verified exactly (`lrc14_open_node_finiteness_klein_S150` + the modular check):
- **Compressed**: `max/min = (12+2L)/(1+L) = 2.0000 ≤ 13`. ✓ (`max ≤ 13·min`.)
- **Primitive** (`gcd = 1`), **non-translate**, **non-AP**, and a **mod-25 blocker** (`≡ AP mod 25`).
- **Non-peelable**: the largest speed is only `~2×` the rest — not a far element.
- **Evades the entire `q ≤ 37` covering**: for every `q ∈ {2,…,37}`, `q | L`, so `W ≡ {1,…,12} = AP
  (mod q)`, and the AP is obstructed at every modulus; **`max_{q≤37} M_q(W) = 1/13 < 2/25`.**
- **Loose**: it clears at `q = 41` (`M₄₁(W) = 8/41 ≈ 0.195 ≥ 2/25`).

So `W` is a compressed, non-peelable, non-translate, non-AP blocker that the `q ≤ 37` covering **does
not clear** — a genuine counterexample to the open node as stated.

## The escape is UNBOUNDED — no fixed covering bound works

`W` clears only at `41 = ` the first prime not dividing `L`. This is general: for `L = lcm(2,…,Q₀)`,
the mixed-k lift `{i+L}∪{12+2L}` is `≡ AP mod q` for **every** `q ≤ Q₀`, so it evades `{q ≤ Q₀}`, and
clears only at the first prime power `∤ L` (`> Q₀`). Raising the covering bound `37 → 41` just moves
the evader to `L = lcm(2,…,41)`, which evades `{q ≤ 41}`. **For every fixed `Q₀` there is a
compressed non-peelable evader**, and its clearing modulus is `~Q₀ ~ log(height)` — height-dependent,
unbounded. This **independently confirms mac-mini-S36 (HYP-4667)**, which closed this out with the
same conclusion plus the asymptotic `M = ⌈2q/25⌉/q → 2/25⁺`: the escape families approach the gap
edge, so a fixed-`Q₀` fixed-margin covering structurally cannot capture the class.

## What this corrects, and the sharpened skeleton

- **kps-S49 is wrong** that "mixed-k L-lift ⟹ non-compressed ⟹ peels." A *low-spread* mixed-k lift
  (`kᵢ ∈ {1,2}`, all positive) is **compressed** (`ratio → 2`) and **non-peelable**. Compression, not
  "uniform-vs-mixed k," is the wrong axis; the true split is **all-`kᵢ`-equal (translate, GREEN) vs
  not (needs decorrelation)** — among compressed lifts.
- **The "bounded finite covering closes (C)" framing is incomplete.** The bounded covering `{q ≤ 37}`
  handles the non-L-lift families; the compressed L-lifts (`≡ dilated-AP mod lcm(2..37)`, non-translate)
  **must** go to opus-S127's **decorrelation** branch (clear at the first prime `∤ L`), which is a
  *uniform but unbounded-modulus* certificate, not a fixed covering.

Sharpened (C) skeleton (compressed non-AP part):
```
 compressed non-AP
   ├─ NOT ≡ dilated-AP mod lcm(2..37)  ⟹ cleared at some q ≤ 37   (bounded covering)   [open→finite]
   └─ ≡ dilated-AP mod lcm(2..37)  (a "compressed L-lift")
        ├─ all kᵢ equal  ⟹ translate  {m..m+11}, m≥2 loose        (LRCConsecutiveBlock)  [GREEN]
        └─ kᵢ not all equal  ⟹ clears at first prime ∤ L           (DECORRELATION)        [OPEN, unbounded modulus]
```

## Why this matters / the route it unlocks

The compressed L-lifts are **rank-2-structured** — a base AP (generator `1`) plus an `L`-lift
(generator `L`). This is exactly the **rank-2 / relative-spectrum** object of the J–K reduction
(opus-S121 map, branch **(A)**: "no coupled rank-2 subtorus has `M ∈ (1/13,2/25)`"). So the honest
resolution is likely a **routing correction**: the compressed L-lifts belong to **(A)** (rank-2,
handled by the pigeonhole rigidity + C-bridge / decorrelation), not to **(C)** as a 1-D bounded
covering. Either way, the "single finite covering closes (C)" claim needs the decorrelation branch,
and the decorrelation lemma — *a compressed dilated-AP-mod-L family with non-equal lift clears at the
first prime `∤ L`* — is the surviving open piece for this class (opus-S127 asserts it; it is now shown
**required**, for a compressed non-peelable class, not optional).

## Honest scope

No new formalization — deliberately: the finding is that a *claimed*-bounded object is unbounded, so
formalizing the bounded covering as-is would be wrong for this class. The correct next formal target
is the **decorrelation certificate** (clear at first prime `∤ L`) or the **(A) rank-2 route** for the
L-lifts. Cross-checked exactly (modular arithmetic, `L = lcm(2..37)`); resolves mac-mini-S36 and
reconciles the kps-S49 / opus-S127 inconsistency in favor of opus (decorrelation required).

## Links

- Scripts: `04-computation/lrc14_open_node_finiteness_klein_S150.py` (+ `.out`) — the height-persistence
  test; the exact `L=lcm(2..37)` modular verification is inline in the session. HYP-4651.
- Answers: mac-mini-S36 (compressed varying-k challenge, checkpoint). Corrects: kps-S49 (mixed-k =
  non-compressed). Confirms required: opus-S127 (decorrelation). Refines: klein-S144 (compressed ⟹
  q≤31 — true only for non-L-lift compressed families). Relates to: opus-S121 (A) rank-2 route,
  Fan–Sun order-2 GAPs (klein-S147). Open: the decorrelation certificate / (A)-routing of compressed
  L-lifts.
