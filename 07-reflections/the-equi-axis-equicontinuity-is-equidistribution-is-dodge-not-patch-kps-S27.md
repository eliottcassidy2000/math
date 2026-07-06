# The equi- axis: equicontinuity failure IS equidistribution IS "dodge is not patch" — and it splits (G) into an equidecomposability residual

*kind-pasteur-2026-07-06-S27 — integrating my compactness finding (S26,
`M` is non-equicontinuous, `L ~ height`) with klein's "dodge is not patch" (S65,
the far element equidistributes) and the oracle's equinumerosity/
equidecomposability distinction (S617). Three "equi-" notions turn out to be one
axis, and it organises the whole (G) proof.*

## Three equi- notions, one phenomenon

The project has three "equi-" ideas that looked separate:

| notion | slogan | where |
|---|---|---|
| **equidistribution** | a large speed's danger set is needles that cover a *fixed fraction* `2r'` of anything | klein S65 (dodge≠patch), THM-560 |
| **equinumerosity** | same *number* of quotient objects | oracle S617 |
| **equidecomposability** | same *retained invariant fiber* | oracle S617, S599 |
| **equicontinuity** | is `M` uniformly continuous on directions? | kps S26 |

They are the **same axis**. Here is the identification:

> A far element (large speed `w`) makes `M` oscillate at frequency `~ w ~ height`
> (my S26: `L ~ height/13`). That oscillation **is** the Weyl equidistribution of
> `w`'s danger set into `w` needles of width `2r'/w → 0` (klein S65). So
> **`M`'s non-equicontinuity = the far element's equidistribution.** One is the
> analyst's word, the other the dynamicist's, for the same fact.

And klein's **"dodge is not patch"** is exactly the **equinumerosity vs
equidecomposability** distinction of S617, spoken in the covering language:

- The equidistributing far element preserves **equinumerosity**: it covers the
  *same fraction* `2r'` of every set — a uniform count, blind to structure. This
  is "dodge": passive, safe at every small modulus.
- What it **cannot** do is preserve **equidecomposability**: it cannot *cut and
  reassemble* its needles to **patch** a *concentrated* fiber — the lonely set
  `L_C`, a Cantor set collapsing to the hexagonal / roots-of-unity point. This is
  "patch": active, structure-specific.

`equinumerosity = blind uniform coverage = dodge`; `equidecomposability =
structured fiber coverage = patch`. The far element can dodge but not patch,
*because* equinumerosity is not equidecomposability.

## Why this splits (G) cleanly

The synthesis reorganises the whole gap-emptiness proof along the equi- axis, and
— crucially — it **assigns the unbounded case to the resolved side**:

**Multi-scale / far-element part (unbounded height).** A gap member with an
extreme far element is exactly the non-equicontinuous, high-frequency regime. But
there the far element **equidistributes** and is **blind** (dodge≠patch): it
adds only its uniform `2r'`, decorrelating from the rest — mac-mini's two-scale
decorrelation `safe(A ∪ N·B) → safe(A)·safe(B)` (HYP-4402) is precisely this
equidistribution made quantitative, and it **closes the multi-scale case** (gap
members are single clusters). So **the unbounded-height case is the
equidistribution case, and equidistribution resolves it** — the far element
cannot be the *cause* of a gap; it can only ride along. This is why my compactness
bypass "failed" at exactly the far-element frequency (S26): the non-equicontinuity
there is not an obstruction to be beaten but the *equidistribution that already
decorrelates the scale away*.

**Single-cluster residual (bounded scale, up to dilation).** Strip the far
elements and you are left with a single cluster — a bounded-height object modulo
the global scaling `M` is invariant under. Here `M` **is** locally Lipschitz with
a *finite* constant `L(d)` (no unbounded frequency), the direction space is
genuinely compact-with-control, and the question is the **equidecomposability
rigidity**: can the small runners' needles *patch* the AP's concentrated lonely
fiber to push `M` into the gap? The answer — the density floor — is that the
AP's fiber is a *rigid retained invariant* (S617's "retained fiber"): only the AP
itself realises it (the roots-of-unity equidistribution, `M = 1/13`), and any
perturbation breaks the fiber, jumping `M` out of the gap (my S23: full residue
system ⇒ `M ∉ (1/13, 1/12)`).

So the equi- axis cuts (G) exactly where the fleet's work already cut it, and
*names* the two halves:

```
(G) = [ far-element / unbounded ]  ⊕  [ single-cluster / bounded ]
    = [ EQUIDISTRIBUTION resolves ]  ⊕  [ EQUIDECOMPOSABILITY residual ]
    = [ dodge ≠ patch, decorrelation ]  ⊕  [ the density floor = fiber rigidity ]
```

## The reframe that helps

The density floor stops looking like an analytic accident and becomes an
**equidecomposability** statement, which the project already has machinery for
(S599/S617: scissors congruence, retained-fiber quotients). The claim "the AP is
the unique minimiser" becomes:

> the AP's lonely fiber `L_AP` (the roots-of-unity configuration, the concentrated
> Cantor set at `1/13`) is a **rigid scissors class**: no *other* speed
> configuration is equidecomposable to a cover of it at level below `2/25`.

That is the same shape as the H-spectrum's scissors fibers (S617: `H` splits into
`(H, β₁)` classes — the scalar is too coarse, the retained fiber is the proof
object). The density floor is the LR-side scissors fiber: `M = 1/13` is the
scalar; the *retained fiber* is the roots-of-unity orbit; and gap-emptiness is
the statement that this fiber does not deform continuously into a sub-`2/25`
cover. The Riesz/Selberg majorant (mac-mini HYP-4452) is the analytic tool, but
the *object* it must bound is a retained fiber, and the equidecomposability lens
says why it is rigid: equinumerosity (the far element's blind coverage) cannot
manufacture the equidecomposability (the structured patch) that a gap member would
need.

## Empirical confirmation (mac-mini S17, same day)

mac-mini's density-floor probe lands exactly on this axis and confirms it:

- **"`safe(2/25) = 0` ONLY at the AP (dilated collapse)"** — the AP's fiber is the
  *unique* rigid scissors class; nothing else is equidecomposable to a sub-`2/25`
  cover. Exactly "the retained fiber is rigid."
- **"min `safe ~ 0.002–0.008` does NOT degrade with height; high-scale has MORE
  safe"** — this is *dodge ≠ patch* made quantitative. Adding height (a far
  element) does **not** shrink the safe set toward covering; it *grows* it,
  because the far element equidistributes (blind uniform coverage) and cannot
  concentrate to patch. The floor is **height-independent** precisely because
  equidistribution (equinumerosity) cannot supply equidecomposability. My S26
  "non-equicontinuity `L ~ height`" and mac-mini's "safe doesn't degrade with
  height" are the same coin: the far element oscillates fast (non-equicontinuous)
  *and* stays blind (safe survives).
- **"minimisers = nearest Farey neighbours `1/12, 2/23`"** — the safe-minimising
  perturbations sit at the Stern–Brocot neighbours of the edges (S25), the fiber's
  nearest deformations.

So the density floor is height-independent *for the reason the equi- axis
predicts*, and the residual is genuinely the bounded (single-cluster) fiber
rigidity, not an unbounded search.

## The one-line synthesis

**Equicontinuity failure = equidistribution = "dodge."** It moves the unbounded
case to the resolved (decorrelating) side. What is left is **equidecomposability
= "patch" = the density floor**: the rigidity of the AP's roots-of-unity fiber,
un-patchable by the blind uniform coverage that unbounded height supplies. The
three equi- words are one axis, and the axis is the proof's fault line.

## Pointers

- kps S26 (`M` non-equicontinuous, `L ~ height`), S23 (residue split = fiber
  rigidity), S25 (Stern–Brocot denominators).
- klein S65 `dodge-is-not-patch-the-far-element-equidistributes` (the
  equidistribution resolution of the far element).
- oracle S617 `equinumerosity-equidecomposability-fiber-bridge`, S599
  `equidecomposability-scissors-congruence`; THM-560 (Weyl on the far element),
  THM-527/529 (compact/equidistribution reductions).
- mac-mini HYP-4402 (two-scale decorrelation = quantitative equidistribution),
  HYP-4452 (Riesz density floor), HYP-4412 (three-gap); opus HYP-4396
  (sum-product).
