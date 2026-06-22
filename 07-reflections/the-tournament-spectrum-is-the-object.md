# The object is the spectrum, not the tournament

*kind-pasteur-2026-06-22-S41. A reflection prompted by the owner's correction: a single winding
tournament fixes one phase and is therefore magnitude-blind — so the right object is the *set* of
iso classes swept across all phases. It is, and it resolves the blindness.*

## The correction that mattered

I had spent a session proving that the LRC tight-locus census is *not* a tournament-iso-class
statement: the apex winding tournament `T(S, a*/14)` sees only residues mod 14, so the tight AP, the
loose `12→26`, and the loose `12→96` all collapse to the *same* regular rotational `R₁₃`. I concluded
"tournament iso class is necessary-only."

The owner's reply was one sentence: *"while any single tournament fixes one phase — so then the object
is a set of tournament iso classes."* That is the whole fix. A single tournament is a snapshot at one
phase `t`; the magnitude information lives in *how the snapshot changes as `t` moves*. The right object
is the **tournament spectrum**:

> `Σ(S) = { iso(T(S, t)) : t ∈ [0,1) }`, each class weighted by the measure of the `t`-interval that
> realizes it.

`T(S, t)` changes only at breakpoints `t = k/(sᵢ−sⱼ)` and `k/(2sᵢ)` — and *those denominators are the
magnitudes*. So the spectrum cannot be magnitude-blind. It isn't: `Σ(AP)` has 14 distinct classes
concentrated on `R₁₃` (measure 0.24), while `Σ(12→26)` has 24 classes with `R₁₃` down at 0.011 —
*different spectra for sets with identical apex tournaments*. The snapshot forgot the magnitude; the
movie remembers it.

## The useful invariant: the binding scale

A spectrum is a lot of data; the creative step is reading the one number out of it that decides
tightness. It is the **binding scale** — the denominator of the phase `t*` at which the spectrum's
*deepest sink* occurs (the optimum, where `M(S)` is achieved):

> **`S` is tight ⟺ the deepest-sink class of `Σ(S)` sits at binding scale 14** (the apex Farey node).

Verified exactly: AP and GW land at 14; every loose set lands elsewhere — `12→26 → 12`, `11→24 → 11`,
`13→26 → 27`, `12→36 → 41` (the Farey neighbor), `12→96 → 101`. This is THM-568 (a tight optimum has
denominator 14) re-read on the spectrum — but now it is a *complete* characterization, not a
necessary-only one, because the spectrum carries **both layers at once**: the deepest-sink *class*
(the residue layer — `R₁₃` or the GW-dipole) and the *scale* it sits at (the magnitude layer). The
single apex tournament had only the first; the spectrum has both.

## Why this is the natural geometry

The phase axis `t ∈ [0,1)` is a copy of the Farey / Stern–Brocot tree, and the spectrum is a labelling
of that tree by tournament iso classes. The deepest sink is a marked node. For a tight set the mark
sits at the **apex node `1/14`**; for a loose set it has **migrated** — either *up* the tree to a
coarser node (`q(S) < 14`, the divisibility threshold) or *down* to a Farey child (`3/41`, the
near-miss neighbor, det `[[1,3],[14,41]] = −1`). Tightness is the statement that *the deepest sink
does not migrate off the apex.* The two ways to lose — coarser (divisibility) and finer (Farey
neighbor) — are exactly the two regimes the divisibility-threshold `q(S)` and the Farey-rigidity
analyses found independently; the spectrum unifies them as "which node is marked."

## What it does and does not buy

It **buys** the conceptual correction: tightness *is* a property of an object made of tournament iso
classes — just the spectrum, not a single class. The owner's home-turf intuition was right; my "not a
tournament statement" was too hasty. It also gives a clean language: the census is "the set of `S`
whose spectral deepest-sink stays pinned to the apex node."

It does **not** dissolve the open core. Proving *which* sets keep the deepest sink at the apex — i.e.
that no Farey-child node ever out-sinks `1/14` except for `{AP, GW}` — is the same three-gap /
consec-maximizes rigidity, now phrased as a non-migration statement on the labelled Farey tree. The
spectrum relocates the difficulty into a more natural geometry; it does not remove it.

## The pointer beyond

If the spectrum is a tournament-labelled Farey tree, the right tool is the one that controls *how the
label changes along a tree edge*. Each breakpoint is a single mutual-pair flip (one antipodal crossing
`(sᵢ−sⱼ)t = 1/2`), so the spectrum is a **walk in the tournament flip-graph indexed by the Farey
tree**. The census conjecture becomes: *the only walks whose deepest-sink node is the apex are the
walk of the AP and its unique Jacobsthal-doubled child.* That is a statement about flip-graph geodesics
constrained by Farey geometry — the kind of object the tournament metagraph `G_n` (this project's
central study) is built to analyze. The Lonely Runner's tight locus may yet be a metagraph statement —
not about one tournament, but about a path through the space of them.

— Related: [[lrc14-thread]], HYP-2928 (the spectrum + binding scale), HYP-2927 (single-tournament
necessary-only), HYP-2917 (q=14), HYP-2919 (Farey-41), THM-568, HYP-2605 (the winding tournament);
`the-lonely-runner-is-a-random-round-tournament.md`, the metagraph `G_n` reflections.
