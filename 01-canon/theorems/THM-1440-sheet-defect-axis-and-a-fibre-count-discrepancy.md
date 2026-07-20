---
id: THM-1440
title: "THE SHEET-LOSS STATEMENT IS RUNG 2 OF FIVE, AND ITS OWN AXIS IS NOT WHERE THE ROOM IS. (0) PLACEMENT: 'a constant-Jacobian polynomial local biholomorphism with no loss of sheets at infinity — e.g. a proper Keller map — is an automorphism' is exactly the classical properness rung of THM-1375's lattice, the second-weakest of five. Ours is considerably more comprehensive: three strictly stronger rungs sit above it (nodal-Jelonek, Galois/Campbell, and Smith's non-self-normalising stabiliser as the maximal element), with the ordering PROVED strict by the D₄ witness. (I) THE QUANTITATIVE WEAKENING (Q_k) — 'Keller + at most k sheets lost ⟹ automorphism' — is FALSE at k = 2, witnessed by the THM-1300 map itself: exact elimination finds targets q = (0,−1,−1) and (0,0,−5/3) whose fibre is a single point, defect 2, with the resultant degree collapsing 3 → 1. (II) A DISCREPANCY I DID NOT RESOLVE AND AM NOT ASSERTING AWAY: THM-1315 states the drop is 3 → 2 'by a sheet leaving through infinity', defect 1. My probes found only defect-2 points and never a defect-1 point. Both can coexist as different strata of the non-properness locus, but (Q_1) is NOT established here and must not be cited as settled. (III) AN INSTRUMENT FAILURE, recorded because its signature is the useful part: my first fibre counter used Groebner standard monomials with an unpinned monomial order and a truncated counting box. It returned the correct 3 on the known triple-collision control while returning the box cap 40 on every generic target — a control passing while generic input fails is the signature of a broken instrument, not of a phenomenon"
status: >
  (0) Placement is exact; the lattice and its strictness are THM-1375, already proved.
  (I) VERIFIED-EXACT by elimination (solve F₃ for z, resultant in y, count distinct
  roots, x = 0 branch handled separately), with the instrument validated on BOTH the
  known control and five generic targets before any drop was reported — all six give
  fibre 3.
  (II) OPEN and flagged.  My sampling found defect 2 only; THM-1315 asserts defect 1
  generically on {K = 0}.  I did not reproduce defect 1 and do not dispute it — random
  rational probes on a hyperplane need not meet a codimension-one stratum, so this is
  inconclusive rather than contradictory.  The minimal non-zero defect is therefore NOT
  pinned, and (Q_1) is undecided by this work.
  (III) Recorded as a method note; the discarded run is not used anywhere.
  Nothing here advances JC₂, DC₁ or DC₂.
source: kind-pasteur-2026-07-20-S128c113 (owner: is the no-loss-of-sheets statement roughly equivalent to our understanding, or is ours more comprehensive?)
depends_on:
  - THM-1375    # the five-rung lattice this places the statement in
  - THM-1315    # the 3 -> 2 drop this partially disagrees with
related: [THM-1330, THM-1430]
script: 04-computation/sheet_defect_sharpness_kps_S128c113.py (+ .out)
---

# THM-1440 — placing the sheet-loss statement, and one number I could not confirm

## 0. The answer to the comparison

> *"A constant-Jacobian polynomial local biholomorphism with no loss of sheets at
> infinity — e.g. a proper Keller map — is an automorphism."*

This is **true, classical, and rung 2 of five** in THM-1375's lattice:

```
injective  ⟹  PROPER  ⟹  nodal Jelonek  ⟹  Galois  ⟹  non-self-normalising stabiliser
(Ax–Groth.)   ^^^^^^^^    (Zariski–Lefschetz   (Campbell    (Smith selection rule,
              the quoted    + Deligne–Fulton)   1973)        the MAXIMAL element)
              statement
```

Each arrow *weakens* the hypothesis, so each rung is a strictly stronger theorem, and
THM-1375 verified the ordering is **strict** — `D₄` at cover degree 4 is killed by Smith
and not by Campbell. "No loss of sheets at infinity" is precisely "the Jelonek set is
empty", and THM-1330 records that a counterexample's Jelonek set is *always* non-empty, so
the quoted statement is the trivial direction of that fact.

**So our picture is considerably more comprehensive**, and not only by having more rungs.
We also have, for the specific counterexample: cover degree 3 with `S₃` monodromy; the
proof that degree 2 is impossible; that `d = 3` is the *unique* degree at which Galois-ness
is a single quadratic character; the Jelonek set identified as Zariski's three-cuspidal
quartic; and an explicit symmetric-Jacobian counterexample on `ℂ⁶` (THM-1430).

## I. Can rung 2 be improved along its own axis?

"Proper" means **zero** sheets lost. The natural weakening is quantitative:

> **(Q_k)** Keller + at most `k` sheets lost at infinity ⟹ automorphism.

`(Q_0)` is the quoted statement. **`(Q_2)` is false**, witnessed by the THM-1300 map.
Exact elimination — solve `F₃` for `z`, take the resultant in `y`, count distinct roots,
handle the `x = 0` branch separately — gives:

| target `q` | fibre | resultant degree |
|---|---|---|
| `(−1/4, 0, 0)` (known triple collision) | **3** | 2 (+ the `x=0` branch) |
| five generic targets | **3** each | 3 |
| `(0, −1, −1)` | **1** | **1** |
| `(0, 0, −5/3)` | **1** | **1** |

Two sheets leave at once: the resultant degree collapses `3 → 1`.

## II. The number I could not confirm

THM-1315 §3 states the drop is **`3 → 2`**, "by a sheet leaving through infinity, never by
two finite points colliding" — defect **1**. **I never found a defect-1 point.** Every drop
my probes located was `3 → 1`.

I am not asserting this away, and it is not a contradiction: the non-properness locus is a
codimension-one subvariety, random rational points on a target hyperplane need not meet its
generic stratum, and a defect-2 stratum can perfectly well sit inside a defect-1 locus.
Suggestively, THM-1315's caustic is `Δ ∝ a²·K` with `a` appearing **squared**, and both of
my drop points have `a = 0` — a double factor is what one would expect over a defect-2
stratum, with `{K = 0}` carrying the defect-1 generic behaviour.

**Consequence for the sharpness question: the minimal non-zero defect is NOT pinned.**
`(Q_1)` is undecided by this work. If THM-1315's `3 → 2` is right then `(Q_1)` is false too
and rung 2 is sharp at the first possible threshold; that should be settled by evaluating
`K` explicitly rather than by sampling.

## III. The instrument that lied, and how it announced itself

My first fibre counter read leading monomials off `Poly.monoms()[0]` without pinning the
monomial order, then counted standard monomials inside a truncated box. It returned the
**correct** answer 3 on the known triple-collision control — and the box cap, 40, on every
generic target.

> **A control that passes while generic input fails is the signature of a broken
> instrument, not of an interesting phenomenon.**

The corrected method needs no monomial-order bookkeeping at all, and was validated on the
control *and* five generic targets before any drop was reported. The discarded run is used
nowhere. This is the same discipline as MISTAKE-196 applied to my own tooling rather than
to a search: an instrument earns its negatives only after it has produced the positives it
should.

## Named next

- **Settle `(Q_1)`** by extracting `K` from THM-1315 and evaluating the fibre on `{K = 0}`
  directly. That resolves the discrepancy in §II and pins the sharpness threshold.
- If the defect-1 stratum exists, the interesting object is the **stratification** of the
  non-properness locus by defect — `{K=0}` at defect 1, a deeper locus at defect 2 — which
  is a finer invariant of the counterexample than the Jelonek set alone.
