---
source: codex-2026-06-01-S455
status: integration note
tags:
  - lonely-runner
  - n14
  - n18
  - tournament-analysis
  - first-doubling
  - bridge-fibers
---

# LRC n=14/n=18 Alternating Noise Session

The user asked for a back-and-forth session: work on `n=14` and `n=18`, get
stuck, search around for forced-random inspiration, then bring that noise back
into computations.

The useful external noise was:

- Tao's LRC remarks: finite checking is possible after bounding velocities,
  so local finite certificates are legitimate proof objects.
- Malikiosis-Santos-Schymura: the zonotopal reinterpretation and improved
  finite-checking bounds say covering geometry is not decoration; it is one of
  the main modern proof languages.
- Jensen's mixed-threshold variant: unequal thresholds suggest separating
  unit, half-gate, and bridge rows rather than flattening all endpoint debt.
- Alcantara-Criado-Santos: shifted LRC via covering radii and dyadic
  fundamental domains supports the repo's habit of tracking endpoint cells and
  2-adic debt together.

S455 metabolizes those as exact repo computations, not imported proof claims.

## Alternating Result

The first pass rechecked the local `n`-gate invoice.

For `n=14`, the old picture survives exactly:

```text
forced = (1,3,5,7,9,11,13)
free bridge = one of 2,4,6,8,10,12
```

For `n=18`, the analogous row is unexpectedly tighter:

```text
forced = (1,5,7,9,11,13,17)
free bridge = one of 6,12
```

That is the main new fact from the session.  `n=18` is not "n=14 but larger."
It is the square-core version of the first-even seam, and the square core
collapses the bridge fiber from six cases to two.

## Debt Ledger

The row-parent and n-gate ladders say the bridge trade is conserved:

```text
n=14 row-parent: gap/th=5/924, debt=84,  product=5/11
n=14 n-gate:     gap/th=5/1848, debt=168, product=5/11

n=18 row-parent: gap/th=1/176, debt=176, product=1
n=18 n-gate:     gap/th=1/352, debt=352, product=1
```

So the gate ladder does not solve the obstruction.  It moves the same
gap-debt product to a deeper owner layer.  This sharpens the user's slogan:
gap is archimedean size, endpoint debt is adic size, and the product is the
thing to try to conserve.

The one-slot repair scan also clarified the earlier `n=14` thought that
tightness wants no multiple of `14` but a counterexample wants one.  In the
small family "drop `n-1`, add one speed," the best nonmultiple can have a
smaller gap than the best multiple:

```text
n=14: add 14 has gap/th=1/22, debt=24; add 54 has gap/th=1/30, debt=34
n=18: add 18 has gap/th=1/30, debt=32; add 70 has gap/th=1/42, debt=42
```

This does not refute the multiple heuristic.  It says the heuristic needs a
debt clause: nonmultiples can tighten the archimedean gap, but they carry more
unprotected endpoint debt.

## Tournament Analysis Reading

At `t=1/n`, the safe-distance switchboard distinguishes initial skeletons from
ladder repairs:

```text
n=14 initial safe cycles = 0
n=14 row-parent ladder safe cycles = 93
n=18 initial safe cycles = 0
n=18 row-parent ladder safe cycles = 207
```

The scalar lonely gap sees a one-dimensional obstruction.  The switchboard sees
pairwise crowding.  The current proof direction should keep both: a marked
stationary bracket for the endpoint proof, and a pairwise safe-switch
tournament for the crowding core.

## Favorite Next Target

For the next proof sprint, I would attack `n=18` first.  The reason is not
that it is easier numerically; it is structurally cleaner:

```text
local bridge cases: 2 instead of 6
row-parent product: exactly 1
forced half-gate: 9
bridges: 6 and 12
```

If a two-case Farkas/Hall dual can be built for `n=18`, the `n=14` proof may
become the same dual summed over a six-choice parity bridge family.

## Files

- `04-computation/lrc_n14_n18_alternating_noise_s455.py`
- `05-knowledge/results/lrc_n14_n18_alternating_noise_s455.out`
- `05-knowledge/hypotheses/HYP-1942-lrc-first-even-bridge-fiber.md`
