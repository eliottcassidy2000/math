---
id: HYP-2947
title: LRC14 local signatures need Paris-Harrington rank and J37 orbit labels
status: COMPUTATIONAL SCOUT / local-global guardrail; not a proof
source: codex-2026-06-24
tags: [lrc14, paris-harrington, johnson-solids, j37, local-global, extension-rank]
related:
  - HYP-2949
  - HYP-2948
  - HYP-2943
  - HYP-2247
  - HYP-2246
  - HYP-2245
results:
  - 04-computation/lrc14_ph_j37_pentagonal_minorant_codex.py
  - 05-knowledge/results/lrc14_ph_j37_pentagonal_minorant_codex.out
external:
  - https://arxiv.org/abs/2410.15880
  - https://en.wikipedia.org/wiki/Elongated_square_gyrobicupola
  - https://mathworld.wolfram.com/ElongatedSquareGyrobicupola.html
---

# HYP-2947: local signatures need rank/orbit labels

This note combines two local/global warning objects:

```text
Paris-Harrington: raw bad-coloring count is weaker than bad-child extension rank.
J37: local vertex figure is weaker than global vertex transitivity.
```

The computation is in
[lrc14_ph_j37_pentagonal_minorant_codex.py](/home/bigo/math/04-computation/lrc14_ph_j37_pentagonal_minorant_codex.py:1),
with output at
[lrc14_ph_j37_pentagonal_minorant_codex.out](/home/bigo/math/05-knowledge/results/lrc14_ph_j37_pentagonal_minorant_codex.out:1).

## Paris-Harrington Anchor

The script recomputes the HYP-2247 baby pair-coloring case:

```text
color pairs of [N] into 2 colors,
avoid a homogeneous 3-set H with |H| >= min(H).
```

Bad-coloring counts:

```text
N=1..6: 1, 2, 6, 18, 12, 0.
```

The derivative profile matters more than the count:

```text
N=4 extension_hist = {0: 6, 1: 12}
N=4 edge shell e=3 extends; shells e=2 and e=4 die.
N=5 all 12 surviving bad nodes die.
```

This is the finite toy version of the Paris-Harrington lesson:

```text
bad branch = side choice plus extension rank,
not just side choice or coloring density.
```

## J37 Anchor

The elongated square gyrobicupola `J37` is the solid-form version of the same
warning.  It has:

```text
8 triangular faces + 18 square faces
24 vertices, 48 edges
local vertex figure 3.4.4.4
```

These are the same local face counts and local vertex figure as the
rhombicuboctahedron.  But the 45-degree cupola twist splits the vertex set:

```text
8 polar vertices + 16 equatorial vertices.
```

So `J37` is locally Archimedean-looking but globally not vertex-transitive.

## LRC14 Use

The LRC14 analogue is:

```text
same local C27/Farey/product signature
does not imply same global proof state.
```

A candidate residual can share all easy local labels and still differ by:

```text
extension rank,
owner orbit,
C27 branch-local chart,
unital block reuse,
twist/orbit split.
```

So the next LRC14 proof carrier should attach a J37/PH-style pair:

```text
local signature + global rank/orbit label.
```

The proof target is:

```text
After AP/Goddyn-Wong and exact C27 labels are attached, any local row signature
that looks rhombicuboctahedral must also pass a J37 test: no hidden twist may
split the owner orbits without decreasing the Paris-Harrington bad-child rank.
```

This is not a theorem yet.  It is a guardrail against treating local product,
face, or block data as globally uniform.
