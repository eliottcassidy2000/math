# LRC14 Spectrum as a Marked Farey Walk in the Tournament Metagraph

Codex session note, 2026-06-23.

The working object is no longer one tournament at one phase.  It is the spectrum:

```text
Sigma(S) = { iso(T(S,t)) : t in [0,1) }.
```

Here `T(S,t)` is the winding tournament on the runner speeds.  As `t` moves, the tournament is piecewise constant.  At a breakpoint, one or more runner pairs lie on a collision or antipodal wall, and crossing the wall flips exactly those arcs.  Thus `Sigma(S)` is a walk in the tournament flip graph / metagraph `G_n`.

The LRC datum is an extra mark on this walk:

```text
marked node = M(S) = max_t min_s ||s t||.
```

So the proof object should be:

```text
(Farey/Stern-Brocot node p/q, tournament iso class at p/q, flip-walk neighborhood)
```

rather than only a tournament iso class.

## Tested Ledger

The script [lrc14_spectrum_farey_flipgraph_codex.py](/home/bigo/math/04-computation/lrc14_spectrum_farey_flipgraph_codex.py:1) computes exact `M(S)`, the marked Farey node, and the phase-walk breakpoints for several LRC14 rows.  A result snapshot is stored at [lrc14_spectrum_farey_flipgraph_codex.out](/home/bigo/math/05-knowledge/results/lrc14_spectrum_farey_flipgraph_codex.out:1).

The key rows:

```text
AP        -> q(S)=14, M=1/14, tight apex
GW 12->24 -> q(S)=14, M=1/14, tight apex
11->24   -> q(S)=11, M=1/11, loose-up divisibility node
12->26   -> q(S)=12, M=1/12, loose-up divisibility node
13->26   -> q(S)=14, M=2/27, loose-down Farey child, det=-1
12->36   -> q(S)=14, M=3/41, loose-down Farey child, det=-1
12->96   -> q(S)=14, M=8/101, loose-down nonunit excess
```

This is exactly the owner's synthesis:

```text
tight = marked at 1/14
loose-up = migrated to q(S)<14
loose-down = migrated to a Farey child below the apex in the Stern-Brocot refinement
```

The two older failure modes are therefore one statement:

```text
which Farey node is marked?
```

## The Child Ladder and the Bipartite Obstruction

The unit-excess children above `1/14` are:

```text
p/q = p/(14p-1),  q - 14p = -1.
```

The first few are:

```text
1/13  -> K_{1,13}, planar
2/27  -> K_{2,27}, planar
3/41  -> K_{3,41}, contains K_{3,3}
4/55  -> K_{4,55}, contains K_{3,3}
```

This folds the previous Farey-bipartite observation into the spectrum story.  The first reduced bipartite obstruction in ordinary Farey level was:

```text
3/4 -> K_{3,4}, with (a+b,ab) = (7,12).
```

The first nonplanar unit child above the LRC14 apex is:

```text
3/41 -> K_{3,41}.
```

So the near-miss `12->36` is not just a numerical det `-1` accident.  It is the first child on the LRC14 Stern-Brocot ladder whose Farey packet has numerator 3, hence the first child carrying a `K_{3,3}` minor in the complete-bipartite interpretation.

The child `2/27` is still planar as `K_{2,27}`.  That matches the `13->26` row: it has `M=2/27`, a unit Farey child, but it has not crossed the bipartite Kuratowski threshold.  The `12->36` row crosses both:

```text
Farey unit child: q - 14p = -1
bipartite obstruction: numerator p = 3
```

## Flip-Walk Verification

The walk data supports the metagraph interpretation directly.  For every tested row, the histogram of actual arc changes across consecutive phase intervals matched the event multiplicity histogram.  That is:

```text
number of arcs flipped at breakpoint = number of runner pairs on that wall.
```

Generically the breakpoint is a single arc flip.  Symmetric phases produce multi-flips.  The large `78` flip at the global symmetry wall is the all-pair reversal event on 13 vertices.

Examples:

```text
AP:
  intervals = 92
  directed iso classes = 23
  flip histogram includes {1:8, ..., 78:1}

GW 12->24:
  intervals = 344
  directed iso classes = 84
  single-flip events = 260

12->96:
  intervals = 2096
  directed iso classes = 149
  single-flip events = 2012
```

The exact counts here use full directed tournament isomorphism on `[0,1)`.  Earlier notes reporting `|Sigma(AP)|=14` appear to be using a folded or otherwise different quotient convention.  The migration result does not depend on that convention.

## Endpoint Classes Matter

One useful warning from the computation: the marked optimum class can be a measure-zero endpoint class.

For `12->26`:

```text
q(S)=12
M=1/12
opt-class-measure=0
```

That means the open-interval weighted spectrum and the marked certificate are related but not identical.  A proof should not say "the largest-measure class decides tightness."  It should say:

```text
the marked Farey node decides tightness,
and the flip-walk supplies the local metagraph neighborhood of that mark.
```

This matches the geometric language: `Sigma(S)` is the path; `M(S)` is the deepest sink marked on the path.

## Proof Target

The emerging LRC14 statement can be phrased as a marked-walk non-migration lemma.

For a primitive covering 13-set `S`, construct:

```text
W(S) = phase walk t -> iso(T(S,t)) in G_13
m(S) = marked Farey node M(S)
q(S) = first uncovered divisibility denominator
```

Then the desired classification is:

```text
S tight  <=>  m(S)=1/14
```

with two certified escape modes:

```text
q(S)<14                 -> loose-up, divisibility migration
q(S)=14 and M(S)>1/14   -> loose-down, Farey/Stern-Brocot migration
```

The sharper local target is:

```text
Every non-AP/GW survivor with q(S)=14 must mark a Farey child or nonunit descendant
outside the tight apex classes.
```

The old AP/GW census, the Farey mediant arithmetic, the summand/multiplicand packet picture, and the tournament metagraph now say the same thing in four projections:

```text
summand graph:          apex denominator 14, pair-sum node 7
multiplicand graph:     seed/product thresholds 12 and 14
Farey tree:             marked node 1/14 vs migrated child
tournament metagraph:   arc-flip walk with a marked sink
```

This still does not prove LRC14.  It gives a cleaner target: prove that the marked walk of any covering residual cannot migrate while remaining in the tight residue shell, except for AP and GW.

