# LRC14 Pi, Flower Turns, and Unital Ledgers

Codex session note, 2026-06-24.

This note merges the previous Farey/spectrum work with three new prompts:

```text
22/7 approximates pi,
31^(1/3) approximates pi,
turning by 1/pi radians gives a 22-family flower effect,
and "unital" has design-theoretic and algebraic meanings.
```

The computation is in [lrc14_pi_unital_flower_codex.py](/home/bigo/math/04-computation/lrc14_pi_unital_flower_codex.py:1), with output at [lrc14_pi_unital_flower_codex.out](/home/bigo/math/05-knowledge/results/lrc14_pi_unital_flower_codex.out:1).

## Pi Approximants

The two approximants have different mathematical characters.

```text
22/7       = 3.142857142857143, error +1.264489e-3
31^(1/3)   = 3.141380652391393, error -2.120012e-4
```

So `31^(1/3)` is about `5.965x` closer to pi than `22/7`.  It comes from:

```text
pi^3 = 31.006276680299816.
```

Thus 31 is the nearest integer to `pi^3`.  In the LRC14 ledger:

```text
31 = 2*14 + 3.
```

That is only a checksum, not a proof hook by itself.  But it is compatible with the cube-root/Eisenstein mode from the earlier `A+B+C-D-E-F+G` sessions: cubic structure naturally asks for a three-phase carrier, and the residue of 31 mod 14 is 3, the same numerator threshold where the complete-bipartite Farey packets first become nonplanar.

The rational approximant `22/7` is a continued-fraction/Farey object.  Its reciprocal packet is:

```text
7/22 -> a+b=29, ab=154=11*14, K_{7,22}.
```

So it carries the apex numerator 7 and the flower denominator 22, but it is already far past the `K_{3,3}` threshold in the complete-bipartite reading.

## The Flower Turn

Let:

```text
theta = 1/pi radians.
```

Then:

```text
theta = 0.318309886183791
7/22 = 0.318181818181818
theta - 7/22 = 1.280680e-4
```

So `22` really appears as the denominator of the continued-fraction convergent:

```text
1/pi = [0;3,7,15,...] = 7/22 + small error.
```

After 22 petals:

```text
22/pi = 7.002817496043395 radians.
```

That is extremely close to 7 radians.  This explains a 22-family effect if the construction is being read in radian-step coordinates.

But full-circle closure uses the normalized turn:

```text
delta = theta/(2*pi) = 1/(2*pi^2).
```

Its early convergents are:

```text
1/19, 1/20, 3/59, 4/79, 23/454, ...
```

So if "families" means circular return/parastichy count, the natural small denominators start at 19, 20, 59, and 79, not 22.  The useful conclusion is:

```text
22-family = radian-denominator phenomenon from 1/pi ~= 7/22
20-family = first strong full-circle closure phenomenon from 1/(2*pi^2) ~= 1/20
```

If a rendered flower visibly shows 22 families, the likely cause is that the visual grouping is sensitive to the radian-step convergent rather than to full-circle closure.  That is worth testing with an image/neighbor-clustering script later, but the arithmetic distinction is clear.

## q=3 Unital: Stronger Than Numerology

The q=3 geometric unital has parameters:

```text
v = q^3 + 1 = 28 points
k = q + 1 = 4 points per block
lambda = 1
b = 63 blocks
r = 9 blocks through each point
```

This lands on the LRC14/Farey ledger almost too well:

```text
28 = sum of distinct products in F_4 = 2*14
63 = number of proper reduced Farey packets 0<a<b<=14
4  = first Farey level with 3/4 -> K_{3,4}
9  = |E(K_{3,3})|
3  = numerator threshold for nonplanar K_{a,b}
```

The previous note found:

```text
F_3 product set = {1,2,3,6}, sum 12, all planar
F_4 product set = {1,2,3,4,6,12}, sum 28
3/4 -> (a+b,ab,K) = (7,12,K_{3,4})
```

The unital adds a pair-uniform design frame exactly at the `F_4` crossing:

```text
2-(28,4,1).
```

Because every pair of points lies in exactly one block, the q=3 unital is a natural averaging frame for pair-slot residuals.  This agrees with the older S105 note: the unital is not an exact tiler for LRC14, but it may be a Hodge/averaging carrier.

Concrete proof target:

```text
Build a ledger from the 63 proper F_14 packets to the 63 unital blocks,
then test whether pair-residual averages over the 28 unital points isolate
the AP/GW tight packets or force nonpositive residual leakage.
```

That would turn the unital from a parameter coincidence into a finite certificate.

## Algebraic "Unital" as Mark Preservation

The algebraic meaning of unital is different: a map preserves identity,

```text
f(1_A) = 1_B.
```

For LRC14, the correct analogue is:

```text
a proof quotient is unital if it preserves the marked apex node 1/14.
```

This is exactly where bare tournament iso classes failed.  The apex winding tournament can identify tight and loose rows because it preserves a residue-level structure but forgets the marked spectrum node.  In the marked-spectrum category, that quotient is non-unital: it does not preserve the identity object "the deepest sink is pinned to 1/14."

So the word "unital" gives a useful proof discipline:

```text
Do not accept a quotient unless it is mark-preserving.
```

In practice:

```text
AP        -> marked at 1/14
GW        -> marked at 1/14
12->26    -> marked at 1/12
13->26    -> marked at 2/27
12->36    -> marked at 3/41
```

Any relation that identifies these without retaining the mark is non-unital for the LRC proof.

There is also a "stably unital" analogy: after adding a far runner, taking a quotient, or stabilizing by a frame, does the apex unit still survive?  The loose rows say this is not automatic.  A stable-unital lemma would say exactly which operations preserve the marked `1/14` sink.

## How This Merges With The Previous Topics

The current proof-search dictionary is now:

```text
Farey tree:
  marked node 1/14 vs migrated nodes 1/12, 2/27, 3/41, ...

Summand graph:
  a+b is the pinch denominator; 3/4 has a+b=7.

Multiplicand graph:
  ab is both product node and |E(K_{a,b})|; 3/4 has ab=12.

Complete bipartite graph:
  K_{3,4} is the first reduced Farey carrier containing K_{3,3}.

Tournament metagraph:
  Sigma(S) is a flip-walk; the LRC datum is the marked sink on that walk.

q=3 unital:
  2-(28,4,1) is the pair-uniform frame sitting exactly at the F_4/LRC14 crossing.

Pi/flower:
  22/7 and 7/22 explain the radian 22-family effect; circular closure uses
  1/(2*pi^2) and gives a different parastichy ledger.

Cube root 31:
  pi^3 ~= 31 gives a cubic shell; 31=2*14+3 keeps the three-phase/cube-root
  filter in the same checksum family.
```

The strongest new item is the q=3 unital alignment.  The flower/pi observations are useful coordinates and sanity checks, but the unital gives an actual finite design frame with the right cardinalities to be tested.

