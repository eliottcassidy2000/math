# LRC invariant stack

Codex 2026-06-19.

The prompt was to look for more Lonely Runner invariants and ask what actually
determines the structure.  I built an exact atlas over 640 primitive 13-speed
rows, mixing the LRC14 tight/near-tight families, S3 cluster shapes, covering
biased random rows, and generic primitive rows.

Namespace/caveat: after this run started, origin/main added HYP-2651 for the
core-gap atlas and opened CASE-thm538, disputing the full zero-padded `K`
support-six floor.  This note is therefore HYP-2652/T899.  Its
`support6plus` statistic is an active equal-subset/additive-relation proxy, not
a claim that short zero-padded Fourier relations vanish.

Files:

- `04-computation/lrc_invariant_atlas_codex_20260619.py`
- `05-knowledge/results/lrc_invariant_atlas_codex_20260619.out`

The test was a fiber test.  If two rows share an invariant signature but have
different exact `M(S)` or different safe-measure at level `1/14`, that invariant
does not determine the structure.  I also tracked group counts so overly fine
signatures do not receive fake credit for making nearly every row a singleton.

## The Stack

No single scalar invariant determines the LRC structure.  The data points to a
stack:

```text
CRT obligation
  -> additive anti-coset / relation density
  -> safe-component boundary-owner geometry
  -> binding denominator readout
```

Each layer resolves ambiguity left by the layer below it.

Residues and q-coverage are necessary gates, but not structure.  They filter the
search space and explain why a row is even eligible to be hard.  They do not
explain which endpoint survives.

Additive invariants are much closer.  The best raw correlations with exact
hardness `M(S)` were:

```text
corr(M, three_term_count)       = -0.919
corr(M, longest_consecutive_run)= -0.896
corr(M, difference_energy)      = -0.890
corr(M, four_sum_collisions)    = -0.889
corr(M, support6plus_subsets)   = -0.722
```

That is a strong signal: near-tight rows are relation-dense.  This agrees with
the anti-coset/relation-density direction from HYP-2612, the corrected coset
quotient HYP-2646, and the older relation-lattice covolume picture from
HYP-2606, subject to CASE-thm538's warning about the full kernel.

Safe-component geometry is the first layer that feels almost structural.  It
records which speeds actually form the surviving witness reservoir at level
`1/14`: number of components, boundary owner pairs, lengths, and owner-pair
weights.  Topology alone helps but is often too specific or still ambiguous;
exact boundary-owner geometry is close to determinative because it is almost the
object itself.

Binding denominators are the cleanest readout, not the cleanest input.  They
nearly determine `M`, but only after solving enough of the problem to know the
active pair.

## Tail Laws

The exact atlas highlighted two near-tight one-tail laws:

```text
S = {1,2,...,12,13a}       has M(S) = a/(13a+1)   for tested a=2..9.
S = {1,2,...,11,13,12a}    has M(S) = a/(12a+5)   for tested a>=3.
```

The second family contains the Goddyn-Wong point at `a=2`, where the formula
would not apply but the row is exactly tight:

```text
{1,2,...,11,13,24}: M=1/14.
```

These laws are binding-denominator invariants.  They are not merely "large tail"
phenomena; the denominator is a repaired endpoint created by the small AP core.
For `{1..12,13a}`, the binding denominator is `13a+1`.  For
`{1..11,13,12a}`, the binding denominator is `12a+5` once `a>=3`.

## Tournament Lesson

I tested a tournament on safe-component boundary owner pairs.  It is useful as a
compression, but not as the whole invariant: the owner-tournament fingerprint
had large mixed fibers.  The LRC structure is not determined by the tournament
shadow alone; it needs the boundary lengths and exact endpoint positions.

This reinforces the repo's repeated lesson: tournaments are best used as a
quotient of a richer object.  For this subproblem the richer object is the
endpoint/safe-component decomposition.

## What I Would Keep

The next proof-facing invariants should be:

1. Additive anti-coset ledger: counts of low-height equal-subset identities,
   especially support at least six.
2. Endpoint owner geometry: for each safe component, the owner speeds at its two
   boundaries and the component length.
3. Binding-denominator remainder: `D=14j-r`, but treated as an output to be
   forced by the first two layers.
4. CRT obligation profile only as a gate, not as a determinant.

This is the structural answer I see: LRC is determined by endpoint geometry
conditioned by arithmetic gates, with additive anti-cosets explaining why the
endpoint geometry becomes sparse near the tight rows.
