# Mahler `3/2` frontier

**Current status (2026-08-23): OPEN.**  This is the routed detail behind the
compact entry in [CURRENT-FRONTIER](CURRENT-FRONTIER.md).  Primary-source
records live in
[CORE-PAPERS-MAHLER-THREE-HALVES](../05-knowledge/reference/CORE-PAPERS-MAHLER-THREE-HALVES.md).

## Closest inherited mechanisms

- **PROVED:** THM-2228 separates strict real carry tails from ordinary
  positive-integer stabilization.  The periodic carry control `(01)^infinity`
  decodes to `Phi=-2/5` and never stabilizes.
- **PROVED:** THM-2352 realizes every finite carry cylinder in its unrestricted
  abscissa universe.  It does not preserve the Haar-null safe follower tree.
- **CITED:** Mahler gives countability and a height-counting bound;
  Dubickas--Mossinghoff force any candidate above `2^57`.
- **CONJECTURAL:** the 2026 rational-base normality preprint would exclude
  `Z_(p/q)` numbers for `1<q<p<q^2`; the normality premise is not proved.

## Exact finite prefix structure

[THM-3848](../01-canon/theorems/THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation.md)
proves, for every coprime `p>q>=2`, that the normalized length-`N` safe-prefix
set has a common wall grid, exactly `p^N` safe open atoms, Haar measure
`q^(-(N+1))`, and a reversible full `p`-ary atom refinement.  Pointwise there
is instead exactly one safe lift among the `q` inverse sheets.  These are two
different maps; identifying their branching factors destroys the wall address.

The adjacent relations `p e_n-q e_(n+1)` generate the full integer kernel.
The mixed-power speed row has exact maximum

```text
M(q^N,q^(N-1)p,...,p^N)=floor((p+q)/2)/(p+q),  N>=1.
```

For `p/q=3/2`, the maximum is `2/5` and the only maximizing times are the two
mod-five phases.  At `N=12` this gives a primitive thirteen-speed, hence
fourteen-runner, positive control with twelve norm-five relations.  It misses
every multiple of five and is sieve-trivial, so abundant short relations and
a perfect recursive address do not characterize the LRC(14) hard core.

## Where LRC and Mahler separate

The LRC evaluator retains centered distance `||x||` and identifies the phases
`2/5` and `3/5`.  Mahler's predicate is oriented: it accepts only the lower
half-circle.  Along the mod-five extremizer the phases alternate, so the exact
LRC maximizer fails the Mahler test.  This phase side is the first destroyed
coordinate in the map from an oriented orbit to centered loneliness.

The distinct closed formal safe-tail shift `K` has greedy equality boundary
`d`, renewal law

```text
a_m=1+sum_(k=1)^m d_k a_(m-k),
A(z)=1/((1-z)(1-D(z))),
```

entropy `log(3/2)`, binary-ultrametric Hausdorff dimension `log_2(3/2)`, and
is nonsofic.  The strict set removes the countable backward orbit of `d` and
has the same finite language.  This does **not** prove that the unknown
Z-language is nonsofic: a Mahler candidate must also satisfy THM-2228's
ordinary positive stabilization condition.

## ABC/IUT sidecar

Assuming ABC, actual odd carry-one packets have asymptotic radical exponent at
least one.  The denominator-19 safe-prefix family at positive multiples of 18 is
also logarithmically radical-saturated.  Neither fact supplies a horizon bound:
arbitrarily long even zero-carry runs have adjacent radical `6`, and the
denominator-19 construction changes its initial integer with `m`.  IUT adds no
unconditional input here; only the contested claimed IUT-to-ABC implication
could feed the explicitly conditional ABC consumer.

## Live finite-exact task

Build the fibre product of the strict safe follower tree with THM-2352's
plateau/stabilization clock.  Preserve four sidecars throughout:

1. oriented phase rather than centered distance;
2. equality-boundary versus strict status;
3. native carry residue and suffix follower state;
4. ordinary terminal/stabilization state.

Mandatory hostile controls are the greedy equality boundary,
`(100)^infinity` (safe but nonterminal), the `A=1` itinerary (terminal but
unsafe), `(01)^infinity` (the nonterminal mod-five maximizer), and the
changing-start denominator-19 tower.  A successful finite computation must
state why its filtered universe still contains every genuine candidate; raw
survival counts alone are not a Mahler or LRC theorem.
