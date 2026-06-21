---
id: HYP-2705
title: LRC14 apex gas has a projective synchronizer route for the true-wide deviation
status: OPEN; exact rational scout and proof fragments stored
source: codex-2026-06-20-S66
tangent: T940
depends_on:
  - HYP-2701
  - HYP-2702
  - HYP-2704
  - HYP-2684
  - HYP-2675
  - THM-557
  - OPEN-Q-108
related:
  - HYP-2700
  - HYP-2698
  - HYP-2681
  - HYP-2627
  - HYP-2637
  - THM-554
  - THM-556
  - THM-558
---

# HYP-2705 - Apex Gibbs-Propagator Synchronizer Route

## Namespace Note

The prose index currently contains a live `HYP-2704` collision between the
mac-mini sector-route summary and the file
`HYP-2704-lrc14-survival-seven-packet-quotient.md`.  This file intentionally
claims the next unused file namespace, `HYP-2705`, for the present synthesis.

## Claim

The new external prompts - Gibbs measure, Arnold cat map, Fubini-Study metric,
road coloring, Hopfield/Hebbian learning, propagator weights, Clifford plus
`T`, and the complete-graph crossing carrier - do not give independent proof
routes by themselves.  They collapse into one usable architecture:

```text
apex-prime phase gas on Z/7
  free/stabilizer layer:    decorrelated cut-space occupancy
  projective phase layer:   signed relation/cycle deviation
  synchronizer layer:       residual-sector automaton collapse
```

In this architecture the true-wide LRC14 crux is:

```text
actual survival currency
  = exact death-chain / Gibbs equilibrium quotient
    + projective signed resonant deviation.
```

The death-chain quotient is already cap-friendly after two decorrelated far
hits.  The open work is to bound the signed deviation by a Fubini-Study angle,
relation-lattice magic rank, or finite synchronizer atlas.

## Exact Scout

Script:

```text
04-computation/lrc14_fs_synchronizer_gibbs_codex_s66.py
```

Output:

```text
05-knowledge/results/lrc14_fs_synchronizer_gibbs_codex_s66.out
```

The scout keeps exact `Fraction` arithmetic.  For the missed-count death-chain
transition law, define

```text
|psi_{t,r}> = sum_s sqrt(Pr[t -> s after r iid color hits]) |s>.
```

Then the squared Fubini-Study overlap with the middle-survival subspace
`M=span{|1>,|2>,|3>,|4>}` is the exact rational

```text
cos^2(theta_M) = ||Proj_M psi_{t,r}||^2
               = sum_{s=1}^4 Pr[t -> s after r hits].
```

For the fully missed state `t=6`:

```text
r   cos^2(theta_M)      Pr(sync to 0)     E[C after r hits]
0   0                   0                 -4
1   0                   0                 -4/7
2   30/49               0                  26/49
3   300/343             0                  296/343
4   330/343             0                  2306/2401
5   16620/16807         0                  16616/16807
6   16650/16807         720/117649        116546/117649
```

The Gibbs obstruction ledger is even sharper:

```text
r=0: negative depths = t=6:-4
r=1: negative depths = t=6:-4/7
r=2: negative depths = none
r=3: negative depths = none
```

Thus the missed-count quotient itself stops being negatively charged after two
decorrelated far hits.  Every hard true-wide row must spend that positive
boundary currency through signed resonant deviation away from the quotient.

## Proof Fragments

### 1. Projective projection identity

Let `psi=sum_s a_s |s>` be a unit vector in a Hilbert space with the depth basis
orthonormal, and let `M` be the span of a subset of basis states.  The
Fubini-Study distance from `[psi]` to the projective subspace `[M]` satisfies

```text
cos^2 d_FS([psi],[M]) = ||Proj_M psi||^2 = sum_{s in M} |a_s|^2.
```

For the death-chain amplitudes above, this is exactly the rational middle mass.
No numeric angle is needed.

### 2. Fubini-Study deviation bound

For normalized phase profiles `psi,phi` and any linear observable `L`, after
choosing the global phase of `phi` so that `<psi,phi>` is nonnegative real,

```text
|L(psi)-L(phi)| <= 2 ||L|| sin(d_FS(psi,phi)/2).
```

Proof: `||psi-phi||^2 = 2-2 cos d_FS = 4 sin^2(d_FS/2)`, then apply
Cauchy-Schwarz.  The LRC content is not this inequality; it is the missing
number-theoretic lemma bounding the Fubini-Study angle between the true
resonant phase profile and the decorrelated/stabilizer profile.

### 3. Road-coloring deletion synchronizer

In the six-inner-sector deletion automaton, a word `w` synchronizes a known
missed set `R` to depth `0` if and only if every color in `R` occurs in `w`.
The shortest synchronizing word length for known `|R|=t` is `t`.  For an
unknown nonempty missed subset of the six inner colors, length `6` is necessary
and sufficient.

This is the road-coloring theorem stripped to the LRC residual automaton: the
automaton is already synchronizing, but LRC only needs partial synchronization
into the middle subspace before the cap comparison.  That partial target is why
two hits can be useful even though full synchronization from `t=6` is still
impossible in two or three hits.

### 4. Low-temperature Gibbs reading

Use energy

```text
E_r(t) = max(0, -E[C after r hits | depth t]).
```

The low-temperature Gibbs measure on depths concentrates on `t=6` for `r=0`
and `r=1`, but has no positive-energy obstruction after `r=2`.  This recovers
the HYP-2701 two-far signal in statistical-mechanics language: the finite cap
problem is not a problem of the depth quotient's free energy, but of relation
phase errors that keep the true row from behaving like the free quotient.

## Dictionary For The User's Prompts

* **Gibbs measure:** use a Gibbs weight on proof states or residual depths.
  Low temperature selects the obstruction.  Here that obstruction disappears
  in the depth quotient after two decorrelated hits, so the real energy must
  live in signed relation phase.

* **Arnold cat map:** not a literal replacement for LRC rotations.  It is a
  good model for a Markov partition/transfer-operator proof of decorrelation:
  split stable/free directions from unstable/cycle corrections, then show a
  spectral gap outside resonant finite atlases.

* **Fubini-Study metric:** the right quotient metric for complex phase profiles
  because global phase and scaling are irrelevant.  The proof target becomes a
  projective-angle bound for the signed resonant profile.

* **Road coloring:** HYP-2698/HYP-2702 are residual-profile automata.  A
  synchronizing coloring is the finite-state certificate that every generated
  residual context collapses into the safe middle/cap cone.

* **Hopfield/Hebbian learning:** the symmetric Hebbian weight matrix is a
  co-firing matrix of missed sectors or relation packets.  Diagonalizing the
  empirical density matrix of obstruction profiles should reveal the dominant
  low-rank modes that spend death-chain currency.

* **Propagator weights / double slit:** far-runner sector ownership acts like a
  class/slit; oscillatory kernels are the signed phase contributions from
  endpoint classes.  Interference is the signed deviation from the
  decorrelated cover.

* **Clifford plus `T`:** Clifford/stabilizer corresponds to the classically
  tractable cut-space/decorrelated quotient.  The `T` gate is the relation
  lattice/cycle-space correction.  The LRC analogue of `T`-count is the number
  and height of non-Clifford relation defects that must be routed to a finite
  atlas.

* **Crossing number of complete graph vs alternating graph:** HYP-2627's K14
  squarefree crossing carrier is the full interaction carrier.  An alternating
  subgraph should be treated as a parity/stabilizer skeleton; the difference
  is the signed magic/cycle correction.

## Next Sharp Targets

1. **FS-angle deviation lemma.**  Define a normalized complex phase profile
   `psi_E` over the mod-7 relation/far-packet basis and a decorrelated profile
   `psi_0`.  Prove an explicit bound

   ```text
   d_FS(psi_E, psi_0) <= theta(relation height, Freiman dimension, far gaps)
   ```

   strong enough that the deviation bound above is below the HYP-2701
   two-far boundary margin.

2. **Synchronizer cone lemma.**  Upgrade the all-singleton deletion automaton
   to the generated residual-profile cone of HYP-2698/HYP-2702.  Show every
   actual decorrelated LRC context enters the middle/cap cone before the
   finite resonant atlas.

3. **Magic-rank finite atlas.**  Treat low-height relation defects as
   non-Clifford resources.  High magic rank should be decorrelated by phase
   mixing; low magic rank should be finite and addressable by the existing
   cube-root, Freiman, AP-tail, and K14 squarefree carriers.

4. **Gibbs pressure comparison.**  Turn the cap margin into a pressure gap:

   ```text
   pressure(decorrelated single-block) + pressure(relation correction) < cap.
   ```

   THM-557 supplies the single-block floor; OPEN-Q-108 is exactly the remaining
   correction aggregation.

## Incoming Mainline Signal After Rebase

While this note was being committed, three relevant mainline additions landed.
They sharpen the route rather than competing with it.

* **HYP-2706 / death-chain band automaton:** singleton death-chain margins stay
  positive on the sparse-coordinate frontier, but first-order hit-count
  dominance and per-band monotonicity fail.  This is exactly the warning
  HYP-2705 needs: the projective phase profile must preserve signed band data
  until the final aggregate comparison.

* **HYP-2707 / tournament Clifford-magic hierarchy:** the Clifford/T analogy is
  now precise on the tournament side.  Score data is degree `1` in tile bits,
  `c3` is degree `2` and Gauss-sum/Clifford tractable, while `c5`, higher OCF
  layers, and `H` are higher-degree magic.  Its proved core goes further:
  `c3 mod 2` is a GF(2) quadratic form whose bilinear rank is
  `2*floor((n-1)/2)`, giving THM-555's `E[(-1)^c3]` by the standard Clifford
  Gauss-sum rank formula.  HYP-2705 imports this as a template for LRC: the
  decorrelated cut-space quotient is stabilizer data; signed relation/cycle
  corrections are the magic resource.  The LRC "magic rank" should therefore
  measure polynomial/phase degree or stabilizer rank of the mod-7 relation
  profile, not raw relation count.

* **Crossing/road-coloring scout:** circular tournament drawings have raw
  crossing number `binom(n,4)` independent of orientation, and the tested
  alternating crossing counts are not `c3`.  Likewise, road coloring fails
  literally on the tiling hypercube because every tile-letter is a bijection.
  This corrects HYP-2705's synchronizer language: the road-coloring object is
  the non-invertible residual-sector deletion automaton, not the invertible
  wiggly/tiling automaton.

* **Formal-group n-fold sum:** the machine-checked Cayley transform identity
  `Q(Fsum)=prod Q(a_i)` and unconditional `E+-O=prod(1+-a_i)` split offer a
  clean algebraic model for diagonalizing repeated update kernels.  The
  interior-pole caveat is a useful guardrail for LRC: a global denominator
  condition is not enough if a suffix/partial transfer hits a pole.

## Tournament Analysis

Vertices are proof lenses, not runners or arcs:

```text
fs_middle_projection
road_coloring_sync
gibbs_currency_gap
hebbian_density_matrix
propagator_kernel
clifford_stabilizer_cutspace
cat_map_markov_partition
crossing_squarefree_carrier
raw_analogy
```

Pairwise observable: how much LRC cap/survival data the lens preserves before
scalarization.  Switch/gauge: orient toward preserving the cap predicate and
signed phase labels.  The scout tournament is transitive:

```text
fs_middle_projection
> road_coloring_sync
> gibbs_currency_gap
> hebbian_density_matrix
> propagator_kernel
> clifford_stabilizer_cutspace
> cat_map_markov_partition
> crossing_squarefree_carrier
> raw_analogy
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
hamiltonian_path_count=1
```

Challenged assumption: a creative analogy is useful only after it names the
quotient preserving the LRC cap predicate and the data it destroys.  For this
session, the useful quotient is the projective residual phase profile plus the
residual-sector synchronizer automaton.

## Source Notes

External prompts inspected:

* Batched quantum state exponentiation and quantum Hebbian learning:
  `https://arxiv.org/abs/1803.07039`
* Road coloring theorem:
  `https://arxiv.org/abs/0709.0099`
* Feynman propagators as neural weights:
  `https://medium.com/@denizaskin/feynman-propagators-as-weights-a-quantum-neural-network-architecture-b65415ff41cc`
* Fubini-Study metric formula:
  `https://en.wikipedia.org/wiki/Fubini%E2%80%93Study_metric`
* Arnold cat map:
  `https://en.wikipedia.org/wiki/Arnold%27s_cat_map`
