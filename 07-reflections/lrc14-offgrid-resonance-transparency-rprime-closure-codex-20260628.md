# LRC14 Off-Grid Resonance Transparency And Rprime Closure

Anchors: HYP-3421, HYP-3418, HYP-3415, HYP-3266, HYP-3310, HYP-3255, HYP-3140,
HYP-3136, HYP-3129, HYP-3125, HYP-3124, HYP-2896, OPEN-Q-108.

HYP-3415 identifies the critical path: q-witness closes non-covering rows,
LRC<=13 closes the `Q` half of covering rows, and one decorrelation floor
`R' > 0` remains.  HYP-3418 then corrects the tempting shortcut: the
coprime-to-14/non-resonant subpacket is not enough, because its witness is
near `t=1/2` and that kills even speeds.  The useful move in HYP-3421 is
therefore narrower and sharper: stop treating "resonant survivor" as its own
analytic beast at the full optimum.  The exact computations say something
more directional:
resonance is a grid-local danger, while the useful witnesses are off-grid.
At those witnesses the resonant speeds are not merely harmless on average;
they are visibly safe by exact rational distances.

The named rows expose the pattern.  The zeta-12 seed has `M=1/12`, the
bounded even covering row has `M=1/8`, the many-`14Q` rows show `M=1/9` and
`M=1/8`, and the canonical `{1..11,13,84}` row has `M=7/89` at `37/89`.
All selected optima are off the 14-grid.  Every speed divisible by `2` or `7`
is safe there; every `14Q` tip in the checked rows is safe there.  Conversely,
the `14Q` tips kill every unit grid point `a/14`, exactly where AP/GW
equioscillation lives.

That turns the geometry into a core/bulk split:

```text
unit grid core  = AP/GW finite equioscillation, three antipodal binding pairs
off-grid bulk   = transparency plus Rprime decorrelation
even-speed floor = 2-adic descent under u=2t
14Q resonance   = grid-local danger, off-grid guard flank
```

The canonical `84m` tower is the cleanest theorem fragment:

```text
t_m = (35m+2)/(84m+5)
M({1..11,13,84m}) = 7m/(84m+5)
active pair = (5,84m)
```

Here resonance does not survive as an enemy.  It becomes one of the binding
flanks that certifies the margin over `1/14`.

The strongest proof-facing obligation is therefore not "handle resonant
survivors analytically."  It is the finite classifier:

```text
every residual packet
  -> q-witness / denominator-floor transparency
  -> canonical binding formula
  -> positive owner packet
  -> Rprime row with signed-SPEC/fiber-PGF certificate
  -> named terminal debt
```

After that classifier, the analytic residue is no longer "non-resonant
decorrelation only."  It is HYP-3418's 2-adic even-speed descent plus the
HYP-3129/HYP-3140 constant chase: prove the closed-form all-packet lower bound
for `Rprime = E[N_R | Q]/E[N_R]`.  The edge-witness language remains necessary
because `Rprime` is not a scalar shadow; it needs the tail `R-safe` packet,
the tip `Q-safe` packet, both deletion children, and the signed covariance
payload before quotienting.

Assumption challenge: tournament vertices were explicitly not runners.  The
candidate vertices considered were runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
witness times, off-grid cells, `14Q` tips, PGF coefficients, and proof
obligations.  The useful vertices in this pass are proof carriers.  The
preserved predicate is existence of a time with all speeds at distance
`>=1/14`.  The destroyed information in the raw resonance quotient is the
off-grid witness, endpoint owners, tail/tip child packets, sheet counts, and
SPEC low/tail constants.

This does not close LRC14.  It makes the next closeable statement sharper:
prove the all-packet full-optimum off-grid transparency lemma, then finish
the 2-adic descent and symbolic `Rprime` constant chase.
