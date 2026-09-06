# Independent norm-five profile and compact-base audit

**Status: INDEPENDENT AUDIT PASS.** This accepts the norm-five portions
(Sections 1--5) of `overnight5_20260906_lrc_norm5_profiles.md`. The sharp
selected and physical conclusions, their equality triples, and the physical
bulk identity are already proved in incoming
`01-canon/theorems/THM-4441-lrc14-signed-122-sharp-ray-closure.md`, read at
`origin/main` after `058a8ded9`. This is independent concurrent confirmation
and a compact proof-base reconstruction, not a new sharp-closure claim.
The producer's separate generic coefficient/H535 route is outside this
numerical audit; incoming THM-4437 already supplies the stronger generic
all-height result.

## Full-ray scope and normalization

For a primitive sorted positive ternary-unit triple, the norm-five relation
has three unit coefficients modulo three. Its integer defect obeys
`|delta|<15/14` and `delta=0 mod3`, so it is zero. Primitivity then makes the
complete carrier line exactly `C=k v`, with unit longitudinal step; the live
word deletes precisely multiples of three. Every strict roof is retained.
This proves full-ray completeness, including empty support, without an
unchanged-prefix or selected-direction assumption.

Positivity and sorting give exactly the four displayed sectors. Independently
enumerating signed permutations of `(1,2,2)` agrees with this classification
on every finite-base row. The normalized parameter in our sector III is
`a/c=t`; THM-4441's F4 uses `b/c=2-2t`, so its crossing at `7/8` becomes
our `t=9/16`. The sector correspondence is I/II/III/IV = F1/F3/F4/F2.

With `r=3/14`, `s=k/(rc)` and normalized speeds `(t,b/c,1)`, let
`g_i(s)=min(2,(S_i-|v_i|s)/P_i)`. Every profile has the **common complete**
cutoff kappa. Integrating one cap and its affine continuation gives

```text
J_i=2*kappa-(max(|v_i|*kappa-S_i+2P_i,0))^2/(2|v_i|P_i),
selected continuum = (3/49)*min_i J_i.
```

The positive part in sector IV's third projection is necessary. Some
individual projections jump down at the full cutoff, because another roof
vanishes first. Dropping that cutoff changes the event being measured.

## Independent exact continuum proof

The verifier reconstructs every cap integral directly from its affine
profile, rather than importing the producer's rational formulas. It checks
all selector inequalities by exact Bernstein coefficients after affine
parameter substitution; subdivision, if needed, remains rational. Every
denominator is separately certified nonzero in the relevant open interval.
This replaces the producer's real-inequality solver path.

The selected normalized integrals are:

| Sector | Complete selector | Maximum or supremum |
|---|---|---:|
| I | `(2t+7)/8` before `1/8`, then `2t^2-t+1` | `29/32` |
| II | `(9t+14)/16` before `1/8`, then `(t^2-t+2)/2` | `121/128` |
| III | `(t+7)/8` before `9/16`, then `2t^2-3t+2` | `121/128` |
| IV | `(16t+7)/(8(2t+1))` throughout | supremum `15/16` |

The global selected continuum is therefore at most `363/6272`. Its two
real-ratio maximizers have a speed divisible by three; the continuum bound
does not claim a discrete equality there. The selected coordinate is fixed
for the sum defining a projection. It is not selected independently for
each carrier.

I also reconstructed the full physical lower envelope. In the following
table the profile is the cap 2 until the first breakpoint, the first named
roof until the second breakpoint, then the second named roof to kappa:

| Sector / parameter | First breakpoint | Second breakpoint | Roofs |
|---|---:|---:|---|
| I | `t/2+1/4` | `(1-t)/2` | a, c |
| II, `t<=1/2` | `t/2` | `(1-t)/2` | a, c |
| II, `t>=1/2` | `(1-t)/2` | `t/2` | b, c |
| III | `t-1/2` | `1-t` | a, c |
| IV | `1/4-t/2` | `1/2` | a, b |

For every piece and every parameter in its whole interval, all other
affine profiles lie above the named one at both endpoints. Exact sign
certificates therefore prove dominance throughout the piece. Coverage,
breakpoint ordering, and the active cutoff are checked explicitly. Direct
integration of every complete arrangement gives

```text
integral_0^kappa min_i g_i(s) ds=7/8,
full raw physical integral=2r^2*(7/8)=9/112,
physical continuum=(2/3)*(9/112)=3/56.
```

Thus the physical result is independently proved without a fixed roof and
without relying on rational sample points. The genuine sector-II switch
at t=1/2 is retained. These formulas agree with THM-4441's existing piecewise
proof, including its differently parameterized F4.

The producer's separate saturated-plane proof also passes: for an integral
Bezout vector q with `q.w=1`, `n1=v cross q` satisfies `w cross n1=v`.
The map `(y,u)->u*n1-y*w` has Jacobian `||v||`, and its y-fiber is the full
physical interval. The primitive cross product rules out an omitted lattice
index. Eliminating the coefficient-one coordinate leaves a square clipped
by `|e_j+-e_k|<=r/2`; subtracting the two corner triangles gives area
`7r^2/4=9/112`. This is consistent with the complete envelope calculation.

## Strict tails and independently reconstructed finite bases

The positive live multiplier count satisfies
`R_<(T)<2T/3+2/3`. Layer-cake integration applies to both physical and
individual profiles despite their strict-cutoff jumps. It gives

```text
min E <363/6272+4/(7c),
physical mass <3/56+4/(7c).
```

These pay the strict `11/140` selected target from c=28, the sharp selected
`46/665` target from c=51, and the sharp physical `51/770` target from c=46.
There is no equality at any discarded tail height.

The finite universe is generated independently by all signed permutations
of the relation magnitudes, not by importing the four-cone classifier.
It gives 174 rows through H50, including 131 through H45 and 48 through H27.
For each row, an explicit integral Bezout lift constructs all effective
intervals along the complete primitive ray. Pair-intersection widths give
all three projections; triple intersections give physical mass. A separate
literal three-sheet breakpoint predicate independently checks every one of
the 174 physical masses. Every projection, mass and support count agrees
with the existing native-audited complete H63 table.

The finite equality loci are exactly

```text
selected 46/665: (2,19,20),
physical 51/770: (1,11,20).
```

The switching-roof hostile `(10,11,16)` is independently recovered with
`E=(17/176,9/140,3/55)`, physical mass `331/6160`, and gap `1/1232` between
the selected projection and physical mass. No fixed-projection identity
was used to evaluate or bound the physical profile.

These compact bases prove the two global constants. They do not prove
THM-4441's finer per-cone sharp constants from too small a head; that
theorem's larger exact bases retain their separate role.

## Incoming theorem comparison and reproduction

THM-4437's proved generic statement is `E_i<=6/77` for each projection and
strict inequality for their minimum, outside the three low circuits. Its
three individual equalities must be preserved; strictness of every generic
projection would be false. THM-4441 closes the norm-five circuit below
`6/77`. Along with the separately audited sharp norm-four closure and
excluding the additive family, these inherited results already yield the
global nonadditive `11/140` bound. No new large raw census is required for
that implication, and none was performed by this audit.

The exact verifier is `overnight5_20260906_lrc_norm5_audit_exact.py` and
the matching transcript has suffix `.out`. From the repository root:

```text
python -B 04-computation/overnight5_20260906_lrc_norm5_audit_exact.py
python -B -O 04-computation/overnight5_20260906_lrc_norm5_audit_exact.py
```

Normal and optimized outputs are byte-identical; **1,960** explicit gates
remain active in both modes. SHA-256 values use LF bytes:

```text
source 81253b43f6f87e7ac9e96878dd31e54a7322cdcffa668ad0a337bd3084c154b0
output 778680adceee3b3697bdb9ce1a0f997d3cd3354da409e450f688f32817b7617b
producer source 38165ecdd4925b4015cf716ae0f741297428b3e06211a7a8403cbb95d4361f52
producer output cb5e46db7a5a6ea0bf84fbc9414be4c6e87eff3a92845169653fcea5fa6a959f
inherited H63 table c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa
```

No arbitrary body Haar floor, actual entry, synchronization, or LRC(14)
closure follows from this local norm-five result. No repository file or
Git state was changed during this fifth-round audit.

**Filing note.** Root filed this completed audit after checkpoint `07b2d91b2`.
Local paths were made repository-relative and output line endings normalized;
normal and optimized verification was rerun with matching output. The proof
and finite universes are unchanged.
