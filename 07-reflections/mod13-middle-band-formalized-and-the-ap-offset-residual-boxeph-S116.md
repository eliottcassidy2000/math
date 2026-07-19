# The mod-13 middle-band witness is formalized; the AP offset residual remains open

*boxeph-2026-07-18-S116, corrected during codex-2026-07-18-S75 review.*

The durable deliverable is `LRCMod13Blocking.lean`: a kernel-checked witness at
times `b/13`.  The original session title and finite script overstated what that
witness proves.  In particular, neither the script nor the Lean module proves
that every tight twelve-speed arithmetic progression has `a=d`, and neither
proves the general n=12 equality classification.  Those are still parts of the
offset-vanishing residual.

## 1. What the Lean module proves

In `namespace LonelyRunner`:

- `mod13_middle_far`: if `2 <= r <= 11`, then `|13k+r| >= 2` for every
  integer `k`.
- `sieve13_middle_witness`: if every residue `(c_i b) mod 13` lies in
  `[2,11]`, then at `t=b/13` every runner is at distance at least `2/13`
  from every integer.
- `no_middle_band_witness_of_tight`: the corresponding contrapositive.

Consequently, if a family has global margin at most `1/13`, then for every
integer `b` some residue lies outside the middle band, hence in
`{0,1,12}`.  If one additionally knows `13` divides no speed and takes
`b` to be a unit modulo 13, the zero case is impossible and this becomes
the familiar conditional pair-blocking conclusion

```text
for every b in (Z/13Z)^*, some c_i b is congruent to +/-1 mod 13.
```

The nonzero-residue hypothesis matters.  A residue `0` makes the runner
close at `b/13`; it does not contradict tightness.

## 2. The valid arithmetic-progression slice

Let

```text
C(a,d)={a+kd : 0<=k<=11}.
```

Suppose the twelve residues of `C(a,d)` modulo 13 are all nonzero.  Then
`d` is a unit modulo 13 and the twelve-term progression is the full residue
cycle with one class omitted.  The omitted class is `a-d`.  Therefore

```text
all C(a,d) residues are nonzero  iff  a is congruent to d mod 13.       (1)
```

This is an exact residue calculation, but it is conditional.  The middle-band
lemma alone does not show that a tight AP has no zero residue, and (1) does not
upgrade `a congruent d (mod 13)` to the integer identity `a=d`.

There is one additional unconditional half that is elementary and was already
available from the translate-block witness.  At

```text
t=1/(2a+11d)
```

all twelve speeds have distance at least `a/(2a+11d)`.  Hence `a>d`
implies a margin strictly larger than `1/13`; only the `a<=d` half remains
for an AP classification.  The S116 computation samples that residual but
does not close it.

## 3. What the finite script says

`lrc14_n12_tight_locus_boxeph_S116.py` enumerates reduced rationals with
denominator at most 300 for eleven displayed pairs `(a,d)`.  Thus its
nonhomogeneous values are certified lower bounds for the true supremum, not
exhaustive upper bounds and not a classification over all `a,d`.  The rows
are useful counterexample/consistency telemetry:

- homogeneous dilates `c*{1,...,12}` reproduce `1/13` (whose exactness is
  known independently by dilation invariance and the classical AP case);
- the displayed shifted APs already have witnesses strictly above `1/13`;
- every displayed tight row satisfies `a=d`, but the finite list cannot prove
  the universal converse.

## 4. Honest frontier

Delivered here: a machine-checked mod-13 middle-band witness, its precise
conditional pair-blocking consequence, the residue identity (1), and the
elementary dispatch `a>d` for twelve-term APs.

Still open here: exclude every `a<d` nonhomogeneous AP uniformly; remove the
nonzero-residue condition; prove `a=d` rather than only a mod-13 congruence;
and classify arbitrary tight twelve-speed families.  Those are all-moduli
or offset-vanishing statements, not consequences of one prime slice.

Cross-links: S115's conditional mod-13 pair blocker, HYP-4382,
THM-635, `LRCMod13Blocking.lean`.
