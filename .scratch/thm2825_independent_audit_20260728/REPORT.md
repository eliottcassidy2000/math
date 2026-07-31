# THM-2825 independent hostile audit

## Verdict

**ACCEPT the proposed collar theorem after keeping its metric, direction,
and type boundary explicit.**  The exact scale-selected collar is a real
new structure, not just a numerical coincidence:

- the independently reconstructed `7 x 9 x 9` bank has `193` nonempty and
  `374` empty cells, `63,308` weighted common pieces, `587` weighted
  pulled-target right pieces, and `195,517` common/right pairs;
- every piece has length `26,444,880=6h/91`, its cell's uniform weight,
  and one common half-step lattice, where `h=401,080,680`;
- the unscaled relation in each cell is complete bipartite and therefore
  never unique (`189` to `526` common choices per right piece);
- centered circular distance selects a unique common piece at `R+h` for
  every one of the `587` right pieces, collision-free and in the
  right-to-common direction;
- delayed content agrees exactly for even half-step offset and reverses
  exactly for odd offset on all `195,517` pairs.  Thus `+h` reverses all
  `587`, while `+2h` preserves all `587`.

The theorem must not be promoted as a typed common-atom cospan or an
iterable endomorphism.  The `+h` and `+2h` images lie in `M`, while the
domain is `R` and `M intersect R` is empty.  Reverse nearest-neighbour
selection has `242` exact ties among the `63,308` common pieces.  The
odometer wrap satisfies

```text
13*tau=14h=7*(2h),
```

but direct `+14h` lands in `M` for only `527/587` right pieces and leaves
the cell support for the other `60`; it preserves delayed content on all
`527` survivors.  Hence the scale identity does not supply a uniform wrap
map.

## Typed sidecars

The full ancestry-event audit succeeds:

```text
+h:   same ancestry chamber 587/587
+2h:  same ancestry chamber 587/587
+14h: same ancestry chamber 527/527 surviving
```

A sampled cell on each live clock gives the same sharp factor/carrier
hostile.  In the first witness

```text
cell=(clock,s,t)=(1,0,3),
R=(142004190428100,142004216872980),
M=R+h=(142004591508780,142004617953660).
```

The right piece's source-native and target-adjacent `E3` bits are false,
whereas every native/pulled bit of the common piece is true.  Its
thirteen-twist carrier mask is `(source empty,target delta_0)`, whereas
the common image is `(source delta_0,target delta_0)`.  Thus translation
repairs the missing source carrier rather than preserving it.  The same
failure occurs for `h`, `2h`, and surviving `14h` on all six sampled right
pieces.

The separate target-twelve wrap audit in
`.scratch/cofiber_wrap_flux_20260728` additionally found no global fixed
or translated endpoint-address preservation.  A selected full-bank endpoint
probe remains a useful secondary audit, but is not needed for the metric and
semantic theorem above.

## Prime-power spectrum sidecar

Splitting each cell's half-step lattice by parity:

- both common colour masses are nonzero modulo `13` in all `193` cells, so
  THM-2839 applies separately to both `C_(13^5)` common-colour masks;
- both right colour classes are nonzero in only `14` cells; in `179` cells
  all right pieces have one parity;
- the two `C2` character augmentations are units for the right mask and the
  combined `M+R` mask in all `193` cells;
- the common mask has full `C2 x C_(13^5)` unit status in `88` cells and
  fails it in `105`, exactly because even-minus-odd is zero there;
- no ordinary total augmentation vanishes.

These are signed group-algebra consequences only.  They do not construct a
positive physical inverse, endpoint current, or LRC(14) exclusion.

## Reproduction

```bash
python3 .scratch/thm2825_independent_audit_20260728/audit.py \
  > .scratch/thm2825_independent_audit_20260728/audit.normal.out
python3 -O .scratch/thm2825_independent_audit_20260728/audit.py \
  > .scratch/thm2825_independent_audit_20260728/audit.optimized.out
cmp .scratch/thm2825_independent_audit_20260728/audit.normal.out \
    .scratch/thm2825_independent_audit_20260728/audit.optimized.out
```

The normal transcript is the stored transcript.  Normal and optimized
outputs are byte-identical.

```text
audit.py                 ad1a6323a735703a5c0a95dd0d692eafd77f7e6495945b12367af35310cb38cf
audit.normal.out         b9846e35c615b6ff0183f6df47ba8e8f99796a195afb6480f922254e35cbe82d
audit.optimized.out      b9846e35c615b6ff0183f6df47ba8e8f99796a195afb6480f922254e35cbe82d
Python assertion nodes   0
```
