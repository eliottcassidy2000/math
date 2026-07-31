# Projected k=3 z246--z244 candidate checkpoint

Status: **SCOUT / EXACT LOCAL COMPUTATION; NOT CANONICAL AND NOT A PROVED
CAP.**  The ordinary/optimized compositor, pinned dependency hashes,
independent replay, and theorem integration do not yet exist.

Current proved dependency at the start of this descent is THM-2941 plus the
pending THM-2979/2981 packages.  Candidate inherited ledger after z247 is
`375251`.

## z246

Atlas body rows: `194`.  Exact denominator quotient:

```text
49578 = 18946 crude + 29116 common-status + 1516 residual.
```

The residual states lie on 26 bodies.  Exact terminal reduction:

```text
zero-high hostile passes = 1402
one-high cases           = 1752
body-distinct low pairs  = 212
two-or-more-high gaps    = strict positive on all 26 residual bodies
high-ray recurrence checks = 1084084
```

Cardinality-forced least-`r` and effective-order histograms agree:

```text
{2:1507, 3:126, 4:34, 5:53, 6:31, 7:1}
```

There are no forced-density failures.  The THM-2984 primitive-unit control
has `2242772` scalar-reachable unit/case incidences, and all `2242772` have an
absolute fixed-safe cell; there are no ray-resolved failures.

The 26 residual bodies are

```text
(1,3,8,10,12,14)       (1,3,10,11,12,14)
(1,4,10,11,12,14)      (1,6,8,9,10,14)
(1,6,8,10,11,14)       (1,6,8,10,12,14)
(1,6,8,10,13,14)       (1,6,10,11,12,14)
(1,8,9,10,12,14)       (1,8,10,11,12,14)
(1,9,10,11,12,14)      (2,3,8,10,11,14)
(2,3,8,10,12,14)       (2,3,8,11,12,14)
(2,3,10,11,12,14)      (2,5,8,11,12,14)
(2,6,8,10,11,14)       (2,6,8,10,12,14)
(2,8,9,10,11,14)       (2,8,9,11,12,14)
(2,8,10,11,12,14)      (2,9,10,11,12,14)
(3,4,10,11,12,14)      (3,8,10,11,12,14)
(3,9,10,11,12,14)      (6,8,9,10,11,14).
```

## z245

One atlas row, `E=(1,8,10,12,13,14)`:

```text
119 = 0 crude + 113 status + 6 residual.
```

The two-high gap is `613634317/162615937120`.  The six residual states give
six one-high cases and one body-distinct low pair.  Cardinality-forced
least/effective histogram:

```text
{2:1, 3:2, 4:2, 6:1}.
```

No forced-density failure occurs.  The THM-2984 control has
`80640/80640` scalar-reachable unit incidences safely located.

## z244

Two atlas rows.  The quotient is

```text
7 = 0 crude + 7 status + 0 residual.
```

Thus the combined **candidate only** consequence is cap `z1<=243` and ledger

```text
375251 - 194 - 1 - 2 = 375054.
```

## z243 startup boundary

The z243 atlas has `154` rows.  Exactly three have `high_floor <= 243`, so
their later-HIGH coordinate comes from strict label ordering rather than the
projected wall: every later label is then automatically high.

```text
E=(1,2,3,8,10,12),       high_floor=166
E=(1,3,4,8,10,12),       high_floor=166
E=(2,4,6,8,12,14),       high_floor=232
```

A direct exact quotient check closes all three *before* the terminal split:

```text
states/crude/status/residual = 1/1/0/0
states/crude/status/residual = 13/13/0/0
states/crude/status/residual = 1/0/1/0
```

Consequently every actual z243 residual, if any, still has
`first < high_floor`, where the strict projected wall forces a later high and
the one-high terminal reduction remains applicable.  This is an exact local
scout result, not yet an ordinary/optimized canonical z243 package.  The
general boundary rule is a dichotomy: below the floor the wall forces a later
high; at or above the floor ordering makes all later labels high.  A literal
`HIGH-TAIL` print token is only one possible maximizer representation and is
not the logical gate.

## Reproduction

The generic local drivers are untracked development files:

```text
.scratch/lrc14_k3_exact_bank_probe.py
.scratch/lrc14_k3_z275_density_torsion_probe.py
```

Representative commands:

```text
python3 .scratch/lrc14_k3_exact_bank_probe.py --first 246 --processes 8 --compact
python3 .scratch/lrc14_k3_z275_density_torsion_probe.py --first 246 --body ... --compact
python3 .scratch/lrc14_k3_exact_bank_probe.py --first 245 --processes 1 --compact
python3 .scratch/lrc14_k3_z275_density_torsion_probe.py --first 245 --body 1,8,10,12,13,14 --compact
python3 .scratch/lrc14_k3_exact_bank_probe.py --first 244 --processes 2 --compact
```

For z246, pass precisely the 26 residual bodies printed by the bank's
`RESIDUAL_SUMMARIES`.  The exact list was sent to the THM-2981 owner in the
same session; the future canonical compositor must derive it rather than
trust this prose.
