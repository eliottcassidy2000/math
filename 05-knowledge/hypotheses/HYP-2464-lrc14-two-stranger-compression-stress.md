# HYP-2464 - LRC14 two-stranger compression stress still opens in Q27

**Status:** OPEN synthesis with exact bounded support lemma.

**Source:** codex-2026-06-13.  Extends HYP-2463, HYP-2444, HYP-2443, and OPEN-Q-082.

## Claim

The HYP-2463 hard-residue hull is not the whole two-stranger story, but its
resource-independence lesson persists.

Let

```text
CORE = 7*{1,...,12},     MAX_R = 13*84 = 1092,
Q27 = {d*m : d in {1,2,7,14}, 1 <= m <= 27} \ {1}.
```

Delete one runner from `CORE` and add two distinct non-core speeds
`r1,r2 in [1, MAX_R]`, keeping only primitive rows.  Among all `6,868,368`
such true two-stranger rows, every row that blocks all plain shells `q<=27`
has a Q27 witness.

This is a bounded lemma, not a proof of LRC14.  Its proof value is that the
plain-shell residuals are highly structured: they all spend a `13`-clock speed
and a visible deleted-core address before Q27 opens.

## Exact Evidence

Script: `04-computation/lrc14_two_stranger_compression_stress_codex.py`

Stored output: `05-knowledge/results/lrc14_two_stranger_compression_stress_codex.out`

The scan covers:

```text
12 * (binom(1080,2) - binom(144,2)) = 6,868,368
```

primitive rows.  Results:

- `877` rows block every plain shell `q<=27`.
- `0` of those `877` miss Q27.
- First-Q27 histogram:

```text
{30:44, 32:88, 34:25, 36:132, 38:267, 40:130,
 44:5, 54:7, 70:93, 84:84, 91:1, 161:1}
```

- `636/877` residuals use zero old hard residues from HYP-2444/HYP-2463.
  Thus the old hard-atom hull is not a sufficient compression target by itself.
- Every residual has at least one added speed divisible by `13`; `858` have
  exactly one, and `19` have two.
- The deleted-core histogram is

```text
{14:8, 28:48, 35:4, 42:139, 56:16, 63:15, 70:240, 77:17, 84:390}.
```

No residual deletes `7`, `21`, or `49`.

The two latest one-off divisor-fiber rescues are:

```text
delete 35, add (702,757),  q=91=7*13,  twist 1
delete 70, add (910,1080), q=161=7*23, twist 16
```

The large residual blocks are also sharply addressed:

- `q=70` occurs only after deleting `70`; `83/93` such rows have mod-13 pair `(0,1)`.
- `q=84` occurs only after deleting `84`; its top shell-27 antipodal pairs involve class `10`.

## Interpretation

HYP-2463 said the old shell-27/13-clock hard packet cannot stack.  HYP-2464
adds that the next broader two-stranger residual class is also resource-typed,
even when neither stranger is one of the old hard residues.

The apparent obstruction has three compulsory coordinates:

1. **A 13-clock debt.**  Every residual uses an added speed in residue `0 mod 13`.
2. **A deleted 7-core address.**  Only nine of the twelve core speeds can be
   deleted in residual rows; `7,21,49` never appear.
3. **A Q27 divisor-fiber escape.**  Every residual opens at `2*m` or at a
   `7`/`14` fiber, with the latest outlier `161=7*23`.

This shifts the compression theorem.  A would-be arbitrary Q27 blocker need not
literally compress to the old hard-residue hull; it should compress to a
resource pattern with `13`-clock debt, deleted-core address, shell-27 pair data,
and divisor-fiber escape.  If that compression fails, the failure should be a
low clock, AP/Vstar/2AP descent, or odd owner/Bprime channel.

## Assumption Challenge / Tournament Analysis

This session again refused to use runners as the default tournament vertices.
Candidate vertex sets considered:

- runners,
- deleted core addresses,
- added-residue pairs,
- denominator fibers,
- unit twists,
- shell-27 classes,
- 13-clock debts,
- owner/Bprime openings,
- AP/Vstar/2AP descent targets,
- proof obligations.

The selected Tournament Analysis uses compression maps as vertices.  It
preserves exact finite blocker laws and proof leverage, while destroying most
raw time geometry.  The stored tournament has score histogram
`{0:1, 2:3, 4:1, 5:1, 6:1}`, one directed 3-cycle, SCC sizes
`[3,1,1,1,1]`, and three Hamiltonian paths.  The leader is
`divisor_fiber_Q27`, followed by `thirteen_clock_debt` and
`deleted_core_address`.

The challenged assumption is that the one-stranger hard residues are the only
atoms worth tracking.  The bounded two-stranger scan refutes that assumption:
most residual rows are non-hard, but the same resource currencies remain
visible.

## Next Proof Route

1. Prove a finite quotient lemma: rows blocking all plain `q<=27` shells must
   contain a `13`-clock debt plus a deleted-core address, after allowable
   compression/descent moves.
2. Upgrade the bounded `MAX_R=1092` result using the same kind of large-speed
   Bprime/Cor-B2 escape used in HYP-2444 for one stranger.
3. Add owner/Bprime certificates for the `877` residuals, but with a faster
   interval engine than the current Fraction-based one.
4. Replace the old target "compress to hard residues" by the broader target
   "compress to resource coordinates"; then use HYP-2463 and HYP-2464 together
   as the finite obstruction atlas.
