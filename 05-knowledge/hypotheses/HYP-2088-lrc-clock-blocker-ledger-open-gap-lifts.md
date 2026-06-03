---
id: HYP-2088
status: SUPPORTED by S573 exact bounded audit; global classification open
source: codex-2026-06-03-S573
related:
  - HYP-2087
  - HYP-2086
  - HYP-2085
  - HYP-2084
  - HYP-2083
  - HYP-2082
  - HYP-2055
  - HYP-2052
---

# HYP-2088: strict sub-edge LRC rows must satisfy a three-clock blocker ledger

**Claim.** The S572 floor-only reading is too strong once lifted rows are
included.  The sharper necessary statement is a three-gate clock-blocker
ledger.  For `k=n-1` speeds, any row with

```text
M(S) < 2/(2n-1)
```

must simultaneously satisfy:

```text
D_q: for every 2 <= q <= n-1, some speed is divisible by q;
U_a: every unit antipodal shell {a,2n-1-a} is hit mod 2n-1;
N_j: every n-clock tick j/n has some runner at distance 0 or 1/n.
```

**Why these are bounds.**
- If `D_q` fails, then `t=1/q` scores at least `1/q >= 2/(2n-1)`.
- If `U_a` fails, S553 gives the inverse-clock witness at the edge.
- If `N_j` fails, then `j/n` scores at least `2/n > 2/(2n-1)`.

So every strict sub-edge row must block all three independent clock families.

**S573 exact audit.** `lrc_second_gap_bounds_s573.py` applies these gates and
then exact-checks the surviving rows in expanded primitive boxes:

```text
k=3,B=60; k=4,B=40; k=5,B=32; k=6,B=26;
k=7,B=22; k=8,B=20; k=9,B=18; k=10,B=17.
```

All-gate survivors exact-checked:

```text
2623, 5676, 3237, 1446, 4654, 238, 145, 111.
```

**New correction.** The open interval `(1/n, 2/(2n-1))` is nonempty in these
expanded boxes:

```text
n=7: (1,5,6,11,16,17), M=5/33, witness t=10/33.
n=8: (1,2,3,4,5,7,18), M=3/23, witness t=4/23.
n=8: (1,3,4,5,7,13,18), M=3/23, witness t=4/23.
```

Thus the repo should stop treating `2/(2n-1)` as a global second value over
all lifted integer representatives.  It is still a powerful edge for unit-shell
and lower-residue reductions, but lifted flip-set/nonunit-hole rows can sit
strictly between the floor and that edge.

**Interpretation.**
- Addition gives the antipodal shell ledger.
- Unit multiplication gives the `U_a` visibility gate.
- Divisibility gives the `D_q` small-denominator gate.
- The even/floor clock gives the `N_j` gate.

The proof target is therefore not "below edge implies floor."  It is:

```text
three-gate ledger => M(S) >= 1/n,
```

plus a classification of the lifted open-gap families.

**Connection to incoming Burnside/Fourier work.** HYP-2085 and HYP-2087 show
that binary lonely time words are too coarse but their Fourier/dihedral shadows
still expose clock families.  The `D_q/U_a/N_j` ledger is the labelled proof
obligation layer that should be attached to those time words: it records which
denominators, unit shells, and floor ticks are actually being blocked before
the binary word forgets ownership.

**Honest scope.** This is a bounded exact correction and a stronger necessary
ledger, not a global theorem.  It refutes HYP-2084 as a global floor-only claim
but preserves its local bounded observation inside the original smaller boxes.

**See:** `04-computation/lrc_second_gap_bounds_s573.py` (+.out),
`07-reflections/lrc-clock-blocker-bounds-open-gap-lifts-s573.md`, HYP-2084,
HYP-2083, HYP-2082, S553/S553b.
