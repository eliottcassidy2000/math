---
id: HYP-2064
status: PROGRESS; exact carryover invariant found, proof consequence open
source: codex-2026-06-02-S560
related:
  - HYP-2063
  - HYP-2061
  - HYP-2059
  - THM-396
  - THM-369
---

# HYP-2064: the n=17 skip-8 packet transfers to n=14 as a predecessor-of-apex skip-6 packet

## Statement

The `n=17` prime-gate attempt isolated a "skip the middle gate" packet:

```text
{r} union {17*q : 1 <= q <= 16, q != 8}.
```

Returning to `n=14` shows the same mechanism with a composite correction.  The
best row-parent, gate, and double-gate packets are:

```text
{r} union {7*q  : 1 <= q <= 13, q != 6}
{r} union {14*q : 1 <= q <= 13, q != 6}
{r} union {28*q : 1 <= q <= 13, q != 6}
```

for primitive breaker choices `r`.  Thus the n=17 half-gate `8` carries back to
the `2*7` row as the predecessor-of-apex gate `6`: the apex/bridge `q=7` is kept
as a shield, and the missing middle packet shifts one step left.

## Evidence

`lrc_n14_n17_gate_skip_transfer_s560.py` scans exact gate-packet rows and finds:

```text
n=14 scale 7  skip 6: gap/th = 5/924
n=14 scale 14 skip 6: gap/th = 5/1848
n=14 scale 28 skip 6: gap/th = 5/3696
n=17 scale 17 skip 8: gap/th = 1/272
```

For `n=14`, the scale chain `7 -> 14 -> 28` preserves `skip=6` and halves the
exact visible gap at each gate-depth step.  The forbidden length remains
`142/143`; the boundary count doubles from `84` to `168` to `336`.  This is the
old endpoint-debt export phenomenon in a sharper gate-skip language.

The top-gate unit-wall lemma also transfers: no `n`-multiple leaves every unit
wall `a/n` as a closed witness.  So the n=17 "prime gate" was not prime-only;
it was the top case of THM-369.  What is special at `n=14` is the extra
`q=7` apex/bridge, which changes the optimal skipped gate.

## Tournament Analysis

Assumption challenge: the useful vertices are not runners.  S560 considered
vertices including runners, gate multiples, skipped labels, endpoint leaves,
small-pinch pair-cells, and whole proof obligations.  It chose whole
gate-packet rows as vertices because the preserved predicate is
counterexample-likeness of a repair packet.

- **Vertices:** best gate-packet rows: `n14` scale `7`, `14`, `28`, and `n17`
  scale `17`.
- **Pairwise observable:** `(gap/th, missing moduli, forbidden length, scale,
  skip)`.
- **Switch/gauge:** smaller gap wins; ties prefer sieve-complete rows and
  longer coverage.
- **Tie Hamiltonian path:** lexicographic row label after the observable.
- **Fingerprint:** transitive, score histogram `{0:1,1:1,2:1,3:1}`, no directed
  3-cycles, singleton SCCs, one Hamiltonian path.

What the quotient preserves: ladder-level gap/debt ordering and the skipped
gate label.  What it destroys: endpoint ownership, small-pinch shield graph
incidence, and runner-level CRT class coupling.  The next proof step must put
those labels back.

## Consequence for the n=14 proof search

The n=17 carryover suggests that the `n=14` near-miss is not an arbitrary
seven-ladder.  It is a stable gate-packet:

```text
missing q=6, retained apex q=7, exported debt under scale 7 -> 14 -> 28.
```

This should be combined with HYP-2061/THM-396.  The retained apex is the natural
sum-multiple shield; the skipped predecessor is the corridor where endpoint
debt leaks.  A plausible next theorem is:

> In the `n=14` primitive gate-packet family, any packet that keeps the apex
> shield and skips a predecessor gate has exact positive gap and exports
> endpoint debt under every further `2`-adic lift.

## Sources

- `04-computation/lrc_n14_n17_gate_skip_transfer_s560.py`
- `05-knowledge/results/lrc_n14_n17_gate_skip_transfer_s560.out`
- `07-reflections/lrc-n14-n17-gate-skip-transfer-s560.md`
- HYP-2063
- HYP-2061
- THM-396
