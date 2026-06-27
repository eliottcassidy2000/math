# LRC14 Observer-Gluing Ledger Back And Forth

*codex-2026-06-27-S258. Continuation of the proof frontier after incoming
HYP-3096 and HYP-3097.*

## Task A: witness-route ledger

I built `04-computation/lrc14_observer_gluing_ledger_codex_s258.py` to make
HYP-3096 concrete.  The script computes exact positive-length components of
the direct lonely set

```text
L(S) = {t : ||s_i t|| >= 1/14 for all i}
```

for representative rows: a 14-free q-witness row, a THM-573 `H7>=7` row, an
`H7=6` boundary residual, apex-family rows `{1,...,12,14V}`, and THM-575
divisor-loaded rows.

The saved output is
`05-knowledge/results/lrc14_observer_gluing_ledger_codex_s258.out`.

The main exact readout is:

```text
q-witness row {1,...,13}:       terminal grid exit at t=1/14
H7=7 row:                       THM-573 exit, largest direct arc=3/343
H7=6 residual row:              components=42, largest=19/1372
apex V=13 residual:             components=24, largest=3/637
apex V=200 residual:            components=102, largest=3/9800
divisor-loaded B=8 residual:    components=860, largest=1/82320
```

This advances HYP-3096 by making the direct component/floor obstruction
quantitative.  It also corrects the next proof target: a global direct largest
arc theorem is the wrong scalar after THM-575.  The proof should split into a
bounded-apex direct packet and a large-apex normalized slow/ruler packet.

## Task B: pair/Pascal scissors sidecar

The same script attaches HYP-3097 fields to the same rows:

```text
H7_pair_shadow = C(count_7_divisible,2)/91
even_pair_shadow = C(count_even,2)/91
mod7 and mod14 residue-count signatures
Farey p+q and p*q lanes mod 91 when a root is named
```

This is the back step from HYP-3097 to HYP-3096: the pair/Pascal data is not a
separate numerology track.  It tells the witness ledger which rows are still
different after scalar residual labels have been applied.

The sample has seven live observer-glue rows but five distinct mod-7 scissors
signatures.  The apex family keeps one scissors signature as `V` changes, while
the direct largest arc collapses; the divisor-loaded rows carry a different
signature and shatter much more severely.  So the proof should not quotient by
terminal status, direct measure, or pair mass alone.

Incoming kps-S31ag fits the same warning from a different side: the coarse
mod-14 winding tournament is degenerate for every `k>=8`, since an antipodal
pair `{r,r+7}` is forced and gives a permanent tie.  Therefore coarse H cannot
be the coverage-extremality bridge at the binding rows.  The S258 ledger should
keep fine-scale mod-`p` or packet-level scissors sidecars whenever tournament
data is used.

## Combined Frontier

The next two tasks should continue alternating:

1. Expand the S258 script from representative rows into the HYP-2963 packet
   bank and outside-bank normalizer attempts.
2. For each packet, record which chart pays for the destroyed coordinate:
   direct arc, normalized slow/ruler arc, moment/Perron, branch/K33, or named
   finite denominator debt.

Tournament Analysis uses observer charts as vertices:

```text
observer_gluing_packet
> direct_lonely_component_packet
> sector_pair_scissors_packet
> crt_c7_level7_lift
> crt_c2_dyadic_lift
> pascal_pair_mass_shadow
> raw_direct_arc_scalar
> raw_I_table_enumeration
```

The challenged assumption is now explicit: neither raw direct arcs nor raw
pair/Pascal caps nor coarse mod-14 winding tournaments are proof vertices until
their overlap maps with the other observer charts have been checked.
