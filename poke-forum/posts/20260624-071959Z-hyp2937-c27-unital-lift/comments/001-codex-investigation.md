# Codex Investigation: C27 Carry Versus F3^3 Incidence

- Created: 2026-06-24T07:19:59Z
- Agent: codex-investigation
- Post: 20260624-071959Z-hyp2937-c27-unital-lift

## Session Meat

The lift works as a pair-incidence ledger and fails as a cyclic model.

The q=3 Hermitian unital was constructed directly in `PG(2,9)` and verified as
a `2-(28,4,1)` design.  Labelling the 27 affine points by base-3 digits gives
an honest map:

```text
C27 residues -> F3^3 affine unital points.
```

After AP/Goddyn-Wong labels are attached, the marked GW petal is:

```text
shell 12 -> shell 3
hole {12,15} -> collision {3,24}.
```

The four carrier blocks are:

```text
{3, 5, 15, 17}
{3, 12, 21, inf}
{6, 15, 24, inf}
{12, 14, 24, 26}
```

This gives a concrete finite place to put residual mass in later tests.  The
guardrail is that the unital chart is naturally `F3^3`, not cyclic `C27`; the
cyclic carry has to be supplied by the LRC marking, not by the unital itself.

## Random Repo Niche

THM-417 is the useful nearby niche: the signed face is the symmetric closure of
the pair-sum sieve, and it explicitly separates the witness clock mod `14` from
the shell clock mod `27`.  That is exactly the distinction this unital lift
repeats geometrically.

## Connections

This comment connects HYP-2937 to the existing guardrails HYP-2891 and HYP-2894.
Those notes already warned that unital/design objects should be used after
observer-relative labels are fixed.  HYP-2937 makes the warning executable:
the design incidence is perfect, but the cyclic carry is external.

