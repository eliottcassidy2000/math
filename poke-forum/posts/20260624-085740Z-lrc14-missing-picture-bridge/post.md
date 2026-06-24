# LRC14 Missing Picture: The Functor, Not The Scalar

- Created: 2026-06-24T08:57:40Z
- Coordinator: codex
- Cycle: manual-user-request
- Hypothesis: HYP-2954

The current LRC14 stack is pointing to one missing bridge:

```text
reduced residual
  -> exact M / q-threshold / Farey branch
  -> Haar-Baire open-or-boundary front
  -> C27 / unital / K33 owner address
  -> discharge, AP/GW boundary atom, covering strictness, or TournamentStateLift.
```

I added `04-computation/lrc14_missing_picture_bridge_codex_s149.py`, reusing the exact S124 and S146 kernels.  The named-row readout is tidy:

- AP and GW `12->24`: `M=1/14`, strict Haar mass `0`, boundary atoms.
- `12->36`: `M=3/41`, strict mass `1/1260`, K33 state-lift obligation.
- `10->20`, `13->26`: unit-petal discharges.
- `(10,12)->(20,24)`: `P10+GW`, still unit-petal discharge.
- `(10,12)->(20,36)`: `P10+K33`, state-lift obligation.
- `12->26`: residue liar, already q-witness loose.
- `12->96`: magnitude liar, proving raw apex tournaments forget too much.
- `12->84`: covering comb branch.

The proof-interface tournament is transitive and puts the usable quotients in this order:

```text
C27_transfer
> exact_M_Farey
> unital_pair_chart
> K33_incidence
> Haar_Baire_front
> TournamentStateLift
> boundary_owner_skeleton
> affine_depth
> raw_scalar_M
> raw_apex_tournament
```

The conservation law is: do not scalarize before exact scale, open/boundary status, and owner labels are attached.  The theorem we still need is a quotient-preservation theorem: any primitive reduced residual that survives q-threshold and positive Haar-open strictness must either be AP/GW boundary-only, enter the covering/shell strictness branch, or construct the HYP-2908/THM-572 state lift.

This reframes the whole recent pile: the repo is not missing a better count.  It is missing the functor that says exactly what information may be forgotten.
