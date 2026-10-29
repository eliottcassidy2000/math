# Incoming work integrated at close-out

The proof checkpoint was rebased onto origin/main at4961bb2f35. The
following new work changes the next questions but does not change this
packet's operation menus or proofs.

## Read-only scans versus persistent attachment

The incoming [moving-selector theorem C1](../moving_selector_atoms_20261028/THEORY.md)
proves the exact minimum for a read-only cursor on Hamilton presentations:

    K_scan(T,i)=(i,(R_j(T),R_j(tau T))_(0<=j<n),B(T)),

where R_j is the tournament rooted at first, last and current position j.
The two rooted models plus their ordered vertex-orbit lists suffice.
This is a minimum for its stated observation/update language. Whether this
profile reconstructs every Hamilton presentation remains open there.

Its C2 constructs a different obstruction for ordered child runs:
individual cursor probes on C3[C3,1,1] cannot distinguish two orderings
that differ by a transposition within the inner triangle. Actual strong-child
assembly and path-dependent complement do distinguish them. C3 identifies
the component-order fiber over an aligned pair as k!/|G| and proves that
individual anchored component probes classify every ordering exactly when
their common pointwise stabilizer is trivial.

These are inherited results, not independently replayed by this packet.
They refine our [J6](PORTS.md): a cursor that cannot leave marks is a
different language from movable substitutions that persist. J1's register
and J3's strong tags supply exactly such persistent joint names.

The remaining smaller-carrier questions are therefore explicit:
faithfulness of K_scan, behavioral minimality under supplied-run assembly,
and fixed-selector or bounded-depth attachment menus. The arbitrary-position
persistent menu is already discrete by J2-J3.

## The missing closure operation can be multiplication

Incoming [Session92 joint-closure algebra](../session92_joint_order_carriers_20261029/algebra/JOINT_CLOSURE_ALGEBRA.md)
proves the standard finite-state correspondence between equitable
partitions and function algebras closed under transitions AND pointwise
products. Its native folded observer has two states whose linear
continuation signatures agree, while a shared-successor squared observation
gives3 versus9.

This supplies a concrete counterpart to J5. A spectrum is a trace-level
summary; a linear response space need not be an algebra; an algebra
containing every arc readout is already the full function space. These
are three distinct sizes, so "minimizing K" must specify whether K is a
partition state, a linear realization, or a faithful relation carrier.

The supplied-run and nonlinear folded examples have their own declared
interfaces. Neither should be silently transplanted to LRC or to ordinary
substitution without path-dependent complement.

## A separately implemented bounded replay of our controls

[Session92's source-phase audit](../session92_joint_order_carriers_20261029/phase_controls/SOURCE_PHASE_CONTROLS.md)
pins our first checkpoint b07ec51b20 and imports no producer modules.
At close-out we ran its checker normally and with python -O. Both gave
the same49 successful checks:

- all four selected-attachment Hamilton counts, each cross-checked by
  direct permutation enumeration as well as subset DP;
- the exact uniform length germ of both phase-hostile rows on
  0<=epsilon<=1/72;
- the speed7 continuation on the entire interval
  0<=epsilon<=1/84, using every rational endpoint-order cell.

The four counts are307/777 for the first core and293/981 for the second,
in exceptional-sink/C3 order. The common length germ is epsilon, and
speed7 survival at the original threshold is false/true.

Reproduction from the repository root:

    python -B 05-knowledge/results/session92_joint_order_carriers_20261029/phase_controls/check_source_phase_controls.py
    python -B -O 05-knowledge/results/session92_joint_order_carriers_20261029/phase_controls/check_source_phase_controls.py

This independently implemented audit checks the four counts and uniform
phase example. It does not audit J1-J3's all-size register proofs, L1-L3's
general germ/Gibbs results, or the full multivariate claim in L5.
Those retain the same-agent proof-review status declared in this packet.
