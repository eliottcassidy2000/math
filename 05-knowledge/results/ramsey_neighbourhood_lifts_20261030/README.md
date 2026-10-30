# Exact Ramsey capacities from joint incidence and fractional colorings

2026-10-30. **PROVED scoped theorems / CITED structural inputs /
FINITE-EXACT controls / VERIFIED external certificate replay.**
Same-agent proof review with separate computational algorithms; no
literature-priority or formalization claim.

The main advance closes both open extension targets of the
[previous Ramsey packet](../ramsey_circuit_interfaces_20261030/README.md):
for EVERY U5-free tournament, the exact maximum order at transitive
budget k is

    G(k)=sum_(0<=j<k) 2^popcount(j)=A006046(k).

A uniform circular rounding rule replaces the earlier T5/T7/T9 finite
chart banks. At most two half-circle-isolable marked points survive,
and they are adjacent; this proves the rule for every odd carousel.
Liu's published prime structure theorem then extends the exact integer
bound to the whole U5-free class. The Pascal-cell tournaments attain
it at every k. This refines a classical exponent, rather than claiming
a new unrestricted Ramsey bound. [C1-C3 proof](CIRCULAR.md).

A second all-size theorem turns a proper b-fold a-coloring into an
integer-weight packing inequality. At ratio5/2 it yields a function w
with w(2h)=(5/2)w(h), w(2h+1)=(3/2)w(h)+w(h+1). For the explicitly
defined nonempty C5/join/union graph class D,

    F_D(A,B)<=min(w(A)B,w(B)A),
    F_D(A,B)=w(A)B if 2^ceil(log_2 A) divides B.

Thus every clique budget has an exact asymptotic slope and infinitely
many exact realizations. The complete first four capacity rows are also
proved, including F_D(3,B)=4B-[4 does not divide B]. The scalar bound
is integral even when unattainable: actual child realization retains an
extra two-level parity guard. [G1-G4 proof](GRAPH_CAPACITY.md).

The supplied Dauth paper supplies an independently useful direction:
common-neighborhood clique signatures and exact flag inequalities.
We replayed its author's17985 hand-arithmetic checks and exact flag
certificates giving R(3,3,6)<=86, R(3,3,7)<=125 and R(3,3,8)<=174.
The preprint's stronger173 and other new bounds are cited, not replayed.
These improvements belong to the author. A small derived defect lemma
prices total unsaturation and forces actual42-edge port types; those
same types decisively block a tempting shorter argument.
[N1-N3 and source audit](NEIGHBOURHOODS.md).

## Verification and scope

- All87376 charts through15 vertices,1250 larger weighted queries,
  and native maximal-support checks through11.
- All33867 labelled tournaments through6 were examined; the19107
  U5-free survivors satisfy the bound. Literal attainers through k=16
  and the full81/27 hostile are checked.
- The exact graph decomposition DP tests311604 candidates;32 literal
  graph witnesses independently recover clique/independence numbers.
  All130388 feasible tested C5 weight vectors obey G1.
- A separate auditor checks the classical maximum recurrence and actual
  recursive realization through clique budget256, without importing the
  producer. Normal and optimized runs agree.
- Incoming Session120's two231-vertex spectral minimizers have their
  exact fourth moments replayed and explicit U5 witnesses. Spectral
  optimality does not silently satisfy the forbidden-pattern hypothesis.

[AUDIT.md](AUDIT.md) gives precise universes and commands;
[RECOVERY.md](RECOVERY.md) retains provenance, portfolio, connection
maps, incoming work and stopping reasons. The source manifest is an
integrity record, not an additional proof.

## Next decisive targets

1. Determine the exact residue law for F_D(A,B) beyond the first four
   clique budgets; G2 already fixes every asymptotic slope.
2. Replace the capacity recurrence after admitting Paley7, whose literal
   transitive supports refute the U5-free law. A failed rounding chart
   alone is not a capacity counterexample.
3. Turn common-neighborhood signatures of actual critical ports into
   new certified exclusions, preserving roots, overlap and sampling
   denominators. The new preprint supersedes the old91 search boundary.

General complement/re-marking closure, A000568 enumeration, LRC(14),
Collatz and unrestricted circuit lower bounds remain open. No count,
spectral statistic or capacity query is promoted to a full carrier.
