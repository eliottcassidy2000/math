# LRC14 Lean Proof-Frontier Ledger, S259

The useful formalization move was not to assert another theorem.  It was to
make the current edge typed enough that the false shortcuts are structurally
excluded.

`TournamentH7.LRCProofFrontier` records HYP-3107/T1184/LTI-245/LTT-143:

```text
solved: q-witness, level-7 lift, pair-Pascal cap RHS plus THM-577
        symbolic dense cap values, Mreach readout
open:   coverage extremality, reflection-Perron, Node-3, finite-ruler glue
false:  coarse mod-14 winding-H as a direct bridge at k>=8
```

The exact new Lean checks are the cap arithmetic:

```text
cap_k = k(k+1)/182 for k=10..13
cap_9 = pair mass - 1/4004
cap_8 = pair mass - 1081/76440
```

After rebasing over THM-577, the `k=10,11` equalities should be read as
symbolic closed-form overlap values, not merely as table-backed dense cap
entries.  The remaining problem is not the rational value; it is the
coverage-extremality / packet-optimality route that turns those values into
proof-bearing residual exits.

The top-level Lean theorem stays conditional:

```text
residual classifier + residual-to-finite-address packets
  -> CuttingEdgeBranchCoverage
  -> LRC14Statement
```

After the second rebase, there is also a parallel closure route:

```text
residual classifier + residual-to-observer-gluing certificates
  -> ObserverGluingCoverage
  -> LRC14Statement
```

That makes `LRCObserverGluingLedger` the producer interface below this
frontier, and `LRCProofFrontier` the map that keeps all open producer
obligations visible.

The final upstream HYP-3099 rebase sharpened the tournament side: cap
optimality is bounded but the improvement tournament is non-transitive, so the
frontier should ask for a finite certificate over local minima instead of a
greedy exchange proof.  It also kills the tempting apex-7/H=7 shortcut; the
surviving theorem currency is order-2 antipodal structure plus retained
observer sidecars, not raw Hamiltonian-path numerology.

The final HYP-3100 rebase adds the contradiction-grammar side.  The
`LRCBleedingEdgeFrontier` wrapper is now the third conservative closure route:
finite-address packet, observer-gluing certificate, or explicit
`BleedingEdgeFrontierCoverage` with encoding/predicate/destroyed-coordinate
sidecars.

The final mainline fetch also added HYP-3101/HYP-3102/HYP-3103/HYP-3104 plus
HYP-3106 signals.  Those do not change the Lean theorem, but they sharpen its
open producer fields: coverage extremality now has a normal-fan/Cech
component route, a PGF-zero test signal, and a maximizer signal atlas, while
observer-gluing needs first-obstruction cocycles and HYP-3106 perspective
sidecars before any quotient is trusted.

There is a live namespace warning in that sentence: current mainline still
overloads `HYP-3101` between the normal-fan component file and the S31ah
tournament certificate-toolkit index entry.  The old HYP-3103 split has been
repaired by assigning miss-count PGF zeros to HYP-3103 and perspective
groupoid controlled forgetting to HYP-3106, so the frontier should cite route
names whenever historical ID meanings matter.

The S31ah toolkit addendum is also a useful negative result for this Lean
frontier.  It validates the H/Omega certificate engine and the `{7,21}` gaps,
but it agrees that the coarse LRC14 bridge is vacuous: the mod-14 obstruction
is order-2 antipodal matching, not an H=7/Omega-K3 event.  That makes H a
terminal certificate only after a fine-scale or packet-preserving encoding
has already carried the LRC predicate.

This is the right boundary.  The formal module knows the cap RHS is solved,
but refuses to collapse the LHS extremality, Node-3, fine-scale tournament, or
finite-ruler obligations into the same scalar.

The creative possibilities are now concrete.  The most promising is to replace
the degenerate coarse winding tournament by a fine mod-`p` or sector-pair
observable, then test whether that observable tracks coverage extremality
without losing magnitude or the HYP-3101 component packet.  The second is to
turn the `k=8,9` cap debt into an Eberlein/Hankel or HYP-3103/HYP-3104
PGF-zero/maximizer-signal certificate.  The third is to make the
S258/HYP-3098 observer-gluing ledger output actual `FiniteAddressBranchPacket`
records for representative residual rows, with HYP-3102 first-obstruction
status attached.
The fourth is to make HYP-3093/HYP-3097's triad into an experiment schema:
equinumerosity supplies the Pascal/base-count shadow, equidecomposability
supplies the sector-pair/component/endpoint-owner scissors fiber, and
equidistribution supplies the Haar/Weyl limit only after the resonance sidecar
has been retained.

The hidden lesson: formalization is functioning as a microscope.  It separates
what is solved, what is merely named, and what is a tempting but degenerate
quotient before those categories get mixed by prose.
