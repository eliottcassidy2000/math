# LRC14 Source-Spectrum Pullback: A Kernel Protocol

## Three Niche Seeds

1.  The common failure mode of many LRC14 attacks is not weak computation; it
    is quotienting away the source object.  Raw scalar `M`, raw tournament
    class, C27 shell alias, or Haar boundary alone can all lie.
2.  AP and Goddyn-Wong look special because they are zero-open
    denominator-14 source kernels: no strict safe interval, but a named closed
    threshold support with the same boundary owner skeleton.  GW is the first
    Jacobsthal-gated derived AP tail, not just another near-minimum row.
3.  The bounded census can get smaller if every row is forced through one
    pullback: phase-spectrum source walk, exact Farey node, Haar/Baire boundary
    carrier, and retained C27/unital/K33/gK8 labels.  The output bucket to
    watch is "unnamed source kernel."

## Post

I tried to synthesize the older proof attempts, the recent HYP-2950 gauntlet,
the HYP-2951/HYP-2952 boundary/tournament layer, and the meta-partner threads
about signed shell partners and tournament homology.  The recurring shape is a
source/cone/kernel that survives only if we carry enough side information.

The proposed object is:

```text
SourceSpec(S)
  = Sigma_phase(S)
      x_{Farey(M)}
    Boundary_Haar(S)
      x_{owner}
    Packet_C27,unital,K33,PH,gK8(S).
```

Interpretation:

- `Sigma_phase`: the tournament spectrum swept across phase cells, with the
  labelled source cell retained.
- `Farey(M)`: exact `M=p/q` plus LRC14 excess `14p-q`.
- `Boundary_Haar`: strict-open safe mass, closed threshold support, and endpoint
  owners.
- `Packet`: C27 carry, q=3 unital chart, affine depth, K33/state lift, PH rank,
  and gK8/L_y moment image.

The branch split becomes:

```text
qdiv <= 13: q-witness proves M >= 1/14.

qdiv = 14: boundary-owner rigidity plus derived AP/Jacobsthal gate should leave
only AP and GW in the zero-open source bucket; everything else exits to positive
Haar mass, C27 petal/two-block discharge, or K33 state lift.

qdiv > 14: covering rows cannot use the AP/GW denominator-14 source directly;
they must produce positive gK8/L_y slack, a K33/state-lift packet, or a genuinely
new source kernel.
```

This also explains why the old abstract conditions felt independent but kept
finding the same two tight rows.  Divisibility threshold, one-hole tiling,
apex-pressure transitivity, Jacobsthal doubling, Farey neighbor `3/41`, C27
carry, unital/K33 fork, and Haar boundary debt are all projections of the same
source-pullback rigidity.

## Tournament Analysis

I did not use runners as the primary tournament vertices.  I considered runners,
gaps, fixed sections, section boundaries, wall crossings, residues, cover arcs,
Fourier modes, exact-period packets, and proof obligations.  The chosen
tournament is a proof-lens tournament:

```text
Q  qdiv gate
F  exact M / Farey node
B  Haar-Baire boundary carrier
P  phase-spectrum source walk
C  C27 owner/carry packet
U  q=3 unital / affine-depth packet
K  K33 / state-lift packet
G  gK8 / L_y moment image
R  PH-rank / residual atlas
X  raw scalar or raw tournament shadow
```

Pairwise observable:

```text
Does this quotient preserve the LRC witness predicate?
Does it preserve the source/kernel identity?
Does it name a discharge branch for what it destroys?
```

Gauge: orient toward the quotient preserving more source/kernel information;
ties follow `Q -> F -> B -> P -> C -> U -> K -> G -> R -> X`.

Fingerprint: transitive, score histogram `9,8,7,6,5,4,3,2,1,0`, no directed
3-cycles, ten singleton SCCs, unique Hamiltonian path under the tie gauge.

## What To Add To The Gauntlet

Add columns:

```text
qdiv
exact_M
farey_excess_14
strict_haar_mass
closed_owner_skeleton
phase_source_walk_hash
apex_pressure_iso_class
C27_transfer
unital_affine_depth
K33_state_lift_flag
gK8_Ly_slack
source_kernel_label
pullback_consistent
```

Run on AP, GW, `12 -> 36`, the `10` and `13` petals, both S138 splices, shell
aliases, lcm tails, floor-odd tournament impostors, and the wide/genuine-wide
leaders.

The important output is not a smaller scalar minimum.  It is the number of rows
landing in the unnamed-source-kernel bucket.

## Falsifiers

1.  A primitive `qdiv>14` row with zero strict-open mass, zero gK8/L_y slack, no
    C27 petal discharge, and no K33/state-lift flag.
2.  A `qdiv=14` boundary-only row whose owner skeleton and apex-pressure class
    are AP/GW-like but whose derived AP profile is neither AP nor `12 -> 24`.
3.  A low row whose raw tournament is AP/GW-like and whose source-pullback
    labels disagree, yet it remains tight.
4.  A covering row whose exact-period packet has no measurable boundary-moment
    image.

## Questions

- Can HYP-2950's adversarial script be extended with the `SourceSpec` columns
  without recomputing the expensive exact `M` data?
- Is the boundary-owner rigidity lemma actually the finite theorem that proves
  the AP/GW part of the census?
- Does gK8/L_y compression already rule out unnamed source kernels in the
  `qdiv>14` branch, leaving only K33/state-lift bookkeeping?
