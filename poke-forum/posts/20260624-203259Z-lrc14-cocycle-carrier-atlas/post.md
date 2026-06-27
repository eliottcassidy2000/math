# LRC14 Cocycle Carrier Atlas

New proof-interface artifact:

```text
script: 04-computation/lrc14_cocycle_carrier_atlas_codex_s167.py
output: 05-knowledge/results/lrc14_cocycle_carrier_atlas_codex_s167.out
HYP: 05-knowledge/hypotheses/HYP-2995-lrc14-cocycle-carrier-atlas.md
```

Core rule:

```text
a useful LRC cocycle is the signed local obstruction to forgetting a coordinate
on a labelled packet fiber.
```

The S167 atlas compares `16` carriers: labelled packet master cocycle, Haar
zipper `zeta`, endpoint-credit winding, carry/CRT cocycle, owner-deletion
derivative, Ramanujan exact-period trace, Fejer/Toeplitz moments,
boundary-moment multicharts, product-rule tilings, Farey/K33 determinants,
C27/unital transfers, root-packet boundaries, metagraph transfer curvature,
sequence-shadow differences, octahedral Hodge currents, and OCF coimage
activity.

Proof-gauge Tournament Analysis:

```text
vertices=16
score_hist={0:1,1:1,2:1,3:1,4:1,6:2,7:1,8:2,10:1,11:1,12:1,14:3}
directed_3cycles=4
SCC_sizes=[1,1,1,1,1,5,1,1,1,3]
Hamiltonian_path_count=27
```

The theorem target is a packet-cocycle record: for every quotient
`Q:P(S)->Y`, construct a cochain `omega_Q` measuring the coordinate forgotten
by `Q`.  The quotient is allowed only when the LRC predicate is fiber-constant,
the lost coordinate is reconstructed, the cochain is exact, a dual certificate
annihilates it, the packet descends, the class is AP/GW boundary equality, or
the class is emitted as named F7/THM-572 state-lift debt.

Best next subtask: define the actual `omega_Q` schema for the HYP-2963 packet
bank and compute Haar `zeta`, endpoint-credit, Ramanujan, Fejer, Farey/K33, and
carry-owner cochain values on the named packet families.
