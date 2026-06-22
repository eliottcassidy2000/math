---
id: HYP-2882
status: SYNTHESIS / proof-target; supported by HYP-+2880, HYP-2878, HYP-2879, HYP-2628, HYP-2630, HYP-2632
source: codex-2026-06-22-S100
tags: [tournaments, lrc14, totient, euler-copy, odd-holes, odd-cliques, strong-atoms, coimage, character-kernel, tournament-analysis]
related:
  - HYP-+2880
  - HYP-2878
  - HYP-2879
  - HYP-2880
  - HYP-2881
  - HYP-2877
  - HYP-2628
  - HYP-2630
  - HYP-2632
  - HYP-2523
  - THM-029
  - THM-079
  - THM-200
---

# HYP-2882: primitive packets mediate the odd clique / odd hole and totient bridges

The common structure behind the recent tournament and LRC threads is not the
number `7` by itself.  It is the order in which primitive packet information is
retained before quotienting.

```text
primitive packets
  -> quotient / coimage / conflict graph
  -> signed or forbidden boundary
  -> scalar value.
```

When that order is respected, the `H=7` obstruction, the `E_7` odd holes, the
strong-ear boundary values `49,75`, and the LRC exact-period/totient packet
machinery line up.  When it is violated, the key obstruction disappears into a
scalar count.

## Tournament side

HYP-+2880 and the KPS S31g recursion note sharpen the `H=7` obstruction.

For the odd-cycle conflict graph `Omega(T)`,

```text
H(T) = I(Omega(T), 2).
```

For a clique `K_r`,

```text
I(K_r,2) = 1 + 2r.
```

The small connected graph preimages of `I(G,2)` show that `H=7` has a unique
connected preimage: `K_3`.  This is not a generic odd-clique law.  It is a
bottom-level uniqueness accident: `3,5,7,9` are forced clique values, and only
the `K_3` preimage is forbidden as a conflict graph.  From `H>=11`, non-clique
preimages already exist.

The "odd clique / odd hole duality" should therefore be phrased carefully:
`K_3` is not a forbidden induced subgraph in the Strong Perfect Graph Theorem
sense.  It is the unique conflict-graph preimage of the scalar value `7`.  The
dual odd hole appears because three pairwise vertex-conflicting triangles force
a directed `C_5`; HYP-+2880 verifies the cycle-space identity

```text
C_5 = triangle_1 xor triangle_2 xor triangle_3
```

with the three triangles pairwise vertex-conflicting.  Thus the tournament
primitive packet is:

```text
three conflicting primitive odd cycles
  -> conflict clique K_3
  -> forced C_5 hole in cycle space
  -> scalar H=7 is skipped.
```

HYP-2878 supplies the metagraph side: `E_7` first gains odd holes, specifically
`1496` chordless `C_5` and `196=14^2` chordless `C_7` holes.  The established
literal bridge is the cycle-space `C_5 = K_3` obstruction above; the still-open
bridge is whether the chordless `C_5` holes in the `E_7` iso-class adjacency
are the same primitive packet after the metagraph quotient.

## Totient side

The user's copy-rule reframe from HYP-2628 is exactly the primitive packet law:

```text
sum_{d|n} c(d) = n  =>  c = mu * id = phi.
```

On a `q`-grid, `phi(d)` is not decorative arithmetic.  It is the number of
residues of exact denominator `d`.  Therefore the LRC proof quotient must be

```text
raw denominator
  -> exact-period phi packets
  -> squarefree mask / unit seam
  -> mod-7 coimage
  -> chi_7 / affine signed kernel.
```

This is exactly what HYP-2630 and HYP-2632 found.  Euler-copy capacity is
uniform on `F_7^*`, so it cannot separate QR from NQR.  The split lives one
layer later in the signed character phase:

```text
2*S_9(1,1,1,1,a,a)/U = -43 - 7*chi_7(a),
```

and in the affine zero lane

```text
a+b = 2 mod 7.
```

The totient packets are therefore analogous to the primitive odd cycles in the
tournament story.  They must be kept until the character/clique-hole boundary
has been evaluated.

## Strong atoms and ears

HYP-2877 and HYP-2879 make the same point inside the `H` spectrum.  The scalar
`H` is multiplicative over strong components, so strong components are the
primitive packets of the Hamiltonian-path count.  Ear insertion refines the
packet law:

```text
H(T+x) = start(sig=1) + end(sig=0) + Q(sig=0,sig=1).
```

Balanced ears cover almost all `n=8` strong values, but miss exactly `49,75`
at the non-isomorphic strong-ear quotient; boundary ears fill them.  HYP-2880
then shows why raw rooted cut size is too coarse: rooted `w=3` can still
realize `49,75`, but only the insertion profile sees the near-doubling defect

```text
49 = 2*25 - 1,
75 = 2*39 - 3.
```

This is the same quotient-order lesson:

```text
strong parent atom
  -> labelled ear insertion profile
  -> strong-ear quotient
  -> scalar H value.
```

## LRC proof reframe

The missing LRC piece may be a primitive-packet theorem rather than another
scalar denominator or tournament-class census.

Candidate theorem shape:

```text
After exact-period phi expansion and wall deletion, every persistent
non-lonely LRC(14) covering packet either
  (a) collapses to the AP boundary D=14, or
  (b) contains a K_3/C_5-type primitive incompatibility whose signed
      chi_7 + affine kernel has positive witness mass.
```

The "one odd clique and odd hole" language then becomes useful, but only after
translation:

- odd clique: a small set of primitive obstruction packets that pairwise
  conflict, like `K_3` in `Omega(T)`;
- odd hole: the forced cycle-space boundary that prevents the obstruction from
  staying scalar, like `C_5`;
- totient packets: the exact-period primitive copies that prevent the LRC
  denominator quotient from deleting the boundary before the `chi_7` phase is
  computed.

This suggests a concrete next test: build the support-six repeated-residue
conflict graph whose vertices are exact-period phi packets or additive
frequency shells, and add an edge when two packets cannot both remain in the
same signed covering certificate.  Then check whether the affine zero lane
`a+b=2`, the Legendre selector `Q(a,b)`, and the `E_7` `C_5/C_7` holes are
visible as odd hole / odd clique fingerprints of that packet graph.

## Assumption challenge

Candidate vertices considered: runners, tournament vertices, free arcs,
directed triangles, odd cycles, `Omega` vertices, `E_7` even-graph classes,
strong components, ears, exact-denominator packets, divisor masks, unit
residues, additive Fourier shells, quadratic characters, affine zero lanes,
and proof obligations.

Chosen quotient: primitive packet before scalar quotient.  This preserves the
information needed by the LRC and tournament predicates: conflict incidence,
exact-period multiplicity, `chi_7` phase, and ear insertion profiles.  It
destroys raw row identity, raw runner identity, and raw divisor counts.  The
challenged assumption is that the right bridge is a direct numeric analogy
between `phi`, `H`, and `7`; the stronger bridge is categorical/operational:
all three are primitive-packet expansions followed by a dangerous quotient.

## Next tests

1. Lift HYP-+2880 from the cycle-space `C_5` identity to the actual `E_7`
   metagraph holes: label the `1496` chordless `C_5` holes by their underlying
   even graph and test for the triangle-fan `K_3` packet.
2. Build the exact-period support-six packet conflict graph from HYP-2632 and
   test whether its odd holes are the affine zero lane and Legendre selector in
   graph form.
3. Compare HYP-2879's exposed-slot matrix `Q_T` to `E_7` holes: look for the
   same `K_3 -> C_5` obstruction in the ear insertion graph rather than in raw
   fixed-path rows.
4. Formalize the quotient-order lemma as a reusable guardrail: do not project
   `raw denominator -> squarefree` or `rooted row -> cut weight` before the
   primitive packet incidence has been evaluated.
