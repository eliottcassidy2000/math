# Even graphs as parity-dual addresses

The even-graph thread is now useful because it stopped pretending to be the
same object as the tournament obstruction.

The newest scar calculation says the key fact plainly.  Tournaments miss the
clique value `7` because the conflict graph packet `K_3` is not realizable:
three pairwise vertex-conflicting triangles force the directed pentagon layer.
Even graphs do the opposite.  `K_3` is degree-even, so the even-graph spectrum
realizes `7`; it misses the even clique values instead.  The projection into
even graphs therefore heals the tournament scar while preserving the
cycle-space address where that scar was born.

The incoming S38 cycle-half note fits this exactly.  With the base Hamiltonian
path fixed, the even graph is the tournament's cycle half; once that cut is
fixed, the cycle half carries `H`.  But the cycle-only quotient `E_n` forgets
the cut and no longer separates `H` values.  So the slogan is not "even graphs
are the obstruction."  It is "even graphs are the address before the quotient
that would erase the obstruction."

That is exactly the warning the LRC proof needs.

For LRC(14), the dangerous quotient chain is

```text
raw speed set
  -> exact-denominator phi packets
  -> support-six relation lattice
  -> mod-7 coimage
  -> chi_7 / affine signed kernel
  -> scalar margin.
```

If we project early, the live obstruction can disappear, just as tournament
`7` disappears in the even-graph shadow.  But the shadow is still valuable if
we read it as an address quotient.  It says to keep support incidence before
we ask for a scalar inequality.

The concrete LRC object should be a signed support-six packet graph.  Its
vertices are not runners.  Better candidate vertices are exact-period packets,
projective coimage classes, low-height wall labels, Fourier residue shells,
affine zero lanes, or proof obligations.  An edge means either compatibility
inside one low-height wall ledger, signed conflict under the HYP-2632
Legendre/affine kernel, or high additive-energy coupling.  The fingerprints to
watch are the ones tournament analysis already taught us to measure: score
histograms, SCCs, directed 3-cycles, chordless odd holes, and edge flips under
the relevant gauge.

The first finite packet graph is already exact.  For the repeated-residue
kernel from HYP-2632, take residue vertices `{0,2,3,4,5,6}`.  The negative
`4+2` loop weights are `-4,-25,-18,-25,-18,-18`; the positive `4+1+1` edge
weights are `0,1,8`, with the zero edges exactly on the affine matching
`a+b=2 mod 7`.  The new script
`04-computation/lrc14_repeated_packet_graph_codex_s101.py` verifies

```text
loop(a) + sum_b edge(a,b) = 0
```

for every vertex.  This is a local conservation law, not just a smaller signed
total.  The analytic lift should preserve this incidence balance through the
reciprocal hyperplane sums before any triangle inequality is applied.

This reframes the genuine-wide branch.  The target is no longer "control every
raw wide cluster."  It is:

```text
delete/certify finite low-height walls;
route relation-depth <= 2 to the finite atlas;
show relation-depth >= 3 has either a parity-null packet graph or enough
spread-cycle decorrelation for the margin.
```

Even graphs do not finish the proof.  They do identify the layer where a proof
should live: exact-period support-six incidence before scalarization.  The next
best computation is to lift the repeated-packet balance to the actual
reciprocal tails, then extend the packet graph from repeated residues to the
full HYP-2617 coimage atlas and compare odd-hole/edge-flip profiles with the
E7 metagraph, not with raw Hamiltonian-path counts.
