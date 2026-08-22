# Private sheets, bidirected tournaments, and the marker polynomial

**Research reflection / provenance, not a truth source.**  Exact claims are
routed to audit-pending THM-3473 and proved THM-3469.  LRC(14), the D5 map,
and the `7 tensor 13` spectral closure remain open.

## Inheritance pass

- **Closest proved mechanism:** THM-3469's nearest-`p` chart, with six affine
  core channels, the `2p` order-three backbone, and the `3p-6` septimal
  repair.
- **Canonical hostile:** when `k==2 (mod 3)`, the global rank is four although
  every owner in the displayed eight-owner presentation has private sheets.
- **Corrected near miss:** a pairwise tournament cannot see singleton support,
  so it cannot certify irredundancy.
- **Least-used sidecar:** the complete Boolean support vector on each sheet,
  rather than its union bit or pairwise coactivity graph.

The resulting theorem candidate is not merely “eight private witnesses.”  It
classifies every sheet support and finds only eleven of the `2^8` possible
Boolean types.

## The whole family compresses to one sparse polynomial

Let `z_i` mark owner `i`.  The complete support enumerator is

```text
H_k=sum_i p_i z_i
    +2k z_1z_4z_5z_6
    +(2k-1)z_2z_3z_5z_8
    +(2k+2e_2)z_5z_7,
```

where `(p_1,...,p_8)` is the exact private-count vector.  This is the most
economical carrier found in this lane:

- its degree-one part is the eight-owner deletion certificate;
- its degree-two monomial is the repair overlap;
- its two degree-four monomials are the backbone collisions;
- it has no degree zero, three, five, six, seven, or eight term.

Under one common marker it becomes

```text
(36k-2-2e_2)t+(2k+2e_2)t^2+(4k-1)t^4.
```

The limiting degree-one mass is `6/7`; each higher support packet has mass
`1/21`.  Thus the septimal danger threshold returns as the exact limiting
overlap mass `1/7`.

## The user's generalized tournament is exact here

Declare two owners related when they co-occur on at least one sheet.  The
observable is symmetric, so encode it by both arcs; unrelated pairs are
missing edges.  The graph is exactly

```text
K4(1,4,5,6) union K4(2,3,5,8) union K2(5,7).
```

It is two bidirected tetrahedra plus one bidirected edge, glued at the `2p`
backbone vertex.  This realizes the user's “tournament of size four with
missing and both-way edges” literally, after changing “tournament” to the
typed bidirected two-section.

The quotient contract is sharp:

```text
vertices:       the eight labelled owners,
observable:     existence of a common active sheet,
orientation:    both directions; no arbitrary gauge,
ties/missing:   both-way edges or no edge,
preserved:      pairwise coactivity existence,
lost:           sheet counts, singleton atoms, k mod 3, ancestry,
needed sidecar: the eleven support coefficients of H_k,
target:         presentation irredundancy.
```

The graph alone fails the target because all private supports are isolated
monomials and therefore invisible to edges.

## A ternary state on the Berggren/Fibonacci branch

The only count anomalies depend on `k mod 3`.  They migrate cyclically:

```text
0: p+1 loses two;
1: p-1 loses two;
2: 2p and 3p-6 each lose two.
```

On THM-3469's U-spine lanes `t mod 21 in {5,8,12,15}`, this ternary state has
minimal ambient period `63`.  Each colour occurs four times per period, so
each is a subset of the harmonic series with coefficient `4/63` ambient and
`1/3` conditioned on the lane.

This is a residue colouring of a distinguished Berggren branch, not a claim
that the Berggren tree itself has become ternary for this reason.  The tree
ancestry chooses the quadratic labels; the three-state automaton decorates
the affine cover presentation on the selected lanes.

## Incoming Rule 30 and factorial work suggest the next operator

THM-3471's three forced Rule 30 strips cancel after unmarked scalarization,
while a transverse slack marker detects the circuit at first order.  THM-3466
expresses factorial face descent by commuting deletion operators
`product(1-partial_i)`.  THM-3473's first homogeneous marker part detects
private sheets that the pairwise quotient loses.

The common grammar is:

```text
labelled object -> scalar quotient -> invisible cancellation
               + transverse/deletion marker -> recovered first response.
```

This does not yet give the D5 map.  A lawful bridge would start from the
multivariate owner marker, attach endpoint labels to every support monomial,
and prove that one deletion derivative maps to a boundary-current coboundary.
The cheapest hostile is an eight-owner cover with positive `p_i` but zero
endpoint current at the common half centre.  If that hostile wins, the
missing coordinate is not Boolean ancestry but endpoint orientation.

## Updated concept board

1. **Private Boolean atlas:** independently audit the eleven support types and
   all residue-boundary counts.
2. **All-modulus rank transport:** re-audit repaired THM-3472; numerical grade
   and presentation ancestry must remain separate.
3. **Marker-current bridge:** enrich `H_k` by endpoint/current coefficients
   and test its first derivatives on hostile zero-current rows.
4. **Bidirected quotient:** determine the smallest sidecar reconstructing
   `H_k`; the graph plus `k mod 3` is still insufficient without coefficients.
5. **Spectral closure:** ask whether a weighted support polynomial can feed a
   nonzero `7 tensor 13` contraction, with a direct zero-current control.

The main conceptual gain is that the family is now an algebraic object, not a
list: a degree-graded sparse support polynomial whose pairwise shadow is two
tetrahedra and whose residue sidecar is ternary.
