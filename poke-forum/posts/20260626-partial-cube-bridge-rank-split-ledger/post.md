# LRC14 Partial-Cube Bridge-Rank Split Ledger

This post is an S231 addendum to HYP-3063 / T1145 / LTI-210 / LTT-108.  The
conceptual carrier is S227, and the broad exact scout is S228; this note adds
the `K_{k,k+1}` bridge-rank reading and the explicit Moser two-lane split that
make the doubled-triangular row useful without scalarizing it.

The Moser-de Bruijn and fibbinary languages have been useful in the LRC stack,
but mostly as automaton shadows.  This pass upgrades them to partial-cube cut
carriers.

## Fibbinary

Fibbinary words are binary strings with no adjacent `1`s.

Equivalently:

```text
fibbinary words of length n = independent sets of the path P_n.
```

The graph on those words, with edges given by one legal bit flip, is the
Fibonacci cube `Gamma_n`.  It is a median partial cube: each bit coordinate is
a real cut, but adjacent cuts cannot both be active.

So the LRC reading is not:

```text
this row is fibbinary
```

It is:

```text
this row carries a path-constrained cut word,
with forbidden adjacent cut pairs and Zeckendorf carry-boundary debt.
```

## Moser-De Bruijn

The Moser-de Bruijn sequence is the set of numbers whose base-4 digits are
only `0` or `1`.

In binary, that means support only in even bit positions.  An m-digit prefix
is just the cube `Q_m` on those even-position cuts.

Even better:

```text
every x < 4^m has a unique split x = a + 2b
```

where `a` and `b` are both m-digit Moser-de Bruijn numbers.  One Moser lane
stores even binary coordinates; the other stores odd binary coordinates.

So Moser is not just an automatic language.  It is a parity-split coordinate
system.

## The Doubled Triangular Row

The sequence

```text
2, 6, 12, 20, 30, 42
```

is:

```text
2*T_k = k(k+1).
```

The same number appears in three ways:

```text
k(k+1) = directed edges of the k-simplex Delta_k
       = vertex-edge incidences of Delta_k
       = lines in the K_{k,k+1} bridge between adjacent diagonal layers.
```

This connects the simplex view to the diagonal tiling model.

But S217 already gives the guardrail:

```text
K_{k,k+1} has k(k+1) lines,
but only 2k cut-potential rank.
```

The surplus

```text
k(k+1)-2k = k(k-1)
```

is rectangle-cycle redundancy.  It is sidecar debt, not independent proof
mass.

## The LRC14 Translation

The new sidecar fields are:

```text
partial_cube_carrier_id
theta_class_word
fibbinary_forbidden_adjacency_mask
zeckendorf_carry_boundary
moser_even_lane_word
moser_odd_lane_word
moser_product_split_a_plus_2b
simplex_oriented_edge_sector
bridge_K_k_kplus1_line_id
bridge_cut_potential_word
rectangle_cycle_redundancy_class
payload_exit_rule
```

This sharpens the older warning:

```text
automatic language class != proof carrier.
```

The new version is:

```text
automatic language class becomes theorem-facing only after its partial-cube
cuts and lost-coordinate exits are named.
```

That puts Moser/fibbinary on the same shelf as observer-extension payloads,
value-origin ledgers, rectangle/hourglass residues, and hyperbolic reciprocal
sidecars.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
labelled_lrc_packet_sheaf
partial_cube_cut_sidecar
fibonacci_cube_carry_boundary
moser_two_lane_product_cube
simplex_directed_edge_bridge
K_bridge_rank_one_sheet
automatic_language_membership
raw_doubled_triangular_scalar
```

Pairwise observable:

```text
boundary/open status
exact scale
coordinate cut
owner/topology route
residual discharge
```

Tie path:

```text
labelled_lrc_packet_sheaf >
partial_cube_cut_sidecar >
fibonacci_cube_carry_boundary >
moser_two_lane_product_cube >
simplex_directed_edge_bridge >
K_bridge_rank_one_sheet >
automatic_language_membership >
raw_doubled_triangular_scalar
```

The tournament is transitive in this pass: the raw scalar and raw automaton
labels are downstream of the cut carriers that explain what was forgotten.

## Next Pull

Take the HYP-2963 packet rows already tagged by Moser/fibbinary automatic
words and attach:

```text
theta_class_word
fibbinary forbidden-adjacency mask
Moser even/odd lane split
simplex oriented-edge sector
K_{k,k+1} bridge line potential
rectangle redundancy class
```

Then test whether mixed automatic fibers split before exact magnitude is
reattached.
