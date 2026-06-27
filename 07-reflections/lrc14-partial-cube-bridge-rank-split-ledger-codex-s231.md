# LRC14 Partial-Cube Bridge-Rank Split Ledger - S231

This is a bridge-rank and two-lane split addendum to
HYP-3063/T1145/LTI-210/LTT-108.  The conceptual carrier is S227 and the broad
finite exact scout is S228; this S231 ledger adds the `K_{k,k+1}` rank/debt
reading and the explicit Moser split `x=a+2b`.

The useful move in this pass is to stop treating Moser-de Bruijn and
fibbinary as just automaton labels.

Fibbinary words are independent sets of a path.  The resulting Fibonacci cube
is a median partial cube: its coordinates are real cuts, but adjacent cuts are
not allowed to fire together.  Moser-de Bruijn words are even-bit base-4
digits, so an m-digit prefix is just `Q_m`; two Moser lanes split every binary
coordinate through `x=a+2b`.

That makes both languages controlled-forgetting objects.  If an LRC quotient
forgets a fibbinary or Moser coordinate, the missing thing is a partial-cube
Theta class or carry boundary, not a vague sequence membership bit.

The doubled triangular row

```text
2, 6, 12, 20, 30, 42
```

is `k(k+1)`.  It is twice the number of edges of the k-simplex, the number of
directed simplex edges, and the number of lines in the `K_{k,k+1}` bridge from
the diagonal-layer model.  S217 then says the bridge has only `2k` cut rank,
so the surplus `k(k-1)` is rectangle-cycle sidecar debt.

This gives a clean synthesis:

```text
Moser/fibbinary = partial-cube cut sidecars
k(k+1)          = ordered simplex/bridge sector sidecar
rectangle debt  = cycle-space residue sidecar
```

For LRC14, the sidecar fields to add are:

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

Assumption challenge: I considered runners, raw sequence names, automaton
states, partial-cube hyperplanes, simplex directed edges, bridge lines,
rectangle cycles, endpoint owners, exact packet routes, and proof obligations.
The selected vertices are proof carriers because raw sequence membership does
not preserve LRC boundary/open status or theorem route.

The next useful computation is a packet-fiber audit: inside the HYP-2963 rows
already tagged by Moser/fibbinary automatic words, attach partial-cube cut
fields and see whether mixed route fibers split before exact magnitude is
reattached.
