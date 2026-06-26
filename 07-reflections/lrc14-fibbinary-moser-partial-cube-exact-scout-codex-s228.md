# LRC14 Fibbinary/Moser Partial-Cube Sector Carrier

codex-2026-06-26-S228

The useful merge is not "Moser plus fibbinary proves something."  The useful
merge is a type discipline for three nearby objects that the repo has been
touching separately:

```text
fibbinary      = path independence = Fibonacci cube = partial cube
Moser          = even-bit Boolean cube = simplex face lattice on even slots
2*T_{k-1} row  = ordered 1-face sectors of the (k-1)-simplex
```

The exact scout `04-computation/lrc14_fibbinary_moser_partial_cube_codex_s228.py`
confirms the count origins.  Fibbinary length `n` has layers
`C(n-r+1,r)`:

```text
Gamma_5: vertices 13, edges 20, layers [1,5,6,1]
Gamma_6: vertices 21, edges 38, layers [1,6,10,4]
```

Moser with `m` base-4 digits has binomial/simplex layers:

```text
M_5: vertices 32, cube_edges 80, layers [1,5,10,10,5,1]
M_6: vertices 64, cube_edges 192, layers [1,6,15,20,15,6,1]
```

And the user's row is exactly:

```text
k=2..7: k(k-1)=2*C(k,2)=2*T_{k-1}=2,6,12,20,30,42.
```

That last row is the ordered-pair / directed-edge sector row.  It belongs next
to HYP-3057/HYP-3059's value-origin and observer-cut ledgers, not as a new
sequence shortcut.

## Main Guardrail

There are tempting collisions:

```text
20 = Gamma_5 edge count = ordered sectors of K_5
12 = M_3 cube edge count = ordered sectors of K_4
```

These are not the same object.  The value is usable only after the origin is
typed:

```text
value_origin_type:
  fibonacci_cube_edges
  moser_even_bit_cube_edges
  ordered_edge_sector_count
  simplex_face_layer
  path_independence_layer
```

This is the same controlled-forgetting pattern as the observer-extension work:
counts are allowed to line up, but only the retained sidecar says what the
count can do in a proof.

## LRC Import

The sidecar I would add to the HYP-2963 stack is:

```text
partial_cube_carrier
language_id
native_transition
bit_phase
independence_complex_layer
simplex_face_layer
ordered_edge_sector_count
exact_M
endpoint_owner
route_label
safe_component_certificate
```

The key transition split:

```text
fibbinary closes under x -> 2x
Moser closes under x -> 4x
Moser does not close under x -> 2x unless bit phase is retained
```

So "automatic sparse language" is too coarse.  A dyadic LRC branch should not
use Moser membership without phase; a base-4 branch should not replace Moser by
generic fibbinary membership without recording that it changed transition.

## Assumption Challenge

I considered vertices as runners, speed gaps, fixed circle sections, residues,
cover arcs, automaton states, partial-cube coordinates, simplex faces,
ordered-pair sectors, and proof obligations.  The selected vertex set is proof
carriers:

```text
labelled_lrc_packet_sheaf
partial_cube_carry_state
fibbinary_fibonacci_cube
moser_even_bit_cube
simplex_face_layer
ordered_edge_sector_pronic
raw_sequence_membership
```

The LRC predicate preserved is not sequence membership.  It is the ability to
retain exact scale, endpoint/topology, route, automaton transition, and sector
origin until the quotient is safe.  Raw sequence membership destroys exactly
those coordinates.

The packet-priority tournament is transitive:

```text
labelled_lrc_packet_sheaf
> partial_cube_carry_state
> fibbinary_fibonacci_cube
> moser_even_bit_cube
> simplex_face_layer
> ordered_edge_sector_pronic
> raw_sequence_membership
```

The fieldwise warning remains cyclic in spirit: fibbinary, Moser, and simplex
views preserve different transitions, so none should be collapsed into a
single "sparse/automatic" scalar.
