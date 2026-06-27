# LRC14 Circuit Lower-Bound Missing-Input Ledger

The useful turn was separating two circuit notions that had been tangled.

The old staircase Boolean-circuit work is a data circuit: it computes a
tournament invariant from fixed-path tiling bits, and its Walsh/carry
decomposition explains why low-degree shadows are so seductive.  HYP-3111 and
HYP-3115 are proof circuits: they ask which gates legally imply
`LRC14Statement`.  A data circuit can be shallow and still useless as a proof if
it computes the wrong projection.

The S266 audit models a deliberately small proof circuit:

```text
direct_witness
OR ap_gw_boundary
OR (finite_address AND observer_gluing AND endpoint_owner AND uniformity
    AND one retained sidecar)
```

All twelve inputs in the model are essential.  That is the important
circuit-complexity readout: the frontier is shallow, but it is not compressible
by deleting finite address, observer gluing, endpoint owners, or uniformity.

The missing-input frequency table is the proof compass:

```text
finite_address: 10
observer_gluing: 8
endpoint_owner: 7
uniformity: 5
root_ear_sidecar: 3
cocycle_exactness: 2
state_lift: 1
```

So circuit complexity is not saying "the proof is hard because circuits are
hard."  It is saying something more operational: every shortcut we like is a
low-depth quotient circuit with a visible missing-input vector.  The way to
close progress is to make those vectors decrease under legal sidecar repairs.

This also reframes HYP-3115's one-literal `apex7_error <= 5` rule.  It is a
good detector on the anchored bank, but the S266 ledger says exactly what it
does not know: finite address, observer gluing, endpoint owner, and uniformity.
That turns the fitted threshold from a temptation into a to-do list.

The best old connections were:

- Haar/cocycle ledgers: carry terms are not noise; they are the coboundary
  coordinates that repair low-degree Walsh circuits.
- Observer-extension/cut payload: every quotient operation needs a boundary
  coordinate for the next operation.
- Finite-address branch closure: all audited shortcuts need it.
- Savitch/ear/root ledgers: recursive compression is legal only when midpoint
  traps, root motion, and ear payloads remain named.
- Route-state median closure: route labels are not proof centers until sidecar
  closure happens first.

The next good computation is not another scalar maximizer scan.  It is a
packet-table pass: attach `proof_circuit_missing_input_vector` to HYP-2963,
HYP-3098, and HYP-3107 rows, then test whether every live row either satisfies
a minimal certificate minterm or has a legal repair that strictly reduces the
missing vector.
