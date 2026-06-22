# Primitive packets before quotients -- Codex S100

The useful synthesis from the tournament archive and the totient-copy work is a
quotient-order rule:

```text
primitive packets first, scalar quotients last.
```

On the tournament side, `H=7` is special because `I(G,2)=7` has the unique
connected preimage `K_3`.  That `K_3` is not a forbidden odd clique in the
perfect-graph theorem sense; it is the unique small conflict-packet preimage of
the scalar value.  The obstruction becomes geometric only one layer earlier:
three pairwise vertex-conflicting triangles force a `C_5`, and HYP-+2880
verified the cycle-space identity `C_5 = triangle xor triangle xor triangle`.

On the LRC side, the totient reframe is the same move.  The copy rule
`sum_{d|n} c(d)=n` gives `c=phi`, and `phi(d)` is exactly the count of
exact-denominator packets on a `q`-grid.  HYP-2630/HYP-2632 show why that
primitive layer must be retained: Euler-copy capacity is uniform on `F_7^*`;
the actual split is the later `chi_7` phase plus the affine zero lane
`a+b=2 mod 7`.

So the abstract bridge is not "`H` is like `phi`" as a number-theoretic slogan.
The bridge is operational:

```text
odd cycles / strong atoms / exact-period denominators
  -> incidence or coimage quotient
  -> signed character or forced-hole boundary
  -> scalar count.
```

This also explains why several older routes stalled.  Raw divisor counts,
squarefree masks, raw rooted cut weight, and raw tournament classes are all
valid addresses, but they are not proof carriers.  They delete the primitive
incidence before the obstruction has been forced.

The next concrete experiment should build an LRC packet-conflict graph from
HYP-2632: vertices are exact-period phi packets or additive frequency shells;
edges mark incompatible signed covering packets.  Then test whether the
affine zero lane and Legendre selector appear as odd holes or odd clique
shadows, paralleling `K_3 -> C_5` on the tournament side.
