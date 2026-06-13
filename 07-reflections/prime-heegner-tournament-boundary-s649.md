# Prime-Heegner Tournament Boundary S649

The useful indexing trick is to stop treating the Euler prime streak as a
flat list.  It is a three-part carrier:

```text
x=0       source/base prime p
x=1..p-2 interior prime window
x=p-1    forced square sink p^2
```

With that decomposition, the user's `n-2` hunch becomes precise.  Standard
zero-based indexing says the lucky polynomial has `p-1` prime values,
`x=0..p-2`.  That is exactly the length of a Hamiltonian spine on `p`
vertices.  If the base term `x=0` is set aside, the remaining interior prime
inputs are `x=1..p-2`, count `p-2`, which matches Moon's strong-tournament
triangle floor.

So there are two adjacent tournament shadows:

```text
p-1 -> fixed Hamiltonian spine length
p-2 -> interior strong floor
```

The endpoint is what keeps this from being numerology.  The failure is not an
accidental composite:

```text
f_p(p-1) = p^2.
```

That is a forced square sink.  The Heegner side channel explains why no earlier
sink appears for the six lucky primes:

```text
d = 4p - 1
f_p(x) = ((2x+1)^2 + d)/4
```

and the class-number-one list projects as

```text
{7,11,19,43,67,163} -> {2,3,5,11,17,41}.
```

The small Heegner rows `{1,2,3}` are not noise.  They are boundary exceptions:
`d=1,2` are outside this quadratic shape, and `d=3` gives `p=1`, the degenerate
edge before primality starts.  This mirrors the tournament habit of keeping
the source and sink instead of erasing them from the count.

The new technique I want to keep is the endpoint-square carrier.  Whenever a
sequence has a long clean prefix and a forced failure at the first endpoint,
try to encode it as:

```text
source/base term
interior proof obligations
endpoint failure
spine length
off-spine deformation fiber
```

For this example, the off-spine fiber is `C(p-1,2)`: after fixing a Hamiltonian
path on `p` vertices, those are the remaining tournament arcs.  This does not
prove Heegner uniqueness, but it gives the right shape of transfer.  The
arithmetic proof lives in the norm/class-number side channel; the tournament
proof language lives in the spine/interior/fiber side channel.  The bridge is
meaningful only if both sides stay attached.

This is also a useful guardrail for LRC and unit distance.  In LRC, `n-2`
should often be read as an interior or residual obligation, while `2n-1`
remains the discriminant/root clock.  In unit distance, the same distinction
appears as a unit spine plus hidden bulk and endpoint repairs.  The raw scalar
is tempting, but the proof is almost always hiding in the side channel.
