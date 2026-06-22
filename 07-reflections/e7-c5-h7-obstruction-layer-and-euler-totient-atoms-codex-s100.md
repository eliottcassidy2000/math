# E7 C5, H7, and Multiplicative Atoms

The pentagon identification has two levels.

S37 verifies the support-level statement: a directed `C_5` is the XOR of three
pairwise vertex-conflicting triangles, hence the H=7/K3 cycle-space support.
S100 separates that support from the E7 metagraph-hole object.

Under the fixed-path cycle-space map, a directed 5-cycle support is one even
graph class in `E_7`.  An `E_7` chordless `C_5` hole is a five-class cycle in
the quotient metagraph.  The S100 audit makes that distinction exact: the
directed pentagon support is class `3`, while E7 has `1496` metagraph C5 holes
touching `48/54` classes.  The K3-forces-pentagon classes hit `835` of those
holes, but no hole is fully made of them.

So the metagraph-hole match is not a single object.  It is support equality
plus a quotient-level obstruction layer:

```text
scalar forgetting -> cycle-space closure -> first odd-cycle obstruction.
```

For `H=7`, the scalar target wants `alpha=(3,0)`, or `Omega=K_3`.  Incidence
closure in true tournaments blocks that point: three mutually-conflicting
triangles force a directed pentagon and therefore extra odd-cycle mass.  For
`E_7`, quotienting the cycle space first creates metagraph pentagons at the
same apex prime.  They are different graphs, but they mark the same failure of
a scalar/quotient shadow to retain incidence.

The Euler-totient thread is the right arithmetic analogy.  In the LRC residue
work, `phi(D)` counts exact-denominator packets before the scalar witness count
is evaluated.  In tournaments, strong components supply the multiplicative
packet law:

```text
H(T)=product H(C_i).
```

This is an Euler product over strong atoms, not ordinary integer factorization.
That distinction matters: `49` and `75` are composite integers but single
irreducible strong atoms, while `21=3*7` is absent because the strong atom `7`
is absent.  The atom-count function for strong `H` values is the tournament
analogue of the exact-period totient packet, and it should be kept before
scalarizing to raw `H`.

The next test should not try to make an E7 metagraph C5 hole equal `Omega=K_3`.
It should compare incidence profiles: which E7 C5 holes pass through the directed
pentagon support class, which pass through the five `k3_forces_pentagon`
classes, and whether the complement-paired C7 heptagon classes have an
analogous strong-atom packet signature.
