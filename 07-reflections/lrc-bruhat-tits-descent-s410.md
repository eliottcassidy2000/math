# LRC Bruhat-Tits Descent

codex-2026-05-31 S410

The existing LRC endpoint program already has a p-adic geometry hiding inside
it.  Every forbidden endpoint is rational:

```text
e = (n*m +/- 1)/(n*v).
```

If we look only at primes dividing `n`, the reduced denominator of `e` gives a
depth vector.  That vector is the shadow of a point on the boundary of a
product of Bruhat-Tits trees.

This makes THM-360 feel less like a one-off divisibility trick.  The unit
endpoints `a/n` are exactly the first visible boundary layer:

```text
d_p(a/n)=v_p(n),  p | n.
```

Only a speed divisible by `n` protects them.  In tree language, only an
`n`-gate translates all unit rays back to the root.

But then the gate creates its own endpoints.  If `v=n*q`, those endpoints are

```text
(n*m +/- 1)/(n^2*q).
```

At every `p | n`, the numerator is a p-adic unit, so their depth is exactly

```text
2*v_p(n)+v_p(q).
```

That is the whole Bruhat-Tits descent slogan:

```text
gate repair is debt export to child vertices.
```

For `n=14`, the initial tight row exposes the six unit endpoints at depth
`{2:1,7:1}`.  A `14`-gate kills that layer, but the S380 14-multiple ladder
still has exposed endpoints at `{2:2,7:2}` and `{2:4,7:2}`.  The anomaly is
real: composite structure changes where the moat lives.  It does not remove
the moat.

For `n=16`, this becomes a clean rank-one model.  The initial row exposes
`{2:4}`.  A `16`-gate exports to `{2:8}`.  THM-367 is then exactly a radial
tree law: owner `u=2^k` has endpoints at depth `4+k`, and a protector
`p=2^j*q` moves upward by `j` edges, with the odd residue deciding branch
selection.

The hard wrinkle is the nine-fan.  From `u=16` onward, lower protectors can
cover all endpoints of a single dyadic owner.  So a proof cannot say "every
branch has a private leaf" locally.  It needs a global product-tree potential:
closing one branch must consume residue budget or primitivity budget needed
to close a descendant branch.

The proof shape I would pursue:

```text
1. Convert every protection arrow into a local tree displacement vector.
2. Add residue branch labels, not just depths.
3. Show every directed labelled endpoint cycle has positive divergence unless
   it uses a super-gate.
4. Show super-gates export to strictly deeper children and cannot occur
   indefinitely with only k speeds and gcd 1.
```

This fits the repo's earlier Bohr-boundary descent, Cayley-Dickson debt, and
dyadic endpoint-flow threads.  Bruhat-Tits is the cleaner language for the
same mechanism: LRC is trying to forbid compact all-protected cores in a
finite shadow of a p-adic building.
