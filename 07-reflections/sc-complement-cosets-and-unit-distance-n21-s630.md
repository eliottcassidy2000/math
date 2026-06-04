# SC Complement Cosets And Unit-Distance n=21 S630

The incoming S629/THM-409 packet sharpened the object.  A self-complementary
tournament does not usually hand us a single canonical "flip every edge and
swap these vertices" map.  It hands us a coset:

```text
Anti(T) = Iso(T, T^op).
```

Once one anti-automorphism `tau` exists, every other one is `aut * tau`.
That is the mathematical object describing the back-and-forth swap.  It is a
torsor for `Aut(T)`, or equivalently the odd half of the signed automorphism
group `Aut^+/- (T)` with a sign map to `C2`.

The node-perspective rule is very clean.  Under any anti-map, a vertex of score
`d` is sent to a vertex of score `n-1-d`.  Under all anti-maps at once, a vertex
is not sent to one distinguished vertex; it is sent to an automorphism orbit.
In this S630 supplement it is exact for every SC class through `n=7`:
`|Anti|=|Aut|`, score
reversal never failed, and the anti-image set of every vertex was an
automorphism orbit.

This makes SC tournaments less mystical and more usable.  The perspective of a
node changes by dualizing its out-neighborhood/in-neighborhood and then
forgetting coordinates up to automorphism.  Rigid SC classes really do swap
vertices one-to-one.  Symmetric SC classes swap orbit-to-orbit.  The regular
cyclic `R_7` is the extreme: every vertex has the same perspective, so the
anti-coset sends every vertex across the whole single automorphism orbit.

That also helps with the H=7/H=21 relation.  Self-complementarity is a duality
mechanism, not a repair mechanism for forbidden Hamiltonian-path counts.  S630
finds `H=7` and `H=21` absent inside exact SC classes through `n=7`, which is
only finite evidence but fits the stronger HYP-2200/THM-115 guardrail.  The
anti-coset explains how a class mirrors itself; it does not make the forbidden
OCF scalar legal.

On the unit-distance side, the useful discovery is the `n=21` split:

```text
57 = 20 + 37.
```

The `20` is the unit Hamiltonian spine.  The `37` is the tile/bulk surplus, and
it is the centered-hex number `3*3*4+1`.  This is much better than treating
`21` as a raw echo of forbidden `H=21`.  The exact Moser row has 21 vertices,
but its edge carrier is a spine plus a centered-hex bulk.  It is not trying to
be a tournament with 21 Hamiltonian paths.

So the current picture is:

- `7 = Phi3(2)` and `21 = 3*Phi3(2)` remain meaningful as Eisenstein/cyclotomic
  scalars;
- tournament `H` forbids those scalars as Hamiltonian-path evaluations;
- unit-distance constructions may show those scalars as decomposed edge
  carriers;
- SC tournaments explain the complement/conjugation side through anti-cosets,
  not through raw H-value matching.

The next sharp computation is not another numerology table.  It is a carrier
audit of the exact `n=21` graph6 cores: do all of them preserve the `20+37`
spine/bulk story, and does the `37` bulk actually behave like a centered-hex
shell under deletion and frontier extension?

The roots-of-unity piece should stay alive but disciplined inside HYP-2206.
`Phi3(2)=7` and `3*Phi3(2)=21` are not arbitrary, but they are also not
licenses to identify carriers.  The better invariant is: when a
primitive-cube-root scalar appears, what side channels must be retained for the
statement to remain true?
