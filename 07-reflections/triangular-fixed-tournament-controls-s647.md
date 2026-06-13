# Triangular Fixed Tournament Controls S647

The useful translation is simple enough to feel almost too obvious:

```text
tournament arcs = C(m,2)
fixed Hamiltonian path = m-1 forced arcs
off-path deformation fiber = C(m-1,2) arcs
```

So a Hamiltonian path is not merely an object counted by `H(T)`.  It is also a
section of the triangular tournament carrier.  Once the section is chosen, the
unresolved orientation choices form one smaller triangular grid.

That makes the even perfect controls sharper than they looked in S644.
Euclid-Euler gives:

```text
6    = C(4,2)
28   = C(8,2)
496  = C(32,2)
8128 = C(128,2)
```

but tournaments add the section/fiber decomposition:

```text
C(m,2) = (m-1) + C(m-1,2).
```

The row `m=8` is the loud one:

```text
28 = 7 + 21.
```

This is exactly the permanent H-gap pair, but in a different role.  `7` is the
spine length of a fixed Hamiltonian path.  `21` is the number of off-path cells
left after the spine is chosen.  They are not path counts there; they are
dimensions of the perfect arc-control.

The finite computation makes the guardrail real.  The fixed-base-path fiber at
`m=8` has size:

```text
2^C(7,2) = 2^21 = 2097152.
```

Using the full labelled `n=8` H-spectrum, the exact fixed-fiber count is:

```text
fixed_count_h = h * labelled_count_h / 8!.
```

This reconstructs the whole fixed-path fiber and finds:

```text
H=7  count=0
H=21 count=0
```

So the perfect control sees `7` and `21` exactly as section/fiber dimensions
while refusing them as `H` values.  That is not a coincidence; it is a role
mismatch.

This is probably the point to carry back into LRC.  We keep seeing integers
that are real but typed differently:

```text
H=21             forbidden Hamiltonian-path count
UD n=21          vertex count with a unit-spine/bulk edge carrier
LRC n=14         row whose pair-count aliquot shadow is 21
C=27             LRC shell clock
28=7+21          perfect tournament arc control
```

The proof mistake would be to demand that these are the same scalar.  The proof
move is to ask what coordinate each number occupies: count, section, fiber,
modulus, bulk, aliquot shadow, or carry mass.

There is also an algorithmic lesson.  Fixing a Hamiltonian path turns a
tournament into a fiber of `2^C(m-1,2)` off-path choices.  That is the same
tiling picture the repo has used forever, but now the perfect controls say
which rows are clean calibrators.  At `m=8`, the off-path budget is `21`; if a
method cannot distinguish "21 choices of deformation" from "21 Hamiltonian
paths," it is quotienting away the wrong side channel.

For LRC `n=14`, the analogy I would keep is cautious:

```text
choose a section first,
then measure the deformation fiber,
then only at the end ask the loneliness predicate.
```

That is the same shape as the current `27` shell/owner/carry program.  The
triangular fixed controls do not prove LRC, but they give a clean miniature
where a tempting scalar is correct in one coordinate and impossible in another.

The next useful tool is a typed scalar ledger.  Every occurrence of `7`, `21`,
`27`, and `28` in these cross-problem packets should be labelled as one of:

```text
count
section length
fiber dimension
shell modulus
bulk mass
aliquot shadow
owner/carry residue
```

That small discipline may prevent a lot of false bridges and expose the real
ones faster.
