# Sequence Shadow Recursions S633

The prompt was a good antidote to brute-force tunnel vision.  The tempting
thing is always to ask for the next value of the hard sequence.  But the repo
has already learned, painfully, that raw next terms can be a mirage: A000568,
self-converse tournaments, LRC source menus, and unit-distance carriers all
become useful only after the side channel is named.

So I treated the self-converse sequence

```text
1, 1, 2, 2, 8, 12, 88, ...
```

as one face of a little family:

```text
T(n)  = all tournament classes
SC(n) = fixed classes under converse
M(n)  = merged complement/converse classes
N(n)  = nonfixed converse pairs
```

This makes the first lesson almost embarrassingly simple.  A hard fixed count
should travel with its merged and nonfixed companions.  The fixed layer is the
special locus; the nonfixed layer is the generic sea; the merged layer is the
quotient the metagraph actually sees.

The second lesson is the bisection.  The even half of SC is exactly `A(m,4)`.
The odd half is not another monster; it is the same base-4 Burnside skeleton
with one fixed vertex charging `2^(#parts)`.  That is the recursive feel I
wanted from this session: the jump from even to odd is not random, it is a
local tax paid by every odd partition part.

The third lesson is that the useful shadows get thinner as they get more
proof-relevant.  A000568 is huge.  SC is the fixed layer.  Round/LRC
self-converse counts are only `2^floor((m-1)/2)`.  Shell transporter orbits are
smaller still, but they remember arithmetic channels like the `C=27` strata
`1,3,9`.  That is exactly the kind of companion sequence we want near an LRC or
unit-distance obstruction: smaller, stranger, and still carrying the hard
predicate.

I do not think the right slogan is "find the recurrence."  It is more like:
for every hard sequence `A`, build its local weather report.

```text
fixed(A), merged(A), nonfixed(A), q-shadow(A),
bisection(A), skinny(A), transporter(A).
```

Then ask which shadow changes when the proof problem changes.  That is a way to
make little bits of progress on famous structures without pretending the next
giant term is the only valid currency.
