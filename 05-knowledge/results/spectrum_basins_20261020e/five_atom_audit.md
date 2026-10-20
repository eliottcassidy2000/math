# Independent audit: a regular five-vertex atom with any source context

**PROVED + independently replayed:** for every nonempty tournament A with
a source, the actual black path-reversal graph contains a four-step path
from `A -> C5` to `A -> 1 -> 1 -> C3`. Here C5 is the regular five-vertex
tournament. No arbitrary interior-context theorem is claimed.

This audits atom5 in the parent's `five_atom_probe.json`, with an independent
literal implementation and a complete larger small universe. The proof
requires the source hypothesis at the two middle moves.

Let `a=(a0,...,a_(m-1))` be any Hamiltonian path of A, with source u=a0.
Choose C5 labels b0,...,b4 with outgoing masks `(6,12,24,17,3)`; thus bi beats
its next two cyclic vertices. Write `a^r` for the reversed list. Use

```text
P1=(a,b0,b1,b2,b3,b4),
P2=(b0,a^r,b4,b3,b2,b1),
P3=(b3,b0,b2,b4,a,b1),
P4=(b2,b0,b3,b4,b1,a^r).                               (1)
```

Every path is directed in its current state. This follows by reading the
five-vertex arcs and the alternating direction of the A path. For m>=2,
only three kinds of cross pairs are ever toggled:

```text
{a_last,b0} in P1/P2,
{u,b4}     in P2/P3,
{a_last,b1} in P3/P4.
```

They cancel pairwise. Every internal A path edge is toggled four times;
all other internal A arcs remain untouched. The same cancellation holds
at m1 when first and last coincide, since the three b labels differ.
The final five-vertex outgoing masks are `(4,29,8,1,13)`, literally

```text
b1 -> b4 -> C3(b0,b2,b3).
```

This proves the target identity uniformly in A and its path.

The first move is black because `A -> C5` is nonselfconverse: its first
SCC is a singleton and its last SCC has five vertices. The last move is
black because the target is nonselfconverse: its first SCC is singleton
and its last SCC is C3. A blue path's positional reversal is an
anti-automorphism, and remains one after R, so nonselfconverseness at
either endpoint rules out blue.

For both P2 and P3, positional reversal maps the ordered pair

```text
(u,a_(m-3))  to the reversed pair  (b4,a_(m-4)).          (2)
```

Here the right side means `(r(a_(m-3)),r(u))`, the pair to which an
anti-automorphism must send its orientation. At P2, u is in position m
and a_(m-3) in position3; at P3 they are in positions4 and m+1. Subtracting
positions from m+4 verifies (2) in both cases.

For every m>=5, the first arc in (2) is true because u is the source and
a_(m-3) is not its path neighbor. It has not been toggled. The second arc
is false because A initially dominates C5; a_(m-4) is not u or a_last,
so no prior cross-boundary toggle reaches it. These contradictory known
orientations exclude anti-reversal and make both middle moves black.

The cutoff m>=5 matters: at m4 the first pair is a path edge and the
second includes u. The remaining m1..4 source/path presentations, modulo
renaming along the path, number1,1,1,2. Direct inspection of (1) makes all
four moves black in all five cases. The independent checker includes
all their labelled versions and every larger source context through m6.

The complete independent finite universe is:

| m | Labelled source contexts | All their Hamiltonian path presentations |
|---:|---:|---:|
|1|1|1|
|2|1|1|
|3|2|2|
|4|8|12|
|5|64|192|
|6|1024|7680|

All1100 contexts and7888 presentations give four directed black steps and
the literal restored target. No inherited canonicalizer, class filter,
blue predicate, move function, or certificate compiler is imported.
Known-arc index and toggle checks for m5..100 supplement the uniform proof
above; they are not being used as an all-size census.

```text
python3 05-knowledge/results/spectrum_basins_20261020e/five_atom_audit.py
python3 -O 05-knowledge/results/spectrum_basins_20261020e/five_atom_audit.py
```

Checks use explicit exceptions under both Python modes. The JSON records
the source hash, exact universe, total presentations, and literal target.
