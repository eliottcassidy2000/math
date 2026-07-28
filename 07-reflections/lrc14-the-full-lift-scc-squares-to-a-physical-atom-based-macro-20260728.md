# The full lift SCC squares to a physical atom-based macro

**Status: FINITE-EXACT + VERIFIED in normal and optimized modes.  This is a
whole-open-cylinder support statement, not a target action, relation-index
identification, semantic endpoint current, row exclusion, or LRC(14)
conclusion.**

## Inheritance pass

The closest proved mechanism is THM-2707: `3346` packet addresses share one
open cylinder `I` and form the complete directed eleven-partite physical-lift
graph.  Its canonical hostile is the frozen THM-2680 `following` atom, absent
at the old terminal address and at ten vertices of the displayed eleven-word.
The corrected near miss is THM-2706: ordinary `D`-arrow factorizations of the
transverse `C_221` macro require two nonphysical ghosts.  The least-used
sidecar is the exact residue part of the THM-2707 lift graph, rather than its
SCC or total edge count.

That sidecar changes the picture.  The old following atom is not sporadically
present at one lucky vertex.  It is exactly one complete residue part of the
physical packet graph.

## 1. Exact atom bank

Retain THM-2707's open cylinder

```text
I=(960117507257/1930018885886,
   324519717452867/652346383429468)
```

and its packet addresses `n in Z/(13^6)`.  The frozen THM-2680 following atom
has metadata

```text
rail=(source,owner,deep)=(1,0,12),
future=1, j=2, h=6, epsilon=1, kappa=1.                  (1)
```

An independent exact scan gives

```text
packet bank P:                    3346 addresses,
following-atom bank A:            304 addresses,
transit bank P\A:                 3042 addresses,
A={n in P:n=0 mod13}.                                    (2)
```

The three descriptions of `A` agree exactly:

```text
following atom at the midpoint
 = following atom on every point of the open I
 = residue-zero part of the packet graph.                (3)
```

Equation `(3)` includes the delayed prefix.  It is not inferred from one
sample: `13^6(z+7n/13^6)` differs from `13^6z` by the integer `7n`, so the
delayed coordinate is constant in `n`, and exact wall containment verifies
the whole open cylinder for every one of the `304` endpoints.  Equality at an
open wall is allowed only at the excluded endpoints of `I`.

## 2. Squaring the lift relation restores the missing part

Let `R_lift` be the intrinsic binary relation on packet addresses

```text
n R_lift n'  iff  n mod13 != n' mod13.                   (4)
```

It is a bidirected complete multipartite relation, not a tournament: equal
residues are genuine ties.  THM-2707 shows that `(4)` is exactly the physical
nonzero-lift condition.  Since `A` is one entire part,

```text
R_lift intersection (A x A)=empty,
(R_lift o R_lift) intersection (A x A)=A x A.            (5)
```

Every ordered pair of atom endpoints has all `3042` transit packets as
physical midpoints.  Numerically,

```text
A -> P\A edges:                  924768,
P\A -> A edges:                  924768,
two-step A-to-A macro paths:  281129472.                 (6)
```

The smallest control is the exact loop

```text
0 -> 1 -> 0,
signed lift numerators (7,-7),
positive representatives (7,13^6-7).                    (7)
```

Both steps are nonzero modulo `13`, their signed sum is zero in `Z`, and the
pure translation phase therefore has trivial holonomy at every integral
frequency.  The whole interval `I` carries the fixed packet at the midpoint
and the full following atom at both endpoints.

This is a useful relational version of transitivity.  A one-step observable
cannot move inside the atom part, but its square is the universal relation on
that part.  The irregularity is concentrated in the first quotient digit:
one residue class carries the atom and the other ten active classes supply
all physical toothpicks through which the macro closes.

## 3. Why this does not contradict the ghost theorem

THM-2706 and `(5)` live in different categories.

```text
THM-2706: ordinary chronological D-arrows between transverse C_221 phases;
here:     physical slope-lift edges inside one fixed packet cylinder.
```

The former has no allowed midpoint and needs two formal ghosts.  The latter
has `3042` physical midpoints for every ordered atom-endpoint pair.  Thus
"midpoints are missing" was never an object-free statement; it depends on
which operation supplies the arrows.  The new macro does not advance the
delayed base and hence cannot replace a chronological `D` edge.

## 4. The remaining type obstruction

The source label `1` in `(1)` is a retained rail label, not a proved exclusive
THM-2305 source.  Likewise the lift residue in `(4)` is a physical deck
coordinate, not automatically either Fourier multiindex in a THM-2334
endpoint pair.  Positive support and zero phase holonomy therefore provide a
based carrier but not its complex amplitude or target character.

The next object should be a fibre product, not an identification:

```text
(physical two-step atom macro)
          x_(fixed endpoint expansion)
(paired blocker--graft Fourier current).                 (8)
```

It must retain the left and right Fourier multiindices `(u,v)`, test the true
difference character `eta.(u-v)`, and sum over the full endpoint constraint
`X`.  Treating a packet address or lift residue as `u` would discard exactly
the coordinate that THM-2563 and the THM-2701 paired-cylinder addendum say is
missing.

## 5. Exact reproduction

Run

```bash
python3 04-computation/lrc14_following_atom_based_macro_probe_20260728.py
python3 -O 04-computation/lrc14_following_atom_based_macro_probe_20260728.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_following_atom_based_macro_probe_20260728.out
```

with LF-normalized SHA-256 values

```text
script  db2aef5fa48611c61939494630e22298e57fe20c7eeb81f08a47ab27869d61dd
output  de0988d683aa1b237a56e5ef73dbef27c136ebf169ddc934f5d040725ca62db6
```

The universe is all `13^6` lift addresses with the full inherited packet
filters.  Positive controls are the complete THM-2707 census and `(7)`;
hostile controls are the `3042` nonzero-residue packets, none of which even
has the following atom at its midpoint.  No scalar obligation is removed.
