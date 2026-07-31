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

### 2.1 The congruence lock descends one scale

Write every atom address as `n=13m`.  The `304` reduced addresses in
`Z/(13^5)` are not an unstructured residue cloud.  They form exactly `35`
consecutive toothpicks: five near each of the seven clock centres

```text
c_j=floor(j*13^5/7),              j=0,...,6.              (8)
```

Here is the complete compressed atlas; each interval is an offset from
`c_j`.

| `j` | `c_j` | five-digit base-13 word | tooth offsets | size |
|---:|---:|:---:|:---|---:|
| 0 | 0 | `00000` | `0..8, 13..32, 37..43, 319..322, 326` | 41 |
| 1 | 53041 | `1B1B1` | `0..9, 13..33, 38..44, 319..323, 327` | 44 |
| 2 | 106083 | `39393` | `0..9, 13..33, 37..44, 319..323, 327` | 45 |
| 3 | 159125 | `57575` | `0..9, 13..33, 37..44, 319..323, 327` | 45 |
| 4 | 212167 | `75757` | `0..9, 13..33, 37..44, 319..322, 327` | 44 |
| 5 | 265209 | `93939` | `0..8, 13..33, 37..44, 319..322, 327` | 43 |
| 6 | 318251 | `B1B1B` | `0..8, 13..32, 37..44, 319..322, 327` | 42 |

The alternating palindromes are forced, not numerology.  Since
`13=-1 mod7`, multiplication by `13` reflects the clock `j/7` to
`(7-j)/7`, while multiplication by `13^2` fixes it.  More explicitly, for
`j=1,...,6`,

```text
j/7=0.overline((2j-1),(13-2j)) in base 13,               (9)
```

and `c_j` is its five-digit truncation.  The odd truncation is therefore a
palindrome.  This is the exact toothpick self-similarity behind the earlier
`35`-component census: seven reflected clocks, each resolved into five fine
teeth.  The small size bias `41,44,45,45,44,43,42` records the oriented
half-open boundary and should not be symmetrized away.

There is also an operation-level descent.  For atom endpoints
`a=13m`, `c=13m'` and any transit midpoint `b`, the signed physical lift
numerators telescope:

```text
7(b-a)+7(c-b)=7(c-a)=13*7(m'-m).                         (10)
```

Thus the two-step macro has a canonical depth-five numerator
`7(m'-m)`, independent of which of the `3042` midpoints realizes it.  The
smallest nonconstant control is

```text
0 -> 1 -> 13:       (7,84) -> total 91 -> descended 7.   (11)
```

One-step physical lifts cannot stay in the atom congruence class; their
two-step composites do, and uniformly admit division by thirteen.  The
valuation is exactly one only when `m'-m` is nonzero modulo `13`, as in
(11); deeper endpoint coincidences can contribute further factors, and a
loop contributes zero.  This is an address coboundary/renormalization law,
not a delayed chronological map.

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
(paired blocker--graft Fourier current).                 (12)
```

It must retain the left and right Fourier multiindices `(u,v)`, test the true
difference character `eta.(u-v)`, and fix one marked frequency `X` while
summing its decomposition fibre.  Summing over all `X` is the recombination
that THM-2568 kills.  Treating a packet address or lift residue as `u` would
discard exactly
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
script  bc81895a029d070d3856a9a9cfeb622a39878947d98aad0768b5c405fdbf6818
output  a636763b55b5471a9a4aa1d6512f3a15d693d6daa1918851d1a62c22dfbcc8f6
```

The universe is all `13^6` lift addresses with the full inherited packet
filters.  Positive controls are the complete THM-2707 census and `(7)`;
hostile controls are the `3042` nonzero-residue packets, none of which even
has the following atom at its midpoint.  No scalar obligation is removed.
