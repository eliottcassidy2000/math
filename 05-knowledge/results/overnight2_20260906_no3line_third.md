# The third factorial moment detects eight-cycles and quadratic four-cycle dependence

**Status: FINITE-EXACT at n=8, with a PROVED contrast reduction, complete
two-path union census, and independent source/profile audit PASS.**

In the saturated no-three-in-line model, independently label the two shores
of a simple bipartite 2-regular skeleton `G` uniformly by `{0,...,n-1}`.
Its edges become grid points. Let `X_3` count unordered nonaxis collinear
triples of those points, and write

```text
M_3(G)=E[(X_3)_3]=E[X_3(X_3-1)(X_3-2)].
```

At `n=8`, the third factorial moment has the form

```text
M_3(G)=A+B*c_4(G)+D*c_6(G)+E*c_8(G)+F*c_4(G)^2,

E=172483/529200 !=0,
F=11881/50400 !=0.                                   (1)
```

Here `c_L` counts connected `L`-cycle components, with `L` the full graph
cycle length. This report computes the complete coefficients `E,F` of
the full third factorial moment. It does not calculate `A,B,D` from the
restricted census, and it makes no assertion that `E,F` have these values
or signs at other grid sizes.

In particular,

```text
M_3(C_8 disjoint-union C_8)-M_3(C_16)
    =172483/264600 >0.                               (2)
```

These two skeletons have the same mean and variance by the proved
short-cycle theorem, but their distributions of `X_3` differ. The new
dependence was allowed by the edge-budget theorem; (1)--(2) show that
the exact geometry-weighted coefficients actually survive cancellation.
No zero-defect probability or asymptotic density conclusion is asserted.

## Inheritance, live board, and exact question

The closest proved mechanism is
[the short-edge homomorphism/copy theorem](overnight_20260906_moments_pairprofile_theorem.md).
The inherited joint-event instrument is
[the full second-moment compiler](../../04-computation/overnight_20260906_no3line_pairprofiles.py).
The hostile is `C_16` versus two `C_8` components: they have identical
`n,c_4,c_6`, so the first two moments cannot distinguish them. The corrected
near miss is replacing a term allowed by a universal polynomial with a
term whose weighted coefficient has been proved nonzero. The sidecar
needed to decide that distinction is the exact multiplicity of geometric
event triples with each complete bipartite union type.

The five-object board is: line triples; non-induced edge copies; short-cycle
weight; signed coefficient contrasts; whole union-incidence multiplicity.
The cheap positive test was an eight-edge cycle copy; the decisive test
must include every other contributing type because cancellation is large.
The census below handles that cancellation in full.

## 1. Why only two coefficients and two edge sizes need computation

An ordered triple of distinct line events uses at most nine grid points.
Its union has a simple bipartite incidence graph `H`. If its degree exceeds
two, it has no copies in a 2-regular skeleton and contributes zero.

The inherited theorem says that the number `C_H(G)` of non-induced edge
copies is a polynomial in cycle counts; every monomial satisfies

```text
sum_L L*exponent(c_L) <= |E(H)|.                      (3)
```

Thus the only cycle monomials permitted for nine edges are
`1,c_4,c_6,c_8,c_4^2`. In particular `c_4*c_6` is excluded because its
weight is ten. The same theorem applies at `n=8` after omitting graph
prototypes with too many vertices in either shore.

For each fixed profile and for its geometry-weighted sum, let

```text
G_0=C_16,
G_8=C_8 disjoint-union C_8,
G_4=C_4 disjoint-union C_12,
G_44=C_4 disjoint-union C_4 disjoint-union C_4 disjoint-union C_4.
```

The coefficient-extraction formulas are exactly

```text
E=(M_3(G_8)-M_3(G_0))/2,
F=(M_3(G_44)-4*M_3(G_4)+3*M_3(G_0))/12.              (4)
```

The denominators are fixed by the cycle counts: `G_8` has two eight-cycles,
while the second contrast cancels its constant and linear four-cycle
terms and leaves `16F-4F=12F`. None of these skeletons has a six-cycle.
Both contrasts annihilate every profile with at most seven edges, by (3).
Therefore enumerating only eight- and nine-edge event unions computes
the full coefficients `E,F`, rather than a truncation of those coefficients.

## 2. Exact incidence and label weights

Let `W_8(H)` count **unordered sets of three distinct** collinear grid
triples whose union is of shore-preserving type `H`. Each corresponds to
exactly `3!=6` ordered event triples in the third factorial moment.

For an `H` with `r` row vertices and `s` column vertices, write

```text
N_8(H)=C_H(K_8,8)=(8)_r*(8)_s/aut_shore(H).           (5)
```

Independent uniform row and column permutations act transitively on the
grid copies of `H`. Hence a fixed grid copy is contained in the random
skeleton with probability `C_H(G)/N_8(H)`, and

```text
M_3(G)=6*sum_H W_8(H)*C_H(G)/N_8(H).                  (6)
```

The contrast is applied to the complete `C_H(G)` counts in (6).
Different event triples can have the same union; all such multiplicities
remain in `W_8(H)`. The copies are non-induced, so additional skeleton
edges between the chosen vertices are allowed. No independence between
the three collinearity events is used.

For degree-at-most-two graphs, components are paths and even cycles.
The sorted list of `(edges,row_vertices,column_vertices,is_cycle)`
therefore retains the whole shore-preserving type. A cycle with `L` edges
has `L` shore-preserving automorphisms; an even-edge path has two; an
odd-edge path has one. Repeated equal components supply their multiplicity
factorials. These are the exact denominators used in (5).

## 3. Complete bounded census and signed cancellation

The independent determinant and slope constructions both recover exactly
648 nonaxis grid triples at `n=8`. The authoritative event-universe size
is computed as the integer

```text
binomial(648,3)=45,139,896.
```

The exhaustive unordered-event census partitions that universe as follows:

| Category | Count |
|---|---:|
| Union has a vertex of degree greater than two | 20,781,112 |
| Degree at most two, but fewer than eight union edges | 2,939,296 |
| Retained eight- or nine-edge union | 21,419,488 |
| Total | 45,139,896 |

There are 150 retained incidence types. Every copy profile on every one
of the seven 2-regular cycle types at `n=8` is counted independently in
Python by enumerating skeleton edge subsets through size nine. The code
checks the full polynomial profile law on all seven types before taking
the geometry-weighted sum. The reduced census is practical: the native
45-million-event pass took about four seconds in the pilot environment.

The exact contributions of the two retained edge sizes are:

| Union edges | Contribution to `E` | Contribution to `F` |
|---|---:|---:|
| 8 | `768463/2116800` | `456371/2116800` |
| 9 | `-26177/705600` | `42631/2116800` |
| Total | `172483/529200` | `11881/50400` |

Summing separately over positive and negative type contributions shows
why an isolated geometric witness was insufficient:

```text
E = 25988161/2116800 - 8432743/705600
  = 172483/529200,

F = 5120467/705600 - 4954133/705600
  = 11881/50400.                                    (7)
```

Besides (2), the resulting raw quadratic contrast is

```text
M_3(G_44)-4*M_3(G_4)+3*M_3(G_0)=11881/4200 >0.       (8)
```

The full certificate manifest retains all 150 union types, their unordered
geometric multiplicities, their grid-copy denominators, their copy counts
in all seven skeletons, their five polynomial profile coefficients, and
both signed ordered-moment contributions.

## 4. Independent and hostile controls

The native program's fast validity test first forms the union of two
event triples. A third event is allowable exactly when each of its new
edges avoids every already-saturated row and column; existing edges are
allowed because the union is a set. This is equivalent to the complete
union having degree at most two, since each event is itself a matching.

A second complete native pass bypasses this shortcut and rebuilds all
literal vertex degrees for **every** one of the 45,139,896 event triples.
It recovers the same partition and all 150 multiplicities. The two native
passes share their component-signature routine; they are a full validity-
filter replay, not a claim of two wholly independent implementations.

At `n=4`, a separate Python set/adjacency census agrees with every native
union multiplicity and every rejected union. A third route enumerates all
`4!^2=576` independent shore labelings for each skeleton and counts the
collinear triples directly on each resulting board. It gives:

| Skeleton | `E[X_3]` | `E[(X_3)_2]` | `E[(X_3)_3]` |
|---|---:|---:|---:|
| `C_8` | `2` | `61/18` | `6` |
| `C_4` disjoint-union `C_4` | `2` | `74/9` | `104/3` |

The third column recovers the inherited variance controls, and the fourth
agrees with formula (6), independently fixing the ordered factor six and
the automorphism/label denominators. Additional exact controls verify that
the two large-skeleton contrasts vanish on every profile with at most
seven edges, that two `C_8` components have exactly two eight-cycle copies,
and that four `C_4` components have six copies of the two-four-cycle graph.

## 5. Reproduction and remaining boundary

```powershell
python 04-computation/overnight2_20260906_no3line_third.py
python -O 04-computation/overnight2_20260906_no3line_third.py
```

The driver requires `g++` on `PATH`, compiles the native enumerator in a
temporary directory, and uses Python integers and `Fraction` for all
counts and weights. The optimized run also compiles with `-DNDEBUG`.
All 357,857 Python gates and all native runtime gates remain active.

- [Native census source](../../04-computation/overnight2_20260906_no3line_third.cpp)
- [Exact weighting and independent controls](../../04-computation/overnight2_20260906_no3line_third.py)
- [Normal transcript](overnight2_20260906_no3line_third.out)
- [Matching optimized transcript](overnight2_20260906_no3line_third_optimized.out)
- [Complete signed certificate manifest](overnight2_20260906_no3line_third_certificates.json)
- [Independent source, copy-profile, and normalization audit](overnight2_20260906_no3line_third_audit.md)

The independent audit constructs skeleton copy profiles by cyclic edge runs
and component convolution, without importing either producer. It recovers
all 1,050 relevant counts in the seven skeletons, both exact coefficient
sums, and the signed cancellations. Its separate size-four label census
also recovers the full histograms and the first three factorial moments.
It reviews the native geometry source and the matching transcripts; its
independent reconstruction takes the frozen size-eight geometric
multiplicities as input rather than claiming a third full geometry census.

Frozen SHA-256 values over LF bytes:

```text
C++ source:   bc6b901088573ab2642cbca8b7cf4856110090869a18a5d6940d2110430ce83e
Python:       4a00e332070e7415cc6fe02999f226a1072dbcbbfbaac95120dceed706eb293e
both outputs: a8c1c4a7ecf28f19699f1b2700b73b72f9ccf0a893b53f0c0e3d1df73f708bc1
certificates: f2c566ac1b2bcedb530af72f3a290479841e7db87844b331558bed68c93ba727
```

The next research question is how the exact geometric coefficient sums
depend on `n`; this result supplies their first decisive nonzero controls.
Their positivity here does not justify a sign conjecture at all sizes,
an independence approximation, or a deduction about `P(X_3=0)`.
