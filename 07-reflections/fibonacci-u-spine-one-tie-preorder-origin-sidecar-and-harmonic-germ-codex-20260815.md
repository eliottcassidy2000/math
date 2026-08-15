# Fibonacci-selected U-spine: one-tie preorder, origin sidecar, and harmonic germ

**Research synthesis, 2026-08-15.**  Truth-bearing claims are in
[THM-3454](../01-canon/theorems/THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order.md)
and its exact companion.  This reflection records the concept board, quotient
anatomy, and next probes.  It does not add an LRC or Jacobian consequence.

## Outcome

The hidden sequence

```text
q_t=4t(t+1)+3=(2t+1)^2+2=Q_(t-1)
```

labels the already-proved parabolic Berggren `U`-spine.  Fibonacci does not
create a second tree branch here.  It selects sparse spine indices `t=F_n`,
at actual rooted depths `F_n-1`, on that branch.

At every pair of selected spine indices, four observables coincide:

```text
leaf-chain distance
= rooted U-tree distance
= Farey cross-determinant magnitude
= sqrt(Lorentz pairing/2).
```

For the four spine indices (whose rooted depths are one lower)

```text
(F_(k-1),F_k,F_(k+1),F_(k+2)),       k>=4,
```

the six distances have levels

```text
01 < 12 < {02,23} < 13 < 03.
```

The correct size-six object is the strict transitive orientation of `K6`
with one missing edge, or the corresponding total preorder.  Calling it a
`T6` would erase the theorem's central equality.

## Inheritance pass

| role | inherited object |
|---|---|
| closest mechanism | THM-3333's Lorentz determinant-square identity and THM-3334's fixed-cusp `U`-fan |
| Fibonacci selector | THM-3339's exact Cassini/Pell descent |
| canonical hostile | THM-3345's source-dependent off-ray Berggren path costs |
| corrected near miss | MISTAKE-394's tied root admitted to a product-`T6` support sentence |
| least-used sidecar | one absolute spine index (equivalently rooted depth plus one), plus adjacent/opposite edge incidence |

The key distinction is between two graphs on the same six labels:

- THM-3339 uses edge **products**, has no tie for `k>=3`, and isolates the
  opposite pair `{03,12}`;
- THM-3454 uses edge **distances**, has one forced tie for `k>=4`, and isolates
  the adjacent pair `{02,23}`.

`S4` preserves incidence, so no relabeling can identify those pairs.

## Concept board

1. **Fixed-cusp leaf chain.**  The fractions `t/(t+1)` are all adjacent to
   `1/1`, but among themselves form an infinite path indexed by `t`.
2. **Parabolic U-spine.**  The same index gives
   `P_t=U^(t-1)(3,4,5)` and exact rooted depth `t-1`.
3. **Fibonacci sparse selector.**  The golden recurrence acts on the selected
   index set, not on the parabolic branch transition.
4. **Six-distance preorder.**  A recurrence equality is visible as one
   missing comparison between adjacent `K4` edges.
5. **Origin kernel.**  Pairwise distances determine spine indices, and hence
   rooted depths, only modulo common translation; the recurrence is not
   invariant under that translation.
6. **Full-germ repair.**  A harmonic total loses a subset, while its full
   Dirichlet germ recovers the coefficients.

Every new object changed the others as follows.  The fan calculation typed
the tree transplant; the tree depth exposed the translation kernel; the
kernel explained why order and Pell data still admit false windows; the gap
view turned the tie into a recurrence equation; and the harmonic view supplied
an independent scalar-versus-germ instance of the same loss.

## The gap transplant is the mechanism

Write an increasing spine-index window in gap coordinates:

```text
x_0<x_1<...<x_m,             g_i=x_(i+1)-x_i.
```

Then the vertex sequence is Fibonacci-recurrent iff

```text
g_1=x_0
and
g_i=g_(i-1)+g_(i-2) for i>=2.
```

The four-point conditions are exactly

```text
g_1=x_0                 <=> d_12=x_0,
g_2=g_0+g_1             <=> d_23=d_02.
```

So the adjacent-edge tie is not numerology.  It is the first recurrence law
of the gap sequence.  The other condition aligns the gap sequence with the
vertex sequence at one marked origin.

If `h_i=x_i-1` denotes actual rooted depth, the same index recurrence becomes
`h_(i+1)=h_i+h_(i-1)+1`, and the marked seam is `d_12=h_0+1`.  Metric
differences do not see this affine shift; absolute recurrence equations do.

This also explains the sharp hostile:

```text
(2,3,5,8)     and     (1,2,4,7)
```

have exactly the same six labelled distances and both initial pairs have
absolute Cassini norm one.  Translation preserved all pairwise distances and
accidentally preserved the absolute Pell gate, but destroyed the marked seam
`d_12=x_0`.

## Typed map ledger

| connection | source | target | preserved | destroyed | sidecar | cheapest test |
|---|---|---|---|---|---|---|
| fixed-cusp transplant | index `t` | fraction leaf and `U`-node | index difference | full Farey graph distance | fixed cusp and branch label | `t=1,4` |
| Lorentz lift | ordered spinor `(t+1,t)` | null triple `P_t` | determinant square | determinant sign | chart order | swap hostile `(t,t+1)` |
| four-window quotient | absolute spine indices (rooted depths plus one) | six distances | all differences | common translation | one index/depth coordinate | true/false pair above |
| weak order | six distances | five-level preorder | comparisons and equality | magnitudes | exact weights for reconstruction | `(2,4,7,12)` |
| harmonic valuation | subset of `N` | one scalar at `s=1` | additive mass | term labels | full Dirichlet germ | `{2}` versus `{3,6}` |

The third row is another exact instance of the controlled-forgetting card in
`META-PATTERNS.md`.  For

```text
q: Z^4 -> Z^6,          q(x)=(x_j-x_i)_(i<j),
```

the kernel is the diagonal translation copy of `Z`.  The recurrence locus is
not saturated under that kernel.  One absolute spine index, equivalently one
rooted depth with the chart shift remembered, separates every seam.
This is structurally parallel, but not mathematically identified, with:

- the divisor seam in the exact-six `k=2` body quotient;
- the missing margin modules in THM-3450's centred D5 carrier; and
- the delayed central digit in THM-3449/3452's Hensel Heisenberg action.

## Branch labels and the order-six recurrence

The full branch labels have second difference eight.  Fibonacci sampling
instead gives

```text
R_n=q_(F_n)=Q_(F_n-1)=3,11,11,27,51,123,291,731,1851,...
```

and the minimal characteristic polynomial

```text
(x^3-2x^2-2x+1)(x^2-x-1)(x-1).
```

The six modes come from the square-Fibonacci, Fibonacci, and constant pieces
of `4F_n^2+4F_n+3`.  Their count happens to be six, but no canonical map to
the six `K4` edge labels has been proved.  Forcing one would be exactly the
kind of cardinality analogy the tournament validity gate forbids.

The leading `3` is only the polynomial extension at the degenerate depth
`-1`; the rooted geometric sample starts at `n=1`.  Independently, the first
`6 x 6` Hankel determinant is `393216`, so no recurrence of smaller order can
fit the sequence.

## Subsets of the harmonic series

The set-valued maps

```text
A |-> {t/(t+1):t in A},
A |-> {P_t:t in A},
A |-> {1/t:t in A}
```

are faithful Boolean embeddings.  Summing the last set is not faithful:
`1/2=1/3+1/6`.  The full function

```text
D_A(s)=sum_(t in A)t^(-s),      s>1,
```

is faithful by the least-differing-coefficient limit as `s->infinity`.
This is a clean scalar/full-germ boundary.  On Fibonacci-indexed subsets even
the convergence bit disappears, because every reciprocal-Fibonacci subseries
converges.

## New probes

### 1. Four-point turnpike preorder atlas

For arbitrary positive gaps `(a,b,c)`, the six line distances are

```text
a,b,c,a+b,b+c,a+b+c.
```

Classify their realizable weak orders by the adjacent/opposite incidence orbit
of every tie.  The Fibonacci chamber is the wall `c=a+b` with `a<b`.  This
could provide a lawful tournament-with-missing-edges grammar for LRC speed
difference data without pretending every chamber is a tournament.

### 2. Origin-restored exact-six transport

The exact-six body graph needed divisor data to stop false composition.
The present metric quotient needs one root coordinate.  Test whether the
smallest exact-six sidecar is likewise a one-coordinate origin in each
divisor chart, rather than the full divisor label.  The decisive object is the
kernel-pair conflict graph, not another SCC census.

### 3. Dirichlet-germ current rather than harmonic mass

Where the repo currently retains only one harmonic coefficient, ask whether a
finite jet or full Dirichlet/Fourier germ separates the live quotient seams.
THM-3450 warns that centred interaction alone misses margins; the first probe
should explicitly include those margin modes and an adversarial equal-value
pair.

### 4. Do not overread the order-six coincidence

Before comparing the six recurrence modes with six tournament edges, require
an explicit equivariant map, preserved predicate, and hostile.  The first
hostile is edge incidence: the mode factorization has no demonstrated
adjacent/opposite `K4` partition.

## Status boundary

THM-3454 is an elementary but exact bridge among already-proved carriers.  It
does not move the LRC(14) live ledger, does not supply the `7 x 13` spectral
closure, and does not classify Keller counterexamples.  Its useful export is
the mechanism: a recurrence may live as a tie in a distance quotient, but one
marked origin is required to lift it back to the absolute sequence.
