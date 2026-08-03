# Composite CRT-even boundary blindness and its interior collapse

**Status: PROVED + INDEPENDENTLY FINITE-EXACT AUDITED.**  The results below
concern the finite right-boundary and full open-cell micro-staircase models.
They do not realize arbitrary speed tuples, retain cell widths or endpoints,
or prove LRC(14).

## Inheritance pass

- Closest proved boundary mechanism: THM-3316, whose prime-field
  interpolation classifies prime right-boundary blockers.
- Closest all-modulus interior mechanism: THM-3317, whose unique-protector
  cells and exact gcd shift costs expose every perturbation of total cost
  strictly below `n`.
- Canonical hostile: at `n=14`, the normalized one-defect vector
  `(0,1,0,...,0)` blocks all 196 right-boundary tests but misses interior
  cells.
- Least-used sidecar made decisive here: the CRT endpoint-owner hypergraph,
  refined by the parity class of bad shifts for the half-turn residue `7`.

## Theorem 1: a large right-boundary blind flat at every `n=2p`

Let `p` be an odd prime and `n=2p`.  In the normalized right-boundary model,
let `E_p` be the `(p-1)`-coordinate flat

```text
w_i=0 unless i is one of 2,4,...,2p-2.
```

Then `w in E_p` blocks every right-boundary test `(a,s)` if and only if at
least one even coordinate of `w` is zero.  Hence the flat contains exactly

```text
(2p)^(p-1)-(2p-1)^(p-1)
```

right-boundary blockers.

### Proof: the zero-coordinate direction

Assume `w_(2j)=0` for some even index `2j`.  For every boundary parameter
`a mod 2p`, a zero coordinate already has endpoint bin `0` or `-1`, so it
blocks all shifts:

1. If `a=0`, use coordinate `1`.
2. If `gcd(a,2p)=1`, use the odd coordinate `i=-a^(-1) mod 2p`; then
   `ai=-1`.
3. If `gcd(a,2p)=2`, use the odd coordinate `i=p`; then `ai=0`.
4. If `a=p`, use the assumed zero even coordinate; then `a(2j)=0`.

These are all gcd classes because `p` is prime.

### Proof: the all-active direction

Assume all `p-1` even coordinates are nonzero and take `a=p`.  Every fixed
odd zero coordinate has bin `p`, which is interior, while every active even
coordinate has bin `0`.  For a nonzero residue `v mod 2p`, its bad shifts at
this boundary are:

```text
gcd(v,2p)=1:  {0,-v^(-1)}, with the second shift odd;
gcd(v,2p)=2:  {0,p};
v=p:          all even shifts.
```

If no coordinate has residue `p`, shift `2` is safe.  If some coordinate has
residue `p`, it covers only the `p` even shifts; the remaining at most `p-2`
coordinates cover at most `p-2` of the `p` odd shifts.  An odd safe shift
remains.  Thus the vector is not a right-boundary blocker.

At `n=14`, this proves that the six-dimensional even flat has exactly

```text
14^6-13^6 = 2,702,727
```

boundary blockers, of which `2,702,726` are non-scalar false positives.  The
hostile `(w_2=1)` is merely the first point of this large CRT-blind family.

## Theorem 2: interior cells collapse every `2p` blind flat

For every odd prime `p`, the zero vector is the only full open-cell
micro-staircase blocker in `E_p`.  Equivalently, every one of the

```text
(2p)^(p-1)-1
```

nonzero vectors in this flat has an explicit safe open cell and shift.  In
particular, at `n=14` all `14^6-1=7,529,535` nonzero vectors are exposed.

### Ordinary residues leave an odd shift

For each even `j`, use the unique-protector cell `alpha_(2p,j)` from
THM-3317.  Its bin at `j` is `2p-1`, and every other coordinate has an
interior bin.  Therefore, when `j` is supported, every zero coordinate is
harmless on this cell.

For every nonzero residue `v != p`, and for every bin `b`, the two endpoint
congruences

```text
b+s v = 0 or 2p-1 mod 2p
```

forbid exactly one odd shift.  For unit `v`, the two unique solutions have
opposite parity.  For even `v`, the soluble congruence has two solutions
separated by `p`, again of opposite parity.  Thus, if no active coordinate
has residue `p`, the at most `p-1` active coordinates forbid at most `p-1` of
the `p` odd shifts.  One odd shift is safe.

### The half-turn activation graph strictly lowers `2`-adic valuation

For `v=p`, an interior bin forbids all odd shifts precisely when its bin is
`p-1` or `p`; otherwise it forbids no odd shift.  Put a directed edge
`j -> h` between distinct even indices when the bin at `h` on the
`j`-protector cell is `p-1` or `p`.

This graph has the exact edge rule

```text
j -> h  iff  h=j/2 or h=3j/2 and h is an even index in [2,2p-2].   (1)
```

To prove (1), write `j=2r`, `h=2q`, and write the THM-3317 cell as

```text
alpha=(1-t/(2p))/j,
t=1/2+1/(16p^2)  if r<=(p-1)/2,
t=1/(4p)          if r>=(p+1)/2.                         (2)
```

Set `x=q/r=m+u/r`, with `0<=u<r`.  Integer ratios, and ratios that wrap after
the small subtraction in (2), have fractional part above `3/4` and are not
dangerous.  Otherwise, with `d=2u-r`, the bin is `p-1` or `p` exactly when

```text
-1 <= p d/r - x t < 1.                                  (3)
```

In the first branch of (2), `r<=(p-1)/2`.  If `d<=-1`, the middle expression
in (3) is below `-p/r<-2`.  If `d>=1`, then

```text
p/r-x t >= (p-(p-1)t)/r > 1.
```

Thus `d=0`; since `t>1/2`, condition `x t<=1` leaves only
`x=1/2,3/2`.  In the second branch, `x<2`; a negative `d` gives a value below
`-1`, while for positive `d`,

```text
p/r-x t > p/(p-1)-1/(2p) > 1.
```

Again `d=0`, hence `x=1/2` or `3/2`.  These two ratios are precisely (1), and
direct substitution places them in bin `p-1`.

Every edge in (1) lowers `nu_2(j)` by exactly one.  Therefore the danger
graph is acyclic.  If `H` is the nonempty set of half-turn coordinates,
choose a member `j` of minimum `2`-adic valuation.  No outgoing danger edge
from `j` lands back in `H`.  On its protector cell, the selected half-turn
has bin `2p-1`, so it forbids every even shift and no odd shift.  Every other
half-turn forbids no odd shift.  At most `p-2` remaining ordinary coordinates
each forbid one odd shift, leaving at least two odd safe shifts.  This proves
Theorem 2.

For `p=7`, the general graph specializes to the following concrete table:

```text
chosen j | bins at coordinates (2,4,6,8,10,12) | other odd-danger sites
---------+---------------------------------------+------------------------
2        | (13,12,12,11,11,10)                  | none
4        | (6,13,6,12,5,12)                     | 2,6
6        | (4,8,13,3,8,12)                      | none
8        | (3,6,10,13,3,6)                      | 4,12
10       | (2,5,8,11,13,2)                      | none
12       | (2,4,6,9,11,13)                      | 6
```

One reverse-topological sink-elimination priority useful for a constructive
`n=14` witness is

```text
2, 6, 10, 4, 12, 8.
```

The mechanism is more informative than a union bound: the half-turn's large
gcd cost `p` is neutralized by choosing a cell whose odd-shift activation
graph strictly descends in a valuation unavailable to an active minimum.

## Exact gain over weighted fragility

Inside the right-boundary blocker family, the strict cost criterion of
THM-3317 already exposes `1,953,510` nonzero vectors.  Exactly `749,216`
nonzero boundary blockers have cost at least `14`; every one contains a
half-turn.  The activation-table argument closes all `749,216` remaining
cases; the general valuation graph explains why the table works.  Thus the
result is a genuine composite/torsion extension of
THM-3317, not a recount of its certified region.

There is also a useful zero/half-turn duality in this flat:

```text
boundary blindness  <=> at least one coordinate equals 0,
new proof pressure  <=> at least one coordinate equals 7.
```

Both hypersurface unions have cardinality `14^6-13^6`.  The boundary
projection loses which side of this duality the interior parity test sees.

## Independent finite-exact audit

The companion script reconstructs all `812` atomic open-cell patterns at
`n=14`.  The fixed odd zero coordinates leave `210` potentially safe cells,
or `2,940` cell-shift candidates.  It forms one safe-candidate bit mask for
each even coordinate and residue, intersects the first and last three
coordinates separately (`14^3=2,744` masks on each side), and checks all

```text
2,744^2 = 7,529,536
```

vector pairs.  Exactly the zero tuple has empty safe-candidate intersection.
An independent boundary-only meet-in-the-middle pass finds exactly
`2,702,727` blockers, agreeing with Theorem 1.  Named controls retain:

```text
w_2=1:  336 full-cell safe candidates but all 196 boundaries blocked;
w_6=7:   56 full-cell safe candidates but all 196 boundaries blocked;
w=0:      0 safe candidates, the scalar positive control.
```

It separately reconstructs the general edge formula (1) for all `1,134`
unique-protector rows at the `25` odd primes `3<=p<=101`; all `940` observed
edges lower the `2`-adic valuation and the maximum out-degree is two.  This is
a finite control for the symbolic proof, not the source of its all-prime
quantifier.

Reproduce from the repository root with

```text
python 04-computation/lrc_composite_even_flat_crt_interior_collapse_20260803.py
```

The frozen output is
`05-knowledge/results/lrc_composite_even_flat_crt_interior_collapse_20260803.out`.

## Connection passport and limits

- **Source:** THM-3316 prime right-boundary interpolation.
- **Target:** composite `n=2p` right-boundary blockers and their full-cell
  interior refinement.
- **Map:** replace field inversion by the endpoint-owner hypergraph
  `a -> {i:ai in {0,-1}}` and split it by CRT gcd class.
- **Preserved predicate:** whether a zero coordinate blocks every shift at a
  fixed boundary parameter.
- **Destroyed information:** all interior-bin values and shift parity away
  from the boundary.
- **Needed sidecar:** the half-turn odd-shift activation table on
  unique-protector cells.
- **Cheapest decisive test:** the two `14^3` safe-mask families above.

The result closes one precise composite blind flat in the finite residue
model.  It does not show that an arbitrary speed set induces a vector in this
flat or transfer the safe cell to physical time.  Those missing realization,
width, and endpoint coordinates remain necessary before any LRC consequence.

This supplies a proved family of instances of the leak grammar in HYP-1832:
after the CRT boundary quotient creates a large false-positive flat, the sole
large-cost residue `p` is forced into an odd/even activation sidecar and then
leaks through a valuation-descending cell choice.  It also resembles the parity split in
HYP-2077, but the labels must not be conflated: `E` uses **even coordinate
indices** in a residue vector, whereas HYP-2077 uses **even physical speeds**.
A real connection would require a speed-to-residue map preserving that parity
label; no such lift is asserted here.
