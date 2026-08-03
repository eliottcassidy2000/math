---
id: THM-3325
title: "Composite even-flat boundary blindness and valuation-descending interior collapse"
status: >
  PROVED + FINITE-EXACT HOSTILE-AUDITED.  For every odd prime p, the
  normalized micro-staircase flat modulo 2p supported on the p-1 nonzero even
  coordinate indices has (2p)^(p-1)-(2p-1)^(p-1) right-boundary blockers, but
  its only full open-cell blocker is zero.  The boundary projection is blind
  exactly when at least one even coordinate vanishes.  Interior unique-protector cells
  recover the missing parity sidecar: half-turn danger arrows are precisely
  j->j/2 and j->3j/2 when those targets are even and in range, so every arrow
  lowers the 2-adic valuation and an active minimum exposes a safe odd shift.
audit: >
  The exact companion reconstructs all 812 n=14 atomic patterns, audits the
  complete 14^6 even flat by two independent 14^3 meet-in-the-middle mask
  families, and finds only zero among 7,529,536 full blockers.  A separate
  boundary pass finds exactly 2,702,727 blockers.  It reconstructs the
  general danger graph on 1,134 protector rows for every odd prime through
  101; all 940 edges lower the 2-adic valuation.  The all-prime quantifier is
  proved symbolically, not inferred from this bank.
source: root/second-creative-jacobian-lrc/2026-08-03
depends_on:
  - THM-3317-unique-protector-cells-and-weighted-scalar-fragility
  - THM-364-lrc-scalar-ramp-cell-blocking
related:
  - THM-3316-prime-right-boundary-interpolation-forces-scalar-rigidity
  - HYP-1832-lrc-torsion-crt-leak-dichotomy
script: 04-computation/lrc_composite_even_flat_crt_interior_collapse_20260803.py
output: 05-knowledge/results/lrc_composite_even_flat_crt_interior_collapse_20260803.out
script_sha256: 848ed70743a97048b18123801f6c4309b7e8d7829b9bc1df31e49088a31bd3c9
output_sha256: b3f46e1abc5bbe0ee02c8731480b3ffb3b9dd23a45d260ebb66622431bf7610d
hash_basis: LF-normalized bytes
---

# THM-3325 -- composite even-flat boundary blindness and valuation-descending interior collapse

**PROVED + FINITE-EXACT HOSTILE-AUDITED.**

## 1. Statement

Let `p` be an odd prime, put `n=2p`, and work in the normalized residue
space `(Z/nZ)^(n-1)`.  Define the even-coordinate flat

```text
E_p={w:w_i=0 for every odd i}.                            (1)
```

For `a,s in Z/nZ`, call the right-boundary candidate `(a,s)` blocked by `w`
when some `1<=i<=n-1` satisfies

```text
a i+s w_i in {0,-1} mod n.                              (2)
```

Then:

1. `w in E_p` blocks all `n^2` right-boundary candidates if and only if at
   least one of `w_2,w_4,...,w_(2p-2)` is zero.  Hence `E_p` contains exactly

   ```text
   (2p)^(p-1)-(2p-1)^(p-1)                              (3)
   ```

   right-boundary blockers.

2. The zero vector is the only member of `E_p` that blocks every shift on
   every open micro-staircase cell.  Every nonzero `w in E_p` has an explicit
   safe cell and an odd shift.

At `p=7`, the boundary count is

```text
14^6-13^6=2,702,727,                                    (4)
```

including `2,702,726` nonscalar false positives, whereas the full-cell count
is one.  This is a finite residue/cell theorem, not a realization theorem for
physical speeds and not a proof of LRC(14).

## 2. Exact right-boundary classification

Suppose first that `w_(2j)=0` for some even coordinate.  Every boundary
parameter `a` has a zero-coordinate endpoint owner:

```text
a=0:                 i=1 gives ai=0;
gcd(a,2p)=1:         i=-a^(-1), which is odd, gives ai=-1;
gcd(a,2p)=2:         i=p gives ai=0;
a=p:                 i=2j gives ai=0.                 (5)
```

These exhaust the gcd classes modulo `2p`, so `(2)` holds for every shift.

Conversely, suppose all `p-1` even coordinates are nonzero and choose `a=p`.
Every odd zero coordinate has interior bin `p`, while every active even
coordinate has bin zero.  A nonzero residue `v` forbids the following shifts:

```text
gcd(v,2p)=1:       {0,-v^(-1)}, with -v^(-1) odd;
gcd(v,2p)=2:       {0,p};
v=p:               all even shifts.                    (6)
```

If no coordinate equals `p`, shift `2` is safe.  If at least one coordinate
equals `p`, those coordinates cover only the even shifts and the remaining
at most `p-2` coordinates cover at most `p-2` of the `p` odd shifts.  An odd
shift remains safe.  This proves the iff statement and `(3)`.

## 3. Odd-shift cost on a unique-protector cell

Fix a nonzero `w in E_p` and an active even coordinate `j`.  Use the explicit
cell `alpha_(2p,j)` of THM-3317.  Its bin at `j` is `2p-1`, and every other
coordinate has an interior bin.  Thus every zero coordinate is harmless.

For a nonzero residue `v!=p` and any bin `b`, the endpoint congruences

```text
b+s v=0 or -1 mod 2p                                   (7)
```

forbid exactly one odd shift.  If `v` is a unit, the two solutions differ by
the odd unit `v^(-1)` and have opposite parity.  If `gcd(v,2p)=2`, exactly one
of the two consecutive right sides is soluble; its two solutions differ by
`p` and again have opposite parity.

If no active coordinate has residue `p`, at most `p-1` ordinary residues
therefore forbid at most `p-1` of the `p` odd shifts.  One odd shift is safe.

## 4. The half-turn danger graph

It remains to treat active half-turn residues `v=p`.  On an interior bin `b`,
such a coordinate forbids every odd shift exactly when

```text
b in {p-1,p};                                           (8)
```

otherwise it forbids no odd shift.  Put an arrow `j->h` between distinct even
indices when the bin at `h` on the `j`-protector cell satisfies `(8)`.

The exact graph law is

```text
j->h  iff  h=j/2 or h=3j/2,
             with h even and 2<=h<=2p-2.                (9)
```

To prove it, write `j=2r`, `h=2q`, `x=q/r=m+u/r`, `0<=u<r`, and use

```text
alpha=(1-t/(2p))/j,
t=1/2+1/(16p^2)   if r<=(p-1)/2,
t=1/(4p)          if r>=(p+1)/2.                        (10)
```

If `u=0`, then the small subtraction in `(10)` wraps an integer ratio to a
bin strictly above `p`; the same is true whenever the subtraction wraps a
positive `u/r`, since `x t<p/2`.  Thus neither case is dangerous.  In the
remaining no-wrap case, the bin lies in `{p-1,p}` exactly when, for
`d=2u-r`,

```text
-1 <= p d/r-x t < 1.                                   (11)
```

In the first branch of `(10)`, a negative `d` puts the middle expression
below `-p/r<-2`, while a positive `d` gives

```text
p/r-x t >= (p-(p-1)t)/r >1.                            (12)
```

In the second branch, `x<2`; a negative `d` gives a value below `-1`, while a
positive `d` gives

```text
p/r-x t > p/(p-1)-1/(2p)>1.                            (13)
```

Thus `d=0`, so the fractional part of `x` is `1/2`.  The bound `xt<=1` in
the first branch, and `x<2` in the second, leave only `x=1/2` and `x=3/2`.
Direct substitution puts both in bin `p-1`, proving `(9)`.

Every arrow in `(9)` lowers `nu_2(j)` by exactly one.  Let `H` be the nonempty
set of active half-turn coordinates and choose `j in H` of minimum `2`-adic
valuation.  No arrow from `j` lands in `H`.  On the `j`-protector cell, the
chosen half-turn has bin `2p-1`, so it forbids all even shifts but no odd
shift; every other half-turn also forbids no odd shift.  At most `p-2`
ordinary coordinates each forbid one odd shift, leaving at least two safe
odd shifts.  This proves that every nonzero vector in `E_p` is safe.

The zero vector is a scalar ramp and is a full blocker by THM-364, completing
the classification.

## 5. Exact audit and failure boundary

At `n=14`, the companion reconstructs all `812` atomic patterns and retains
`2,940` cell-shift candidates after the fixed odd coordinates are filtered.
Two independent `14^3` mask banks inspect all

```text
14^6=7,529,536                                          (14)
```

vectors and find only zero as a full blocker.  A separate boundary-only pass
finds exactly `(4)`.  The named hostile controls are

```text
w_2=1:   336 safe full-cell candidates, 196/196 boundaries blocked;
w_6=7:    56 safe full-cell candidates, 196/196 boundaries blocked;
w=0:       0 safe candidates.                            (15)
```

The script also reconstructs `(9)` for all `25` odd primes through `101`,
checking `1,134` protector rows and `940` arrows.  This is a hostile audit of
the symbolic proof, not the source of its quantifier.

Primality and the restriction to the even-coordinate flat are load-bearing.
For a general odd cofactor, extra gcd classes can forbid more than one odd
shift.  For vectors with active odd indices, the right-boundary owner
hypergraph and the valuation graph change.  Finally, the map from physical
speeds to this residue flat does not preserve cell width, endpoint origin, or
physical time; none is claimed here.

The mechanism gives a proved family instance of HYP-1832's CRT-leak grammar:
the boundary quotient destroys parity, while interior protector cells restore
it through a valuation-descending sidecar.
