---
id: THM-3317
title: "Unique-protector cells and weighted scalar fragility"
status: >
  PROVED + FINITE-EXACT AUDITED.  For every n>=3 and every coordinate, an
  explicit open micro-staircase cell has that coordinate as its unique
  zero-ramp endpoint protector.  A perturbation of any scalar ramp is
  therefore not a full blocker whenever the sum of its exact gcd shift-costs
  is less than n.
source: root/creative-jacobian-lrc/2026-08-03
depends_on:
  - THM-363-lrc-scalar-gauge-reindexing
related:
  - THM-364-lrc-scalar-ramp-cell-blocking
  - THM-3316-prime-right-boundary-interpolation-forces-scalar-rigidity
scripts:
  - 04-computation/lrc_microstaircase_unique_protector_sparse_fragility_opus_20260803.py
outputs:
  - 05-knowledge/results/lrc_microstaircase_unique_protector_sparse_fragility_opus_20260803.out
---

# THM-3317 -- unique-protector cells and weighted scalar fragility

**PROVED + FINITE-EXACT AUDITED.**

## Statement

Fix `n>=3`.  For every `1<=j<=n-1`, define

```text
             1/2+1/(4n^2)     if 2j<=n-1,
t_(n,j) =
             1/(2n)           if 2j>n-1,

alpha_(n,j)=(1-t_(n,j)/n)/j.                              (1)
```

Then `alpha_(n,j)` lies in an open cell of the full arrangement and

```text
floor(n {i alpha_(n,j)}) in {0,n-1}    iff    i=j.        (2)
```

Thus `j` is the unique endpoint protector of the zero ramp on this cell.

For a nonzero residue `a mod n`, put

```text
kappa_n(a)=2                  if gcd(a,n)=1,
           gcd(a,n)           otherwise.                 (3)
```

Let `v in (Z/nZ)^(n-1)`, choose any scalar ramp `mi`, and set

```text
w_i=v_i-mi,               J={i:w_i!=0}.                  (4)
```

If `J` is nonempty and

```text
sum_(i in J) kappa_n(w_i)<n,                              (5)
```

then `v` is not a full micro-staircase blocker.

## Proof of the explicit cell

Write `i=qj+r`, where `0<=r<j`.  Equation `(1)` gives

```text
i alpha_(n,j)=q+r/j-i t_(n,j)/(nj).                      (6)
```

At `i=j`, the fractional part is `1-t/n`, strictly inside the last bin.

First suppose `2j<=n-1`.  If `i=qj` is a different multiple, then `q>=2`,
`qt>1`, and `qt<n-1`; hence `n{ialpha}=n-qt` is in an interior bin.  If
`r>=1`, first set `t=1/2`.  The margin above the lower endpoint `1/n` has
numerator

```text
2nr-i-2j >= 2n-(n-1)-(n-1)=2.                           (7)
```

Increasing `t` by `1/(4n^2)` moves the fractional part by less than this
strict margin.  The upper bound follows from
`r/j<=1-1/j<1-1/n`.  Thus every `i!=j` is interior.

Now suppose `2j>n-1`.  There is no different multiple of `j` among
`1,...,n-1`.  At `t=0`, every other fractional part has distance at least
`1/j>1/n` from an integer.  The perturbation `t=1/(2n)` moves it downward by
less than half its smallest lower margin and away from the upper endpoint.
All inequalities are strict.  This proves `(2)` and also shows that `(1)` is
not a breakpoint.

## Proof of weighted fragility

By THM-363, subtracting a scalar ramp preserves the full-blocker predicate,
so it is enough to treat `w`.  Choose any `j in J` and the cell `(1)`.
Coordinates outside `J` have residue zero and, by `(2)`, an interior bin; they
block no shift.

For `i in J`, if its cell bin is `b_i`, its bad shifts solve

```text
s w_i=-b_i             or             s w_i=-b_i-1 mod n. (8)
```

Let `g=gcd(w_i,n)`.  If `g=1`, each congruence has one solution and their
solutions are distinct.  If `g>1`, at most one of two consecutive right-hand
sides is divisible by `g`; that congruence, when soluble, has exactly `g`
solutions.  Coordinate `i` therefore forbids at most `kappa_n(w_i)` shifts.
The union bound and `(5)` leave a shift forbidden by no coordinate.  Together
with the cell `(1)`, it is a safe candidate, proving the theorem.

## Consequences and failure boundary

- Every one-coordinate deformation of every scalar ramp is safe for every
  `n>=3`.
- For prime `p`, the theorem covers every perturbation supported on at most
  `(p-1)/2` coordinates.  THM-3316 is stronger at primes and closes all
  support layers by a different boundary-interpolation mechanism.
- At `n=14`, all nonzero residues except the half-turn `7` have cost `2`,
  while `kappa_14(7)=7`.  Thus `(5)` covers six ordinary defects or one
  half-turn plus three ordinary defects.

The strict inequality is not necessary.  Two half-turn defects at coordinates
`2,3` for `n=14` have cost exactly `14`, yet exact checking finds the safe
candidate at cell `254`, shift `1`, with pattern

```text
(4,8,13,3,7,12,2,7,11,1,6,10,1).                       (9)
```

This hostile prevents interpreting the weight sum as an iff criterion.

## Verification record

The companion checks every explicit cell for `3<=n<=30`, every congruence
cardinality for `3<=n<=60`, and every vector satisfying `(5)` after the exact
gauge `w_1=0` for `3<=n<=8`.  It audits scalar positive controls, the equality
hostile `(9)`, the independent prime critical-layer moment obstruction, and
dynamic-programming class counts at `n=13,14,15`.  Normal and optimized runs
equal the frozen transcript.

```text
script sha256 (LF)  4ab6544d589f90263ecbeb4fc39b175d9cb788d406de911c76deeb45d87b96f0
output sha256 (LF)  f47f93daae5e9ebef79098e133fe00cc2ce6677ddcaa0d8a320166d01950eb9f
```

This theorem retains only residue support and gcd shift-costs.  It discards
cell width, prime-grid realization, endpoints, divisibility branches, and
physical time, so it is not an LRC lift.
