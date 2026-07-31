# J6 recursive flag and parity-ladder audit

**Status:** scratch finite-exact/structural audit.  The first audit covered
the three rank-one hostile carriers; the second, independent full-pipeline
audit now covers the complete fixed-order four-root battery: `73` suffix
branches, `25` scalar-open branches, and `784` H4-pair residuals.  This is
still not a uniform `j=6` theorem and not LRC(14).

## 1. General singleton-complement corollary

Apply `THM-2893-complement-cap-finite-core-flag-lemma.md` to a current
nonempty carrier `C` of mass `h`, with `p` remaining covering labels, using

```text
k=p,       s=p-1,       ell=2,       B=B1,
```

where `B1` is a uniform cap on every allowed singleton coverage.  The
complement block has size `k-s=1`.  THM-2893's strict finite-core condition
is exactly

```text
B1 < (7-s)h/7 = (8-p)h/7.                              (1)
```

Put

```text
theta=h-B1,
H_(p-1)={w:c_C(w)>=theta/(p-1)}.                        (2)
```

If a `p`-set covers `C`, THM-2893 gives

```text
|K intersect H_(p-1)| >= k-s+1 = 2.                    (3)
```

Therefore every hypothetical `p`-cover contains a pair `L` in the finite
core whose literal residual must be covered by the remaining `p-2` labels.
For `p>=4`, `ell=2<s=p-1`, so the flag condition on `L` is vacuous: **all**
pairs in `H_(p-1)` must be checked.  For `p=3`, `ell=s=2`, so only heavy
pairs are admissible.

This is the exact one-step descent

```text
p -> p-2.                                               (4)
```

Iteration gives an odd/even proof ladder only if `(1)` is re-established
on every new nonempty literal residual, with its new mass, allowed-label
set, singleton cap, and high core.  One may not reuse the parent `B1` or
`H`.  Empty residuals are cover witnesses, not successful recursive rows.
The even terminal `p=2` needs its own direct pair gate (or the separate
`(k,s)=(2,1)` specialization); the odd terminal `p=1` is singleton
containment.

### Equality and the `p=7` wall

In abstract THM-2893, `(1)` is a sufficient finite-core condition and is
uniformly sharp.  For the rational-interval LRC carriers it is also an exact
finiteness boundary.

Let `D` be a common denominator of all interval endpoints of `C`.  For every
positive multiple `w` of `D`, each scaled endpoint is integral.  The exact
periodic tooth primitive then gives

```text
c_C(w)=h/7.                                             (5)
```

Hence:

- if the threshold in `(2)` is greater than `h/7`, the discrepancy estimate
  makes `H_(p-1)` finite;
- if it equals `h/7`, infinitely many common-denominator multiples lie in
  `H_(p-1)` by `(5)`;
- if it is below `h/7`, convergence `c_C(w)->h/7` makes the core infinite.

Thus the strict sign in `(1)` is load-bearing, not a cutoff convention.
At `p=7`, `(1)` becomes

```text
B1<h/7,
```

which is impossible because `(5)` forces every uniform singleton cap to be
at least `h/7`.  For `p>=8` the right side is nonpositive.  The
singleton-complement parity ladder therefore has a structural ceiling at
six remaining labels.

## 2. First j6 flag: `(k,s,ell)=(5,4,2)`

On each of the three rank-one hostile carriers the exact global singleton
maximum satisfies

```text
B1=q1<3h/7.                                             (6)
```

This is precisely `(1)` for `p=5`.  The core

```text
H4={w:c_C(w)>=(h-q1)/4}
```

has sizes `16,12,16`.  Because `ell=2<s=4`, all

```text
C(16,2)+C(12,2)+C(16,2)=306
```

pairs must be retained; calling them “heavy pairs” would be a quantifier
error.

For each pair, literal subtraction leaves a three-slot residual.  Two
incomparable sufficient tests are applied:

```text
q1+q2+q3 < residual mass,
B2+q1     < residual mass,                              (7)
```

where the singleton ranks and pair cap are globally tail-sealed.  Equality
in either line is not a closure.  Their union closes `301/306` pairs.  Five
literal residuals survive both tests.

## 3. Second flag: `(k,s,ell)=(3,2,2)`

For each of the five survivors, its own singleton maximum satisfies

```text
q1<5m/7,                                                (8)
```

which is `(1)` for `p=3`.  The finite `H2` cores have sizes

```text
(3,3,3,2,2).
```

Now `ell=s=2`, so the pair must be heavy.  Exact literal union tests find
one heavy edge in each core.  Subtracting it leaves a nonempty residual
which a final singleton would have to cover.

The longest-interval tooth lemma gives exact geometric horizons

```text
(41,38,32,31,33).
```

Equality at the horizon is scanned.  After excluding the inherited prefix,
the first `H4` pair, and the heavy `H2` edge, there are respectively

```text
(23,20,15,13,15)
```

finite singleton checks, `86` total.  Exact interval/tooth containment
agrees with exact scalar full coverage on all `86`; there are zero covers.
All five recursive rows close.

The independent verifier is

```text
.scratch/lrc_thm2892_structural_compression_20260729/
  j6_recursive_geometric_audit.py
```

and gives byte-identical normal and optimized output.

## 4. Audit of the concurrent implementation

The abstract quantifiers and equality conventions above are correct:

- first core: `q1<3h/7`, high membership `>=`;
- first pair flag vacuous, so all 306 pairs;
- tests `(7)` strict;
- second core: `q1<5m/7`;
- second pair flag heavy, with union equality included;
- recursive residuals required nonempty;
- geometric horizon equality scanned;
- the full excluded prefix retained at every descent.

In the initial discovery version,

```text
.scratch/lrc_sparse_j6_gate_audit_20260729/
  suffix_h4_closure_probe.py
```

computed the five final horizons from

```text
(99/70)r/(6m),
```

producing `(160,226,223,139,143)` and `796` checks.  It did not yet implement
the advertised longest-interval horizons.  The current corrected source and
the independent verifier above both implement the geometric version,
reproducing `(41,38,32,31,33)` and `86` checks.

## 5. Full fixed-order four-root independent audit

A fresh standalone verifier imports no repository LRC implementation.  It
rebuilds the danger sets by rational interval sweeps and evaluates coverage
with an independently written periodic antiderivative.  It reconstructs all
four root rankings, all suffixes, all high cores, all pair caps, and the
recursive geometry.  Its exact census is

```text
root gate sizes              (19,20,21,13)
suffix branches              73
scalar-closed / scalar-open  48 / 25
H4 pair residuals            784
top-three closed             771
B2+q1 closed                 773
adaptive union closed        779
adaptive failures            5
finite pair evaluations      7551
```

The rank stratification is decisive:

```text
rank one:      pairs 361, top3 348, B2+q1 350, union 356, failures 5
non-rank-one:  pairs 423, top3 423, B2+q1 423, union 423, failures 0.
```

Thus every non-rank-one branch closes under each certificate separately,
and the only five recursive obligations are exactly the earlier rank-one
rows.  Their H2 cores, unique heavy edges, geometric horizons, and `86`
noncontainment checks all reproduce exactly, with zero singleton covers.

The verifier is

```text
.scratch/lrc_thm2892_structural_compression_20260729/
  full_j6_four_root_independent_audit.py
```

Normal and optimized runs are byte-identical.  It has no Python `assert`
nodes.  Frozen hashes at audit time are

```text
script       827eae32b8db64af3833ea8c2536194e809af8a5d52161a954e9ed3c0a8e762e
audit ledger eae3f0a54ea86dfb61a63e7340f5199689bbf8714f63da28c9276f289611bac4
```

## 6. Exact scope and structural object

The finite-exact result closes the four selected roots and all of their `73`
fixed-order suffix branches.  It does not close the other `3428` seven-body
roots, the uniform `j=6` rung, or LRC(14).

As in THM-2892, the intrinsic object is an undirected heavy
hypergraph/flag recursion decorated by literal residuals.  There is no
tournament orientation: every defining union observable is symmetric, and
the pair flags at `p>=4` are vacuous rather than oriented.
