# Zero-mode-cochain ancestry is divisor pullback plus a Boolean gcd gate

**Date:** 2026-08-15  
**Status:** PROVED-ANALYTIC reduction, conditional only on proved
[THM-3405](../01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md);
FINITE-EXACT unrestricted positive-transverse rank table for `15<=q<=28`,
with independent union-state and exhaustive-combination solvers and literal
witness replay.  No LRC(14) ledger decrement.

## 1. Inheritance and the corrected target

The closest mechanism is THM-3405: a vanishing complete mode cochain has only
the fixed-zero and half-twist centre gauges after the active owner gcd is
known.  The canonical hostile is MISTAKE-389's `q=15,c=1/150` physical
partition: it satisfies the half-grid equations but fails the mode-gcd gate.
The corrected near miss is therefore the old divisor-chart rank table, now
properly typed as a larger physical slice.  The least-used sidecar is the
**primitive owner gcd after quotienting by the active common scale**.

The live concept board is:

| object | representation | invariant | operation | lost coordinate |
|---|---|---|---|---|
| zero mode cochain | THM-3405 scalar word | Boolean twist `epsilon` | common sheet translation | absolute sheet origin |
| active owners `U` | `U=dV`, `gcd(V)=1` | `g=gcd(q,d)` | divisor pullback | primitive quotient `Q=q/g` |
| primitive cover | sheet bitmask plus prime breakers | gcd-one realizability | set union | literal positive lifts |
| physical half-grid cover | divisor chart | block rank | affine normalization | selected-mode centre/cochain |
| rank-four witness | four-block clutter | full union plus gcd gate | dilation | no intrinsic pair orientation |

The exact connection contract is:

| field | value |
|---|---|
| source | distinct positive transverse owners with zero complete THM-3398 cochain |
| target | primitive fixed/half-centre covers on divisor quotients |
| map | divide owners by their gcd and reduce sheets modulo `Q=q/gcd(q,d)` |
| preserved | strict danger sets, full cover, transversality, Boolean twist, and owner count |
| destroyed | common fibre multiplicity and literal owner labels |
| required sidecar | `gcd(V)=1`, encoded by prime-breaker bits |
| cheapest decisive test | reproduce q15--28 by two finite solvers and replay literal owners |

## 2. Exact divisor-ancestry reduction

Let an active zero-cochain family be

```text
U=dV,       gcd(V)=1,
g=gcd(q,d), q=gQ, d=g d_0, gcd(Q,d_0)=1.              (1)
```

THM-3405 gives, after one common sheet translation, only

```text
c_epsilon=epsilon g/(2qd),       epsilon in {0,1}.    (2)
```

For `u=dv`, substitution gives

```text
u(c_epsilon+ell/q)
 =epsilon v/(2Q)+d_0 v ell/Q.                         (3)
```

Multiplication by `d_0` permutes `Z/QZ`.  Hence each danger set in `(3)` is
the `g`-fold inverse image of the primitive block

```text
B^epsilon_(Q,r)
 ={ell in Z/QZ:
   14 min(x,2Q-x)<2Q,
   x == r(2ell+epsilon) mod 2Q}.                       (4)
```

At zero twist `r` is read modulo `M_0=Q`; at half twist it is read modulo
`M_1=2Q`.  Transversality is exactly `Q` not dividing `r`.  Therefore every
zero-cochain certificate descends to a primitive full cover on `Q`, and every
primitive full cover pulls back along `Z/qZ -> Z/QZ`.

## 3. The gcd-one condition is a finite Boolean realization gate

For selected residue types `r_1,...,r_s`, positive lifts
`v_i == r_i mod M_epsilon` with `gcd(v_1,...,v_s)=1` exist exactly when

```text
gcd(M_epsilon,r_1,...,r_s)=1.                         (5)
```

Necessity is immediate.  For sufficiency, a full strict danger cover uses at
least two nonfull blocks.  Put `R=gcd(r_2,...,r_s)`.  For each prime `p|R`:

- if `p|M_epsilon`, condition `(5)` forces `p` not to divide `r_1`, so every
  lift of `r_1` avoids `p`;
- if `p` does not divide `M_epsilon`, exactly one residue class of the lift
  parameter `k` makes `p | r_1+kM_epsilon`.

Choose one allowed class modulo every latter prime and combine them by CRT.
Then replacing `r_1` by `r_1+kM_epsilon` gives gcd one.  Thus no unbounded
owner search remains.

Computationally, `(5)` is ordinary set cover on an augmented universe:

```text
Q sheet bits
+ one breaker bit beta_p for every prime p|M_epsilon,

r covers beta_p  iff  p does not divide r.             (6)
```

This is the exact Boolean realization promised by the ancestry language.
The prime bits are not decorative: deleting them admits quotient covers whose
every positive lift retains a forbidden common factor.

Define `rho_epsilon^prim(Q)` as the minimum augmented cover rank, or infinity
if none exists.  Equations `(1)--(6)` prove the divisor-minimum formula

```text
rho_ZMC(q)
 =min_(Q|q,Q>=2) min(rho_0^prim(Q),rho_1^prim(Q)).     (7)
```

This is an exact min-convolution on the divisibility poset.  A primitive
quotient `Q` is an ancestry node; multiplying the fibre degree `g=q/Q` moves
within its dilation ray without changing owner count.

## 4. Primitive and global classifications

The exact primitive pairs `(zero,half)` are

| `Q` | `2..7` | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rank pair | `(inf,inf)` | `(inf,4)` | `(inf,4)` | `(inf,5)` | `(inf,6)` | `(inf,5)` | `(inf,7)` | `(inf,7)` |

and

| `Q` | 15 | 16 | 17 | 18 | 19 | 20 | 21 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| rank pair | `(6,6)` | `(5,5)` | `(8,8)` | `(5,5)` | `(9,9)` | `(6,6)` | `(8,8)` |

| `Q` | 22 | 23 | 24 | 25 | 26 | 27 | 28 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| rank pair | `(7,6)` | `(11,6)` | `(6,5)` | `(11,6)` | `(8,7)` | `(10,5)` | `(8,7)` |

Applying `(7)` gives

```text
q:          15 16 17 18 19 20 21 22 23 24 25 26 27 28
rho_ZMC(q):  6  4  8  4  9  5  8  6  6  4  6  7  4  7. (8)
```

The ancestry explains every nonprimitive drop:

```text
16 <- 8,   18 <- 9,   20 <- 10,   22 <- 11,
24 <- 8,   26 <- 13,  27 <- 9,    28 <- 14.           (9)
```

At `q=22,26,28`, the displayed ancestor ties the primitive half layer; it is
not the unique minimizer.  This is precisely the degree-graded monoid view:
the family is not a list of isolated witnesses, but primitive Boolean atoms
plus multiplicative fibre degrees.

One literal half-centre witness in each grade is:

| `q` | ancestor `Q` | owners `U` |
|---:|---:|---|
| 15 | 15 | `(1,4,6,7,8,10)` |
| 16 | 8 | `(2,6,10,14)` |
| 17 | 17 | `(1,3,4,5,7,8,9,11)` |
| 18 | 9 | `(2,10,12,14)` |
| 19 | 19 | `(1,3,4,5,7,8,9,11,12)` |
| 20 | 10 | `(2,6,8,14,18)` |
| 21 | 21 | `(1,4,5,6,8,11,13,14)` |
| 22 | 11 | `(2,4,6,10,14,18)` |
| 23 | 23 | `(1,4,5,7,9,11)` |
| 24 | 8 | `(3,9,15,21)` |
| 25 | 25 | `(1,9,10,11,19,21)` |
| 26 | 13 | `(2,4,6,10,14,18,22)` |
| 27 | 9 | `(3,15,18,21)` |
| 28 | 14 | `(2,6,8,10,18,22,26)` |

All use `c=1/(2q)`.  Their active gcd is `g=q/Q`, the THM-3405 scalar is
`a=g`, every owner centre word is `H_i=u_i`, and
`gcd(q,u_i)|H_i`.  Direct strict danger masks cover all sheets.

## 5. Two elementary lower-bound controls

The finite solver is decisive for the table, but two formerly cap-sensitive
entries also admit short analytic lower bounds.

For `q=25`, the only possible ancestors are `Q=5,25`, and `Q=5` has no
primitive cover.  A primitive `Q=25` family must contain a residue breaking
the prime five.  Such a block has size at most three at zero and four at half
twist; every five-divisible transverse block has size at most five.  Five
owners therefore cover at most `3+4*5=23` or `4+4*5=24` sheets.  Rank six is
necessary and attained.

For `q=27`, the ancestors are `Q=3,9,27`.  Quotient three has no primitive
cover.  A primitive family must contain a three-breaking residue.  With
three owners the mass bounds are

```text
Q=9:  2+3+3=8<9,
Q=27: 4+9+9=22<27.                                  (10)
```

Thus rank four is necessary, and the `Q=9` half-twist pulls back to the
displayed witness.

## 6. Three nested ranks, now cleanly separated

The corrected comparison is

```text
fixed c=0:       6,5,8,5,9,6,8,7,11,6,11,8,10,8
zero cochain:    6,4,8,4,9,5,8,6, 6,4, 6,7, 4,7
physical HG:     3,2,8,2,9,2,3,2, 6,2, 5,2, 3,2.     (11)
```

The inclusions of feasible loci run in the opposite direction to the ranks:

```text
fixed-zero subset zero-cochain subset synchronized half-grid physical.
```

At `q=17,19,23`, the zero-cochain and half-grid minima happen to coincide.
Every other gap is genuine lost mode-centre information.  The infinite even
parity ladder `P=-a q^2/2` supplies a uniform positive-drift hostile, not a
candidate zero-cochain shortcut.

## 7. Why the rank-four object is not a tournament

The primitive reduction really is Boolean, but it is a **cover clutter on
sheet and prime-breaker vertices**.  In the zero-cochain slice every pair
centre gap vanishes, so a pairwise observable has only ties.  Orienting those
ties would add gauge, not information.  XOR of sheet masks also loses the
target: the cover predicate is OR, and XOR cancels overlaps.

There is a useful method-level echo of
[THM-3407](../01-canon/theorems/THM-3407-hadamard-core-multitoggle-response-plaquette-shells-and-trade-distance.md).
There, pair plaquettes lose an oriented triple-minor sidecar.  Here, sheet
blocks lose prime divisibility and mode-centre sidecars.  In both cases the
repair is to augment the quotient by the smallest coordinate that restores
realizability.  There is no object-level theorem transfer.

Thus the recurring “tournament of size four” intuition should be retyped:
four is often the first realizable owner count, but the predicate-preserving
carrier is a four-edge augmented clutter.  Any tournament is at most a
scheduling or visualization layer over that clutter.

## 8. Verification and new frontiers

Run

```text
python 04-computation/lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.py
python -O 04-computation/lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.py
```

The companion pins THM-3405, THM-3401, and the MISTAKE-389 repaired
half-grid artifact; checks all 54 primitive twist banks for `2<=Q<=28`; and
uses two distinct solvers.  The union-state route visits `36,580` states and
`565,480` transitions.  The exhaustive route checks `394,418` subsets.  All
fourteen displayed witnesses are reconstructed as literal strict covers and
as exact divisor pullbacks.  There is no floating point or `assert`-dependent
truth gate.

The highest-value continuations are:

1. prove or refute the apparent universal primitive rank-four floor, then
   turn the `Q=8,9` ancestors into exact infinite rank-four families rather
   than upper bounds;
2. compute the primitive spectrum beyond `Q=28` and classify which new
   ancestors beat all proper divisors;
3. intersect the q23 rank-six primitive half-twist with the reserved
   exceptional-edge leakage problem, keeping cover rank distinct from LRC
   row exclusion;
4. transport one zero-cochain certificate through the actual LRC body/core
   sidecars.  Formula `(7)` alone closes no row and leaves LRC(14) open.
