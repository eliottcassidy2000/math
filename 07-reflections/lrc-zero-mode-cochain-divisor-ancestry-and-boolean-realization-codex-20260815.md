# Zero-mode-cochain ancestry is divisor pullback plus a Boolean gcd gate

**Date:** 2026-08-15
**Status:** PROVED-ANALYTIC reduction, conditional only on proved
[THM-3405](../01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md);
PROVED-ELEMENTARY universal rank floor four, primitive rank four exactly at
half-twist quotients eight and nine, and global rank four iff the sheet degree
is divisible by eight or nine; [THM-3415](../01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md)
proves, downstream of THM-3408's COMPUTER-ASSISTED PROVED cutoff, global rank
five iff `(10|q or 12|q)` while `8` and `9` do not divide `q`;
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
| rank-five atom | anchor plus four disjoint petals | OR=XOR on sheets | divisor pullback | prime-breaker sidecar |

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

For at least two selected residue types `r_1,...,r_s`--in particular, for
every possible full transverse cover--positive lifts
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

### Universal atom floor: three owners never suffice

The rank-four phenomenon is not merely the first value found in the table.
It is universal.

For a primitive family let

```text
m_i=Q/gcd(Q,v_i)                                      (7a)
```

be the quotient order of owner `v_i`.  Primitive gcd one implies

```text
lcm(m_1,...,m_s)=Q.                                   (7b)
```

The converse for fixed literal lifts is false and is not used; for example,
`Q=5,V=(2,4)` has lcm five but owner gcd two.  Exact residue-class
liftability is supplied instead by the augmented gate `(5)`.

At zero twist an order-`m` owner sees the complete grid `j/m`, so its exact
block size is

```text
(Q/m) z(m),       z(m)=1+2 floor((m-1)/14).           (7c)
```

Every such block contains sheet zero.  For `m>=3`,
`z(m)/m<=1/3`.  If no owner has order two, three block masses sum to at most
`Q`, and their forced common sheet makes their union at most `Q-2`.

All order-two zero blocks are the same half-sized subgroup.  If at least two
owners have order two, duplicates leave at most that subgroup plus one block
of relative size `1/3`, still short of a cover.  It remains to consider
exactly one order-two owner.  If both remaining orders are at least four,
`z(m)/m<=1/4`, with equality only at four.  If one order is three, the other
must satisfy `z(n)/n>=1/6`.  Writing `z(n)=2k+1` on
`14k+1<=n<=14k+14`, this inequality gives `n<=12k+6`, hence `k<=2`.
The complete critical list is therefore finite:

| remaining orders `(m,n)` | `Q=lcm(2,m,n)` | excess of summed capacities over `Q` |
|---|---:|---:|
| `(3,3)` | 6 | 1 |
| `(3,4)` | 12 | 1 |
| `(3,5)` | 30 | 1 |
| `(3,6)` | 6 | 0 |
| `(3,15)` | 30 | 1 |
| `(3,16)` | 48 | 1 |
| `(3,17)` | 102 | 1 |
| `(3,18)` | 18 | 0 |
| `(3,29)` | 174 | 1 |
| `(3,30)` | 30 | 0 |
| `(4,4)` | 4 | 0 |

All three blocks contain zero, so inclusion--exclusion removes at least two
from these summed capacities.  Even the excess-one cases cover at most
`Q-1` sheets.

At half twist an order-two block is empty.  An order-`m>=3` block has at most
`ceil(m/7)` quotient phases, and

```text
ceil(m/7)/m<=1/3,                                     (7d)
```

with equality only at `m=3`.  Thus three owners could meet the mass bound
only if all have order three.  Equation `(7b)` then gives `Q=3`, where every
nonempty half-twist order-three block is the same singleton.  It cannot
cover.  We have proved

```text
rho_ZMC(q)>=4             for every q>=2.              (7e)
```

This is the exact strong-atom floor suggested by the degree-spectrum grammar:
the first zero-cochain cover is a four-owner object.  Its mechanism is forced
overlap plus quotient order, not tournament orientation.

### Complete classification of the rank-four atoms

The same order calculus closes the equality case.  First consider zero twist.
In a minimum cover duplicate blocks are useless.  The order-three and
order-four zero blocks are each unique.  If there is no order-two owner and
no order-three owner, all four densities are at most `1/4`; their common zero
prevents equality.  If there is one order-three owner, the other three must
contribute at least `2/3`.  But apart from the unique order-four block every
order `m>=5` has `z(m)/m<=1/5`, and

```text
1/4+1/5+1/5 < 2/3.
```

Thus two order-four blocks would be required, but they are duplicates.

If there is one order-two owner, its block is exactly the even sheets.  On
the odd-sheet coset an order-`m` zero block occupies the relative fraction

```text
o_0(m)=z(m)/m,                                  m odd;
o_0(m)=4 ceil(floor((m-1)/14)/2)/m,             m even. (7f)
```

One has `o_0(m)<=1/3`, with equality only at `m=3`.  Three remaining blocks
could cover all odd sheets only if all have order three; those blocks are the
same kernel coset.  Hence zero twist has no primitive rank-four atom.

At half twist put

```text
a(m)=2 ceil(floor((m-1)/7)/2),
h(m)=a(m),                         m even;
h(m)=max(a(m),z(m)),               m odd.             (7g)
```

This is the exact maximum number of dangerous phases on an order-`m`
quotient.  Apart from order three,

```text
h(m)/m<=1/4, equality only at m=8;
h(m)/m<=2/9 away from m=8, equality only at m=9.      (7h)
```

The nonempty order-three half block is unique.  With no order-three owner,
four-cover mass forces four order-eight blocks, and they give the primitive
`Q=8` partition.  With one order-three owner, the other three need total
density `2/3`.  Two order-eight blocks force the third order into

```text
{5,8,9,10,11,12,15,17,22,23,24,29,36};
```

one order-eight block leaves only the pairs

```text
(5,9), (9,9), (9,10), (9,15);
```

and no order-eight block forces `(9,9,9)`.  These are exactly eighteen order
profiles.  An exact augmented-mask check of all `978` realizations in them,
including every prime-breaker bit, leaves one positive profile:

```text
(3,9,9,9) at Q=9, with residues (1,5,6,7).            (7i)
```

Together with the separately checked all-order-eight partition `(1,3,5,7)`,
for `979` tests in total, this proves

```text
rho_epsilon^prim(Q)=4
 iff (Q,epsilon)=(8,1) or (9,1).                      (7j)
```

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

### Exact rank-four classification and harmonic pullbacks

Equations `(7)`, `(7e)`, and `(7j)` give the complete global classification:

```text
rho_ZMC(q)=4                iff 8|q or 9|q.             (9a)
```

Its full support is the union of sixteen residue classes modulo 72, so

```text
#{q<=N:8|q or 9|q}=(2/9)N+O(1),
sum_(q<=N, 8|q or 9|q) 1/q=(2/9)log N+O(1).           (9b)
```

Restoring ancestry turns two primitive four-owner atoms into a positive
density subset of the harmonic series.

For the Berggren `U`-spine denominator

```text
Q_n=4n^2+12n+11,
```

direct reduction modulo nine gives `9|Q_n` exactly for
`n==1,5 mod 9`.  Hence those two ninths of the spine have exact zero-cochain
rank four.  They simultaneously lie in the previously proved half-grid
rank-three branch, making the one-unit sidecar tariff explicit.

The Fibonacci sequence modulo eight has period twelve with zeros at indices
`0,6`; modulo nine it has period twenty-four with zeros at `0,12`.  Therefore

```text
rho_ZMC(F_n)=4              whenever 6|n, n>=6.        (9c)
```

Here the order-nine pullback is already contained in the order-eight index
support.  The recurrence tree and divisor ancestry agree because Fibonacci
divisibility transports a primitive Boolean cover along a multiplicative
fibre grade.

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

Concurrent [THM-3408](../01-canon/theorems/THM-3408-fixed-zero-additive-order-duality-and-six-core-corridor.md)
independently proves lcm descent for fixed-zero quotient modes and adds an
additive-order fractional dual.  Its result and `(7)` meet exactly on the
zero-twist layer.  Neither subsumes the other: THM-3408 retains stratum
density and a six-core obstruction through `q=20000`, while the present
reduction retains both Boolean twists and the primitive gcd breaker needed
for exact mobile ranks.

Concurrent [THM-3409](../01-canon/theorems/THM-3409-q15-exceptional-edge-positive-cochain-rigidity-and-leakage-tariff.md)
closes the orthogonal positive-drift control.  The exceptional capped q15
edge has exact complete-pair tariffs `(L1,Linf,L2^2)=(50,6,206)` and faithful
tree tariffs `(10,3)`; no realization has zero cochain.  Thus q15 now has an
exact zero-cochain ancestry classification and an exact nonzero-cochain
tariff without identifying their owner edges.

[THM-3410](../01-canon/theorems/THM-3410-projective-cochain-wedge-ray-tree-tariff-and-residue-scalar-hubs.md)
is now **PROVED analytic + PROVED-ELEMENTARY + VERIFIED-EXACT, with independent
audit requested**.  It writes every realized affine cochain as the integral
wedge `P=A wedge u`, contracts equal primitive rays to an exact tariff MST, and
extends the parity ladder to scalar fibres over the ternary and five-colour
half-grid partitions.  Combining its `a=1` hubs with `(9a)` prices the exact
rank loss on the entire zero-cochain rank-four support:

```text
8|q or 9|q, q even:  zero-cochain rank 4 -> half-grid rank 2,
                      (tau_1,tau_infinity)=(q^2/2,q^2/2);
9|q, q odd:           zero-cochain rank 4 -> half-grid rank 3,
                      (tau_1,tau_infinity)=(4q^2/9,2q^2/9).          (12)
```

These are comparisons of two nested cover predicates, not claims that one
owner edge realizes both minima.  The projective lift identifies the exact
integer sidecar erased by the half-grid quotient; it still gives no LRC row
embedding or ledger decrement.

[THM-3415](../01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md)
closes the next global grade:

```text
rho_ZMC(q)=5 iff (10|q or 12|q) and 8 not|q and 9 not|q.
```

Its atoms at `Q=10,12` are disjoint pointed partitions, so OR and XOR agree
on sheet masks after retaining one anchor plus four petals; no augmented
prime-breaker XOR statement is made.  The support is
32 classes modulo 360, with natural and harmonic density `4/45`; among Farey
vertices its denominator density is `43/864`.  Fibonacci indices are exactly
`15 mod 30`, while the odd Berggren `2c+1` labels never have rank five.  This
is the precise ancestry/Boolean/recurrence transport sought here; it remains
a cover-clutter theorem, not a tournament or an LRC exit.

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
as exact divisor pullbacks.  The floor audit checks the zero-layer order
formula through order 500, freezes all eleven capacity-critical order pairs,
checks all 979 primitive rank-four profile realizations, replays 111 exact
rank-four degrees through `q=500`, and independently recovers the
mod-eight/mod-nine Pisano zeros.  There is no floating point or
`assert`-dependent truth gate.

An independently written rare-coordinate branch-and-bound verifier then
checked both twists for every `2<=Q<=500`, tracking literal gcd transitions
instead of using the union-state BFS or subset enumeration above.  It examined
`374,250` raw types (`371,760` nonempty; `184,338` distinct mask/gcd-transition
types), found no cover of rank at most three, and found exactly the two
primitive rank-four half-twist positives `Q=8` with residues `(1,3,5,7)` and
`Q=9` with residues `(1,5,6,7)`.  The frozen rank-three and rank-four digests
are respectively
`e5b660100f95d7c41c9a6460a42b2283435c1020d8b74e0e0a36a48e0d79f82b`
and
`c55882988ee546925d84844b55980a93e676541257077c48ac2af8083d92fd08`;
endpoint controls `Q=7` half twist and `Q=14` zero twist fail, while adjacent
`Q=8` half twist succeeds.  The audit also exposed and motivated the two
scope repairs in MISTAKE-390; neither repair changes the cover theorem or its
rank consequences.

After THM-3415 the highest-value continuations are:

1. classify the primitive rank-five pattern beyond the finite `Q<=500`
   audit, where the only observed families are `Q=16,18` at zero twist and
   the `8/9` breaker lifts plus atoms `10,12` at half twist;
2. combine the PROVED THM-3414 fixed-zero six-owner bases
   `{15,16,18,20,24}` with a half-twist six-owner classification to decide
   global rank six; THM-3414 inherits THM-3408's computer-assisted cutoff and
   is not an independent proof of that arithmetic gate;
3. intersect the q23 rank-six primitive half-twist with the reserved
   exceptional-edge leakage problem, keeping cover rank distinct from LRC
   row exclusion;
4. transport one zero-cochain certificate through the actual LRC body/core
   sidecars.  Formula `(7)` alone closes no row and leaves LRC(14) open.
