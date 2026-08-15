# LRC(14): the one-coset affine-cochain theorem and the q=7/q=8 transition

**Date:** 2026-08-14  
**Status:** SUPERSEDED AS CURRENT STRUCTURAL SOURCE BY THM-3395.  The
independently audited canonical theorem now proves the q<=7 criterion and the
q=5/q=6 literal classifications.  This note remains a synthesis of the sharp
q=7/q=8 boundary, tournament sign shadow, quotient branch transplants, and
harmonic monoid actions.  Nothing here gives a refined-ledger decrement or
proves LRC(14).

## 1. The actual threshold is phase spacing `1/7`

Fix a sheet degree `q` and a transverse speed `u`, so `q` does not divide
`u`.  Put

```text
g=gcd(u,q),                 m=q/g.                     (1)
```

At source time `t`, the phases on the q sheets are

```text
u(t+k/q),                   k in Z/qZ.                 (2)
```

There are `m` distinct values modulo one, equally spaced by `1/m`, and each
is repeated on a coset of size `g`.  The danger arc
`{x: ||x||<1/14}` has length `1/7`.  Therefore, whenever

```text
m=q/gcd(u,q) <= 7,                                   (3)
```

at most one phase value in `(2)` is dangerous.  Strictness matters at
`m=7`: an open arc of length exactly `1/7` still cannot contain both endpoints
of a spacing interval.

Thus one owner blocks either no sheet or exactly one coset

```text
B_i(k_i)={k in Z/qZ: k == k_i (mod m_i)},             (4)
```

of size `g_i`.  This is automatic for every transverse speed when `2<=q<=7`.
At `q=8`, an odd speed has `m=8`, and one danger arc can first meet two
adjacent phase classes.  The q=7/q=8 boundary is therefore the structural
transition between one-coset clutters and multi-coset blockers.

## 2. The complete affine-coset criterion

Let `U={u_1,...,u_r}` be transverse speeds satisfying `(3)`.  Choose for each
owner a coset representative `k_i`, and suppose the chosen cosets cover the
sheet group:

```text
Union_i B_i(k_i)=Z/qZ.                                 (5)
```

For every oriented pair define the finite affine gap set

```text
P_ij(k_i,k_j)={p in Z:
  p == (k_j-k_i)u_i u_j  (mod q gcd(u_i,u_j)),
  14|p| < q(u_i+u_j)}.                                 (6)
```

Take `p_ji=-p_ij`.  The proposed common-phase certificate is one choice
`p_ij in P_ij` on every edge of `K_r` satisfying, for all distinct `i,j,h`,

```text
u_h p_ij + u_i p_jh + u_j p_hi = 0.                  (7)
```

The criterion is

```text
the selected coset cover fires at one source time
iff the complete affine cochain (6) has zero circulation (7). (8)
```

This simultaneously contains THM-3388's q=3 triangle, THM-3389's typed
q=4 complete cochain, and the q=5/q=6 probe criteria.

## 3. Proof of necessity

Assume one source time `t` fires every selected blocker.  Choose integers
`a_i` so that the owner interval containing `t` has centre

```text
x_i=a_i/u_i-k_i/q.                                    (9)
```

Its radius is `1/(14u_i)`.  Define

```text
p_ij=q u_i u_j(x_i-x_j)
    =q(u_j a_i-u_i a_j)+(k_j-k_i)u_i u_j.             (10)
```

The first form gives the congruence in `(6)`.  Since the two open intervals
both contain `t`, their centres satisfy

```text
|x_i-x_j| < 1/(14u_i)+1/(14u_j),                     (11)
```

which is exactly the strict inequality in `(6)`.  Finally the normalized
gaps

```text
delta_ij=p_ij/(q u_i u_j)=x_i-x_j                    (12)
```

sum to zero around every triangle; clearing denominators gives `(7)`.

## 4. Proof of sufficiency

Suppose `(5)`--`(7)` hold.  Equation `(7)` says the complete rational
one-cochain `(12)` is closed.  On a complete graph, closedness on every
triangle is equivalent to exactness, so there are rational potentials `z_i`
with

```text
z_i-z_j=delta_ij.                                     (13)
```

We need one common shift `s` such that every shifted potential is a lawful
owner centre:

```text
z_i+s in (1/u_i)Z-k_i/q.                              (14)
```

For a pair `i,j`, the two affine lattices in `(14)` intersect exactly when

```text
p_ij == (k_j-k_i)u_i u_j (mod q gcd(u_i,u_j)),        (15)
```

which is already in `(6)`.  Choose a common multiple of all denominators and
of all `u_i`; after scaling, `(14)` is an ordinary system of integer
congruences.  Pairwise compatibility is sufficient by the generalized CRT.
Hence one `s` gives lawful real centres

```text
x_i=z_i+s.                                             (16)
```

The second part of `(6)` makes the corresponding open real intervals
pairwise intersect.  Ordinary one-dimensional Helly now gives a common
`t` in all of them.  Owner `i` fires on the representative `k_i`; every
sheet in its coset `(4)` has the same phase modulo one.  The coset cover `(5)`
therefore fires every sheet.  This proves `(8)`.

Notice what is and is not used.  The CRT/cochain converse works for any finite
number of owners.  The bound `(3)` is used only to make one owner's sheet set
a single coset.  In an irredundant cover each owner has a private sheet, so
there are automatically at most `q` selected owners.

## 5. The whole q=2 through q=7 atlas has one grammar

For the literal six-speed LRC body, the one-coset theorem produces the
following hierarchy.

| q | blocker sizes | minimal-edge ranks | literal status | exact atlas rows | core rescues |
|---:|---|---|---|---:|---:|
| 2 | 1 | 2 | PROVED, THM-3387 gcd graph | 252 | 0 |
| 3 | 1 | 3 | PROVED, THM-3388 | 588 | 3 |
| 4 | 1,2 | 2,3,4 | PROVED, THM-3389 | 619 | 0 |
| 5 | 1 | 5 | PROVED, THM-3395; independent q=5 convergence audit | 1,619 | 2 |
| 6 | 1,2,3 | 3,4,5 | PROVED, THM-3395; independent q=6 convergence audit | 1,478 | 7 |
| 7 | 1 | 7 | PROVED, THM-3395 | 2,079 | 0 |

At prime `q<=7`, every blocker is a singleton and a minimal cover has rank
`q`.  At q=7 a body has at most six transverse owners, so no cover edge can
exist; every candidate row survives this sheet obstruction.  Composite q=4
and q=6 lower the rank by allowing subgroup fibres of sizes two and three.

The twelve core rescues through q=7 occur only at q=3,5,6.  They are not
errors in the global clutter: the transverse speeds do cover every sheet, but
their common-phase cells lie entirely inside the divided core danger union.

## 6. Tournaments with missing and bidirected edges are the sign shadow

For a fixed pair `(i,j)`, project the affine set `(6)` to its available signs.
Ignoring a possible zero value gives four coarse states:

```text
empty:       no pair overlap;
i -> j:      only positive gaps;
j -> i:      only negative gaps;
both:        both signs occur.                         (17)
```

This is exactly a directed graph with missing and bidirected edges.  Its
four-state pair alphabet is a Boolean square, the same cardinality as XOR
two-bit data.  It is a lawful quotient of the affine problem, but not a
faithful one.

If a nonzero complete cochain exists, `sign(p_ij)` orders the real centres
`x_i`; its tournament is transitive.  A directed cycle is therefore an
immediate obstruction.  The converse fails because transitivity remembers
only signs, while `(6)`--`(7)` require exact integral magnitudes and affine
residue classes.

The hostile ladder makes the information loss explicit:

```text
q=3: (1,4,7);
q=4: (2,7,11) and (1,3,11,5);
q=5: (1,4,2,3,9);
q=6: (2,8,14)=2(1,4,7).                               (18)
```

Each displayed assignment has every required pair gap, but no global exact
cochain.  Thus a tournament is a useful obstruction and serialization
sidecar, never an iff certificate by itself.

The independent q=5 audit sharpened its hostile: the star gaps are forced to
`(-1,-1,-1,1)`, which force `p_14=13` although the legal fibre is only
`{-2,3}`.  The robust order `(13,16,18,19,17)` has pair-fibre cardinalities
`(4,4,5,4,3,5,5,6,5,5)` and still has no complete cochain.  It also proves a
q=5 boundary simplification: distinct sheet labels can never be tangent,
because tangency forces `5|p` while the affine residue is nonzero modulo five.

This also explains why q=4 sheets and a four-state directed pair must not be
identified merely by cardinality.  The sheet carrier is the cyclic group
`Z/4Z`; the pair-state carrier in `(17)` is a Boolean square.  Their subgroup
and quotient structures differ.

## 7. An exact q=3 to q=6 branch transplant

THM-3388 uses the q=3 cover edge

```text
E_3={1,4,5}.                                           (19)
```

Doubling gives the q=6 pure pair-blocker edge

```text
E_6={2,8,10}=2E_3.                                    (20)
```

Under `Z/6Z -> Z/3Z`, every size-two blocker collapses to one q=3 sheet.
The source-time map is `s=2t`:

```text
(2u)(t+k/6)=u(s+k/3).                                 (21)
```

Hence q=6 cover witnesses for `(20)` are exactly lifts of q=3 witnesses for
`(19)`.  This is not an analogy; it is an iff transport preserving common
phase after the explicit source rescaling.  It destroys the binary coordinate
inside each antipodal pair and needs no further sidecar on this pure stratum.

Now apply the commuting multipliers `7,11,13`.  At every depth, the q=6
ternary ancestry is exactly twice the q=3 ancestry.  The free word counts and
exponent support agree, while harmonic mass is halved:

```text
mass(E_3 orbit)=29029/14400,
mass(E_6 orbit)=29029/28800.                           (22)
```

Both distinct-support shell masses satisfy

```text
1001H_d=311H_(d-1)-31H_(d-2)+H_(d-3),  d>=3,          (23)
```

because the recurrence belongs to the multiplier monoid, not to the root
edge.  The q=5 positive edge `{1,2,3,4,8}` supplies another orbit with root
mass `53/24` and total harmonic mass

```text
(53/24)(1-1/7)^(-1)(1-1/11)^(-1)(1-1/13)^(-1)
=53053/17280.                                          (24)
```

Unique prime factorization identifies each collision-quotiented support with
a subset of the harmonic series.  Every subfamily has finite mass because the
full product in `(22)` or `(24)` converges.  The free ternary word tree is a
different representation: commuting words with the same exponent triple
collide, so word subsets require multiplicity or a Parikh-vector quotient.

## 8. Consequences and next hostile tests

The one-coset theorem replaces five isolated fibre stories by one carrier:

```text
coset cover on Z/qZ
  + complete affine integral 1-cochain
  + zero triangle circulation
  + core-cell sidecar.                                 (25)
```

The next decisive work is therefore concrete.

1. Use the completed THM-3395 and q=5/q=6 audits as the pinned base for the
   q=8 mode-expanded theorem.
2. Freeze the exact common-phase cells of all twelve core rescues.  Preliminary
   exact geometry shows repeated translated two-collar atoms rather than
   twelve unrelated accidents.
3. Formulate the q=8 replacement carrier: one odd owner may block two adjacent
   singleton sheets, so a blocker is now a short cyclic interval of cosets,
   not one coset.
4. Transport only the pure quotient strata, such as `(19)`--`(21)`, to the
   ancestry/current frontier.  Mixed q=6 patterns lose row/column incidence
   under separate q=2 and q=3 projections and require the full q=6 cochain.
5. Keep the Jacobian degree-monoid analogy at grammar level until an explicit
   map is supplied.  Common dilation acts within one sheet grade; polynomial
   composition multiplies Keller degrees.  Those are structurally parallel
   operations, not yet the same theorem.
