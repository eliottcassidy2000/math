---
id: THM-4375
title: "Bilateral source-cone reciprocal orbits and fibre imbalance"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4368 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED; JC(2) OPEN. Reciprocal-realizable packet types form exactly a looped
  finite-power path band, with explicit orbit/address generating functions,
  primitive-ray scale intervals, and sharp source-fibre imbalance. This is not
  a tournament and gives no bracket, seam, JC(2), or DC(2) consequence.
source: root + boundary_orbit_count + clean-room referee / next-sharp session, 2026-09-03
depends_on:
  - THM-4368-diagonal-boundary-valuation-triangular-address-and-simplex-stream-rank
related:
  - THM-4369-source-packet-pascal-circuit-kernel-and-boundary-basis
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
  - THM-2352-q-adic-prefix-residue-collision-spectrum
mistake_firewall:
  - MISTAKE-222
primary_script: 04-computation/jc2_bilateral_source_cone_reciprocal_orbits_thm4375.py
primary_output: 05-knowledge/results/jc2_bilateral_source_cone_reciprocal_orbits_thm4375.out
primary_script_sha256: e1467c3f5f9ebe6973c8a93981393ce6a3214456f3fda4530ae757d2ef1fcd5a
primary_output_sha256: 54c6250c18b50fb7df5a90123f85a0d10411ccbb9f43a73613c2f3d1c2b9f5a9
independent_referee_script: 04-computation/jc2_bilateral_source_cone_reciprocal_orbits_independent_referee_thm4375.py
independent_referee_output: 05-knowledge/results/jc2_bilateral_source_cone_reciprocal_orbits_independent_referee_thm4375.out
independent_referee_script_sha256: 256b0193a048a1db43a33aaf78724ef750e30cb64650e25c81f2f6648d2de0e0
independent_referee_output_sha256: 1938da0bcbbadf2a5ea2da29b224494c73c2abb4ea8fd6136edfe61c796ee76e
hash_basis: raw LF bytes
audit: >
  PASS WITH NONEMPTY-BLOCK, POSITIVE-SEGMENT, CUMULATIVE-ASYMPTOTIC, AND
  EMPTY-RAY-INTERVAL WORDING REPAIRS. The 6,375,878-assertion primary and
  3,338,648-check import-free referee independently verify the cone, orbit and
  address laws, generating functions, ray scales, fibre formulas, equality
  boundaries, asymmetry maximum, and all named low-order hostiles.
---

# THM-4375 -- Bilateral source-cone reciprocal orbits and fibre imbalance

**PROVED ELEMENTARY RELATIVE TO THM-4368 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THIS IS AN INTERNAL PACKET-GRAPH AND FIBRE THEOREM. IT HAS NO
`JC(2)`, `DC(2)`, LRC(14), TOURNAMENT, OR STERN--BROCOT DYNAMICAL
CONSEQUENCE.**

## 1. Setup and result

Fix `ell>=2` and put

```text
s=ceil(ell/2),                 rho=ceil(ell/3),
R(u,v)=(v,u).
```

Following THM-4368, the positive boundary pair `(u,v)` denotes the normalized
packet type

```text
N=v,                           n0=s+u-1.
```

Let `S_ell` be the set of source-realizable boundary pairs, and let

```text
B_ell=S_ell intersect R(S_ell)
```

be the **bilateral source cone**: both a packet type and its reciprocal
boundary reflection are source-realizable.

### Theorem A (the bilateral cone is exactly a path-power band)

```text
B_ell={(u,v) in Z_{>0}^2 : u,v>=rho and |u-v|<=s-1}.       (1)
```

Thus, if positive integers `rho,rho+1,...` are vertices and `(u,v)` is the
arc `u -> v`, reciprocal reflection reverses an arc in the looped
`(s-1)`-st power of the infinite path.  The quotient by `R` is the set of
its loops and unoriented edges.  This is not a tournament: pairs at distance
at least `s` are absent, and the loops `u=v` are honest fixed ties.

Every reciprocal orbit has a unique representative

```text
(u,v)=(rho+j+d,rho+j),          j>=0, 0<=d<=s-1.           (2)
```

The orbit is fixed exactly when `d=0`.

### Theorem B (triangular blocks and exact generating functions)

Let

```text
tau=u+v-1,
Addr(u,v)=T(tau-1)+u.
```

At fixed `tau`, the bilateral addresses form either the empty set or one
centered consecutive interval

```text
T(tau-1)+L_tau, ..., T(tau-1)+U_tau,                       (3)

L_tau=max(rho, ceil((tau+1-(s-1))/2)),
U_tau=min(tau+1-rho, floor((tau+1+(s-1))/2)).              (4)
```

When the interval in (3) is nonempty, its endpoints are reciprocal and obey

```text
L_tau+U_tau=tau+1,
[T(tau-1)+L_tau]+[T(tau-1)+U_tau]=tau^2+1.                (5)
```

Put `tau0=2rho-1`.  If `O_ell(tau)`, `M_ell(tau)`, and
`F_ell(tau)` count reciprocal orbits, oriented packet types, and fixed types
in block `tau`, respectively, then

```text
sum_(tau>=1) O_ell(tau) z^tau
    =z^(2rho-1)(1-z^s)/((1-z)(1-z^2)),                    (6)

sum_(tau>=1) F_ell(tau) z^tau
    =z^(2rho-1)/(1-z^2),                                  (7)

sum_(tau>=1) M_ell(tau) z^tau
    =z^(2rho-1)(1+2 sum_(d=1)^(s-1) z^d)/(1-z^2)
    =z^(2rho-1)(1+z-2z^s)/((1-z)(1-z^2)).                 (8)
```

Equivalently, for `k=tau-tau0>=0`,

```text
O_ell(tau)=#{d : 0<=d<=min(s-1,k), d congruent k (mod 2)}, (9)
M_ell(tau)=2 O_ell(tau)-1_[k even],                        (10)
F_ell(tau)=1_[k even].                                    (11)
```

There is therefore exactly one source-realizable reciprocal fixed point in
every odd block

```text
tau=2rho-1, 2rho+1, 2rho+3, ...,
```

and none in an even block.  The first is

```text
(u,v)=(rho,rho),
Addr=2rho^2-2rho+1.                                       (12)
```

More explicitly, the orbit `(j,d)` in `(2)` occupies the address pair

```text
A_-=T(tau-1)+rho+j,
A_+=A_-+d,                    tau=2rho-1+2j+d.            (12a)
```

Thus the address separation is exactly the boundary-coordinate separation
`d`, while the midpoint is forced by `A_-+A_+=tau^2+1`.  Equivalently,

```text
{A_-,A_+}={(tau^2+1-d)/2,(tau^2+1+d)/2}.                 (12b)
```

At `d=0` these two endpoints genuinely coalesce; there is no missing
orientation bit.  The fixed addresses are exactly

```text
Addr(w,w)=w^2+(w-1)^2,                  w>=rho,           (12c)
```

the integer centers of their odd triangular blocks.  This identity is
derived here; no external sequence property is being imported.

For a finite initial segment of `n>=1` vertices, put
`D=min(s-1,n-1)`.  Its number of reciprocal orbits, including loops, is

```text
(D+1)n-D(D+1)/2.                                         (13)
```

Choose one canonical address for each reciprocal orbit.  Both the oriented
and orbit address sets are infinite but have natural density zero.  More
precisely, as `K -> infinity`, their cumulative counts through blocks
`tau<=K` are

```text
sum_(tau<=K) M_ell(tau)=(s-1/2)K+O_ell(1),
sum_(tau<=K) O_ell(tau)=(s/2)K+O_ell(1),                 (14)
```

while `T(K)` natural-number addresses have appeared.

The low-order boundaries are included without exceptions.  At `ell=2`,
`(s,rho)=(1,1)`, so `B_2={(w,w):w>=1}`: all three generating functions are
`z/(1-z^2)`, every fibre has size one, and there are no nonfixed arcs.  At
`ell=3`, `(s,rho)=(2,1)`: there is one reciprocal orbit in every block
`tau>=1`, the oriented counts alternate `1,2,1,2,...`, the first nonfixed
orbit `(1,2)<->(2,1)` has fibres `1<->1`, and the first sharp imbalance is
`(2,3)<->(3,2)` with fibres `1<->3` in block `tau=4`.

### Theorem C (exact Stern--Brocot scale interval)

Let `gcd(p,q)=1`, `p,q>0`, and consider the primitive rational ray
`(u,v)=g(p,q)`.  If `p=q`, necessarily `(p,q)=(1,1)`, and the bilateral
scales are exactly

```text
g>=rho.                                                    (15)
```

If `p!=q`, the bilateral scales are exactly the integer interval

```text
ceil(rho/min(p,q)) <= g <= floor((s-1)/|p-q|).            (16)
```

The interval is empty when its lower endpoint exceeds its upper endpoint.
Its size is the positive part of upper minus lower plus one.  Thus every
non-diagonal primitive ray has only finitely many source-realizable reciprocal
scales, whereas the fixed ray `1/1` has infinitely many.  Reducing `(u,v)` by
its gcd forgets the scale `g`, and that scale is not recoverable from the
primitive Stern--Brocot address.  In particular, every fixed point lies over
the single primitive ray `1/1`; this quotient collapses the infinite fixed
sequence `(rho+j,rho+j)` to one rational node.

### Theorem D (source fibres before and after reflection)

For `(u,v) in S_ell`, write `mu_ell(u,v)` for THM-4368's number of source
monomials expanding to that packet type.  For a reciprocal orbit in `B_ell`,
use its canonical representative

```text
(u,v)=(w+d,w),             w>=rho, 0<=d<=s-1,
```

the two endpoint fibre sizes are exactly

```text
mu_+(w,d)=mu_ell(w+d,w)
 =min(w,s-1+d)-max(0,ell-2w)+1,                           (17)

mu_-(w,d)=mu_ell(w,w+d)
 =min(w+d,s-1-d)-max(0,ell-2w-2d)+1.                     (18)
```

Equivalently, with `[x]_+=max(x,0)`,

```text
mu_+(w,d)=s+d-[s-1+d-w]_+-[ell-2w]_+,                    (19)
mu_-(w,d)=s-d-[s-1-w-2d]_+-[ell-2w-2d]_+.                (20)
```

Both are positive, and they satisfy the sharp conservation ceiling

```text
mu_+(w,d)+mu_-(w,d)<=2s.                                 (21)
```

Equality holds if and only if

```text
w>=max(s,s-1+d),                                          (22)
```

and then the endpoint fibres stabilize separately to

```text
(mu_+,mu_-)=(s+d,s-d).                                    (23)
```

Consequently

```text
|mu_+-mu_-|<=2s-2.                                       (24)
```

For `s>=2`, equality in `(24)` holds if and only if, after choosing the
representative with first coordinate at least the second,

```text
d=s-1,                  w>=2s-2.                          (25)
```

Then `(mu_+,mu_-)=(2s-1,1)`.  The first sharp-asymmetry orbit occurs in
block `tau=5s-6`.  For `s=1`, every bilateral orbit is fixed and the
asymmetry is zero.

The exact fibres themselves are the exponent tuples `(a,b,c,e)` obtained by

```text
max(0,ell-2N)<=e<=min(N,n0-N),
c=N-e,        a=2N+e-ell,        b=n0-N-e,                (26)
```

with `(N,n0)=(v,s+u-1)`.  Equations `(17)` and `(18)` count these two
intervals; reciprocal boundary reflection does not in general biject them.

## 2. Proof

THM-4368 gives source realizability of `(u,v)` precisely when

```text
v>=rho,                  u-v>=1-s,
u+v>=ell-s+1.                                             (27)
```

Applying `(27)` also to `(v,u)` gives

```text
u,v>=rho,                |u-v|<=s-1,
u+v>=ell-s+1.                                             (28)
```

For `ell>=6`,

```text
2ceil(ell/3)>=2ell/3>=ell/2+1>=floor(ell/2)+1=ell-s+1;
```

the same inequality is immediate for `ell=2,3,4,5`.  Hence the last wall in
`(28)` follows from `u,v>=rho`, proving `(1)`.  Taking
`w=min(u,v)` and `d=|u-v|` proves the unique parametrization `(2)`.

At fixed `tau`, write `S=tau+1=u+v`.  Equation `(1)` is equivalent to

```text
rho<=u<=S-rho,            ceil((S-(s-1))/2)<=u
                          <=floor((S+(s-1))/2),           (29)
```

which proves the interval `(3)`--`(4)`.  Both intervals in `(29)` are
symmetric under `u -> S-u`, so `L_tau+U_tau=S`.  The triangular-address
identity follows from

```text
2T(tau-1)+(tau+1)=tau^2+1.
```

For the canonical representative `(2)`,

```text
tau=2rho-1+2j+d.                                         (30)
```

Summing `z^tau` over `j>=0` and `0<=d<=s-1` proves `(6)`.
Keeping only `d=0` proves `(7)`, and giving weight two to `d>0` proves
`(8)`.  Coefficient extraction gives `(9)`--`(11)`.  Formula `(13)` follows
by summing the `n-d` pairs at each distance `0<=d<=D`.  The eventual
two-block coefficients of `(8)` sum to `2s-1`, and those of `(6)` sum to
`s`; this proves `(14)` and hence zero natural density.

On a primitive ray, `(1)` becomes

```text
g min(p,q)>=rho,             g|p-q|<=s-1.                 (31)
```

Solving these two integer inequalities gives `(15)`--`(16)`.

Finally, substitute

```text
(N,n0)=(w,s+w+d-1)          and
(N',n0')=(w+d,s+w-1)
```

into THM-4368's exact formula

```text
mu(N,n0)=min(N,n0-N)-max(0,ell-2N)+1.
```

This is `(17)`--`(18)`.  Replacing each minimum by
`min(x,A)=A-[A-x]_+` proves `(19)`--`(20)`.  Dropping their four nonnegative
defects proves `(21)`.  All four defects vanish exactly when
`w>=max(s,s-1+d)`, proving `(22)`--`(23)`.  Positivity plus `(21)` gives
`(24)`.  To attain `2s-2`, the two positive fibre sizes must be `2s-1` and
`1`; since `mu_-<=s-d<=s`, the larger one must be `mu_+`, and
`mu_+<=s+d` forces `d=s-1`.  Its first defect then forces `w>=2s-2`.
Conversely these conditions give `(2s-1,1)`.  Since `rho<=s<=2s-2` for
`s>=2`, the least such `w` belongs to the bilateral cone; substituting
`w=2s-2,d=s-1` into `tau=2w+d-1` gives the claimed first block `5s-6`.
This proves `(25)` and completes the proof.

## 3. Earliest `ell=10` controls and hostile

Here `(s,rho)=(5,4)`.  The first bilateral block is the honest tie

```text
(u,v)=(4,4),       tau=7,       Addr=25,       mu=3.      (32)
```

Its three source monomials, written as exponent tuples `(a,b,c,e)` for
`x^a u^b p^c y^e`, are

```text
(0,2,2,2), (1,1,1,3), (2,0,0,4).
```

The very next block is already the minimal reciprocal-fibre hostile:

```text
(u,v)=(4,5) <-> (5,4),
tau=8,       Addr=32 <-> 33,       mu=4 <-> 3.            (33)
```

The `Addr=32` fibre is

```text
(0,3,5,0), (1,2,4,1), (2,1,3,2), (3,0,2,3),
```

whereas the `Addr=33` fibre is

```text
(0,3,2,2), (1,2,1,3), (2,1,0,4).
```

Thus the ambient trace reciprocity from THM-4368 survives, but it has no
cardinality-preserving source-monomial lift even at the first nontrivial
bilateral orbit.  The fixed point in `(32)` must remain a tie, and `(1)`
shows why forcing the band into a tournament would be doubly lossy: it would
both orient loops cosmetically and invent all missing long edges.

There is also an exact finite scale ladder on the primitive ray `4/5`.  By
`(16)` its bilateral scales are precisely `g=1,2,3,4`; the fifth scale exits
the width-four band.  The source-fibre counts in the displayed orientation
and its reciprocal are

```text
g       (u,v)       Addr(u,v) <-> Addr(v,u)       mu(u,v) <-> mu(v,u)
1       (4,5)          32     <-> 33                  4    <-> 3
2       (8,10)        144     <-> 146                 3    <-> 7
3       (12,15)       337     <-> 340                 2    <-> 8
4       (16,20)       611     <-> 615                 1    <-> 9
```

This extends THM-4368's same-ray hostile: primitive reduction preserves the
rational direction but loses the exact scale interval, triangular address,
trace, and source-fibre size.  It is an internal packet calculation, not a
continued-fraction dynamical theorem.

## 4. Connection ledger and scope

```text
source:      source-realizable THM-4368 boundary packet types
target:      looped path-power arcs / reciprocal orbits / natural addresses
map:         (u,v) -> u->v; R(u,v)=(v,u); Addr=T(u+v-2)+u
preserved:   exact ordered boundary pair; reciprocal involution; trace type
destroyed:   individual source exponent e after passing to (u,v)
ray quotient destroyed: gcd scale g; trace/fibre data can vary with g
needed sidecar: g for a primitive rational ray; full exponent tuple for a
                source monomial rather than only its packet trace
hostile:     ell=10, Addr 32<->33 has fibre sizes 4<->3
```

The path-power language is an exact graph model internal to this packet
family.  The Stern--Brocot statement is only the exact classification of gcd
scales on rational rays.  The tournament observation is a validity warning:
the intrinsic relation is a partial looped graph with ties, not a tournament.
None of these statements transfers a theorem between Jacobian, lonely-runner,
tournament, or continued-fraction problems.

## 5. Reproduction

The primary derives the bilateral cone independently from direct source
exponent enumeration and checks every displayed consequence on its printed
finite universe.  The clean-room referee imports no repository code or Python
modules; it uses a separate forward/inverse exponent reconstruction, direct
block decoding, and rational generating-function recurrences.

```text
python -B 04-computation/jc2_bilateral_source_cone_reciprocal_orbits_thm4375.py
python -B 04-computation/jc2_bilateral_source_cone_reciprocal_orbits_independent_referee_thm4375.py
```

The primary passes `6,375,878` assertions and the independent referee passes
`3,338,648` checks.  Each agrees byte-for-byte with its frozen raw-LF output
under normal, optimized, isolated, and fixed-hash-seed runs.

```text
primary script:    e1467c3f5f9ebe6973c8a93981393ce6a3214456f3fda4530ae757d2ef1fcd5a
primary output:    54c6250c18b50fb7df5a90123f85a0d10411ccbb9f43a73613c2f3d1c2b9f5a9
referee script:    256b0193a048a1db43a33aaf78724ef750e30cb64650e25c81f2f6648d2de0e0
referee output:    1938da0bcbbadf2a5ea2da29b224494c73c2abb4ea8fd6136edfe61c796ee76e
hash basis:        raw LF bytes
```

The finite universes include `ell=2..60, u,v=1..80`, blocks through
`tau=240`, primitive coprime rays through `p,q=30` with scales through 90,
and fibre tests through `ell=100`, beyond every asserted equality threshold.
They verify consequences of the proved formulas, not an open-conjecture claim.
