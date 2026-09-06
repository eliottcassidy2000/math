# Virtual Hadamard doubling, midpoint defects, and the remaining signed transport

Status: **PROVED / INDEPENDENTLY AUDITED** for the new path
inclusion and fixed-family transport certificate (root and observer).
The virtual sign is a
direct consumer of proved THM-4440 and the incoming complete-path factors.
The universal rootwise transport inequality is **OPEN**; the declared
finite banks below are **FINITE-EXACT**, not an unbounded proof.

## 1. Inheritance and the exact objects

The source is the incoming
[complete binomial-path / full-support Hadamard factorization](overnight3_20260906_moments_root_geometry.md),
read in full. The new consumer is
[THM-4440, signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md),
whose [full proof](nc2_signed_duplication_overnight_hexagon_sep05.md)
proves strict negativity of a square coefficient at an interior zero
coefficient of a real-rooted carrier. The complete first-channel roots are
supplied by
[THM-4436, complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).

The canonical hostile is the actual-mass word `PQQPQPQ` in the
[contiguous-operation note](nc2_channel_contiguous_overnight_hexagon_sep05.md):
adjacent masses need not interlace. Its mixed-sign individual Fourier
pairs prohibit replacing a signed sum by a sum of individually negative
terms. The corrected near miss in the incoming proof is deleting the
negative-index auxiliary support; the concrete prefix `21+20t+5t^2`
fails real-rootedness. We retain that full support here. The least-used
sidecar is the **specific midpoint lattice section hit or skipped by a
complete path**, not just the endpoint count.

The live board is: full auxiliary paths; ordinary carrier roots;
Hadamard versus ordinary squares; midpoint residues and skipped levels;
negative phase evaluation; complete affine carries; quotient trace/norm.
The main new connection has an exact preserved coefficient map and an
equally exact failure boundary.

Let integers satisfy

```text
0<A<B, delta=B-A, h>=1, x>=0,
0<=r<B, 0<=z<A, y=B*h+r, m=x+y+z, L=floor(x/delta).
```

No coprimality or positive-charge interpretation is required for the
next two lemmas. Define finite two-sided sequences, zero outside the
specified nonnegative-count domains,

```text
alpha_j=binom(m,z+A*j),                 j>=0, z+A*j<=m,
beta_j=binom(x+y-A*j,x+delta*j),        x+delta*j>=0, y-B*j>=0.
```

Thus the full beta support is `-L<=j<=h`. Put

```text
F(t)=t^L sum_(j>=0) alpha_j*t^j,
G(t)=sum_(j=-L)^h beta_j*t^(j+L),
P(t)=sum_(j=0)^h alpha_j*beta_j*t^j
    =sum_(j=0)^h m!*t^j/[(x+delta*j)!(y-B*j)!(z+A*j)!].       (1)
```

For ordinary coefficientwise multiplication `star`, the full incoming
identity is `F star G=t^L P`. Both `F` and `G` are real-rooted;
`F` has exactly its declared zero factor `t^L`, and `G(0)>0`.

For precision, this real-root input follows from the actual path
construction: `binom(U+V+(p-q)j,V+p*j)` is the path matrix entry from
`(-q*i,p*i)` to `(U-q*j,V+p*j)`. Every finite minor is a signed sum of
path families. Swapping tails at the first intersection cancels the
intersecting families, and source/sink order excludes a nonidentity
disjoint pairing. Thus its Toeplitz matrix is totally nonnegative.
Translating the **entire** finite support preserves this construction.
The only cited characterization is that a finite PF sequence has a
generating polynomial with real nonpositive roots; see
[Branden, arXiv:math/0403364v2, Section2, Theorem2.1 and its finite consequence](https://arxiv.org/pdf/math/0403364#page=2),
whose statement was checked directly. Apply it with `(p,q)=(A,A)`
for alpha and `(delta,B)` for beta. Their positive translated constant
terms exclude zero; multiplication by `t^L` accounts for F's zeros.
No simplicity theorem or infinite Hadamard closure is needed here.

Define the virtual row

```text
V(t)=t^(-2L) [F(t)^2 star G(t)^2]
    =sum_(j=0)^(2h) V_j*t^j,
V_j=(sum_i alpha_i*alpha_(j-i))*(sum_l beta_l*beta_(j-l)).      (2)
```

The actual doubled Laurent row is

```text
Q_raw(t)=sum_j q_j*t^j,
q_j=(2m)!/[(2x+delta*j)!(2y-B*j)!(2z+A*j)!],                 (3)
```

with `q_j=0` if any factorial argument is negative. Only
`-1<=j<=2h+1` can occur. The possible lower carry is also removed if
`2x-delta<0`; it must not be inserted outside its genuine count domain.
For positive-charge first fibres both usual carry endpoints are retained.

## 2. The virtual row is strictly negative at every first root

Here is the useful general lemma. Let real-rooted real polynomials
`F,G` have degrees `p,q` and orders at zero `r_F,r_G`. Suppose

```text
r_F<q,       r_G<p.
```

For any nonzero real `t`, set

```text
H_t(u)=u^q F(u)G(t/u).
```

This is a real-rooted polynomial: its nonzero roots are those of F
and `t` divided by the nonzero roots of G. Its order at zero is
`r_F` and its degree is `p+q-r_G`. The indicated coefficient is
strictly interior, and literal multiplication gives

```text
[u^q]H_t=(F star G)(t),
[u^(2q)]H_t^2=(F^2 star G^2)(t).                            (4)
```

THM-4440 therefore proves

```text
(F star G)(t)=0  ==>  (F^2 star G^2)(t)<0                    (5)
```

for every such real `t!=0`. Repeated and zero roots are allowed with
the stated order inequalities. If the entire Hadamard overlap is
absent, the interior condition fails and a strict conclusion is false.

For the factors in (1), `r_F=L`, `q=L+h>L`, `r_G=0`, and `p>=L+h>=1`.
All roots of P are simple and negative by THM-4436. At any such root
`rho`, equation (5) and the positive even factor `rho^(2L)` give

```text
V(rho)<0.                                                   (6)
```

This is an honest rational-polynomial carrier reduction for the
**virtual** row. It does not identify V with the true doubled moment.

## 3. Actual coefficients dominate virtual coefficients by midpoint inclusion

For every integer `j`, under the weak parameters of Section1,

```text
q_j>=V_j>=0,                                                (7)
```

where `V_j=0` outside `0<=j<=2h`. In particular the entire Laurent
defect `D_raw=Q_raw-V` has nonnegative coefficients.

To prove this, factor the actual coefficient as

```text
q_j=binom(2m,2z+A*j)
       *binom(2x+2y-A*j,2x+delta*j),                         (8)
```

with the nonnegative-count convention. The first factor counts all
northeast unit paths to `(2(m-z)-A*j,2z+A*j)`. The alpha convolution
in (2) counts those paths passing through a selected midpoint

```text
(m-z-A*i, z+A*i).
```

Every such point has coordinate sum m. A northeast path reaches that
level once, so splitting and concatenating paths is injective. The
selected vertices retain the residue `Y=z mod A`; when `A>1`, other
midpoint residues can be missing. Hence the alpha convolution is at
most the first binomial in (8).

For beta, all selected split vertices are

```text
(y-B*l, x+delta*l).
```

They lie on the linear level

```text
delta*X+B*Y=delta*y+B*x.
```

Each northeast edge increases this functional by the positive amount
delta or B. A path can hit at most one selected vertex. Splitting
there and translating its second half gives exactly
`beta_l*beta_(j-l)`, with **all** admissible negative l retained.
Thus the beta convolution injects into the paths to
`(2y-B*j,2x+delta*j)`, counted by the second binomial in (8).
Some paths jump across the level without a vertex on it; with
nonprimitive steps the selected lattice class can also be a proper
subset of its integer vertices. Neither phenomenon invalidates the
injection. Multiplying the two inequalities proves (7).

More precisely, `D_raw` counts pairs of full paths for which at least
one component misses its selected midpoint section. It is the
disjoint union "alpha misses" and "alpha hits but beta misses."
This gives a faithful combinatorial object for the extra coefficients.
No coprimality is needed for this description.

Crucially, (7) does **not** imply an inequality between values at a
negative argument. A nonnegative-coefficient polynomial such as 1 is
positive at the negative root of `1+t`. The needed signed transport is
an additional assertion, not a consequence of path counting alone.

## 4. The first failed equality and the surviving signed target

Already for `(A,B,h,r,z,x)=(1,3,1,0,0,1)`, the actual primitive
support is `(-3,1,9)` and the first mass is four. The full factors are

```text
F=(1+t)^4,     G=4+t,     P=4+4t,
V=16+64t+28t^2,
Q_raw=28+280t+28t^2.
```

At `rho=-1`, `V=-20`, whereas `Q_raw=-224`. The two rows are neither
equal nor proportional, despite both being negative. The first failed
implication is identifying parameter doubling of the path endpoints
with squaring both auxiliary generating polynomials.

The natural repaired target is

```text
P(rho)=0  ==>  D_raw(rho)=Q_raw(rho)-V(rho)<0.                (9)
```

Together with (6), this would prove the desired actual doubled sign.
The universal version of (9) on positive-charge collided first fibres
remains **OPEN**. The midpoint-count description preserves all extra
weights, but their negative-phase signed sum has not been evaluated
in general. In particular, no termwise negative frequency-product
claim is being revived.

A three-channel exact control is `(A,B,h,r,z,x)=(1,2,2,0,0,1)`,
corresponding to `(-4,1,6)`:

```text
P=5+30t+10t^2,
V=37+620t+2070t^2+1440t^3+210t^4,
Q_raw=45+840t+3150t^2+2520t^3+210t^4,
D_raw=8+220t+1080t^2+1080t^3.
```

Modulo P, the defect is `1088+6160t`. At the two roots, its trace
is `-16304` and its norm is `50304`, so both values are negative.
This is an exact small control for the stronger target, not its proof.

## 5. An unbounded transport control on the incoming lower-carry family

The full incoming
[width-fifteen family proof](trinomial_width15_empty_core_returns_sep06.md)
was read, including its certificate audit and next boundary. It already
proves actual negativity for `(-15,2g-15,3g-15)`, every integer
`g>=8` with `gcd(g,15)=1`. Its normalized first polynomial is

```text
P_g(t)=6t^2+20(g-5)t+(g-5)(g-6),
(A,B,h,r,z,x)=(2,3,2,0,1,g-7).
```

We do not count that return theorem again. The stronger new control is
that **the actual row is below the virtual row at both first roots**,
uniformly over this same unbounded family. This includes the lower
carry and lies mostly outside the cubic real-rooted-core phase sector.

Define fixed-degree binomial polynomials

```text
alpha_i=binom(g,1+2i),
beta_i=binom(g-1-2i,6-3i),
V_g(t)=sum_(j=0)^4 [sum_(i=0)^j alpha_i alpha_(j-i)]
                     [sum_(l=j-2)^2 beta_l beta_(j-l)] t^j,
Q_g(t)=sum_(j=0)^5 (2g)_(15-j)*t^j/[(15-3j)!(2j)!].         (10)
```

At integral `g>=8`, these are the exact rows above, including the
vanishing boundary binomials when a selected split is unavailable.
The true anchored second row is `t^(-1)Q_g(t)`. Put

```text
D_g(t)=t^(-1)Q_g(t)-V_g(t),
D_g mod P_g=C(g)+E(g)*t.
```

The inverse is exact in the first-root quotient:
`t^(-1)=-(6t+20(g-5))/[(g-5)(g-6)]`. Write
`K(g)=g(g-1)(g-2)(g-3)(g-4)`. Exact polynomial identities are

```text
trace(D_g)=-K(g)(g-5)*J(g)/294226732800,
norm(D_g)=K(g)^2(g-5)*H(g)/961881892157362176000000,           (11)
```

where

```text
J(g)=52280181899g^8-1435119555457g^7+16676493669661g^6
 -106396392380355g^5+403767244795376g^4-920771860729028g^3
 +1207291516251904g^2-803077287018720g+190335031286400,

H(g)=1738714564090842281g^17-108547937091182402339g^16
 +3143807684397515907261g^15-56040143508356807289959g^14
 +687752861501332059667787g^13-6158087029819152968307873g^12
 +41588689690171421410824927g^11-215967719730019617541515213g^10
 +871039120856486066254008288g^9-2736437794988899479308430712g^8
 +6673222243886469185321657488g^7-12505944664059081956718237072g^6
 +17689194160735783031309784064g^5-18348193363506055389859038976g^4
 +13337602056379737902528845824g^3-6310852212346447791291872256g^2
 +1703694415249151173406638080g-195068578693303718052249600.
```

Every coefficient of `J(s+8)` and `H(s+8)` is strictly positive; the
complete coefficient lists are printed in the companion exact output.
All prefactors in (11) are positive for real `g>=8`. The first roots
are real and distinct, with discriminant `8(g-5)(47g-232)>0`.
A negative trace and positive norm force both defect values negative.
Thus (9) holds on this entire fixed family, not just the sampled g.

Here is the algebraic validity gate behind (11). The sole possible
denominator in the inverse-carry term divides `(g-5)(g-6)`, and the
constant term of Q contains both factors in `(2g)_15`. It cancels
exactly. Reducing positive powers by
`t^2=-(10/3)(g-5)t-(g-5)(g-6)/6` shows
`deg C<=14`, `deg E<=13`, hence trace degree at most14 and norm at
most28. The checker obtains (11) over the exact field `Q(g)` and
independently reproduces the quotient with standard-library Fraction
recurrences at all 29 distinct points `g=8,...,36`. These degree bounds
make the latter an exact polynomial identity certificate, not merely
a specialization census. Nonprimitive g values are used only to
certify these identities, not as first-return claims.

This stronger transport certificate is more complicated than the
incoming actual-sign certificate: normalized degrees `8,17` here
versus `3,7` there. The virtual row supplies a new structural object;
it has not automatically simplified the symbolic elimination.

## 6. Finite controls and the next precise obligation

The [standalone checker](../../04-computation/nc2_hadamard_transport_overnight_hexagon_sep05.py)
uses exact SymPy polynomial arithmetic and rational root isolation,
literal midpoint path dynamic programming, and independent Fraction
quotient recurrences. It imports no repository mathematical producer.

Its declared root-sign head contains all primitive positive-charge first
fibres with `A<B<=6`, `h in {1,2,3}`, every `r,z`, and `x in {1,2,3}`:
638 rows. Its wide bank contains arithmetic-progression heights
`{2,3,4,5,6,8,10,12}`, both r residues, x values
`{5,7,11,101,503,997}`, followed by 300 seeded raw proposals with
`B=3..12`, `h=1..8`, all residue ranges, and x drawn from
`{1,2,3,5,7,11,23,101,997}`. The same exact eligibility filter retains
229 indexed wide controls. Every root of every first row is covered.
There are 1,271 head roots and 1,179 wide roots. The banks contain
865 distinct parameter tuples; indexed bank totals are not added as
though they were disjoint.

The eligibility test is `gcd(A,B)=1`, `A*x-(B-A)z>0`, and
`gcd(A*(B*h+r)+B*z,m)=1`. It makes the support
`(-a,m*A-a,m*B-a)` primitive with first return at m. For every target
root, exact interval evaluation of the polynomial remainder, with
rational refinement, certifies both `V(rho)<0` and
`Q_raw(rho)-V(rho)<0`. Multiplication by
`rho^(-floor(2z/A))` restores the lower-carry sign before it is tested.
Neither midpoint evaluation nor a floating root estimate is used.

The independent path-DP universe includes noncoprime A,B and x=0,
which are legal for (7), even when they lack a positive-charge first
interpretation. It separately counts all paths and those hitting the
selected midpoint. The fixed-family algebraic certificate retains
every shifted coefficient and all 29 interpolation controls.

The exact next obligation is to evaluate the signed generating row of
**missed midpoint path pairs** at a root of the one-block Hadamard row.
Candidate operations include pairing complementary midpoint residues
or skipped-level crossing classes and using the proved contiguous
Euler identities on whole grouped responses. The earlier positive
individual-pair hostile prevents a termwise sign shortcut. An
equivalence between endpoint doubling and auxiliary squaring has
already been refuted; any proposed transport must preserve (2)--(3).

No universal trinomial two-rung closure, no general sharp Laurent
width theorem, and no general actual-mass interlacing is claimed here.
Both root and observer independently audited the full written argument,
including the fixed-family quotient, degree bounds, signs and boundaries.
The observer separately replayed the path, symbolic-family and hostile
controls: 4,361 explicit gates, 144 indexed path cases, 864 endpoints,
and all 29 independent Fraction certificate points. The full producer
passes 20,908 explicit gates; none is a Python assertion.

Reproduce from the repository root:

```bash
python3 -B 04-computation/nc2_hadamard_transport_overnight_hexagon_sep05.py
python3 -B -O 04-computation/nc2_hadamard_transport_overnight_hexagon_sep05.py
```

Normal, optimized and [stored output](nc2_hadamard_transport_overnight_hexagon_sep05.out)
agree byte for byte. Raw LF SHA256 values are:

```text
source 6f37cd9b1702f616f234e9cf032b0b4ae286768ab6bc622ff9820f69d1813ce2
output 4c8f58cfa2aa9c63af2c8a8341ceb5e7a0fb654e79e60695a722a8d27b937ca6
trace  bbe0f30470377ea72ee78a8b80f2ae152dece0f1b2ba32b195f0cfb18f33ab88
```
