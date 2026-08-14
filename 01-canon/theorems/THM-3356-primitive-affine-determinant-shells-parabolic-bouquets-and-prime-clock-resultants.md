---
id: THM-3356
title: "Primitive affine determinant shells, parabolic bouquets, and prime-clock resultants"
status: RESERVED / UNPROVED PROVISIONAL CANDIDATE UNDER AUDIT
source: codex-2026-08-14-affine-determinant-shells
depends_on: []
related: []
---

# THM-3356 -- primitive affine determinant shells and prime-clock resultants

**RESERVED / UNPROVED PROVISIONAL CANDIDATE UNDER AUDIT.**

This proof candidate is quarantined from the proved graph while its exact
companion, scopes, and independent hostiles are audited.  In particular it
does not prove the carrier inequalities or physical census reserved for
THM-3355, and it does not decrement the open LRC(14) affine-tail ledger.

## 1. Inheritance, portfolio, and the two views

The closest proved mechanism is
[THM-3333](THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md):
for the spinor square

```text
Phi(m,n)=(m^2-n^2,2mn,m^2+n^2)                          (1)
```

and Lorentz pairing one has

```text
<Phi(u),Phi(v)>_L=2 det(u,v)^2.                         (2)
```

[THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
identifies the consecutive-parameter U-spine and its parabolic horocycle;
[THM-3336](THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md)
proves primitive Gaussian content curvature; and
[THM-3346](THM-3346-u-spine-prime-toggle-root-atlas-and-conjugation-monodromy.md)
splits common U-spine prime powers into difference and reflected-sum channels.
[THM-3353](THM-3353-split-prime-parabolic-branch-transplant-and-unary-transducer-compiler.md)
compiles a selected split-prime toggle on special U-spine subsequences.

The canonical hostile is [MISTAKE-371](../MISTAKES.md): a selected primary
channel is not automatically the full Gaussian content.  The least-used
sidecar here is the residue of an affine-shell point inside its parabolic
stabilizer orbit.

The session portfolio became the following three-item board.

| lane | object / representation | question | lost coordinate / cheapest test |
|---|---|---|---|
| anchor | the open disconnected-low rays `dq-(d+a)p=c` | can the incoming primitive quotient expose a reusable invariant? | physical context and mass; test the literal `22,890 -> 14,168` quotient |
| niche | `S_c(u)={x:det(x,u)=c}` | how many integral parabolic orbits and contents occur? | the Lorentz scalar loses sign and orbit residue; test `c=2` |
| wildcard | the two speeds as U-spine indices `p,q` | which shared Gaussian prime clocks can occur along one ray? | a gcd loses its two local branches; test fixed polynomial resultants |

Both non-anchor views survive, but they preserve different predicates.

| Item | Determinant-shell view | Dual U-spine-index view |
|---|---|---|
| source | one primitive affine ray | its two integer coordinates `(p,q)` |
| map | `(q,p) -> Phi(q,p)` | `(p,q) -> (z_p,z_q)`, `z_t=(t+1)+it` |
| preserved | transverse determinant, torus intersection, Lorentz shell | common hypotenuse factors and same/opposite Hensel branch |
| destroyed | determinant sign, parabolic-orbit residue, spinor content | LRC context, carrier chamber, owner, phase, and physical overlap |
| sidecar | oriented shell charge and residue modulo `|c|` | row label `(d,a,c)` and the two channel labels |
| decisive hostile | `c=2` is two horocycles, not one | the norm-`85` toggle changes the LRC carrier coefficient |

## 2. A primitive affine shell is a finite parabolic bouquet

Let

```text
u=(D,d) in Z^2,             gcd(D,d)=1,
S_c(u)={x in Z^2:det(x,u)=c},          c!=0.             (3)
```

Choose `v` with `det(v,u)=1`.  Every shell point has a unique expression

```text
x=c v+k u,                 k in Z.                      (4)
```

Indeed `(v,u)` is a unimodular basis, and taking determinant with `u`
forces its first coordinate to be `c`.  Replacing `v` by `v+j u` changes
`k` by `-jc`; hence the class

```text
[k]_c in Z/cZ                                                (5)
```

is intrinsic.

Define the integral parabolic transvection

```text
T_u(x)=x+det(x,u)u.                                      (6)
```

For `u=(D,d)` its matrix is

```text
T_u=[1+Dd  -D^2]
    [ d^2  1-Dd],        det(T_u)=1,       T_u u=u.      (7)
```

On (4),

```text
T_u(c v+k u)=c v+(k+c)u.                                (8)
```

Therefore `S_c(u)` is the disjoint union of exactly `|c|` integral
parabolic orbits, indexed by (5).  The affine enumeration

```text
x_n=x_0+n u                                               (9)
```

interleaves those `|c|` orbits; it is one parabolic orbit exactly when
`|c|=1`.  This is the first boundary: an arbitrary affine ray must not be
called a single horocycle.

Unimodularity of `(v,u)` also gives the exact content formula

```text
gcd(x_1,x_2)=gcd(c,k).                                  (10)
```

Consequently the content-`g` stratum, for `g` dividing `|c|`, consists of

```text
phi(|c|/g)                                               (11)
```

parabolic orbits.  After dividing by `g`, each becomes a primitive orbit in
`S_(c/g)(u)`.  This is a canonical divisor stratification, not an empirical
moving-gcd filter.

There is one further parity sidecar.  If `x=g y` with `y` primitive, then

```text
gcd(Phi(x))=g^2 epsilon(y),
epsilon(y)=2  iff both coordinates of y are odd,
epsilon(y)=1  otherwise.                                (12)
```

Thus the shell residue and divisor `g` recover the spinor content only after
the retained mod-two colour is supplied.  This is the same raw-versus-
primitive parity boundary as THM-3333, now resolved orbitwise.

## 3. Lorentz shells and sum-of-two-squares quadratics

Equations (2)--(3) give, for every `x in S_c(u)`,

```text
<Phi(x),Phi(u)>_L=2c^2.                                 (13)
```

For `x=g y`, the primitive normalized shell is

```text
<Phi(y),Phi(u)>_L=2(c/g)^2.                             (14)
```

In the torus-curve language of THM-3333, `|c|` is the geometric
intersection number.  The case `|c|=1` is precisely one oriented Farey-star
orbit about `[u]`; larger `|c|` are higher intersection shells.  Passing to
projective spinors identifies the two oriented signs, while retaining only
the residues coprime to `c` as primitive determinant-`|c|` vertices.

Now fix `x_0 in S_c(u)` and put

```text
A=||u||^2,        h=x_0.u,        C=||x_0||^2,
Q(n)=||x_0+n u||^2=A n^2+2h n+C.                        (15)
```

Lagrange's two-square identity gives

```text
AC-h^2=c^2,
disc(Q)=(2h)^2-4AC=-4c^2,                              (16)
A Q(n)=(An+h)^2+c^2.                                   (17)
```

Thus every primitive affine determinant shell carries a quadratic
hypotenuse sequence of square discriminant `-4c^2`.  The U-spine polynomial
`2n^2+2n+1` is the specialization

```text
x_0=(1,0),        u=(1,1),        A=2, h=c=1.           (18)
```

For a positive point `x_n=(q_n,p_n)` with `q_n>p_n`, the corresponding
right triangle has inradius

```text
r_n=p_n(q_n-p_n).                                       (19)
```

This is THM-3333's triangular-number defect on the same shell; it is not an
LRC margin.

### 3A. The tame affine prime-toggle atlas

Let

```text
M=product_j p_j^(e_j),
p_j=1 mod 4 distinct odd primes,       gcd(M,Ac)=1.      (20)
```

For each prime power choose a square root `iota_j^2=-1`.  Equation (17)
shows that the roots of `Q(n)=0 mod M` are exactly

```text
n=A^(-1)(-h+epsilon_j c iota_j) mod p_j^(e_j),
epsilon_j in {+1,-1}.                                  (21)
```

The derivative `2(An+h)` is a unit at every root, so the two branches are
distinct Hensel lifts.  CRT makes the root set a free transitive
`F_2^omega(M)` torsor.  Toggling every prime sends

```text
n |--> -2h A^(-1)-n mod M.                              (22)
```

To see the Gaussian allocation literally, identify vectors with Gaussian
integers and put

```text
z_n=(x_0+n u)_1+i(x_0+n u)_2,
L_n=(An+h)+ic=u_C conjugate(z_n),       N(L_n)=A Q(n),   (23)
```

where `u_C=u_1+i u_2`.  Since `A` and `c` are units modulo `M`, a root is a
nonzero isotropic vector at every `p_j`.  Hence

```text
G_M(n)=gcd_(Z[i])(M,L_n)                                (24)
```

has norm `M`, up to a Gaussian unit, and (21) chooses one of the two factors
above every `p_j`.  The involution (22) conjugates all factor choices.  Its
quotient is the usual fixed-norm parent torsor `F_2^omega(M)/<1>`.

For roots `r,s`, define selected channels

```text
delta_-=gcd(M,s-r),
delta_+=gcd(M,A(r+s)+2h).                               (25)
```

At a prime power, equal signs in (21) put that full power in `delta_-`, and
opposite signs put it in `delta_+`.  Therefore

```text
gcd(delta_-,delta_+)=1,          delta_- delta_+=M.      (26)
```

Moreover, writing `cont` for coordinate gcd,

```text
delta_-=gcd(M,cont(z_s conjugate(z_r))),
delta_+=gcd(M,cont(z_r z_s)).                            (26a)
```

This follows by multiplying (23):
the extra factors `A` and `u_C^2` are units at `M`, while the imaginary
coordinates distinguish `s-r` from `A(r+s)+2h`.

Equations (20)--(26) extend the U-spine root cube of THM-3346 to every tame
primitive affine shell.  They do **not** identify selected channels with
full contents when extra common primes occur; MISTAKE-371 remains
load-bearing.  Primes dividing `c` are discriminant-ramified, primes dividing
`A` destroy the inverse in (21), and `2` is parity-ramified.  No Boolean
toggle is claimed at those primes.

## 4. Ramification is exactly the coherent Gaussian-toggle gate

There is a complementary law at a split prime dividing the anchor norm.  Let
`p=1 mod 4` be prime and let `H` be a primitive integral Gaussian-rotation
matrix satisfying

```text
H^T H=p^2 I,             det(H)=p^2.                    (27)
```

Let `u` be primitive with

```text
H u=0 mod p.                                                (28)
```

Since `H mod p` has rank one, its kernel is exactly the line spanned by
`u mod p`.  For every `x in S_c(u)`,

```text
(H/p)x is integral
iff x mod p lies in <u mod p>
iff det(x,u)=0 mod p
iff p|c.                                                  (29)
```

Thus a chosen local prime toggle acts coherently on either the entire
integral shell or none of it.  When it acts, with

```text
u'=(H/p)u,       x'=(H/p)x,                              (30)
```

one has

```text
||u'||=||u||,       ||x'||=||x||,       det(x',u')=c.    (31)
```

If additionally `v_p(||u||^2)=1`, then `u'` is primitive and the conjugate
matrix gives the inverse shell map.  A coordinate reflection or the
unit/global-conjugation gauge of THM-3353 may reverse the displayed sign of
`c`, but not its magnitude.

The valuation-one hypothesis is necessary only for primitive output, not for
the integrality equivalence.  With `u=(3,4)`, `p=5`, and the displayed `H`,
one gets `(H/5)u=(5,0)`: a coherent integral toggle may create content when
the anchor norm contains `p^2`.

The gate (29) is the ramified counterpart of Section 3A.  Away from `Ac`, a
split prime gives two simple clock branches along the shell.  When it divides
both the anchor norm and `c`, every shell point lies on the anchor's isotropic
line modulo `p`, so the same prime allocation can be toggled coherently.

The minimal positive control is

```text
p=5,       H=[ 3 4;-4 3],       u=(7,6),       u'=(9,-2).
```

For `x=(2,1)`, `det(x,u)=5` and `(H/5)x=(2,-1)` is integral.  The Farey
neighbour `x=(6,5)` has determinant one and maps to `(38,-9)/5`, not an
integral spinor.  Prime toggles therefore do not preserve a Farey star.

## 5. An affine LRC ray has two fixed U-spine resultants

Return to integers

```text
D=d+a,
p=p_0+dn,          q=q_0+Dn,
dq-Dp=c.                                                (32)
```

The determinant-shell view uses `u=(D,d)` and `x_n=(q,p)`.  A second view
treats `p` and `q` separately as indices of the consecutive spinors

```text
z_t=(t+1)+it,             C_t=N(z_t)=2t^2+2t+1.         (33)
```

THM-3346 gives their full content channels

```text
g_-(p,q)=gcd(C_p,q-p),
g_+(p,q)=gcd(C_p,p+q+1),                                (34)
gcd(C_p,C_q)=g_- g_+,        gcd(g_-,g_+)=1.            (35)
```

Define the two row fingerprints

```text
R_-=c^2+(a-c)^2,
R_+=(c+d)^2+(d+a-c)^2.                                 (36)
```

They are Gaussian norms.  More strongly, put

```text
L_-=q-p,             L_+=p+q+1,          b=2d+a.
```

The affine equation gives the two exact resultant identities

```text
a^2 C_p
 =R_-+2d L_-(dL_-+a-2c),                               (37)
b^2 C_p
 =R_++2d L_+(dL_++a-2c).                               (38)
```

Since `g_-` divides both `C_p,L_-` and `g_+` divides both `C_p,L_+`,

```text
g_-|R_-,                  g_+|R_+,
gcd(C_p,C_q)|R_-R_+.                                      (39)
```

Every prime divisor of `C_t` is `1 mod 4`: a prime `3 mod 4` dividing the
norm of the coprime coordinates `(t+1,t)` would divide both.  Therefore the
number of shared Gaussian prime-toggle coordinates along the entire infinite
ray is bounded by the fixed quantity

```text
rho(d,a,c)=# {ell:ell=1 mod 4 prime, ell|R_-R_+}.        (40)
```

The identities claim divisibility, not equality.  Primes dividing the
coefficients `a` or `2d+a`, higher valuations, and common factors of the two
fingerprints are precisely why one must retain (34), rather than replacing
it by `gcd(L_-,R_-)` and `gcd(L_+,R_+)` without qualification.

## 6. The incoming affine quotient and its exact finite fingerprint ledger

The disconnected-low Dirichlet reduction has

```text
1<=d<=8,       0<=a<=7d,       1<=|c|<=46,              (41)
```

with `c>0` at `a=0`, `c<0` at `a=7d`, and both signs in the interior.
For a formal solvable row, divide

```text
g=gcd(d,a,c),       (d,a,c)=g(d',a',c').                (42)
```

Then `gcd(d',a')=1`, the geometric affine equation is unchanged, and the old
enumeration is one residue subsequence of the primitive step
`(d',d'+a')`.  Thus the `22,890` formal rows map onto primitive rows.

For fixed primitive denominator `d`, the number of admissible directions is
`7 phi(d)` after endpoint weighting, and each contributes `2*46` signed
charges.  Hence the exact primitive count is

```text
N_d=2*46*7 phi(d)=644 phi(d),                            (43)

(N_1,...,N_8)=(644,644,1288,1288,2576,1288,3864,2576),
sum_(d=1)^8 N_d=644*22=14,168.                          (44)
```

The entire normalization fibre is also explicit.  For one primitive row put

```text
M=min(floor(8/d),floor(46/|c|)).                         (44a)
```

Its formal preimages are

```text
disjoint_union_(m=1)^M Z/mZ,                            (44b)
```

because scale `m` permits the `m` residue representatives

```text
(md,ma,mc,p_0+jd,q_0+j(d+a)),       0<=j<m.             (44c)
```

Thus one primitive row has exactly

```text
T_M=M(M+1)/2                                             (44d)
```

formal lifts.  The exact depth distribution is

```text
M : primitive rows
1 : 12,236       2 : 1,512       3 : 112       4 : 182
5 :     28       6 :    14       7 :  14       8 :  70, (44e)
```

and weighting (44e) by `T_M` gives `22,890`.  This triangular phase stack has
its own U-spine identity

```text
4T_M+1=C_M=2M^2+2M+1.                                  (44f)
```

For `M<=8` these grades are `5,13,25,41,61,85,113,145`, with at most two
distinct split primes.  The normalization-fibre grade therefore cannot
encode the rank-five endpoint fingerprint in (49); the two prime-toggle
scales are genuinely different.

The `a=1` slice has `8*92=736` rays and anchor direction

```text
u=(d+1,d)=z_d;                                          (45)
```

it is literally parallel to the U-spine.  The `17` low-unit carrier rows are
exactly the positive-cone pieces with `d+(d+a)<8` and `|c|=1`; by Section 2
each is a single Farey-star orbit.

Factoring (36) on all `14,168` primitive rows gives the FINITE-EXACT
distribution

```text
rho : number of rows
  0 :   25
  1 :  229
  2 : 4214
  3 : 7168
  4 : 2371
  5 :  161.                                             (46)
```

Thus every incoming affine tail has shared U-spine prime-toggle rank at most
five.  This is sharp, not merely a fingerprint upper bound.  On

```text
(d,a,c)=(1,1,-44),       q=2p-44,
R_-=3961=17*233,          R_+=3965=5*13*61,             (47)
```

the CRT choice

```text
(p,q)=(14,426,006,28,851,968)                           (48)
```

satisfies

```text
g_-=3961,       g_+=3965,
gcd(C_p,C_q)=15,705,365=5*13*17*61*233.                (49)
```

The coherent shell gates of Section 4 form a much smaller and differently
typed ledger.  A gate prime must divide both `c` and the anchor norm
`D^2+d^2`.  Exactly `1,166` primitive rows have such a prime, distributed as

```text
p : rows
5 : 909,   13 : 135,   17 : 68,   29 : 24,   37 : 14,   41 : 16. (49a)
```

No row has two coherent gate primes, because the two smallest possible ones
already have product `5*13>46>=|c|`.  Of these rows, `978` have anchor
valuation one and `188` have deeper valuation, where primitive content can
collapse.  This matching/isolates statement concerns simultaneous integral
shell transport only; it does not bound the endpoint U-spine clock rank (46).

## 7. Exact boundary: arithmetic organization is not LRC transport

The theorem supplies three useful structures:

1. a canonical divisor/orbit split of every affine ray;
2. a tame CRT prime-clock atlas for every determinant shell; and
3. a fixed rank-at-most-five U-spine fingerprint on the current LRC tail.

None controls the physical overlap integral.  Even a coherent Gaussian
toggle changes the coordinate data used by the carrier inequality.  For
example the norm-`85` factor toggle relates the primitive directions

```text
(D,d)=(7,6)   and   (9,2).                              (50)
```

At `|c|=5` it is integral by (29), but the high-chamber error coefficients
are respectively

```text
12*6/35+10=422/35,       12*2/35+10=374/35.             (51)
```

For the current strict target these give different cutoffs, `1898` and
`1682`.  Norm, Lorentz shell, and prime allocation therefore do not preserve
the LRC carrier certificate.  The physical context, translated cells,
component labels, and exact mass engine remain load-bearing.

Likewise, (46) does not make the open tail finite by itself: a fixed finite
prime fingerprint can occur along infinitely many parameter values, and a
failed sufficient certificate is not an unsafe LRC row.  THM-3355 remains
reserved and unproved.

## 8. Exact evidence awaiting freeze

The provisional companion

```text
04-computation/primitive_affine_determinant_shells_thm3356.py
```

currently checks:

- all `22,890` formal rows and their `14,168` primitive images;
- `566,720` shell, Lorentz, content, resultant, and two-channel evaluations;
- the totient and triangular-fibre counts, `17` Farey rows, `736`
  U-spine-parallel rows, the `1,166` coherent-gate rows, (46), and the sharp
  witness (47)--(49);
- `5,022` tame root atlases, `13,236` Gaussian allocation checks, and
  `39,240` ordered root-pair channel checks; and
- `82,400` coherent-toggle integrality checks for every split prime below
  `100`.

Normal and optimized replays must match a frozen stored transcript before
promotion.  The finite sweeps certify their declared universes; equations
(3)--(40) carry the universal proof.
