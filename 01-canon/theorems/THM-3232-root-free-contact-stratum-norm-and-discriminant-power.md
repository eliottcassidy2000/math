---
id: THM-3232
title: "Root-free contact-stratum norm and discriminant power"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Every nonempty first-contact stratum of the THM-3229 squarefree
  Hasse--Pluecker gcd flag carries a canonical etale unit c_m.  Its norm is
  the nonzero root-free resultant quotient
  (-1)^r Res(S_m,Omega_m)/Res(S_m,f'g'), equal after splitting to the
  unordered product of all selected-root first-live coefficients.  The norm
  has exact affine weight r(m-1); combining it with Disc(S_m) produces an
  affine- and row-swap-invariant nonzero scalar.  For f=P and g=P+P^m this
  norm is exactly Res(P,P')^(m-1), a signed discriminant power.  The norm
  removes root choice but forgets every branch label and additive prime carry.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent rational companion checks two split degree-two
  contact strata branchwise, five irreducible quadratic power-family norms,
  the resultant sign, trace cancellation, cubic affine and discriminant
  weights, the balanced invariant, diagonal/swap/general-GL2 frame laws, an
  irreducible F5 Frobenius norm, and the coarse-stratum and repeated-root
  boundaries.  Normal and optimized runs byte-match the stored transcript
  and the LF-normalized hashes below.  An independent hostile audit rederived
  the etale-unit typing, resultant sign, split/Galois product, target-frame
  law, affine and discriminant weights, balanced invariant, finite-field
  norm, and both sharp boundaries, and found no defect.
depends_on:
  - THM-3229-hasse-pluecker-simple-root-contact-gcd-flag-and-degree-termination
related:
  - THM-3221-selected-root-osculating-separation-and-minimal-jet-prime-carry
  - THM-3178-squarefree-resultant-tangent-cone-and-first-witt-norm
  - THM-2067-a-galois-orbit-product-proof-of-one-variable-constant-term-nonvanishing
  - THM-2022-gmc2-frobenius-lowest-balanced-face
script: 04-computation/contact_stratum_norm_discriminant_thm3232.py
output: 05-knowledge/results/contact_stratum_norm_discriminant_thm3232.out
script_sha256: 350f6ea1adfa4cabe79d1e8b1d09b8d8982363c1ab66dc3304ea4296b2e572ff
output_sha256: 0690f8cfa91548929be01842dd3585a818bd9394a1f1aff05f49620f9990a256
hash_basis: LF-normalized bytes
---

# THM-3232 -- root-free contact-stratum norm and discriminant power

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3229 replaces a list of selected common roots by a squarefree divisor
flag, but its first-live coefficient is still stated one root at a time.  A
finite etale norm performs the next quotient: it retains exact simultaneous
nonvanishing on one contact stratum while removing every root ordering and
Galois choice.

The resulting scalar is multiplicative, not additive.  That distinction is
load bearing: traces of the live root coefficients can cancel exactly while
their norm cannot.

## 1. The contact algebra and its unit

Let `K` be a perfect field and let `f,g in K[x]` be nonzero and
nonproportional.  Use THM-3229's squarefree simple-common-root flag

```text
G_1=G_0,                 G_m=gcd(G_(m-1),Omega_m),
Omega_m=f^[m]g'-g^[m]f'.                                (1)
```

Fix `m>=2` for which the first-contact quotient

```text
S_m=G_(m-1)/G_m                                           (2)
```

has degree `r>=1`, and make `S_m` monic.  It is squarefree, so

```text
A_m=K[x]/(S_m)                                            (3)
```

is a finite etale `K`-algebra of dimension `r`.  Both `f'g'` and `Omega_m`
are units in `A_m`: the first assertion is the simple-root saturation in
`G_0`, and the second is exactly the quotient by `G_m`.

Define the normalized contact class

```text
c_m=-Omega_m/(f'g') in A_m^*.                            (4)
```

At every geometric root `a` of `S_m`, THM-3229 gives

```text
c_m(a)
 =[u^m] g(a+u)/g'(a)-[u^m] f(a+u)/f'(a),                 (5)
```

and all lower normalized differences vanish.  Thus `(4)` is not a new
choice of coefficient: it is the simultaneous etale class of the already
proved selected-root coefficients.

## 2. Root-free norm and exact resultant formula

Let

```text
N_m=Norm_(A_m/K)(c_m).                                    (6)
```

The norm is nonzero because `c_m` is a unit.  With the resultant convention

```text
Res(P,Q)=lc(P)^deg(Q) product_(P(a)=0) Q(a),              (7)
```

monicity of `S_m` gives `Norm(h)=Res(S_m,h)` for every polynomial class
`h`.  Multiplicativity and `Norm(-1)=(-1)^r` therefore give the exact sign

```text
N_m
 =(-1)^r Res(S_m,Omega_m)/Res(S_m,f'g') !=0.             (8)
```

The denominator in `(8)` is nonzero by simple-root saturation.  After a
splitting extension,

```text
N_m=product_(S_m(a)=0)c_m(a).                             (9)
```

Equation `(9)` proves simultaneously that the answer is independent of root
ordering and fixed by the Galois group.  Conversely it shows the exact loss:
`N_m` says that every factor in the product is live but does not recover any
one factor, branch label, or additive relation among the factors.

There is also a polynomial Pluecker semi-invariant

```text
W_m=Norm_(A_m/K)(-Omega_m)=(-1)^r Res(S_m,Omega_m),
N_m=W_m/Res(S_m,f'g').                                    (10)
```

Thus denominator clearing never changes the nonvanishing locus on the
simple first-contact stratum.

## 3. Target-frame covariance

Constant independent row rescalings leave `c_m` and `N_m` unchanged.  More
generally, let

```text
F=alpha f+beta g,        G=gamma f+delta g,
Delta=alpha delta-beta gamma !=0,                         (11)
```

and assume `F'G'` remains a unit on `A_m`.  Since

```text
Omega_m(F,G)=Delta Omega_m(f,g),                          (12)
```

the normalized first-contact class transforms by

```text
c_m(G<-F)=Delta f'g'/(F'G') c_m(g<-f).                    (13)
```

Taking norms gives the full target-frame law.  In particular `W_m` has pure
determinant weight `r`:

```text
W_m(F,G)=Delta^r W_m(f,g).                                (14)
```

A row swap sends `c_m` to `-c_m` and hence sends `N_m` to `(-1)^r N_m`.
The derivative quotient in `(13)` is load bearing: the normalized norm is
not a general `GL_2` invariant.

## 4. Affine source weight and a balanced invariant

Let `T(y)=rho y+sigma` with `rho!=0`, put

```text
F(y)=f(T(y)),            G(y)=g(T(y)),
S~_m(y)=rho^(-r)S_m(T(y)).                                (15)
```

Then `S~_m` is monic and defines the transported stratum.  Hasse derivatives
give

```text
Omega_m(F,G)=rho^(m+1)(Omega_m(f,g) composed T),
F'G'=rho^2(f'g' composed T),                              (16)
```

so each transported root coefficient and its norm have weights

```text
c~_m=rho^(m-1)c_m,
N~_m=rho^[r(m-1)]N_m.                                    (17)
```

The monic discriminant has the opposite weight

```text
Disc(S~_m)=rho^[-r(r-1)]Disc(S_m).                        (18)
```

Consequently

```text
I_m=N_m^(r-1) Disc(S_m)^(m-1) in K^*                    (19)
```

is invariant under affine source changes.  It is also unchanged by constant
independent row rescalings and by swapping `f,g`, because `r(r-1)` is even.
Thus `(19)` is a root-free, affine-source, unordered-row scalar.

For `r=1`, `(19)` is identically one.  This is the sharp single-root
boundary: without a chosen tangent scale, a nonzero coefficient of weight
`m-1` has no nonconstant affine scalar invariant.

## 5. Exact discriminant-power family

Let `P in K[x]` be monic and squarefree and take

```text
f=P,                  g=P+P^m=P(1+P^(m-1)).              (20)
```

Every root of `P` is simple and has first normalized contact exactly `m`, so

```text
S_m=P,                 c_m=(P')^(m-1) in K[x]/(P).        (21)
```

Therefore

```text
N_m=Res(P,P')^(m-1)
   =(-1)^[r(r-1)(m-1)/2] Disc(P)^(m-1).                  (22)
```

This family proves that the norm is not merely a formal nonzero quotient: it
recovers an exact signed power of the discriminant.  For the irreducible
quadratic `P=x^2-2` over `Q`, the norms for `m=2,...,6` are

```text
-8, 64, -512, 4096, -32768.                              (23)
```

At `m=2`, the class is `c_2=2x`.  It has trace zero but norm `-8`.  Hence an
additive Galois average can erase a live stratum even in the smallest
irreducible example; multiplicative orbit aggregation is essential.

## 6. Frobenius orbit product

If `K=F_q` and `S_m` is irreducible of degree `r`, then

```text
A_m=F_(q^r),
N_m=c_m^[1+q+...+q^(r-1)] in F_q^*.                      (24)
```

For a reducible squarefree stratum, `(24)` applies on every irreducible
factor and their field norms multiply.  The exact companion checks
`P=x^2+2` over `F_5`: for the order-two power family, `c_2=2x` and both the
resultant norm and the Frobenius power equal `3`.

This is the precise bridge to the orbit-product move in THM-2067 and the
Frobenius language surrounding THM-2022.  The bridge stops one step short of
either application: here no competing valuation identity is supplied, and
no Gaussian-moment channel is proved to realize `c_m`.

## 7. Sharp boundaries

1. **Exact stratum is load bearing.**  If roots of different contact depths
   are combined before quotienting by `G_m`, `Omega_m` vanishes on the deeper
   roots and the coarse norm is zero.  The companion realizes contacts two
   and three on two roots: the norm over `G_1` vanishes at depth two, while
   the norm over `S_2=G_1/G_2` is nonzero.
2. **Multiple roots.**  A repeated `S` has zero discriminant and a non-etale
   quotient.  Formula `(9)`, the Frobenius interpretation, and `(19)` then
   fail.  THM-3178's conormal norm survives repeated resultant collisions at
   first Witt order; it is not a replacement for the higher normalized
   contact class here.
3. **Proportional rows.**  The THM-3229 flag does not terminate and there is
   no first-contact stratum.
4. **General target frames.**  Equation `(13)`, not invariance, is the exact
   law.  Zeros of `F'G'` leave the simple normalized chart.
5. **No branch or additive carry.**  The norm cannot recover which root
   carried a label, the THM-3227 trie edge, or its additive prime carry.
6. **No physical survivor.**  No Wick channel, LRC packet, Jacobian owner,
   or resonant PRS row is proved to realize this etale unit.

## 8. Cross-frontier connection contract

```text
source:      one nonempty THM-3229 first-contact divisor S_m;
map:         c_m=-Omega_m/(f'g') in A_m, then finite etale norm;
target:      nonzero resultant quotient N_m and balanced scalar I_m;
preserved:   simultaneous nonvanishing, contact depth, Galois orbit,
             constant row scaling, exact affine weight;
destroyed:   root order, branch identity, trace/additive relations,
             THM-3227 affine labels and prime carry;
sidecar:     a physical carrier plus a branch/reference selector if an
             additive or owner-labelled conclusion is required.
```

For Jacobian/PRS work, `(8)` is a bounded-degree elimination certificate for
one entire first-contact stratum.  For Gaussian moments, it packages all
Galois-conjugate root coefficients without cancellation, but NC2/GMC(2) is
already proved by THM-2022 and no new arbitrary-radial channel theorem is
claimed.  For LRC, the norm explains why an orbit product can retain
nonvanishing while still being unable to supply the semantic root/owner
intertwiner isolated by the current frontier.

## 9. Exact companion

The assertion-independent companion

```text
04-computation/contact_stratum_norm_discriminant_thm3232.py
```

uses exact rational and finite-field arithmetic and two independent norm
presentations: multiplication determinants and Sylvester resultants.  It
checks exactly:

```text
two split degree-two contact strata and four branch coefficients;
five irreducible quadratic discriminant-power norms;
trace-zero versus norm-nonzero cancellation;
cubic affine weights and the balanced scalar;
diagonal, swap, and general GL2 target-frame laws;
one irreducible F5 Frobenius orbit norm;
coarse mixed-depth and repeated-root boundaries.             (25)
```

Normal and optimized runs byte-match

```text
05-knowledge/results/contact_stratum_norm_discriminant_thm3232.out
```

and the LF-normalized hashes are pinned in the frontmatter.

QED.
