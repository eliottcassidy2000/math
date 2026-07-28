---
id: THM-2826
title: "Uniform triple-pole Faber obstruction and simple-pole chamber exclusion"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every normalized reduced Faber degree m=4R-2 with R>=7,
  a nonsplit polynomial exact-square-prefix chart cannot contain a finite
  response point with (ord V,ord M)=(3,1).  The all-degree proof replaces
  THM-2823's finite resultants by the exact carrier
  H_R in Q[T,s,Td], a Phi/H/Psi mod-three pure-face cycle, and THM-2760's
  coprime exact-prefix pair.  In THM-2796's balanced normal form every
  simple pole part creates this local pair.  Thus every balanced passport
  containing a simple pole, including the full h>N/2 chamber, is excluded
  from every such normalized degree.  All-parts-at-least-two passports,
  other charts, JC(2), and DC(2) remain open.
source: root/sextic-uniform-triple-pole-faber-obstruction-2026-07-28
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2823-degree-twenty-six-triple-pole-faber-valuation-obstruction
  - THM-2817-sextic-e3-maximal-pole-power-chebyshev-accessory-classification
script: 04-computation/jc_uniform_triple_pole_faber_obstruction_thm2826.py
output: 05-knowledge/results/jc_uniform_triple_pole_faber_obstruction_thm2826.out
script_sha256: 6e8ba0307fb6a494641a4d876469dc01c32400efb495510785d8dfdbfb489c19
output_sha256: 7b8f50528328647f08a4feb9816017ae4baa25a495d57bd3fb1078fcad1e9173
hash_basis: LF-normalized bytes
---

# THM-2826 -- the triple-pole obstruction is uniform in Faber degree

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2823 proved the first case `m=26` by a complete Newton split, two
cubic resultants, and one exceptional two-scale unit ideal.  Those finite
eliminations conceal an all-degree mechanism.  The third response always
lands in a three-generator valuation cone, the regular boundary rotates
through three different response channels, and the polar boundary is
exactly the coprime prefix pair already proved in THM-2760.

## 1. Uniformly normalized Faber bank

Fix `R>=7`, and put

```text
m=4R-2,                         alpha=R-1/2.           (1)
```

After scaling the top coefficient, a full reduced representative in this
congruence class starts as

```text
Q=E_m+a_(R-1)E_(m-4)+sum_(j=1)^(R-2) a_j E_(4j-2).  (2)
```

The legal constant target translation `P_c=P+c` has triangular law

```text
E_m(P)=E_m(P_c)-alpha c E_(m-4)(P_c)+lower rows.      (3)
```

Since `alpha!=0`, choosing `c=a_(R-1)/alpha` kills the next row.  We work
in the resulting normalized bank

```text
Q=E_(4R-2)+sum_(j=1)^(R-2) a_j E_(4j-2).             (4)
```

This missing `R-1` row is load-bearing only for strict separation between
the top regular face and the lower bank.

Use the nonsplit polynomial exact-square-prefix chart

```text
H=Vz^2+B_src z+C_0,              L=A_src z+E_0,
P=H^2+L,                         U^2=V,                (5)
```

and its standard invariants

```text
q=A_src/U,                       T=q^2=A_src^2/V,
d=C_0-B_src^2/(4V),
s=A_src B_src/(2V)-E_0.                                  (6)
```

The inherited Hamiltonian identities are

```text
Phi_Q=0,                         Psi_Q in C,
K_Q=R_Q/q,                       A_src K_Q=lambda M,   (7)
```

where `lambda!=0` and `M` is the response carrier.

## 2. The all-degree third-response carrier

For the row indexed by `j>=1`, let

```text
A(t)=1+2dt^2+qt^3+(d^2-s)t^4
    =(1+dt^2)^2+t^3(q-st),

A(t)^(j-1/2)=sum_(n>=0)c_n^(j)t^n,                    (8)

phi_j=Phi_(4j-2)/q=4c_(4j-1)^(j)/q,
psi_j=Psi_(4j-2)=4c_(4j)^(j),
K_j=(4c_(4j+1)^(j)+2d c_(4j-1)^(j))/q.               (9)
```

Coefficient extraction gives

```text
K_j=-(d/2)phi_j+T H_j,                                (10)

H_j=4q^(-3)[t^(4j+1)](1+dt^2)A(t)^(j-1/2).           (11)
```

The crucial point is not merely that `(11)` is polynomial.  Expand

```text
(1+dt^2)A(t)^(j-1/2)
 =sum_(u>=0) binom(j-1/2,u)t^(3u)(q-st)^u
              (1+dt^2)^(2j-2u).                      (12)
```

Every summand with `u<=j` has degree at most `4j`, so only
`u=j+1+h` contributes to `(11)`.  If `ell` copies of `-st` and `k`
copies of `dt^2` are chosen, the coefficient condition is

```text
ell+2k=j-2-3h.                                        (13)
```

Moreover

```text
q^(u-ell-3)=q^(2k+4h)=T^(k+2h).                       (14)
```

Therefore the exact finite formula is

```text
H_j =
 4 sum_(h,k,ell>=0; ell+2k=j-2-3h)
   binom(j-1/2,j+1+h) binom(j+1+h,ell)(-1)^ell
   binom(-2-2h,k) T^(2h)(Td)^k s^ell.                 (15)
```

In particular,

```text
H_j in Q[T,s,Td]                                     (16)
```

for every `j`; this is an identity, not a finite-degree observation.  Put

```text
H_Q=H_R+sum_(j=1)^(R-2)a_jH_j.                        (17)
```

On `Phi_Q=0`, equations `(10)` and `(17)` become

```text
K_Q=T H_Q.                                            (18)
```

## 3. Local valuation ledger at a triple pole

Let `beta` be a finite point with

```text
ord_beta(V)=3,                    ord_beta(M)=1.       (19)
```

Write

```text
a=ord_beta(A_src),                b=ord_beta(B_src),  (20)
```

allowing `b=infinity` when `B_src=0`.  Polynomiality makes `a,b`
nonnegative, and `(7)` makes `A_src` nonzero.  With `v=ord_beta`,

```text
v(T)=2a-3,                        v(K_Q)=1-a,
v(d)>=min(0,2b-3),               v(s)>=min(0,a+b-3). (21)
```

Whenever a displayed minimum is negative it is an equality, because a
regular source term cannot cancel a pole.

## 4. The two regular `H` regions

If `a>=3`, then

```text
v(T)>=3,                          v(d)>=-3,
v(s)>=0,                          v(Td)>=0.            (22)
```

Every monomial in `(15)` is regular, hence `v(H_Q)>=0`.  Equation `(18)`
gives `v(K_Q)>=3`, contrary to `v(K_Q)=1-a<=-2`.

If `a=2,b>=1`, then

```text
v(T)=1,                           v(d)>=-1,
v(s)>=0,                          v(Td)>=0.            (23)
```

Again `H_Q` is regular, so `(18)` gives `v(K_Q)>=1`, contrary to
`v(K_Q)=-1`.

This is the first uniform replacement for THM-2823's degree-specific
polynomial `(8)`.

## 5. The regular mod-three face cycle

It remains in this section to take `a<=1,b>=2`.  Then `d` is regular.
Only `(a,b)=(0,2)` permits a pole in `s`, and there `v(s)=-1`.

Exactly one of the following pure-`q` terms occurs in the top row:

```text
R=3k+1:  phi_R has
          4 binom(R-1/2,4k+1) T^(2k);

R=3k+2:  H_R has
          4 binom(R-1/2,4k+3) T^(2k);

R=3k:    psi_R has
          4 binom(R-1/2,4k) T^(2k).                  (24)
```

All three half-integer binomial coefficients are nonzero.

These terms are the unique least-valuation faces.  For the hostile boundary
`(a,b)=(0,2)`, assign

```text
v(q)=-3/2,                        v(s)=-1,
v(d)>=0.                                               (25)
```

A monomial contributing to `t^n` with `j_q` copies of `qt^3` and
`ell` copies of `-st^4` has valuation at least

```text
-3j_q/2-ell >= -n/2,                                (26)
```

with equality only for the pure-`q` choice.  The other regular subcases
have `s` regular, so the same uniqueness is stronger.  Every retained row
has index at most `R-2`; its relevant coefficient occurs at least eight
`t`-degrees earlier and therefore has strictly larger valuation.  Thus no
lower coefficient in `(4)` can cancel `(24)`.

For `R=3k+1`, the first term in `(24)` contradicts `Phi_Q=0`.  For
`R=3k`, the third has a pole and contradicts `Psi_Q in C`.  For
`R=3k+2`, the middle term gives a pole in `H_Q`, whereas `(18)` and
`(21)` require

```text
v(H_Q)=v(K_Q)-v(T)=4-3a>=1.                          (27)
```

This proves the regular region in every congruence class.  The rotating
`Phi/H/Psi` certificate, rather than one fixed response, is why the
mod-three pattern is structural.

## 6. Polar regions are the THM-2760 prefix

The only remaining pairs are

```text
(a,b)=(0,0),(1,0),(2,0),(1,1),(0,1).                 (28)
```

Pass to the harmless local square-root extension and put

```text
omega=B_src/(2U),                 q=A_src/U.           (29)
```

The polar initial forms of `(6)` are exactly

```text
in(d)=-omega^2,                   in(s)=q omega.       (30)
```

Set

```text
rho=q/omega^3.                                        (31)
```

Equations `(30)` are precisely the exact-prefix specialization of
THM-2760.  Its all-degree polynomials `P_R,Q_R` give

```text
in(phi_R)=4 in(s)^(R-1) P_R(rho),
in(psi_R)=4 in(s)^R Q_R(rho),                         (32)

gcd(P_R,Q_R)=1,                    P_R(0)Q_R(0)!=0.    (33)
```

Since `v(s)<0` throughout `(28)` and the next retained index is `R-2`,
the top row in `(32)` strictly dominates every lower row.

For `b=0`, `a=0,1,2`, one has `v(rho)=a+3>0`.  For `(a,b)=(1,1)`,
one has `v(rho)=1`.  Thus `(33)` makes both top initial fluxes nonzero.
The first contradicts `Phi_Q=0`, and the second has negative valuation
and contradicts `Psi_Q in C`.

For the exceptional pair `(a,b)=(0,1)`, `rho` is a nonzero unit.
If both top initial fluxes vanished, their residue would be a common root
of `P_R,Q_R`, contrary to `(33)`.  Hence at least one of the two flux
equations fails.  This replaces THM-2823's exceptional Groebner basis and
quintic resultant by one cited all-degree coprimality theorem.

Sections 4--6 exhaust every nonnegative pair `(a,b)`, including
`b=infinity`.  We have proved:

> **Uniform local obstruction.**  For every `R>=7`, no normalized
> nonsplit polynomial exact-square-prefix degree `4R-2` chart contains a
> finite point with `(ord V,ord M)=(3,1)`.

## 7. The simple-pole chamber in every balanced passport

THM-2796 writes every balanced response carrier as

```text
V=v S D T^2,                     M=S E T,              (34)

D=product_(j=1)^h(x-beta_j)^(p_j),
T=product_(j=1)^h(x-beta_j),
sum_j p_j=N,                     p_j>=1,               (35)
```

with `S,E,T` pairwise disjoint and squarefree.  At a pole `beta_j`,

```text
(ord_beta_j V,ord_beta_j M)=(p_j+2,1).                (36)
```

Consequently every simple part `p_j=1` supplies `(19)` and is excluded by
the local theorem.  A simple part is odd, so THM-2796's squareclass gate
also guarantees that this application is genuinely nonsplit.

There is an immediate passport-level consequence.  If `h>N/2` and every
positive part were at least two, then

```text
N=sum_j p_j >=2h>N,                                  (37)
```

which is impossible.  Hence:

> **High-pole chamber exclusion.**  Every balanced passport with
> `h>N/2` is excluded from every normalized reduced degree `4R-2`,
> `R>=7`.

In particular, THM-2796 has `N=s+2e` and `h<=e+1`.  Its extremal
zero-simple layer

```text
s=0,                            N=2e,
h=e+1                                             (38)
```

lies in `(37)` and is empty in every such degree without any accessory
classification.  The power/Chebyshev sextic examples of THM-2817 and the
two- and three-pole `e=2` charts are special cases of this general
simple-pole statement, not separate inputs.

## 8. Scope and sharp residual

This theorem is an entry obstruction for the inherited nonsplit polynomial
exact-square-prefix chart.  It does not say that the abstract rational
response maps of THM-2796 do not exist.  It does not enter an arbitrary
Keller pair into this chart, treat nonpolynomial prefixes, or prove
`JC(2)` or `DC(2)`.

The local multiplicity is also sharp for the present proof.  Balanced
passports whose pole parts all satisfy `p_j>=2` have
`ord_beta(V)>=4`; the valuation cone changes, and neither `(22)` nor the
pure-face comparison has been proved there.  Thus the next local target is

```text
(ord V,ord M)=(4,1),                                  (39)
```

not another enumeration of the already-eliminated simple-pole chamber.

## 9. Exact companion and proof/computation boundary

The companion is finite verification of the symbolic proof, not the source
of its all-degree quantifier.  It:

1. reconstructs the three relevant coefficients for `R=1,...,18` both by
   the Faber recurrence and by an independent generalized multinomial sum;
2. checks `(15)` through `R=18` and verifies every monomial lies in
   `Q[T,s,Td]`;
3. checks the two regular `H` bounds and the strict `Phi/H/Psi` face cycle
   against every retained row for `R=7,...,18`;
4. reconstructs the THM-2760 `P_R,Q_R` formulas directly and checks their
   finite gcd controls for `R=7,...,18`;
5. checks exhaustiveness of the local valuation partition and all
   high-pole integer partitions through `N=24`; and
6. contains neither a Python `assert` node nor a float literal.

The all-degree proof comes from formulas `(12)--(15)`, the general
weighted-face argument `(24)--(27)`, and proved THM-2760—not from those
finite ranges.  Run

```text
python 04-computation/jc_uniform_triple_pole_faber_obstruction_thm2826.py
python -O 04-computation/jc_uniform_triple_pole_faber_obstruction_thm2826.py
```

Normal, optimized, and stored transcripts agree exactly.
