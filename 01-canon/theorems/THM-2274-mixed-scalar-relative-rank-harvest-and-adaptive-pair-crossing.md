---
id: THM-2274
title: "Mixed scalar relative-rank harvest and adaptive pair crossing"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every live scalar
  five-unit/three-blocker cover has two independent scalar relations of
  coefficient height at most 2116. Starting from THM-2275's primitive
  height-20 relation, support at least three is handled by a degree-35
  mixed relative Selberg packet. If the relation has support two, smoothing
  across the adaptive cut formed by those two coordinates gives an
  independent crossing relation: height 708 when the pair contains the
  guard and height 2116 otherwise. All 36 pair types are covered. The fixed
  THM-2203 section lifts rank two to the original row at height 4232. This
  removes THM-2275's dependency-ray exception for rank, but excludes no
  scalar profile and does not prove LRC(14).
source: codex-2026-07-25-mixed-scalar-relative-rank
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-2085-explicit-height-57-rank-seven-selberg-gate
  - THM-2137-deep-scalar-tail-boundary-complexity
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
  - THM-2164-relative-packet-rank-harvesting
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2275-mixed-scalar-relation-and-guard-blocker-crossing
related:
  - THM-2199-effective-positive-subspace-rank-lift
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation
script: 04-computation/lrc14_mixed_scalar_relative_rank_thm2274.py
output: 05-knowledge/results/lrc14_mixed_scalar_relative_rank_thm2274.out
script_sha256: c0e4faeabd2fc1cdaef671662a5885abc6742272b2bfa310fece529ed6caedf6
output_sha256: 689787c7a78a1b5599f1e024de22ddd70ba944b3e4067aebef36d87686748490
hash_basis: LF-normalized working-tree bytes
---

# THM-2274 -- mixed scalar relative-rank harvest

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2275 supplies one height-20 relation directly on the nine scalar
coordinates. This theorem forces a second relation on the same fixed
coordinate section:

```text
support at least three
  -> retain the complete known-relation Fourier face
  -> independent relation of height at most 35;

support two
  -> use the support as an adaptive 2+7 cut
  -> independent relation of height at most 708 or 2116.        (1)
```

The support-two branch is not treated as an exceptional dependency ray.
It supplies the cut on which a second relation must cross.

## 1. Scalar statement

Use the live scalar hypotheses of THM-2275. Thus

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3) in Z_(>0)^9,         (2)
```

where

```text
H,q_1,...,q_5 are thirteen-units;
H is odd;
q_1,...,q_5 are pairwise distinct;
c_1,c_2,c_3 are distinct positive multiples of thirteen, (3)
```

and, outside a null set,

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j).             (4)
```

Here

```text
C_H={t:||Ht||>1/7},          D_a={t:||at||<1/14}.     (5)
```

Let

```text
Lambda_*={r in Z^9:r.w_*=0},

W_K^*=span_Q{r in Lambda_*:||r||_infinity<=K}.        (6)
```

Then

```text
dim_Q W_2116^*>=2.                                    (7)
```

More precisely, THM-2275 gives a primitive

```text
0!=p in Lambda_*,              ||p||_infinity<=20.    (8)
```

The second relation can be chosen with height

```text
35,       if |supp(p)|>=3;

708,      if |supp(p)|=2 and the H-coordinate is in supp(p);

2116,     if |supp(p)|=2 and the H-coordinate is not in supp(p). (9)
```

Replacing THM-2275's relation by its primitive generator preserves (8).

## 2. The mixed signed tensor

The null safe event complementary to (4) has one guard coordinate of
interval length `5/7` and eight ordinary coordinates of interval length
`6/7`. Write their indicators as

```text
chi_0=1_{||x||>1/7},
chi_i=1_{||x||>=1/14},             1<=i<=8.           (10)
```

Endpoint choices do not affect Haar measure.

For degree `K`, choose the Vaaler lower and upper polynomials from
THM-2085:

```text
L_i<=chi_i<=U_i,          D_i=U_i-L_i,
epsilon=1/(K+1),          d=2epsilon.                 (11)
```

Their zero coefficients are

```text
u_0=5/7+epsilon,          u=6/7+epsilon,
Fourier(D_i,0)=d.                                      (12)
```

Use the signed tensor minorant

```text
P_K(x)
 =product_(i=0)^8 U_i(x_i)
  -sum_(i=0)^8 D_i(x_i) product_(j!=i)U_j(x_j).       (13)
```

THM-2085's product telescope gives

```text
P_K(x)<=product_(i=0)^8 chi_i(x_i)                    (14)
```

away from endpoints. Its zero coefficient is

```text
B_K
 =u_0 u^8(1-d/u_0-8d/u).                              (15)
```

This tensor is signed. Positivity is required only after the complete
known-relation face is summed; no termwise domination is asserted.

## 3. A coefficient-free face majorant

For every nonzero integer `n`,

```text
|Fourier(chi_0,n)|
 =|sin(2pi n/7)|/(pi|n|),

|Fourier(chi_i,n)|
 =|sin(pi n/7)|/(pi|n|),             i>=1.            (16)
```

THM-2164's elementary residue bounds imply in both cases

```text
|Fourier(chi_i,n)|<5/(16|n|).                         (17)
```

The Vaaler order relations also give

```text
|Fourier(U_i,n)|<=|Fourier(chi_i,n)|+epsilon,
|Fourier(D_i,n)|<=d.                                  (18)
```

Suppose that the primitive relation `p` in (8) has support size

```text
s>=3,                                                 (19)
```

and let

```text
g=1 if the guard coordinate belongs to supp(p),
g=0 otherwise.                                        (20)
```

At the frequency `ell p`, every active coordinate has magnitude at least
`ell`. Put

```text
A_ell=5/(16ell)+epsilon,
Z_(s,g)=u_0^(1-g) u^(8-s+g).                          (21)
```

Form the entire coefficient of (13) first, then take absolute values.
Equations (18)--(21) give

```text
|Fourier(P_K,ell p)| <= C_ell(s,g),                   (22)

C_ell(s,g)
 =Z_(s,g)[
      A_ell^s
     +s d A_ell^(s-1)
     +d A_ell^s(
          (1-g)/u_0+(8-s+g)/u
       )
   ].                                                  (23)
```

The three terms respectively price the all-upper product, a defect on an
active coordinate, and a defect on a zero coordinate.

Assume for contradiction that every relation in the degree-`K` tensor
box is proportional to `p`. The orbit average of `P_K` then sees only the
frequencies

```text
0, +/-p, +/-2p, ...
```

that remain inside the tensor box. Overcounting through `ell=K` yields

```text
integral_(R/Z) P_K(w_*t)dt
 >=B_K-2sum_(ell=1)^K C_ell(s,g).                     (24)
```

This is the mixed-interval analogue of THM-2164's whole-face packet. Any
additional large relation outside the tensor box is irrelevant because it
is not a Fourier frequency of `P_K`.

## 4. Exact degree-35 positivity

At

```text
K=35,             epsilon=1/36,                       (25)
```

the exact companion evaluates (24) for all thirteen valid types

```text
3<=s<=9,        g in {0,1},        s-g<=8.            (26)
```

The unique worst type is

```text
(s,g)=(3,1),                                          (27)
```

and its lower bound is

```text
10912708836373489079295440740131626642463500500311420793
----------------------------------------------------------------
1804923419549521964583407381958492581322740462946091008000

>1/166>0.                                             (28)
```

Therefore the left side of (24) is positive. By (14), the exact scalar
safe event would have positive measure, contradicting (4). Hence an
independent relation of height at most `35` exists.

The adjacent degree is an honest boundary only for the displayed
coefficient-free majorant. At `K=34`, the same worst type gives

```text
-2764922237396438494134094765517911935722912519
---------------------------------------------------
3044225083073044036025020827143413039104000000000

<0.                                                    (29)
```

Equation (29) does not say that the actual row or a sharper
residue-sensitive packet fails at degree `34`.

## 5. Support two becomes the adaptive cut

It remains to treat

```text
|supp(p)|=2.                                          (30)
```

Let `A` be the two labelled coordinates in `supp(p)` and let `B` be the
seven-coordinate complement. Since `p` is supported on `A` and is a
relation,

```text
L_A(p):=sum_(i in A)p_i(w_*)_i=p.w_*=0.               (31)
```

Thus any relation `r` satisfying

```text
L_A(r)!=0                                             (32)
```

is automatically independent of `p`.

Use THM-2275's mixed form of the positive crossing lemma. Smooth every
safe interval with the normalized squared-Fejer polynomial at bandwidth
`N`. Each coordinate has

```text
L1 error <3/(2N),          Fourier degree 2N-2.        (33)
```

For a two-coordinate block and a seven-coordinate block, the total
telescope errors are

```text
eta_A=3/N,              eta_B=21/(2N).                (34)
```

If the exact block-safe masses have lower bounds `alpha,beta`, a crossing
relation follows whenever

```text
(alpha-eta_A)(beta-eta_B)>eta_A+eta_B.                (35)
```

The extracted relation satisfies (32).

### 5.1 The pair contains the guard

The other coordinate of `A` is ordinary. The union bound gives

```text
mu(A-safe)>=5/7+6/7-1=4/7.                            (36)
```

The seven coordinates in `B` are distinct ordinary speeds. Indeed, the
five `q_i` are pairwise distinct, the three `c_j` are pairwise distinct,
and a thirteen-unit cannot equal a positive multiple of thirteen.
THM-1221 therefore gives

```text
mu(B-safe)>=15/154.                                   (37)
```

At `N=355`, the exact margin in (35) is

```text
(4/7-3/355)(15/154-21/710)-27/710
 =21177/135854950
 >0.                                                  (38)
```

The resulting crossing relation has height

```text
2N-2=708.                                             (39)
```

At `N=354`, the same ledger has margin

```text
-3/15010072<0.                                        (40)
```

This is the adjacent boundary of (35), not an impossibility claim.

### 5.2 The pair omits the guard

Now `A` consists of two distinct ordinary speeds. THM-1166 gives danger
overlap at least `1/91`, so

```text
mu(A-safe)
 =1-2/7+mu(D_a intersection D_b)
 >=66/91.                                             (41)
```

The complement `B` consists of the guard and six distinct ordinary
speeds. THM-2137's exact odd-guard six-comb floor gives

```text
mu(B-safe)>=delta_6=191/6930.                         (42)
```

At `N=1059`, equation (35) has exact margin

```text
(66/91-3/1059)(191/6930-21/2118)-27/2118
 =44021/78582173670
 >0.                                                  (43)
```

The crossing relation has height

```text
2N-2=2116.                                            (44)
```

At `N=1058`, the same ledger gives

```text
-861505/47060301288<0.                                (45)
```

Again, this is certificate-relative.

There are exactly

```text
8
```

two-coordinate supports containing the guard and

```text
binomial(8,2)=28
```

supports omitting it. Thus Sections 5.1--5.2 cover all

```text
8+28=36=binomial(9,2)                                 (46)
```

scalar pair types. Possible equality between `H` and an ordinary scalar
coefficient creates no loophole in Section 5.1: the union-bound floor
`4/7` remains valid.

## 6. Fixed-section original-row lift

THM-2203 identifies the scalar labels with the original-row coordinate
section

```text
S_I=(8H,16q_1,...,16q_5,16c_1,16c_2,16c_3).         (47)
```

The relation map is

```text
(r_H,r_1,...,r_8)
 |->(2r_H,r_1,...,r_8)                               (48)
```

on that section, with zeros on the four forgotten coordinates. It is
injective and sends scalar relations to original-row relations. Therefore
it preserves linear independence.

The height-20 first relation lifts to height at most `40`. The second
relation heights in (9) lift respectively to

```text
70,          1416,          4232.                    (49)
```

Hence the original thirteen-speed row has two independent relations
supported on the fixed nine-coordinate section by uniform height

```text
4232.                                                 (50)
```

The factor two in (48) is load-bearing. No original-row height `2116`
claim is made.

## 7. What this improves, and what it does not

THM-2275 combined its height-20 relation with THM-2266's small pair atlas.
That gave rank two by height `9841` on `120` interior profiles except when
the two relations were proportional along explicit guard-owner or
owner-unit rays.

The present theorem removes that rank exception and the interior-profile
restriction:

```text
every live scalar cover
  -> scalar relation rank at least two by height 2116
  -> fixed-section original relation rank at least two by height 4232.
                                                               (51)
```

This is much smaller and more localized than THM-2199's universal ambient
rank-twelve bound. It is also logically different from THM-2164's ambient
height-105 rank-two result: both relations here live on the nine-coordinate
scalar section, so they survive the quotient and fixed-section lift.

Nevertheless, relation rank two is not an LRC contradiction. The theorem
does not:

1. bound the scalar speeds;
2. select a blocker owner or root sheet;
3. retain THM-2263's expiration target;
4. force a nonzero minor on a prescribed pair of labels;
5. exclude any of the `165` remaining scalar profiles.

The connection and loss ledger is

```text
source:
  one quotient-faithful height-20 scalar relation and a null mixed
  safe event;

target:
  two independent relations on the same nine-coordinate section;

map:
  retain every multiple of a support-at-least-three relation as one
  Fourier face, or turn a support-two relation into the cut that a
  new positive spectral relation must cross;

preserved:
  scalar labels, the fixed coordinate section, linear independence,
  and the distinction between guard and ordinary safe intervals;

destroyed:
  root digits, owner labels, exact crossing frequency, and all
  post-expiration geometry;

cheapest hostile probes:
  support three containing the guard, a guard pair at N=354/355,
  and an ordinary pair at N=1058/1059;

needed sidecar:
  an anchored minor, bounded cofactor/carry termination, or a labelled
  owner transition that converts fixed-section rank two into a phase
  contradiction.                                             (52)
```

## 8. Reproduction

Run

```bash
python 04-computation/lrc14_mixed_scalar_relative_rank_thm2274.py
python -O 04-computation/lrc14_mixed_scalar_relative_rank_thm2274.py
```

Both modes reproduce the stored transcript after platform newline
normalization. The companion uses exact `Fraction` arithmetic and explicit
raising checks. It verifies:

```text
all thirteen support/guard types at degrees 34 and 35;
the exact worst type and fractions (28)--(29);
the complete 8+28 scalar pair partition;
the four adjacent crossing margins (38), (40), (43), (45);
the scalar heights 708 and 2116;
the fixed-section height 4232.                             (53)
```

Independent audit separately reconstructed the relative coefficient
majorant, checked all pair types and distinctness hypotheses, recomputed
the four crossing margins, and verified that (48) preserves rank. QED.
