---
id: THM-3592
title: "Universal exponent-two three-by-three weight Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT / REPAIRED AFTER FAILED AUDIT; REAUDIT PENDING.  For every
  squarefree polynomial Sigma over C of degree at least two, no polynomial
  Darboux pair on c^2 e=Sigma(b), after scalar removal, can have exactly three
  nonzero homogeneous weight pieces in each output.  The proof classifies all
  three-point sumset collisions, deletes disconnected support components, and
  closes the equal-AP, diagonal, hooked, reflected-triple, and two Euclidean
  connected families.  Together with THM-3569 and THM-3583 this gives a
  universal seven-piece nonconstant support floor.  No Darboux pair, planar
  Jacobian counterexample, or solution of JC(2) is claimed.
source: kps-s188 / unequal-step three-by-three Darboux hostile lane, 2026-08-21
audit: >
  The standard-library exact companion passes ordinary and optimized Python
  with byte-identical output.  The first independent hostile audit correctly
  rejected the pushed version because simultaneous support reversal is not a
  regularity-preserving symmetry: the reflected hook and two reflected
  Euclidean representatives had been omitted.  They are now proved below and
  included in the companion.  A second hostile audit is pending.
depends_on:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3579-equal-step-three-by-three-danielewski-darboux-nonentry
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
related:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-3584-all-exponent-equal-step-three-by-three-danielewski-darboux-nonentry
script: 04-computation/jc2_universal_exponent_two_three_by_three_nonentry_thm3592.py
output: 05-knowledge/results/jc2_universal_exponent_two_three_by_three_nonentry_thm3592.out
script_sha256: 863795c6254ba1ab2e246ba078885a625724978d53834dcf3f27dfdf19753539
output_sha256: effd1c69f39035adf8f8c59d663496f96f7db2247f982d2025b6de7b9b85b894
hash_basis: LF-normalized bytes
---

# THM-3592 -- universal exponent-two three-by-three weight Darboux nonentry

**PROVED + VERIFIED-EXACT / REPAIRED AFTER FAILED AUDIT; REAUDIT PENDING.**  This closes
the unequal-step `3 x 3` cell left by THM-3579 and the last six-piece cell left
by THM-3583.  The first hostile audit caught a reversal-sensitive omission;
Sections 5.5 and 7 now close the three missing representatives, and the
independent-audit qualifier is deliberately pending until that repair is
reaudited.

All rings are over `C`.  Let `Sigma in C[b]` be squarefree with
`deg Sigma>=2`, and put

```text
A_Sigma=C[b,c,e]/(c^2 e-Sigma(b)).                       (1)
```

Use the Poisson bracket and grading

```text
{b,c}=c^2,          {c,e}=-Sigma'(b),          {b,e}=-2ce,
wt(b,c,e)=(0,1,-2).                                      (2)
```

Subtract constant weight-zero pieces before counting support.  The theorem is:

> If `P,Q in A_Sigma` each have exactly three retained homogeneous weight
> pieces, then `{P,Q}!=1`.

The word *retained* permits a nonconstant weight-zero coefficient but excludes
a scalar, which can be subtracted without changing the bracket.

## 1. Homogeneous algebra and the three rigidity gates

On `c!=0`, a weight-`u` element has the unique form `c^u f(b)`.  It extends
regularly across `c=0` exactly when

```text
Sigma^ceil(-u/2) divides f                    for u<0.   (3)
```

For homogeneous pieces define

```text
W_(u,v)(f,g)=v f'g-u f g'.                              (4)
```

Then

```text
{c^u f,c^v g}=c^(u+v+1) W_(u,v)(f,g).                  (5)
```

We use the following three elementary gates throughout.

### 1.1 Singleton zero rows

For nonzero polynomials and nonzero weights of the same sign, unique
factorization gives

```text
W_(u,v)(f,g)=0
 iff
f=A h^(|u|/delta),   g=B h^(|v|/delta),
delta=gcd(|u|,|v|),                                  (6)
```

with nonzero constants `A,B`.  With strict opposite signs, the logarithmic
derivative equation would make `f^|v|g^|u|` constant, so no retained negative
coefficient is possible.  If exactly one weight is zero, `(4)` forces its
weight-zero coefficient to be constant and hence removable.  If both weights
are zero, the row vanishes tautologically.  We call these alternatives the
singleton gate.

### 1.2 The scalar-arm gate

Let `beta` be a simple root of `Sigma`.  A complementary pair has weights
`(-R,R-1)`.  If the two coefficient orders at `beta` are `m,n`, its first
possible local term has order and multiplier

```text
m+n-1,             (R-1)m+Rn.                         (7)
```

By `(3)`, `m>=ceil(R/2)`.  A nonzero constant at the arm therefore occurs
only for

```text
R=2,              (m,n)=(1,0).                        (8)
```

The `R=1` candidate has zero multiplier.  Thus at each arm every scalar fibre
must contain at least one address with weights

```text
(-2,1)  or  (1,-2),                                    (9)
```

whose negative coefficient is simple and whose positive coefficient is a
unit at that arm.  When two such weight addresses coexist, this statement
does not by itself select one address uniformly over all arms; a separate
compatibility or no-alternation argument is then required.  This is the
scalar-arm gate.

### 1.3 Euler factors

Whenever `Sigma|h`, every nonzero polynomial `K` and every integer `a>=1`
satisfy

```text
E_a(h,K)=K h'+a h K',
deg E_a=deg h+deg K-1>=1,
lc(E_a)=(deg h+a deg K)lc(h)lc(K) !=0.                 (10)
```

Consequently `E_a(h,K)` cannot divide the unit `1`.  THM-3569 supplies the
homogeneous, one-by-arbitrary, `2 x 2`, and `2 x 3` exclusions used below.

## 2. Complete affine classification of two three-point supports

Translate the supports independently and write their positive gap vectors as

```text
A={0,x,x+y},                  B={0,u,u+v}.              (11)
```

Common dilation, exchange of `A,B`, and simultaneous reversal of the two gap
vectors are used only to name collision types.  The later arm analysis uses
the actual absolute weights and does not assume that reflection preserves
Danielewski regularity.

Let `F_s={(i,j):a_i+b_j=s}` and let `rho_A(d)` count increasing pairs in `A`
of difference `d`.  Put

```text
C(A,B)=sum_d rho_A(d)rho_B(d)
      =sum_s binom(|F_s|,2).                            (12)
```

Every fibre has size at most three.  A triple fibre is necessarily

```text
F={02,11,20},                                           (13)
```

and exists exactly when `(u,v)=(y,x)`.  Hence there is at most one triple;
writing its indicator as `tau`, double-counting lost addresses gives

```text
|A+B|=9-C(A,B)+tau.                                    (14)
```

The difference multiset of `(x,y)` is `{x,y,x+y}`, with `x` repeated exactly
when `x=y`.  Comparing these two three-term Schur multisets gives the complete
classification for `|A+B|<=7`:

| size | type | normalized gap vectors | nontrivial fibres |
|---:|---|---|---|
| 5 | equal AP | `(d,d ; d,d)` | `01=10`, `02=11=20`, `12=21` |
| 6 | diagonal `D` | `(x,y ; x,y)`, `x!=y` | `01=10`, `02=20`, `12=21` |
| 6 | hook `H`, `bar H` | `(d,d ; d,2d)`; `(d,d ; 2d,d)` | `01=10`, `11=20`, `02=21`; `01=20`, `02=11`, `12=21` |
| 7 | reflected `R` | `(x,y ; y,x)`, `x!=y` | `02=11=20` |
| 7 | Euclidean `E+`, `bar E+` | `(p,q ; p,q-p)`; `(q,p ; q-p,p)`, `0<p<q`, `q!=2p` | `01=10`, `12=20`; `02=10`, `12=21` |
| 7 | Euclidean `E-`, `bar E-` | `(p,q ; q-p,p)`; `(q,p ; p,q-p)`, same range | `02=11`, `12=20`; `02=10`, `11=20` |
| 7 | AP extension | `(d,d ; d,t)`, `t notin {d,2d}` | two doubles |
| 7 | AP contained | `(d,d ; p,d-p)`, `p!=d/2` | two doubles |
| 7 | AP dyadic | `(d,d ; 2d,2d)` | two doubles |

Here and below a semicolon separates the two gap vectors.  For completeness,
the multiset proof is short.  If neither support is an AP, all three
differences have multiplicity one.  Three common differences force equal or
reversed gap vectors.  Two common differences cannot be the two gaps in both
Schur triples, and cannot be gap-plus-total in both without forcing the third;
therefore they are the two gaps in one triple and a gap plus the total in the
other.  The missing gap is `q-p`, giving `E+` or `E-`.  If exactly one support
is an AP, multiplicities `2,1` show that `C=3` is precisely the hook and that
`C=2` is precisely extension or containment.  If both are APs, equality gives
size five and ratio two gives the dyadic size-seven case.  Thus no
triple-plus-double size-six pattern was omitted.

## 3. Component deletion and the five connected families

Define `Gamma_B(A)` on the three vertices of `B` by joining `b_j,b_l` when
`|b_j-b_l|` is a difference in `A`; define `Gamma_A(B)` symmetrically.  If

```text
a_i+b_j=a_k+b_l,                                        (15)
```

then `j,l` lie in one component of `Gamma_B(A)`.  Thus bracket rows coming
from different components are disjoint.  Exactly one component can contain
the scalar row.  Every other component has total bracket zero and can be
deleted.  If either graph is disconnected, the retained component on that
side has at most two vertices, reducing to the `3 x 2` theorem THM-3569.

For `|A+B|>=8`, there is at most one collision pair, so at least one component
graph is disconnected.  The three AP size-seven forms are also disconnected:
extension has the third non-AP vertex isolated, containment has the middle
vertex isolated, and the dyadic form has the middle vertex isolated on the
fine-scale side.  The connected list is therefore exactly

```text
equal AP, D, H/bar H, R, E+/bar E+, E-/bar E-.          (16)
```

The equal-AP case is THM-3579.  It remains to close the other five types.

## 4. The diagonal family reduces to `2 x 2`

In `D`, each scalar candidate is an off-diagonal of a `2 x 2` rectangle:

```text
01/10 has corners 00,11;
02/20 has corners 00,22;
12/21 has corners 11,22.                                (17)
```

Both other corners are globally singleton fibres.  Their brackets vanish,
while the chosen diagonal is the whole scalar fibre.  Restricting to those
two pieces of each output would therefore give a `2 x 2` Darboux pair,
contrary to THM-3569.  So `D` is empty.

## 5. The hooked six-row family

Write the relative supports as

```text
A={0,d,2d},                    B={0,d,3d}.              (18)
```

Its fibres, in increasing offset, are

```text
00;       01+10;       11+20;       02+21;       12;       22. (19)
```

A scalar cannot be a singleton by the homogeneous theorem.  Testing the two
addresses and two orientations of `(9)` in each of the three double fibres
gives twelve cases.  The singleton gate, including its zero-weight boundary,
leaves exactly these four ladders:

| scalar address | gate | step | `supp(P)` | `supp(Q)` |
|---|---|---|---|---|
| `01` | `(-2,1)` | `d=2k+1` | `{-2,2k-1,4k}` | `{-2k,1,4k+3}` |
| `10` | `(1,-2)` | `d=2k+1` | `{-2k,1,2k+2}` | `{-2,2k-1,6k+1}` |
| `11` | `(1,-2)` | `d=3k+1` | `{-3k,1,3k+2}` | `{-3k-3,-2,6k}` |
| `02` | `(-2,1)` | `d=2k+1` | `{-2,2k-1,4k}` | `{-6k-2,-4k-1,1}` |

Here `k>=1`.  The other eight gates die as follows:

| scalar address and gate | singleton obstruction |
|---|---|
| `01:(1,-2)`; `10:(-2,1)` | `00` has strict opposite signs |
| `11:(-2,1)` | `12` has strict opposite signs |
| `20:(-2,1)` | `00` has strict opposite signs |
| `20:(1,-2)` | `12` is mixed, or zero/nonzero at `d=1` |
| `02:(1,-2)` | `00` has strict opposite signs |
| `21:(-2,1)` | `12` has strict opposite signs |
| `21:(1,-2)` | `12` is mixed for `d>1`; `22` is zero/nonzero for `d=1` |

At `d=3`, the first two surviving descriptions coincide and both scalar
addresses have gate weights.  The low singleton is then a `(-2,-2)` row, so
`(6)` ties the two negative coefficients to the same base and gives them the
same arm orders.  Scalar contribution cannot alternate between arms to evade
the simple-order conclusion.

For the first, second, and fourth survivors, the singleton containing the
simple weight `-2` forces `d` odd; `d=1` has a mixed or removable zero-weight
endpoint.  For the third survivor put `a=d-1` and
`delta=gcd(a,a+3)=gcd(a,3)`.  The low singleton has

```text
f_0=A h^(a/delta),       g_0=B h^((a+3)/delta).         (20)
```

If `ell` is the arm order of `h`, cancellation in the adjacent row requires

```text
((a+3)-a)ell/delta=1.                                  (21)
```

Indeed the two initial multipliers are `a-2aell/delta` and
`-(a+3)ell/delta`; the second is nonzero, and the first could vanish only if
`delta=2ell`, impossible for `delta in {1,3}`.  Thus `delta=3`, `ell=1`, and
`a=3k`, proving `d=3k+1` rather than merely guessing it from a finite census.

Use `f_i,g_j` for the coefficients at the three ordered weights of `P,Q`.
We now close the four ladders.

### 5.1 First ladder, primary orientation

Put

```text
p=2k-1,                     q=4k+3.                    (22)
```

The low singleton and scalar gate give

```text
f_0=A h,                    g_0=B h^k,       Sigma|h.  (23)
```

The two upper singleton rows are

```text
q f_1'g_2=p f_1g_2',
q f_2'g_2=4k f_2g_2'.                                  (24)
```

Let `L_2,L_3` denote the coefficients of the next two double rows:

```text
L_2=W_(p,1)(f_1,g_1)+W_(4k,-2k)(f_2,g_0),
L_3=W_(-2,q)(f_0,g_2)+W_(4k,1)(f_2,g_1).               (25)
```

Using `(24)` to eliminate `g_1` from `L_2=L_3=0` gives

```text
(q h'g_2+2h g_2')
[A p q f_1g_2+16k^3 B f_2^2 h^(k-1)]=0.               (26)
```

The first factor cannot vanish, since it would make
`(h^q g_2^2)'=0`.  If `k>1`, logarithmically differentiating the second
factor and using `(24)` gives

```text
(k-1)(h'/h+2g_2'/(qg_2))=0,                            (27)
```

which makes the first factor vanish after all.  Thus `k>1` is impossible.

For `k=1`, unique factorization in `(24)` gives

```text
f_1=L K,             f_2=U K^4,             g_2=M K^7. (28)
```

The first bridge integrates to

```text
g_1=K[C_0-(4BU/L)hK^2].                                (29)
```

The other bridge forces

```text
7ALM+16BU^2=0,                                         (30)
```

and the scalar row factors exactly as

```text
(K h'+2hK')
[AC_0-BL-(12ABU/L)hK^2]=1.                             (31)
```

This contradicts `(10)`.

### 5.2 First ladder, mirror orientation

The low and upper singleton solutions are

```text
f_0=A h^k,        g_0=B h,
f_1=L K,          f_2=U K^(2k+2),        g_2=M K^(6k+1). (32)
```

Substitute the first middle row into the second.  With `q=6k+1`, the result is

```text
(K h'+2hK')K^(4k+2)
[kqAM(hK^2)^(k-1)+4(k+1)^2U^2B/L]=0.                  (33)
```

For `k>1`, the bracket is a nonzero constant plus a nonzero positive-degree
power of `hK^2`; it cannot vanish.  At `k=1`, `(33)` is exactly `(30)`, and
the exceptional calculation `(28)--(31)` closes it.  The first hooked ladder
is empty in both orientations.

### 5.3 Central hooked ladder

The singleton rows give

```text
f_0=A h^k,             g_0=B h^(k+1),
f_1=L K,               f_2=U K^(3k+2),
g_2=M K^(6k),          Sigma|h.                        (34)
```

The lower bridge has the particular solution

```text
g_1=lambda hK,              lambda=(k+1)LB/(kA).       (35)
```

The difference of any two solutions, say `w`, obeys

```text
3h w'-2h'w=0,
(w^3/h^2)'=w^2(3h w'-2h'w)/h^3=0.                     (36)
```

If a nonzero polynomial homogeneous mode exists, arm integrality requires
`3 ord(w)=2 ord(h)`, hence `ord(h)` is a multiple of three and
`ord(w)>=2`.  The particular term also then has arm order at least three.
But `g_1` is the unique weight-`-2` scalar channel and must be simple by
`(8)`.  Therefore the homogeneous mode is absent and `(35)` is exhaustive.

The scalar row is

```text
-(K h'+3hK')
[lambda L K+(k+1)(3k+2)UB h^k K^(3k+1)]=1,            (37)
```

contradicting `(10)` with `a=3`.

### 5.4 Upper hooked ladder

Put `p=2k-1`.  The endpoint equations give

```text
f_0=A h,              g_0=B h^(3k+1),
f_1=L K^p,            f_2=U K^(4k),           g_2=M K,
Sigma|h.                                                 (38)
```

Let

```text
R_2=W_(p,-4k-1)(f_1,g_1)+W_(4k,-6k-2)(f_2,g_0).       (39)
```

If `S` denotes the scalar row, direct elimination gives the polynomial
identity

```text
S-[4kU/(pL)]K^(2k+1)R_2
=(K h'+2hK')
 [AM+16k^2(3k+1)U^2B/(pL) h^(3k)K^(6k)].              (40)
```

A Darboux pair has `R_2=0` and `S=1`, contradicting `(10)`.  This closes all
four ladders in the displayed orientation `H`.

Two nondeletion checks are worth recording.  In the primary ladder, the two
summands in each zero row must cancel nontrivially:

```text
ell_2=W_(p,1)(f_1,g_1)=-W_(4k,-2k)(f_2,g_0) !=0;
otherwise the surviving 01 scalar rectangle reduces to 2 x 2,

ell_3=W_(-2,q)(f_0,g_2)=-W_(4k,1)(f_2,g_1) !=0;
otherwise deleting g_2 reduces to 3 x 2.                    (41)
```

Equations `(25)--(26)`, rather than either deletion, are therefore the minimal
genuinely connected leak obstruction.  This is also why a `2 x 4` theorem
cannot simply replace the hooked calculation: deleting one vertex here gives
`3 x 2`, while the full cell has three pieces on both sides.

### 5.5 The simultaneously reversed hook

The additive classification also has the reversal-sensitive representative

```text
A={0,d,2d},                    B={0,2d,3d},
00; 10; 01+20; 02+11; 12+21; 22.                      (41a)
```

It cannot be discarded by simultaneous reversal, because `(3)` distinguishes
positive from negative absolute weights.  Testing the twelve scalar-arm
placements against the singleton rows `00,10,22` leaves, for `d>=3`, exactly
the address `21` with weights `(1,-2)`.  The only same-sign cases not killed
immediately by a mixed or zero/nonzero singleton are:

```text
20:(1,-2):  ord(g_0)=1 forces ord(f_0)=(2d-1)/2;
11:(-2,1): ord(f_1)=1 forces both (2d-1)/2 and (d+2)/2 integral;
12:(-2,1): ord(f_1)=1 forces both (3d-1)/2 and (d+2)/2 integral.
```

The first two are impossible because `2d-1` is odd, and the last asks
opposite parities of `d`.  All remaining placements have a strict mixed-sign
singleton (with the direct zero/nonzero boundary at `d=1,2`).

For the survivor, the absolute supports are

```text
supp(P)=(1-2d,1-d,1),          supp(Q)=(-2d-2,-2,d-2).
```

The singleton equations give, for nonzero constants `A,B,L,U,M`,

```text
f_0=A h^(2d-1),   f_1=L h^(d-1),   g_0=B h^(2d+2),
f_2=U K,          g_2=M K^(d-2),              Sigma|h. (41b)
```

At an arm write `ell=ord(h)>=1`.  Scalar effectiveness at `21` gives
`ord(g_1)=1` and `ord(f_2)=0`.  In the lower double row

```text
W_(1-2d,-2)(f_0,g_1)+W_(1,-2d-2)(f_2,g_0)=0,          (41c)
```

the two summands have respective orders

```text
(2d-1)ell,                    (2d+2)ell-1,
```

and nonzero initial multipliers `(2d-1)(1-2ell)` and
`-(2d+2)ell`.  Cancellation would require `3ell=1`, impossible.  When
`d=3`, the other address `12` also has weight pair `(-2,1)`, but `(41b)`
gives `ord(f_1)=2ell>=2`, so it cannot supply a simple arm or alternate with
`21`.  Hence `bar H` is empty as well.

## 6. The reflected triple

Normalize its gaps as

```text
(x,y ; y,x),                    0<x<y.                 (42)
```

The only nonsingleton fibre is `02+11+20`, so it must be scalar.  The six
arm gates in that fibre all fail:

| gated address | gate | singleton obstruction |
|---|---|---|
| `02` | `(-2,1)` | `10=(x-2,1-x-y)`; `x>=3` is mixed, `x=2` is zero/nonzero, and `x=1` leaves `01=(-2,0)` |
| `02` | `(1,-2)` | `00=(1,-2-x-y)` is mixed |
| `11` | `(-2,1)` | `01=(-2-x,1)` is mixed |
| `11` | `(1,-2)` | `10=(1,-2-y)` is mixed |
| `20` | `(-2,1)` | `00=(-2-x-y,1)` is mixed |
| `20` | `(1,-2)` | `01=(1-x-y,y-2)` is mixed for `y>2` and zero/nonzero for `y=2` |

Thus the reflected family is empty.

## 7. The four oriented Euclidean families

For `E+`, the double `01+10` is a `2 x 2` diagonal whose other corners
`00,11` are singleton, so it cannot be scalar.  The other double is
`12+20`.  Its four gates die respectively on

```text
12:(-2,1)  -> 02 mixed,
12:(1,-2)  -> 22 mixed,
20:(-2,1)  -> 00 mixed,
20:(1,-2)  -> 02=(1-p-q,q-2).                          (43)
```

In the last row `q>2`; equality `q=2,p=1` is exactly the excluded hooked
boundary `q=2p`.

For `E-`, the double `12+20` is a `2 x 2` diagonal with singleton corners
`10,22`, so it cannot be scalar.  The other double is `02+11`.  Its gates are

```text
02:(-2,1)  -> 10=(p-2,1-q); p>=3 mixed, p=2 zero/nonzero,
                 p=1 leaves 01=(-2,0),
02:(1,-2)  -> 00 mixed,
11:(-2,1)  -> 01 mixed,
11:(1,-2)  -> 10 mixed.                                (44)
```

For the simultaneously reversed representative `bar E+`, the doubles are
`02+10` and `12+21`.  The latter is a `2 x 2` diagonal with singleton corners
`11,22`.  In the former, the four gates

```text
02:(-2,1), 02:(1,-2), 10:(-2,1), 10:(1,-2)
```

die respectively on the singleton cells `20,00,00,20` by strict mixed signs.

For `bar E-`, the double `02+10` is the rectangle with singleton corners
`00,12`.  The four gates in `11+20` die as follows:

```text
11:(-2,1) -> 01 mixed,       11:(1,-2) -> 21 mixed,
20:(-2,1) -> 00 mixed,
20:(1,-2) -> 01 mixed for p>=3, zero/nonzero for p=2,
             and 21 mixed for p=1.                    (44a)
```

Thus all four oriented Euclidean representatives are empty.  Sections 3--7
exhaust every possible three-by-three collision pattern without quotienting
by a false regularity symmetry, proving the theorem.

## 8. Exact support floor and scope

THM-3569 excludes one retained piece against arbitrary width and all cells
with one side at most two and the other at most three.  THM-3583 strengthens
the two-piece sector to

```text
#supp(P)=2  implies  #supp(Q)>=5,                       (45)
```

and the present theorem excludes `3+3`.  Therefore every polynomial Darboux
pair in the stated universal squarefree exponent-two scope has

```text
#supp(P)+#supp(Q)>=7.                                   (46)
```

At total seven the only unexcluded partitions are

```text
(2,5), (3,4), (4,3), (5,2).                            (47)
```

This is a nonconstant homogeneous-support floor, not an existence theorem.
The cells `2 x 5`, `3 x 4`, and larger remain open, as do existence of any
polynomial Darboux pair on `Y_Sigma` and the planar Jacobian conjecture.

The boundaries are sharp and intentionally outside the statement:

1. On the Laurent chart,

   ```text
   {b,-1/c}=1,                                          (48)
   ```

   so retaining one pole destroys polynomiality and defeats the support gate.

2. If `Sigma=ub+v` has degree one, then `A_Sigma=C[c,e]` and

   ```text
   {c,-e/u}=1.                                          (49)
   ```

3. If `Sigma` has a repeated root `beta`, the Poisson tensor vanishes on
   `c=0,b=beta`; evaluating there makes every polynomial bracket zero.  Thus
   the singular repeated-root target has a stronger no-go for a different,
   Poisson-degenerate reason.

4. Positive characteristic is excluded because it can annihilate the
   derivative, arm, coprimality, and leading-degree multipliers.

## 9. Optimization-safe exact verification

Reproduce with

```bash
python3 04-computation/jc2_universal_exponent_two_three_by_three_nonentry_thm3592.py
python3 -O 04-computation/jc2_universal_exponent_two_three_by_three_nonentry_thm3592.py
```

The companion uses only the Python standard library and unconditional
`require` gates.  It checks:

- the identity `(14)`, triple criterion, and complete normalized collision
  catalogue on every positive gap quadruple in `[1,24]^4`, with primitive
  dilation normalization;
- every component/deletion verdict in that census and the exact fibre words;
- all diagonal and all four oriented Euclidean `2 x 2` rectangles in the
  hostile ranges;
- all twelve hook gate placements for `1<=d<=64`, including the `H2` bridge
  valuation and every zero-weight boundary;
- all twelve reversed-hook gates for `1<=d<=64`, its unique candidate ladder,
  and 3,968 exact arm-valuation mismatches;
- all reflected-triple and four-orientation Euclidean arm gates through gap
  24;
- every displayed hooked differential identity for `1<=k<=8`, including the
  `k=1` compatibility, cube first integral, and H3 row subtraction;
- the Euler leading-degree gate and all three sharp hostiles; and
- the exact support-floor partitions.

Those finite ranges are exact identity and hostile controls.  The all-support,
all-degree conclusion is the Schur-multiset classification, component lemma,
arm valuation, UFD, logarithmic-derivative, and Euler-factor proof above; it is
not extrapolated from the finite ranges.

**QED.**
