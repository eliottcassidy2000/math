---
id: THM-4058
title: "Exceptional affine triangle period and simple-zero monomial obstruction ladder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with the strict
  formal/local and fixed-first-output boundary stated below.  For the
  exceptional quartic of THM-4054, the affine three-branch fixed-a cokernel
  is generated in every retained degree by one normalized triangular
  period.  Its response on w^r first occurs in degree r+4 and is
  -(r/2)rho.  Consequently, for every m>=2 and every nonzero scalar gamma,
  the fixed-a constant packet for H(t)=t+gamma t^m exists through cutoff
  m+2 and fails first at cutoff m+3, with response
  gamma*binom(m,2)*rho.  This does not obstruct a mixed pair, construct a
  global polynomial pair, or imply anything about JC(2).
source: jc2-double-zero-rebuild-20260824 / fixed-a period lane, 2026-08-24
audit: >
  PASS.  The production companion reconstructs the three branch graphs and
  pairwise intersections exactly over K=Q(alpha), proves that their common
  intersection jet first opens at A^5, checks delta=(26/15)rho, the exact
  edge vectors, and the signed triangle area -(15/52)delta^2.  A separate
  implementation over the good split prime 137 reconstructs all four roots,
  the triangle, and complete retained fixed-a matrices.  For every root and
  every m=2,...,12 it independently finds zero response at cutoff m+2 and
  the sole response binom(m,2)*gamma*rho at cutoff m+3, with rank
  rows-(N+1).  Normal and optimized transcripts agree byte for byte; both
  companions have zero Python assert nodes and zero float literals.  The
  split-prime calculation is an independent hostile audit, not the
  characteristic-zero proof.
depends_on:
  - THM-4054-exceptional-affine-simple-zero-retained-packet
related:
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
script: 04-computation/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058.py
independent_script: 04-computation/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058_independent_audit.py
output: 05-knowledge/results/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058.out
independent_output: 05-knowledge/results/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058_independent_audit.out
script_sha256: 99d6e27e93129636faba409750fbd1b71fd3f60a7b65117ea24f906cd78888e1
independent_script_sha256: 9cf33bbd5d022893a2dc31bec837ad708d1d9ef501015fb38da54feaea7a94d2
output_sha256: 8bf827144b5da0e14dab9fc43eae59f2ec695d6bf60e4298e6b1d28b0369e26b
independent_output_sha256: 12037bf83bf0754677a6f80d1b89830e2e3d5cd9acb66f3201a50e623fbce939
hash_basis: raw LF bytes
---

# THM-4058 -- the exceptional fixed-`a` obstruction is a triangular period

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with a strict
formal/local boundary.**  THM-4054 found the first fixed-`a` failure for the
single displacement `H=t+t^2`.  The mechanism is uniform: the three affine
branch surfaces form a formal triangle whose normalized period is the
complete retained cokernel.  The triangle first opens in degree five, so
stable powers form a sharp one-step obstruction ladder.  A simple-zero
source reparametrization then turns that ladder into an all-`m` theorem.

All characteristic-zero calculations below are over the exceptional field

```text
K=Q[alpha]/(72783360 alpha^4-77822208 alpha^3-28419741 alpha^2
                                      +7849770 alpha-1276420).       (1)
```

The notation `Q_alpha`, `a=q/D^2`, `c=xD(D+2)`, `A=a+3/4`, and `rho` is as
in THM-4054.  A retained cutoff is total maximal-ideal degree on each of the
three source germs.

## 1. A parity-free normalization coordinate

Put

```text
s=xD.                                                        (2)
```

The identities of THM-4054 give the exact inverse-coordinate formulas

```text
D=1+a s^2,             x=s/(1+a s^2),
q=a(1+a s^2)^2,        c=3s+a s^3.                         (3)
```

At `a=-3/4` and `c=0`, the three solutions `s=2,0,-2` correspond in this
order to the retained source points `x=-1,0,1`.  Define

```text
R(a,s)=a(1+a s^2)^2-Q_alpha(s/(1+a s^2)).                  (4)
```

For `i=0,1,2`, let `s_i(A,c)` be the unique formal solution of the last
equation in `(3)` near `2,0,-2`, and put

```text
R_i(A,c)=R(-3/4+A,s_i(A,c)).                               (5)
```

On the affine source `q=Q_alpha(x)+w`, the three image surfaces are the
graphs `w=R_i(A,c)`.  Exact differentiation gives

```text
R_i=A+lambda_i c+O((A,c)^2),
(lambda_0,lambda_1,lambda_2)=(3/4,-1/3,-3/4).             (6)
```

The three tangent planes are distinct.  No evenness of `Q_alpha` is used.

Let `G_i(A,c)=G(A,c,R_i(A,c))`.  Since
`Jac_(x,w)(a,c)=-3` on the affine source, inversion of the `(a,c)` Jacobian
gives `a_x=-3R_(i,c)`.  Therefore the fixed-`a` density operator is

```text
L_i G=-3G_c+a_xG_w
     =-3(G_c+R_(i,c)G_w)
     =-3 dG_i/dc.                                         (7)
```

Thus the fixed-`a` problem is simultaneous exact differentiation along the
three graph branches.

## 2. The exceptional intersection triangle opens in degree five

For `ij=01,12,20`, let

```text
v_ij(A)=(c_ij(A),w_ij(A))                                 (8)
```

be the unique pairwise intersection of graphs `i` and `j`.  Exact implicit
series arithmetic shows that all three vertices share one jet through
degree four:

```text
c_ij=c_*(A)+C_ij A^5+O(A^6),
w_ij=w_*(A)+W_ij A^5+O(A^6),                              (9)
```

where

```text
c_*=
 (-64/27+64alpha/3)A^2
 +(-74752/6561+87296alpha/729+8192alpha^2/81)A^3
 +(-36628480/1594323+9671680alpha/59049
      +15073280alpha^2/6561-4014080alpha^3/2187)A^4,      (10)

w_*=A+(64/81-64alpha/9)A^2
 +(74752/19683-87296alpha/2187-8192alpha^2/243)A^3
 +(-576256/19683+88960alpha/6561
      -709568alpha^2/729-81920alpha^3/81)A^4.             (11)
```

Orient branch `0` from `v_20` to `v_01`, branch `1` from `v_01` to `v_12`,
and branch `2` from `v_12` to `v_20`.  Recall the THM-4054 response

```text
rho=
 -2073506706944/1678822119
 +(372679949312/62178597)alpha
 -(184159683584/6908733)alpha^2
 -(73442787328/2302911)alpha^3.                           (12)
```

Set `delta=[A^5](c_20-c_12)`.  Exact reduction in `(1)` gives

```text
delta=(26/15)rho
=-4147013413888/1937102445
 +(745359898624/71744535)alpha
 -(368319367168/7971615)alpha^2
 -(146885574656/2657205)alpha^3.                          (13)
```

The oriented leading edge vectors are

```text
Delta c=delta(5/13,-18/13,1),
Delta w=delta(15/52,6/13,-3/4)
       =(lambda_0 Delta c_0,lambda_1 Delta c_1,
                                      lambda_2 Delta c_2). (14)
```

Both coordinate sums vanish.  Since the norm of `rho` in THM-4054 is
nonzero, `delta` is nonzero under every complex embedding of `K`.

## 3. The normalized triangular period is the complete cokernel

For an arbitrary triple `f=(f_0,f_1,f_2)` in `K[[A,c]]^3`, define

```text
Pi(f)= integral_(c_20)^(c_01) f_0 dc
      +integral_(c_01)^(c_12) f_1 dc
      +integral_(c_12)^(c_20) f_2 dc,

Lambda(f)=Pi(f)/(delta A^5).                              (15)
```

Every endpoint difference is divisible by `A^5`, so `Lambda(f)` is a
well-defined formal series in `A`.  For every target germ `G`, `(7)` gives

```text
Pi(dG/dc)=
 [G(v_01)-G(v_20)]+[G(v_12)-G(v_01)]
                         +[G(v_20)-G(v_12)]=0.            (16)
```

This proves annihilation.  The following graded count proves completeness.
In total degree `d=n+1`, the restrictions of homogeneous target functions
to the tangent planes `(6)` have dimension

```text
binom(d+2,2)-binom(d-1,2)=3d.                             (17)
```

Here the second term is zero for `d<3`; for `d>=3`, the common vanishing
ideal is generated by the product of the three distinct plane equations.
The kernel of tangential `c`-differentiation is exactly the common line
`K A^d`: a branchwise constant tuple `(h_i A^d)` is a target restriction
only when all three `h_i` agree on `c=0,w=A`.  Hence the degree-`n` rank is

```text
3d-1=3n+2                                                (18)
```

in the `3(n+1)=3n+3` dimensional branch-triple space.  The graded cokernel
is therefore exactly one-dimensional in every degree.

If the first homogeneous part of `f` has

```text
f_i(A,0)=q_i A^n+O(A^(n+1)),                             (19)
```

then `(9),(14),(15)` give

```text
Lambda(f)=mu(q)A^n+O(A^(n+1)),
mu(q)=5q_0/13-18q_1/13+q_2.                              (20)
```

The identities

```text
sum mu_i=0,                 sum mu_i lambda_i=0           (21)
```

show directly that `mu` spans the one-dimensional graded cokernel in
`(18)`.  Consequently the coefficients

```text
Lambda_0,...,Lambda_N,       Lambda_j(f)=[A^j]Lambda(f), (22)
```

are `N+1` independent relations on the cutoff-`N` fixed-`a` map.  The
filtered matrix has the graded diagonal ranks `(18)`, so its cokernel has
dimension at most `N+1`.  Therefore `(22)` is the complete cokernel:

```text
f lies in the cutoff-N fixed-a image
       iff Lambda_0(f)=...=Lambda_N(f)=0.                 (23)
```

On each affine branch, `(x-p,w)` and `(A,c)` are formal coordinate systems
because their Jacobian is `-3`.  They preserve maximal-ideal powers.  Thus
`(23)` applies to exactly the retained cutoff and complete degree-at-most
`N+1` target-jet domain used in THM-4054, not to a different filtration.

## 4. The all-power affine obstruction ladder

Translate `v_20` to the origin and divide the leading triangle by
`delta A^5`.  Equations `(6),(14)` give its vertices

```text
(0,0),              (5/13,15/52),
(-1,3/4),           (0,0).                              (24)
```

The exact oriented line integral is

```text
Pi(w)=-(15/52)delta^2 A^10+O(A^11).                      (25)
```

The common vertex value satisfies `w_*=A+O(A^2)`, whereas both coordinates
along every small edge differ from the common jet by `O(A^5)`.  Expanding
around `w_*`, the constant term integrates to zero, the linear term gives
`r w_*^(r-1)Pi(w)`, and every term of binomial degree at least two has
strictly larger valuation.  Hence, for every integer `r>=1`,

```text
Lambda(w^r)=-(r/2)rho A^(r+4)+O(A^(r+5)).                (26)
```

Completeness `(23)` and nonvanishing of `rho` now prove the sharp ladder:

```text
w^r is in the fixed-a retained image through cutoff r+3,
w^r is not in that image at cutoff r+4.                  (27)
```

This is an all-`r` characteristic-zero statement, not an extrapolation from
the finite audit range.

## 5. Every nonlinear monomial simple zero hits the ladder

Let

```text
H(t)=t+gamma t^m,             m>=2,        gamma!=0,     (28)
```

and let `u=H(t)` with formal inverse `t=h(u)`.  Given a target germ
`G(A,c,w)`, set

```text
G_tilde(A,c,u)=G(A,c,h(u)).                              (29)
```

The nonlinear fixed-`a` operator, expressed on the affine graph `u=R_i`, is

```text
L_H G=-3(H'(t)G_c+R_(i,c)G_w)
     =-3(dG_tilde/dc)/h'(u).                             (30)
```

Thus `L_H G=f` is equivalent, under a filtered formal target automorphism,
to the affine equation

```text
-3 dG_tilde/dc=h'(u)f.                                  (31)
```

In particular, the nonlinear constant-density relation is
`Lambda(h'(u))`.  Formal inversion gives

```text
h'(u)=1-m gamma u^(m-1)+O(u^(2m-2)).                    (32)
```

Now `Lambda(1)=0` exactly.  Apply `(26)` with `r=m-1`.  Every term hidden in
`(32)` has stable exponent at least `2m-2` and, by `(26)`, cannot contribute
before degree `2m+2>m+3`.  The first response is therefore

```text
Lambda(h')=
 gamma binom(m,2)rho A^(m+3)+O(A^(m+4)).                (33)
```

Combining `(23)` and `(33)` proves, for every `m>=2` and every nonzero
`gamma`,

```text
the fixed-a constant packet exists through cutoff m+2,
it fails first at cutoff m+3,
its first normalized response is gamma binom(m,2)rho.    (34)
```

The result holds after every complex embedding of `(1)`.  It also holds for
arbitrary nonzero complex `gamma` after scalar extension.  At `gamma=0`,
the affine pair `G=-c/3` gives unit density exactly, so the nonlinearity in
`(28)` is load-bearing.

## 6. Exact and hostile replay

The production companion uses an exact four-coordinate implementation of
`K`.  It Newton-solves `(3)` and the three pairwise intersections through
degree six, verifies `(6),(9)--(14),(25)`, and independently composes the
formal inverses for `m=2,...,12`.  It uses no retained matrix.

The independent companion uses a separate finite-field implementation.  At
the good split prime

```text
p=137,        alpha=44,82,92,134,                         (35)
```

it reconstructs the triangle and obtains

```text
rho mod 137=(8,85,12,135).                               (36)
```

It then builds the complete fixed-`a` retained matrix at every cutoff
`4,...,15`.  For all four roots and every `m=2,...,12`, it finds

```text
rank at cutoff N = rows-(N+1),
response at N=m+2 =0,
sole response at N=m+3 =gamma binom(m,2)rho.              (37)
```

The serialized audit-table hash is

```text
1641250c558fbdeea9f2e2ccff1a515761d4af596f95ae5bce4587d12edc72c9. (38)
```

This modular path is an independent hostile control.  The proof of the
all-`m` characteristic-zero statement is Sections 1--5.

## 7. Connection contract and strict boundary

| field | contract |
|---|---|
| source | complete target function jets `G(A,c,w)` with first output frozen to `a`; after reparametrization, the affine right side `h'(u)` |
| target | three retained source-density germs through total cutoff `N` |
| map | fixed-`a` Jacobian density, changed to simultaneous branch differentiation by `(7)` and `(30)` |
| preserved predicate | finite-jet membership, first failure cutoff, and normalized response under the filtered simple-zero reparametrization |
| destroyed information | coherence of chosen primitives across all cutoffs, convergence, global regularity of `a`, freedom to change the first output, behavior at infinity, injectivity, and properness |
| needed sidecar | a mixed-pair nonlinear obstruction or construction, followed by global denominator and algebraization control |
| cheapest decisive tests | `delta=(26/15)rho`; edge closure and signed area; the graded rank `(18)`; direct retained ranks and response polynomials in `(37)` |

The theorem obstructs only the gauge `F=a`.  THM-4054 already proves that
the full mixed tangent bank absorbs the quadratic perturbation through its
first fixed-`a` failure.  Nothing here contradicts that escape or promotes
the fixed-coordinate period to an obstruction for arbitrary pairs.  No
global polynomial pair, Keller map, counterexample, or consequence for
`JC(2)` is constructed.

## 8. Reproduction

```bash
python3 04-computation/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058.py
python3 -O 04-computation/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058.py
python3 04-computation/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058_independent_audit.py
python3 -O 04-computation/jc2_exceptional_affine_triangle_period_monomial_ladder_thm4058_independent_audit.py
```

Both normal/optimized pairs are byte-identical to their stored transcripts.
The production and independent script/output hashes are pinned in the
frontmatter.
