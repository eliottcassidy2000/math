---
id: THM-2238
title: "Terminal-cover moment, hinge, Gram, and fibrewise convex hierarchy"
status: >
  PROVED + HOSTILE-AUDITED. For every finite residual target and finite
  terminal-mask family, selection of p distinct masks has a monotone
  fixed-normalization factorial-moment lower hierarchy which is exact at
  degree p. Its degree-one outer selection problem is exactly the Ky--Fan
  hinge dual used by THM-2216, after the singleton-meet loss and optimization
  over the threshold. Degree two is an exact binary quadratic program over
  the terminal-mask intersection Gram matrix. Unaggregated fibre capacities
  give a stronger mixed-integer hinge and a forced convex pair-overlap
  correction psi_H(S)=sum_(ell>=1)(S-ell H)_+. Continuous hypersimplex
  relaxations are rigorous lower bounds but need not be exact; explicit
  finite witnesses have strict integrality and prove that the global Gram
  and fibrewise count programs are numerically incomparable. This is an
  abstract finite cover theorem, not a uniform LRC estimate or a proof of
  LRC(14).
source: codex-2026-07-24-terminal-cover-hierarchy
depends_on:
  - THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy
  - THM-2236-pointwise-nested-binomial-minorants-and-cubic-vertex-fan
related:
  - THM-2209-sharp-quadratic-reversed-peel-and-joint-fourier-ledger
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
---

# THM-2238 -- the terminal-cover hierarchy

This theorem separates three optimizations which share similar syntax but
retain different information:

```text
mask selection,
factorial-moment minorants,
and blocker-capacity hinge envelopes.
```

Their degree-one boundary agrees exactly. Beyond degree one, the required
Gram matrices and incidence sidecars are different.

## 1. Finite terminal-cover packet

Let `R` be a finite residual target set, let `Qset` be a finite terminal-mask
index set, and let

```text
N_q subset R                         for q in Qset.
```

Fix

```text
1<=p<=|Qset|.
```

For every `p`-element set `T subset Qset`, put

```text
K_T(y)=sum_(q in T)1_(N_q)(y),
I_s(T)=sum_(y in R) binom(K_T(y),s),
U(T)=|R\union_(q in T)N_q|,
W=|R|.                                                   (1)
```

Thus

```text
U(T)=sum_(y in R)1_(K_T(y)=0),             I_0(T)=W,
```

and `I_s(T)` is the sum of all `s`-fold intersection sizes among the
selected masks.

For `1<=r<=p`, let `C_(r,p)` be THM-2236's complete coefficient polyhedron

```text
1-K+sum_(s=2)^r a_s binom(K,s)<=1_(K=0)
                                      for K=0,...,p.     (2)
```

Define the fixed-normalization adaptive value

```text
B_r(T)
 =W-I_1(T)+max_(a in C_(r,p))
      sum_(s=2)^r a_s I_s(T),                           (3)

Phi_r=min_(T subset Qset, |T|=p) B_r(T).               (4)
```

The maximum exists whenever it is finite. Indeed, `C_(r,p)` is a nonempty
closed polyhedron and a linear functional bounded above on a nonempty
polyhedron attains its optimum. It is finite on every realizable packet
because (2) makes every objective at most `U(T)`.

## 2. Monotonicity, exactness, and the unrestricted LP

The embedding

```text
C_(r,p) -> C_(r+1,p),             a_(r+1)=0
```

shows

```text
B_1(T)<=...<=B_p(T)<=U(T).                            (5)
```

At `r=p`, the exact inclusion--exclusion polynomial belongs to
`C_(p,p)`, so

```text
B_p(T)=U(T),                 Phi_p=min_(|T|=p)U(T).   (6)
```

Taking minima in (5) preserves the order:

```text
Phi_1<=...<=Phi_p.                                    (7)
```

Consequently `Phi_p>0` if and only if no set of at most `p` masks covers
`R`: if fewer than `p` masks cover, adjoining arbitrary distinct masks
gives a `p`-mask cover.

THM-2236's canonical pointwise-nested polynomial gives the smaller
structured value

```text
J_r(T)=sum_(y in R)P_(r,p)(K_T(y))<=B_r(T).           (8)
```

Let `Lambda_r(T)` be THM-2210's unrestricted adaptive factorial-moment LP
value for the same packet. Its dual contains both `C_(r,p)` and the zero
polynomial, hence

```text
Lambda_r(T)>=max(0,B_r(T)).                           (9)
```

After the same outer mask selection,

```text
min_T Lambda_r(T)>=max(0,Phi_r).                     (10)
```

The min and max in (3)--(4) may not be interchanged. For `r>=3`, the
normalized maximizing polynomial can depend on `T`; THM-2236's hostile
cubic already changes vertex with the moment packet.

## 3. Degree one is the Ky--Fan hinge

Put

```text
c_q=|N_q|.
```

Since `B_1(T)=W-sum_(q in T)c_q`,

```text
Phi_1=W-Top_p(c).                                     (11)
```

The exact integer Ky--Fan identity gives

```text
Top_p(c)
 =min_(theta in Z_{\ge0})
    [p theta+sum_(q in Qset)(c_q-theta)_+].
```

Therefore

```text
Phi_1
 =max_(theta in Z_{\ge0})
    [W-p theta-sum_q(c_q-theta)_+].                  (12)
```

This is the precise selection dual behind THM-2216. In its two-blocker
application the actual capacity row is first replaced by the singleton meet
envelope. Optimizing THM-2216's exact hinge over `theta` then gives (12) for
that envelope. A certificate at one prescribed threshold is a further
sufficient relaxation, not another equality statement.

## 4. Degree two is the terminal-mask Gram program

Assume `p>=2`. Let

```text
A[y,q]=1_(y in N_q),
G=A^T A,
c=diag(G).                                            (13)
```

Thus `G` is positive semidefinite and records every pairwise terminal-mask
intersection. If `z in {0,1}^Qset` selects `p` masks, then

```text
I_1=c^T z,
I_2=(z^T G z-c^T z)/2.                               (14)
```

THM-2209's sharp coefficient `2/p` yields

```text
Phi_2
 =min_(z in {0,1}^Qset, 1^T z=p)
    [W+(1/p)z^T Gz-(1+1/p)c^T z].                   (15)
```

This `mask x mask` Gram matrix is not THM-2216's
`blocker x blocker` superlevel Gram. The former retains intersections among
chosen terminal masks; the latter bounds a meet of singleton blocker rows.

Relaxing `z` to the hypersimplex

```text
0<=z<=1,                  1^T z=p
```

enlarges a minimization domain. Its optimum is therefore a rigorous lower
bound for `Phi_2`, so a positive relaxed optimum proves noncoverage. The
relaxation can lose positivity.

### Strict integrality-gap witness

Let `Qset={1,2,3,4}`, let `p=2`, and let the six targets `y_E` be indexed by
the two-subsets `E` of `Qset`. Define

```text
N_q={y_E:q notin E}.                                  (16)
```

Every mask has size three. Every pair `{i,j}` misses exactly `y_{i,j}`, so

```text
Phi_1=0,                     Phi_2=1.
```

But `z=(1/2,1/2,1/2,1/2)` gives load one on every target and value zero in
(15). Thus PSD and continuous optimization do not remove the integral
selection sidecar.

The same capacity vector can also arise from a coverable set family, so no
criterion using only the `c_q` can recover the quadratic improvement.

## 5. Fibrewise hinge before aggregation

Now partition

```text
R=disjoint_union_(x in X)Y_x,
H_x=|Y_x|>=1,
h_q(x)=|N_q intersection Y_x|.                      (17)
```

For a selected set `T`, put

```text
S_x(T)=sum_(q in T)h_q(x).
```

The union in one fibre has size at most `S_x(T)`. Hence

```text
U(T)>=sum_x (H_x-S_x(T))_+.                          (18)
```

After minimization over one common `p`-set of labels,

```text
Psi_1
 =min_(|T|=p)sum_x(H_x-S_x(T))_+                    (19)
```

is a mixed-integer hinge certificate. It dominates the aggregated
degree-one gap:

```text
Psi_1>=max(0,W-Top_p(c)),                            (20)
```

because the sum of positive parts dominates the positive part of the sum.
This retains the common-winner incompatibility lost by optimizing labels
independently in each fibre.

The hypersimplex relaxation of (19) remains a safe lower bound on the
integer optimum. Positivity of that relaxation still proves noncoverage.

## 6. Forced pair overlap from fibre counts

Assume `p>=2`. For integers

```text
H>=1,             0<=S<=pH,
S=aH+b,           0<=b<H,
```

define

```text
psi_H(S)
 =H binom(a,2)+ab
 =sum_(ell>=1)(S-ell H)_+.                           (21)
```

If integers `k_1,...,k_H` lie in `{0,...,p}` and sum to `S`, discrete
convexity of `binom(k,2)` shows

```text
sum_(i=1)^H binom(k_i,2)>=psi_H(S).                  (22)
```

Indeed, if two entries differ by at least two, moving one unit from the
larger to the smaller strictly decreases the sum. The minimum is therefore
the balanced vector with `b` entries `a+1` and `H-b` entries `a`, giving
(21).

Apply (22) to the multiplicities `K_T(y)` in each fibre. Then

```text
I_(2,x)(T)>=psi_(H_x)(S_x(T)).
```

The pointwise quadratic minorant and nonnegativity of uncovered mass give

```text
U(T)
 >=sum_x max(
      0,
      H_x-S_x(T)+(2/p)psi_(H_x)(S_x(T))
    ).                                               (23)
```

This count-only program is a coarsening of the full fibre incidence tensor,
but it is stronger than global capacities. It is not generally ordered
against the un-clipped global Gram value (15). A fibre-resolved, separately
clipped exact Gram program retains at least as much information as either
coarsening.

For a nondegenerate boundary witness take four masks, `p=3`, and five
one-point fibres with incidence rows

```text
e_1,e_2,e_3,e_4,e_1+e_2.
```

Every integral three-mask choice misses one singleton, so the fibrewise
count program has minimum one. The global quadratic value has minimum
`2/3`, attained by a triple containing masks one and two. Thus the
fibrewise count program can be strictly larger.

The reverse strict inequality already occurs with `p=2`, one fibre of
size two, and two masks which both contain the same one target. The only
selection has `S=2` and `psi_2(2)=0`, so (23) gives zero. Its multiplicity
row is `(2,0)`, however, and the exact global quadratic value (15) is one.
Thus the two rigorous quadratic lower programs are numerically
incomparable when their clipping and aggregation differ.

## 7. Consequence and boundary

The hierarchy identifies an information **fork**, not a linear ladder:

```text
                         /-> per-fibre capacity rows -> convex count bounds
global capacities c ---<
                         \-> global mask Gram G -> global moment tensors

{per-fibre rows, global G}
  + finer incidence data
  -> fibre-resolved Gram family (G_x)_x
  + higher incidence data
  -> higher fibre-resolved incidence tensors.

Global factorial moments through degree p -> the exact Boolean union.
```

The two branches are incomparable refinements of `c`: per-fibre rows retain
where singleton capacity is spent, while `G` retains pair intersections
after forgetting their fibres. The paired sidecar `{h,G}` retains both
kinds only coarsely; the family `(G_x)_x` is a further common refinement,
not data determined by that pair. Higher global factorial moments recover
the full multiplicity histogram and hence the exact uncovered count at
degree `p`; higher fibre-resolved tensors also retain its localization.
None of these objects by itself transports between thirteen-adic depths or
proves a uniform positive margin.

For a fixed rational interval packet, one may atomize at every endpoint
and clear atom denominators, converting a measured cover into the finite
multiset model above. Equations (19) and (23) are then the cheapest
count-preserving noncoverage certificates, while (15) is the first exact
label-intersection upgrade. This applies inside an atomized owner-label
search for THM-2222, but it is **not** by itself an upper bound for
THM-2222's fourfold intersection `S_4`: that problem uses the separate
upper moment LP in THM-2222. Any combined search must keep the coefficient
triple common across all checkpoints. Proving the required extremal
inequality remains OPEN. QED.
