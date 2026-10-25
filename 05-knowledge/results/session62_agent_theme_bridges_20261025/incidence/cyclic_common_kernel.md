# Full cyclic incidence kernel: wrap modes, path generators, and coefficient characteristic

2026-10-25. **PROVED** statements below; **FINITE-EXACT** verification is
separately scoped in §8. This closes the modular-wrap boundary of
`05-knowledge/results/trijection_hourglass_ports_20261025/incidence/common_kernel.md`
at commit `fde0dc637d81397de2b857e61d1e993ba11231bc`.
The earlier theorem is correct in its stated nonwrapping range; no correction
or retraction is needed. These are linear-incidence results, not assertions
that the arrays are binary errors, realizable sections, or norm-graph bounds.

## 1. Exact objects and the interval answer

Let `p>3` be prime, let `K` be **any coefficient field**, and initially let
`I={0,...,L−1}⊆F_p`, `1≤L≤p`. The geometry is over `F_p`, not necessarily K.
Put

```
u=(1,0),  v=(0,1),  w=(3,9),
T_I={(-a−3b,j−9b): a,b∈{0,1}, j∈I}.
```

The four x-coordinates `0,−3,−1,−4` are distinct because `p>3`; consequently
`T_I` has `4|I|` points even at full cyclic length. For `d∈{u,v,w}`, define
`Λ_d:K^{T_I}→K^{F_p}` by summing on each physical affine line parallel to d,
using cross-product labels. Empty physical lines give zero output coordinates.
Set

```
J=ker Λ_u ∩ ker Λ_v ∩ ker Λ_w,
W_u=ker Λ_v ∩ ker Λ_w,
W_v=ker Λ_u ∩ ker Λ_w,
W_w=ker Λ_u ∩ ker Λ_v.
```

Write `(x)_+=max(0,x)` and `δ=1` if `char K=p`, otherwise `δ=0`.

**Interval theorem, p≥17.** For every `1≤L≤p`,

```
dim_K J = (L−1)+(L−10)_+ +(L−p+6)_+ +(L−p+3)_+
          + δ·1_{L≥p−6}.                                      (1)
```

In particular the physical coefficient field `K=F_p` has a new cyclic mode
at `L=p−6`; the ordinary bridge-boundary corrections begin at `L=p−5`
and `L=p−2`. The previous formula survives `L=p−8,p−7`. At `L=p`,

```
dim_K J = 2p−2+δ.                                             (2)
```

**Small physical primes.** In the following table the ordinary row lists
`dim J` for `L=1,...,p` when `char K≠p`, including characteristic zero.
In characteristic p, add one exactly at the listed lengths.

| p | Ordinary dimensions, indexed by L=1,...,p | Add one when |
|---|---|---|
| 5 | 0,1,2,4,8 | L≥4 |
| 7 | 0,1,2,3,4,8,12 | L≥5 |
| 11 | 0,1,2,3,4,6,8,10,12,16,20 | L≥9 |
| 13 | 0,1,2,3,4,5,6,8,10,12,16,20,24 | L≥10 |

These rows are consequences of the general path theorem below, not extrapolated
rank tables. Small-prime exceptions concern cyclic shift overlap, not collisions
of the four physical columns. No coefficient-characteristic exceptions other
than p occur for J.

## 2. Lossless cyclic reduction

The following derivation works for **every nonempty subset** `I⊆F_p`, not only
an interval. Extend every column by zero off I, and use indices modulo p.
Write the four columns as `A=f_(0,0,·)`, `B=f_(0,1,·)`,
`C=f_(1,0,·)`, `D=f_(1,1,·)`. Define

```
h=A,  P=A+B,  R=A+C,  S=B+D,
(τ_d r)_j=r_(j−d),  Q=τ_3+τ_6+τ_9.
```

The physical u-line label is `j−9b`; the physical w-line label, after dividing
its cross product by 3 in **F_p**, is `j+3a`; the physical v-lines are the four
columns. Therefore the exact equations, with no no-wrap assumption, are

```
R_j+S_(j+9)=0,
P_j+(R+S−P)_(j−3)=0,
ΣA=ΣB=ΣC=ΣD=0.                                                (3)
```

The first gives `R=r`, `S=−τ_9r`. The second is

```
(1−τ_3)(P+Qr)=0,
```

because `(1−τ_3)Q=τ_3(1−τ_9)`. Since adding 3 permutes `F_p` in one cycle,
its invariant vectors are precisely the constants over **every** field K.
Thus for a unique scalar c,

```
R=r,  S=−τ_9r,  P=−Qr+c·1,
A=h,  B=P−h,  C=r−h,  D=S−P+h.                               (4)
```

Let `G=F_p\I` and `H=I∩(I−9)`. Support and the column sums become exactly

```
supp(h)⊆I,  supp(r)⊆H,
(Qr)_j=c  for j∈G,
Σh=0,  Σr=0,  p c=0 in K.                                    (5)
```

Indeed `ΣP=−3Σr+p c`. This is not division by 3 in K. Conversely (5) makes
all four arrays in (4) supported on I and makes every equation in (3) hold.
The construction is bijective, so neither dimension nor generators are lost.

## 3. The all-prime, arbitrary-puncturing path theorem

Construct a graph Γ on **all p residues**, with one undirected edge

```
{j−6,j−3}  for each j∈G.
```

Mark each vertex outside H as an **anchor**, required to have value zero.
The step-3 edges form one p-cycle. Since I is nonempty, Γ is a proper subgraph
of that cycle and hence a disjoint union of paths, including isolated vertices.
On a gap index j, `r_(j−9)=0`; the boundary equation is therefore simply

```
r_(j−6)+r_(j−3)=c.                                           (6)
```

A component is **free** if it has no anchor. Define:

* `k`: number of free components;
* `o∈{0,1}`: whether at least one free component has odd vertex count;
* `b∈{0,1}`: whether, in each anchored component, all anchors lie in the same
  bipartition class. Empty anchor sets impose no condition.

**Path theorem.** For every prime `p>3`, every field K, and every nonempty
`I⊆F_p`,

```
dim_K J = |I|−1 + k−o + δ b o.                               (7)
```

This includes `I=F_p`: then Γ has p isolated free vertices, so `k=p,o=b=1`,
and (2) follows immediately.

**Proof and constructive basis.** When `c=0`, an anchored component is zero.
Each free path C supports its alternating vector `ψ_C`, with entries
`+1,−1,+1,...` in endpoint order. These vectors have sum zero for even-sized
paths and sum one for odd-sized paths, in every characteristic, including 2.
If no odd free path exists, all k vectors obey `Σr=0`. Otherwise choose one
odd free path O: the even-path vectors and `ψ_C−ψ_O` for the other odd paths
form a basis of dimension `k−1`. This proves the homogeneous term `k−o`.

If `char K≠p`, (5) forces `c=0`. Suppose `char K=p`. A path equation propagates
values `t,c−t,t,c−t,...`. Opposite-parity anchors force `c=0`; same-parity
anchors permit every c. Hence a nonzero-c boundary solution exists iff `b=1`.
For `b=1`, construct an integral `c=1` solution `r^0` by placing alternating
0/1 values on every component, with all anchors zero (choose either parity on
free paths). Put `β=Σr^0`. If `o=1`, then

```
(r,c)=(r^0−βψ_O,1)                                           (8)
```

obeys all of (5), giving exactly one new dimension. If `o=0`, all homogeneous
solutions have sum zero, so the sum condition would require `βc=0`. Here I
must be proper: full I has odd isolated free paths. There is at least one edge,
so `β≥1`; all ones lie in H, so `β≤|H|≤p−1`. Thus β is nonzero in characteristic
p, forcing `c=0`. This proves (7).

To obtain a complete basis of J, lift these parameter pairs using (4), and add
`|I|−1` h-generators `h=e_(i_t)−e_(i_(t+1))` for any ordering of I. The h and
(r,c) parameters are independent by the bijection. This also proves the
small-prime table by finite path counting, with no field-by-field elimination.

## 4. Evaluating intervals for p≥17

The nonwrap calculation from the source gives `L−1+(L−10)_+` for `L≤p−9`.
It can also be read directly from (5): `c=0`,
`r` is supported on `{0,...,L−10}`, and only its sum must vanish.

For `m=p−L` with `1≤m≤8`, write

```
H={0,...,L−10} ⊔ {p−9,...,L−1},
B=L−9,  x_t=r_(p−9+t),  0≤t≤8−m.                             (9)
```

The first block has B vertices, interpreted as empty if B=0. For p≥17 it is
never negative in this wrap range. Gap equations involve only x. Their complete
reduced forms and homogeneous dimensions are:

| m | Independent/forcing equations | d_m=dim(x-space at c=0) |
|---|---|---|
| 1 | x_2+x_5=c | 7 |
| 2 | x_1+x_4=c, x_2+x_5=c | 5 |
| 3 | x_0+x_3=c, x_1+x_4=c, x_2+x_5=c | 3 |
| 4 | x_2=c, x_0+x_3=c, x_1+x_4=c | 2 |
| 5 | x_1=c, x_2=c, x_0+x_3=c | 1 |
| 6 | x_0=x_1=x_2=c | 0 |
| 7 | c=0, x_0=x_1=0 | 0 |
| 8 | c=0, x_0=0 | 0 |

Thus `d_m=(6−m)_++(3−m)_+`. For m≤6, the particular boundary solution
`x_0=x_1=x_2=c`, all other x zero, always works. The low block is nonempty
there, so `Σr=0` removes one independent constraint and the characteristic-p
constant supplies precisely one extra parameter. For m=7 the low block has
at least one vertex, again imposing one sum constraint. For m=8 the only
zero-low-block case is `p=17,L=9`; the lone x is already killed, so no spurious
negative dimension occurs. This gives

```
dim J=(L−1)+(L−10)_+ + d_m + δ·1_{m≤6}
```

throughout `1≤m≤8`. At m=0 use (2), which equals (1). This proves (1).

### Explicit interval generators

The h-generators are adjacent cube differences. In the nonwrap range, the
r-generators are consecutive differences in the low block: these are the
source's twelve-point bridges. In the wrap range use the following x-bases:

* m=1: `e_0,e_1,e_3,e_4,e_6,e_7,e_2−e_5`;
* m=2: `e_0,e_3,e_6,e_1−e_4,e_2−e_5`;
* m=3: `e_0−e_3,e_1−e_4,e_2−e_5`;
* m=4: `e_0−e_3,e_1−e_4`;
* m=5: `e_0−e_3`;
* m=6,7,8: none.

Take the low-block adjacent differences. For every listed x-vector v, extend
it by zero and subtract `(Σv)e_0` in the low block. This enforces `Σr=0` and
provides all ordinary wrap generators. The low anchor always exists when a
listed x-vector is present. In characteristic p and m≤6 add

```
r=e_(p−9)+e_(p−8)+e_(p−7)−3e_0,   c=1.                        (10)
```

At L=p instead take `r=e_i−e_(i+1)` for `0≤i≤p−2`, together with `(r,c)=(0,1)`
only in characteristic p. Lift all these vectors by (4). They are independent
and complete. The arbitrary-puncturing path basis covers the small primes.

## 5. All three pairwise incidence spaces

For any nonempty I put `n=|I|`,

```
η_s=|I∩(I+s)|,  c_s=2 if η_s=0, otherwise 1,
κ=# components on I with edges of differences ±3 or ±9 modulo p.
```

Then, over every field,

```
dim W_u=2n+η_3−4+c_3,
dim W_v=η_3+η_9+κ,
dim W_w=2n+η_9−4+c_9.                                       (11)
```

Proof: for two line partitions form the bipartite incidence graph with support
points as edges. Multiplying one part's rows by −1 identifies its matrix with
an oriented graph incidence matrix; its rank is `#vertices−#components` in
every characteristic, by a spanning forest. For v versus w the two pairs of
columns have overlap η_3, giving the first expression; u versus v gives η_9
and the third. For u versus w, each fixed j is a connected four-edge square.
Squares j and k meet exactly when `k−j∈{±3,±9}`, yielding κ; their total vertex
count is `(2n−η_3)+(2n−η_9)`, which gives the middle expression.

For an interval and `t=s mod p∈{1,...,p−1}`,

```
η_s=(L−t)_++(L−(p−t))_+.
```

For p≥17, additionally

```
κ=min(3,L)−1_{L≥p−8}−1_{L≥p−7}.
```

Before wrap the squares have the three ordinary residue classes mod3 (or L
isolated classes for L<3). The first wrap-9 join merges two classes and the
second merges the third. The general κ definition handles every small prime.

## 6. Coupling, and the necessary zero-weight boundary

For **nonzero** `s_u,s_v,s_w∈K`, define

```
M:W_u⊕W_v⊕W_w→K^{T_I},
M(F_u,F_v,F_w)=s_uF_u+s_vF_v+s_wF_w.
```

Applying Λ_u to a zero weighted sum forces `F_u∈J`; cyclically all three
components belong to J. Therefore

```
ker M={(F,G,−(s_uF+s_vG)/s_w):F,G∈J},
dim ker M=2 dim J,
rank M=dim W_u+dim W_v+dim W_w−2dim J.                        (12)
```

At full length, the three W dimensions are `(3p−3,2p+1,3p−3)`, so
nonzero-weight M has rank `4p−1−2δ`: image codimension3 in characteristic p
and codimension1 otherwise.

A basis consists of the two channels of every J-generator. For the source's
physical weights `(-3,-9,1)` over `F_p`, these are `(F,0,3F)` and `(0,F,9F)`.

If those numerical weights are reused over characteristic 3, the first two
weights are **zero**. Then instead `ker M=W_u⊕W_v⊕{0}`, `rank M=dim W_w`.
For example geometry p=17,L=2 over K=F_3 has `dimJ=1` and each W dimension2,
but this degenerate coupling has nullity4, not2. Formula (7) for J is unaffected.

## 7. Why the characteristic-p term is not a Fourier shortcut

No diagonalization was used. The cyclic shift always has a one-dimensional
fixed space, over every K. If `char K≠p`, its minimal polynomial `X^p−1` is
separable and a splitting-field Fourier view is available. If `char K=p`,
`X^p−1=(X−1)^p`; the p-cycle is one size-p Jordan block (it has a cyclic vector),
so a semisimple Fourier argument would be invalid. The actual change in (5)
is that the invariant constant has total sum `pc=0`. The path parity and sum
sidecar determine whether this potential mode survives the support boundary.
Counting constants alone would wrongly predict it at every wrap length.

The information retained by (4) is the full labelled physical line partition,
not only point count or the abstract box coordinates. The path reduction
preserves every equation and restores omitted support through anchored vertices.
Forgetting those anchors changes J; forgetting `pc=0` erases its field dependence.

## 8. Reproducibility and hostile controls

`verify_wrap.py` is a standard-library-only checker. Its physical matrices are
constructed directly from the point coordinates and physical cross products;
the candidate dimensions and basis use the separate path construction. All
checks raise explicit errors rather than using disabled Python assertions.
It verifies support, every physical line equation, independence, dimensions of
J and all three W spaces, and selected literal weighted coupling ranks. A
second, separately constructed sparse interval basis checks the explicit table
and generators in §4, without invoking the forest.

```
python verify_wrap.py --output verification.json
python -O verify_wrap.py --output verification_optimized.json
```

The declared universe is every interval `1≤L≤p` for
`p∈{5,7,11,13,17,19,23,29,31,37,41,43}`, over coefficient prime fields with
characteristics in `{2,3,5,p,101}` (duplicates removed), plus **every nonempty
subset** of F_p for p=5,7,11 over characteristics `{2,3,p,101}`. There is no
symmetry quotient, support filter, or genericity assumption. The full run has
10,195 direct field cases and 550 direct weighted coupling-rank checks.

Explicit hostiles include p=29,L=23: dimension36 over F_29 but35 off that
characteristic; p=29,L=24: ordinary dimension38 rather than the unwrapped37;
and full p=29: dimensions57/56 rather than the old expression47. The p=17,L=2,
characteristic3 weight example above blocks misuse of the nonzero-weight lemma.

Independent auditing also checked 3,525 direct physical-matrix cases and all
235 full physical incidence matrices by exact integer Smith normal form for
p∈{5,7,11,13,17,19,23,29,31,37,43}, every 1≤L≤p. In that **finite universe** all
nonzero Smith factors were1 except one p precisely when predicted. These finite
controls are supplementary; the proofs in §§2–6 establish the stated all-prime,
all-field and arbitrary-puncturing claims.

## 9. Inheritance and scope ledger

* Closest proved mechanism: the source's common-kernel lemma and cube-plus-bridge
  parametrization; recovered directly at the stated commit.
* Canonical hostile: disjoint support need not give independent physical
  incidence, and cubes alone fail from the first bridge length.
* Corrected near miss: integer finite-support uniqueness cannot be carried into
  a cyclic group without retaining the constant mode and its total sum.
* Least-used sidecar: actual modular line labels, now compressed losslessly to
  path components together with anchors and their parity.
* New operation: boundary wrap, then arbitrary puncturing of the vertical index
  set. Preserved predicate: all three physical zero-line-sum equations.
* No norm-graph failure estimate, binary-error realizability, LRC implication,
  or identification with a different geometrical object is asserted.
