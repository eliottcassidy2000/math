---
id: THM-4115
title: "Uniform ear-cut Walsh variance and sharp growth refinement"
status: >
  PROVED ELEMENTARY WALSH/PARSEVAL + VERIFIED-EXACT + INDEPENDENTLY
  VERIFIED-EXACT. The complete labelled cut fibre has an exact degree-two
  Walsh expansion and variance. Its lower support at H(T) gives a
  second-moment maximum bound and improves every full-cut selected-bank
  recurrence from (n+3)/4 to (n+1)(n+2)/(4n). Equality holds for C3, so the
  order-three coefficient is sharp. For strong n>=4, a cyclic-order argument
  also forces F_1>=n-1 in even order and F_1>=2 in odd order. Cut-image
  intervals, overlap, and global H-spectrum completeness remain OPEN.
source: codex-padic-zeta-irrationality-20260825
depends_on:
  - THM-001-redei
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
related:
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4102-selected-order-ten-strong-ear-solid-interval
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-4118-ear-response-lattice-and-stateful-unit-component-intervals
  - HYP-2879-strong-ear-atom-calculus
  - HYP-9029-strong-interval-tiling-law
  - MISTAKE-402
script: 04-computation/tournament_ear_cut_walsh_variance_thm4115.py
output: 05-knowledge/results/tournament_ear_cut_walsh_variance_thm4115.out
independent_audit_script: 04-computation/tournament_ear_cut_walsh_variance_thm4115_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_ear_cut_walsh_variance_thm4115_independent_audit.out
script_sha256: 27546bc169cf07ce95260afa56322ae8e6b44432e379e549605bf28f85225faa
output_sha256: 25c0005e8743b57ba849e5bc42e40c9ae55fa7b8106bebe8bef4e8af37dcc8b8
independent_audit_script_sha256: 0007b1b52ee39bc169668580959bac770f83fcd0604a41bdf9097446a23dc1bd
independent_audit_output_sha256: 287560799d8be6550d2e409449661900885e8cf52333c8e42e68e39c084dd46d
semantic_sha256: 7580aaf0fd320afad386d6401b3e95dea31c91fc2cb3fc67b085773dcdf0836e
independent_semantic_sha256: d805f598ce72128f22bd46bf4c8ba696b471be4934455af00b02837fa508cc0c
hash_basis: raw LF bytes for files; canonical compact JSON for semantic ledgers
audit: >
  PASS. The primary subset-DP path exhausts all 1,098 labelled parents and
  33,864 cut ears through order five and freezes selected order-nine and
  order-ten controls. The independent path imports no primary code and
  literally scans 124,468 parent orderings and 23,717,424 child orderings.
  Both normal/-O streams byte-match their LF frozen outputs and agree on the
  Walsh reconstruction, variance, support inequality, full F_1 refinement,
  sharp C3 equality, boundary hostiles, the maximum-neighbor local-step
  obstruction, and recursive product.
---

# THM-4115 -- uniform ear-cut Walsh variance and sharp growth

**PROVED ELEMENTARY WALSH/PARSEVAL + VERIFIED-EXACT + INDEPENDENTLY
VERIFIED-EXACT.**

THM-4111 averaged the complete cut fibre and thereby retained only its
constant Walsh coefficient. This theorem restores the complete degree-one and
degree-two coefficients. Their squared norm gives enough dispersion to force
a larger cut value at every order. It still does not locate that value inside
a solid interval.

## 1. The cut field as a degree-two Walsh polynomial

Let `T` be a tournament on `V={1,...,n}`, with `n>=2`, and write

```text
H=H(T),
X_S=H(T+x_S)                 (S subseteq V).               (1)
```

Use THM-4097's convention `x->i` exactly when `i in S`. Its integral cut
field is

```text
X_S=H+cut_w(S)+sum_(i in S)h_i,
w_ij=(Q_ij+Q_ji)/2 >= 0,     sum_i h_i=0.                 (2)
```

Put

```text
epsilon_i(S)=+1 if i in S and -1 otherwise,
W=sum_(i<j)w_ij.                                              (3)
```

The indicator of membership is `(1+epsilon_i)/2`, while the indicator that
an edge crosses the cut is `(1-epsilon_i epsilon_j)/2`. Equation `(2)` is
therefore the exact Walsh expansion

```text
X_S=H+W/2
       +(1/2)sum_i h_i epsilon_i
       -(1/2)sum_(i<j)w_ij epsilon_i epsilon_j.             (4)
```

Thus the cut response has no Walsh degree above two. This is an identity for
every tournament and every labelled cut, not a finite-order observation.

## 2. Exact mean and variance

Let expectation be uniform over all `2^n` subsets. Distinct Walsh characters
are orthogonal. Equation `(4)` gives

```text
mu=E[X_S]=H+W/2,                                             (5)
Var(X_S)=1/4 (sum_i h_i^2+sum_(i<j)w_ij^2).                 (6)
```

Comparing `(5)` with THM-4111's exact all-cut mean yields

```text
W=((n-1)H+F_1(T))/2,                                        (7)
```

where `F_1(T)` counts old-vertex orderings with exactly one bad adjacency.
In particular `W>0` for `n>=2`.

## 3. Every cut has lower support H

For every Hamiltonian path `v_1,...,v_n` of `T`, record the binary cut word
along that path. There is always a legal insertion slot for `x`: either the
word starts with `1`, ends with `0`, or contains an adjacent transition
`0,1`. Choose the first legal slot in a fixed priority order. Inserting `x`
there produces a Hamiltonian path of `T+x_S`, and deleting `x` recovers the
old path. Hence this construction is injective and

```text
X_S>=H for every S subseteq V.                              (8)
```

The constant cuts attain equality: `X_empty=X_V=H`.

## 4. The support-sensitive second-moment bound

Let `M=max_S X_S`. Pointwise on the cut cube,

```text
(X_S-H)(M-X_S)>=0.                                          (9)
```

Averaging and using `(5)--(6)` gives

```text
(mu-H)(M-mu)>=Var(X_S).
```

Because `mu-H=W/2>0`,

```text
M >= mu+Var(X_S)/(mu-H)
  = H+W/2+(sum_i h_i^2+sum_(i<j)w_ij^2)/(2W).              (10)
```

This is the exact variance-sensitive floor. It uses both the Walsh norm and
the lower-support sidecar; a mean or variance in isolation would not imply it.

Put `b=binom(n,2)`. Cauchy's inequality on the nonnegative edge weights and
deletion of the nonnegative field energy give

```text
sum_(i<j)w_ij^2 >= W^2/b.                                  (11)
```

Substitution of `(7)` into `(10)--(11)` proves the full one-defect refinement

```text
M >= ((n+1)(n+2)/(4n)) H
     +((n(n-1)+2)/(4n(n-1))) F_1(T)                        (12)
  >= ((n+1)(n+2)/(4n)) H.                                  (13)
```

### Strong-parent additive surplus

There is a further all-order consequence when `T` is strong and `n>=4`.
Let `C_k(T)` count cyclic orderings of `V(T)`, modulo rotation, with exactly
`k` bad cyclic adjacencies, and let `F_r(T)` count linear orderings with
exactly `r` bad adjacencies. Breaking a cyclic order at a good edge preserves
its defect count, while breaking it at a bad edge lowers that count by one.
Therefore

```text
F_r=(n-r)C_r+(r+1)C_(r+1),
H=nC_0+C_1,             F_1=(n-1)C_1+2C_2.                (13a)
```

We need one elementary extension lemma. Suppose a cyclic order
`v_1,...,v_m` has `k in {1,2}` bad edges and a new vertex `x` is added. Put
`b_i=1[x->v_i]` and let `d_i` indicate that the cyclic edge
`v_i,v_(i+1)` is bad. Inserting `x` in gap `i` changes the defect count to

```text
k'=k-d_i+b_i+(1-b_(i+1)).                                  (13b)
```

For `k=2`, a nonconstant cyclic bit word has a `0->1` transition, which gives
`k'=2-d_i in {1,2}`; if the word is constant, inserting at either bad gap
gives `k'=2`. For `k=1`, insert at the unique bad gap unless it is a `0->1`
transition. In that exceptional case, an equal-bit gap gives `k'=2`; if no
such gap exists, the bits alternate, and a second (good) `0->1` gap gives
`k'=1`. The latter case is used only from `m>=4` onward.

Every strong tournament contains a directed triangle `a->b->c->a`. Add any
fourth vertex `x`. The cyclic order `(a,b,c,x)` already has one or two bad
edges unless both `c->x` and `x->a` hold. In that remaining case
`(a,b,x,c)` has the bad edge `x,c` and at most one other bad edge. The
extension lemma then adds every remaining vertex. Hence `C_1+C_2>0` for every
strong `T` of order at least four.

If `n` is even, THM-001 and `H=nC_0+C_1` force `C_1` to be positive odd, so
`F_1>=n-1`. If `n` is odd, `(13a)` makes the now-positive `F_1` even, so
`F_1>=2`. With

```text
eta_n=n-1 if n is even, and eta_n=2 if n is odd,            (13c)
```

equation `(12)` has the strict strong-parent refinement

```text
M >= ((n+1)(n+2)/(4n))H
     +((n(n-1)+2)/(4n(n-1)))eta_n,          n>=4 strong.   (13d)
```

Equality in the `F_1` floor would require `(C_1,C_2)=(1,0)` in even order or
`(0,1)` in odd order. No proportional bound such as `F_1>=H` is asserted
here.

The factor in `(13)` is greater than one. Therefore a maximizing cut is
nonconstant. If `T` is strong, that child is strong by THM-4097. By THM-001
every `X_S` is odd, so `M` is at least the least odd integer above either
rational floor.

For `T=C_3`,

```text
H=3, F_1=0, W=3, h=0, w_12=w_13=w_23=1,
mu=9/2, Var=3/4, M=5.                                      (14)
```

Equality holds in `(10)`, `(12)`, and `(13)`. Thus the order-three universal
coefficient `5/3` is sharp. The hypothesis `n>=2` is necessary for the
displayed quotient: at `n=1`, `W=mu-H=0`.

## 5. Refined selection-robust recursive growth

Use exactly the full-cut selected-bank construction of THM-4111. Start with a
nonempty bank `B_n` of strong order-`n` tournaments; expand every nonconstant
cut from every retained parent; then keep an arbitrary strong witness for each
attained scalar value. If

```text
M_n=max{H(T):T in B_n},                                     (15)
```

then `(12)--(13d)` applied to a maximizing parent gives, for `n>=4`,

```text
M_(n+1) >= oddceil(
  ((n+1)(n+2)/(4n))M_n
  +((n(n-1)+2)/(4n(n-1)))eta_n)
 >= oddceil(((n+1)(n+2)/(4n))M_n).                         (16)
```

At the initial strong order `n=3`, the second line remains valid and is exact
for `C_3`; the additive `eta_n` statement begins at order four.

Ignoring the helpful odd rounding and iterating from `n` to `m>n` yields

```text
M_m >= M_n * m(m+1)! / (n(n+1)! 4^(m-n)).                  (17)
```

This strictly improves THM-4111's `(n+3)/4` recurrence and still tends to
infinity. It is selection-robust because the full cut fibre is expanded before
one representative per scalar value is retained.

For the selected frontier parents already used by THM-4097/4102/4104:

| parent order | `H` | exact variance floor, oddceil | universal oddceil | actual cut maximum |
|---:|---:|---:|---:|---:|
| `9` | `3,357` | `14,481` | `10,259` | `15,487` |
| `10` | `15,621` | `76,225` | `51,551` | `93,751` |

The exact variance floors use the full right side of `(10)`. The observed
maxima are finite-exact controls, not consequences of the inequalities.

## 6. What the restored moment does and does not recover

The strong order-five codes `1015` and `759` from THM-4111 have the same

```text
(H,F_1,mu)=(9,30,51/2),                                    (18)
```

but variances `305/4` and `315/4`, respectively. Thus the second moment
separates the known equal-mean hostile. Their exact variance floors are
respectively `994/33` and `333/11`, both with odd ceiling `31`, while their
actual maxima are `41` and `43`. The first image omits the internal odd value
`21`. Thus even the exact variance floor does not imply a solid image. The
exhaustive primary audit finds no equal-`(H,F_1,Var)` pair with different full
cut image through order five; that finite absence is not a determinacy
theorem.

A separate **FINITE-EXACT scout** enlarges the order-nine selected bank before
comparison. It retains the first at most 32 distinct labelled witnesses in
each of all `1,482` exact `H`-fibres arising from the complete
`6,008*254=1,526,032` strong order-eight ear stream. Among `46,314` parents it
finds `11,938` repeated `(H,F_1,Var)` groups and `88,527` equal-triple parent
pairs, with zero different-image pairs. This is a high-multiplicity hostile
search, not an all-order injectivity result. The canonical one-per-`H` bank by
itself supplies no such evidence because its `H` coordinate is already
injective. Reproduce with

```bash
python -B 04-computation/order9_ear_variance_image_determinacy_cap32_scout.py
python -B -O 04-computation/order9_ear_variance_image_determinacy_cap32_scout.py
```

There is also a direct hostile to descending from the maximum by one-bit
`-2` moves. At code `1015`, the positive drops from a maximizing cut to a
one-bit neighbor are exactly `{4,8,18}`; at code `759` they are
`{6,10,12,14}`. Exhaustively among labelled strong parents, the numbers having
no maximum with an `M-2` neighbor are

```text
order 3: 0/2,       order 4: 0/24,       order 5: 400/544.  (19)
```

This is **FINITE-EXACT**, not an all-order theorem. It already refutes the
naive strategy “find the large value by moments, then walk downward through
single-bit `-2` flips.” Any interval mechanism must instead use overlap among
multiple cuts/parents, paired flips, or another labelled local sidecar.

The complete `(H,h,w)` Walsh carrier in `(4)` determines every cut value.
Passing that carrier to the scalar variance in `(6)` retains only the squared
norm of its degree-one and degree-two coefficients. This quotient destroys
their labelled signs, incidence, higher moments, small-ball profile, and the
locations of attained values. Consequently this theorem proves none of:

- a solid interval around or below the maximum;
- overlap of recursively selected solid intervals;
- a complete order-`n` spectrum at any new order; or
- global H-spectrum completeness.

All remain **OPEN**. THM-4114 identifies the symmetric pair curvature
`Delta_i Delta_j X=-(Q_ij+Q_ji)=-2w_ij`. If
`K=Q+Q^T=2w`, then the degree-two Walsh coefficient is `-K_ij/4` and

```text
Var(X)=1/4 sum_i h_i^2+||K||_F^2/32.                        (20)
```

Its full directed boundary network retains `Start/End/Q`, whereas
`(6)` contracts both `h` and `w` to one squared norm. Nothing here treats that
curvature or variance as determining the cut image.

## 7. Exact controls and reproduction

The primary referee uses subset path DP and the THM-4097 boundary matrix. It
checks every labelled parent through order five:

```text
1,098 parents, 33,864 cuts, 570 strong parents,
16,668 nonconstant ears over strong parents.                (21)
```

It verifies every individual Walsh value, both moments, lower support,
second-moment floors, parity, strong-ear inheritance, the `n=1` denominator
hostile, the transitive `n=2` factor control, C3 equality, the equal-mean
hostile, the two large selected rows, and the recursive product algebra.

The independent referee imports no primary code. It literally enumerates
`124,468` parent orderings and `23,717,424` child orderings over the same cut
universe, reconstructs the Boolean polynomial by singleton and second
differences, and reproduces every structural and boundary check.

Run

```bash
python -B 04-computation/tournament_ear_cut_walsh_variance_thm4115.py
python -B -O 04-computation/tournament_ear_cut_walsh_variance_thm4115.py
python -B 04-computation/tournament_ear_cut_walsh_variance_thm4115_independent_audit.py
python -B -O 04-computation/tournament_ear_cut_walsh_variance_thm4115_independent_audit.py
```

Both normal/-O streams must byte-match their LF frozen outputs. The finite
audits are hostile controls for the elementary identities and inequalities;
they are not the source of their all-order quantifiers.
