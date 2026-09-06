# Injective labels and a genuine conditional cycle-cumulant expansion

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent referee](overnight6_20260906_no3line_palettes_audit.md)
accepts the full proof and exact controls, with a third board-generation path.
This is a structural continuation of the cycle-defect work, not a new
evaluation of its n=8 coefficients or a no-three-in-line existence theorem.
The general conditioning and cumulance identities are classical. The new
repository synthesis identifies their exact geometric input, local cycle
kernel, and decisive small failures of discarding the conditioning variable.

## 1. Inheritance and the coordinate restored

The closest proved mechanism is
`05-knowledge/results/overnight4_20260906_no3line_cycle_defect.md`, with its
independent audit. Its universal type ring has an exact connected formal
logarithm, but the geometric weighting map is not a positive moment state
for disjoint-union multiplication. Its homogeneous forest square has value
`-397/529200`. Its n=8 third-cumulant contrast is `11881/6300`, refuting
addition of scalar contributions from individual skeleton cycles.

The older `overnight_20260906_no3line.md` supplies the correct uniform
shore-label model and the n=4 cycle-blind-mean hostile. The short-cycle
copy theorem in `overnight_20260906_moments_pairprofile_theorem.md` supplies
the local non-induced embedding normalization. The corrected near miss is
to confuse a formal type product with multiplication of actual random
observables. The least-used relevant sidecar is the **actual set of row
and column labels assigned to each named cycle**, before those labels
are permuted within the cycle.

Board: formal path defects; complete event unions; injective labels;
component palettes; conditional independence; ordinary/factorial cumulants.
The precise repair below retains the label coordinates until after genuine
products and conditional cumulants have been formed.

## 2. The random variable and the all-size palette theorem

Fix a simple bipartite 2-regular graph `G` with n vertices on each shore.
Name its connected components `G_i=C_(2k_i)`, where `k_i>=2` and
`sum k_i=n`. Choose independent uniform bijections `rho` from the left
shore to `{0,...,n-1}` and `sigma` from the right shore to the same set.
The board is

```text
B={(rho(u),sigma(v)): uv is an edge of G}.
```

Let `T_n` be all unordered triples of grid points with distinct row and
column coordinates and **zero integer determinant**. Define
`I_T=1[T subset B]` and `X=sum_(T in T_n) I_T`. Thus X counts nonaxis
collinear triples; axes cannot contribute because B has two points in
every row and column. An event T can use edges of one, two, or three
skeleton components. Nothing here reduces collinearity modulo a prime.

Write `R_i=rho(left(G_i))` and `C_i=sigma(right(G_i))`. The palette variable
`P=((R_i,C_i))_i` retains the association to each named component, including
when several cycle lengths agree. The row palettes are uniform ordered
partitions of their prescribed sizes; column palettes are an independent
such partition. Every palette pair has probability

```text
(product_i k_i! / n!)^2.
```

**Theorem (whole-event local factorization).** Conditional on P, the
internal left/right bijections in different components are independent and
uniform. For any finite family of grid events, let F be their complete
union of grid edges. If an edge `(r,c)` has r in `R_i` and c in `C_j` with
`i!=j`, its conditional presence probability is zero. Otherwise put
`F_i=F intersect (R_i x C_i)`, and let `r_i,s_i` be its numbers of incident
rows and columns. Then

```text
Pr(F subset B | P)
 = product_i inj_sh(F_i,C_(2k_i)) / ((k_i)_(r_i)*(k_i)_(s_i)).       (1)
```

Here `inj_sh` counts edge-preserving, shore-preserving injective maps from
the vertices of F_i to those of the fixed abstract cycle; empty factors
are one. There is no induced-subgraph condition or shore-exchange quotient.
The equivalent expression is `copies(F_i,C_(2k_i))/N_(k_i)(F_i)`, where
`N_k(H)=(k)_(r(H))*(k)_(s(H))/aut_sh(H)`. Both numerator and denominator
divide by the same complete shore-preserving automorphism group.

**Proof.** Conditioning specifies the images of each component's vertices,
but each of its `k_i!` left bijections and `k_i!` right bijections remains
possible. The conditional sample space is their Cartesian product, of
size `product_i (k_i!)^2`, with constant weight. The inverse images of the
r_i specified row labels and s_i specified column labels are independently
uniform injections, so precisely the numerator in (1) of the
`(k_i)_(r_i)*(k_i)_(s_i)` injections contain all edges. Independence across
components proves the product. This also automatically rejects every local
degree-three union and every incompatible cycle or size. Repeated events
are harmless: their complete union, rather than a list with multiplicity,
controls the joint presence probability. QED.

For a single triple, this gives an especially small positive conditional
mean formula. If T is palette-compatible, let `s_i` be how many of its
three points lie in rectangle `R_i x C_i`. Then

```text
E[X | P] = sum_(compatible T in T_n) product_i q_(k_i,s_i),          (2)
q_(k,0)=1,
q_(k,1)=2/k,
q_(k,2)=2(2k-3)/(k(k-1)^2),
q_(k,3)=4(2k-5)/(k(k-1)^2(k-2)) for k>=3;
q_(k,s)=0 for s>k.
```

Indeed `q_(k,s)=m_s(C_(2k))*s!/(k)_s^2`, where m_s is the number of
s-edge matchings. For s=1,2,3 the counts are `2k`, `k(2k-3)`, and
`2k(k-2)(2k-5)/3`; the last is used only for k>=3. Distinct coordinates
of a line triple are essential for this matching reduction.

## 3. The genuine connected expansion, and the outer term it needs

For a fixed P, an event depends only on the internal permutation variables
of the cycles containing its points. In fact four-cycles may be removed
from this dependency set: `C4=K_(2,2)`, so its entire palette rectangle is
already present deterministically. Let `S_P(T)` retain only components
with `k_i>=3` in that event; impossible events are identically zero.

For any r>=2 and any list `T_1,...,T_r`, **with repetitions allowed**, form
the graph on the r occurrences by joining two when their sets S_P intersect.
If this graph is disconnected, then

```text
kappa(I_(T_1),...,I_(T_r) | P)=0.                                  (3)
```

An occurrence with empty support is conditionally constant and also has
zero joint cumulant of order>=2. Independence of disjoint collections of
component variables proves (3): their multivariate moment generating
function factors, so its logarithm splits, and a coefficient involving
variables from both factors is zero. Connectedness is necessary, not
sufficient, for a nonzero conditional joint cumulant. By multilinearity,
`kappa_r(X|P)` is the sum of these kernels over all ordered r-tuples of
events. This is an actual cumulant expansion, not the logarithm of the
formal component-type polynomial.

For three events, put `p_J=Pr(union_(j in J) T_j subset B | P)`. Formula
(1) gives the exact connected kernel

```text
p_123-p_12*p_3-p_13*p_2-p_23*p_1+2*p_1*p_2*p_3.                  (4)
```

Unconditioning requires more than averaging (3). With
`mu_P=E[X|P]`, `v_P=Var(X|P)`, the exact third identity is

```text
kappa3(X) = E[kappa3(X|P)] + kappa3(mu_P) + 3 Cov(mu_P,v_P).       (5)
```

For a direct proof, write `X-E[X]=(X-mu_P)+(mu_P-E[X])`, cube, and
condition. The term linear in `X-mu_P` vanishes; its square gives v_P,
and its cube gives the conditional third cumulant. This proves (5) with
no independence premise. At all orders the corresponding partition
identity is the classical law of total cumulance: each partition block
first supplies a conditional joint cumulant, and an outer joint cumulant
combines those block variables. The checked primary source is
[Brillinger, *The calculation of cumulants via conditioning*, AISM 21
(1969), 215-218](https://www.ism.ac.jp/editsec/aism/pdf/021_1_0215.pdf),
equation (4) and its theorem on p.215; its proof and univariate corollary
are on p.216. Only this finite-moment identity is imported.

## 4. Two minimal controls and one nondegenerate decomposition

For n=4 and `G=2C4`, P already fixes the entire board. The 36 equiprobable
pairs of palettes for the first named component have histogram

| X | 0 | 2 | 4 | 6 | 8 |
|---|---:|---:|---:|---:|---:|
| Number of palettes | 18 | 8 | 4 | 4 | 2 |

They describe 18 distinct boards, each counted twice. Hence

```text
E[X]=2, Var(X)=56/9, kappa3(X)=16,
E[kappa3(X|P)]=0, kappa3(mu_P)=16, Cov(mu_P,v_P)=0.               (6)
```

Each four-cycle alone has no nonaxis line triple; every event in (6)
spans components. This is the smallest disconnected simple bipartite
2-regular skeleton: each component needs at least two vertices per shore.
Thus a claim that the genuine cumulant is just the average of the
conditional connected internal-cycle terms fails at the first admissible
disconnected size. The issue is the palette mixture, not residual local
dependence overlooked by (3).

There is also a sharp typing trap. For deterministic X=x, the ordinary
conditional cumulants of order>=2 vanish, but
`log E[(1+t)^X | P]=x log(1+t)`. Its factorial cumulants are
`(-1)^(r-1)*(r-1)!*x`; in particular its third is 2x. In (6) the average
conditional third factorial cumulant is 4, whereas the unconditional third
factorial cumulant is `4/3`. Ordinary and factorial conditional cumulants
cannot be interchanged. Event repetitions in the logarithm cannot be
silently dropped.

At n=5, `G=C4+C6`, there are 100 equally likely palettes and six equally
likely boards per palette. Direct exact enumeration gives

```text
E[X]=13/3, Var(X)=1769/225, kappa3(X)=26159/675,
E[kappa3(X|P)] = -871/675,
kappa3(mu_P)   = 6467/360,
3 Cov(mu_P,v_P)= 7949/360.                                      (7)
```

All three terms in (5) matter; averaging conditional cumulants alone even
gives the wrong sign. The board-generation path is just a complete 2x2
palette rectangle plus a 3x3 rectangle with one perfect matching removed.
An independent path using all `5!^2=14400` shore labelings of a fixed
`C4+C6` reproduces every conditional histogram and the 600 distinct boards.

Conditional connected interactions across cycles also genuinely occur.
For n=6 with two C6, choose row and column palettes `{0,1,2}` and
`{3,4,5}`. Each component is its 3x3 rectangle minus a uniform one of six
perfect matchings, independently. Let A be the diagonal triple on 0,1,2,
B the diagonal triple on 3,4,5, and C the diagonal triple on 1,2,3. Then

```text
(p_A,p_B,p_C,p_AB,p_AC,p_BC,p_ABC)
  = (1/3,1/3,1/3,1/9,2/9,1/6,1/9).
kappa(A,A,B | P)=0,           kappa(A,B,C | P)=1/54.              (8)
```

A and B depend on different independent components; C meets both. Formula
(4) proves (8) directly. The 36 conditional boards give an independent
literal determinant/indicator check. This is the smallest size allowing
two **stochastic** components after palette conditioning: both must have
length>=6. It is a specific minimum in that model, not a global minimality
claim about arbitrary cumulant examples.

## 5. A general all-four-cycle reduction

For n=2m and `G=m C4`, every conditional board is deterministic. Uniform
boards of this cycle type are therefore equivalent to the following exact
palette model: choose an unordered perfect matching of the row labels,
an independent perfect matching of the column labels, and a bijection
between their m pairs. Put a complete 2x2 rectangle on each matched pair.
The number of distinct boards, each equiprobable, is

```text
((n)!/(2^m*m!))^2 * m! = (n!)^2/(2^(2m)*m!).                    (9)
```

Every board has a unique such triple of choices, proving the count and
uniformity. Equivalently, named palette assignments represent each board
m! times, and each has `2^(2m)` internal labelings. Thus every ordinary
cumulant of order>=2 in this entire all-size sector comes from the palette
distribution. Formula (2) is pointwise X here because `q_(2,1)=q_(2,2)=1`.

Every nonaxis line triple has component occupancy `(2,1)` or `(1,1,1)`;
it cannot have occupancy `(3)`. The two points in a `(2,1)` event must be
opposite corners of their rectangle. Consequently X splits exactly into
the sum, over ordered distinct blocks A,B and each of A's two diagonals,
of points of B on that diagonal line, and the sum over three distinct
blocks of collinear
choices of one point per block. This is an inter-block geometric statistic
on random paired labels. Its multiplication still has the usual label
exclusion and shared-block information; (9) supplies no sign or tail bound.
A point at the intersection of A's diagonals contributes twice to the
first sum: it produces two different triples. Merging the two diagonal
lines into an unweighted union would lose that multiplicity.

## 6. A positive algebra that actually retains the missing labels

There is a small exact way to state what changes in the averaging map.
For injective partial assignments `alpha:A subset left(G)->[n]` and
`beta:B subset right(G)->[n]`, let `u_(alpha,beta)` be the indicator that
the full random labels extend both assignments. Their product is zero if
either union conflicts on a vertex or repeats a label on different
vertices of the same shore; otherwise it is the cylinder of the union.

The state is

```text
phi_n(u_(alpha,beta)) = 1 / ((n)_(|A|)*(n)_(|B|)).                (10)
```

This is an actual positive state: `phi_n(Y^2)>=0` for every real linear
combination Y of cylinders. To see the whole algebra explicitly, impose
the extension relations that a partial cylinder is the sum of its one-
vertex extensions over unused labels. Iteration expresses it as a sum of
full-assignment cylinders, which are pairwise orthogonal idempotents and
sum to one. Thus (10) is precisely the uniform state on `(S_n)^2`.

Conditional on P, a cylinder incompatible with its component palettes has
state zero; a compatible cylinder has state

```text
product_i 1 / ((k_i)_(|A intersect left(G_i)|)
               *(k_i)_(|B intersect right(G_i)|)).                (11)
```

Hence the genuine state is a finite mixture of the product states (11).
Geometry is carried by summing cylinders whose label coordinates satisfy
the required determinants; products use actual union and collision rules.
This realizes (1)-(5) in a positive algebra. It does **not** assign the old
formal path-defect D or homogeneous Q an unjustified random-observable
meaning. Disjoint union of unlabelled types is a different multiplication.

Even disjoint domains do not make unconditional cylinders independent:
two different left vertices both demanding label zero have separate
expectations 1/n but product zero. This is the elementary lost coordinate
behind the global falling-factorial denominator. Component palettes move
that exclusion into an explicit outer random variable; they do not erase it.

## 7. Verification, maps, and stopping boundary

The standalone companion `overnight6_20260906_no3line_palettes.py` imports
no repository modules or event census. It checks all 34 one-shore partial
injections at n=3 and every product on full permutations; positive two-shore
squares; local matching kernels through three edges; every n=4/n=5 palette
mean and singleton line-event kernel; bounded complete union controls;
all 576 and 14,400 independent shore-label assignments; (6)-(8); and the
all-C4 counting identity. Its geometry paths are literal integer
determinants and canonical pair-line counts. The historical n=4 output is
used only as a comparison value, not as computational input.

```text
python -B 04-computation/overnight6_20260906_no3line_palettes.py
python -B -O 04-computation/overnight6_20260906_no3line_palettes.py
```

Normal and optimized LF outputs agree: **15,676** active exact gates.
No third-event census is repeated or enlarged. The local k<=7 checks do
not replace the all-size conditioning proof.

```text
source SHA-256 cc0fc2774d8475766da8f90ab4a5bb12d070170bfc5ba2b4d3bbf0c74e05f80f
output SHA-256 006b6f15076863925e168088e7aae264f68893b34231403548772cc07a5cd5b9
semantic SHA-256 e0c664137c0c816a91d2946e2e8a11826918c0b0f4ddf97297e0038bf42c16fe
```

| Source -> target | Preserved predicate | Lost information / sidecar / decisive check |
|---|---|---|
| Full labels -> named component palettes | Which coordinate labels belong to each cycle | Internal bijections; retain their local injection kernels (1) |
| Event family -> complete local unions | Simultaneous presence of every event | Event list is restored for cumulants; use all partition terms and repetitions |
| Conditional kernels -> scalar unconditional cumulant | Exact probability after total cumulance | Averaging conditional cumulants alone loses outer palette terms; (6),(7) refute it |
| Formal type ring -> scalar geometric weight | Whole union-type incidence counts | Multiplication/collision data; actual cylinders (10) repair positivity |
| All C4 skeletons -> random paired palettes | The entire board, including Euclidean coordinates | No loss; (9) is an exact alternative carrier |

The new stopping point is concrete: residual internal-cycle connectedness
is now rigorous, and the required outer coordinate is explicitly finite.
The global palette distribution still couples geometric events, including
events spanning cycles. No component-additive unconditional formula,
all-n sign theorem, convergence bound, or positive probability of X=0
follows merely from this conditioning. The original LRC/no-three-line
existence problems remain open at their existing scope.

**Filing:** root integrated this independently audited report after
`f5f0f7f75`. The proof and transcript are frozen in the sixth manifest.
