# Session128: Pisot substitution attractors and exact feedback transport

**Status: PROVED scoped constructions; FINITE-EXACT certificates; CITED context.**
This is a proposed connection of the Session127 arithmetic source to a standard
substitution construction. It is not attributed to Amanda's Loop, does not
settle a Pisot spectral conjecture, and does not claim research priority for
classical Pisot/Rauzy machinery.

## Results first

1. The k-bonacci field, ordered substitution, and contracting conjugate
   representation fit an explicit commuting diagram, with an integral
   unimodular change of basis. The field degree is not discarded.
2. The Tribonacci prefix graph has a unique compact internal-space attractor;
   its contraction factor is strictly below `3/4`. A depth12 approximation has
   Hausdorff error below `531441/4194304` in the stated complex norm.
3. An arbitrary-size family of static linear feedback matrices is exactly
   absorbed by a single common linear shear. Its digit graph, all labeled
   compositions, global address-equality relation, and Fourier cocycle transport
   exactly. There is no perturbative small-gain hypothesis.
4. **Arithmetic rigidity:** every nonzero member of that shear family fails
   even to send the original embedded order into itself. Transporting the
   arithmetic module and metric is mandatory, not optional.
5. Identical incidence matrices and eigenvalues do not determine the ordered
   digit system. Standard and twisted Tribonacci give an exact hostile.
6. The previous sparse unit `beta^(-16)=56 beta−103` has internal modulus
   `beta^8>100`; it lies outside this attractor, whose radius is below4.

## 1. Inheritance, scope, and concept board

Repository snapshot: `03a9a2c8839e62880cfe84c74ccfb0a2413d2f8d`.
The exact inherited Session127 packet is
`05-knowledge/results/session127_gsma_norm_collapse_20261031/README.md`, with
its `gsma/INHERITANCE_AND_KERNEL.md` and `tribonacci/REPORT.md`.

- **Closest proved mechanism:** multiplication by a multinacci unit, exact
  sparse-unit coordinates, and the companion relation from `THM-227`,
  `01-canon/theorems/THM-227-knacci-mersenne.md`.
- **Canonical hostile:** `MISTAKE-028` in `01-canon/MISTAKES.md` rejects using
  recurrence/trace coincidences as tournament-spectrum theorems. This packet
  likewise distinguishes matrix eigenvalues from dynamical spectrum.
- **Corrected near miss:** Session127 replaces the historical floating
  continued-fraction tail by exact arithmetic; no floating tail is imported.
- **Least-used sidecar:** prefix order and translations, absent from the
  incidence matrix, are the missing coordinates recovered here.

The live board is (i) arithmetic unit/norm, (ii) ordered substitution,
(iii) contracting embedding, (iv) labeled affine feedback, (v) Fourier
cocycle. The anchor is the exact source-to-attractor adapter; the niche is
feedback compatibility; the wildcard is sparse-unit escape. None asserts
progress on LRC(14), which remains outside this packet.

Amanda's supplied paper, as independently audited in the sibling packet,
has a nonlinear threshold update for `b` and `sigma`, not this substitution
or a specified matrix feedback operator. Its positive scalar gain can
separately disappear under a sign observer. Neither that fact nor a
Tribonacci ratio limit identifies Amanda dynamics with the construction
below. The construction here is independent and new **to this session**.

## 2. A typed field–substitution–conjugate diagram

For any integer `k>=2`, define

    f_k(t)=t^k−sum_(j=0)^(k−1)t^j,
    beta = its positive real root,
    K=Q(beta), O=Z[beta].

The standard ordered substitution is

    sigma(j)=1(j+1)  (1<=j<k),     sigma(k)=1.

Our column-source incidence convention is `M_ij=#i in sigma(j)`. Thus M has
first row all1 and subdiagonal1. The power-basis multiplication matrix C for
beta has subdiagonal1 and last column all1. These matrices are not silently
identified.

Let

    ell_j(t)=t^(j−1)−sum_(r=0)^(j−2)t^r,
    eta(v)=sum_j v_j ell_j(beta),

where the sum is empty at j=1. The coefficient matrix U whose columns are
`ell_j` is upper triangular with diagonal1, hence `det U=1`. Directly,

    U M=C U,             eta(Mv)=beta eta(v).

Consequently `eta:Z^k -> O` is an additive lattice isomorphism that
intertwines these operators. The `ell_j(beta)` are positive: starting with
`ell_k(beta)=beta^(-1)>0`, use
`beta ell_j=1+ell_(j+1)` backwards. They are the natural tile lengths with
`ell_1=1`.

### Why these are Pisot units (classical elementary proof)

Descartes' rule and `f_k(1)<0<f_k(2)` give one positive root beta in `(1,2)`.
The telescoping identity is

    (t−1)f_k(t)=t^(k+1)−2t^k+1.

For every radius r sufficiently close to1 from above and below beta,
`r^k(2−r)>1`: its derivative at1 is `k−1>0`. On `|z|=r`, the term
`z^k(z−2)` has modulus at least `r^k(2−r)>1`. Rouché's theorem therefore
puts exactly k zeros of `(z−1)f_k(z)` inside that circle: the same number
as `z^k(z−2)`. Besides the zero1, these are all k−1 roots other than beta.
Letting r decrease to1 puts those roots in the closed unit disk. A root z
on the unit circle would satisfy `|2−z|=1`, forcing z=1; but `f_k(1)!=0`.
Thus all other roots have modulus strictly below1.

This also proves irreducibility over Q. If a monic integer factorization
existed, a nonconstant factor not containing beta would have all its roots
strictly inside the unit circle, so its nonzero integer constant term would
have absolute value below1, impossible. Hence `[K:Q]=k`; the constant term
also makes beta a unit. No maximal-order claim is needed.

### Where the compact attractor lives

Write the signature of K as `(r1,r2)`, retaining one from each complex pair,
and omit only the Perron real embedding. The internal space is

    H=R^(r1−1) x C^r2,           dim_R H=k−1.

Let `star:K -> H` be this tuple of embeddings. Let A multiply each
coordinate by the corresponding non-Perron conjugate. In the maximum of
absolute coordinate norms, `||A||=c<1`, and

    star(beta x)=A star(x),       star(eta(Mv))=A star(eta(v)).

This is a real-linear quotient of the full Minkowski representation after
extending scalars, not a real-linear equivalence: its expanding dimension
has been removed. On K itself each field embedding is injective; that fact
does not make the projected additive module a discrete lattice in H or
supply a continuous inverse. The compact attractor lives in the real
completion H; its translations come from `star(O)`. A limit point is not
asserted to belong to the countable field image.

The fields of degrees3,4,5 are related by this common construction, not by
field embeddings between them. Degree divisibility forbids such embeddings.

## 3. Ordered prefix graph and compact attractor

For each occurrence `sigma(j)=p i s`, give the edge `e:j->i` the digit
`t_e=star(eta(abelianization(p)))`. Retain edge identity and multiplicity.
For the standard substitution all digits are0 or `d=star(1)`:

    R_1=union_(j=1)^k A R_j,
    R_i=A R_(i−1)+d                 (2<=i<=k).

On the product of nonempty compact subsets of H, use maximum Hausdorff
distance. The displayed set operator contracts by c, so the contraction
mapping theorem gives a unique nonempty compact tuple `(R_i)` and
convergence from any nonempty compact starting tuple. This proves neither
pairwise disjoint interiors nor a translation tiling nor pure dynamical
spectrum. Those are distinct statements with further proof obligations.

For k=3 choose the conjugate alpha with positive imaginary part. Then
`H=C`, `A z=alpha z`, `d=1`, and

    R1=alpha(R1 union R2 union R3),
    R2=alpha R1+1,
    R3=alpha R2+1.

Equivalently,

    R1=alpha R1 union (alpha^2 R1+alpha)
                    union (alpha^3 R1+alpha^2+alpha).

The cubic discriminant is−44; its other roots are a complex pair. Vieta
gives `|alpha|^2=beta^(-1)`. Since `f_3(9/5)=−26/125` and f is increasing
above1, beta>9/5, so `|alpha|^2<5/9<9/16` and `|alpha|<3/4`.
Every path expansion is a convergent digit sum, hence
`|z|<=1/(1−|alpha|)<4`. From the starting sets `{0}`, the depth-n Hausdorff
error is at most `c^n/(1−c)<4(3/4)^n`.

**FINITE-EXACT:** `verify.py` compares every graph-path translation,
including its coefficient multiplicity and terminal/source labels, with
every ordered prefix in `sigma^n(j)` for all j and `0<=n<=12`. At depth12
the word lengths are1705,1431,927; component support counts are927,504,274.
No floating coordinates or plotted resemblance are used in this equality.

## 4. Exact nonperturbative linear-feedback absorption

Consider a fixed finite labeled graph with branch maps

    f_e(x)=A x+t_e,      e:j->i,
    g_e(y)=B y+t'_e.

A componentwise affine map `H_i(x)=S x+b_i`, with one invertible linear S,
intertwines **every labeled branch on its entire ambient space** iff

    B S=S A,
    t'_e=S t_e+b_i−B b_j                 for every edge e:j->i.       (F)

Proof: expand `H_i(f_e(x))=g_e(H_j(x))` and compare linear and constant
coefficients. These conditions are necessary for ambient branch conjugacy,
not claimed necessary for arbitrary homeomorphisms of only the attractor.
They transport all finite compositions and infinite path codings. If A
contracts in a norm, B contracts in `||y||'=||S^(-1)y||`; uniqueness gives
`R'_i=S R_i+b_i`. Injectivity preserves address coincidences within each
fixed target component, or in the disjoint labeled union. Different b_i
can change cross-component spatial intersections: for example, two
self-loop attractors initially both `{0}` can be translated to `{0}` and
`{1}`. When all b_i agree, as in the shear family below, global spatial
coincidences are preserved too. No claims
about a single-valued map on an overlapping attractor are needed.

A feedback law `y_next=A y+G y+t'_e` has `B=A+G`. Calling its feedback
"absorbed" requires (F), not just stability, matching eigenvalues, or a
formal resolvent. If an actuator matrix D is specified, a realizable gain
also needs the separate factorization `G=D L`. We use full-state additive
feedback and do not silently claim single-input realizability.

### An arbitrary-size positive family

Set `s=1−beta`, `p=beta^(-1)` and use real coordinates `z=x+y alpha`, so
`alpha^2=s alpha−p`. Then

    A = [[0,−p],[1,s]],       d=(1,0).
    S_h = [[1,h],[0,1]],     h any real number.
    B_h=S_h A S_h^(-1)
       = A + G_h,
    G_h=[[h,hs−h^2],[0,−h]].

Since `S_h d=d`, the same graph and same literal digit vectors0,d work
with `b_i=0`. Thus the feedback-modified compact tuple is exactly `S_h R_i`
for every h, with the same contraction factor in the transported norm.
Gains can be arbitrarily large. The finite checker verifies the algebraic
identities for `h=0,1,−2,10,7/3`, on all1705 distinct depth12 points; the
all-real-h assertion follows from the displayed symbolic calculation.

**Hostile:** replacing A by `A+gI` changes its trace by2g. Similarity is
therefore impossible unless g=0. Even a stable modified system need not be
conjugate to the original one. Same eigenvalues would still not replace the
edge-digit equations in (F).

### Arithmetic rigidity: only h=0 preserves the original order image

The abstract field is represented in two different embedded copies:
`K_beta=Q(beta)` is real, while `K_alpha=Q(alpha)` is not. The planar
coefficient scalars p,s belong to K_beta; they are not being identified
with their nonreal conjugates. Let `Lambda=Z+Z alpha+Z alpha^2`, a rank3
additive module, not a planar lattice, represented in `(1,alpha)` coordinates as

    Lambda={(a−cp,b+cs):a,b,c in Z}.

Since s is irrational, `Lambda intersect (R x {0})=Z x {0}`. If
`S_h Lambda subset Lambda`, subtract alpha from `S_h alpha` to get
`(h,0) in Lambda`, hence h is an integer. Likewise, subtract alpha^2 from
`S_h alpha^2` to get `(hs,0) in Lambda`, hence hs is an integer. Irrationality
of s forces h=0. Conversely h=0 plainly preserves Lambda.

The same argument with rational rather than integer coefficients proves
the stronger rigidity `S_h K_alpha subset K_alpha iff h=0`.

Consequently every nonzero shear needs the transported module `S_h Lambda`.
It is a faithful conjugate real representation of the arithmetic source,
not the original literal complex field embedding or order. The transported
metric is also essential; large shears can have Euclidean operator norm
above1 while their spectral radius remains below1. Explicitly, the Gram
matrix of the complex modulus in these coordinates is

    Q=[[1,s/2],[s/2,p]],       A^T Q A=p Q.
    Q_h=S_h^(-T) Q S_h^(-1),  B_h^T Q_h B_h=p Q_h.

Q is positive definite because alpha is nonreal; its determinant is
`(Im alpha)^2>0`. These equations give an exact quadratic Lyapunov metric
for every h. This rigidity argument
was supplied by the parent lane and independently checked here.

## 5. Fourier cocycle covariance: a genuine spectral-data statement

For any labeled digit graph, define

    H_ij(y)=sum_(e:j->i) exp(2 pi i <y,t_e>),
    H_n(y)=H(y) H(A^T y) ... H((A^T)^(n−1)y).

The matrix entries are finite Fourier sums of all length-n path digits,
with multiplicities. In particular `H(0)=M`; the value at0 forgets all
prefix positions. This matrix need not itself be a unitary dynamical
operator.

Under (F), put `D_b(y)=diag(exp(2 pi i <y,b_i>))`. Direct substitution gives

    H'(y)=D_b(y) H(S^T y) D_b(B^T y)^(-1).

Using `S^T B^T=A^T S^T`, intermediate diagonal factors cancel:

    H'_n(y)=D_b(y) H_n(S^T y) D_b((B^T)^n y)^(-1).

This is an exact all-depth Fourier-cocycle covariance law. For the shear
family b=0 it reduces to `H'_n(y)=H_n(S_h^T y)`. Although the one-step
literal digit matrix is unchanged, its frequency iteration changes with
B; the dual-coordinate transport must not be omitted. This proves a
spectral-data transport statement, not a purity theorem for any symbolic
or tiling dynamical spectrum.

## 6. Same matrix, different ordered digit law

The twisted substitution `(21,13,1)` shares M, f, beta and every matrix
trace with `(12,13,1)`. But its edge `1<-1` has digit `alpha−1` rather than0,
and its edge `2<-1` has digit0 rather than1; all other edge digits agree.

There is no ambient componentwise real-affine conjugacy preserving these
edge identities with a common linear part. Indeed both systems have A,
so (F) forces S to commute with complex multiplication alpha. A real
linear map commuting with a nonreal complex multiplication is itself
complex multiplication by a scalar c. The self-edge equation forces
`b1=−1`. The edges `1<-2` and `1<-3` force
`b2=b3=−1/alpha`. The remaining edges respectively force

    c=alpha^(-1)−alpha,       c=alpha^(-1),

a contradiction. This does not rule out differently labeled, non-affine,
or merely measurable identifications. It pinpoints the exact failure of
the attempted matrix-only equivalence.

## 7. The sparse-unit/Lucas coordinate junction and escape

There is a direct inherited coordinate connection, not merely a repeated
integer. For a sparse element `x=p−q beta`, its substitution-basis vector is

    U^(-1)(p,−q,0,...,0)^T=(p−q,−q,0,...,0)^T.

Therefore the Session127 Lucas47 shear `r=p−q` is literally the first
substitution coordinate. The inverse unit beta^(-16) has coordinates
`(-47,56,0)`, while its sign-normalized negative `103−56 beta` has
coordinates `(47,−56,0)`. This identifies two coordinate descriptions of
the same source element. It does not make the vector a nonnegative letter
count, legal prefix, or attractor address. The triadic lane recovered this
connection; the present checker verifies the coordinates exactly.


Session127's exact identity transports under the field embedding to

    alpha^(-16)=56 alpha−103.

Its internal modulus is

    |alpha^(-16)|=beta^8>(9/5)^8=43046721/390625>100.

It therefore lies outside every R_i, all of which lie in `|z|<4`.
The Perron coordinate tends to zero for inverse unit powers while the
contracting-conjugate coordinates of those **inverse** powers expand.
Norm±1, sparse coordinates and a large continued-fraction digit are not
attractor-membership conditions. Any successful bridge needs an explicit
normalization/address operation and must show what it preserves.

## 8. Primary-source context and boundaries

The construction belongs to established Pisot substitution/Rauzy-fractal
machinery. Baake–Grimm's [*Fourier transform of Rauzy fractals and point
spectrum of 1D Pisot inflation tilings*](https://arxiv.org/abs/1907.11012),
v2, §§3,4,6, treats internal-space contraction, Fourier cocycles and the
Pisa family. Their §6.2 compares standard and twisted Tribonacci with the
same incidence matrix but different windows and diffraction intensities.
Our normalization has `ell_1=1`; theirs uses `ell_1=beta`, so geometric
coordinates are related by multiplication by the corresponding conjugate.
Thuswaldner's [*Unimodular Pisot substitutions and their associated
tiles*](https://www.numdam.org/articles/10.5802/jtnb.556/) distinguishes
contact-graph boundary constructions and an additional coincidence
condition. These sources provide context, not an imported proof of our
feedback theorem or an audited statement of every current Pisot conjecture.

No conclusion here derives pure spectrum, a translation tiling, measurable
conjugacy of Amanda's threshold loop, or equivalences between degrees3,4,5
from the Pisot property alone.

## Reproduce

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json

The program is dependency-free. Explicit exceptions retain every check in
optimized mode. It records the exact finite universe, preserved labels,
positive shear controls, scalar-gain hostile, twisted digit hostile, and
the escaping inverse unit. `SOURCE_PINS.json` records inherited file
hashes; `sources.json` records primary citations; `SHA256SUMS` fixes this
packet's files. General mathematical claims rely on the proofs above, not
on finite depth12 or the five sampled gains.
