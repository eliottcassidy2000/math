# Node-coloured defects, needles, Radon data, and the typed metagraph

**Status:** exact synthesis of proved identities and finite audits, followed by
explicit conjecture targets. This is not a proof of LRC(14). The standard
fourteen-runner case remains open.

The recurring failure in this part of the repository is not a shortage of
statistics. It is a type error. A tournament class, a fixed-path tiling, a
complement line, a factorization witness, a path realizer, and an LRC safe
component answer different questions. Equal cardinalities between two of
these sets do not supply a map, and a map that is adequate for one operation
need not be a congruence for the next operation.

The useful replacement is a **typed, operation-indexed metagraph**. It retains
literal witnesses until a minimization theorem proves that a field is
redundant. The node-coloured defect algebra, the Hamiltonian-path Walsh stalk,
the path-holonomy law, and the finite-field Radon proposal are four shadows of
that one rule.

## 1. The object has several sorts

At size `n`, put

```text
S_n = {(a,b): 1 <= b < a <= n, a-b >= 2},
m_n = |S_n| = binom(n-1,2),
X_n = F_2^(S_n).
```

The fixed path is `n -> ... -> 1`. A point of `X_n` is a literal marked-path
tiling and hence a labelled tournament with that Hamiltonian path. The maps
already proved in THM-781, THM-796, and THM-830 form the span

```text
                         q_n
                 X_n ----------> M_n
                  |
                  | orbit by kappa
                  v
                 L_n = X_n/<kappa> ----> Sym^2(M_n).
```

Here `M_n` is the converse-merged tournament-node set and `L_n` is the set of
literal complement lines. The endpoint map records a node loop or an unordered
pair of nodes. It also carries a blue/black colour. These are three sorts, not
three names for one graph.

THM-830 makes the incidence numerical. If `Phi_n(G)` is the merged node of the
gap tournament `G`, then

```text
K_n^c(u,v) = #{G : Phi_n(G)=u, Phi_n(G^op)=v, colour(G)=c}.       (1.1)
```

For `u != v`, this is the literal multiplicity of the coloured node edge. For
`u=v`, the loop multiplicity is `K_n^c(u,u)/2`. Its row sums recover the
coloured tiling fibres. Thus node adjacency, tiling incidence, and line
multiplicity are all pushforwards of one literal kernel, but none determines
the other two without the fibre weights.

The exact weighted handshake dictionary is

```text
N_n^c(u)=sum_v K_n^c(u,v),
ell_n^c(u,u)=K_n^c(u,u)/2,
ell_n^c(u,v)=K_n^c(u,v)=K_n^c(v,u),              u!=v,
|L_n^c|=(1/2)sum_(u,v)K_n^c(u,v).                         (1.2)
```

Thus a node loop consumes two tilings in the same fibre, whereas a cross-node
line consumes one tiling at each endpoint. Forgetting multiplicity converts
this weighted multigraph to a simple graph and destroys both fibre size and
the number of literal continuation witnesses.

There are at least four further witness sorts:

```text
G_n       defect-groupoid arrows (x,delta),
F_n       two-arrow factorizations g = g_2 o g_1,
R_n       direct self-line realizers (T,P,sigma),
O_n       owner-labelled LRC wall, gap, and proof-obligation packets.
```

The first three are purely tournament-theoretic. The fourth is where the LRC
metric lives. No result presently gives a lossless map from `O_n` to a bare
tournament node.

## 2. The node-coloured defect algebra

Staircase reflection has `f` fixed tile coordinates and `r` transposed pairs.
Write

```text
h=f+r,                 m=h+r,
x in F_2^h,            delta in F_2^r.
```

The arrow `(x,delta)` has target `x+i(delta)`, where `i` inserts the defect in
the paired coordinates, and composition is

```text
(x+i(delta_1),delta_2) o (x,delta_1)
    = (x,delta_1+delta_2).                                  (2.1)
```

This is an elementary abelian pair groupoid. Reflection is arrow inversion.
Blue arrows have `delta=0`; black arrows have `delta != 0`. A seed colour may
retain the full defect and the merged node of the decoded endpoint. It is
already more informative than a node colour or a blue/black bit.

For an arrow colour `c`, define the ordered factor tensor

```text
p_ab(g) = #{(g_1,g_2): g=g_2 o g_1, c(g_1)=a, c(g_2)=b}.        (2.2)
```

The recursive refinement

```text
c_(s+1)(g) = (c_s(g), (p_ab^s(g))_(a,b))                       (2.3)
```

is the coarsest refinement of the seed for which all composition structure
constants are cellwise constant. Before quotienting by inversion, the Jordan
version retains `p_ab+p_ba`. Its symmetric-closure interpretation requires
the input colour to be inversion-invariant; omitting this hypothesis would be
false in general.

At the binary blue/black level those constants have a closed all-size form.
With `B` denoting defect zero and `K` a nonzero defect, the factor tensor is

```text
output B:  p_BB=1, p_BK=p_KB=0, p_KK=2^r-1,
output K:  p_BB=0, p_BK=p_KB=1, p_KK=2^r-2.               (2.3a)
```

This follows by taking the first factor defect `alpha`: the second is
`alpha+gamma`. It is an exact equitable two-colour algebra, but precisely
because the row depends only on whether `gamma` vanishes, it cannot distinguish
two blue outputs or two black outputs. The useful refinement is not more
binary counting; it is the attachment of the exact defect to both factor-node
incidences.

The size evolution makes the same information loss quantitative. Writing
`m_n=h_n+r_n` as above,

```text
m_(n+1)-m_n=n-1,
h_(n+1)-h_n=floor(n/2),
r_(n+1)-r_n=floor((n-1)/2).                               (2.3b)
```

There are `2^(m_n-1)` literal complement lines and `2^(h_n-1)` blue lines,
so the blue fraction is exactly `2^(-r_n)`. Every size step adds
`floor((n-1)/2)` new defect directions; the binary colour collapses their
entire nonzero space to one symbol, while the factor deck retains it as the
new recursive fibre.

THM-851 gives the exact first-step census:

```text
n                 3    4     5      6       7
arrows            2    8    64   1,024  32,768
seed cells         2    4    23     255   7,926
ordered cells      2    6    64   1,024  32,768
Jordan cells       2    5    40     544  16,640
reflection orbits  2    6    40     544  16,640.              (2.4)
```

Thus ordered factor data reconstructs every arrow at `n=3` and `5<=n<=7`,
while Jordan data reconstructs every reflection orbit there. Both fail in the
small exceptional chart `n=4`. This exception is useful: even an exact algebra
can have a low-size accidental kernel, so an all-size claim needs a proof, not
an interpolated table.

The crucial information is **which defect carries which factor pair**. If the
midpoint defect labels are erased and only the unlabelled multiset of
intermediate node pairs is retained, the deck misses `1,1,36` distinctions at
`n=5,6,7`. A node-coloured algebra is therefore a labelled incidence tensor,
not merely a histogram of adjacent node colours.

For fixed output `(x,gamma)`, the full deck can be written as the function

```text
alpha |-> {C(x,alpha), C(x+i(alpha),alpha+gamma)},
alpha in F_2^r.                                             (2.5)
```

Because `C` retains the exact defects, the factor colours recover the index
`alpha`; the multiset presentation and this indexed function are equivalent.
Its Walsh transform in `alpha` is therefore lossless. The unlabelled midpoint
deck is its node-pair population shadow. A sharp compression problem is to
find the least Walsh sector of (2.5) that reconstructs reflection orbits and
remains closed under A/B/C and the THM-853 return semigroup.

## 3. The H drift is a second operation algebra

Let `H(T)` be the number of directed Hamiltonian paths and let `T^e` denote
one arc flip. Put `M=binom(n,2)` and

```text
b_H(T) = (1/M) sum_e (H(T^e)-H(T)).                          (3.1)
```

THM-848 decomposes `H=sum_r H_(2r)` into even Walsh levels supported on
vertex-disjoint unions of positive even paths. The exact drift is

```text
b_H = -(4/M) sum_r r H_(2r).                                (3.2)
```

If `U_1` counts the Hamiltonian needles with one boundary defect and `K` is
the odd-cycle-forest partition function, then

```text
sum_e Delta_e H = U_1-(n-1)H,
U_1 = 2K-(n+1)H,
K = nH + (M/2)b_H.                                         (3.3)
```

So the H-drift functional-form question has a sharp answer: one step needs
the degree-Euler image `K=(n-D)H`, equivalently the length-weighted Walsh
stalk. `H` alone is insufficient. Iterated means require the whole Krylov
stalk generated by the Walsh levels; even the full first-moment stalk need
not determine diffusion.

The live THM-848 addendum now computes the next Krylov step. If `a_j(T)` counts
vertex orders with exactly `j` forward consecutive arcs, put
`A_T(x)=sum_j a_jx^j` and let `S` be the unnormalized sum over all arc flips.
Then

```text
S a_j=(j+1)a_(j+1)+(n-j)a_(j-1)-(n-1)a_j,
S A_T=(1-x)[(1+x)A_T'-(n-1)A_T].                          (3.3a)
```

Writing `b_r=a_(n-1-r)` gives

```text
K=(b_1+(n+1)b_0)/2,             S K=b_2+b_1-Mb_0.          (3.3b)
```

Thus `K` is the second Krylov needle, and its drift requires the next top
coefficient. Direction by direction there is an exact odd-cycle-forest current:

```text
K(T^e)-K(T)=4[sum_(C in C_e^+)K(T-V(C))
                  -sum_(C in C_e^-)K(T-V(C))].             (3.3c)
```

The pair `(H,K)` determines the averaged `K` drift on every merged node through
`n=7`, but exactly one of its 1,727 fibres splits at `n=8`; the split is the
first realized ambient level-kernel direction. The unordered `K`-gradient
multiset separates all 272 size-seven nodes. This is the same lesson as the
Kakeya covariance algebra below: a low-order moment pair can close accidentally
on a finite locus, while the operation-stable object is the full tridiagonal
coefficient/Walsh stalk. There is no all-size two-scalar `(H,K)` closure.

The energy of level `2r` is

```text
E_(2r)/mu_n^2
 = [w^r] ((1+w)/(1-w))^(n-2r) / (n)_(2r),
mu_n = n!/2^(n-1).                                         (3.4)
```

At `n=8`, (3.4) gives `E_6/mu_8^2=1/1680` and total variance ratio
`131/560`; the former formula first fails there. This is another next-size
guardrail against pattern extrapolation.

The multiset

```text
Grad_H(T) = multiset_e (H(T^e)-H(T))                        (3.5)
```

classifies every converse-merged node through `n=8` in the exact atlases:
`3,10,34,272,3528` cells at `n=4,...,8`. This is finite evidence, not an
all-size theorem. It is nevertheless a strong candidate node stalk because it
records all one-step directional responses without naming arbitrary labels.

Equations (2.2) and (3.5) are the same recursive idea at different operation
alphabets: retain the multiset of coloured outcomes of every legal local
factor or flip. Their common abstraction is the operation profile

```text
Prof_O(x) = multiset_(o in O) (o, colour(o x)).             (3.6)
```

## 4. Self-lines require path holonomy

Let `kappa_P` flip all pairs off the marked Hamiltonian path `P`. The black
quasi-fixed endpoint set is stable under the Klein four generated by full
flip `kappa` and staircase reflection `sigma`, and its orbits are free. Hence

```text
4 divides Q_black(n),             2 divides selfK(n).       (4.1)
```

The finite endpoint counts are

```text
n                    5    6    7     8
Q_black(n)           8   12   88   404
selfK(n)             4    6   44   202
Klein-four orbits    2    3   22   101
SC(n)                8   12   88   176.                    (4.2)
```

Therefore `Q_black(n)=SC(n)` at `n=5,6,7` is a small-size coincidence, not an
all-size self-line law. There can be no all-size two-to-one bijection with
self-converse classes in the proposed form.

The surviving all-size law is witness-level. If `pi` realizes
`pi.T=kappa_P T`, and `L_pi` is the linear relabelling action on pair masks,
then

```text
pi^k.T = T + [k odd]1 + sum_(j<k) L_pi^j p_0,              (4.3)
N_(L_pi) p_0 = [ord(pi) odd]1.                             (4.4)
```

Thus odd witness order requires a mod-2 decomposition of `K_n` into the
translates of one Hamiltonian path, and is possible only for
`n=1,2 (mod 4)`.

Most importantly, relabelling does not commute with a fixed path flip:

```text
sigma kappa_P = kappa_(sigma P) sigma,
sigma^2 T = Flip_(E(P) triangle E(sigma P)) T.             (4.5)
```

The square is an automorphism only when the two path edge sets agree. Exact
search finds no such direct black realizer for `n=5,...,8`. The quotient node
forgets precisely the obstruction it was asked to decide. The minimal stalk
is the marked pair `(P,sigma P)`, or at least its symmetric-difference mask,
not just the permutation cycle type.

In general the realizer-pair mass is

```text
sum_(t in Q_n^K) |Aut(T_t)|,                               (4.6)
```

not `2 selfK(n)`. Equality in the audited range comes from asymmetric
endpoints. This is the correct groupoid weighting behind any future counting
identity.

## 5. The three signed recursions are charts, not scalar laws

The full staircase has the exact `B3` Cech cover

```text
A+B+C-D-E-F+G.                                             (5.1)
```

THM-801 proves literal gluing, while THM-830 identifies every face as a gap-
tournament operation. With address variables `X,Y,Z`, the typed enumerator is

```text
T_n=h_(n-3)(X,Y,Z),
T_n=(X+Y+Z)T_(n-1)
    -(XY+XZ+YZ)T_(n-2)+XYZ T_(n-3).                        (5.2)
```

Setting `X=Y=Z=1` retains the cardinality and destroys the direction that
owned each coordinate. That direction is exactly what deletion, gap
contraction, and continuation use.

There is an exact invariant-ring form of this loss. Reflection exchanges
`X,Y`, and

```text
Q[X,Y,Z] = Q[s=X+Y,p=XY,Z]
           direct-sum (X-Y) Q[s,p,Z].                       (5.3)
```

In homogeneous degree `n-3`, the two summands have dimensions `h_n` and
`h_(n-1)`, so `m_n=h_n+h_(n-1)`. The half chart is the Reynolds-invariant
summand. What it destroys is the entire anti-invariant converse-current
module, not one unstructured sign bit.

The even half chart is the `B2` fold

```text
A+B-C.                                                     (5.4)
```

Its scalar recurrence is valid, but the fixed reflection line, shortcut seam,
and endpoint coupling are sidecars. THM-830 shows that the middle `B` term is
a tournament on cuts rather than deletion of one fixed original vertex.

The prompt-order odd word is

```text
A+B-C+D-E-F+G.                                             (5.5)
```

Its actual Venn chart has corners `A,D,B`, edges `A+B-C`, `A+D-E`, and
`B+D-F`, and center

```text
A+B+D-C-E-F+G.                                             (5.6)
```

The slots `C` and `D` have the same scalar size in the relevant count, so
they cancel after forgetting their roles. They are geometrically different:
one is an overlap and one is a corner. Any recursive carrier that cancels
them before transport is too small.

Algebraically, in odd homogeneous degree the invariant piece is covered by
the ideals `(s)` and `(Z)`, yielding the two-generator fold. In even degree,
the pure monomial `p^((n-3)/2)` is not covered and forces the extra corner.
This is the invariant-ring reason for the parity change and for retaining the
equal-size `C,D` roles separately.

The natural arithmetic companion is not another scalar sequence but the
product ledger

```text
Div(D) x B_r.                                              (5.7)
```

The divisor coordinate retains primitive period capacity `phi(d)` and CRT
labels; the Boolean coordinate retains which face, overlap, or far packet
owns the term. Mobius inversion is legal on each typed axis. Multiplying the
two only after their labels are retained prevents a cancellation in one axis
from erasing an obstruction in the other.

There is a sharper common algebra behind (5.1), (5.4), and the Hunter moments.
For three Boolean observables `1_A,1_B,1_C`, the seven signed slots in (5.1)
are exactly the nonempty raw moments

```text
m_S=integral product_(i in S)1_(A_i),
mu(A union B union C)=sum_(empty!=S subset {A,B,C})
                         (-1)^(|S|+1)m_S.                  (5.8)
```

This is Mobius inversion on the Boolean subset lattice. The `A+B-C` chart is
its two-generator restriction. The prompt-order word (5.5) is not a different
scalar identity: it is a role-sensitive traversal of the same degree-three
moment cube, in which equal numerical coefficients do not license swapping an
overlap with a corner.

Centering is an invertible triangular change of coordinates only while all
lower moments remain attached. If `p_i=mu(A_i)` and `F_i=1_(A_i)-p_i`, then

```text
integral F_1 F_2 F_3
 =m_123-p_1m_23-p_2m_13-p_3m_12+2p_1p_2p_3.              (5.9)
```

Thus the Hunter edge stalk `theta_(Eij)` is the centered image of the `G`
corner of this same `B3` algebra. Erasing `G`, erasing the A/B/C carry
contrast, and erasing `theta_(Eij)` are three versions of one operation: a
degree-three Boolean moment has been projected away before the next recursive
action. More generally, an operation-stable carrier should be graded by typed
joint moments and should truncate only after a closure theorem shows that the
next degree is determined. The n=8 `(H,K)` split is the arc-cube analogue of
the same failure of premature truncation.

## 6. Three meanings of a needle

The repository now has three exact needle constructions. They should be
related by their incidence algebra, not identified literally.

1. A **Hamiltonian needle** is a directed vertex ordering in the tournament
   arc cube. Its Walsh boundary is an even path forest, and arc flips give
   the gradient stalk (3.5).
2. A **Farey needle** is a threshold interval moving through a rational gap.
   THM-841 shows the violation ladder has exact breakpoint profiles but not a
   homogeneous toothpick recursion.
3. An **affine finite-field needle** is a line of `F_7^2`. Its line sums give
   a discrete Radon transform.

For `v` in the eight directions `P^1(F_7)`, define

```text
R_v f(s) = sum_(x: v dot x=s) f(x),        s in F_7.       (6.1)
```

For any `d` distinct directions, over characteristic zero,

```text
rank(R_D)=1+6d,                 dim ker(R_D)=48-6d.        (6.2)
```

Proof: Fourier transform in `s` reveals `f-hat(a v)` for `a in F_7`. Every
direction supplies its six nonzero frequencies, while all directions share
the zero frequency. The eight projective frequency lines partition the 48
nonzero frequencies. Hence all eight directions are necessary and sufficient
to reconstruct an arbitrary scalar function on `F_7^2`.

The integer address simplex has an additional odometer fibre over this plane.
For

```text
S_r={(alpha,beta,gamma)>=0: alpha+beta+gamma=r},
```

choose canonical residues `xi=(a,b,c) in {0,...,6}^3` and put
`L_xi=(r-a-b-c)/7`. Then

```text
S_r ~= disjoint-union_(xi in H_(r mod 7)) Comp_3(L_xi),
(alpha,beta,gamma)=xi+7(i,j,k).                            (6.3)
```

An A/B/C increment advances the corresponding residue odometer; a `6->0`
wrap advances its carry coordinate. Thus the carry fibre recursively carries
the same B3 face structure.

At the `n=14` depth `r=11`, 33 residue fibres are singletons, 15 are triples,
and `(6,6,6)` is absent. Scalar residue pushforward has rank 48 and kernel
dimension

```text
30=15(3-1).                                                (6.4)
```

The full eight-direction Radon map has rank 49 on the abstract plane but only
rank 48 after composition with the 78-address pushforward. It reconstructs
residue sums, not the two independent A/B/C carry contrasts inside each
triple fibre. The three direction orbits under face triality have sizes
`3+3+2`; the chi charge lies in the two-direction orbit.

With no-carry and A/B/C carry treated as four channels, one common
reflection-invariant six-pencil deck has ranks `(33,15,15,15)` and total rank
78. Five pencils cannot reconstruct the 33-point no-carry channel. The exact
audit finds this full-rank six-set unique under true endpoint reflection; the
formal-triality six-set has total rank 77. Thus full triality and minimal
tomography are competing gauges, not interchangeable symmetries.

This supplies the precise Kakeya warning. One line in every direction is a
geometric hitting object; all parallel line sums in every direction are a
tomographic reconstruction object. A single chi-seven pencil has rank seven
and a 42-dimensional kernel. Seven directions still leave a six-dimensional
kernel. A proof that samples one attractive pencil cannot claim to have
preserved the plane.

The j=4 flood audit provides a structured, non-generic function. Its 21 floods
are indexed by the edges of `K_7`; Fano triples organize those edges, and the
negative masks form one Fano line. But only the identity in `GL(3,2)` preserves
the full exact numerical flood data. The Fano plane is therefore an address
system, not yet a numerical quotient. The one-sided endpoint needles visit at
most four masks and at most two masks on the negative line, so they do not
close a local chi-seven sweep.

The flood observables also occupy different edge modules. For edge `{a,b}`,

```text
r_(a,b)=x_a+x_b,       x=(16,16,14,14,10,14,6),            (6.6)
```

so `r` has orbit-span rank 7. The exact `m` and `V1` fields each have full
orbit-span rank 21. Point-star and Fano-triangle rows together span only 13
dimensions, leaving an eight-dimensional invisible edge field, and both `m`
and `V1` have nonzero curl there. A Fano transform may compress `r`; it cannot
compress the other two without retaining their residual edge cocycle.

The seven-comb Hunter wall supplies a fourth needle algebra, and it makes the
same preservation failure nonlinear. Let `E` be the prefix safe set,
`e=mu(E)`, and let `x_1,...,x_m` be its remaining danger combs. Define vertex
and edge anomalies

```text
s_i = mu(E intersect D_(x_i))-2e/13,
x_i=g_ij a_ij, x_j=g_ij b_ij,  gcd(a_ij,b_ij)=1,
h_ij=[Q(a_ij+b_ij)-Q(b_ij-a_ij)]/(169a_ijb_ij),
eta_ij=mu(E intersect D_(x_i) intersect D_(x_j))
       -e mu(D_(a_ij) intersect D_(b_ij)),
c_ij=e h_ij+eta_ij.                                      (6.7)
```

The corrected pair law gives `|eta_ij|<=2c_E/g_ij`. Hunter--Kounias then
becomes the exact defect identity

```text
L_H=(165-22m)e/169-sum_i s_i+MST(c),
MST(c)=max_T sum_(ij in T)c_ij.                           (6.8)
```

At `m=7`, the baseline term is the positive margin `11e/169`. The residual
state is therefore a seven-node coloured graph: node colours are restricted
single-comb anomalies, while edge colours retain reduced projective ratio,
common scale, the mod-13 sawtooth, endpoint discrepancy, and incidence. Raw
speed, raw difference, and pair-density histograms all forget part of this
state.

The maximum-tree evaluator is the tropical Kirchhoff polynomial. If
`lambda_1>...>lambda_q` are the distinct edge credits, `F_l` contains the
edges of credit at least `lambda_l`, and `r_l=m-kappa(F_l)`, then Kruskal gives

```text
MST(c)=sum_l lambda_l(r_l-r_(l-1)).                       (6.9)
```

Thus only connectivity-changing levels survive this particular evaluator,
but their incidence cannot be abelianized to a multiset. This is the same
structural lesson as the node-coloured factor deck: factor defects must remain
coupled to factor-node pairs; H derivatives must remain coupled to targets;
Hunter pair defects must remain coupled to the two comb vertices and the
prefix endpoints. In each case, a marginal contains every scalar value and
still loses the operation.

The node and edge defects are, more precisely, moments of one centered
indicator algebra. Put `p=2/13`, `F_E=1_E-e`, and
`F_i=1_(D_(x_i))-p`. Then

```text
s_i=<F_E,F_i>,                    h_ij=<F_i,F_j>,
theta_(Eij)=integral F_E F_i F_j,
c_ij=e h_ij+p(s_i+s_j)+theta_(Eij).                         (6.9a)
```

The completed restricted edge matrix is a Gram matrix:

```text
G^E_ij=c_ij-p(s_i+s_j),                         i!=j,
G^E_ii=22e/169+9s_i/13.                                   (6.9b)
```

Therefore it is positive semidefinite, coupling every node colour to every
incident edge through all principal-minor inequalities. The previously named
endpoint error is
`eta_ij=p(s_i+s_j)+theta_(Eij)`: its irreducible part is a third centered
moment. This is the exact sense in which the prefix supplies a new stalk. It
also aligns the Kakeya algebra with the tournament Walsh viewpoint: both are
graded moment algebras of centered Boolean observables, although one lives on
the circle and the other on the arc cube.

The projective part nevertheless has an exact compression. On `Z/13`, its
numerator kernel is

```text
H_(r,s)=Q(r+s)-Q(r-s)=C_Q(R-I).                             (6.9c)
```

It annihilates the seven-dimensional reflection-even sector and is positive
definite on the six-dimensional odd sector. In the basis
`u_i=e_i-e_(-i)`, its matrix is

```text
u_i^T H u_j=8 min(i,j)(13-2 max(i,j)).                      (6.9d)
```

Thus `H` is a rank-six positive-semidefinite Gram kernel. For a common-scale,
pairwise-coprime packet the global `h_ij` values are inner products of six
residue sine channels. Pairwise gcd changes break that direct common
embedding, and the endpoint terms `eta_ij` remain outside it. The correct
split is a low-rank reflection-odd residue current plus a prefix-coupled
metric correlation.

For arbitrary gcd patterns the `L^2` covariance Gram from (6.9a) survives,
with diagonal `22/169`, though its rank can grow. The exact scalar floor is
`h_ij>=-11/1014`, attained at reduced ratio `1:12`. Consequently every
seven-vertex tree has projective sum at least `-11/169`, exactly the negative
of the ideal Hunter margin if its edges are treated independently.

Complete-graph ratio consistency supplies the missing strictness. At the
threshold `lambda_0=-3/676`, the edges with `h_ij<lambda_0` have only the ten
reduced ratio types

```text
{1:9,2:9,4:9,1:10,3:10,1:11,2:11,3:11,1:12,1:25}.
```

Including reciprocals, this ratio set has no multiplicative triangle. Hence
the bad-edge graph is triangle-free and its good complement has at most two
components. A good spanning tree gives six edges at least `lambda_0`; in the
two-component case, five good edges plus a cross edge suffice. The unique
minimum edge is `-11/1014` at `1:12`, whose ratio graph has maximum degree two,
so a seven-vertex cross cut cannot consist entirely of minimum edges. The
second edge level is `-18/1859`. Therefore, for every seven-speed packet,

```text
MST(h)>=-237/7436,
11/169+MST(h)>=19/572.                                    (6.9e)
```

This is a strict projective theorem for arbitrary gcd patterns. It does not
make every edge or even the maximum tree positive: all 21 projective edges of
`{4,9,21,32,70,170,189}` are negative, with exact maximum-tree value
`-505/1447992`.

An exact low packet found by search has a more suggestive shape:

```text
{5,8,22,36,64,176,288}={5} union A union 8A,
A={8,22,36},
MST(h)=-941/334620,       Kruskal rank word=(1,3,1,1).      (6.9e1)
```

Its value and two maximizing trees are certified, but global minimality is
not claimed. The form is nevertheless a concrete cross-thread clue: the same
copy-plus-source grammar seen in the violation ladder appears inside a low
projective Hunter packet, now with scale eight as the copied channel and one
primitive speed as the source.

The tropical evaluator is itself recursively describable. If
`tau(G)=MST(c)` and a new comb vertex `v` has incident credits `d_i`, then

```text
tau(G+v)=max_(pi partition of V)
 [sum_(B in pi)tau(G[B])+sum_(B in pi)max_(i in B)d_i].      (6.9f)
```

Deleting `v` from a tree proves the formula. It shows exactly what insertion
must preserve: induced subtree values and best incident credits on every
partition block. The current tree value or threshold word alone is not closed
under the next insertion.

Equivalently, put graph edges at tournament vertices and order them by credit.
The resulting tournament is transitive away from its declared tie gauge. With
the graphic-matroid rank sidecar, its binary derivative

```text
k_j=r({e_1,...,e_j})-r({e_1,...,e_(j-1)})
```

has six ones and gives `MST(c)=sum_j c_(e_j)k_j`. Without endpoint incidence,
the same transitive tournament cannot determine the word. At each credit
level, contracting higher levels gives a three-way edge classification: loops
occur in no maximum tree, bridges in every maximum tree, and other nonloops in
some but not all. This is a literal tournament/metagraph correspondence: edge
order is the tournament shadow and graphic circuits are the information it
destroys.

The strict margin also gives an all-packet endpoint criterion. For a fixed
prefix `E` with `c_E` components and `e=mu(E)`, let `Gamma_h(x)` be the least
`sum_(ij in T)1/gcd(x_i,x_j)` among the trees maximizing the projective sum.
Evaluating the restricted credits on such a tree gives

```text
L_H >=19e/572-2c_E[sum_i 1/x_i+Gamma_h(x)].                (6.9g)
```

In particular, if `x_i=g a_i` for arbitrary distinct positive `a_i`, then
sorting gives `sum_i1/a_i<=H_7=363/140`, while the six tree gcd terms sum to
at most `6`. With no pairwise-coprime hypothesis,

```text
L_H>=19e/572-(1203/70)c_E/g.                              (6.9h)
```

Every fixed common-dilate ray is Hunter-positive beyond an explicit common
scale. This is new asymptotic territory, but not the full radius-seven chart:
incoherent small-gcd packets can keep the endpoint term in (6.9g) large.

The exact six-packet replay makes the two layers visible. All four positive
Hunter pilots have connected positive-credit graphs. Both consecutive pilots
have all 21 `eta_ij<0` and all 21 `c_ij<0`, although 14 and 18 of their 21
global projective defects are positive. Their conditional-overlap tournaments
remain transitive with the same coarse fingerprint as the successful pilots.
The signed-credit connectivity grade, not that tournament fingerprint,
detects the change.

This also distinguishes linear tomography from tropical selection. The
Radon transform first adds valued atoms along affine lines. Hunter first
forms restricted pair atoms and then selects a maximum tree. Averaging edge
credits before taking the tree is not equivalent to taking the tree in every
stalk. A useful combined object is consequently a **stalkwise tropical Radon
algebra**: linear line-sum reconstruction on the typed carry module, followed
only at the declared proof endpoint by the tropical tree character. This is a
concrete bridge between the Fano/chi-seven probe and Kakeya comb needles, not
an identification of their underlying geometries.

The live primitive-H6 pull supplies the dual tropical graph operation. For a
missing-label row `R` with retained core `P=[12] minus R`, a full antipodal
pair `{r,13-r}` creates an oriented AP cusp at `a/13`, `a=r^(-1) mod 13`.
The exact core-safe germ gives a vertex cap

```text
B_r(P)=2/ell_r(P),
ell_r(P)=min_(p in P) c([ap])/p,                         (6.10)
```

and a possible provider `s` gives the directed edge weight

```text
w_(r->s)=c([as])/2.                                     (6.11)
```

Tightness forces the disjunction

```text
u_r<=B_r(P)  or  u_s>=w_(r->s)u_r for some s!=r.        (6.12)
```

If every vertex avoids its cap, choose one provider edge from each vertex.
The resulting functional digraph contains a directed cycle, and multiplying
(6.12) around it forces

```text
product_(r->s in C) w_(r->s)<=1.                        (6.13)
```

After taking logarithms this is a nonpositive directed-cycle condition. It is
not a tournament condition: a pair can carry two arrows or none because its
two oriented germs use different gauges. The twenty primitive-core rows with
three full antipodal pairs have minimum-cycle-product census
`1/16:6, 1:2, 3/2:12`. The twelve expanding rows and two equality rows force a
cap; exact fixed-coordinate recursions then close the remaining forced slices.
That local algebra closes fourteen rows. A later exact longest-component
recursion closes the six product-`1/16` rows in 2,653,600 states. Thus all
twenty three-pair rows are loose, leaving at that checkpoint 903 rows with at
most two full pairs and the exceptional mixed-parity branch. THM-857's later
all-root component recursion closes that entire residual at scale one. This
shows exactly what the cusp quotient forgot: the literal residual component
word, remaining labelled progressions, last speed, and shortcut witness.

This produces a useful tree/cycle duality:

```text
Hunter:   undirected edge credits -> max-plus spanning tree,
AP pins:  directed handoff ratios -> min-plus cycle obstruction.  (6.14)
```

Both evaluations require vertex-edge incidence and both become simple only
after the correct local stalk has been formed. The antipodal count `f(R)` is
only a routing grade: it preserves the number of AP cusps but forgets their
signs, owners, exact germ lengths, and height parity. This is another concrete
case where a symmetry class predicts which algebra to use without being a
predicate-preserving quotient itself.

The next experiment should replace the scalar `f` in (6.1) by a value in the
free module generated by

```text
(node, exact defect, factor colour, carry role, seam, owner, metric obligation). (6.5)
```

One then computes ranks only after applying declared forgetful functors. The
scalar plane-rank theorem is a ceiling after residue aggregation, not a
certificate for the integer address tensor. A structured image may need fewer
directions or carry contrasts, but that must follow from a proved restriction
on the valued function, not from the Fano picture alone. Radon acts on Parikh
vectors; ordered j=4 peel chronology remains a separate history stalk until a
commutation theorem proves otherwise.

## 7. The ladder has a source term

The violation ladder's dyadic breakpoint sets obey the exact relation

```text
D_(2p) = D_p union (1/2)D_p
         union {1/m : p < m <= 2p, m odd}.                 (7.1)
```

The last set is a birth term. It is the information that a pure rescaling
picture drops. The corresponding cumulative counts match the classical
toothpick totals at `1,2,4,8`, giving `1,3,11,43`, and fail at `16`, where the
ladder has `159` rather than `171`.

This is not a disappointing analogy. It identifies the exact recursive datum
that must be preserved: inherited breakpoints plus newly primitive odd
denominators. The same pattern occurs in (5.7), where divisor transport alone
does not create the new primitive-period packets. Recursion is generally
**copy plus source**, not exact self-similarity.

There is nevertheless an exact self-affine subchannel. If
`W_a^H(s)=sum_(j=0)^H zeta^(2^j) 1_(s<a/2^j)`, then

```text
W_a^H(s)=zeta 1_(s<a)+zeta^2 1_(s<a/2)+zeta^4 1_(s<a/4)
         +W_a^(H-3)(8s).                                  (7.2)
```

The full ladder is therefore a self-similar chi-labelled stalk plus the
non-self-similar primitive-odd source in (7.1). This is the same odometer
pattern as the residue-plane base plus carry-simplex stalk in (6.3).

## 8. What the four-dimensional LRC object is

There are two exact uses of "four-dimensional" in the current frontier, and
they must remain distinct.

First, after the `f<=3` branch, fix one of the `binom(14,9)=2002` nine-speed
cores. The four ordered far speeds form a rank-four semilinear cone chart. It
has radial and transverse noncompact directions. A point must carry an
owner-coloured threshold loop, rational gaps, marked relations, residue
addresses, and the next proof action. The bare `Z^4` coordinate does not decide
the metric predicate.

Second, an affine family `V_i(c)=c a_i+r_i` is the integer-slope fibre of

```text
Phi_(A,R)(u,t)=min_i ||a_i u+r_i t||,
X_(A,R)={(u,t,c,lambda): u=ct, Phi_(A,R)(u,t)>=lambda}.     (8.1)
```

This is a mixed two-torus, integer-slope, clearance suspension. For fixed `c`
it has only two continuous coordinates. It is not a second four-manifold and
does not turn the global LRC14 moduli space into a four-dimensional space.

The tournament metagraph can serve as a **base atlas** for combinatorial
types of these fibres. The metric and owner data form a stalk over it. The
target is not a map from all LRC packets to one node; it is a stratified
pullback in which a node change transports the stalk and every wall crossing
updates its owner and clearance.

At `n=14`, one literal tiling has 78 bits. Reflection rewrites these as

```text
fixed bits f=6, paired objects r=36, h=f+r=42,
defect bits delta in F_2^36,              h+r=78.          (8.2)
```

The chi-seven grading rewrites the same 78 bits as one neutral block of size
12 and six nonzero charge blocks of size 11. These are two coordinate systems
on the same bits, not 156 independent coordinates. Other necessary views are

```text
six nonconstant Walsh levels H_2,...,H_12,
91 directional H derivatives,
the 2^36 factor-defect deck at each object,
marked path holonomy,
owner/wall/carry/clearance data.                            (8.3)
```

The apparent size is a reason to seek transforms and coherent algebras, not a
license to average away labels.

THM-856--858 sharpen this base/stalk picture into four interacting layers:

```text
projective ratio skeleton,
arithmetic ramification base,
literal metric-component fibre,
legal operation history.                                      (8.4)
```

THM-856 gives a positive horizontal projective tree margin `19/572`.
THM-857 solves the complete H6 metric fibre over scale one, with `2[12]` as
the only covering equality. THM-858 makes the H5 ramification base finite. What is
still missing is not another scale-one node census: it is the vertical
transport law joining fibres. Analytically, that coupling is the centered
prefix--comb--comb third moment on small-gcd maximizing trees; combinatorially,
it is the endpoint/progression cocycle and correlated AP-window word. The four
layers are a carrier decomposition, not four Euclidean coordinates.

## 9. The universal recursive definition

Let `O` be a declared alphabet of legal operations: face restrictions,
internal deletion with seam repair, inverse lift, CF return, arc flip,
reflection, path-relative complement, factorization, or an LRC observer
update. Let `P` be the declared terminal predicate or valued obligation.
Define

```text
x ==_(O,P) y
  iff P(w x)=P(w y) for every legal operation word w.      (9.1)
```

This is the operation-indexed Nerode equivalence. It is the smallest quotient
that is simultaneously a congruence for the chosen operations and sufficient
for the chosen task. Changing `O` or `P` changes the quotient. Static
injectivity, one-step lumpability, and LRC metric closure are three different
claims.

Across sizes, the literal objects and face maps form a presheaf on the
operation category. A proposed statistic is safe exactly when it is a natural
quotient of that presheaf for the declared subcategory. A failure of
naturality is not noise; its kernel pair is the next defect stalk.

This reframes several historical results uniformly:

```text
THM-801: literal B3 face descent is functorial.
THM-813: projected node/edge CF descent can fail functoriality.
THM-840: a kernel may be a congruence for insertion but not deletion.
THM-848: H requires its Walsh/Krylov operation stalk.
THM-851: node-defect colours close under factor composition through n=7.
THM-854: a fixed-path quotient drops the square holonomy.
THM-841: dyadic scaling needs a primitive-denominator source term.
```

## 10. Formula ledger for future tournament transfers

The following small formulas are proved or exact finite targets worth keeping
together.

```text
tiling dimension:       m_n=binom(n-1,2)
half-object dimension:  h_n=floor((n-1)^2/4)
reflection defect rank: r_n=floor((n-2)^2/4)
full reconstruction:    m_n=h_n+r_n
size increments:         Delta(m,h,r)=(n-1,floor(n/2),floor((n-1)/2))
line count:             |L_n|=2^(m_n-1)                 (m_n>0)
blue-line fraction:      |L_n^B|/|L_n|=2^(-r_n)
cyclic triangles:       C3=binom(n,3)-sum_i binom(d_i,2)
mean Hamilton paths:    mu_n=n!/2^(n-1)
H drift:                b_H=2(K-nH)/binom(n,2)
K drift:                S K=a_(n-3)+a_(n-2)-binom(n,2)H
projective H13 kernel:  rank=6, kernel=reflection-even sector
projective tree floor:  MST(h)>=-237/7436; margin>=19/572
restricted covariance: G^E_ij=c_ij-(2/13)(s_i+s_j) is PSD
Hunter defect:          L_m=(165-22m)e/169-sum_i s_i+MST(c)
common-dilate Hunter:   L_7>=19e/572-(1203/70)c_E/g
Kruskal derivative:     MST(c)=sum_j c_(e_j) Delta r_j
score OU law:           E[Delta x|T]=-8(x-n(n-1))/(n(n-1))
defect composition:     delta_total=delta_1+delta_2
binary factor rows:      B -> (BB,KK)=(1,2^r-1)
                         K -> (BK,KB,KK)=(1,1,2^r-2)
black endpoint parity:  4 | Q_black, 2 | selfK
witness norm:            N_(L_pi)p_0=[ord(pi) odd]1
Radon rank:              rank R_D=1+6|D|
B3 enumerator:           T_n=h_(n-3)(X,Y,Z)
deletion ownership:      (n-a)+(a-b-1)+(b-1)=n-2.
```

These formulas do not imply one another. They identify compatible coordinates
that can be joined and then subjected to (9.1).

## 11. Exact next experiments

1. Compute THM-851's factor closure at `n=8` without materializing the roughly
   billion naive factor incidences. Use defect translation/convolution and
   test whether the ordered closure reconstructs all arrows and whether the
   Jordan closure reconstructs all reflection orbits.
2. Join the `n=7,8` factor-defect colours with `Grad_H`. Determine whether
   gradient cells are already factor-coherent, or record the first crossed
   pair of fibres. This directly compares composition and flip operation
   algebras.
3. Extend the black self-line census to odd `n=9`, retaining witness cycle
   type, automorphism weight, and path-holonomy module. Test (4.4), not the
   refuted SC count.
4. Build the full eight-direction `F_7^2` valued Radon deck for the 21 j=4
   floods, with the A/B/C carry simplex attached before residue pushforward.
   Compute scalar, node-valued, defect-valued, and future-refinement ranks
   separately; compare each with the exact `48+30` address split.
5. Search for a natural role map `tau` on the 21 floods that commutes with
   subtree restriction and reflection. Fano incidence alone is not such a
   map because the exact numerical stabilizer is trivial.
6. Attach owner-event cocycles and clearance margins to the four-far scale
   action. Test whether the resulting packet is a congruence under peel and
   insertion, rather than only under static residue classification.
7. Use THM-857's closed 924-root scale-one certificate language as boundary
   data for arbitrary-scale/common-sheet transport. Lift its component word,
   remaining progressions, last speed, and shortcut witnesses to each
   ramification fibre. Attach the seven-comb projective edge types, rational-
   period stalks, Kruskal rank word, and reciprocal-gcd tree cost. Remove every
   common-dilate ray with (6.9h), then test whether the remaining small-gcd
   endpoint packets share a finite operation-stable third-moment language
   before invoking an AP-window argument.
8. Join the Hunter tree word to the AP-pin cycle word. Test whether every
   negative tree certificate either has an expanding directed handoff cycle
   or descends to a residual-endpoint grammar of the type that closed the six
   product-`1/16` rows. Retain exact endpoints when that implication fails.
9. Finish THM-853's operation-semigroup census and minimize the joined
   `(node, line, defect, factor, gradient, holonomy)` state under the actual
   multiword language. Compare it with each marginal before adding LRC data.

## 12. Assumption challenge

Tournament Analysis need not put runners or tournament arcs at its vertices.
Depending on the question, legitimate vertices include:

```text
gap positions, face directions, factor colours, defect sectors,
Hamiltonian needles, path masks, affine-line directions,
owner events, Fourier levels, and proof obligations.
```

For each choice, the pairwise observable and gauge must be declared, followed
by the tie Hamiltonian path and the standard fingerprint. More importantly,
the analysis must state which LRC predicate its quotient preserves and what it
destroys. A transitive carrier tournament can still lose higher overlap,
absolute metric, or continuation truth.

The most promising recursive object is therefore not a larger untyped graph.
It is the smallest natural, operation-closed quotient of the literal typed
incidence system for the terminal proof obligation. The computations above
show several of its stalks exactly. They also show, with unusual consistency,
where every smaller description breaks.
