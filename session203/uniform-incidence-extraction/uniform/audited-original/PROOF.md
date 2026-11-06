# No direct affine extraction by a uniform binary-incidence gadget

Status: **PROOF CANDIDATE, pending independent audit.** This is an all-order
obstruction for a specified physical constructor, not an obstruction to the
unrestricted tournament spectrum. No eventual odd-spectrum threshold is proved.

## 1. Physical constructor and statement

Let A be any tournament, with no bound imposed on its order. For every vertex
a of A choose two bits e_a(0), e_a(1). A binary word b=(b_1,...,b_m) supplies
new vertices x_1,...,x_m. Orient x_i -> x_j for i<j, and orient a -> x_i
exactly when e_a(b_i)=1. All edges within A are fixed. Call the resulting
actual tournament T_A(b). The same pair of incidence bits is used at every
position; this is the **uniform binary-incidence** restriction.

A vertex a is sensitive if e_a(0) differs from e_a(1). The two sensitive
polarities, e_a(b)=b and e_a(b)=1-b, are both permitted. A fixed leading digit,
fixed prefix or suffix, or any fixed additional gadget can be absorbed into A,
provided its incidences to each new digit vertex have this stated form.

Put

    u_m(b)=1+2^m+sum_(i=1)^m b_i 2^(m-i).

Thus u_m(b)-1 has binary word 1b, and u_m runs through the full interval
[2^m+1,2^(m+1)]. The inherited E(u) family is included: put its root, its
sensitive apex, and its leading chain vertex in A. The forward order sector
of E(u) equals u, but its total count has a correlated extra term.

**Theorem 1 (quantitative uniform-incidence obstruction).** Fix real numbers
alpha>0 and beta. If m>=2 and

    N_m/16 > alpha*(3*2^m+1)/2 + beta,                 (1)
    N_m = 3^m-(3+m/2)2^m+2m+3,

then no choice of A or its two incidence types satisfies

    H(T_A(b)) = alpha*u_m(b)+beta for every b in {0,1}^m.  (2)

The size, internal arcs and incidence types of A can be chosen separately at
each length m. There is no fixed-size hypothesis. For any fixed alpha,beta,
(1) holds for all sufficiently large m, since (3/2)^m dominates every linear
function of m.

**Corollary 2.** No such constructor realizes either 6u-1 or 6u+1 uniformly
on the single full block 8193<=u<=16384. In particular, no uniform
binary-incidence gadget extraction of either expression works for all u.
This excludes an arbitrary-size version of the previously tested fixed
small grafts; it does not exclude arbitrary surgery on E(u), changing the
incidence rule with the position, target-dependent gadget selection, or a
constructor with a different encoding of u.

For m=13 one has N_m=1516528 and N_m/16=94783. The larger target average is
6*(3*8192+1)/2+1=73732. These are exact integers.

## 2. Inheritance, novelty and live concepts

The closest mechanism is the positive weighted-cut formula of Session192,
`05-knowledge/results/session192_surprising_approaches_20261105/weighted_cuts/PROOF.md`,
Section1. It already gives strict quadratic interactions when one apex is
varied. That formula is not claimed as new here. Session199 source-cut work
specializes related strict interactions to physical O/Z words. Session202's
`positive-overlap/EDGE_SECTOR.md` supplies the free edge-order sector and the
open 6u+-1 extraction target; its bounded graft probes were not all-order
obstructions.

The new combination is a dichotomy: two sensitive apexes force a three-slot
exponential average, while exactly one sensitive apex has unavoidable strict
quadratic interactions. Together they cover every fixed uniform binary
incidence table, irrespective of gadget size, and give the explicit dyadic
block obstruction above. No claim of literature priority is made.

Anchor: the eventual odd spectrum remains open. Niche: test whether free
sector control can be converted directly to a scalar count. Wildcard:
complementing digit incidences does not avoid either half of the dichotomy.
Live concepts are physical bits, strict cut curvature, three-slot entropy,
induced path injection, and the distinction between pointwise template choice
and one uniform template on a whole digit cube.

## 3. Induced monotonicity, with its injection

For any tournament T and induced subtournament U,

    H(T) >= H(U).                                    (3)

To see it, fix an ordering of vertices outside U. Insert them one at a time
in any Hamilton path of U. For a new vertex v and old path p_1,...,p_t,
prepend v if v->p_1. Otherwise, if p_t->v append v. In the remaining case
the indicators p_i->v begin with1 and end with0, so there is a transition
p_i->v->p_(i+1); insert v there. Resolve any choice deterministically.
Deleting all inserted vertices recovers the original order on U, proving
injectivity. This is an actual path injection, not a count heuristic.

## 4. Two sensitive apexes force a three-slot average

Choose any two sensitive vertices a,c of A. Work only with their induced
tournament and the m digit vertices. For each assignment of digit vertices
to three ordered bins L,M,R, put the vertices within each bin in increasing
index order and form

    L, a, M, c, R.                                  (4)

Require L and R to be nonempty and |M|>=2. All internal chain edges are
correct. Exactly four other conditions are needed: last(L)->a,
a->first(M), last(M)->c and c->first(R). They involve four distinct digit
vertices. Under independent uniform random digits, each prescribed edge has
probability1/2, regardless of either sensitive polarity. The four events are
independent. Therefore (4) is a Hamilton path with probability1/16.
Different assignments give different vertex words.

There are exactly

    N_m = [3^m-3*2^m+3] - m[2^(m-1)-2]
        = 3^m-(3+m/2)2^m+2m+3                     (5)

such assignments for m>=2. The bracket in the first line counts assignments
with all bins nonempty; subtract those with the middle bin a singleton.
For m=2,3 both sides are0. Linearity of expectation and (3) give

    average_b H(T_A(b)) >= N_m/16.                  (6)

No independence between different assignments is claimed or needed.
Formula (6) is valid however large the rest of A is and whatever its arcs.

## 5. One sensitive apex has strictly negative pair coefficients

Here is a self-contained version of the inherited cut fact. Let U be any
nonempty tournament, and adjoin v with v->i iff z_i=1. Write Q_ij for the
number of vertex words of U in which i is immediately followed by j and
all steps except that exposed step are directed correctly. Set

    c_ij=Q_ij+Q_ji.

Then the exact Hamilton-path count is the multilinear polynomial

    H(U+v)=H(U)+sum_(i->j in U) c_ij z_i(1-z_j).     (7)

To verify (7), remove v from a Hamilton word. If v is first or last, its
contribution is sum_i z_i start_i + (1-z_i) end_i. If v is internal, its
contribution is sum_(i!=j) Q_ij(1-z_i)z_j. The divergence identity
start_i-end_i=sum_(i->j)Q_ij-sum_(j->i)Q_ji follows by summing directed
steps along every old Hamilton path. Cancel crossing terms to obtain (7).

For distinct i,j the capacity is strictly positive. More precisely,

    c_ij >= 2 H(U-{i,j}) >= 2,                       (8)

with H(empty)=1. Take any Hamilton word z_1,...,z_k of U-{i,j}, and insert
the block (i,j) or (j,i), leaving its internal step exposed. If k>=1 put
A_l=1_(z_l->i), B_l=1_(z_l->j). The number of valid insertions of both
block orientations over all gaps is

    2+sum_(l=1)^(k-1) [(1-B_l)A_(l+1)+(1-A_l)B_(l+1)] >=2.

At k=0 both exposed words work. Deleting the block recovers the old word,
so there is no overcount across old words. This proves (8).

Now assume A has exactly one sensitive vertex v. Removing v leaves a
fixed tournament U independent of b. At every digit vertex, z_i=b_i for
one polarity, or z_i=1-b_i for the other. All remaining z-coordinates are
fixed. By (7), for every pair i<j of digit positions the coefficient of
b_i b_j is exactly -c_(x_i,x_j), for either common polarity, and is strictly
negative by (8). Thus H(T_A(b)) has a nonzero quadratic coefficient whenever
m>=2. It cannot equal an affine-linear function of b on the entire Boolean
cube: multilinear polynomial representation on {0,1}^m is unique.

This last uniqueness can be checked without an algebraic theorem. The
mixed finite difference at positions i,j is the displayed nonzero constant;
an affine-linear function has zero mixed finite difference.

## 6. Complete dichotomy and scope

If no vertex of A is sensitive, H(T_A(b)) is independent of b, while
alpha*u_m(b)+beta is not. If exactly one is sensitive, Section5 contradicts
(2) for m>=2. If at least two are sensitive, averaging (2) would give
alpha*(3*2^m+1)/2+beta, contradicting (1) and Section4. This proves Theorem1
and the stated numeric corollary.

The map from digit words to actual tournaments preserves every vertex and
orientation. It forgets no correlated coordinate by declaring it free.
The obstruction is to **the same incidence table on every word of a full
cube**. It says nothing about choosing among many tables after inspecting b,
changing old digit-to-digit arcs, feeding derived digits, or changing digit
position rules. It is not a density-zero theorem for all outputs of such
gadgets, and is not a nonexistence statement for any individual target in
the dyadic block. All those target counts are already covered by finite
physical certificates in the inherited spectrum work.

## 7. Literal validation

`verify.py` uses subset Hamilton-path dynamic programming, explicit exposed
path permutations, and independent bin assignments. Its finite controls
check (7),(8) for every labeled parent through order4 and all apex cuts;
check all sensitive polarity pairs and both apex arcs for m=2..8;
check the three-bin inventory and its four independent edge conditions;
and check the exact m=13 bound. These validate implementation and constants.
The all-order obstruction follows from Sections3-6, not from the census.

Run:

    python3 verify.py > VALIDATION_NORMAL.json
    python3 -O verify.py > VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json
