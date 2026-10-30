# Proofs with exact interfaces

**PROVED for the quantified models below.** These are elementary arguments
assembled from inherited signed-incidence/odd-cover mechanisms and classical
constraint-preserving operations. They have received same-agent review and
finite exact controls, not a Lean build or an independent-agent audit. P3–P4
are a specialization of an established LP/affine viewpoint, not a claim to a
new general CSP classification.

## P1. Odd descent and the canonical two-cover

Let a finite graph have a totally ordered domain D_v at every named vertex v.
An oriented edge e:u->v requires x_v=f_e(x_u), where f_e is a fixed strictly
increasing or strictly decreasing function on its domain of definition.
Partial functions are allowed. Every used value must be in their domains.
Allow arbitrary unary allowed subsets at the vertices.

A d-sheet matching cover has variables (v,i), 1<=i<=d. For each edge choose
a permutation pi_e of the sheets and impose

    x_(v,pi_e(i)) = f_e(x_(u,i)).

The same edge function and the same unary restrictions are used on every
sheet. Sheet permutations for different edges need not agree.

**Theorem.** If d is odd, taking the median of each fibre produces a solution
of the original network. Consequently base feasibility is equivalent to
feasibility on each odd-sheet matching cover.

**Proof.** Permuting a list does not change its median. A strictly increasing
function takes its middle order statistic to the middle order statistic of
its image; a strictly decreasing function reverses the order and still takes
the middle element to the middle element. Thus each edge equality descends.
The median is an actual fibre entry, so every unary restriction and partial
function domain remains satisfied. Conversely, repeat a base assignment on
every sheet. This works for every matching and every d. QED.

**Finite-cover compression.** Define the canonical two-cover by preserving
sheet labels on increasing edges and swapping them on decreasing edges.
From any finite-cover solution assign the fibre minimum and maximum to its
two copies. Increasing maps send minimum to minimum and maximum to maximum;
decreasing maps exchange them. These extrema are actual entries, so all
unary restrictions survive. Thus existence on some finite cover is equivalent
to existence on this particular two-cover.

If all edge functions increase, the fibre minimum already descends to the
base for every finite degree. More generally this holds when the edge signs
are balanced: reverse the domain orders at suitable vertices to make every
edge increasing. Balance means that the product of monotonicity signs is +1
around every cycle. There is no assertion that every unbalanced network has
an obstruction.

**Sharp scope.** Extrema do not preserve future unary-support queries:
[0,1,2] and [0,0,2] have identical extrema, but the query asking whether a
sheet can equal 1 separates them. Multiplicities and weighted sums are also
lost. The theorem fixes the constraint network before compression. A pin on
just one sheet is not a base unary pin.

An odd function need not qualify: f(x)=x^3-3x satisfies f(-x)=-f(x), but
median(f(0),f(1),f(2))=0 differs from f(1)=-2. Functional oddness and order
compatibility are different hypotheses.

**Fixed-clock arithmetic corollary.** For a specified valuation a, the
accelerated Collatz branch U_a(n)=(3n+1)/2^a is strictly increasing on the
ordered set of positive odd n with v_2(3n+1)=a. Therefore a network with
one fixed valuation at each edge has base feasibility if and only if it
has a solution on any finite matching cover: fibre minima suffice, with
no odd-degree restriction. The same argument applies to fixed branches
of 5n+1. It cannot prove convergence, since the actual branch choices
are variable and the 13->33->83->13 control is a valid fixed-clock cycle.
The corollary concerns solutions of a specified network, not the existence
of a uniform descending potential for all starting integers.

## P2. Signed-affine solvability, all-modulus counts and cover degrees

Let a finite graph carry equations

    x_v = s_e x_u + b_e,       s_e in {+1,-1}, b_e in Z.       (1)

Multiple edges are allowed. Unary conditions are absent unless explicitly
added below. In each connected component choose a root and spanning tree.
Tree propagation expresses every coordinate as

    x_v = a_v t + c_v,         a_v in {+1,-1}, c_v in Z.       (2)

Substituting in every remaining edge gives either 0=h_e or 2t=r_e,
after changing the sign of a row with coefficient -2. All quantities remain
integral; we have not divided a modular congruence by two.

**Theorem P2a.**

    (1) has an integer solution
      iff it has a real solution and a solution modulo 2.

**Proof.** In a component with no 2t row, real solvability says all h_e vanish,
and any integer t works. In a component with such rows, real solvability says
all h_e vanish and all r_e agree with one integer r. The unique root is r/2.
Modulo-two solvability is exactly r even. Equation (2) proves the conclusion
component by component. It also proves that every real-solvable system has a
half-integral solution. QED.

This is the integral mechanism behind an odd cycle of edge sums. When all
s_e=-1, a cycle is unbalanced exactly when its length is odd. With mixed
signs the relevant invariant is the sign product, not its length alone.

**Theorem P2b, exact all-modulus formula.** For a component let Delta be the
nonnegative gcd of all h_e and, when negative pins exist, all differences
r_e-r_0 for one chosen pin r_0. The gcd of an empty/all-zero list is zero.
Let g_q=gcd(2,q). For every integer q>=2 the number of solutions modulo q is

    balanced component:
      q,     if q divides Delta;
      0,     otherwise;

    unbalanced component:
      g_q,   if q divides Delta and g_q divides r_0;
      0,     otherwise.                                      (3)

Multiply these numbers over connected components.

**Proof.** Equation (2) is a bijection between root residues and assignments
satisfying the tree. All zero rows and pin differences hold iff q divides
Delta. The remaining equation 2t=r_0 has exactly g_q solutions iff g_q divides
r_0. The root is free if there is no negative pin. QED.

Consequently, writing Delta_* for the gcd of the component contents, the
system is solvable modulo an odd prime p exactly when p divides Delta_*.
If Delta_*=0 every odd prime works; otherwise precisely the odd prime
divisors of this explicit integer work. When consistent, the count is
p^b for b balanced components. This includes inconsistent real systems that
nevertheless have solutions at selected primes. For q=2^h, a nonzero
content also bounds the admissible precision by h<=v_2(Delta_*), with the
additional pin parity condition in (3). All zero contents mean real
compatibility, but do not decide the evenness of a negative pin.

The classical signed-incidence Smith form has only 1 and possibly 2 as
nonzero invariant factors: tree elimination leaves unit rows and, in an
unbalanced component, one row with entry 2. Formula (3) additionally retains
the translated right-hand sides and their compatibility content.

**Theorem P2c, finite-cover classification.** The minimum degree of a finite
matching cover on which (1) has an integral solution is:

- no such degree, if the base system is real-infeasible;
- one, if the base system is integer-feasible;
- two, otherwise.

In the third case no odd degree works.

**Proof.** Average a solution over each fibre. Affinity and the bijective edge
matchings yield a real base solution, excluding the first case. In a
real-compatible system choose a half-integral base solution x. On the
canonical signed double cover put floor(x_v),ceil(x_v). Increasing integral
translations preserve these two values, and

    floor(-x+b)=b-ceil(x),     ceil(-x+b)=b-floor(x).

Thus the double cover has an integer solution. P1 forbids any odd-degree
repair of an integrally infeasible base. QED.

Covers need not be connected. The theorem says that a particular double
cover suffices; it does not say every even cover does.

**Interval extension.** P2c continues to hold with closed unary intervals
whose finite endpoints are integers, allowing infinite endpoints. A
real-feasible free tree component has a nonempty root interval with integral
endpoints, so it has an integer root. A pinned component retains its
half-integral root. Floor and ceiling preserve these intervals. For P2a
one may use real feasibility with the intervals and modulo-two feasibility
of the equations; narrow intervals need not realize every residue, but
the preceding root argument proves existence of an integer solution.
Arbitrary unary holes are allowed in P1 but are not covered by this
real-plus-mod-two criterion.

**Minimal obstruction.** The triangle

    x+y=y+z=z+x=1

has unique real solution (1/2,1/2,1/2), no integral solution, and solutions
modulo every odd integer. Its signed double cover is a six-cycle with
alternating labels 0 and 1. Odd moduli cannot certify integral feasibility.
This does not contradict recovering a known bounded integer from enough
odd-modulus observations, as in the prior counting packet: that is a
different target with an a priori integral unknown.

## P3. Consistent odd-denominator marginals round exactly

Take a finite constraint problem with finite totally ordered variable
domains. Every constraint is a relation on a specified ordered tuple of
variables. Assume that each relation is closed under coordinatewise
median of every odd number of its tuples. Unary constraints/pins are
included in the language.

For each constraint c give a probability distribution lambda_c on its
actual allowed tuples. The unary marginals at every occurrence of a shared
variable must agree. Suppose all these probabilities are rational with
a common odd denominator D.

**Theorem.** The coordinatewise medians specified by these unary marginals
form a genuine solution.

**Proof.** List each tuple of relation c exactly D*lambda_c(tuple) times.
There are D tuples. The list of values at a variable has multiplicities
depending only on its unary marginal, hence the same median in every
constraint containing that variable. Closure puts the resulting tuple in
each relation. Variables in no constraint may be assigned any allowed
value; empty domains are immediately infeasible. QED.

There is no alignment assumption between different tuple lists, and no
claim that unary marginals determine joint statistics. Closure is what
makes the common median sufficient for feasibility.

Conversely an actual solution supplies such distributions with D=1.
Thus for this class, feasibility is equivalent to the existence of
consistent local marginals with odd common denominator.

An ordinary rational local relaxation is insufficient: Boolean disequality
is odd-median closed, but a triangle has consistent edge measures assigning
1/2 to (0,1) and (1,0) and has no two-coloring. Its local equations force
one-half marginals. Odd denominators eliminate precisely this witness.

### P3b. The operation need not be a median

The same proof works if, for every odd D, each variable domain has a
symmetric D-input operation and the coordinatewise operations preserve
every constraint relation. Symmetry means dependence only on the
input multiset. Shared marginals therefore give the same output at
every occurrence. This is the classical symmetric-polymorphism mechanism.

For affine binary constraints Mx=b, take XOR of an odd number D of
assignments. Linearity gives M(sum x_i)=D*b=b in F_2. Constant unary
pins are preserved too. For affine constraints over F_p and D coprime
to p, the symmetric operation D^(-1)*sum x_i works instead. Combining
with P4 over Z_(p) gives the corresponding local denominator criterion.
Gaussian elimination is already simpler for an explicitly supplied
linear system; this is a bridge between representations, not a new
linear-solving algorithm.

Median and XOR are not interchangeable. The even-parity tuples
110,101,011 have median111, which has odd parity, whereas their XOR
is000. This is why a binary code can have odd-copy descent without
having ordered-median descent. The same failure occurs in the repo's
seven-coordinate Hamming carrier, tested on its actual syndrome map.
No convex median-cost estimate is transferred to XOR.

## P4. Positive-face test for odd denominators

Let A and b be integral, P={x in R^n: Ax=b, x>=0}, and put

    R_(2) = Z_(2) = {a/d : a in Z, d a positive odd integer}.

This ring retains powers of two as nonunits. It is not Z[1/2].

Suppose P is nonempty. Define J to be the indices positive in at least
one point of P. Let A_J contain just those columns.

**Theorem.**

    P intersects Z_(2)^n
      iff A_J z=b has a solution over Z_(2).                 (4)

To test the right side, take an integer Smith normal form

    U A_J V = diag(d_1,...,d_r,0),       d_i>0.

Writing c=Ub, the necessary and sufficient conditions are

    c_i=0 for every zero row;
    v_2(c_i)>=v_2(d_i) for i=1,...,r,                       (5)

with v_2(0)=infinity.

**Proof.** Necessity in (4) is immediate because every feasible point
vanishes off J. Smith reduction consists of integral invertible operations,
so it works over Z_(2). Divisibility by d_i in that ring is exactly the
two-adic test (5).

For sufficiency, choose a solution z_0 over Z_(2). There exists a feasible
point x_* positive in all coordinates of J: average, over the finitely
many j in J, a feasible point with j positive. The real kernel of A_J
has a basis of integral columns H, obtained by clearing denominators in a
rational nullspace basis. Therefore x_*-z_0=Ht for some real vector t.
The ring Z_(2) is dense in the real line: approximate by rationals with a
large odd denominator. Approximate every coordinate of t closely enough
that z_0+Ht' stays positive. This is an exact Z_(2)-solution on J;
extend it by zero outside J. If J is empty the assertion is simply
b=0 with x=0. QED.

All coordinates of this solution have a common odd denominator, by taking
the least common multiple of their denominators.

**The face is essential.** Let

    A = [[1,2,0],[1,0,1]],     b=(1,0).

There is an unrestricted integer solution (1,0,-1). Nonnegativity forces
x_1=x_3=0 and x_2=1/2, so P has no odd-denominator point. Testing the full
unrestricted matrix instead of A_J would give a false answer.

**Finite exact algorithm.** Decide real feasibility. Find J by testing
whether each coordinate can be strictly positive (equivalently maximize
a capped auxiliary variable bounded by that coordinate and by 1).
Compute Smith data on A_J and test (5). A constructive rational-nullspace
approximation as in the proof produces a certificate when the test passes.
This packet establishes correctness, not a new bit-complexity bound.

Combining P3 and P4 yields a complete finite feasibility procedure for
the specified odd-median-closed constraint systems. Its equality matrix
encodes normalization and shared unary marginals, and its nonnegative
coordinates encode masses of actual allowed tuples. It need not enumerate
graph covers, and it rejects the triangle of disequalities.

**Generalization.** The same proof works with Z_(p), denominators coprime
to any fixed prime p, replacing v_2 by v_p. The special significance of
p=2 here is that odd arities admit a middle order statistic commuting
with both increasing and decreasing maps. This does not furnish a
p-ary analogue for arbitrary constraints when p is odd.

### P4b. Exact denominator exponent and cover divisibility

For the same nonempty P, let e(P) be the minimum v_2(D) among positive
integers D for which some x in P satisfies D*x in Z^n. Then

    e(P) = max(0, max_(1<=i<=r) [v_2(d_i)-v_2(c_i)]),      (8)

where a term with c_i=0 contributes zero. No zero-row incompatibility is
possible because P is real-feasible.

**Proof.** A rational x has a denominator with v_2(D)<=e exactly when
2^e*x belongs to Z_(2)^n. Apply P4 to the scaled polyhedron
{y>=0: A*y=2^e*b}; its support face is still J. Smith divisibility becomes
e+v_2(c_i)>=v_2(d_i). The smallest such e is (8). A smaller exponent is
impossible, and the positive-face construction attains this exponent. QED.

Each obstruction has a modular dual witness: choose a row y of U and h
with y*A_J divisible entrywise by 2^h but y*b not divisible by 2^h.
Every feasible denominator D must satisfy
v_2(D)>=h-v_2(y*b). Retaining the positive face is essential here too.

For an arbitrary finite constraint problem, form the usual local marginal
polytope: masses on actual allowed tuples, a marginal at each named variable,
normalization and occurrence-by-occurrence consistency. There is no median
hypothesis in the following statement.

**Cover-divisibility theorem.** If the local polytope is nonempty, (8)
is the minimum v_2(D) of a satisfying D-sheet factor-graph cover. If the
local polytope is empty, no such finite cover exists.

Here a factor-graph cover duplicates each variable and each constraint D
times, with a permutation matching at every named occurrence/port. Constraints
and their port order are unchanged. The cover may be disconnected.

**Proof.** A cover assignment produces normalized tuple and variable
frequencies with denominator D. Conversely, clear denominators of a local
solution. At each constraint create D copies carrying its allowed tuples
with the specified multiplicities; at each variable create D copies carrying
its marginal values. For each occurrence, equal marginal counts allow a
value-preserving bijection between its constraint copies and variable
copies. These bijections are exactly the cover matchings. This proves
the equivalence for each D, so (8) gives the minimum two-adic exponent. QED.

This does not determine the smallest total cover degree, its connectedness,
or its symmetry group. Odd factors may still be needed. With the median
hypotheses of P3, base feasibility is equivalent to e(P)=0. Thus every
infeasible problem in that class has either an empty real local polytope
or an explicit two-primary modular witness on its positive face.

**Higher-order control.** Give three variables the domain Z/q, with two
identity relations around a triangle and the third relation x->x+1 mod q.
There is no base solution. Marginal consistency forces every residue to
have mass 1/q, so a satisfying cover has degree divisible by q. A q-sheet
cover is obtained by giving each residue one sheet and rotating the last
matching. Its minimum degree is q and its denominator exponent is v_2(q).
For q=4 this exceeds the signed-affine exponent one. For q=3 even an odd
cover repairs the obstruction: the cyclic permutation is not monotone on
a linearly ordered domain and the relations are not odd-median closed.
This identifies the exact failed hypothesis in an attempted universal
odd-cover theorem.

## P5. Sharp cost and marked-point boundaries

Let d=2k+1, m=median(x_1,...,x_d), and phi be a nonnegative convex
real-valued function on an interval containing the x_i. Then

    phi(m) <= sum_i phi(x_i)/(k+1)
           = [2d/(d+1)] average_i phi(x_i).                 (6)

**Proof.** If a point left of m and a point right of m both had cost
strictly smaller than phi(m), convexity would contradict the value at m.
Consequently on at least one of those closed half-lines all costs are
at least phi(m). At least k+1 of the sample entries lie on that side.
Their sum proves (6), using nonnegativity for the remaining entries. QED.

The factor is sharp: take k zeros and k+1 ones with phi(x)=max(0,x).
Summing (6) with fixed nonnegative vertex weights also controls additive
unary convex costs. It does not cover arbitrary joint-port potentials.
The proof even applies to nonnegative quasiconvex costs, but no such
extension is needed for the computational claim.

For the point indicator phi(x)=1[x=0], take k copies of -1, one 0, and
k copies of 1. Median cost is 1 and average cost is 1/d. No uniform factor
exists for arbitrary point-supported weights. This is why a typical-source
or average-cost argument does not automatically control a designated
Collatz source or a particular runner boundary.

If at least k+1 sheets carry a designated value n, their median is n.
This strict-majority threshold is sharp. Merely giving n positive mass
is insufficient for median descent to retain it; requiring majority is
a new constraint that must itself be proved feasible.

## P6. The integral interface behind odd functions and transfer

Let X have an involution tau with f fixed points and r exchanged pairs.
Over Z/q the number of functions satisfying h(tau x)=-h(x) is

    gcd(q,2)^f q^r.                                        (7)

Each pair has one free value; each fixed point requires 2h(x)=0.
This is inherited from the October odd-function packet and checked here.

For M=Z^X let M^+ and M^- be its invariant and anti-invariant subgroups.
At a pair with values (a,b), the natural coordinates are S=a+b and D=a-b.
An integral pair is recovered exactly when

    S = D (mod 2),       a=(S+D)/2, b=(S-D)/2.

Thus M/(M^++M^-) is (Z/2)^r. Over a ring where 2 is invertible the parity
gluing disappears. Over Z_(2) it remains. Fixed points have no nonzero
anti-invariant integer value.

There is also a classical transfer statement. If homomorphisms res and Tr
on an abelian two-primary torsion group satisfy Tr(res(a))=d*a, with d odd,
then res is injective: res(a)=0 gives d*a=0, and multiplication by an odd
number is injective on every element of two-power order. With exponent
bounded by 2^h, choose u with ud=1 mod 2^h; then u*Tr is a left inverse.

These are exact mechanisms, not a claim that every appearance of an odd
number has a common action. An application to square classes, cohomology,
or Selmer groups must construct its actual restriction and transfer.
Injectivity under odd degree says nothing by itself about Sha finiteness.
