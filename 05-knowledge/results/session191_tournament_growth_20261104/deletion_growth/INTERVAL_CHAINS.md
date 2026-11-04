# Exact deletion growth inside Moon's interval chains

**PROVED derivation from inherited chain and one-arc-join laws; finite checks
are retained separately. Independently audited.**

This explains why the unrestricted deletion minima begin like Fibonacci
numbers but cease doing so. It proves an all-order formula **inside the
interval-chain family**, not for arbitrary strong tournaments.

## 1. Definitions and credited inputs

For a composition p=(p1,...,pk), all pi>=2, let IC(p) have ordered vertices
x0,...,xN, N=sum pi. Every arc points forward except the block endpoint arcs
x_(li+pi)->x_li, where li=sum_(j<i) pj. Put f(p)=2^(p-1)+1.
The empty composition denotes a singleton tournament.

Moon's known product formula, written out in repository THM-4744,
`01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, is

    H(IC(p1,...,pk)) = product_i f(pi).                       (1)

All nonempty chains here are strong. The second input is the previously
proved one-arc-join law (Session29
`05-knowledge/results/session29_framed_profiles_20261022/session29-profiles/ADDENDUM.md`, sections 2 and 5;
restated with attribution in Session190's imported
`single_cycle_floor_crevasses_lorentzian_joins_20261104.md`, section 4).
For disjoint rooted tournaments (A,a),(B,b), form J from A->B by reversing
only a->b. If e_a is the number of paths ending at a, s_b the number starting
at b, R_a the number of ordered spanning two-bin path covers whose second
bin is nonempty and starts at a, and R'_b the converse end-bin response, then

    H(J) = H(A)H(B) + R_a(A)R'_b(B) - e_a(A)s_b(B).           (2)

The bins may otherwise be empty. Formula (2) follows by distinguishing paths
that cross back on b->a from paths that do not. It is not new here.

## 2. The endpoint response of a chain

Let X=IC(sigma,p) and r its last spine vertex. Let A=IC(sigma), with last
vertex u. The p-1 new interior vertices form a transitive set L, and r has
only one out-neighbor, u. Write a=H(A). Then

    H(X)=a f(p),    e_r(X)=a,    R_r(X)=H(X),    s_r(X)=1.     (3)

The H identity is (1). Removing r leaves A->L, so every Hamiltonian path
of X-r ends in L and can be followed by r; this gives e_r=a.

For R_r, the second path is either the singleton r, giving a choices for
the first path, or begins r,u. In the latter case its A-part starts at u
and its complementary A-part is the first path. Every vertex of L can be
assigned independently to either path, in forced increasing order. Thus

    R_r(X) = a + 2^(p-1) R_u(A).                             (4)

Starting with R=H=1 for the singleton, induction gives R_r=H(X).
Similarly, a Hamiltonian path starting at r must continue with u, then the
remaining A-vertices and then L. Its number is s_u(A)=1 by induction.
Taking the converse yields the first-endpoint versions of (3).

## 3. Every strong-deletable vertex and its exact increment

Deleting an internal spine boundary between two blocks makes all remaining
vertices on its left dominate all those on its right, so the deletion is
not strong. Deleting the first or last spine endpoint creates a source or
sink. Therefore the only possible strong-deletable vertices are interiors
of blocks. All those interiors are strong-deletable when the chain has at
least four vertices, as the following constructions also show.

**A block of length p>=3.** Removing any of its p-1 interior vertices
replaces that part by p-1. Consequently

    Delta = 2^(p-2) product_{other blocks q} f(q).            (5)

**A triangle block, p=2.** Write the chain as IC(sigma,2,tau), and let v be
the middle of this triangle. Put X=IC(sigma), Y=IC(tau), with distinct
endpoint roots a=last(X), b=first(Y). On deleting v, the remaining tournament
is precisely J(X,a;Y,b). Using (2),(3),

    H(T)=3 H(X)H(Y),
    H(T-v)=2 H(X)H(Y)-e_a(X)s_b(Y),
    Delta=H(X)H(Y)+e_a(X)s_b(Y).                             (6)

This deletion is strong whenever X or Y is nontrivial, hence whenever
there is any other block. For example, if X is a singleton and Y is strong,
the reversed arc b->a and the other forward arcs connect both ways; the
other cases follow identically. The one-block C3 has no strong deletion and
is excluded from this family of marked pairs.

If the triangle is internal and its neighboring block lengths are p,q,
let c be the product of f over all other blocks. Formula (3) gives

    Delta = c (f(p) f(q)+1).                                (7)

If the triangle is at one end and its neighboring block length is q, the
corresponding formula is

    Delta = c (f(q)+1) = c (2^(q-1)+2).                      (8)

These formulas isolate the neighboring blocks; the rest of the chain enters
only through a multiplicative factor. They retain the marked vertex and are
not an assertion that a scalar H determines arbitrary insertion response.

## 4. Minimize exactly over every chain order

Let M(t) be the minimum of product f(pi) over compositions of t into parts
at least two. Set M(0)=1 and M(1)=M(t)=infinity for t<0. For t>=2,

    M(3j)=5^j,
    M(3j+1)=9*5^(j-1)       (j>=1),
    M(3j+2)=3*5^j.                                         (9)

Indeed a part p>=5 can be split as 2+(p-2), strictly lowering its product:
f(p)-3f(p-2)=2^(p-3)-2>0. A part4 has the same cost as two part2s, while
three part2s cost27 and can be replaced by two part3s costing25. Thus at
most two 2s are needed, giving (9).

Let g(n) be the minimum Delta over all strong deletion pairs with T an
n-vertex interval chain. Put N=n-1. Then, for every n>=4,

    g(n) = min {
       2 M(N-3), 4 M(N-4), 6 M(N-5), 10 M(N-6),
       16 M(N-7), 26 M(N-8), 46 M(N-9), 82 M(N-10)
    }.                                                      (10)

Here is a complete reduction to the finite list, rather than an observed
recurrence. Unmarked parts at least5 can be split as above. In (5), a marked
part p>=5 can also be split into an unmarked2 and a marked p-2; its
coefficient falls from 2^(p-2) to 3*2^(p-4). So marked parts3 and4 suffice.
In (8), replacing a neighbor q>=5 by an outer2 and neighbor q-2 changes the
coefficient by

    (f(q)+1)-3(f(q-2)+1)=2^(q-3)-4 >=0.

In (7), the same replacement at a neighbor p>=5 changes it by

    f(p)f(q)+1 - 3(f(p-2)f(q)+1)
       = (2^(p-3)-2)f(q)-2 >0.

Thus neighbors2,3,4 suffice. Their coefficient/cost pairs are
(10,6),(16,7),(28,8),(26,8),(46,9),(82,10); the pair(28,8) is dominated by
(26,8). End triangles add (4,4),(6,5),(10,6), and marked parts3,4 add
(2,3),(4,4). All outside blocks have minimum product M at the remaining
cost, and that minimum can be placed on either permitted outside side.
Every term in (10) is therefore attainable whenever it is finite. This
proves (10) in both directions.

Substituting (9) in (10) gives

    n       4  5  6   7   8   9  10  11
    g(n)    2  4  6  10  16  26  46  78.                  (11)

For every k>=0 the exact all-order formulas are

    g(9+3k)=26*5^k,
    g(10+3k)=46*5^k,
    g(11+3k)=78*5^k.                                       (12)

For completeness, (9) satisfies M(t+3)=5M(t) for t>=2. Therefore every
term in (10) scales by5 when N increases by3 once N>=12. The finite
boundary values N=8,...,14 are 26,46,78,130,230,390,650, which connect the
three residue classes to that uniform recurrence and prove (12).

Simple attaining families are

    IC(3,2,3, 3,...,3),       with k trailing 3s;
    IC(3,2,4, 3,...,3),       with k trailing 3s;
    IC(2,3,2,3, 3,...,3),     with k trailing 3s,

respectively, always deleting the middle vertex of the displayed internal
triangle block. The third line's marked triangle is its second part2.

For example IC(3,2,4) has H=135; deleting its marked middle vertex gives89,
so the increment is46. The Busch floor at order10 is125: minimizing H and
minimizing a strong deletion increment are distinct optimization problems.

## 5. Boundary and significance

The first six minima, n=4,...,9, match 2F_(n-2), but at n=10 the chain
minimum is46 instead of42. The separate complete rooted-increment census
proves that46 is also the unrestricted strong-deletion minimum at order10.
Equation (12) is proved only for interval chains. Its extension to all
strong deletion pairs is an open conjecture until separately proved or
falsified; finite agreement is not such a proof.

This is a growth law for induced vertex arrival, not a +/-2 arc patch and
not a spectrum-completeness theorem. It identifies a useful distinction:
count-neutral growth is singleton SCC padding, and the cheapest currently
certified cyclic growth is controlled by a local marked block plus a
multiplicative outside chain.

Independent checks use `python3 verify_chains.py` (optional first argument sets the maximum order; the default is 12). The all-order result rests on the proof above, not the finite checks.
