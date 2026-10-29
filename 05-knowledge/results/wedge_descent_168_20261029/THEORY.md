# Proofs and failure boundaries

**Status:** PROVED elementary statements W1--W7, subject to their explicit
hypotheses; W1 imports the independently audited Session103 triangularity
lemma. This packet has same-agent proof and independent-algorithm checks,
not an independent audit or Lean verification. Classical connections are
identified as such; no literature-first claim is made.

## W1. Fixed-port nilpotence and exact recurrences

Let T be a loopless finite digraph on n vertices with a marked directed
Hamilton cycle C and successor s. Define D_C by putting s(u)->v for every
T arc u->v outside C. Suppose D_C is acyclic and fix 1<=k<=n.
Let M_k count spanning covers
by k ordered nonempty vertex-disjoint directed paths, with rows their ordered
starts and columns their ordered ends. Singleton paths are permitted.

For an ordered start tuple a, cut C before every member of a. Order the
resulting paths by their starts and call their ordered end tuple Phi(a).
The inverse construction cuts C after a supplied end set and matches each
end to the start of its segment, so Phi is a bijection.

Set B[a,b]=M_k[a,Phi(b)]. Choose a topological enumeration f:V->{0,...,n-1}
of D_C and put w(a)=sum_i f(a_i). For every counted cover with ends t,

    sum_(used u->v) (f(v)-f(s(u)))
        = sum_i f(s(t_i)) - sum_i f(a_i).

Cycle arcs contribute zero and all other arcs contribute a positive integer.
The unique all-cycle cover has t=Phi(a) and count1. Also
sum_i f(s(Phi(b)_i))=w(b). Thus B[a,a]=1 and every off-diagonal positive
entry satisfies w(b)>w(a). This is the inherited Session103 mechanism.

The least and greatest possible tuple weights differ by

    [(n-k)+...+(n-1)] - [0+...+(k-1)] = k(n-k) = h.

Consequently N=B-I has N^(h+1)=0: a contributing chain of h+1 off-diagonal
entries would raise weight by more than the full available range. It follows
over the integers that, for every nonnegative integer d,

    B^d = sum_(j=0)^h binom(d,j) N^j,
    B^-1 = sum_(j=0)^h (-N)^j,
    (I-zB)^-1 = sum_(j=0)^h z^j N^j/(1-z)^(j+1).       (W1)

Every fixed entry of B^d is therefore polynomial in d of degree at most h;
its (h+1)-st forward difference vanishes. Actual degree can be lower. Because
N is nonnegative, its nilpotence index equals one plus the longest path in
its support. This supplies an independent combinatorial test. For the
grammar word `ca`, n=6,k=3, the computed index is10, so degree9 is attained.
We do not claim equality for every n,k or every word.

This saves work when repeatedly composing one already constructed table:
precompute its nonzero N powers, then use binomial coefficients for any
depth. It does not remove the cost of constructing the table or yield a
fixed-size stationary recurrence for all A000568 tournaments. The indexed
carrier still has n!/(n-k)! ordered port states.

If an arc is weighted by q^(f(v)-f(s(u))), every full cover has the same
total exponent w(b)-w(a) for fixed boundary data. Hence the resulting B(q)
is D(q)^-1 B(1) D(q), where D(q)[a,a]=q^w(a). Such a potential weight is
only a change of coordinates. It cannot recover an erased internal path
clock. Counting actual physical steps or noncycle arcs needs a separate
variable or sidecar.

The all-order `be`/`ee` classification applies to the specific classical
unique-Hamiltonian-cycle grammar in Session103. It is not a classification
of all tournaments. Its singular direction is different: under the
bridge-sink hypotheses, the integral vector `(e_a-e_s) wedge e_r` is in
ker M2. In our `ee` control, `(a,s,r)=(0,2,1)` works directly in the native
ordered-pair coordinates. Potential telescoping gives the positive case;
an antisymmetric boundary relation gives this negative case.

## W2. A decisive test for closure under a menu of port operations

Let B_i=I+N_i be any finite family of m-by-m nonnegative real matrices,
with N_i nonnegative and diagonal zero. Let D be the union of the supports
of the N_i. The following are equivalent:

1. D is acyclic.
2. Some ordering makes every B_i unit upper triangular.
3. Every finite product of the B_i is unipotent.

Here unipotent means all eigenvalues equal1, equivalently a suitable power
of the matrix minus I is zero. These conditions concern the **same labelled
port carrier**, with the same identification at every seam.

Proof: a topological order gives2 from1, and triangular products give3.
Conversely, if D has a directed cycle v0->v1->...->vr=v0, choose a matrix
for each of its r edges. In their ordered product the all-diagonal choice
contributes1 to every diagonal entry; the chosen cycle contributes strictly
positive extra weight to the v0 diagonal. The trace is greater than m,
contradicting unipotence. A simple cycle has r<=m. Thus failure has a finite
constructible word witness, not just an asymptotic obstruction. QED.

If the common support has longest path length h, every product of more
than h matrices selected from the N_i is zero. For arbitrary depth d,

    B_(i1)...B_(id)
      = I + sum_(r=1)^h sum_(1<=j1<...<jr<=d) N_(ij1)...N_(ijr). (W2)

Order is retained: the N_i need not commute. A bounded finite matrix menu
then has polynomial growth in d of degree at most h. Independent nonnegative
weights on off-diagonal transitions preserve the same support argument.
Negative weights can cancel cycle contributions, so the converse cannot be
inherited after replacing native counts by arbitrary signed sectors.

**Native hostile.** For the strong four-vertex seed and k=2, B has size12.
Let P be the permutation of ordered pairs induced by swapping vertices1,3.
Then tr(B P^-1 B P)=22. Each factor separately is unipotent with determinant1.
All diagonal entries of the product are at least1 and at least one exceeds1;
repeated multiplication therefore has an exponentially growing diagonal.
Exhausting all24 vertex re-markings gives product trace12 in17 cases,16 in6,
and22 in1. Trace12 alone does not certify an acyclic union. Our first probe,
swapping0,1, had trace12 and did not supply a counterexample; completing the
24-case universe found the stated witness.

Even the 2-by-2 pair of opposite shears illustrates the same mechanism:
[[1,1],[0,1]] times [[1,0],[1,1]] equals [[2,1],[1,1]]. This small example
is explanatory; the12-by-12 witness is the actual tournament-observer test.

This is a concrete test for a proposed small attachment carrier. First
authenticate the attachment as an operator on the named native ports.
Then form the support union and either produce a common integer rank or
return a cycle word. A summary of determinants, traces of individual
operators, or current automorphism orders does not perform that test.

## W3. The actual K4--Fano--168 map and the different six-edge lattices

Let E={x in F2^4: sum x_i=0}, and p=1111. Its dimension is3; its nonzero
vectors are p and the six e_i+e_j indexed by K4 edges. Every permutation
of the four coordinates fixes p and acts faithfully on E. On the other
hand |GL(E)|=(8-1)(8-2)(8-4)=168; it is transitive on its seven nonzero
vectors. Thus its p stabilizer has order24, exactly this S4. The Fano lines
through p are {p,e_i+e_j,e_k+e_l}, the three opposite-edge matchings.

The classical exceptional isomorphism GL(3,2) ~= PSL(2,7) is CITED from
Brown--Loehr, *Why is PSL(2,7) isomorphic to GL(3,2)?* (2009), DOI
10.4169/193009709X460859. Session102 independently checks the source-adjoined
Paley7 **switching** group as PSL(2,7). These are different permutation
actions of isomorphic groups. The ordinary Paley7 automorphism group is21;
the stabilizer of the adjoined projective point in the8-point action is21,
whereas the stabilizer of p in the7-point Fano action is24. They cannot be
identified as the same marked action.

Over Q, oriented edges are the basis e_i wedge e_j of Lambda^2 Q^4,
of dimension6. The boundary map sends e_i wedge e_j to e_j-e_i. Its image
is the3-dimensional zero-sum subspace W, and its kernel is Lambda^2 W,
also dimension3. With the standard inner product, gradients give the
orthogonal cut complement. A cross product in oriented Euclidean3-space
is the additional identification *:Lambda^2 W -> W; wedge products
themselves need no such metric/orientation identification. A bivector
field is a varying section of this space; no arithmetic dynamics follows
merely from its dimension.

The integral cut and cycle sublattices of the oriented K4 edge lattice have
index16 and quotient (Z/4)^2. Use the three nonroot vertex gradients and
the three triangles through vertex0 as bases; their combined Smith diagonal
is (1,1,1,1,4,4), reproduced here. This is the classical K4 critical-group
lattice calculation. In contrast THM-3045's **unsigned permutation lattice**
has rational decomposition1+[22]+[31] and integral index24, with Smith
(1,1,1,2,2,6). They are not conflicting answers to one integral problem.
Neither gives the four spurious Z/6 factors of the separately defined T5
three-run gluing problem simply by counting edges.

Likewise the 3/4/6 minima of the grammar are different observations:
linear counting, accept/reject, and acceptance plus legality. Holding its
five port names and three output values fixed, the six-state machine has
exactly one automorphism. Each named letter has image size3, so it is not
a permutation of six K4 edges. This rules out an operation-intertwining
identification with the full S4 edge action. A bare six-element bijection
would preserve none of this structure.

## W4. Reflection unifies the two arithmetic rows, with an exact cubic boundary

For any real c and offsets d_1,...,d_k,

    c + sum_i [(c-2d_i)-(c+2d_i)] = c-4 sum_i d_i,
    c^2 + sum_i [(c-d_i)^2-(c+d_i)^2] = c(c-4 sum_i d_i). (W4)

Expansion proves both identities. For positive offsets, the positive
balancing center is c=4 sum d_i. All square bases are then positive.
The doubled linear bases also remain positive. The extra solution c=0
in the square equation is the symmetric signed-base degeneracy.

With d_i=i, c=2k(k+1). The linear row starts at2k^2 and advances by2;
the square bases start at k(2k+1) and advance by1. In each row there are
k+1 terms on the left and k on the right. Their common totals are

    linear: L_k=k(k+1)(2k+1),
    square: S_k=k(k+1)(2k+1)[12k(k+1)+1]/6.              (W5)

Thus Delta^4 L_k=0 and Delta^6 S_k=0, exact all-order recurrences. The first
square totals are25,365,2030,7230,19855,45955. The user's rows are the first
three instances. These square identities are classical; Klinger's2008
page explicitly lists the first two, so no discovery claim is made.

For cubes put s=k(k+1)/2. Consecutive offsets give

    C_k(c)=c^3+sum_(j=1)^k [(c-j)^3-(c+j)^3]
          =c^2(c-6s)-2s^2.                              (W6)

If c>0 is a root then c>6s. At6s the expression is negative; at6s+1 it is
(6s+1)^2-2s^2>0. Its derivative3c(c-4s) is positive for c>6s. Hence the
unique positive root lies strictly between6s and6s+1, so is not an integer.
There are **no** two adjacent blocks of k+1 and k consecutive positive cubes
with equal sum, for any k>=1. Other block-size differences are not excluded:
3^3+4^3+5^3=6^3 is a3-versus1 identity. The square proof closes on the
first offset moment; cubes introduce the third moment. This identifies
the first failed implication in extending the same balance law.

## W5. 169 states: finite-field norm data versus a split-ring trace clock

For every odd prime power q, on F_(q^2)^* consider

    x -> (a,b)=(x^(q+1), x^(q-1)).

The targets lie in F_q^* and the norm-one group respectively. Every image
pair obeys a^((q-1)/2)=b^((q+1)/2). The kernel is exactly {1,-1}, since
gcd(q-1,q+1)=2. Each target group is cyclic of even order, and its half-power
sign is surjective onto{1,-1}; exactly half of the target pairs have equal
signs. Their number is(q^2-1)/2, equal to the image size, so the compatibility
condition is also sufficient. One additional choice within{+x,-x} recovers
the source. Without such a choice this is not an injective chart.

At q=13, the two clocks have orders12 and14, their compatible pair count
is84, and the source count is168. We exhaust all elements in
F13[t]/(t^2-2); 1+2t has exact order168. This cyclic group is abelian,
whereas PSL(2,7) is nonabelian. Their common order alone gives no group
isomorphism or arithmetic action on the Paley switching class.

The older g58 example has a **different**169-element source:
O/13 with O=Z[u], u^2=5u+1 and u=(5+sqrt29)/2. Here29 is square mod13,
so O/13 is split F13 x F13, not F169. The inherited trace-carrier theorem
collapses169 states to13 under the sampled u^12 clock; adding the primitive
u clock recovers169. At11 the same unit is inert and u^12=-1; at13 it is
split and u^12=1. The weighted-CRT packet derives the exact lift orders.

This addresses the Zagier/12/13 annotation: g58^24=u^12 is the actual
sample supplied by the level-two modular formula. Its prime-dependent
trace loss and the new norm/ratio loss both require retained clocks, but
they have different kernels and different source rings. The prior packet
also separates Zagier's singular-moduli formula from his three-term
Apéry-type recurrences; neither is being identified by the integer12 alone.

## W6. Doubling, squaring, and an exact ABC prime budget

The summand shadow x->z with distinct positive summands is all x<z except
z=2x. Its missing arcs form doubling rays. In the proper-divisor shadow,
deleting equal factors removes x->x^2, giving squaring rays. This inherited
equal-parent deletion theorem is a precise additive/multiplicative duality.
Full two-parent generation still needs both parents; a shadow edge alone
does not certify that the companion has been generated.

For a primitive ABC triple a+b=c, 0<a<b and gcd(a,b)=1, define

    g=gcd((b-a)^2,4ab),
    (A,B,C)=((b-a)^2,4ab,c^2)/g.

Then g=1 if c is odd and g=4 if c is even. This is another primitive positive
ABC triple, and

    rad(ABC)=rad(abc) rad(oddpart(b-a)).                  (W7)

Proof: a and b are coprime; every odd prime of b-a is coprime to abc,
since gcd(b-a,c) divides2. All odd primes dividing abc remain present in
ABC, and no others enter except those of b-a. In the even-c case, a,b are
odd and (c/2),(b-a)/2 have opposite parity, so2 remains in the new radical
after division by4. This proves the formula and the gcd normalization.

The height becomes c^2 or c^2/4. It strictly increases except at the sole
case c=4,a=1,b=3, where the triple is fixed. Examples:

    1+2=3 -> 1+8=9 -> 49+32=81;
    radicals: 6 -> 6 -> 42;
    1+3=4 -> 1+3=4.

This exact prime budget explains both the initial radical coincidence and
the next new prime7. It is not an ABC inequality and gives no universal
bound on the new prime budget. A descending inverse must have an explicit
domain and a proof that every target enters it. Neither is supplied by a
formal reversal of the identity.

## W7. What a natural-number descent proof must still provide

For a deterministic transition f and a genuine absorbing target set R,
universal finite hitting is equivalent to the existence of a rank
rho:X->N with rho=0 on R and rho(f(x))<rho(x) off R. Sufficiency is
well-ordering; necessity uses the hitting time. The latter construction
is circular if termination has not already been proved.

A useful sufficient arithmetic form instead supplies, for every source
n outside a checked finite floor, a **legal finite return word** ending at
some 0<m<n. Strong induction then proves termination when the floor's
targets are independently grounded. For odd Collatz this is exactly the
inherited first-descent receipt interface. Our999 finite receipts preserve
each valuation exponent and both source and endpoint; they are regression
controls, not a proof of all-source coverage.

For a fixed authenticated affine return n'=(A n+C)/D with D>A and C>=0,
descent holds precisely when n>C/(D-A). This turns a proposed rule into
a decidable threshold plus a residue/valuation guard. A finite or recursive
rule library still needs a coverage theorem. Guard satisfaction, source
identity, physical clock, and arrival in the actual rooted component are
separate proof obligations.

The exact expanding prefixes are a hostile to bounded corrections. For
L>=1, n=2^(L+1)-1 has L odd accelerated Collatz steps of valuation1:

    U^j(n)=3^j 2^(L+1-j)-1,        0<=j<=L.

The endpoint ratio grows without bound as L grows. Thus log n plus a
globally bounded correction cannot decrease on every step. The much
stronger THM-4507 excludes arbitrary corrections built from finitely many
polynomial valuations, even with finite colours. That theorem leaves
selected-return ranks and unbounded fibre information open. Different L
use different positive sources: these are not one divergent positive orbit.

For5x+1, the odd cycle13->33->83->13 rules out any strictly decreasing
rank for every nonterminal native step. Passing to bijective coordinates
cannot remove this cycle. Any claimed transfer to3x+1 must exhibit an
arithmetic hypothesis that this cycle fails, not merely a new state name.

Finally, a positive atomic source law nu(n)>0 detects each nonterminating
source. If S_t is the set surviving to time t, then nu(S_t)->0 iff every
source terminates (continuity from above). This repairs the possibility of
discarding individual integers in an average; it does not prove that the
survivor mass tends to zero. For LRC the analogous endpoint-atomic repair
detects isolated target points but still needs a uniform positivity proof.
Topology, dimension counts, and an average over the wrong source measure
cannot substitute for these missing inequalities.
