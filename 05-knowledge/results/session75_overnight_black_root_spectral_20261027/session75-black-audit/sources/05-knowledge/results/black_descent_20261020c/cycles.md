# Cycle-changing black descent

Status: PROVED for Sections 2–7; FINITE-EXACT controls in Section 8.
Strong black descent is proved at every order. Black connectivity remains OPEN.
Scope: converse-merged complement-line metagraph, n>=3.
This is the Hamiltonian-path reversal operation, not single-arc reversal.

## 1. Inheritance and target

The inherited [black closed-path theorem and blue-cut coset lemma](../blue_quotient_20261020b/structure.md)
show that every strong black-active tournament has a black closed
Hamiltonian path. The exceptions among strong tournaments are exactly C3
and the regular tournament R5, which are pure blue. The same note records
that endpoint-current averages and equal-endpoint existence alone do not
give the desired descent.

The [endpoint transport theorem](../merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md)
gives the exact current below. The [uniform prefix rerouting proof](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md)
is a dependency of the inherited black closed-path theorem. These are
proved dependencies, not an inference from the finite census.

Concept board: (i) black components and score descent; (ii) cycle cuts and
anti-automorphism cosets; (iii) prime rotation orbits; (iv) additive
boundary compatibility; (v) ordered SCC normalization. The new operation
uses (iii) and (iv) to repair the fixed-cycle obstruction in (ii).
It does not itself finish (v), hence does not prove black connectivity.

For a tournament T write s(v) for outdegree and c3(T) for the number of
directed triangles. A Hamiltonian path P=(a,...,b) is **closed** when b->a.
It is **blue** when its positional reversal is an anti-automorphism of T,
and **black** otherwise. Let T^P reverse all arcs of P. Then

    c3(T^P)-c3(T) = s(a)-s(b)-1.                         (1)

Only a and b change outdegree, respectively by -1 and +1; (1) also follows
immediately from c3(T)=binom(n,3)-sum_v binom(s(v),2).
Thus a black closed P with s(b)>=s(a) gives a strict black descent.
Path reversal represents the same merged line as complementing nonpath
arcs, by [THM-4729, Proposition K](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md).
The directed-triangle count is invariant under converse, so this descent
survives the converse quotient.

## 2. An additive boundary lemma, with a self-contained proof

We use the classical Cauchy–Davenport inequality over the prime field F_p:
for nonempty A,B,

    |A+B| >= min(p, |A|+|B|-1).

Here is a proof so that the argument does not depend on an unexamined
literature claim. Put a=|A|, b=|B|. If a+b>p, every translate x-B meets A,
so A+B=F_p. Otherwise, if |A+B|<=a+b-2, extend A+B to a set C of size
a+b-2 and form

    f(X,Y)=product_(c in C)(X+Y-c).

This polynomial vanishes on A x B. The coefficient of X^(a-1)Y^(b-1)
is binom(a+b-2,a-1), nonzero in F_p because a+b-2<p.
On the other hand, for every polynomial of total degree at most a+b-2,
that coefficient is

    sum_(x in A,y in B)
       f(x,y) /
       [product_(x' in A\{x})(x-x') product_(y' in B\{y})(y-y')].

To verify the identity, apply the two univariate Lagrange coefficient
functionals to each monomial. Each functional annihilates powers below
its set size minus one and takes that critical power to 1. The total
degree bound leaves only the stated monomial. The displayed sum is zero,
a contradiction. The case a=b=1 is immediate.

Induction gives the consequence we need: for nonempty S_0,...,S_(ell-1),

    |S_0+...+S_(ell-1)|
        >= min(p, sum_i |S_i|-ell+1).                    (2)

## 3. Prime-orbit block normal form

**Lemma (PROVED).** Let C=(v_0,...,v_(n-1)) be a directed Hamiltonian
cycle of T, where n=p ell, p is an odd prime, and ell>=3. Suppose the
rotation v_j -> v_(j+ell) is an automorphism. Then T has a directed
Hamiltonian cycle in which each of this automorphism's ell vertex orbits
occurs as one consecutive block.

Write

    u_(i,j)=v_(i+ell j),  i in Z_ell, j in F_p.

The i-th orbit is a circulant tournament on p vertices. Choose
epsilon_i in {+1,-1} so that u_(i,j)->u_(i,j+epsilon_i).
Because p is prime, for every starting coordinate t_i the word

    u_(i,t_i), u_(i,t_i+epsilon_i), ...,
    u_(i,t_i-epsilon_i)

is a Hamiltonian path through that orbit.

For each pair of consecutive orbit labels in cyclic order, define

    S_i={d in F_p: u_(i,j)->u_(i+1,j+d)}.

Rotation invariance makes this independent of j. The given Hamiltonian
cycle makes every S_i nonempty: for i<ell-1 it contains 0, and the last
set contains 1.

If any S_i=F_p, the sum of these nonempty sets is already F_p. Otherwise
all S_i are proper, and we may reverse the cyclic order of the orbit
blocks. Its cross sets are, in reverse order,

    -(F_p \ S_i).

They are all nonempty, and the total cardinalities in the two directions
add to p ell. Choose the direction whose total is at least ceil(p ell/2).
For p>=3 and ell>=3,

    ceil(p ell/2) >= p+ell-1.

For the smallest pair (3,3) both sides are 5. Otherwise
(p-2)(ell-2)>=2 gives the inequality even without rounding.
Therefore (2) says the sum of the chosen cross sets is all of F_p.

Choose d_i in the cross set for the chosen cyclic block order so that

    sum_i d_i = sum_i epsilon_i   in F_p.               (3)

Start with t_0=0 and recursively set

    t_(i+1)=t_i+d_i-epsilon_i.

Equation (3) makes the coordinates close cyclically. The endpoint of
block i has coordinate t_i-epsilon_i, so the cross difference from this
endpoint to the next block's start is

    t_(i+1)-(t_i-epsilon_i)=d_i.

Every boundary arc, including the final closing arc, is therefore
present. Concatenating the orbit paths proves the lemma.

This proof is constructive: compute cross sets, choose the better block
direction, solve the finite sumset equation (3), and recover all starts.
It requires no assertion that every tournament edge lies on a Hamiltonian
cycle; that stronger assertion is false.

**Score consequence.** If one rotation orbit has score a with
2a!=n-1, the resulting cycle has an internal orbit edge whose cut path P
has equal endpoint scores a and is black. Indeed a blue path reverses
outdegrees, requiring the endpoint scores to sum to n-1. Equation (1)
then gives c3(T^P)=c3(T)-1.

## 4. All odd strong black-active tournaments descend

**Theorem (PROVED).** Let n>=3 be odd. Every strong black-active
tournament on n vertices has a black closed Hamiltonian path P satisfying
s(last P)>=s(first P). Its merged black neighbor T^P has strictly fewer
directed triangles.

If T is regular, the inherited black closed-path theorem supplies P, and
its endpoint scores are equal. Assume henceforth T is nonregular.
Choose any Hamiltonian cycle C, and let rho be its one-step rotation.
For its cut beginning at v_j the reversal is

    r_j(v_i)=v_(2j-1-i).

Let B be the indices of blue cuts. If B is nonempty, the inherited coset
lemma, or the identity r_j r_i=rho^(2(j-i)), gives

    B=i+H,  H={h in Z_n: rho^(2h) is an automorphism}.    (4)

In particular B is a coset of a subgroup. If B is empty, the endpoint
score differences s(v_(j-1))-s(v_j) sum to zero, so some cut is a black
descent. If B is all of Z_n, rho^2 is an automorphism. Since n is odd,
rho is an automorphism, making T regular, contrary to assumption.

Suppose C has no descending black cut. Every cut outside B then obeys

    s(v_(j-1)) < s(v_j).                                (5)

If |B|=1, (5) gives n distinct strictly increasing scores along the
Hamiltonian word obtained by cutting the unique blue boundary. This is
impossible: strong tournament scores lie in {1,...,n-2}.

Thus b=|B| is an odd divisor of n with b>=3. Choose an odd prime p|b.
Equation (4) gives rho^(2n/p) as an automorphism. Its order is p, so taking
a suitable power gives rho^(n/p) as an automorphism. Put ell=n/p. If
ell=1 this rotation is transitive and T is regular. Otherwise oddness
gives ell>=3, and the prime-orbit block lemma applies.

Scores are constant on its rotation orbits. Nonregularity gives an orbit
whose score a differs from (n-1)/2. The score consequence of Section 3
then supplies an equal-endpoint black closed path, a contradiction.
This proves the theorem.

The theorem asserts existence of a descending black move from each
strong black-active node. It does not assert that the descending neighbor
is strong, nor that different nonstrong terminal basins have already
been connected.

## 5. The even-order obstruction reduces to two blue cuts

**Proposition (PROVED).** For even n, if a Hamiltonian cycle C has either
no blue cuts, all blue cuts, or more than two but fewer than n blue cuts,
then T has a black closed path with nonincreasing endpoint current, hence
strict triangle descent.

With no blue cuts the averaging argument applies. With all cuts blue,
the explicit two-orbit construction in
[the inherited theorem, even-order proof](../blue_quotient_20261020b/structure.md)
gives a black closed path with equal endpoint scores.

For a nonempty proper B, write |B|=2d, where d is an odd divisor of n/2;
this follows from (4) because tournament automorphisms have odd order.
If |B|>2, choose an odd prime p|d. As above rho^(n/p) is an automorphism,
and ell=n/p is even. If ell=2 then d=p and |B|=n, excluded here.
Thus ell>=4, and Section 3 applies. Every orbit score is noncentral,
because 2a cannot equal odd n-1. Its internal edge gives the descent.

Consequently a counterexample to the all-order strong descent assertion
would have even order, and **every Hamiltonian cycle would have exactly
two blue cuts**, with strict score increases at every other cut. This
is a substantially narrower obstruction than arbitrary self-converseness.
The next section excludes it.

## 6. Terminal-block insertion closes the even-order case

**Theorem (PROVED).** Every strong black-active tournament of even order
has a black closed Hamiltonian path P with s(last P)>=s(first P).

Suppose otherwise, and write n=2m. By Section 5, every Hamiltonian cycle
has exactly two opposite blue cuts and every other cut has strictly
increasing scores. Choose a Hamiltonian cycle

    C=(a_0,a_1,...,a_(n-1))

whose first vertex v=a_0 has minimum score delta. Its cut at v has
nondecreasing endpoint score, so under the contrary assumption it must
be blue. Put u=a_(n-1). Blue reversal gives

    s(u)=Delta=n-1-delta,

the maximum score. The two blue cuts are at 0 and m. The scores
strictly increase within each of the two blocks of length m. Therefore

    Delta-delta >= m-1,
    delta <= floor(m/2).                                (7)

Assume first m>=3. Among the vertices a_1,...,a_(n-3), choose the last
one y=a_(k-1) that beats v. Such a vertex exists: v has n-1-delta
inneighbors, and excluding u and a_(n-2) still leaves at least one,
by (7). Let

    Q=(x=a_k,...,u),  t=|Q|.

The choice of y excludes both final positions, so t>=2. Every vertex
of Q except u and possibly a_(n-2) is beaten by v. Also v beats a_1,
which lies outside Q because y cannot be a_1. Counting outneighbors
of v gives

    t-2 <= delta-1,
    2 <= t <= delta+1 < m.                              (8)

The final strict inequality follows from (7) and m>=3. In particular
y and all of Q lie inside the second increasing-score block, so

    s(y)<s(x),  Delta-s(x)>=t-1.                         (9)

The indegree of x is

    n-1-s(x)=delta+Delta-s(x).

At most t-2 vertices of Q beat x, because x beats the next vertex of Q.
The vertex y accounts for one further inneighbor. Hence at least

    delta+Delta-s(x)-(t-2)-1 >= delta

vertices a_j with j<=k-2 beat x. Their successors a_(j+1) are distinct.
The indegree of u is delta, and its predecessor a_(n-2) belongs to Q
and already beats u. Consequently at most delta-1 of these distinct
prefix successors beat u. There is therefore an index j<=k-2 with

    a_j -> x,  u -> a_(j+1).                            (10)

Insert the intact terminal block Q between a_j and a_(j+1). The new word

    v,...,a_j, x,...,u, a_(j+1),...,y

is a Hamiltonian path by (10), and it closes because y->v.
Its first vertex has minimum score delta; its last has
s(y)<Delta by (9). Thus

    s(y)>=delta,   s(y)+delta<n-1.

The first inequality gives strict triangle descent. The second excludes
blue reversal, which would force the endpoint score sum to be n-1.
This contradiction settles m>=3.

For n=4, strong scores lie in {1,2}; their sum is 6, so there are two
of each. A contrary two-block cycle must have score word (1,2,1,2).
Each score-1 vertex already beats its following score-2 vertex on the
cycle, while one of the two score-1 vertices must beat the other.
That vertex would have score at least 2, a contradiction.

Together with Section 4 this proves:

**All-order strong descent theorem.** For every n>=3 and every strong
black-active tournament T, some black closed Hamiltonian path P satisfies
s(last P)>=s(first P), and therefore c3(T^P)<c3(T).

The even proof is constructive after choosing a hypothetical obstructing
cycle: the terminal block is found from the last eligible inneighbor of
the minimum-score start, and the two boundary arcs follow by pigeonhole.
The extra coordinate that the original fixed-cycle method lacked is
the movable terminal block together with its two endpoint degree counts.

## 7. Deleting an extremal-score vertex gives a separate short criterion

**Proposition (PROVED).** Let T be strong of order n>=4. If a
minimum-outdegree vertex v has T-v nonstrong, then every Hamiltonian
cycle cut beginning at v is a descending black path.

Non-self-converse T has only black paths, and its minimum-score start
already guarantees descent. Suppose T is self-converse, so

    minimum score + maximum score = n-1.                (6)

Write the strong components of T-v as C_1->...->C_r, where r>=2.
Every Hamiltonian word of T-v visits these components in that order.
Thus the predecessor u of v on any Hamiltonian cycle belongs to C_r.
Put t=|C_r|. If t>1 then C_r is strong of order at least 3, so u has at
most t-2 internal outneighbors and at most one additional outneighbor v.
Consequently s(u)<=t-1. Any vertex of C_1 beats all t vertices of C_r,
so the maximum score is at least t. If t=1, u has precisely one
outneighbor v, whereas the maximum score is at least ceil((n-1)/2)>=2.
In both cases s(u) is strictly below the maximum score.

Equation (6) now gives s(u)+s(v)<n-1, excluding a blue cut. Since v is a
minimum-score vertex, s(u)>=s(v), giving the required descent.

Applying this to the converse proves the analogous sufficient condition
for a maximum-score vertex whose deletion leaves a nonstrong tournament.
This gives a shorter constructive route whenever an extremal-score
deletion is nonstrong. Section 6 also covers the other even-order cases.

The following stronger statement has only finite support here:

> Every nonregular strong self-converse tournament has a Hamiltonian
> cycle with an edge u->v where v has minimum score and u has less than
> maximum score.

This statement is stronger than the descent theorem and remains OPEN;
neither arbitrary prescribed edges nor arbitrary equal-score edges are
known here to extend to Hamiltonian cycles.

## 8. Exact controls and explicit failed shortcuts

Reproduce with:

    python 05-knowledge/results/black_descent_20261020c/cycles.py

The [script](cycles.py) uses explicit exceptions rather than Python
assertions. Its [output](cycles.out) records:

| p | number of orbits | all nonempty cross-set systems | signed boundary certificates | literal tournament checks |
|---|---:|---:|---:|---:|
| 3 | 3 | 343 | 2,744 | 28 |
| 3 | 4 | 2,401 | 38,416 | 381 |
| 5 | 3 | 29,791 | 238,328 | 2,360 |

Every sign vector in {+1,-1}^ell is checked for each cross-set system.
The certificate solver independently uses finite sumset dynamic
programming and verifies each boundary equation. The sampled literal
tournaments additionally check every constructed cycle arc and the
rotation automorphism. Cross sets are allowed to be arbitrary nonempty
subsets: this tests a larger universe than the original-cycle coordinates
with specified 0/1 entries.

The complete class controls read the inherited authenticated n8 class
table and the [independently authenticated n9 census](../blue_quotient_20261020b/census.md).
They test all strong self-converse classes; non-self-converse classes are
covered by the elementary all-black cycle-average proof.

| n | self-converse classes | strong self-converse classes | a descending black cut appears by |
|---|---:|---:|---|
| 8 | 176 | 162 | first rooted cycle: 161; second: 1 |
| 9 | 2,752 | 2,661 | first rooted cycle: all 2,661 |

The [literal certificate file](cycles_certificates.json) contains one
Hamiltonian path and triangle change for each of these 2,823 classes.
The cycle ranks depend on the canonical labeling and DFS order, and are
not asserted to be invariant or to persist at larger orders.
Each certificate is checked by actually reversing its path arcs and
independently counting directed triangles before and after.

The terminal-block construction in Section 6 is also implemented
literally. The checker examines all 30 rooted Hamiltonian cycles in the
strong self-converse n6 classes and all 4,284 in the n8 classes. There
is respectively one and nine two-blue-cut cycles with no descent among
their existing cuts. The construction repairs all ten, with every new
arc, black reversal predicate, and actual triangle change checked.
One n8 witness uses a block of length 3; the others use length 2.
The two additional bad n6 cycles have all six cuts blue, and are handled
by the inherited two-orbit construction rather than fed to this routine.

There are fixed-cycle failures already at n=6. Its adjacency bitsets are

    [2,12,25,17,35,15],

with scores [1,2,3,2,3,4]. The natural cycle (0,1,2,3,4,5) has blue cuts
exactly at positions 0 and 3, and none of its black cuts decreases
triangles: all four have triangle current zero.
The new cycle (0,1,2,4,5,3) has no blue cuts; its cut beginning at 0
decreases triangles by 2. Analogous n8 and n9 fixed-cycle failures and
their repaired paths are included in the output. This refutes a
fixed-cycle proof, not the cycle-changing theorem.

The shortcut "an even strong self-converse tournament always has a
closed path with equal endpoint scores" is also false. At n8, five
complete-table classes have no such path at all:

    68182161, 68176401, 8921188, 34624274, 33702274.

These integers use the n8 table's upper-triangle, most-significant-bit
orientation encoding; they must not be confused with the different
tiling canonical codes in the older structural note. The script
exhausts every rooted Hamiltonian cycle in each negative case.
All five nevertheless have strict black descent. Thus the surviving
target is s(last)>=s(first), not compulsory equality.

The extremal-predecessor statement at the end of Section 7 was separately
checked on all strong self-converse classes at n=4,6,8 and all 2,646
nonregular strong self-converse classes at n=9, with no failure. This is
only a finite signal for that stronger statement; the all-order descent
proof uses Sections 2–6.
