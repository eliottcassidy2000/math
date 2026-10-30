# Integer Ramsey capacity from finite parity charts

2026-10-30. **PROVED scoped deductions / CITED background / FINITE-EXACT
certificate banks.** Same-agent review, with separate exact computations.
No claim of literature priority. In particular, the exponent below is
already implied by Liu's stronger U5-free theorem; the contribution here
is an integer refinement, a finite proof compiler, and explicit algorithms.

## R1. The complete five-slot weighted polytope

Let T5 have vertices Z/5 and i -> j iff j-i is 1 or 2 modulo5. Write
tau(T) for the maximum order of an induced transitive subtournament.
The maximal transitive supports of T5 are

    012, 123, 014, 034, 234.                         (1)

The inherited substitution identity is

    tau(Q[T_i]) = max_{S transitive in Q} sum_{i in S} tau(T_i).

This is P3 of
[tournament_open_problem_routes_20261029/PROOFS.md](../tournament_open_problem_routes_20261029/PROOFS.md),
not a new identity. It holds with empty slots: a transitive set meets only
a transitive quotient support, and maximum child sets on such a support
can be united.

Define P={w>=0: sum_{i in S}w_i<=1 for each S in(1)}. Its vertices are
exactly zero, the five coordinate unit vectors, the five vectors equal
to 1/2 on a cyclic triple and zero elsewhere, and (1/3,...,1/3).

**Proof of completeness.** A vertex with one positive coordinate is a
unit vector. Two positive coordinates lie in a transitive support, so
the corresponding face is a simplex. For three positive coordinates,
a transitive support again gives a simplex; a cyclic support gives the
three pair constraints, whose sole fully positive vertex is (1/2)^3.
For four positive coordinates, rotate so the missing coordinate is4.
The only potentially tight nonredundant positive-support inequalities
are w0+w1+w2<=1, w1+w2+w3<=1, and w0+w3<=1. Three hyperplanes cannot
isolate four positive coordinates. With all five coordinates positive,
all five inequalities must be tight. Subtracting consecutive cyclic
three-term sums gives w_i=w_(i+3), hence all weights are1/3. Each listed
point is a vertex. Exact rational enumeration of all252 active bases
provides a separate check.

Set p=log_2(3). Convexity of sum w_i^p on P gives

    max_{S in(1)} sum_{i in S} w_i >= (sum_i w_i^p)^(1/p).       (2)

The half-triple and unit vertices have p-mass1. The uniform vertex has
mass5/3^p<1: p>3/2 follows from9>8, and 3^(3/2)>5 from27>25.
Strict convexity gives equality for nonzero w exactly on rays through
unit vertices or half cyclic triples. A cyclic triple with equal weights
forces p>=log_2(3), so this is the sharp all-weight exponent.

For every tournament built recursively from induced subtemplates of T5,
induction yields tau(T)>=|T|^(log_3(2)). The balanced five-child power
alone suggests the stronger exponent log_5(3); unequal cyclic-triple
children refute that inference. Even requiring every five-child slot to
be nonempty does not improve the asymptotic exponent: put three copies
of the old tree on the cyclic triple013 and singleton leaves on2,4.
Starting at a singleton, its order and transitive number are exactly

    N_d=2*3^d-1,     tau_d=2^(d+1)-1.                (3)

**Literature boundary.** Liu's Theorem6.1 proves the same exponent for
ALL U5-free tournaments, a larger class, and uses cyclic-triangle powers
for sharpness. We do not claim a new Erdős--Hajnal exponent.
[Primary source, 1208.0398v2](https://arxiv.org/pdf/1208.0398v2).

## R2. A parity certificate for every integer weight configuration

Let Q be a fixed tournament, and let M be its maximal transitive supports.
For every subset O of vertices and epsilon in{0,1}, seek a map
c:O->{0,1,2} such that, for every S in M and j in{0,1,2},

    |S intersect c^-1(j)| <= ceil((|S intersect O|-epsilon)/2)
                              + epsilon*[j=2].                (4)

This is a finite constraint system on the actual joint incidence of M.
The empty-subset case has all counts zero. A bank satisfying(4) is a
sufficient certificate, not a necessary characterization of all possible
Ramsey capacity inequalities.

Define G(0)=0,G(1)=1 and, equivalently,

    G(k)=sum_{0<=j<k} 2^popcount(j),
    G(2h)=3G(h),    G(2h+1)=2G(h)+G(h+1).            (5)

For every nonnegative integer vector a with transitive support sums<=k,
a complete bank(4) proves

    sum_i G(a_i) <= G(k).                            (6)

**Proof.** Write k=2h+epsilon and a_i=2b_i+e_i, and take O={i:e_i=1}.
For a support S with r=|S intersect O|, feasibility gives

    sum_{i in S}b_i <= h-ceil((r-epsilon)/2).

Using the certified coloring, define three integer layers
a_i^(j)=b_i+[i in O and c(i)=j]. Equation(4) makes their support
budgets respectively h,h,h+epsilon. Coordinatewise(5) gives the exact
mass identity

    sum_i G(a_i) = sum_{j=0}^2 sum_i G(a_i^(j)).       (7)

Induct on k. At k=0 all weights vanish; at k=1 every pair of vertices
is transitive, so at most one coordinate is positive and it is1.
For k>=2 all three layer budgets are strictly smaller than k.
Apply induction and(5) to(7), giving(6). This argument quantifies over
ALL integer weights; it is not an extrapolation from tested budgets.

### A human-readable complete certificate for T5

For epsilon=0, a triple containing r odd vertices allows at most
ceil(r/2) of each color. With |O|<=3, give its vertices distinct colors.
With |O|=4, the three transitive triples containing the unique even
vertex impose a three-edge PATH on the odd vertices; two-color that
path. This also prevents a monochromatic triple. With |O|=5, distribute
the colors in group sizes2,2,1.

For epsilon=1, a triple with exactly one odd vertex forces that vertex
to color2. If |O|<=2, use color2 everywhere. If |O|=3, a forced vertex
must complete a transitive triple with the same two even vertices.
Each pair of T5 belongs to at most two supports in(1), so at most two
odd vertices are forced. Put them in color2, put at most two vertices
there altogether, and distribute the remainder singly into colors0,1.
For |O|=4 use sizes1,1,2; no triple has just one odd vertex. For |O|=5,
put a cyclic triple in color2 and the remaining two vertices separately
in colors0,1. No transitive triple is wholly color2. All conditions(4)
follow, including odd triples of size2 or3.

The checker independently searches all64 charts. It also supplies
complete256-chart and1024-chart banks for the odd carousel tournaments
T7 and T9 (i -> j for 0<j-i<=(n-1)/2 modulo n). All maximal transitive
supports were separately recovered from the literal relations. Thus
(6) is **PROVED by a finite certificate** for each of T5,T7,T9, and for
their induced subtemplates by zero extension.

**Hostiles and boundary.** For U5 obtained by reversing edge0->1 of T5,
epsilon=1,O={2,3,4} has no chart. This refutes this sufficient compiler
for that template, not inequality(6). Paley7 gives a genuine obstruction
to(6): unit weights have k=3 but sum G(1)=7>G(3)=5. Its chart also fails
at epsilon=0,O={0,1,2,3}. No all-tournament assertion is made.

## R3. Exact capacity and a logarithmic arithmetic algorithm

Let C be the hereditary substitution closure of the templates T5,T7,T9
(or T5 alone). Include the empty tournament as a bookkeeping object.
Let F(k)=max{|T|:T in C and tau(T)<=k}. Then

    F(k)=G(k) for EVERY k>=0.                         (8)

**Upper bound.** Induct on a substitution tree. A child with budget a_i
has order at most G(a_i); the parent's actual budgets satisfy the
joint transitive-support constraints. Apply R2.

**Attainment by literal Pascal cells.** Let

    V_k={(j,t):0<=j<k, 0<=t<=j, binom(j,t) is odd}.

In characteristic2, (1+X)^j is the product of (1+X^(2^r)) over set bits
r of j. Thus a cell is present exactly when t is a bit-submask of j;
|V_k|=G(k). At each bit position its pair (j_r,t_r) is one of
A=(0,0),B=(1,0),C=(1,1). Orient between two distinct cells by their
LEAST significant differing bit-pair, using A->B->C->A. There are no
ties, and this defines the actual tournament P_k, not just a count.

Splitting by the lowest bit gives three cyclic children with budgets
(h,h,h) at k=2h, or (h+1,h,h) at k=2h+1. Induction using the exact
substitution rule proves tau(P_k)=k. These are permitted C3 substitutions
inside T5, so P_k lies in C and attains(8). The inclusions V_k into
V_(k+1) preserve every old arc. Adding Pascal row k adds exactly
2^popcount(k) new vertices. The bit-pair address is the alignment
coordinate; the integer k is the source cutoff, not a field phase.

Consequently the transitive-k Ramsey threshold RESTRICTED TO C is
G(k-1)+1. This is not a new unrestricted Ramsey number. At powers of
two the familiar values F(2^d)=3^d are recovered, while(8) settles
all intervening integer budgets too.

**Algorithm.** Read the binary digits of k from most to least significant,
starting with (a,b)=(G(0),G(1))=(0,1). On bit0 update
(a,b)=(3a,2a+b); on bit1 update (a,b)=(2a+b,3b). The final a is F(k).
This takes O(log(k+1)) arithmetic operations and O(log(k+1)) storage
bits, or O(log^2(k+1)) elementary bit operations with ordinary arithmetic.
An adjacency query on Pascal labels takes O(log(k+1)) bit steps; writing
the whole adjacency matrix still costs its output size.

**Independent optimization check.** Before deriving R2, we solved the
T5 integer attachment optimization directly:

    F(k)=max sum_i F(a_i),  0<=a_i<k,
                       a0+a1+a2,a1+a2+a3,a0+a1+a4,
                       a0+a3+a4,a2+a3+a4 <=k.        (9)

Unary roots can be contracted. At any nonunary root, every positive
child budget is <k because every pair is transitive. Attaching a
singleton proves F(k+1)>=F(k)+1. Thus the fifth coordinate can always
be set to min(k-a0-a1,k-a0-a3,k-a2-a3,k-1). The resulting four-loop
algorithm takes O(K^5) arithmetic operations for all budgets through K,
and stores attaining trees. It tested18,855,887 candidates through64,
with an independent five-loop optimization through6 and native witnesses.
R2 improves exact evaluation from that search to the bit algorithm.

**Existing sequence and recurrence literature.** G is A006046, not a
new integer sequence. Hwang--Janson--Tsai, Proposition C.3, with
(alpha,beta)=(2,1), solves the binary maximum recurrence that also arises
from arbitrary C3-only substitutions; they credit Knuth--Lossers.
The five-slot integer inequality(6) is the additional step here.
[Primary manuscript](https://arxiv.org/pdf/2210.10968v1),
[sequence record](https://oeis.org/A006046).

## R4. A two-sided ordinary-graph companion

Let D be the closure of singletons under disjoint union, join and
substitution into the undirected five-cycle C5. Write alpha,omega for
independence and clique numbers. For every G in D,

    |G| <= (alpha(G)*omega(G))^gamma, gamma=log_4(5). (10)

The exponent is sharp: balanced C5 powers have order5^d and
alpha=omega=2^d. This is a theorem about D, not arbitrary C5-free graphs.

**Proof.** The polytope x>=0, x_i+x_(i+1)<=1 on C5 has vertices the
eleven independent-set indicators (including empty) and (1/2)^5.
Indeed, a fractional vertex with a zero coordinate lies on path
components, whose alternating edge constraints force integral vertices.
With no zero, all five edges must be tight, giving the half vector.

For a C5 substitution let a_i=omega(G_i)/omega(G) and
b_i=alpha(G_i)/alpha(G). Then a lies in this edge polytope and b in its
complementary-cycle version. The separately convex function
sum_i(a_i*b_i)^gamma takes a maximum at a pair of vertices. Two
indicator vertices meet in at most one coordinate (a clique and an
independent set); one half vector and one indicator give at most
2*(1/2)^gamma<1; two half vectors give5*(1/4)^gamma=1. Therefore the
sum is<=1, and induction proves(10) for substitution. Union and join
follow from sum x_i^gamma <= (sum x_i)^gamma. All144 vertex pairs were
checked with exact rational coordinate products.

This companion uses two complementary incidence systems. The tournament
calculation uses one transitive-support system. Equal fivefold counts
do not identify the graphs, orientation gauges, or their operations.

## R5. Exactly what the Ramsey observer forgets

A labelled vertex set of a tournament is transitive iff it contains
no cyclic triple. Thus its labelled cyclic-triple hypergraph determines
ALL weighted transitive-subset queries, although it does not determine
the tournament's automorphism group or literal path complements.
The retained twelve-vertex81/27 pair is therefore a successful common
Ramsey observer and a hostile control for any claim of full reconstruction.
Its two automorphism orders are replayed; no new metagraph closure or
A000568 term follows from(8).

The Pascal map retains named bit-pair addresses and nested source cutoffs.
Forgetting them retains the scalar capacity but loses the actual arcs
and the transport needed for subsequent metagraph operations. Conversely,
retaining all arc data is unnecessary for the restricted Ramsey query.

## R6. Full symmetry and decoration recurrences on the extremal family

Let |Aut(P_k)|=3^a(k), and let O_k(q) count vertex q-colorings of P_k
modulo its FULL tournament automorphism group, with colors named. Then

    a(1)=0,
    a(2h)=1+3a(h),       a(2h+1)=2a(h)+a(h+1),

    O_1(q)=q,
    O_(2h)(q)=(O_h(q)^3+2O_h(q))/3,
    O_(2h+1)(q)=O_h(q)^2*O_(h+1)(q).                (11)

For clarity, the full cycle index, with formal cycle variables x_j, is

    Z_1=x_1,
    Z_(2h)=(Z_h^3+2 Z_h[x_j -> x_(3j)])/3,
    Z_(2h+1)=Z_h^2 Z_(h+1).                         (12)

**Proof of full-group scope.** In C3[A,B,C], any proper module is
contained in one child: a module meeting two children must contain the
third (an outside vertex there sees the two in opposite directions),
and then must contain the other children as well. Thus these three
nonempty children are intrinsically the maximal proper modules.
Every automorphism permutes them through a rotation of C3. At an even
budget all three children are identical, giving Aut(P_h) wreath C3.
At an odd budget their orders are G(h+1),G(h),G(h), with the first
strictly larger; no nonidentity rotation survives, giving the direct
product of their full groups. This proves the group order recurrence.
Burnside on the permitted rotations proves(11) and(12).

The group orders at budgets1 through7 are

    1, 3, 3, 81, 27, 81, 729,

on orders1,3,5,9,11,15,19. Full native isomorphism enumeration checks
these groups, and direct cycle counts check(11) for q=1,2,3,4. The
nine-vertex order81 member is the familiar C3[C3,C3,C3]. The
eleven-vertex order27 member is not the different twelve-vertex81/27
control. The latter still requires its own literal relations.

Here the binary budget decides where cyclic rotations survive at an
LCA. It is not Frobenius on a finite field. These are exact decorated
orbit recurrences for a Ramsey-extremal family, not the unrestricted
tournament-isomorphism sequence A000568.

## Precisely scoped next questions

1. Does every odd carousel T_(2r+1) admit(4)? The complete banks here
   stop at9. This asks for a uniform coloring rule, not more size samples.
2. Does |T|<=G(tau(T)) hold for all U5-free tournaments? Liu gives its
   sharp power-law relaxation. The present packet proves the integer
   statement only for the declared palette. A separate literature and
   proof audit is needed before calling the broader formulation new.
3. What replacement recurrence is optimal after adjoining Paley7?
   Its seven vertices at budget3 are a decisive obstruction to G.
4. Which actual path complements remain in a certified palette, and
   what full target relation must be retained when they leave? Scalar
   capacities alone cannot answer this question.
