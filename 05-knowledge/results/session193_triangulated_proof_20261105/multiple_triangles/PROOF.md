# Root-cut path parity closes every extreme-arc sector

Session 193, 5 November 2026. **PROVED and independently reviewed.**
This packet proves an all-order lower bound and reduces the still-open rooted
deletion minimum to a substantially smaller structural sector. It does not
prove the unrestricted minimum or the Hamiltonian-path spectrum conjecture.

## 1. Inheritance, scope, and notation

Write H(T) for directed Hamiltonian vertex words and H(empty)=1. The inherited
marked odd-cycle formula is

    Delta(T,v) = H(T)-H(T-v)
               = sum_{Gamma: v covered} 2^|Gamma|
               = 2 sum_{odd directed C containing v} H(T-V(C)).

Collections Gamma consist of vertex-disjoint directed odd cycles. Each cycle
is counted modulo cyclic rotation, with its directed cyclic order retained.
The inputs are [Session 191 rooted growth](../../session191_tournament_growth_20261104/deletion_growth/PROOF.md) and [Session 192 even-path continuation](../../session192_surprising_approaches_20261105/weighted_cuts/EVEN_PATH_CONTINUATION.md). The new sharp single-triangle bound g(n) is [Session 193 single-triangle proof](../single_triangle/PROOF.md), Theorem 6, independently reviewed.
Thus

    g(4),...,g(8) = 2,4,6,10,16;
    g(9+3k)=26*5^k, g(10+3k)=46*5^k,
    g(11+3k)=78*5^k, k>=0.

For a root a of a tournament A, let R_a(A) count ordered spanning two-path
covers whose second path starts at a; the first path may be empty. Let e_a(A)
count Hamiltonian paths ending at a. Dually L_b(B) has its first path ending
at b and s_b(B) counts Hamiltonian paths starting at b.

For a tournament T rooted at v, put A=N^-(v), B=N^+(v). When T is strong of order at least three, these sets are nonempty. Reverse cross arcs are the arcs b->a with b in B, a in A;
their number t is exactly the number of triangles through v.

The live concepts are marked odd cycles, parity-resolved rooted paths,
positive cut capacities, SCC extreme blocks, and one-triangle replacement.
The canonical hostile is S191's H=33 tournament: one-at-a-time removal of
reverse cross arcs raises Delta. We do not assume orientation monotonicity.

## 2. Elementary rooted-path parity identity

For each nonempty directed simple path P starting at a, give P the weight
H(A-V(P)). Let O_a and V_a be the total weights when |V(P)| is odd and even,
respectively. Then

    O_a+V_a=R_a(A),     O_a-V_a=e_a(A).                 (1)

The first identity is just the definition of the ordered path cover.
Here is a complete induction for the second. Define the signed sum F_a(A)
with sign (-1)^(|V(P)|-1). Separate P=(a) from paths whose next vertex is b:

    F_a(A)=H(A-a)-sum_{b:a->b} F_b(A-a).

For |A|=1 it is 1=e_a(A). Induction and the last-vertex recurrence give,
when |A|>=2,

    F_a(A)=H(A-a)-sum_{b:a->b}e_b(A-a)
          =sum_{b:b->a} e_b(A-a)=e_a(A).

This proves (1) for all tournaments, including weak ones. By reversing arcs,
for paths ending at b the odd/even sums are (L_b+s_b)/2 and (L_b-s_b)/2.
In particular these are nonnegative integers. No asymptotic or formal signed
weight is used in the lower bound below.

## 3. A universal single-crossing lower bound

**Theorem 1.** For every rooted tournament (T,v), with the empty-side sum
understood as zero,

    Delta(T,v) >=
       sum_{b->a, b in B,a in A}
           [ R_a(A)L_b(B)+e_a(A)s_b(B) ].              (2)

Each bracket is a positive even integer, at least two.

**Proof.** Retain only the following marked odd-cycle collections:

* the unique cycle through v has the form

      v, P_B, P_A, v,

  where P_B is a nonempty directed path entirely in B ending at b, P_A is a
  nonempty directed path entirely in A starting at a, and b->a is a reverse
  cross arc;
* every other cycle is wholly inside A-V(P_A) or B-V(P_B).

All displayed edges at v are valid by the definitions of A and B. The root
cycle is odd exactly when |V(P_A)| and |V(P_B)| have the same parity. Its
weight is two. The remaining cycles contribute independently the exact
factor H(A-V(P_A)) H(B-V(P_B)), by the inherited odd-cycle formula.

This is an injection into the full marked collection set. The root cycle
and its direction determine P_A, P_B and the sole B->A transition uniquely.
Its other cycles retain their side, so no collection is counted twice,
including collections assigned to different reverse cross arcs.

For a fixed b->a the retained weight, using (1), is

    2[ ((R_a+e_a)/2)((L_b+s_b)/2)
       +((R_a-e_a)/2)((L_b-s_b)/2) ]
      =R_a L_b+e_a s_b.

Summing and discarding the other positive objects proves (2). A fixed edge
has the triangle v,b,a with arbitrary internal-side cycle collections, so
its bracket is at least 2H(A-a)H(B-b)>=2. QED.

The omitted quantity is itself a nonnegative even integer. Its objects
are exactly those whose root cycle crosses B->A at least twice, or whose
root cycle crosses once but some other odd cycle meets both A and B.
Thus (2) is a structural projection with an explicit, positive remainder,
not an unjustified independent sum of cut capacities.

**Corollary 2 (exact star formula).** If all reverse cross arcs share one
head a, then (2) is equality and becomes

    Delta = R_a(A) sum_{b:b->a}L_b(B)
               +e_a(A) sum_{b:b->a}s_b(B).             (3)

The dual statement holds if they share one tail b.

Indeed every directed cycle crossing the two sides must use that common
head (or tail). The root cycle cannot cross twice, and no disjoint additional
cycle can cross at all. Hence there are no omitted objects. This includes,
but is strictly broader than, the inherited single-cross-arc formula.

## 4. Extreme arcs give an all-order sector closure and strict gap

Let A_0 be the source SCC of A, and B_1 the sink SCC of B. Call a reverse
cross arc **extreme** if it goes from B_1 to A_0. Let k be their number.

**Theorem 3.** If T is strong, has n>=4 vertices and k>=1, then

    Delta(T,v) >= k*g(n)+2(t-k).                       (4)

In particular Delta>=g(n); if t>=2 then Delta>=g(n)+2. A root attaining
g(n) in this sector must lie in exactly one triangle.

**Proof.** Fix an extreme arc b->a. A vertex in the source SCC of A reaches
all of A, and every B vertex reaches b in its sink SCC. Construct T_e on the
same vertices with the same internal sides and arcs at v, keeping b->a and
orienting every other cross arc A->B. Those two reachability properties make
T_e strong. It has exactly one root triangle. The exact one-triangle identity
(or Corollary 2) gives

    Delta(T_e,v)=R_a(A)L_b(B)+e_a(A)s_b(B)>=g(n).

This compares one summand, not a sequence of arc changes; there is no claim
that intermediate tournaments have decreasing Delta. Apply this bound to
each of the k extreme summands in (2), and the lower bound two to each of
the t-k others. QED.

**Corollary 4.** The same closure applies whenever at least one of A and B
is strong (singletons included), or whenever the reverse cross graph is a
star. A possible counterexample to Delta>=g(n), or a multi-triangle equality
case, must have both sides weak and no arc B_1->A_0.

To see the first assertion, strongness forces some reverse arc into A_0 and
some reverse arc out of B_1. If A is strong, the latter is extreme; if B is
strong, the former is extreme. For a common-head star the head must belong
to A_0 and some tail must belong to B_1. The common-tail case is dual.

More explicitly, when A and B are nonempty, strongness of T is equivalent to the two conditions

    some B vertex points to A_0;
    some B_1 vertex points to A.

They suffice because v reaches all B and then all A through A_0, while all
B reaches v through B_1 and A; every A already points to v. Necessity follows
from the source and sink SCC obstructions. The two required arcs need not
be the same; the remaining sector is precisely that separation.

## 5. Hostiles and what remains open

S191's greedy-smoothing hostile has rows (48,5,9,3,14,30), root 0.
Its A={1,2,3} is strong C3; B={4,5} is transitive, 5->4. There are six
root triangles and three extreme arcs from vertex 4. Theorem 1's lower
bound is exactly 30, equal to its actual Delta=33-3. Reversing an individual
reverse arc raises Delta to 32 or 34, as inherited. This is consistent with
(2): that theorem retains combinatorial objects from the original T, not
an orientation-monotone path through intermediate tournaments.

The pure single-crossing bound need not itself reach g(n) in the residual
sector. Let A=(a_0->...->a_r), B=(b_0->...->b_s), r,s>=1, with A->v->B.
Set all cross arcs A->B except b_0->a_0 and b_s->a_r. The result is strong,
both sides are weak, and the extreme pair is a_0->b_s. The right side of
(2) is exactly 2^r+2^s: on a transitive side the first-root R and
last-root L are the appropriate powers of two, the opposite responses
are one, and the relevant endpoint correction products vanish. For r=s=2, this is 8<g(7)=10, while the actual
Delta is 24. For r=s=3 it is 16<g(9)=26, while actual Delta is 80.
Thus the positive remainder cannot be dropped in the remaining proof.
These are obstructions to sufficiency of this lower-bound statistic, not
counterexamples to the desired global inequality.

There is a guaranteed first remainder term. For every reverse arc b->a
with a in A_0 and every reverse arc d->c with d in B_1, absence of extreme
arcs forces b outside B_1 and c outside A_0. Then

    v -> b -> a -> d -> c -> v

is a directed 5-cycle, with five distinct vertices. Different arc pairs give
different directed cycles. Therefore the remainder in (2) is at least

    2 sum_{b->a, a in A_0; d->c,d in B_1}
               H(T-{v,b,a,d,c}).                       (5)

This extra charge is proved, but is not claimed sufficient at all orders.
The unresolved task is to control longer repeated-crossing root cycles and
crossing disjoint cycles. Neither the lower bound nor sector closure implies
an interval-filling or full Hamiltonian-path spectrum theorem.


## 6. Exact controls

The standalone verifier checks all 1,099 labelled tournaments through order
five, all 5,405 rooted presentations, 55,043 literal rooted simple paths,
and 10,758 literal single-crossing root cycles. It also independently computes
the inequality and sector conclusions for all 50,779 rooted strong
presentations of orders four through eight in the inherited complete
7,412-class bank. It checks the H33 hostile and all six bad greedy steps,
star equality, and the forced five-cycle remainder. Regenerate the bank with `g++ -O3 -std=c++17 generate_bank.cpp -o generate_bank` and `./generate_bank > inherited_tournaments_through8.txt`; then run `python3 verify.py` and `python3 -O verify.py`. These finite controls are
not substitutes for the all-order arguments above.
