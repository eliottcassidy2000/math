# Obstruction lifetimes, endpoint conservation, and unavoidable catalogue growth

**PROVED structural framework and counting bounds, conditional only on the
explicitly inherited OCF/strong-deletion facts; FINITE-EXACT checks through
order seven. The full odd H-spectrum remains OPEN. No external priority claim.**

## 1. Inheritance and the question this answers

For a tournament T, H(T) counts ordered directed Hamiltonian vertex words;
H(empty)=1. Let C(T) delete all singleton strongly connected components
(SCCs). A nonempty tournament with C(T)=T is called a core. Isomorphism is
ordinary directed isomorphism: do not merge converse pairs.

Inherited at repository commit d11ba3658601:

* S133 `spectrum/OBSTRUCTION_PROOF.md`: the minimal forbidden basis F_B of
  H<=B is finite, has order at most 3r where r=min{j:3^j>B}, and its equality
  at two thresholds certifies an empty spectral interval. The SCC envelope
  formula below was already used in its `forbidden_basis.py` implementation.
* S133 `core_reduction/THEORY.md`: H(T[S])=H(T) iff S contains C(T).
  Consequently every single deletion from a nonempty core strictly lowers H.
  It also gives the exact K(B) core-order budget. None of these are new here.
* S134 `descent/THEORY.md`: each fixed positive induced gap has bounded
  endpoints; this rules out bounded-increment nested universality. It does
  not bound the number of different incomparable core witnesses.

The new quantitative question is: **how much explicit structural information
must the finite forbidden-basis hierarchy contain as its threshold grows?**
The answer below separates finite recognizability from uniformly small lists.

Concept board: (i) exact-H core fibers, (ii) threshold forbidden bases,
(iii) deletion-envelope lifetimes, (iv) SCC words, (v) catalogue complexity.

## 2. The exact lifetime transform

For a nonempty core F define

    h(F)=H(F),    d(F)=max_{v in V(F)} H(F-v).

The inherited strict deletion theorem gives 1<=d(F)<h(F), with both integers
odd. A tournament F is in the minimal forbidden basis F_B exactly when it
is a nonempty core and

    d(F) <= B < h(F).                                      (1)

Indeed minimality is equivalent to all one-vertex deletions having H<=B,
by induced monotonicity. Any singleton SCC could be deleted without lowering
H, so a minimal bad tournament has none. Conversely (1) verifies every
proper induced subtournament by placing it inside a one-vertex deletion.
Thus each core contributes one half-open integer interval [d,h). This is an
interval encoding of forbidden-basis membership, not a claim about persistent
homology. It retains the core class, not merely the multiset of intervals.

Let

    b(B)=|F_B|,
    a(m)=#{nonempty core classes F : h(F)=m},
    beta(m)=#{nonempty core classes F : d(F)=m}.

All are finite: a(m) uses K(m), and beta(m) is a subset of F_m by (1).
For every integer m>=2, counting intervals entering and leaving gives

    b(m)-b(m-1) = beta(m)-a(m).                            (2)

This endpoint conservation law identifies the exact core-fiber multiplicity
from the basis change **and its births**, not from b alone. In particular,

    m is attained  iff  a(m)>0,
    m is missing   iff  beta(m)=b(m)-b(m-1),              (3)

for m>=2. A missing m also forces beta(m)=0 because d(F) is itself an attained
count. These are exact encodings, not a uniform method for proving a(m)>0.
No scalar-cardinality plateau criterion is claimed; (2) retains births.
The special value 1 is realized by the empty/transitive core convention and
must not be inserted into the nonempty-core count a(1).

### Equivalent locally finite descent graph

There is also an exact graph encoding: use attained H values as vertices and,
for every nonempty core class F, put one class-labelled directed edge

    d(F) -> h(F).

Parallel edges remain distinct. This is an acyclic, locally finite directed
graph: indegree at m is a(m), outdegree is beta(m), and (1) identifies b(B)
with the number of edges crossing the cut {1,...,B}. Every attained h>1 is
reachable from 1. To see this, choose a core of count h, choose a deletion
attaining d, replace that deletion by its H-preserving core, and repeat;
positive counts strictly decrease until 1. For this particular recursive
choice the witnesses really are nested induced cores.

However an arbitrary numerical path through the graph need not admit one
compatible sequence of tournament witnesses: two incident edges may use
unrelated representatives of their shared scalar H. This is precisely the
sidecar warning from S134, not an embedding theorem. The graph packages the
full spectrum and its finite threshold edge cuts; it does not manufacture a
vertex at a prescribed odd integer.

### SCC composition and the exact information sidecar

If F=C_1 ordinal-sum ... ordinal-sum C_s is its unique SCC word, set
h_i=H(C_i), d_i=max_v H(C_i-v), and P=product_i h_i. Then

    h(F)=P,    d(F)=max_i (P/h_i)d_i,
    d(F)/h(F)=max_i d_i/h_i.                              (4)

This follows from SCC multiplicativity after deleting a vertex in one block.
The pair (h,d) therefore suffices to propagate lifetimes along SCC words.
The actual ordered SCC identities must still be retained for isomorphism
counts. Equal scalar pairs do not identify tournaments or their embeddings.
The algebraic rule (4) was already used in S133; its role here is to make the
new interval transform and catalogue lower bounds compositional.

**Sharp loss witness at fixed order and strong status.** The strong six-vertex
rows (48,1,3,7,14,30) and (34,8,3,5,15,30) both have H=23. Their full deletion
decks are respectively

    (1,11,11,11,9,9),    (3,9,9,5,13,5).

Their lifetimes are [11,23) and [13,23). The first is a minimal obstruction
at threshold 11; the second is not. Thus even (order,H,strong) loses the
necessary deletion-envelope coordinate. Source padding gives the other
boundary: it preserves H but makes d=h, so the lifetime becomes empty.

Connection contract: tournament cores -> class-labelled intervals; preserved
predicate is membership in every threshold basis; lost data from bare
intervals includes adjacency, embeddings and the deletion deck; restore
actual core/SCC identities when those targets matter. The H23 pair is the
cheapest decisive test of scalar H compression.

## 3. Superpolynomial spikes in every explicit full basis hierarchy

**Theorem.** For every n>=4 there is an integer B_n satisfying

    f(n-1) <= B_n <= (n-1)!/2,
    b(B_n) >= 2^{n(n-1)/2-n+1}/[n!(n-1)!],               (5)

where f is the inherited strong minimum (Busch). Moreover the witnesses
counted in this lower bound are all strong and have exactly n vertices.

**Proof.** Fix the labelled directed Hamilton cycle
0->1->...->n-1->0. Its n arcs are distinct for n>=3. Every orientation of
all remaining binom(n,2)-n arcs is strong, yielding 2^{binom(n,2)-n}
labelled strong tournaments. Each ordinary isomorphism class has at most n!
labelled copies, hence this universe contains at least
2^{binom(n,2)-n}/n! distinct strong classes.

For each such class assign its integer birth d(F). Since every (n-1)-vertex
Hamilton word and its reversal cannot both be directed paths,

    d(F) <= (n-1)!/2.

Strict deletion gives d(F)<H(F), so F lies in F_{d(F)}. Pigeonholing the
classes into at most (n-1)!/2 possible birth values proves the lower bound
in (5). For n>=4, inherited Moon strong deletion supplies a vertex v with
F-v strong, so d(F)>=H(F-v)>=f(n-1). This proves the lower endpoint and
ensures B_n tends to infinity. QED.

With natural logarithms, (5) implies along these thresholds

    log b(B_n) >= (log(2)/2-o(1))
                    (log B_n / log log B_n)^2.           (6)

For clarity, log of the counting lower bound in (5) is
(log(2)/2)n^2-O(n log n), while log B_n<=n log n and B_n tends to infinity.
Monotonicity of x/log x for large x gives
(log B_n)/(log log B_n) <= (1+o(1))n, which yields (6).
In particular b(B) is not O(B^c) for any fixed c.

The inherited size cutoff provides the complementary universal upper bound

    b(B) <= sum_{j=1}^{3r} 2^{j(j-1)/2}
          = exp(O((log B)^2)),    r=min{j:3^j>B}.         (7)

**Scope.** This is a lower bound on the number of entries in an explicit
complete minimal-forbidden-tournament catalogue. It is not a lower bound on
compressed certificates, arithmetic descriptions, membership algorithms,
or proofs that a single proposed H is missing. It neither proves nor refutes
universality of the positive odd spectrum except 7 and 21. Counting many
witnesses at one H does not create additional H values.

## 4. A uniform constructive lower bound at every threshold

There are two nonisomorphic strong five-vertex tournaments with (h,d)=(15,5):

    rows (16,9,3,5,14),   deletion deck (3,5,5,5,3);
    rows (24,17,3,6,12),  deletion deck (5,5,5,5,5).

Their degree multisets (1,2,2,2,3) and (2,2,2,2,2) already distinguish them.
The counts are verified by direct permutation enumeration in `normal.json`;
all later steps are symbolic.

**Theorem.** For every integer B>=5,

    b(B) >= 2^{floor(log_15(3B))}.                         (8)

The logarithmic floor is the largest r with 15^r<=3B; no floating-point
arithmetic is needed.

**Proof.** Put r=floor(log_15(3B))>=1. Any ordered word of r copies chosen from
the two five-vertex atoms has H=15^r and d/H=1/3 by (4). There are exactly
2^r distinct ordinary core classes because the ordered SCC decomposition is
unique. Append k triangles, where k>=0 is the smallest integer for which
P=15^r3^k>B. The starting inequality 15^r<=3B and minimality of k imply
B<P<=3B. Triangles also have d/H=1/3, so every resulting core has d=P/3<=B<H.
All 2^r words therefore belong to F_B, and their different five-vertex SCC
words remain distinguishable after the common triangle suffix. QED.

This is an all-threshold lower bound, not only a subsequence count. Its
polynomial exponent log_15 2 is much weaker than the spikes (6), but it has
actual uniformly constructible witnesses.

## 5. Exact tests and hostile controls

`check.py` has no external packages or imported catalogue. It enumerates all
strong tournaments through order seven by fixing one Hamilton cycle and
orienting every remaining arc; Camion/Moon ensures every strong class is
represented. Canonicalization checks every permutation within equal-degree
cells. The numbers of ordinary strong classes are

    1,1,6,35,353       at orders 3,4,5,6,7.

Orbit masses independently match the exact first-SCC recurrence for labelled
strong tournaments. Literal permutation counts independently verify H and
every deletion for each representative. The 399 nonempty singleton-free
cores of total order <=7 are obtained from all eligible SCC words. Equation
(2) is checked over their **restricted finite universe**, thresholds 1..189;
this finite receipt does not claim to compute the complete bases that far.

The exact strong order-seven birth histogram is

    d:    15 17 19 23 25 27 29 31 33  37 41 43 45
    count: 1  7  5 13 10  5 27 14 41 104 19 68 39.

Thus b(37)>=104 already from strong seven-vertex classes alone. Their actual
adjacencies and deletion decks are retained in the receipt; larger-order
obstructions are not excluded by this lower bound.

Controls: same-(n,H,strong) H23 collision; source-padding nonminimality;
exact SCC max-ratio multiplication; endpoint-flow closure; robust H15 atom
nonisomorphism; word injectivity; integer triangle-scaling at every B from
5 through 100000. Normal and Python -O output must agree. The arbitrary-n
counting bounds (5) are output as exact integers/rationals for n=4..20;
small lower bounds below one are not falsely treated as strong numerics.

Reproduce:

    python3 check.py > normal.new.json
    python3 -O check.py > optimized.new.json
    cmp normal.new.json optimized.new.json
    cmp normal.json normal.new.json

Remaining frontier: give a substantially smaller compressed hierarchy,
characterize which (h,d) pairs are actually realizable, or use genuine
structural atoms to settle the open full spectrum. Neither interval notation
nor large forbidden bases discharges any of those realization obligations.
