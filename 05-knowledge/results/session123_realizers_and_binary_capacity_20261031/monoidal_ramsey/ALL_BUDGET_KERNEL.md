# Exact all-budget capacity of the full pentagon-composition class

2026-10-31. **Proved, scoped. No novelty-priority claim.**

Let D consist of nonempty finite graphs generated from a singleton by binary
join, disjoint union, and full C5 substitution (all five children nonempty).
Let F(A,B) be the maximum order in D with clique number at most A and
independence number at most B. This is not the class of all graphs and this
quantity is not Shannon capacity.

## Main theorem

For every positive integer A,B,

    F(A,B) = G(A,B) = sum_{0<=i<A, 0<=j<B} 2^popcount(i & j).

Thus the formerly eventual dyadic translation holds from the origin, and
balanced pentagon substitutions supply optimal witnesses for every budget.
The proof below is uniform, not an extrapolation of computed rows.

Put G(0,B)=G(A,0)=0. Splitting each index into its low bit and remaining bits
proves, for a,b>=0,

    G(2a,2b)       = 5 G(a,b),
    G(2a+1,2b)     = 3 G(a,b)+2 G(a+1,b),
    G(2a,2b+1)     = 3 G(a,b)+2 G(a,b+1),
    G(2a+1,2b+1)   = 2 G(a,b)+G(a+1,b)+G(a,b+1)+G(a+1,b+1).

Also G(1,B)=B, G(A,1)=A, G is coordinatewise increasing, and its rectangular
mixed difference at (x,y) is 2^popcount(x & y)>0.

## Weighted pentagon lemma

For nonnegative integer ports a_i,b_i, indexed modulo5, suppose

    a_i+a_(i+1)<=A,       b_i+b_(i+2)<=B.

Then sum_i G(a_i,b_i)<=G(A,B).

Proof by strong induction on A+B. For a zero parent budget every contribution
vanishes. If A=1, positive a_i equal1 and form an independent set of C5,
hence occupy at most two nonadjacent slots. Their b-values sum to at most B,
so the assertion follows from G(1,b)=b. The case B=1 is symmetric.
Assume A,B>=2 and write a_i=2x_i+e_i, b_i=2y_i+f_i, with bits e_i,f_i.
We form five layers of rounded half-ports. Each odd port is rounded upward
in exactly two layers; even ports always remain their half-value.

Use layer labels Z/5. The sets P_i={i,i+2} are disjoint for adjacent i;
the sets Q_i={i,i+1} are disjoint for distance-two i; and |P_i intersect Q_i|=1.
For an odd parent budget 2h+1, rounding every odd port upward on a common
set P of two layers gives layer budgets h+1 on P and h otherwise:
for two odd endpoints the floor sum is at most h-1; for one odd endpoint
it is at most h; for two even endpoints it is at most h.
For an even parent budget 2h, proper two-fold sets P_i give budget h in
every layer: two odd endpoints have floor sum at most h-1 and receive
at most one upward rounding, and one odd endpoint has floor sum at most h-1.
These observations apply to either constraint cycle.

* Both parents even: round a on P_i and b on Q_i. Each layer has budgets
  (A/2,B/2). Since their intersection has size1, the four binary recurrences
  say the sum of five contributions for each child is exactly G(a_i,b_i).
  Apply induction in each layer and the even/even recurrence.
* Both parents odd: choose common size-two sets P,Q with intersection1.
  Round a on P and b on Q. The five budget pairs comprise two low/low,
  one high/low, one low/high, one high/high. The child contribution is exact,
  and induction gives the odd/odd recurrence.
* A=2h+1 odd, B=2k even: round b on Q_i. Let S be the indices with e_i=f_i=1.
  If S is proper, the edges Q_i for i in S form a proper edge subgraph of
  a five-cycle. It has a vertex cover P of size at most2: it is a subgraph
  of the four-edge path, which has such a cover. Pad P to size2 and round
  every odd a on P. For i in S, r_i=|P intersect Q_i|>=1. Its layered total
  minus G(a_i,b_i) is

      (r_i-1)[G(x_i+1,y_i+1)-G(x_i+1,y_i)-G(x_i,y_i+1)+G(x_i,y_i)]>=0.

  All other children have exact layered totals. Thus the original sum is
  at most the layered sum, which induction bounds by
  3G(h,k)+2G(h+1,k)=G(A,B).
  If S is all five indices, use P_i,Q_i instead. All a_i are odd, so their
  adjacent floor sums are at most h-1, and proper rounding adds at most1.
  Every layer therefore has the **low** budgets (h,k). Child totals are
  exact. The result is at most5G(h,k)<=3G(h,k)+2G(h+1,k).
* The remaining parity is symmetric, interchanging the two constraint
  cycles (their isomorphism is multiplication of indices by2 modulo5).

Every layer parent sum is strictly smaller than A+B. This completes induction.

## Upper bound for actual graphs

Specialize the lemma to two adjacent nonzero slots, with a-ports x,y and
both b-ports B, obtaining G(x,B)+G(y,B)<=G(x+y,B). Two nonadjacent nonzero
slots similarly prove superadditivity in the second coordinate. Zeros here
are merely analytical ports, not permission to use empty substitution
children in D. Structural induction now handles joins, disjoint unions,
and full C5 substitutions using respectively these two superadditivities
and the weighted lemma. A singleton is bounded by G(A,B)>=1. Hence every
allowed graph has order at most G(omega,alpha), proving the upper bound.

## Actual optimal witnesses

Construct W(A,B) recursively. W(1,B) is the edgeless graph on B vertices;
W(A,1) is the clique on A vertices. For A,B>=2 set h=floor(A/2), k=floor(B/2).
In C5 slots0,...,4 choose child W(h+e_i,k+f_i), where e_i is1 exactly when
A is odd and i in {0,2}, and f_i is1 exactly when B is odd and i in {0,1}.
Clique-heavy slots are nonadjacent; independence-heavy slots are adjacent.
Their intersection is one slot. All children have positive budgets, and
all recursive parent sums decrease. Exact substitution rules give clique
number A, independence number B, and the binary recurrences give order
G(A,B). Base cliques/edgeless graphs are built by joins/unions of singletons.
This proves attainability **within D**, without assuming hereditary closure.

## Exact translation and fast evaluation

Let L=2^ceil(log2 A), with L(1)=1. Since i<A<=L, the summand is periodic in
j modulo L. Summing a full period gives

    G(A,B+L)-G(A,B) = sum_{i<A} L (3/2)^popcount(i) = L w(A),  B>=0.

The displayed digit sum for w agrees with the inherited slope recurrence.
This proves threshold1 for positive graph budgets (and a formal origin0).
Least period/sharpness results in Session122 remain compatible; no separate
minimal-period claim is needed here.

A four-state binary digit DP evaluates G in O(log max(A,B)) arithmetic
operations and constant many big-integer registers: scan bits of A-1,B-1,
keep whether each prefix is tight, and weight a chosen bit pair (u,v) by
2^(uv). This is a bit-complexity claim only after charging integer arithmetic;
orders can be exponentially large in input bit length. The accompanying
checker implements this DP and a memoized symbolic witness DAG. Expanding
an actual graph necessarily costs at least its output size.

## Source, target, preserved and lost data

Source: graph construction trees in D. Target: integer budget/order profiles.
Map: a graph goes to (omega,alpha,order); optimization gives F=G.
Preserved: exact extrema, explicit attaining constructors, budget composition
inequalities. Lost: graph isomorphism, adjacency, number and classification of
all optimizers. This is an exact numerical theorem, not an equivalence of
graph categories, not a realizability topos, and not an unrestricted Ramsey
breakthrough. The binary matrix [[1,1],[1,2]] and its tensor powers provide a
numerical factorization only; the parity-cover argument is the essential
additional graph-composition upper-bound proof.

## Prior art and verification boundaries

Session122 established exact rows through8, finite generation and eventual
translation; probes through16 and a million random port tuples suggested this
formula but cannot prove it. The proof above removes those finite restrictions.
An initial primary-source search found the directly relevant 1988 publisher
record **Ramsey numbers based on C5-decompositions**,
https://www.sciencedirect.com/science/article/pii/0012365X88900659 . Its full text
was unavailable (HTTP403). Thus prior-art equivalence remains unassessed;
no novelty claim is made. This citation is a search lead, not proof evidence.
Independent parity auditing resides in ../audit_ramsey/AUDIT.md.
