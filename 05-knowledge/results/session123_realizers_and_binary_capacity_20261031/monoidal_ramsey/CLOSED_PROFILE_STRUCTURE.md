# A genuine closed profile category—and what its internal hom forgets

2026-10-31. Proven formalization and explicit obstruction; no novelty claim.

## Typed closed structure

Let Q be the complete lattice of downward-closed subsets of the ordered
abelian group Z. Thus elements are bottom=empty, top=Z, and principal ideals
 down(n)={z:z<=n}. Tensor is Minkowski sum, with empty absorbing; its unit
is down(0). The scalar residual is

    D => E = {z in Z : z+D subset E}.

It is downward closed and satisfies H+D subset E iff H subset (D=>E).
This definition handles top/bottom without any undefined infinity subtraction.
Minkowski sum distributes over arbitrary unions.

Take profiles f:N->Q, ordered pointwise. Define

    (f tensor g)(n) = union_{i+j=n} (f(i)+g(j)),
    unit(0)=down(0), unit(n)=empty for n>0,
    [f,g](n) = intersection_{m>=0} (f(m) => g(n+m)).

Associativity and symmetry follow by reindexing finite decompositions of n;
the unit is immediate. Most importantly,

    h tensor f <= g  iff  h <= [f,g].

Indeed, the left condition says for every n,m that h(n)+f(m) is included
in g(n+m); apply the scalar adjunction and intersect over m. Hence this is
a genuine thin closed symmetric monoidal category, not just an analogy.
For finite integer-valued functions, read f(n) as down(f(n)); then convolution
is max-plus and [f,g](n) is the possibly negative-infinite infimum
of g(n+m)-f(m). No extended-infinity subtraction is being presumed.
This is the discrete convolution pattern of Brian Day's classical construction;
the specialized proof here is self-contained. Primary background:
https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/DayReport.pdf

## Exact relation to graphs, and an internal-hom obstruction

Represent the class-D capacity row by f_A(n)=down(F(A,n)), formally setting
F(A,0)=0. Principal ideals mean guaranteed lower order thresholds, not that
every smaller graph order is attainable: D was not assumed hereditary.
Disjoint union witnesses the appropriate profile tensor inequalities.
This numeric abstraction forgets graph witnesses, port incidence and optimizer
classification; it is not an equivalence of graph categories or a realizability
topos. Negative thresholds have only formal ordered-group meaning.

The ALL_BUDGET_KERNEL theorem proves the complete residual classification:

    [f_A,f_C] = f_C             if A<=C,
    [f_A,f_C] = bottom profile  if A>C.

For A<=C, superadditivity of F(C,.) and monotonicity in C give
F(C,n+m)>=F(C,n)+F(C,m)>=F(C,n)+F(A,m); m=0 attains equality in the infimum.
For A>C, the exact translation formula yields F(A,m)=w(A)m+O(1),
F(C,n+m)=w(C)m+O(1), with w(A)>w(C) because the digit-sum formula for w
has positive increments. The differences tend to minus infinity.
Thus internal hom between whole capacity rows is completely degenerate.
Merely placing these rows in a closed category does not extract additional
graph constraints: witness/arrival/atomic information must be retained.

## A hostile test: forgetting arrival lengths is not strongly monoidal

For a superadditive knapsack profile f with a sharp block of length5 and
profit50, set

    bar(f)(n)=sup_{t>=0} [f(n+5t)-50t].

Consider finite coin sets, where a coin is (length,profit):

    f: (1,0),(2,19),(5,50),
    g: (1,0),(4,39),(5,50).

For a combination define deficit=10*length-profit. The nonsharp coins have
costs10 and1. At residue1 modulo5, f's least deficit is3, realized by three
length2 coins (arrival length6); g's is4, realized by four length4 coins
(arrival length16). Hence bar(f)(1)=7 and bar(g)(1)=6, while both bars at0
are0. Consequently (bar(f) tensor bar(g))(1)=7.
In the mixed coin system, a length2 and a length4 coin have length6, profit58,
and deficit2. No cost0 or1 combination has residue1; any length1 coin costs10.
Thus bar(f tensor g)(1)=8>7. The displayed saturation is not strongly
monoidal for the unchanged N-indexed convolution. Borrowing a sharp block
can permit mixed resource allocations unavailable in nonnegative splits of1.
A localization with a modified tensor or integer-indexed resources may repair
this; the counterexample does not refute such constructions. Recording both
residue deficit and first feasible arrival avoids this specific loss.

## Canonical minimal generators (general lemma)

Suppose an integer-valued superadditive profile f, f(0)=0, is the knapsack
closure of its optimal coins at budgets1 through M. Call b irreducible if

    f(b)>max_{1<=i<b} [f(i)+f(b-i)],

with budget1 irreducible by convention. All irreducibles are<=M: a representation
of a budget>M uses at least two coins and can be split into two nonempty
subcombinations; superadditivity forces equality in the displayed maximum.
Conversely recursively splitting any reducible budget strictly decreases its
positive parts and terminates at irreducibles. Thus these are exactly the
unique minimal budget-indexed optimal coin generators; every generating set
of optimal coins must include each irreducible, as omitting one would supply
a forbidden split. This statement does not classify all graph optimizers.
For class D, irreducibility means every optimal graph at that budget is
connected. An optimal graph at budget b has actual independence number b:
otherwise adjoining an isolated singleton would increase its order while
preserving both budgets (the clique budget is positive). Every connected component of a D graph lies in D: full C5 substitutions
and nontrivial joins are connected, and unions inherit components recursively.
A disconnected optimal graph therefore supplies an equality split by grouping
components; conversely an equality split constructs a disconnected optimizer.
The implication does not say every optimizer at a reducible budget is disconnected.
