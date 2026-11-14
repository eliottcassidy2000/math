# Exact factor selection with prime quotient tournaments

This is a constructive alternative to the
[singleton cyclic factor pair](../chanpak_absolute_reflection_20261114/README.md).
The module quotient supplies a small word grammar; every factor
below comes from an actual assignment of labelled vertices to runs.
It gives exact witnesses including 63, 97 and 881. It does not prove
that some quotient or factor family succeeds at every odd target.

## Count runs before choosing their physical sizes

Let Q be a tournament and let its vertex i be substituted by a
nonempty tournament A_i. Orient each cross edge as Q prescribes.
Compress each maximal consecutive run in the same block in a full
Hamiltonian path. The resulting word uses every label, adjacent
labels follow Q, and j_i occurrences of i can be read as an ordered
cover of A_i by j_i nonempty paths. Conversely every such word and
choice of covers uniquely expands back to a full path. Thus, with
G_j counting ordered nonempty covers,

    H = Σ_valid run multiplicities c_Q(j) · Π_i G_(j_i)(A_i).

There is no cancellation and the run slots are distinguished by
their positions. For TT_n, G_j=j! S(n,j), the number of surjections
of its labelled vertices into j named slots. In particular
G₁=1, G₂=2^n−2 and G₃=3^n−3·2^n+3.

If Q is strong, a directed Hamiltonian cycle of Q read from block
i and back to i gives one word visiting i twice and the others
once. This word alone contributes at least 2^|A_i|−2, since each
nonempty two-slot assignment supplies paths. At least one different
word visits every block once, so |A_i| ≤ floor(log₂(N+1)) at
prescribed H=N. This makes any fixed quotient bank a finite exact
decoder, also when the modules themselves are arbitrary. For
transitive modules use the surjection formulas instead of
enumerating any internal orientations.

## Two distinct product equations on a five-vertex quotient

Use the strong quotient with outgoing masks Q=(28,25,10,16,4),
vertices indexed from zero. It is prime in the modular-decomposition
sense: it has no proper nontrivial set to which every outside vertex
has a uniform orientation. The replay checks this directly.

First put arbitrary blocks A,B,C at positions 0,3,4; positions 1,2
are singletons. Put h=G₁ and p=G₂ for each variable block, so
p=f−2h with f the ordinary two-slot count. The complete word
coefficients, listed by occurrences of (0,3,4), are

| Occurrences | 111 | 112 | 121 | 122 | 211 | 212 | 221 | 222 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Coefficient | 9 | 5 | 3 | 3 | 3 | 2 | 1 | 1 |

These are finite exact coefficients even for unbounded input
orders: the three variable labels induce a transitive tournament,
so each can appear at most once between consecutive singleton
markers, and at most three times overall. Checking those finite
words gives the displayed table and no coefficient at three visits.
It factors to

    N + h_A h_B (f_C−2h_C)
      = (f_A+h_A) [ f_B(f_C−h_C) + h_B h_C ].

For transitive blocks of widths a,b,c≥1 this is the particularly
simple prescribed-target test

    N + 2^c−2 = (2^a+1) [ 2^b(2^c−1)+1 ].

Try the bounded c and every divisor of the form 1+2^a;
for its cofactor v require (v−1)/(2^c−1) to be an integral
power of two at least two. The converse is exact. One can
enumerate the exponents a instead of carrying out a general
factorization. For **N=63**, c=2 gives 65=5·13 and forces
b=2 for the factor choice 5,13. Taking a=b=c=2 gives an
eight-vertex strong witness with outgoing masks
(250,248,243,52,224,192,136,8).

On the **same quotient**, put arbitrary blocks A,B at positions
1,4 instead, leaving the other positions singletons. The possible
run multiplicities of (1,4) have coefficients 9,5,5,4 at
11,12,21,22 respectively. The same finite-word proof applies,
with four chunks around the three singleton markers. This gives

    4N − 11h_Ah_B = (4f_A−3h_A)(4f_B−3h_B).

With transitive widths a,b this is an ordinary Fermat-type factor
test for 4N−11: its factors must be 2^(a+2)−3 and 2^(b+2)−3.
For the prime **N=97**, the exact adjusted product is
377=13·29, hence a=2,b=3. The eight-vertex literal witness
has outgoing masks (248,245,241,22,224,200,136,8).

## Why 63 needs a larger strong quotient in this hierarchy

The [cyclic decoder](CYCLIC_RUNS.md) rules out 63 with *any*
three-block cyclic join, allowing arbitrary orders and arbitrary
tournament modules. Allowing a strong quotient of order four
does not help: every strong order-four tournament is a transitive
order with only its first-last arc reversed. Indeed a directed
Hamiltonian cycle leaves two chord orientations, and its cyclic
relabelings transitively permute the four chord patterns. The two
interior vertices of the transitive order form a module. They
combine into one ordinal block even after arbitrary substitution,
so the blow-up is still a three-block cyclic join.

For the singleton-branch exclusion at 63, first suppose its
other two input orders are at least two.
M=f_A f_B is a multiple of 16 in (63,94.5), hence 64 or 80.
The shifted factor product 189−2M is prime (61 or 29), forcing
one TT₂ (shifted factor 1), but the other f would be 16 or 20
and cannot even exceed that prime. With a singleton input the
test is f+h=63. The size bound gives order at most five; through
four the maximum is 37. At order five, if t counts directed
triangles and c counts spanning directed cycles, odd cycle
packing gives f+h=33+18t+6c. Equality would require 3t+c=5,
so a spanning cycle but at most one triangle. That is impossible:
a strong five-vertex tournament has outdegrees in {1,2,3}, and
the number of 1s equals the number of 3s, at most two. Thus
its triangle count 10−Σ_i binom(outdeg(i),2) is at least three.

Thus among strong quotients on at least three vertices, the minimum
quotient order at 63 is **exactly five**. This does not claim that
seven and twenty-one fail only because of small quotients; their
exclusion is global. It separates 63 from them using an actual
different constructive branch. A nontrivial scalar factorization
does not rescue 63 either, since every nontrivial factorization
uses 7 or 21 as one of its factors.

## A different exact source line reaching 881

Take the strong prime quotient (60,9,42,16,6,26). Leave all
positions singleton except position 1, fixed at TT₂, and position
3, an arbitrary tournament B. If positions 1 and 3 are allowed
unboundedly many visits, their respective multiplicities have
coefficients 31,14,33,21,8,8 at 11,12,21,22,31,32 and zero
elsewhere. Indeed those two labels induce a transitive pair,
so each occurs at most five times among the four singleton
markers; the replay enumerates all these words. The terms
at three visits to TT₂ vanish, so the actual reader is

    N = 97 G₁(B) + 56 G₂(B) = 56 f_B−15 h_B.

At B=TT_n this gives the exact ray 56·2^n−15: at n=1 it
returns 97, and at n=4 it reaches **881**, with outgoing masks
(1016,245,241,758,480,448,384,256,14,502). The method changes
both the geometry and the absolute coefficient, not just the phase.

Run **node replay_prime_quotient.js 63** to get actual widths
on the first five-vertex chart, or test any other target.
The independent literal replay checks all three named examples,
the finite multiplicity tables, and 27 three-width controls.
A rejection on these charts is not an exclusion from the scalar
spectrum. With a fixed quotient and transitive modules the
logarithmic width bound leaves only polylogarithmically many
outputs through X; universal coverage needs a genuinely
growing family or freely varying nontransitive modules.

### A depth requirement for bounded quotient arity

More generally, take a fixed finite bank of strong quotients
of maximum arity q, composed recursively with transitive leaves.
For any nontrivial output H≤X, every leaf width is O(log X):
its immediate strong parent is already at least 2^width−1,
and a single full quotient walk through one run in each parent
block shows counts do not decrease up the tree. At depth d there
are at most q^d leaves and exp(O(q^d)) bank-labelled tree shapes.
Hence such trees give at most exp(O(q^d log log X)) counts
through X. In particular allowing depths up to
(1−ε) log_q log X for fixed ε>0 still gives only X^o(1)
distinct counts. The bound does not apply to unbounded quotient
arity or to arbitrary, nontransitive leaves.
