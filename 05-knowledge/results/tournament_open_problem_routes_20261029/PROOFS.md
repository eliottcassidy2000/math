# Two recovered gaps closed, and a repaired Erdős–Hajnal interface

2026-10-29. **PROVED from explicitly CITED classical inputs / elementary
PROVED deductions; same-agent audit. No literature-priority claim.**
Finite tests are guards, not proofs of the all-size or asymptotic statements.

## P1. The prime-prefix Legendre tournament is quasirandom

This closes the precise repository question
`HYP-9375-the-legendre-tournament-on-primes-three-mod-four-is-quasirandom.md`.
The October25 packet suggested a large-sieve route but did not execute it.
The October26 [Legendre sidecar](../root27_machine_bridges_20261026/anchor/LEGENDRE_SIDECAR.md)
correctly supplied the known CRT/Dirichlet extension theorem and warned that
abstract universality alone does not imply prime-prefix frequencies.

Put P(X)={p<=X: p prime, p=3 mod4}, N=|P(X)|, and orient p->q iff
(p/q)=1. Let K have entries K_pq=(p/q) off the diagonal and zero on it.
Quadratic reciprocity makes K skew-symmetric. These are prime vertices with
varying moduli, not elements of a single Paley field.

**CITED input:** Heath-Brown, *A mean value estimate for real character
sums*, Acta Arithmetica72 (1995),235–275, Theorem1, printed p.237:

    sum*_{m<=M} |sum*_{n<=L} a_n (n/m)|^2
       <<_epsilon (ML)^epsilon (M+L) sum*_{n<=L}|a_n|^2.       (1)

Stars restrict BOTH indices to positive odd squarefree integers. Coefficients
are arbitrary complex numbers. [Publisher](https://doi.org/10.4064/aa-72-3-235-275),
[author's Oxford manuscript](https://ora.ox.ac.uk/objects/uuid%3Ab188b365-d8d1-4267-8548-49f01e6cb6a7/files/m6af40476f1b860392271a2d8ecae593b).
The complete displayed statement was checked in indexed primary text; direct
PDF fetches returned403. The same inequality is independently quoted in
Soundararajan–Young, [JEMS12 (2010), §2.3](https://ems.press/content/serial-article-files/31734).
We use this published theorem, not a numerical estimate of its constant.

Set M=L=floor(X), support a only on P(X), and restrict the outer sum to
P(X). The resulting matrix is K^t. Restriction can only decrease a sum of
squares, and transposition preserves singular values. Thus for every epsilon>0,

    ||K||_2 <<_epsilon X^(1/2+epsilon).                       (2)

The [prime number theorem in a fixed progression](https://terrytao.wordpress.com/2009/09/24/the-prime-number-theorem-in-arithmetic-progressions-and-dueling-conspiracies/)
gives N~X/(2 log X). Consequently

    delta_X := ||K||_2/N <<_epsilon X^(-1/2+epsilon) log X ->0 (3)

on choosing epsilon<1/2. In particular every two subsets A,B of P(X), even
overlapping or chosen adversarially, satisfy

    |sum_{p in A,q in B}K_pq| <= ||K||_2 sqrt(|A||B|).         (4)

Here is a direct counting proof, avoiding an unstated quasirandomness
equivalence. Put W=(J+K)/2, including W_pp=1/2. For a fixed labelled
k-vertex tournament F with m=binom(k,2) arcs, expand

    E_{x_1,...,x_k in P(X)} product_{i->j in F}W(x_i,x_j).

Replace factors one at a time by1/2. In each difference fix the other k-2
variables. All remaining factors incident with x_i depend only on x_i;
those incident with x_j depend only on x_j; there is only one edge between
i and j. They are functions f,g in [0,1], with an additional constant in
[0,1]. The absolute difference is at most

    |f^t K g|/(2N^2) <= ||K||_2/(2N).

Summing the m replacements proves

    |t(F,W)-2^(-m)| <= m delta_X/2.                          (5)

Repeated sampled vertices have probability at most binom(k,2)/N. Removing
them and normalizing by the falling factorial (N)_k changes the answer by
O_k(1/N). Therefore the proportion of ordered DISTINCT k-tuples inducing
the specified labelled F is

    2^(-m) + O_{k,epsilon}(X^(-1/2+epsilon) log X).           (6)

For an unlabelled class [F], multiply by k!/|Aut(F)|:

    #{k-subsets inducing F}/binom(N,k)
      -> k!/(|Aut(F)| 2^binom(k,2)).                         (7)

This is exactly HYP-9375. It is unconditional; no GRH or independence of
prime Legendre symbols was assumed. Constants are not explicitly evaluated.
The same conclusion holds for the first N primes3 mod4 by taking X at the
last prime. The proof is a standard spectral/counting consequence of the
classical large sieve; the contribution here is closing the recorded gap.

**Robust extension.** For any fixed eta>0, arbitrary subsets S_X of P(X)
with |S_X|>=X^(1/2+eta) are quasirandom. Apply (2) to the restricted matrix
and choose epsilon<eta. This includes arithmetic or adversarial selections
provided that cardinality hypothesis actually holds. It supplies no claim
at or below the square-root scale.

**81/27 control.** Inherit the two actual twelve-vertex classes and their
automorphism orders from [DEFECT_TRANSFER](../overnight_rank_return_20261029/DEFECT_TRANSFER.md).
Equation(7) gives their positive limiting frequencies

    12!/(81*2^66),       12!/(27*2^66),

so the second occurs three times as often. The same labelled cyclic-triple
field does not justify merging these classes. This is a deduction from
their inherited automorphism certificates, not a fresh twelve-vertex census.

The inverse-automorphism law is also the stationary class law of the incoming
[fair refresh process](../session105_metagraph_recursion_20261029/README.md).
Both arise by pushing forward the uniform law on labelled tournaments.
This equality of measures does not identify their clocks or transition
operators; legal-Hamilton-path conditioning instead introduces an H factor.

**Scope and loss.** The source is the prime-indexed character matrix; the
map is its normalized kernel; the preserved predicates are averaged fixed
finite induced-pattern counts. Numerical height X is essential. The limit
forgets individual primes, exact incidences and rare named-port fibers.
It cannot replace the full relation in actual complement/re-marking work,
nor imply a fixed-root extension estimate, a finite-field Frobenius action,
or any A000568 recurrence speedup. All finite tournaments occur: this family
is not H-free and supplies no new general Erdős–Hajnal theorem.

**Möbius/Liouville connection with an actual operation.** The squarefree
restriction in(1) is the selector mu(n)^2, not the sign lambda(n). Away from
the prime3, the Jacobi columns for1 and9 are identical, whereas at3 their
values differ because of divisibility. Liouville gives both1 and9 sign+1
and does not remove this duplication. A proposed extension from prime
vertices to composite inputs must retain squarefree kernels, multiplicities
and coprimality exceptions. Dropping these guards is not justified by(1).
This is an arithmetic analogue of retaining actual port multiplicities,
not a theorem identifying the two carriers.

## P2. Close THM-1450's all-size out-neighbourhood VC=1 gap

Use the exact slug
`THM-1450-linearity-trivialises-vc-and-the-two-vc-one-tournaments.md`;
another theorem has the same numeric ID. Its part(C) was finite-exact
through seven and explicitly left one strong-tournament lemma open.

For a tournament T consider the set system {N^+(v):v in V(T)} on V(T),
using OPEN out-neighbourhoods, including v itself as an eligible witness
row. A set is shattered if all its subsets occur as row intersections.

**Theorem.** For every n>=3, this system has VC dimension at most one iff

    T is TT_n or C3 => TT_(n-3).                            (8)

Both have dimension exactly one. At n=1 the dimension is zero; at n=2
only TT_2 occurs and has dimension one. The orientation of the ordered
sum matters: the triangle is at the source end.

**Pair test.** For a->b, row a supplies {b} and row b supplies the empty
set. The pair is shattered iff there exist rows u,w with

    u->a, b->u;                 w->a, w->b.                 (9)

Equivalently, the arc a->b lies in a cyclic triangle and a vertex beats
both its endpoints. The two additional rows are necessarily distinct and
outside {a,b}.

**Missing lemma, proved.** Fix a directed triangle a->b->c->a and a vertex
w outside it. If w beats at least two triangle vertices, their arc together
with the third triangle vertex satisfies(9). If w beats exactly one, rename
it a. Then w->a->b->w is a triangle and c beats both w and a, again giving
(9). Thus in a VC<=1 tournament every outside vertex LOSES to every vertex
of every directed triangle.

If T has no triangle it is transitive. Otherwise take one triangle C.
It dominates the complement. That complement has no triangle D, since
the same lemma applied to D would force D to dominate C. It is therefore
transitive, proving necessity in(8). Conversely, in either family every
arc on a directed triangle is an arc of the unique top C3, and no vertex
beats both its endpoints. The pair test proves sufficiency. This also
proves every strong tournament of order>=4 has VC>=2. QED.

**Enumeration and transport.** For n>=3 there are exactly

    n! + 2*binom(n,3)*(n-3)! = 4n!/3                       (10)

labelled members and two unlabelled classes. Their automorphism groups
are respectively trivial and C3. Their ordinary Hamilton-path counts are
1 and3; this follows directly from the ordered-sum path decomposition.
Their largest transitive sizes are n and n-1, respectively. Equivalently,
VC<=1 is avoidance of the two four-vertex types source=>C3 and the strong
four-tournament: a shattered pair has a witness on four vertices, and these
are exactly the two disallowed types among the four four-vertex classes.
For n>=4 the converse of the second family has VC=2: a tail vertex now
dominates the triangle and shatters a pair. No pair within the transitive
tail or crossing from tail to triangle can be shattered: its arc lies in
no triangle. Every three-set meeting the tail contains such a pair. The
triangle itself has only the three singleton patterns and the full pattern,
so it too is not shattered. Do not silently identify out- and in-VC dimensions.

This class is not closed under arbitrary substitution: C3[C3,C3,C3] is
strong of order9, hence has VC>=2. Nor is it closed under converse for
n>=4. Low VC of one source does not close a metagraph carrier.

## P3. Repair the Hamilton-count obstruction, then keep the right observable

THM-1390's finite nonmonotonicity examples do NOT refute every possible
bound tr(T)>=phi(n,H(T)), and do not show that every EH route is impossible.
Indeed the inherited odd-cycle formula gives H(T)>=1+2c3(T), while
THM-1845's elementary triangle-deletion argument gives tr(T)>=n-c3(T).
Together they prove the perfectly valid, usually weak bound

    tr(T) >= n-(H(T)-1)/2.                                 (11)

No EH power bound follows: the right side becomes negative for large H.
What failed was a particular monotone proxy, not the existence of bounds.

For substitution T=Q[T_1,...,T_s] the correct exact extremal interface is

    tr(T)=max_{S: Q[S] transitive} sum_{i in S}tr(T_i).       (12)

A transitive set in T uses a transitive set of quotient vertices (choose
one representative in each nonempty intersection), and meets each child
in a transitive subset. Conversely their maximum subsets concatenate in
the quotient order. This proves both inequalities in(12).

For p>=1, a sufficient template certificate is

    max_{S transitive in Q}sum_{i in S}w_i >= (sum_i w_i^p)^(1/p)
       for every nonnegative weight vector w.              (13)

If every allowed template satisfies(13), induction yields
tr(T)>=|T|^(1/p) throughout its substitution closure, with singleton leaves.
This is a scoped construction-class theorem, not a theorem about every
H-free tournament. A missing prime-atom decomposition/certificate is still
needed to reach the latter. The classical EH substitution reduction and
the newer bounded-VC theorem are credited in SOURCES.md.

**Finite certificate formulation.** Normalize max_S sum_{i in S}w_i<=1.
The resulting rational polytope is compact (singletons give w_i<=1).
Since sum_i w_i^p is convex, checking it at every vertex suffices for(13).
Completeness of the vertex list must be certified; floating samples alone
do not suffice. This packet did not enumerate that polytope. The successor
[R1-R3 / Ramsey circuit interfaces](../ramsey_circuit_interfaces_20261030/RAMSEY.md)
classifies all12 vertices for regular T5, credits Liu's already-known
sharp exponent, and proves the stronger all-size integer capacity
F(k)=A006046(k) for the T5/T7/T9 hereditary substitution palette by
complete joint parity charts. Its scope is not arbitrary tournaments.

**Cheap hostile.** For the regular five-tournament Q, tr(Q)=3. The uniform
vector suggests p=log(5)/log(3), but an induced C3 with weights1,1,1 and
all other weights0 requires p>=log(3)/log(2), which is strictly larger.
Balanced block counts therefore do not certify unequal attachments.
This is the same reason to retain joint ports in the path-counting lane,
with a different observable and a different composition law.

These elementary deductions carry no claim of novelty in the EH literature.
