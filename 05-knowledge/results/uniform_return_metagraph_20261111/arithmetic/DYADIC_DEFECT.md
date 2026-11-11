# Every fixed branching defect has finitely many dyadic count families

**PROVED elementary classification, independently reviewed.** The root and
arithmetic lanes derived the prefix-tail strengthening independently; the
metagraph lane refereed its mechanism and the intercept convention. This
strengthens [BOUNDED_TEMPLATES](../defect/BOUNDED_TEMPLATES.md). It does not
classify the d=1 stratum, nor prove source-return closure or Chan–Pak.

## 1. Canonical prefix and a bounded tail

Let(T,r) be strong. Set I=N^-(r), d=|I|>=2, R=T-({r} union I), and
b=F_r(2)-H(T). The inherited positive formula and exponential bound give

    b=(q-s)+delta(T-r),
    delta(U)=sum_{a<c in I} c_ac(U),
    b >= (2^d-2) H(R).                                    (1)

Each tagged odd-cycle atom in a capacity has positive weight at least2.
In particular every c_ac, and hence delta, is even. Write
K=floor(b/(2^d-2)). The nontrivial strong SCCs of R have a finite ordered
list: there are at most floor(log_3 K), each has order at most K. Their
total order is therefore at most K floor(log_3 K), using the previously
proved elementary kernel bound. The sharper known bounds are unnecessary
here. Also d<=floor(log_2(b+2)).

Call a residual vertex ordinary if it beats every inport. Define P to be
the maximal INITIAL sequence of SCCs of R all of whose vertices are
ordinary. Thus P is an actual induced ordinal prefix; it may be empty.
Let W be the induced UNROOTED tail on I together with R-P. The marked
root r is NOT part of W. Then P dominates W whenever P is nonempty.

The older exceptional-singleton atom argument bounds the nonordinary
singleton SCCs. For d>=3 there are at most floor(b/2). For d=2, if the
three nonordinary letters occur m0,m10,m01 times, then

    m0+binom(m10,2)+binom(m01,2)<=floor(b/2).

Since m<=binom(m,2)+1 for every nonnegative integer m, the number of
nonordinary singleton SCCs is at most floor(b/2)+2.

**Late-ordinary lemma.** There are at most floor(b/2) ordinary singleton
SCCs in R-P. If R-P is nonempty, its first SCC is nonordinary, so choose
a vertex y there and a in I with a->y. Every ordinary singleton x after
that SCC satisfies y->x, and x->c for any c in I other than a. Hence

    a -> y -> x -> c

closes by a tagged two-edge ear to an odd five-cycle. Distinct x give
distinct tagged atoms, even if the same y,a,c are used. Each costs at
least2 in delta. An ordinary singleton cannot itself be the first
nonordinary SCC, so these are all ordinary singleton SCCs of R-P. QED.

Consequently the TOTAL number of singleton SCCs in R-P is at most
b+2 for d=2, and at most b for d>=3. Together with the bounded strong
core list and bounded d, this bounds the entire tail W in terms of b.
The prefix P contains only a bounded list of nontrivial strong cores
separated by arbitrary transitive gaps. This is a canonical decomposition,
not a choice of a favorable representative in the same scalar H fibre.

## 2. Nonempty prefix is exactly neutrality

The inherited q-s identity counts the Hamilton paths of U=T-r starting
in I. Every Hamilton path starts in the first SCC of U, and every vertex
of a strong tournament can start a Hamilton path. Thus q=s exactly when
the first SCC of U avoids I. If it avoids I, it must be the first SCC of
R and must dominate I, so P is nonempty. Conversely the first SCC of a
nonempty P dominates all other vertices of U and avoids I. Therefore

    q=s  iff  P is nonempty.                              (2)

In particular every nonneutral source(q!=s) has P empty and hence bounded
order at fixed b. This holds for even as well as odd b. The unbounded
part of any fixed branching-defect stratum is entirely neutral.

For nonempty P, let h=H(P) and C=delta(W), using the same marked inport
set I inside W. In an ordered two-path cover counted by W_ac(U), the
second path starts at c in W and cannot enter P. Thus all of P lies in
the first path, as one initial Hamilton segment, and deletion of that
segment is a bijection onto a cover counted by W_ac(W). Consequently

    b=h C,     h divides b.                               (3)

The capacity bound on I gives C>=2^d-2. Since h is odd and C is even,
b is even. Every ODD branching-defect stratum is therefore finite.
This explains the finite b5 sector and the absence of an infinite odd-b
tail without extrapolating those small cases.

## 3. The exact finite grammar and its Hamilton counts

For a bounded coloured tail(W,I), add a new root r that loses to I and
beats W-I; call this finite rooted kernel V. Define

    q0=H(W), H0=H(V), s0=number of Hamilton paths of V starting at r,
    F0=sum_(S containing r) e_V(S,r)H(V\S),
    C=F0-H0-q0+s0,     E=F0-q0.                           (4)

Thus C is delta(W), by the positive defect identity. The distinction
between W and V matters: the intercept below uses H(W), NOT H(V).
For example tail word1 has H0=5,q0=3,s0=1,F0=9, giving C2,E6 and the
ray6*2^k+1; replacing q0 by H0 would incorrectly give+3.

The exact prefix-transfer bijection in
[TAIL_SQUARES section2](TAIL_SQUARES.md#2-the-source-inventory-and-exact-prefix-transfer)
applies to every nonempty P:

    H(T)=E Q_P(2)+(q0-C)H(P).                            (5)

It counts zero, one or two initial P segments separated by the root;
there are no other entries into P. Here E>0: in F0-q0, each two-vertex
path a->r with a in I supplies a positive term. No sign assumption on
the intercept q0-C is needed.

Fix the ordered nontrivial prefix core list P1,...,Pc. Put

    h=product H(Pj),     Q*=product Q_Pj(2),

using empty products1. Insert TT_g0,...,TT_gc into its c+1 gaps, with
each gj>=0, and write m=sum gj. The ordinal product identities give

    H(P)=h,     Q_P(2)=Q* 2^m.

Thus every such family has the exact scalar law

    H(T)=A 2^m+B,
    A=E Q*>0,     B=h(q0-C),     b=h C.                   (6)

There are finitely many possible bounded tails and prefix core lists for
each b. Enumerate them as actual finite tournaments and retain exactly:

* the inport orientation and residual SCC structure;
* a nonordinary first SCC of W-I, if that set is nonempty, ensuring P is
  the MAXIMAL ordinary prefix;
* the equality h C=b;
* the condition that the last SCC of W meets I.

The last condition is exactly strongness when P is nonempty: the root
reaches all W through P, and every vertex of W reaches I and hence r
exactly when its last SCC meets I. If c=0, require m>=1 so P is nonempty.
If c>0, allow m>=0. Every listed graph is a valid source of the stated
defect, and the canonical decomposition proves every neutral source is
listed. Nonneutral sources form the separately enumerated finite part.
No capped-DFA search is required for this stronger finite grammar, though
the inherited automata remain useful small-budget enumerators and controls.

## 4. Class counts, symmetries, and scalar sparsity

For each fixed tail and ordered prefix core list, a total gap length m has
exactly binom(m+c,c) distributions. These give distinct rooted ordinary
classes: the residual SCC order and every gap position are intrinsic.
Choosing one coloured-tail isomorphism representative and one ordinary
representative for each prefix core therefore gives an exact class count,
not a count of labelled descriptions. Different core occurrences at
different positions cannot be permuted, even if isomorphic.

The rooted automorphism group is

    Aut(W,I) times product_j Aut(Pj),                    (7)

independent of gap lengths. The root is retained, the tail colouring I is
intrinsic, every transitive gap is rigid, and all prefix-to-tail incidences
are uniform. Thus the order-generating functions have denominator powers
of1-x and the class counts are eventually polynomial, now with a direct
canonical explanation. The all-neutral c0,m0 case is excluded; it is
already in the finite P-empty part if it is a valid source at all.

**Scalar sparsity corollary.** For any fixed defect budget B0, the set of
Hamilton values H<=X admitting a strong root with d>=2 and0<b<=B0 has
size O_B0(log(X+2)). Indeed there are finitely many rays(6) and finitely
many nonneutral values, and A>0 allows only O(log X) exponents per ray.
Consequently such a fixed bounded branching-defect bank has zero density
among odd targets and cannot prove all-odd target coverage by itself.

This is a scalar-spectrum obstruction, distinct from missing almost all
tournament isomorphism classes and from the finite-state witness-selector
obstruction. It does NOT exclude d1 roots, a defect growing with the
target, or a source-return operation that leaves every fixed budget.

For a fixed strong-source family, order is n=n0+m, so its exact count
can equally be written H=(A/2^n0)2^n+B. Consequently at each sufficiently
large order there are only O_b(1) distinct H values in the fixed-defect
stratum, although their rooted-class multiplicities can grow polynomially.
Every unbounded such family has n=log2 H+O_b(1). Given its finite atlas,
target inversion is exact: test whether(N-B)/A is a permitted power of2,
then enumerate the binom(m+c,c) gap distributions. This concerns the
classified defect stratum, not unrestricted tournament enumeration.

**Weak scalar extension.** The independently reviewed root-SCC factorization
in [RESPONSE_FAMILIES section4](../metagraph/RESPONSE_FAMILIES.md#4-weak-graphs-have-a-canonical-composition-law)
extends the finite union of dyadic VALUE rays, and the O_B0(log X) bound,
to all rooted tournaments with d>=2 and defect at most B0. If the root's
SCC has a nonempty upstream component, H<=b<=B0. Otherwise its strong
rooted component has defect b/t and count H/t, where t is the Hamilton
count of the downstream suffix and divides b. There are only finitely
many such t, so only finitely many multiples of the strong dyadic rays.
This does NOT extend bounded-order nonneutrality or the fixed-order
O_b(1) statement to weak sources: transitive sink padding already supplies
counterexamples to those extensions.

## 5. Arithmetic and construction boundary

At a fixed b, square-plus-one questions reduce to finitely many exact
equations n^2+1=A*2^m+B and a finite exception list. The complete solution
through b6 in TAIL_SQUARES is one instance; the general exponential
Diophantine equations are not all solved by this classification. A prime
guard remains a prime guard. In particular no Landau infinitude follows.

The map retains the bounded tail, ordered cores, gap vector and marked
root. Quotienting to(A,B) loses graph identity, automorphisms, and legal
future operations. The exact dyadic response is useful for enumeration
and exclusions, but does not identify two sources with the same H.
The uniform reset can raise defect24 to268 and beyond; that behavior
is compatible with the corollary and explains why unbounded interface
data remain relevant to the Chan–Pak construction problem.

`dyadic_defect.py` checks canonical decompositions and formulas by independent
literal path tables on all strong tournaments through five vertices with
every branching root, plus nontransitive-prefix controls. These finite
checks supplement the atom, SCC and path-bijection proofs above.

## 6. The complete rooted interface, and an actual neutral core switch

The prefix simplification preserves more than the scalar count. Use the
seven rooted coordinates(H,S,E,d,U,V,K) of the inherited
[root-amalgamation law](../../h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md):
S,E count full paths starting/ending at r; d=H(T-r); U,V count ordered
two-path covers whose distinguished nonempty path starts/ends at r,
allowing the other path empty; and K=Q_(T-r)(2).

Keep the nonempty P and W above, with d>=2 inports. For the finite root
kernel G=W+r write its seven coordinates with subscript0. Thus d0=q0=H(W)
and K0=Q_W(2). Write a0=V0-q0 for the earlier coefficient called E in(4),
so the symbol E below unambiguously denotes the endpoint coordinate.
Let h=H(P), Q=Q_P(2), and C=V0-H0-q0+S0. Then the ENTIRE response is

    H = a0 Q+(q0-C)h,
    S = d = h q0,
    E = h E0,
    U = Q K0+h(U0-K0),
    V = H+h C,
    K = Q K0.                                             (8)

The H formula is(5). A path starting at r must traverse P immediately
and then W, proving S=d=hq0. A path ending at r must traverse P first
and then a kernel path ending at r, proving E=hE0. The nonempty W is
essential to this endpoint argument and is automatic because it contains
the at least two inports. The V formula is the defect identity, and K
is the ordinal product for the deletion P->W.

For U, distinguish the two ordered path bins: the first starts at r.
If all of P lies in the other bin, its contribution is h U0. If all of
P lies in the root bin, inserting it after r removes the old first-edge
constraint and contributes h K0. If both P segments are nonempty, the
contribution is(Q-2h)K0. Their sum is the displayed U formula. These
three cases are disjoint and cover every root-start two-cover.

Therefore any prefix replacements with the same pair(H(P),Q_P(2))
preserve all seven original-root coordinates. Redistributing transitive
gap lengths with fixed total, or permuting the ordered prefix cores,
preserves that pair by ordinal multiplication. The resulting rooted
graphs generally remain distinct ordinary classes and must still be
counted separately. Equality of the seven coordinates is preserved by
every fixed iterated ordered root-amalgamation context; equality of the
four-coordinate substate also preserves every common O/Z/J continuation.
No neutral physical path implementing an arbitrary core permutation is
asserted.

There IS a physical example within the five-core bank. The
inherited [cap9 arc dictionary](../metagraph/cap9.json) contains the named
one-arc reversal

    S5a=(2,4,9,19,7) -> (6,4,8,19,7),
    reversing the old arc2->0.

The target is isomorphic to S5b=(2,4,25,3,11). Both cores have
(H,Q(2))=(9,84). Apply this single reversal INSIDE any occurrence of
that core in the ordinary prefix P. All vertices, the marked root and
every other arc remain unchanged, including every outside incidence.
Equation(8) proves that the operation preserves the full seven-coordinate
rooted response in every such prefix/tail context. It also preserves
strongness and the branching defect. This is an actual source operation
with its future sidecar retained, not merely substitution of equal scalar
counts. The uniform-prefix hypothesis is essential: arbitrary mixed
apex cuts of the two cores need not have equal responses.

`prefix_interface.py` checks eight literal contextual flip pairs, all six
core-order/gap distributions of one two-core example, and literal common
root amalgamation on both sides. The graphs and full seven-tuples are
saved in `prefix_interface.json`. The all-order statement follows from
the bijections above and the inherited exact amalgamation law.
