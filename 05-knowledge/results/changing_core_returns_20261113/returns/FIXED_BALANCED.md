# Fixed balanced deletion count has finitely many exact dyadic families

**PROVED IN NOTE, using inherited birth, SCC, and strong-order results.**
For every fixed odd integer t>=3, the full Hamilton counts H of all
strong marked tournaments with profile(H,t,t,H) form a finite union
of dyadic rays. There is an exact finite family of source templates,
including all placements of transitive gaps, and an effective membership
test once its finite core bank is enumerated.

This fixes the deletion/start coordinate t, not merely the zero
readiness defect. Allowing t to grow is outside the finite-union
conclusion. No all-target balanced selector or Chan--Pak proof follows.

## 1. Inheritance and the precise finite bank

The mechanism is inherited from the
[newborn inverse and weighted-cut identity](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md#3-the-defect-detects-indegree-and-the-legal-inverse)
and the [SCC factorization O1](../../ready_return_pentagon_20261103/OBSTRUCTION.md).
The closest construction is the existing zero-ready prefix clock in
[component inversion](../../zero_ready_component_returns_20261103/component/COMPONENT.md).
The contribution here is their complete fixed-t consequence, not a new
birth formula or a new factorization law. The different
[balanced-controller sparsity theorem](../../source_selection_closure_20261103/BALANCED_BOUNDARY.md)
fixes an operation architecture and allows t to change; it is not the
same quantified statement.

Let G_t consist of every ordered list of nontrivial strong tournaments

    (S1,...,Sj),     j>=1,     product_i H(Si)=t.          (1)

Each H(Si)>=3, so j<=floor(log_3 t). Each H(Si)<=t, and the inherited
published strong-order bound bounds the order of Si in terms of t.
Thus G_t is finite up to isomorphism. It may be empty. For each list
retain actual representatives and every possible mark r in its last
component. Put

    P_i=P_Si(2),
    A_(G,r)=(product_(i<j) P_i) E_(Sj,r)(2),
    B_G=product_i P_i,
    epsilon_j=1 if j1, and0 if j>=2.

All these are positive integers obtained from actual finite graphs.

## 2. Exact support theorem and physical constructors

The claimed set of full counts is exactly

    union_(G,r) {t+2^a A_(G,r): a>=epsilon_j}
       union_G {t+2^c B_G: c>=0}.                       (2)

**Necessity.** In a strong balanced source, F=H and q=s=t. The
weighted-cut law implies root indegree at most1; strongness excludes0.
Delete that root y and let r be its sole predecessor. The parent U has

    H(U)=t,       S_U(r)=0,       F_U(r)=H-t.             (3)

The ordered nontrivial SCCs of U form one of the lists(1). Every other
SCC has count1 and is a singleton. Conversely each such SCC list and
any finite transitive singleton gaps produce a parent with H(U)=t.

The newborn graph is strong exactly when r belongs to the final SCC
of U: every old vertex must reach r before reaching y; when this holds,
y beats every old vertex except r and can also reach r through U.
Here U has a nontrivial SCC, so the singleton boundary cannot invalidate
that last path.

If r is in Sj, no singleton gap follows it. Let a be the total number
of singleton SCCs before Sj. SCC factorization gives

    F_U(r)=2^a (product_(i<j)P_i)E_(Sj,r)(2).

When j>=2, another component already precedes r, so S_U(r)=0 for all
a>=0. When j1, it is necessary that a>=1: with no preceding component
U is strong and every vertex starts a Hamilton path. Thus this is
exactly the first union in(2).

If r is a singleton SCC, it is the final singleton. Let L>=1 be the
total number of singleton SCCs in U. Every one except r contributes
P_singleton(2)=2 before the marked component. Hence

    F_U(r)=2^(L-1) product_i P_i,

which is exactly the second union in(2).

**Sufficiency.** For a term in the first union, arrange the selected
strong components in their prescribed order, insert the chosen number
of singleton SCCs before the last strong component, and mark the
selected root in that last component. The exponent guard guarantees
another component before the mark, hence S_U(r)=0. For a term in the
second union, append c+1 singleton SCCs and mark the last. In both
cases every old vertex reaches r. Add y with r->y and y beating every
other old vertex. The newborn is strong and its profile is precisely
(t+F_U(r),t,t,t+F_U(r)), proving(2). QED.

The full source-template classification permits distributing singleton
gaps anywhere before the marked final component; their total count
alone controls H. Count support therefore forgets gap placement and
the marked core's isomorphism class. The SCC list, gap lengths, and
actual marked vertex are the sidecar for reconstructing the source.

## 3. Consequences and hostile boundaries

Normalize every amplitude in(2) to its odd part and its minimum allowed
2-adic exponent, retaining the smallest threshold for repeated odd
parts. This gives disjoint exact count rays

    H=t+o*2^e,       e>=e_min(o),                       (4)

for a finite set of positive odd o. A proposed H is accepted exactly
when its amplitude H-t has an admitted odd part and meets its exponent
threshold. Every accepted value has a physical source template from
the proof. The finite bank can be expensive to enumerate; no practical
uniform complexity bound is asserted here.

The t1 boundary has a separate exact description. Its deletion parent
has H1 and is therefore transitive. Strongness forces the selected
predecessor to be its final vertex, and start0 forces parent order m>=2.
The parent end-cover count there is2^(m-1). Hence all strong balanced
(H,1,1,H) sources have, and physically attain, exactly

    H=2^a+1,       a>=1.                               (5)

This also makes one elementary number-theoretic connection precise:
the prime members of the t1 balanced count support are exactly the
Fermat primes. Write a=2^v u with u odd. If u>1, then
(2^(2^v))^u+1 is divisible by2^(2^v)+1 and is composite. Thus primality
requires a power-of-two exponent, and every prime of that Fermat form
is in(5). This is an exact support connection, not a claim that only
the five known Fermat primes exist or that more must exist.

For fixed t, the number of distinct admitted H<=X is O_t(log X).
Including(5), the same holds with a constant depending on D when t
ranges over odd1<=t<=D. Thus keeping this deletion/start coordinate bounded cannot
supply all odd targets, all primes, or all7p targets. This is not an
obstruction for arbitrary balanced sources with unbounded t.

The keep-root birth R acts on every template by

    (H,t,t,H) -> (2H-t,t,t,2H-t),

so it advances one exponent along its exact ray. This explains the
closed continuation without implying coverage of missing odd parts.

For t3, the only nontrivial core is C3 with P2=12 and every root F3.
Formula(2) reduces to H=3+3*2^a, a>=1; the second union is already
contained in the first. For t71 the complete instance is
[BIRTH_BUDGET.md](BIRTH_BUDGET.md): eleven core presentations yield62
distinct odd-part rays. Odd part405 is absent, while odd part101 with
minimum exponent3 supplies the actual neighboring ray404*2^k+71.
The literals and normal/optimized controls for that instance are in
the linked note. No analogous finite census is claimed here for every t.
