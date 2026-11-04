# Fixed construction complexity remains sparse after arbitrary multiplication

**PROVED, all orders.** Fixed feedback-arc budgets, fixed finite transitive
clock templates, and fixed-run two-apex languages cannot generate a cofinite
Hamilton-path spectrum, even if arbitrary order joins and source/sink padding
are allowed. The obstruction survives multiplication of arbitrarily many
strong atoms. This is a construction-class limitation, not a new forbidden
value in the unrestricted tournament spectrum.

## 1. The general multiplicative theorem

Let A be a set of integers at least two, and suppose for some fixed C,d that

    #{a in A:a<=X} <= C(1+log X)^d       for all X>=2.      (1)

Let <A> be all finite products of elements of A, including the empty product
1. Then, for every epsilon>0,

    #{n in <A>:n<=X} = O_(A,epsilon)(X^epsilon).          (2)

In particular this multiplicative semigroup has natural density zero and
cannot contain all sufficiently large odd integers.

**Proof.** For every s>0, (1) implies sum_(a in A) a^(-s)<infinity, for example
by splitting into intervals [exp(j),exp(j+1)). Moreover

    Z_A(s) = product_(a in A) (1-a^(-s))^(-1) < infinity,

since -log(1-a^(-s)) <= a^(-s)/(1-2^(-s)). Expanding finite partial products
and then taking a monotone limit counts every finite exponent vector. Distinct
exponent vectors may have the same integer product, but this only overcounts:

    sum_(n in <A>) n^(-s) <= Z_A(s).

For n<=X, n^(-s)>=X^(-s), so the left side gives the bound Z_A(s) X^s.
Take s=epsilon. QED. This requires no prime-distribution theorem or unique
factorization of the generators.

The exact tournament connection is the SCC factorization: if the strong
components of T are C_1->...->C_t, every Hamilton path visits them in this
order, whence H(T)=product H(C_i). Singleton components contribute one.
Thus (2) applies to arbitrary order joins of permitted strong atoms, however
many factors are used. This is stronger than saying that a finite seed bank
misses a new prime: the bank may contain infinitely many variable-clock atoms.

## 2. Fixed feedback-arc budget, including unlimited strong factors

Let A_k be the distinct H-values of nontrivial strong tournaments whose
minimum feedback arc number is at most a fixed k>=1. Then

    #{h in A_k:h<=X} = O_k((1+log X)^(2k-1)).             (3)

Consequently arbitrary tournaments all of whose nontrivial strong components
have feedback arc number at most k have only X^{o(1)} possible H-values below
X. Their total feedback arc number may be arbitrarily large.

**Proof.** Use the cited Busch bound H(C)>=5^((n-1)/3) for a strong tournament
of order n. H(C)<=X therefore forces

    n <= L(X) = floor(1+3 log(X)/log(5)).                 (4)

Choose an ordering with at most k backward arcs. Their endpoint set has
s<=2k names. Since the tournament is strong, the first name is not a source
and the last is not a sink. Hence both belong to the backward-arc endpoint
set. For fixed n,s, the other endpoints can be selected in at most
binom(n-2,s-2) ways. On those s ordered endpoints, the backward pairs can
be chosen in at most

    B_(s,k) = sum_(j=0)^k binom(binom(s,2),j)

ways (binomial coefficients outside their range mean zero). Every such
choice and the endpoint positions specify a tournament; most choices in
this upper bound need not actually be strong or use all endpoints. No n!
factor is needed: an isomorphism class with a suitable ordering has at least
one such ordered presentation, and we only need an upper bound.

Summing n=3,...,L and s=2,...,2k gives at most

    sum_(s=2)^(2k) B_(s,k) binom(L-1,s-1),               (5)

using the hockey-stick identity with an optional extra n=2 term. This is
O_k(L^(2k-1)), proving (3). Apply section 1. QED.

For k=1 the exponent is one, consistent with the classical single-back-arc
counts 1+2^j. The known three-back-arc Mersenne family does not conflict with
this theorem: an infinite family can still have only logarithmically many
values below X. No bound on feedback *vertex* number is inferred here.

## 3. Arbitrary nonuniform clocks on finitely many fixed templates

Take finitely many fixed tournament templates. Each has a fixed finite core
and at most r designated vertices replaced by nonempty transitive modules,
with arbitrary and independent positive module lengths; all other template
names have fixed sizes. Restrict to strong outputs. Their distinct H-values
below X number O((1+log X)^r).

Indeed the total output order is bounded by (4). Each of the r lengths is
at most L(X), so there are O(L(X)^r) length vectors per template. This proof
is independent of the precise exponential response formula, of cancellations
in that formula, and of whether the clocks are equal, shifted, or unrelated.
For empty allowed modules, split into the finitely many templates remaining
after their deletion. The same bound follows.

By section 1, even arbitrary order joins of these infinitely many strong
clock outputs have X^{o(1)} possible counts below X. Allowing a finite union
of such template banks does not help.

This applies to any fixed finite bank of the S190/S191 rooted wrappers and
to any fixed number of independently varied transitive runs. It does **not**
apply when the number of runs or the template itself grows with the target,
or when newly built nontransitive objects are substituted recursively into
modules. Those operations leave the stated class.

## 4. The two-apex language must use unbounded attachment alternation

In `TWO_APEX.md`, restrict the four-letter attachment word to at most r
maximal constant runs. There are 4*3^(j-1) patterns with exactly j>=1 runs.
For chain length at most L-2, the number of positive run-length choices is
binom(L-2,j), after summing over chain lengths. Thus the number of strong
outputs with H<=X is at most

    sum_(j=1)^r 4*3^(j-1) binom(L(X)-2,j),               (6)

plus harmless empty-chain cases. This is O_r((1+log X)^r). Consequently its
arbitrary order-join closure has X^{o(1)} possible H-values.

Each fixed word pattern is indeed a transitive-module clock template:
consecutive equal letters have identical incidences to both apexes, and all
arcs to other chain runs follow their chain order. Unbounded changes of
letters are therefore a mathematically necessary escape from this particular
obstruction if the two-apex grammar is to prove cofinite spectrum. Two apexes
alone are not excluded: with an unrestricted word the number of patterns is
exponential, and hypothesis (1) need not hold.

## 5. Dependencies, credit, and scope

The analytic semigroup argument is elementary Dirichlet-series comparison;
no priority claim is made. Its role here is to close the multiplication
loophole in bounded-complexity tournament constructions. The only external
bound used is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in
Strong Tournaments*, EJC 13 (2006), N3,
https://doi.org/10.37236/1141. Its published abstract's bound was checked in
the inherited S190 packet; this note does not claim a new independent check
of the paper's full proof. SCC multiplicativity is proved directly above.

No arbitrary-target realization, global additional gap, or universal lower
bound on feedback vertex number follows. These are support-counting theorems
for exactly the named construction classes, not counts of labelled objects
or isomorphism classes.
