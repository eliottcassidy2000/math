# The exact Hamilton-count spectrum of paired-black leaves

Status: PROVED consequences of the inherited all-order leaf classification;
FINITE-EXACT independent small-order controls. No all-tournament spectrum,
global minimum-order, bridge classification, or literature-priority claim.

## 1. What is being merged

Here 14, 15 and 27 mean the prior research sessions, not tournament orders.
The request's bare numbers were ambiguous; this interpretation follows the
actual linked mechanisms in their published packets.

- Session 14 classifies non-self-converse leaves by nonpalindromic words in
  singleton and directed-triangle strong components. Its mixed leaves are
  rigid mirrored words around an unknown strong deficit-two center.
- Session 15 proves that the center is uniquely the strong four-vertex
  tournament S4. Its proof uses an all-order cutoff followed by an exhaustive
  certified finite exclusion, rather than extrapolating a small census.
- Session 27's repunit report retains automorphism weights in ordered-SCC
  enumeration. This supplies the distinction between uniform classes and
  uniform labelled tournaments that becomes essential below. Its marked
  response theorem separately explains how context restores information lost
  by bare H; the spectrum here is not a reconstruction invariant.

These are inherited inputs, not new classifications. The new deductions are
the joint count spectrum and multiplicities, its exact inverse/extremal
boundaries, and the measure-dependent asymptotic law. Session 24, recovered
in Session 27, and Session 207 already show that almost every fixed-H node
for H a power of three is a leaf. We do not claim that statement as new.

Pinned primary sources and their hashes are in DEPENDENCIES.json. The older
Session 14 text deliberately leaves the center open; the final Session 15
README and completed audits supersede that historical boundary.

## 2. Object, convention and inherited classification

Let Q_n be the fundamental paired-black, converse-merged tournament
multigraph. Every literal black-line multiplicity is divided by two through
reflection pairing; loops contribute two to degree. A Q-leaf means a
degree-one vertex of THIS quotient, not a leaf of the labelled native graph,
not a single-arc-flip graph, and not merely one distinct neighbor.

For an ordinary tournament T let H(T) count directed Hamilton vertex words
and let f(T)=H(T)/|Aut(T)|. For a non-self-converse Q-node, degree=f(T).
For a self-converse node its degree is (f(T)-B(T))/2, where B counts the
reflection-fixed Hamilton-path orbits. H is unchanged by converse and so is
well-defined on Q-nodes.

Write 1 for a singleton SCC and 3 for a directed triangle SCC, using
concatenation for dominance join. Every Q-leaf is uniquely of one of these
types:

  NS: a nonpalindromic word W in {1,3}, modulo W~reverse(W);
  mixed: A -> S4 -> A^op, where A is an arbitrary {1,3} word.

For mixed leaves, the half-word A is unique and is NOT taken modulo reversal.
SCC order is canonical. The singleton has H=Aut=1, C3 has H=Aut=3, and
S4 has H=5, Aut=1, B=3. In particular NS leaves have f=1, B=0;
mixed leaves have f=5, B=3, regardless of their order.

The source-to-target map is a tournament to its ordered SCC word, then to
(n,H). The second step forgets word order and labels. The whole word is the
sidecar restoring the isomorphism class here; marked response data restores
the labelled tournament in Session 27's different reconstruction problem.

## 3. Exact joint generating function and coefficients

Let L(z,x,y) enumerate Q-leaf vertices, marking order by z, v_3(H) by x,
and v_5(H) by y. Define

  A(z,x)=1/(1-z-x z^3),
  P(z,x)=(1+z+x z^3)/(1-z^2-x^2 z^6).

Then, as a formal power series,

  L(z,x,y) = (A(z,x)-P(z,x))/2
             + y z^4/(1-z^2-x^2 z^6).                         (1)

Proof. A enumerates all rigid SCC words. A palindrome consists of a half
word, its reverse, and either no central letter, one singleton, or one C3.
This gives P, with the half-word's order and triangle count doubled.
Deleting palindromes and dividing the remaining reversal orbits by two
gives the NS term. A mixed leaf has an S4 center and a free half-word,
giving the last term. These cases are disjoint and exhaustive. QED.

At x=y=1 equation (1) specializes to Session 15's leaf series

  (2z^4-z^5-z^7)/[(1-z-z^3)(1-z^2-z^6)].

For k,s>=0 set

  p(s,k)=0 if s and k are both odd;
  p(s,k)=binom(floor(s/2)+floor(k/2),floor(k/2)) otherwise.

The number of Q-leaves with order n=3k+s and H=3^k is exactly

  [binom(s+k,k)-p(s,k)]/2.                                    (2)

There are no such leaves for k=0 or s=0. The formula is zero there, as it
must be. The parity condition on p says that a palindrome has at most one
letter type with odd multiplicity; its first half then determines it.

The number of Q-leaves with H=5*9^k and order n=6k+4+2s is exactly

  binom(s+k,k),                                               (3)

and is zero at all other orders. Equations (2) and (3) count quotient
vertices, not labelled tournaments and not edges.

## 4. Exact spectrum, inverse boundary and largest leaf count

The all-order Q-leaf Hamilton-count spectrum is exactly

  {3^k : k>=1} union {5*9^k : k>=0}.                          (4)

Both necessity and sufficiency follow from the classification. For the
first family, one singleton and k triangles arranged as 1,3,...,3 give a
nonpalindromic word. For the second, take the half-word to be k triangles.
More precisely, a power 3^k occurs in a leaf exactly at orders n>=3k+1;
5*9^k occurs exactly at even orders n>=6k+4.

Thus the minimum order AMONG Q-LEAF REALIZATIONS is

  ell(3^k)=3k+1 (k>=1),
  ell(5*9^k)=6k+4 (k>=0),
  ell(h)=infinity for every other positive h.                 (5)

This ell is not the global inverse m(h). For example C3 has H=3 at order
three but is pure blue, not a Q-leaf, so m(3)=3 while ell(3)=4.
Likewise the entire balanced family H=15^k, k>=1, is absent from the leaf
spectrum: the first branch has no factor5; the second has v5=1 and even
v3, incompatible with v3=v5=k. Consequently Session 208's order-fourteen
H=15^6 atom cannot lie in the Q-leaf stratum. This is a topology-of-the-
quotient obstruction, not nonexistence of those counts as tournaments.

For n>=4 let H_max^leaf(n) be the largest H on a Q-leaf. Then

  n odd: H_max^leaf(n)=3^floor((n-1)/3);
  n=2 mod6: H_max^leaf(n)=3^((n-2)/3);
  n=0 or4 mod6: H_max^leaf(n)=5*9^floor((n-4)/6).              (6)

Proof. In the NS family the largest admissible k is floor((n-1)/3).
In the mixed family n must be even and k<=floor((n-4)/6).
Comparing these maxima gives mixed/NS ratio 5/3 for n=0,4 mod6 and
5/9 for n=2 mod6. Odd orders have no mixed leaves. QED.

An apparently plausible extension, 'almost every leaf is NS, so an NS leaf
maximizes H for all sufficiently large n', is therefore false at every
even order n=0 or4 mod6. Exponential rarity does not remove extremizers.

## 5. The automorphism-weighted transfer

Let E(z,x,y) have coefficient at z^n x^k y^j equal to 1/n! times the number
of labelled n-vertex tournaments with Q-image a leaf and H=3^k5^j. Then

  E(z,x,y) = A(z,x/3)-P(z,x/3)
              + y z^4/(1-z^2-(x/3)^2 z^6).                 (7)

Proof. An NS Q-node represents TWO ordinary isomorphism classes, each with
Aut order3^k, hence contributes 2n!/3^k labelled tournaments. A mixed node
represents one ordinary class with Aut order9^k and contributes n!/9^k.
Apply these weights to (1). There is no permutation of equal SCCs because
an automorphism fixes their total order. This is exactly Session 27's
weighted SCC mechanism. QED.

For example the two Q-leaves at n=4 consist of an NS pair with H=3,
containing16 labelled tournaments in total, and S4 with H=5, containing24.
They share one paired-black edge. Treating Q-node count, ordinary-class
count, labelled count, or bridge count as the same measure fails already
at this first nonempty order.

## 6. Measure-dependent exponential moments and Gaussian law

Consider either of these probability distributions for each n>=4:

  U: uniformly choose a Q_n-leaf;
  V: uniformly choose a labelled n-vertex tournament whose Q-image is a leaf.

Let K_n=v3(H). Put a=1 for U and a=1/3 for V, and let rho(a) be the unique
positive solution of rho+a rho^3=1. Define

  mu(a)=(1-rho(a))/(3-2rho(a)),
  v(a)=rho(a)(1-rho(a))/(3-2rho(a))^3.                       (8)

For either measure, as n tends to infinity,

  E[K_n]=mu(a)n+O(1),       Var(K_n)=v(a)n+O(1),
  (K_n-mu(a)n)/sqrt(v(a)n) converges to N(0,1).              (9)

The same statement for log H multiplies the mean by log3 and the variance
by (log3)^2. Write log H=K_n log3+J_n log5, where J_n indicates a mixed
leaf. The subdominant-pole proof below gives Pr(J_n=1)=O(eta^n). Since
0<=K_n<=n/3, Cov(K_n,J_n)=O(n eta^n) and Var(J_n)=O(eta^n), so this
correction preserves the O(1) moment errors. Its boundedness suffices for
the stated normal limit.
In particular the typical exponential H-growth rates are DIFFERENT:

  U: K_n/n -> 0.19425400402459; log H/n -> 0.21340983594440;
  V: K_n/n -> 0.13357525121498; log H/n -> 0.14674741244670.

The constants displayed are numerical evaluations; (8) defines them exactly.
Convergence of the displayed densities holds in probability. These are not
unconditioned random tournaments or fixed-H ensembles.

Here is also the full fixed-real moment asymptotic. For each real t, under
the applicable measure,

  E[H^t] = [(3-2rho(a))/(3-2rho(a*3^t))]
           * [rho(a)/rho(a*3^t)]^n * (1+O(eta_t^n)),        (10)

for some eta_t in (0,1). The exponentially small error can depend on t;
no uniform assertion as |t| tends to infinity is made.

Proof without a probability limit theorem as a black box. For u>0,
A(z,u) has its unique smallest-modulus pole at rho(u). Indeed if
1=z+u z^3 and |z|<=rho, the triangle inequality forces equality, |z|=rho,
and both z and u z^3 to be positive real; therefore z=rho. The pole is
simple. Its coefficient asymptotic is

  [z^n] A(z,u)=rho(u)^(-n)/(3-2rho(u)) * (1+O(eta^n)).

The other denominator in (1) or (7) is analytic on a strictly larger disk:
at z=rho(u),

  rho(u)^2+u^2 rho(u)^6 = 1-2u rho(u)^4 < 1.

This also bounds its complex roots away by the triangle inequality. Thus
reflection terms and mixed leaves are exponentially smaller than A for
every fixed positive u. The multiplier of A is1/2 for U and1 for V; it
cancels on normalization. Substitution u=a*3^t and y=5^t proves (10).

For a rigorous neighborhood argument, the simple root depends analytically
on complex u near any fixed positive value. Compactness and the strict root
separation give the coefficient estimates uniformly in a sufficiently
small complex neighborhood. In the moment generating function for K_n,
replace u by a exp(t), with y=1. Its logarithm is

  n f(t)+g(t)+O(eta^n),
  f(t)=log rho(a)-log rho(a exp(t)),
  g(t)=log(3-2rho(a))-log(3-2rho(a exp(t))).

The error is analytic and locally uniform, so derivatives are controlled
by Cauchy's estimate. Implicit differentiation gives

  d rho(a exp(t))/dt at0 = -rho(1-rho)/(3-2rho),
  f'(0)=mu(a),   f''(0)=v(a)>0.

This proves the moment assertions. For fixed real s, set t=is/sqrt(vn),
multiply by exp(-is mu n/sqrt(vn)), and expand f to second order. The
linear terms cancel, the quadratic term is -s^2/2, and the higher terms
are O(n^(-1/2)); g contributes o(1). The characteristic functions therefore
converge to exp(-s^2/2), proving the normal limit. QED.

## 7. Scope, checks and next question

The deterministic checks enumerate rigid words through order22, compare
every (n,H) coefficient with (2)/(3), test all-order formulas through order60,
and independently enumerate ALL labelled tournaments through order5. The
last check constructs automorphism orbits of literal Hamilton paths and
their anti-automorphism reversal action to measure Q-degree; it does not
assume leaf classification to decide which small tournaments are leaves.
Explicit larger positive and hostile examples use independent subset-DP H
counts. Exact integer/rational moment recurrences through n=600 check the
asymptotic constants numerically; they do not prove the limits.

The strongest live distinction is now precise: fixed powers of3 concentrate
on leaves as n grows, but among all leaves at growing n the typical H itself
grows exponentially, and its rate changes under automorphism weighting.
Neither statement is a distribution law on strong tournaments. The mixed
family is exponentially rare under both measures, yet dominates maximum H
on two even-order residue classes. The remaining question is what similar
joint laws hold at paired-black degrees2 and3, where Session14's rigid
outer-word reduction begins to change.
