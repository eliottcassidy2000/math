# Critical source energy, the integer clock, and a precise Collatz boundary

**PROVED scoped / FINITE-EXACT; same-agent audit.** The critical tree is the
fixed-point tree of uniform recursive ternary automorphisms. The arithmetic
factor below is inherited from incoming Session94. These statements do not
identify that measure with the native ROOT source law or prove convergence
of a prescribed integer. The new targets are source energy and clock regression.

There are two distinct measures: mu below is a deterministic mass on tree
addresses, while expectation samples the group or the inherited Haar-coded
tree. Neither is automatically the native ROOT source distribution.

## 1. Every positive conservative source flow has growing energy

At every node of the infinite ternary tree, independently choose a uniform
local C3 rotation. A vertex is fixed iff all rotations on its path are zero.
Conditional on being fixed, a node has all three children fixed with
probability1/3 and none with probability2/3.

Let mu be a deterministic nonnegative flow, mu(root)=1 and
mu(v)=sum_(children u of v) mu(u). Set

    M_k=3^k sum_(|v|=k, v fixed) mu(v),
    e_k=3^k sum_(|v|=k) mu(v)^2.

**Theorem E1 (exact source energy).** M_k is a mean-one martingale, and

    E[(M_(k+1)-M_k)^2 | F_k]
      =2·3^(2k) sum_(|v|=k, v fixed) mu(v)^2,
    E M_K² = 1+2 sum_(k=0)^(K-1) e_k >= 1+2K.              (E1)

Proof. For each fixed v, its next contribution is3^(k+1)mu(v) with
probability1/3 and zero otherwise. Its conditional mean is3^k mu(v) and
variance2·3^(2k)mu(v)^2. Conditional contributions at different nodes are
independent. Every depth-k vertex is fixed with probability3^-k. Taking
expectations and adding the orthogonal martingale increments proves the
identity. Cauchy-Schwarz on the3^k masses summing to1 gives e_k>=1. QED.

There is also an exact LCA kernel for this energy. If v,w are depth-K
vertices and lcp(v,w) is their common-prefix length, then

    K_K(v,w)=3^(min(K,1+lcp(v,w))),
    E M_K²=sum_(v,w) mu(v)mu(w) K_K(v,w).                 (E1a)

For distinct vertices, fixing both requires zero rotations at precisely
2K-lcp(v,w)-1 independent nodes; the diagonal requires K. Multiplication
by3^(2K) proves the kernel. Equivalently,

    K_K(v,w)=1+2 sum_(j=0)^(K-1) 3^j 1[prefix_j(v)=prefix_j(w)].

This recovers(E1) by grouping source pairs at their common ancestors. It
connects the **joint** port coordinate to energy, rather than just comparing
their marginal distributions.

**A refinement-activated null mode.** At depth K, a signed mass+1,-1 on
two sibling vertices has M_K=0 for every automorphism: siblings are fixed
together. Extend each signed mass conservatively to the next generation.
Then its mean remains zero, but

    E M_(K+1)²=4·3^K.                                   (E1b)

All earlier prefix masses vanish, while the two depth-K masses have square
sum2; the same variance calculation gives(E1b). Thus a nullspace discarded
at the current depth can carry the entire next increment. Keeping its signed
lift is necessary if refinement is an allowed future operation. This is an
exact kernel-lift mechanism, not a claim that all signed currents amplify.

Equality at a given finite K means uniform mass at every contributing
level0,...,K-1. Thus positive reweighting cannot produce uniform L2 control
for this model. If mu has an atom of mass w on an infinite path, then

    E M_K² >= 1+w²(3^K-1).                                (E2)

For product branch probabilities pi, e_k=(3 sum pi_i²)^k. Uniform pi gives
1+2K; pi=(1/2,1/3,1/6) gives1+12[(7/6)^K-1]; an atomic flow gives3^K.
These examples and a flow choosing children proportionally to reciprocal
Berggren hypotenuses are checked with the independent fixed-weight recurrence.

Extinction of the unweighted fixed tree is almost sure: the least fixed
point in[0,1] of q=(2+q³)/3 is1. Hence M_k→0 almost surely while its mean
stays1. In particular this positive martingale is not uniformly integrable.
This is a statement about the displayed model, not a universal obstruction
to signed native Collatz currents.

## 2. A constructive discounted alternative, with its debt exposed

For deterministic coefficients a_k define

    T_K=1+sum_(k<K) a_k(M_(k+1)-M_k).

Orthogonality gives the exact criterion

    sup_K E T_K² < infinity  iff  sum_k a_k² e_k < infinity. (E3)

The criterion uniform over *all* probability flows is
sum_k a_k²3^k<infinity: the upper bound uses sum mu(v)^2<=1, and a single
atomic flow attains it. For a_k=rho^k, universal control holds exactly for
|rho|<1/sqrt(3); for the uniform flow the threshold is |rho|<1.

This is an explicit controlled system to compare with the target, not a
proof that discounting preserves the desired arithmetic observable. If the
target is M_K, the removed tail M_K-T_K must remain a named debt. More
generally, any signed correction Z_K with
||M_K-1+Z_K||_2<=sqrt(C) must satisfy

    ||Z_K||_2 >= sqrt(2 sum_(k<K) e_k)-sqrt(C).             (E4)

That follows from the reverse triangle inequality. A genuine ROOT transfer
therefore needs an exact map plus cancellation paying this source-specific
energy; renaming a mean or discarding its tail does not pay it.

## 3. Incoming arithmetic factor, and the best terminal-clock correction

Session94's [valuation-tree theorem](../session94_recursive_ports_crt_clocks_20261029/dynamics/VALUATION_TREE_CLOCKS.md)
works on odd 2-adic sources with singular points and their backward iterates
removed. Under conditional odd Haar, successive odd-return valuations a
are iid with P(a=j)=2^-j. Put

    b=1[a even], G=floor((a-1)/2), a=1+b+2G.

Then b and G are independent, P(b=1)=1/3, and P(G=g)=3·4^(-g-1).
Breadth-first offspring3b gives precisely the critical tree. Keeping the
ordered G-marks and unused valuation tail restores the arithmetic source.
An unmarked tree alone loses it, and the valuation clock is correlated with
offspring: E[a b]=8/9, not E[a]E[b]. These are inherited results.

Let Z_k be the frontier population and D_k the sum of valuations read at
levels0,...,k-1. Incoming Session94 proves

    E Z_k=1, E D_k=2k, Var Z_k=2k,
    Var D_k=(8/3)k²(k-1)+2k,
    Cov(Z_k,D_k)=(2/3)k(3k-2).                            (E5)

**Theorem E2 (clock-only linear compensation).** For k>=1, the unique best
real coefficient in Z_k-1-c(D_k-2k) is

    c_k=(3k-2)/(4k²-4k+3),
    min_c Var[Z_k-1-c(D_k-2k)]
      =2k(3k²+5)/[3(4k²-4k+3)] > k/2.                   (E6)

Proof. Minimize Var Z-2c Cov(Z,D)+c² Var D. Its minimizer is Cov/Var D;
substitution of(E5) gives(E6). All coefficients are exact rational numbers.
The residual is asymptotic to k/2. Thus even a depth-dependent optimal
linear use of the **actual accumulated division clock** cannot supply bounded
energy. This excludes neither nonlinear clock functions nor additional
ordered carry/port observations.

The checker independently expands the root recurrence: with z=E Z, d=E D,
z2=E Z², d2=E D², zd=E ZD, the next moments are

    z'=z, d'=d+2, z2'=z2+2z²,
    d2'=d2+2d²+(16/3)d+6, zd'=zd+2zd_mean+(8/3)z,

where zd_mean=z*d, not E ZD. This verifies(E5)-(E6) through depth12.

### A bit recovered from the clock, and its sharp limit

For any processed valuation prefix of length J, let Q be its auxiliary queue
and D its division clock. Since a_i=1+b_i+2G_i,

    Q=1+3 sum b_i-J =1+3D-4J-6 sum G_i.                   (E9)

Thus **Q=D+1 modulo2**, and retaining J as well recovers Q modulo6.
At a complete generation boundary with k>=1, J=1+3 times an integer,
so Z_k=3(D_k-1) modulo6. This is an exact shared-clock bit, not an independent
extra observation. Equivalently the joint generating function obeys
H_k(-s,-z)=-H_k(s,z). It gives a nonlinear refinement that linear regression
does not use.

It does not recover the whole frontier. At height two, the valuation words
(4,1,1,1) and(2,2,2,1) both have J=4 and D=7, but frontiers0 and6. They
are native positive-source cylinders, respectively

    x=165 mod256, endpoint=(81x+331)/128;
    x=129 mod256, endpoint=(81x+175)/128.

Their smallest positive representatives follow165→31→47→71→107 and
129→97→73→55→83. Same clock, same slope and recovered parity, different
frontier and carry. The checker verifies both literal paths and(E9) on all
4095 positive-valuation compositions of total clock at most12. Prefixes in
that identity test are algebraic words; a stopped-tree application additionally
keeps the queue-positivity guard until its first zero.

For standard Collatz the executed clock is D plus the number of processed
odd steps, not D alone. The incoming example41→31→47→71→107 has an extinct
auxiliary tree yet grows arithmetically; source1 has an infinite auxiliary
tree although its arithmetic orbit is already certified. These two hostiles
prevent replacing arithmetic entry by this auxiliary stopping event.

## 4. Lucas word pressure: an exact six, with the source law declared

For rooted cyclic no-11 words put

    P_n(t)=sum_(w in W_n) t^(number of ones in w).

The transfer matrix [[1,t],[1,0]] gives the formal trace initial values
P_0=2, P_1=1 (the word definition above uses n>=1), and
P_n=P_(n-1)+tP_(n-2). Its two eigenvalues have product -t, so

    P_(3n)(t)=P_n(t)^3-3(-t)^n P_n(t).                   (E7)

At t=1, P_n=L_n. At t=6, the eigenvalues are3 and-2. A legal standard
Collatz branch word of length n and k odd branches has affine slope
lambda=3^k/2^(n-k)=2^-n6^k. Thus

    sum_w lambda(w)=2^-n[3^n+(-2)^n],
    mean_uniform_words lambda
      =[3^n+(-2)^n]/[2^n L_n] ~ [3/(2phi)]^n.            (E8)

The fraction with lambda>=1 is bounded above by this mean. This is a finite
uniform-word statement; it supplies no law for a fixed positive integer.
The inherited rational fixed source is B_w/(2^(n-k)-3^k), with B_w>0 for
a nonzero word. Thus its positivity is exactly lambda<1, but integrality
and the actual numerator remain separate conditions.

The proper-clock hostile in(C1) joins the two lanes exactly. At period3m,
the words{2} and{0,m+1} have the same period-m field label but one versus
two odd branches. Their slopes differ by a factor **six**, while their
reflections differ by -2phi^-1. A quotient can therefore erase both an
orientation observable and the odd-charge coordinate. Restoring the field
phase alone does not authorize a change to the source measure.

## 5. A reusable research proposal

The concrete candidate is a **source-decorated refinement system**: retain
the arithmetic source and ordered carry, a typed local chart, the joint
port orbit, the accumulated transport permutation, and the executed clock.
Its transition certificate has three separate obligations:

1. An exact map/intertwiner for each operation, such as(C1)-(C2) or the
   root-framed Berggren identity(W6), including a certificate of any lost data.
2. Source mass pushed through the same joint fibres, evaluated by(W4)-(W5)
   or an explicitly retained correlated law.
3. A signed innovation estimate satisfying a summable version of(E3), with
   the compensating tail still represented and paid.

This is a proposed architecture, not an established universal theorem.
The incoming [affine torsion obstruction](../weighted_crt_clocks_20261028/AFFINE.md),
A3, rules out a literal faithful homomorphism from a nontrivial finite
marking group into positive rational affine maps. Our proposed arrows are
therefore guarded source-to-target maps, not an identification of the groups.
The sharp new helper question is whether the **native ROOT innovation** admits
an exact decomposition with summable squared source energy while its remaining
debt is certified by an independent arithmetic rank. E1 excludes positive
tree reweighting alone; E2 excludes one scalar terminal-clock correction.
They identify what extra structure a successful construction must use.
