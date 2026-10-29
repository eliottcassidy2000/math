# Critical trees compiled from periodic valuation clocks

Status: PROVED in the iid geometric source model; arithmetic application uses the inherited punctured odd 2-adic coding, not a convergence theorem. This extends the session94 ternary calculation to all odd arities and identifies the least periodic source clock.

## Inheritance and the niche recovered

The nearby source is [session94 valuation trees](../../session94_recursive_ports_crt_clocks_20261029/dynamics/VALUATION_TREE_CLOCKS.md). A less-used July source, [THM-2966, spine normal form for critical-run fair extractors](../../../../01-canon/theorems/THM-2966-spine-normal-form-for-critical-run-fair-extractors.md), records why a stopping law must retain actual source-cylinder probabilities and unconsumed boundary tails. Its coin bias is unknown and its deadline is pathwise. Neither hypothesis is silently transferred here: our valuation law is specifically P(a=j)=2^(-j), and our tree stopping is probabilistic. Classical Galton–Watson and plane-tree coefficient extraction supply the combinatorics; the new extension here is the exact periodic-clock compiler and its retained arithmetic clock law. No external-priority claim is made.

## 1. Which arities admit a periodic valuation test?

Fix h>=1. Write a=hG+R uniquely with 1<=R<=h and G>=0. Then

    P(R=r)=2^(h-r)/(2^h-1),
    P(G=g)=(1-2^(-h))2^(-hg),

and R,G are independent. These equations follow by multiplying the two displayed probabilities.

Choose S subset {1,...,h}; give a vertex d children when R belongs to S, and zero children otherwise. The branching process is critical exactly when

    sum_(r in S) 2^(h-r)=(2^h-1)/d.

For integer d>=2 this is possible if and only if d divides 2^h-1. Necessity follows because the left side is an integer; sufficiency follows from the unique h-bit binary expansion of (2^h-1)/d. The event S is unique for this h and d.

Consequently:

* No even arity is realizable by this periodic-residue test.
* Every odd d>=3 is realizable.
* Its least residue clock is h=ord_d(2).
* At a larger admissible h, the chosen periodic event is the repetition of the least-clock event. A smaller period would itself make d divide 2^h-1, proving minimality.

This is a classification of this operation menu, not of all ways to extract a Bernoulli(1/d) variable. A nonperiodic or variable-length source algorithm is outside the claim.

Examples: d=3 has h=2,S={2}; d=5 has h=4,S={3,4}; d=7 has h=3,S={3}; d=9 has h=6,S={4,5,6}; d=15 has h=4,S={4}. The appearances of 3,7,15 are precisely the singleton-residue cases d=2^h-1, not a claimed identification with tournaments or projective spaces.

## 2. Full recursive law with the division clock

Read the iid valuations in breadth-first order, starting with one active vertex. A visited vertex consumes exactly one valuation a, adds a to the division clock D, and has d children or no children by the test above. Let N be the number of visited vertices before extinction, and let

    A(z)=sum_(r not in S) (z/2)^r,
    B(z)=sum_(r in S) (z/2)^r,
    H(z)=1-(z/2)^h.

The exact bivariate probability generating series F(u,z)=E[u^N z^D] satisfies

    F=u [A(z)+B(z) F^d]/H(z).

The unique solution with zero u-constant term defines a full recurrence for its coefficients. At z=1 this is the ordinary critical offspring law with probabilities (d-1)/d and 1/d. Extinction has probability one; N has infinite mean. The recursion does not bound an individual arithmetic trajectory.

For n=dm+1 and L=(d-1)m+1, plane-tree/Lagrange extraction gives

    P(N=n,D=t)
      = binom(n,m)/n * [z^t] B(z)^m A(z)^L H(z)^(-n).

Other n have probability zero. There are binom(n,m)/n ordered full d-ary tree shapes. Each shape has m internal vertices and L leaves; multiplying their independent source series proves the same formula without formal inversion.

The coefficient formula yields an efficient finite-degree recurrence without enumerating valuation words. To compute through clock t, truncate polynomial products and use

    [z^(hg)] H(z)^(-n)=binom(n+g-1,g) 2^(-hg).

For a fixed tree shape, the internal residues are iid with law R conditioned on S, the leaf residues are iid with law conditioned on its complement, and all G variables remain independent of those residues. Thus

    E[D|N=n]=m E[R|S]+L E[R|S^c]+n h/(2^h-1),
    Var(D|N=n)=m Var(R|S)+L Var(R|S^c)
                    +n h^2 2^h/(2^h-1)^2.

Tree shape introduces no extra variance because these quantities depend only on its fixed counts m,L.

There is a simpler universal form. Write mu_I=E[R|S], mu_L=E[R|S^c],
v_I=Var(R|S), v_L=Var(R|S^c), and Delta=mu_I-mu_L. Then

    E[D|N=n]=2n-Delta/d,
    Var(D|N=n)=n [2-(d-1)Delta^2/d^2]+(v_L-v_I)/d.

To prove this, substitute m=(n-1)/d,L=((d-1)n+1)/d above. The
unconditioned valuation mean and variance are both two. The law of total
variance splits the latter into the within-type term and the between-type
term (d-1)Delta^2/d^2. Conditioning on a tree fixes the counts of the two
types, removing exactly that between-type contribution to the linear
variance coefficient. The constant boundary correction remains. For d=5,
h=4, the two conditional residue means are10/3 and4/3, and both variances
are2/9; hence E[D|N=n]=2n-2/5 and Var(D|N=n)=34n/25.


## 3. Singleton-residue specialization: 3 -> 7 -> 15

For d=2^h-1, S={h}, the integer-coefficient formula simplifies to

    P(N=dm+1,D=t)
      = binom(n,m)/n * 2^(-t)
        * [z^(t-hm)] (z+...+z^(h-1))^L (1-z^h)^(-n).

The exact conditional mean is

    E[D|N=n]=2n + [d(2-h)-h]/[d(d-1)].

For h=2,3,4 these offsets are -1/3,-5/21,-17/105. The mean is near 2n because the unconditioned geometric valuation has mean two; conditioning leaves a nonzero boundary correction, not a free clock.

At h=2 every leaf residue is one, so D is necessarily odd at extinction. At h>=3 a one-vertex tree can terminate with residue one or two. Therefore no universal terminal parity survives. A parity shortcut from the ternary case cannot be transported through 7 and 15.

## 4. Source restoration and arithmetic limits

A tree alone retains only whether R is in S. Lossless reconstruction requires R at every visited vertex, G at every visited vertex, the breadth-first order, and the untouched source tail. Knowing only N and D cannot recover the sequence, its affine carry, or its guard cylinder. Even when S is a singleton, the leaf residues must be retained once h>=3.

For accelerated 3x+1, after valuations (a_1,...,a_N), maintain the native affine numerator and clock by

    B_0=0, D_0=0;
    B_(i+1)=3 B_i+2^(D_i), D_(i+1)=D_i+a_(i+1),
    x_N=(3^N x_0+B_N)/2^D_N.

The exact source cylinder is an odd class modulo 2^(D_N+1), not merely the endpoint rational expression. The standard unaccelerated clock counts N+D_N operations. For 3x-1, change the sign of each additive term; other twisted maps need their own typed recurrence.

For each fixed odd d, finite extinction-cylinder unions have natural densities equal to their Haar probabilities. Their masses increase to one, so auxiliary extinction has natural density one among positive odd 3x+1 sources. This uses finite-cylinder lower bounds, not countable additivity of natural density. It does not establish native descent or convergence. At d=3, the inherited source41 and valuation word(2,1,1,1) stop the auxiliary tree while41 ->31 ->47 ->71 ->107 grows. Changing d changes this auxiliary observation; it is not an orbit invariant.

## 5. Shared source clocks do not become independent

For several arities d_j, use a common residue clock H=lcm(h_j), with
h_j=ord_(d_j)(2). The exact one-source-step generating kernel is

    K(z;x_1,...,x_s)
      = [sum_(r=1)^H (z/2)^r product_j x_j^(d_j I_j(r))]
        / [1-(z/2)^H],

where I_j(r) tests the appropriate residue subset modulo h_j. With queue
coordinate Q_j starting at one, its unabsorbed update is Q_j <- Q_j-1+d_j I_j.
A stopped coordinate is held at zero in an absorbed joint calculation.
The common source index and division increment a remain explicit. This
kernel gives an exact multivariate transfer recurrence; it does not make
the separate breadth-first trees independent or equate their generation
clocks.

In the singleton-residue case d_h=2^h-1, d_k=2^k-1,

    P(h divides a and k divides a)=1/(2^lcm(h,k)-1),
    Cov(d_h I_h,d_k I_k)=d_h d_k/(2^lcm(h,k)-1)-1.

The arities3 and7 have covariance -2/3, despite coprime residue clocks2
and3. Arity3 with15 instead has covariance2 because the events nest.
CRT is a bijection of residue addresses, not a factorization of this
geometrically weighted source measure. Replacing the shared clock by
independent marginal branching laws changes the model.

## 6. Verification and stopping boundary

`check_residue_branching.py` independently compares the coefficient formula with a frontier-state dynamic program retaining both visited count and division clock. It checks the unique residue selector, minimal period, conditional mean algebra, the ternary parity specialization and its failure from h=3, and native source guards. These checks instantiate the proofs; they do not certify every integer orbit.

The older unknown-bias critical-run deadline problem remains a different problem. A proof for this fixed Haar law cannot be advertised as a uniformly fair extractor or a pathwise deadline result.

Classical total-progeny provenance: Meyer Dwass, [The total progeny in a branching process and a related random walk](https://doi.org/10.1017/S0021900200026711), Journal of Applied Probability6(3),1969. The tree-shape argument above supplies the specialized coefficient calculation directly.
