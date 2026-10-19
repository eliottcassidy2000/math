# A symbolic Gold family with a golden-ratio boundary

Date: 2026-10-19. Status: **PROVED**. This is a classification of one fixed strategy over a structural family and all its specified contexts, not a classification of all Gold posets. The general Gold Partition Conjecture remains **OPEN**.

## Main theorem

Take two disjoint chains of lengths m,n>=2. Write them as

    a_1<...<a_m,  c_1<...<c_n.

Optionally put any finite poset H above both chains, with no other restrictions on H. Add vertices z,p with only the extra relations

    a_1<p, z<p.

In particular p,z have no prescribed relations to H. Let a=a_1 and b=p be the boundary ports. The earlier nine-element core is the case m=n=3 and |H|=1.

Use the fixed comparisons

    first: c_n versus a_m;
    second: c_(n-1) versus a_m.                              (1)

Give every core extension sigma a nonnegative flag weight

    alpha+beta 1{a first}+gamma 1{p last}
          +delta 1{a first, p last},

where the four coefficients are nonnegative and not all zero.

**Theorem.** Strategy (1) is Gold for every such weight if and only if

    m <= n,
    Q(m,n):=n^2-mn-(m-1)(m-2) <= 0.                        (2)

Under (2), the first comparison is at least (m-1)/(m+n-1) balanced. Its probability of the increasing answer lies in the sharp interval

    [(m-1)/(m+n-1), m/(m+n)].                              (3)

More generally, for arbitrary m,n>=2 the minimum normalized Gold deficit over the whole flag cone is exactly

    Delta(m,n) = min{
        (n-m)/(m+n),
        (m-1)(m-2)/[(m+n-1)(m+n-2)],
        -Q(m,n)/[(m+n-1)(m+n-2)]
    }.                                                    (4)

These formulas are independent of the size and structure of H.

By the arbitrary terminal-upset/initial-ideal fibre theorem in `gold_unbounded_flag_context_20261019.md`, every result transfers to arbitrary exterior posets attached through such boundary sets. Indeed their core weights are proportional to (1+r I)(1+s J), with r,s>=0. Conversely every failure in (2) is approached by genuine finite chain collars and eventually causes this strategy to fail. Thus (2) is also the exact parameter range for the strategy to work in every allowed genuine context.

Section 7 replaces the fork by an arbitrary side poset and proves the same universal lower margin by a flag-preserving interleaving injection. Thus the structural result extends well beyond fork attachments.

A concrete interior case is m=4,n=5,H empty: an eleven-element core has

    t_0-t_1-t_2 >= t_0/56

under every nonnegative flag weight and every allowed exterior context. Its first comparison is at least 3/8 balanced. This is a strict positive normalized margin even when the exceptional boundary sets have unbounded size.

## 1. Triangular cells and conditional shuffles

The first increasing answer c_n<a_m implies the second increasing answer c_(n-1)<a_m. Therefore the joint table has form

    [[A,0],[B,C]],

where A means zero c-elements after a_m, B means exactly one, and C means at least two. The three realized Gold inequalities are

    B+C-A >= 0, A-B >= 0, A-C >= 0.                        (5)

Delete z,p. The remaining extensions consist of a uniform shuffle of the two chains, followed by an arbitrary extension of H. Put T=m+n and h=|H|. Every shuffle has e(H) choices for its final cap order; this common factor will cancel.

Let t be the position of a_1 in the T-letter shuffle. Necessarily 1<=t<=n+1. Conditional on t, the prefix consists of t-1 c-elements, then a_1, followed by a uniform shuffle of chains of lengths

    r=m-1, s=n-t+1.                                       (6)

For any two-chain shuffle with lengths r>=1,s>=0, the probabilities of zero and exactly one c-element following its last a-element are

    p_A=r/(r+s),
    p_B=rs/[(r+s)(r+s-1)],

with p_B=0 if s=0. Write p_C=1-p_A-p_B.

As s decreases, both p_A-p_B and p_A-p_C are nondecreasing. For the first difference this follows from

    p_A-p_B=r(r-1)/[(r+s)(r+s-1)],

with the one-letter case interpreted directly. For the second, couple uniform orders using independent continuous labels on r a-elements and s c-elements. Adding a c-element can only increase the number of c-labels above the largest a-label. The function equal to 1 for count zero, 0 for count one, and -1 for count at least two is decreasing. Its expectation p_A-p_C therefore decreases with s.

Consequently both lower inequalities A>=B and A>=C are weakest in the conditional measure t=1, where s=n.

## 2. How the fork changes the measure

Put N=T+h. Fix one backbone extension after deleting z,p, and insert z in gap j in {0,...,N}. The number of allowed slots for p is

    W(t,j)=N+1-t if j<t, and N+1-j if j>=t.

Hence the total insertion weight is

    K_h(t)=(N+1)(N+2)/2 - t(t+1)/2.                        (7)

This is positive and strictly decreasing in t.

The four flag measures have particularly simple descriptions:

- No flags: weight the uniform shuffle by K_h(t).
- a first: condition on t=1 and forbid z before a. The remaining insertion weight N(N+1)/2 is constant.
- p last: replace W(t,j) by one. All N+1 gaps are allowed, so the shuffle is uniform.
- Both flags: condition on t=1, and use its N allowed gaps. The resulting shuffle measure is the same as the a-first measure.

This is the structural reason a four-flag proof is possible. Two flag measures coincide after normalization, a third is uniform, and the fourth is a decreasing one-coordinate reweighting.

## 3. Proof of the exact margin

Normalize each table to total one. Its realized deficits are 1-2p_A, p_A-p_B, and p_A-p_C.

### The two lower deficits

Section 1 shows that any mixture over t has the lower bounds attained by the a-first measure:

    p_A-p_B >= (m-1)(m-2)/[(T-1)(T-2)],
    p_A-p_C >= [mn+(m-1)(m-2)-n^2]/[(T-1)(T-2)].           (8)

This covers the unflagged and p-last measures, and equality holds for each a-first measure. No assumption (2) was needed to establish these bounds.

### The upper deficit and reverse covariance

Conditional on t, the first answer probability is

    f(t)=(m-1)/(T-t),

which increases with t. Under the uniform shuffle its mean is m/T. The unflagged measure reweights by decreasing K_h(t), so

    E[K_h f]/E[K_h] <= E[f] = m/T.                         (9)

For completeness, the covariance inequality follows by taking two independent copies t,t': the product

    (K_h(t)-K_h(t'))(f(t)-f(t'))

is nonpositive pointwise, and its expectation equals twice the covariance. The a-first probability (m-1)/(T-1) is also at most m/T. Thus every basis measure has

    1-2p_A >= (n-m)/T,                                    (10)

with equality in the p-last basis measure.

Every normalized nonnegative flag mixture is a convex combination of the normalized basis measures. Equations (8) and (10) therefore prove the lower bound (4). Each of its three candidate values is attained in a basis measure, so it is exact. Since the middle candidate is nonnegative for m>=2, condition Delta>=0 is exactly (2).

The same reasoning gives (3): the conditional probability f(t) is minimized at t=1, while (9) controls the unflagged upper endpoint. Both endpoints are attained by basis measures. Under (2), m/T<=1/2, proving the balance statement. In particular (2) implies n<=2m-2, so every member is at least 1/3 balanced directly from this one fixed comparison.

For the original fork construction in genuine finite contexts, the unflagged coefficient is positive. Its three realized deficits are strictly positive whenever (2) holds: the covariance is strict when m=n, and conditioning includes values t>1 that strictly improve each possibly zero lower bound. Thus the fixed strategy is strict even at the boundary Delta=0, although a positive uniform normalized margin need not survive there.

If m>n, the p-last basis has negative upper deficit; long right chain collars with no left collar approach this basis. If Q>0, the a-first basis has negative A-C; long left chain collars with no right collar approach it. This proves the claimed necessity over actual exterior contexts as well as over abstract weights.

## 4. Closed cell counts and an eleven-element certificate

The proof above is probabilistic and does not require closed counts. For an independently checkable formula, use binomial coefficients equal to zero outside their usual range, and put

    K=(N+1)(N+2)/2.

After removing the common factor e(H), the unflagged table has

    Z_0=K binom(T,m)-binom(T+2,m+2),
    A_0=K binom(T-1,m-1)-binom(T+1,m+1),
    B_0=K binom(T-2,m-1)-binom(T,m+1),
    C_0=Z_0-A_0-B_0.                                      (11)

The first subtracted binomial follows from the elementary convolution

    sum_t binom(t+1,2) binom(T-t,m-1)=binom(T+2,m+2).

Conditioning on the last a-element being last or penultimate gives the other two instances of the same identity.

Define two raw shuffle tables by

    Y=[[binom(T-2,m-2),0],
       [binom(T-3,m-2),binom(T-1,m-1)-binom(T-2,m-2)-binom(T-3,m-2)]],

    U=[[binom(T-1,m-1),0],
       [binom(T-2,m-1),binom(T,m)-binom(T-1,m-1)-binom(T-2,m-1)]].

Then

    F_a=N(N+1)Y/2, F_p=(N+1)U, F_ap=NY.                    (12)

For m=4,n=5,h=0, the four exact tables are

    F_0  = [[2828,0],[1799,1841]], total 6468;
    F_a  = [[945,0],[675,900]],    total 2520;
    F_p  = [[560,0],[350,350]],    total 1260;
    F_ap = [[189,0],[135,180]],    total 504.

The sharp margin is 1/56, attained by A-C in either a-first basis measure. Arbitrary genuine contexts approach that bound through long left collars. The construction gives infinitely many core sizes and arbitrary exterior posets, rather than only one finite certificate.

## 5. A genuine golden-ratio and Fibonacci boundary

Solving (2) for n gives

    m <= n <= [m+sqrt(5m^2-12m+8)]/2.                      (13)

Thus the limiting admissible aspect ratio n/m is the golden ratio phi=(1+sqrt(5))/2. This occurrence comes from an exact extension-count inequality: too many trailing c-elements make the a-first flag measure violate A>=C. It is not an analogy inferred from the name of the conjecture.

Let F_k be Fibonacci numbers with F_0=0,F_1=1. Cassini's identity gives

    Q(F_k+1,F_(k+1)) = (-1)^k-F_(k-1) <= 0                (k>=3).

These parameters satisfy m<=n and approach the golden-ratio edge. Hence Fibonacci data supply an explicit infinite sequence of permitted structural cores.

The curved equality boundary Q(m,n)=0 has the Pell form

    (5m-6)^2 - 5(2n-m)^2 = -4.                            (14)

For Lucas numbers L_0=2,L_1=1, the identity L_k^2-5F_k^2=4(-1)^k gives the explicit boundary sequence

    m=(L_(4j+3)+6)/5,
    n=(F_(4j+3)+m)/2,      j>=0,

beginning (m,n)=(2,2),(7,10),(41,65),(274,442). Integrality and parity follow immediately from the Fibonacci/Lucas recurrence modulo 5 and 2. These are exactly the integer equality cases satisfying m,n>=2 and m<=n; the reduction proof follows below.

### Exact boundary recurrence and its budget-matrix link

The diagonal n=m also has Delta=0, usually with Q<0. The following recurrence and exhaustion apply specifically to Q=0.

The affine map

    R(m,n)=(2m+3n-3, 3m+5n-6)                              (15)

preserves Q. It generates the displayed sequence from (2,2). Its inverse is

    R^(-1)(m,n)=(5m-3n-3, -3m+2n+3).

This generates every admissible integer equality case. To prove exhaustion, suppose Q=0 and m>=7. As a polynomial in n, Q is increasing for n>=m. Evaluating it at the two indicated endpoints gives

    Q(m,(8m-6)/5)=-(m-7)(m-2)/25 <= 0,
    Q(m,(5m-5)/3)=(m-1)(m-7)/9 >= 0.

Therefore (8m-6)/5<=n<=(5m-5)/3. The inverse coordinates satisfy 2<=m'<=n', preserve Q, and strictly decrease both coordinates; in particular m'<=(m+3)/5<m. Repeating reaches 2<=m<=6. The five discriminants 5m^2-12m+8 are 4,17,40,73,116, so only (2,2) survives. This proves exhaustion by (15), without appealing to a general Pell-solution theorem.

There is also an exact link to the inherited Gold equality-obligation map

    g(x)=(1-2x)/(1-x),   M=[[-2,1],[-1,1]].

Let

    F=[[0,1],[1,1]], B=[[0,1],[-1,1]],
    S=[[2,3],[3,5]]=F^4.

Direct multiplication gives B^(-1) M B=-F, hence B^(-1) M^4 B=S. The affine recurrence (15) has fixed point c=(6/5,3/5), and is exactly

    R(v)=c+S(v-c).

Thus after translation by c and the stated linear change of basis, equality-boundary evolution is the fourth power of the Gold-budget matrix. The Lucas/Fibonacci formula also follows from this recurrence. This is a statement about equality parameters and a shared exact matrix action; it does not say that conditional extension probabilities evolve by g.

These integer identities and the golden asymptote are mathematically linked to the Gold certificate by (2). They do not identify the comparison process with a Fibonacci torus action or a Collatz orbit, and no finite-field reduction is used to transport inequalities.

## 6. Independent exact controls and scope

`gold_two_chain_fork.py` compares (11)–(12) against actual ideal-DP extension counts for every 2<=m,n<=8 and six choices of cap H, including nonchains. It checks the classification, balance interval, and sharp margin for every 2<=m,n<=80. It also tests genuine unbounded chain collars around the eleven-element core, including 71-vertex examples, the arbitrary-side-poset extension below, and exact Pell-boundary parameters. All arithmetic is integral or rational. Normal and `python3 -O` outputs agree; the published JSON receipt is in `05-knowledge/results`. The script prints a fresh receipt without changing files.

The conceptual proof is the conditional-shuffle law, the monotonicity of the number of trailing chain elements, and reverse covariance under a decreasing insertion weight. These yield an extendable local certificate schema. They do not prove every poset has such a backbone, nor Gold for all width-three or bounded-incomparability-degree posets.

Primary conjecture and sufficient-condition reference: Marcin Peczarski, *Comments on the Golden Partition Conjecture*, Contributions to Discrete Mathematics 12 (2017), 106–109:
https://cdm.ucalgary.ca/article/download/62349/46785/178291


## 7. An arbitrary side-poset extension via a flag-preserving injection

The fork is not essential. Keep the two chains and arbitrary upper cap H. Take any nonempty finite side poset D, any upset W of D, and any maximal element p of D. Add exactly the cross-relations

    a_1<w for w in W,

and no other relations between D and the backbone. The upset assumption makes this transitive. The boundary ports remain a=a_1 and p. They are distinct, minimal and maximal respectively. The original fork is D={z<p}, W={p}.

**Side-poset theorem.** For this entire class, every nonnegative mixture of the four a-first/p-last flag measures satisfies the normalized Gold lower bound Delta(m,n) from (4), and the probability interval (3). In particular (2) implies Gold under every arbitrary exterior upset/ideal collar. The bound and admissible parameter range are sharp uniformly over this structural class because the original fork is a member. Boundary strictness is not claimed for general side posets: when W is empty the fibre can be constant, and m=n can realize zero slack. Interior cases such as m=4,n=5 retain their positive uniform 1/56 bound.

**Proof.** Fix a backbone extension, an extension tau of D, and a binary word specifying their interleaving. All dependence on the backbone order is through the position t of a_1. The only cross-constraints say that the occurrence of a_1 must precede each tau-element in W.

Now replace a backbone order whose first a occurs at t by one whose first a occurs earlier at t'. Keep tau and the binary word unchanged. The physical position of a_1 can only move earlier, so every previously satisfied cross-constraint remains satisfied. This is an injection of the allowed interleaving-word sets at t into those at t'. It preserves the p-last flag: that flag merely requires tau to end in p and the binary word to end with a D-symbol. Summing over tau proves that both the unflagged restriction weight w(t) and the p-last weight w_p(t) are nonnegative, positive, and nonincreasing in t.

For the a-first flag, necessarily t=1 and the binary word begins with a backbone symbol. The remaining interleaving count is independent of the order of the other backbone vertices, so the induced shuffle is uniform conditional on t=1. The same holds when both flags are imposed.

The lower-deficit bounds (8) apply to any mixture over t. The reverse-covariance argument (9) applies to every nonincreasing weight, not only to the quadratic insertion weight (7). Therefore all four measures obey exactly the lower bounds used to prove (4), and their mixtures do too. The earlier boundary fibre theorem then transfers the certificate through arbitrary exterior collars. ∎

The side poset may have arbitrarily large width, arbitrary internal comparabilities, and a proper nontrivial upset W. The proof needs its attachment to be through this one first-arrival constraint; additional cross-relations involving later chain vertices need not preserve the decreasing-weight property. This is a concrete structural limitation rather than an assertion about every poset.


Reproduce from the repository root with `python3 05-knowledge/experiments/gold_two_chain_fork.py`. The helper scripts are in the same experiment directory.
