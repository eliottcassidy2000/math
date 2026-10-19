# Height-controlled entry into a prescribed Collatz bank word

19 October 2026. Status: PROVED, as an elementary combination of classical parity encoding, large deviations and affine carry control. This strengthens the earlier bank-entry calculation by controlling the height relative to the original source and the available time. It is not a new generic Collatz density theorem.

## 1. The result and its boundary

Let

    T(n)=n/2 if n is even, and T(n)=(3n+1)/2 if n is odd.

Use the explicit 32-bit bank word W=KF from the earlier note:

    W=11011011011011010101101110101010.

Its complete lawful source cylinder is

    C_W={y>0: y=2962948091 (mod 2^32)},

and on it

    T^32(y)=(3486784401 y+15810005461)/4294967296 < y.    (1)

The strict inequality holds for every positive lift: the affine fixed point is less than 20, whereas the least positive source is 2962948091.

**Concrete theorem.** Apart from O(X^(1-delta)) positive integers x<=X, there is a time j for which

    T^j(x) belongs to C_W,
    T^(j+32)(x) < T^j(x) <= x^(15/16),
    j+32 <= (7/8) log_2 x,                            (2)

where the explicit positive saving is

    delta= -log_2(1-2^(-32))/256
           = 0.000000000001312123496... .              (3)

The implied constant is effective. Section 4 gives a fully explicit dyadic exceptional-count bound, rather than concealing it in O notation.

Thus this particular bank entry has occurred below a power of the original source, and its certified endpoint is also below the original source. The earlier theorem only gave eventual entry and a drop relative to the entry height. The extra conclusion here comes from simultaneous height and time estimates, not from changing the quantifier “density one” to “every integer.”

The saving in (3) is extremely small because one specified 32-bit word is rare. The theorem is an asymptotic statement, not a useful numerical waiting-time guarantee for an ordinary-sized individual source. No positive integer is assumed to have infinitely many occurrences of W.

**General exponent.** For every fixed nonempty binary word w and every

    log_4 3 < beta <1,

all but O(X^(1-epsilon)) sources x<=X enter w's lawful cylinder at a height at most x^beta within log_2 x Terras steps, for some effective epsilon>0 depending on w and beta. The clock can in fact be c log_2 x for a fixed c<1 chosen with the parameters. If w has a uniform all-lift descent property, its ending point also lies below x^beta.

The exponent threshold log_4 3 is classical Korec territory. Korec also gives a log_2 y bound on the generic descent time (page 86). The result here is the prescribed-word requirement together with an explicit joint time/height exceptional count; a logarithmic descent clock by itself is not claimed as new. No improvement of Korec's exponent, no strengthening of Tao's almost-bounded logarithmic-density theorem, and no general priority claim are intended.

## 2. Classical finite-word facts, with the carry retained

For each length N, initial residues modulo 2^N are in bijection with the N-bit parity vectors. To see the induction step, two lifts of one source residue modulo 2^N differ by 2^N. After the prescribed N steps their endpoints differ by the odd number 3^h, so their next parities are opposite. This is the standard Terras/Everett encoding.

For a prefix of length j write its exact affine branch as

    T^j(x)=A_j x+B_j,
    A_j=3^(H_j)/2^j, B_j>=0,                       (4)

where H_j is its odd-bit count. Under uniform N-bit parity vectors, H_j is Binomial(j,1/2), for j<=N.

The additive term is essential. Its updates are

    B_(j+1)=B_j/2              after a 0 bit;
    B_(j+1)=(3B_j+1)/2        after a 1 bit.

Averaging over the fresh next bit gives

    E[B_(j+1) | first j bits]=B_j+1/4,
    E[B_j]=j/4.                                      (5)

The last identity is also Lemma 2.5 in Rozier–Terracol's paper cited below. Equations (4)-(5) are being applied to a finite uniform parity space, not to a stochastic model imposed on one integer's infinite orbit.

## 3. A completely specified dyadic experiment

Fix an integer m>=512 and put

    X=2^m,  N=floor(7m/8),  J=ceil(3m/4),
    L=32,   k=floor((N-J)/L),
    p=2^(-32),  beta=15/16,
    q=[10^20/(9^9 11^11)]^(1/20) <1.                (6)

The sources are exactly the X integers in [X,2X). Since N<=m, every parity vector of length N occurs exactly 2^(m-N) times in this interval. There is no residue-rounding error.

Inspect the k disjoint candidate blocks of length L starting at

    j=J, J+L, ..., J+(k-1)L.

Every candidate word finishes by time N. Define three bad events on this finite parity space:

1. W is absent from all k candidate blocks.
2. Some prefix with J<=j<=N has H_j>11j/20.
3. Some prefix with J<=j<=N has B_j>X^beta/2.

We bound these events separately on the original uniform space, then take a union bound. In particular, neither parity independence nor the mean carry formula is asserted after conditioning on entry into W or on any selected stopping time.

The first event has probability exactly

    (1-p)^k.                                          (7)

For the second event, use z=11/9 in the elementary exponential Markov bound:

    P(H_j>11j/20)
      <= E[z^H_j] z^(-11j/20)
       = [(1+z)/2]^j z^(-11j/20)
       = q^j.

Consequently its probability is at most

    q^J/(1-q).                                       (8)

The exact integer inequality 10^20<9^9 11^11 establishes q<1. One may also certify q<999/1000 by raising both sides to the twentieth power.

For the third event, (5), Markov's inequality and a union bound give

    sum_(j=J)^N P(B_j>X^beta/2)
       <= sum_(j=J)^N j/(2X^beta)
       <= N(N+1)/(4X^beta).                          (9)

No independence between the different B_j is used.

## 4. Why a good source enters below its original height

Suppose none of the three bad events occurs, and choose any candidate start j at which W occurs. By H_j<=11j/20 and 3^5<2^8,

    A_j=3^H_j/2^j <=2^(-3j/25)
         <=2^(-9m/100)=X^(-9/100).

As x<2X and B_j<=X^beta/2, equation (4) gives

    T^j(x) < 2X^(91/100)+X^(15/16)/2
             <= X^(15/16) <= x^(15/16).             (10)

The middle inequality is valid for m>=73: its required comparison is

    2 <= (15/16-91/100)m =11m/400.

Our declared m>=512 is more than sufficient and guarantees a nonempty candidate-block set. The actual occurrence of W makes its complete branch lawful, so (1) now gives the strict endpoint inequality in (2). Finally

    j+32<=N<=7m/8<= (7/8)log_2 x.

Thus the number of sources in [2^m,2^(m+1)) failing (2) is at most

    2^m [(1-2^(-32))^k
           +q^J/(1-q)
           +N(N+1)/(4*2^(15m/16))].                (11)

The bound can exceed 2^m at small m, in which case the trivial bound is stronger. Formula (11) is still valid; it is not evidence of a practical success rate at that scale.

To deduce (3), note k>=m/256-2. Therefore the first term in brackets is at most

    (1-2^(-32))^(-2) X^(-delta).

The exponent from (8) is -(3/4)log_2 q>0.0054, much larger than delta. The carry error in (9) has exponent 15/16, with only a quadratic factor in m. Hence (11) is O(X^(1-delta)). Summing dyadic intervals gives the stated bound for all sources up to an arbitrary real X; the bounded initial intervals are absorbed into the constant.

For a purely rational comparison of these exponents, q<999/1000 gives -(3/4)log_2 q>3/4000, while -log_2(1-2^(-32))/256<2^(-38)<3/4000. Thus their ordering does not depend on displayed decimal approximations.

## 5. The general prescribed-word power-height theorem

Fix a word w of length L and log_4 3<beta<1. Choose parameters

    1/2<p0<log_3 2,    0<a<c<1,
    eta=1-a(1-p0 log_2 3)<beta.                    (12)

Such parameters exist because eta tends to log_4 3 as a tends to 1 and p0 decreases to 1/2. Work in [2^m,2^(m+1)), use N=floor(cm), J=ceil(am), and inspect disjoint L-bit blocks between J and N.

The three preceding estimates have the same forms:

- Pattern avoidance is (1-2^(-L))^floor((N-J)/L).
- The probability of some H_j>p0 j for j>=J is O(q0^J), where q0<1 is the exponential Markov value obtained from z=p0/(1-p0).
- The probability of some B_j>2^(beta m)/2 is O(m^2 2^(-beta m)).

Outside these events,

    T^j(x) <=2*2^(eta m)+2^(beta m)/2 <=2^(beta m)

for all sufficiently large m. All exponents in the bad-event estimates are strictly positive. Dyadic summation yields the stated effective power-saving exceptional bound. This proves the general theorem without asserting that two unrelated density-one sets may be conditioned on one another.

A finite list of prescribed words, or any fixed finite occurrence quota, can be incorporated with the corresponding finite union or binomial-tail bound. Nothing permits an infinite quota for a fixed integer.

## 6. An optional bound on the whole pre-entry excursion

The same finite-prefix argument also controls growth before entry if desired. The multipliers A_j form a nonnegative mean-one martingale on the finite parity space, because the two possible multiplier factors are 1/2 and 3/2. The elementary first-hit proof of the maximal inequality gives

    P(max_(j<=N) A_j>X^epsilon)<=X^(-epsilon).

Also (5) and a union bound give

    P(max_(j<=N) B_j>X^(1+epsilon))
       <=N(N+1)/(8X^(1+epsilon)).

Outside these two extra events, (4) gives

    max_(j<=N) T^j(x) <=3X^(1+epsilon)<=3x^(1+epsilon).

For example epsilon=1/100 adds a 3x^(101/100) ceiling for the entire inspected trajectory, with the same final O(X^(1-delta)) exception count. This does not assert that the pre-entry trajectory stays below x; an odd source already violates that at its first step.

## 7. Scope, inheritance, and primary sources

The proof does not use THM-4581's almost-sure affine coalescence or THM-4590's basin equidistribution. Their density statements do not supply pointwise control of an integer in an exceptional set. The present proof also leaves arbitrary pointwise entry, recurrence, and Collatz convergence open.

It complements the deterministic bank theorem in collatz_original_height_repair_20261019/ORIGINAL_HEIGHT_REPAIR.md. That theorem handles arbitrarily long contracting bank words once they start at the source, using a 92-step repair or a smaller periodic obstruction. This note handles an explicit statistical set of entries below the original source's height. The statements have different quantifiers and neither turns the other's exceptional set into an empty set.

The basic encoding is due to Riho Terras, “A stopping time problem on the positive integers,” Acta Arithmetica 30 (1976), 241–252; and C. J. Everett, “Iteration of the number theoretic function f(2n)=n, f(2n+1)=3n+2,” Advances in Mathematics 25 (1977), 42–45. The finite encoding proof is included here.

Ivan Korec, “A density estimate for the 3x+1 problem,” Mathematica Slovaca 44 (1994), 85–89, Theorem 1 proves generic power descent for every exponent above log_4 3: https://dml.cz/handle/10338.dmlcz/133225 . That classical threshold is not claimed as new here.

Olivier Rozier and Claude Terracol, “Paradoxical behavior in Collatz sequences,” Discrete Mathematics 349 (2026), 115167, arXiv:2502.00948v5, Lemma 2.5 records E[B_j]=j/4: https://arxiv.org/html/2502.00948v5 . The proof is reproduced in section 2 and needs no computation.

Terence Tao, “Almost all orbits of the Collatz map attain almost bounded values,” https://arxiv.org/abs/1909.03562 , concerns a stronger height conclusion in logarithmic density. No comparison here changes that distinction of density or claims an improvement on that theorem.

The companion checker authenticates W, the rational parameter inequalities, the finite carry identities and small independent parity-count controls. The all-size counting proof is the argument above, not extrapolation from those controls.
