# A000568: growth obstructions and a convergent normalization

Status: PROVED below, with classical provenance; FINITE-EXACT checks are separate.
Source checkout: `5862c94555eafea4e9a0638a17f7c661f8e45d37`.
No discovery-priority, polynomial-time counting, or finite-state impossibility claim.

**Extension:** `ALL_POLES_CORE_EXPANSION.md` proves global meromorphic
continuation with exact poles 4^k of order 3k+1, strengthening Section 2.

## 1. What kind of generating function can exist?

Write N(n)=n(n-1)/2 and a_n for unlabelled tournaments, with a_0=1.
Since every orbit of labelled tournaments has at most n! elements,

    2^N(n)/n! <= a_n <= 2^N(n).

Consequently both Σ a_n z^n and Σ a_n z^n/n! have radius ZERO:
their nth-root lower bounds are respectively 2^((n-1)/2)/(n!)^(1/n)
and 2^((n-1)/2)/(n!)^(2/n), and n!<=n^n suffices.
More generally, division by any fixed real power of n! cannot give positive
radius. These series nonetheless exist perfectly well as FORMAL series.

There is no fixed-order linear recurrence with polynomial coefficients
(over C), valid for all sufficiently large n, for a_n. Indeed, write a
hypothetical recurrence with last nonzero coefficient p_r(n) as
p_r(n)a_(n+r)=-Σ_(j<r)p_j(n)a_(n+j). Beyond finitely many roots of p_r,
all ratios |p_j(n)/p_r(n)| are at most C n^D for some integer D>=0.
The running maximum therefore grows at most by C' n^D per step, giving
|a_n|<=K C''^n (n!)^D. This contradicts the lower bound above, whose
logarithm is (log 2)n²/2-O(n log n). The order-zero case would require
a_n=0 eventually and is also impossible. Thus A is not D-finite; in
particular it is neither rational nor algebraic (characteristic zero).
The same direct bound excludes P-recursiveness of a_n/(n!)^c for fixed c.

This does NOT exclude nonlinear recurrences, unbounded-order recurrences,
recurrences with exponential coefficients, catalytic/infinite-state systems,
or formal plethystic equations. Growth alone does not refute a finite
combinatorial description with n-dependent branching weights.

## 2. Symmetry-normalized counts are bounded and tend to one

Define b_n=n!a_n/2^N(n). Burnside gives the exact identity

    b_n = Σ_(sigma in S_n) Prob[sigma is an automorphism of uniform labelled T]
        = E[|Aut(T)|].

It is an EXPECTED GROUP ORDER, not the probability of asymmetry, and not
a uniform-unlabelled expectation. Only all-odd-cycle permutations contribute.
For such a permutation moving m vertices, no unordered pair is fixed unless
both vertices are fixed: a fixed non-pointwise pair would be a transposition.
There are M=N(n)-N(n-m)=m(2n-m-1)/2 affected unordered pairs. Each affected
pair orbit has odd length >=3, so the loss in independent orientation bits
is at least 2M/3. Thus its fixed probability is at most

    2^[-m(2n-m-1)/3] <= 2^[-m(n-1)/3].

The number of permutations with support m is at most n^m, and nonidentity
contributors have m>=3. Setting q_n=n*2^[-(n-1)/3], for all n with q_n<1,

    0 <= b_n-1 <= Σ_(m=3)^n q_n^m <= q_n³/(1-q_n) -> 0.

Hence (b_n) is bounded, b_n->1, and B(z)=Σ b_n z^n has radius EXACTLY 1.
The bound also implies P(T has nontrivial automorphism)<= (b_n-1)/2,
because tournament automorphism groups have odd order, so a nontrivial
one has order at least three. The weaker bound b_n-1 needs no parity fact.

### First symmetry correction and a genuinely analytic remainder

The support-three permutations are the 2*C(n,3) three-cycles. Each loses
exactly 2n-4 bits. Their contribution is

    L_n = C(n,3)*2^(5-2n).

There are no support-one, support-two, or support-four all-odd contributors.
Let R_n=b_n-1-L_n>=0. For supports 5<=m<=n/4, the preceding loss estimate
is at least m(7n/4-1)/3. Their total is bounded by q'_n^5/(1-q'_n), where
q'_n=n*2^[-7n/12+1/3], once q'_n<1. For m>n/4, the loss is at least
n(n-1)/12 and there are at most n! permutations, giving a bound
n!*2^[-n(n-1)/12]. Therefore

    limsup R_n^(1/n) <= 2^(-35/12) < 1/4,
    b_n-1 ~ C(n,3)*2^(5-2n).

The initially empty ranges cause no issue; this is an eventual bound.
Consequently

    B(z) = 1/(1-z) + z³/[2(1-z/4)^4] + R(z),

where R(z)=ΣR_n z^n is analytic for |z|<2^(35/12). Thus B has a
meromorphic continuation past 4 with a simple pole at 1 and an order-four
pole at 4. This is a modest consequence of classical asymptotic enumeration,
not a new unrestricted enumeration formula. In particular B-1/(1-z) has
radius exactly 4. This elementary argument alone does not classify all later singularities;
the companion all-poles note now does so via the inherited core kernel.

## 3. SCC composition survives as a deformed product

Let s_n count strongly connected unlabelled tournaments, with s_1=1,s_2=0.
Strong components are canonically linearly ordered, and all cross edges are
forced. Thus, formally,

    a_n = Σ_(k=1)^n s_k a_(n-k),   A(z)=1/(1-S(z)).

This identity is classical and gives s_n from a_n; it is not by itself a
closed recurrence determining a_n without additional strong/prime data.
Define d_0=0, d_n=n!s_n/2^N(n). Transforming the recurrence gives

    b_n = Σ_(k=1)^n C(n,k) 2^[-k(n-k)] d_k b_(n-k).

For arbitrary sequences set

    (f star g)_n = Σ_(k=0)^n C(n,k) 2^[-k(n-k)] f_k g_(n-k).

This product is associative and commutative. Its three-factor coefficient
is n!/(i!j!k!)*2^[-ij-ik-jk], independent of bracketing. Its identity is
(1,0,0,...). Equivalently it is ordinary convolution conjugated by
f_n -> 2^N(n)f_n/n!. Hence B=1+D star B and B=(1-D)^(-star).
The displayed inverse is a coefficientwise formal inverse; it is NOT the
ordinary analytic identity B=1/(1-D). Here 0<=d_n<=b_n, so D also converges
at least on the unit disk. No isotype information is created by normalization.

For contrast, labelled tournament probabilities have total sequence 1 and
strong probability rho_n. They obey
1=Σ C(n,k)2^[-k(n-k)]rho_k, since the remaining normalized labelled count
is 1. The same kernel appears, but d_n is not rho_n: d_n includes the
unlabelled symmetry weighting. This distinction is an important hostile.

## 4. Provenance and corrected near miss

* E. M. Wright, *The number of irreducible tournaments*, Glasgow Math. J.
  11 (1970), 97–101, DOI https://doi.org/10.1017/S0017089500000914 .
  Page 97 equation (2) is the unlabelled SCC recurrence; page 98 records
  a_n~2^N/n!; page 100 develops symmetry corrections. Full primary PDF:
  https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C7FD99353A0E91D4433B31D18FC87344/S0017089500000914a.pdf/the-number-of-irreducible-tournaments.pdf
  Retrieved 2026-10-29. Terminology: Wright's irreducible means strongly
  connected, not indecomposable/prime under modular substitution.
* J. W. Moon, *Topics on Tournaments* (1968), cited by Wright for the
  unrestricted asymptotic enumeration; no new priority is asserted here.
* Repo `05-knowledge/results/session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md`
  already supplies Burnside and full prime-cycle-index composition.
* Repo `04-computation/a000568_asymptotic_exact_s26.py` has a maintained
  2026-10-28 correction: fixed relative precision does NOT yield exact integer
  counts. Its old output has wrong float-rounded n=14 and n=15 values.
  `05-knowledge/results/strong_attachment_tail_compiler_20261028/HISTORY.md`
  records the repair lineage. The convergent normalization here does not
  repair or imply that polynomial-time exact-counting proposal.
* `05-knowledge/results/sc_a000568_q_deformed.out` uses SC for another count
  (its n=2 value is 1), NOT the strongly connected count (whose n=2 is 0).

## 5. Reproducibility

Run `python3 session98_growth/limits/verify_limits.py` from the workspace.
This uses exact rational Burnside arithmetic, reconstructs through n=30,
checks the initial OEIS values as external controls, SCC recurrences, the
first symmetry term, all support-loss inequalities through n=30, and
associativity of the twisted kernel through total degree 30. It does not
pretend finite checks prove asymptotic or non-P-recursive statements.

Both checkers use explicit `require`/`RuntimeError` gates (not removable
assertions). Normal and `python3 -O` executions passed with byte-identical
outputs on 2026-10-29.
