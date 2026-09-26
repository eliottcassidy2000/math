---
id: THM-4493
title: "Effective gcd and source-height certificates for chronological Collatz rank"
status: "PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT controls"
source: "crossroads233-20260926 carry lane; root and geometry proof audits"
depends_on: []
proofs:
  - 05-knowledge/results/crossroads233_20260926_carry.md
scripts:
  - 04-computation/experiments/crossroads233_20260926_carry.py
---

# THM-4493: height, gcd overlap, and chronological rank

**PROVED + independently audited.** Reserved in checkpoint bb630ddb8,
then promoted after root and geometry audits. The argument is elementary.
It gives an explicit sufficient rank certificate and a descent-or-rank
alternative. It does not prove Collatz or independence of every injective
chronological prefix (G2).

## The valuation certificate

If integers a_0,...,a_(N-1)>1 are multiplicatively dependent, then some i
satisfies

    a_i divides product_(j!=i) gcd(a_i,a_j).                    (1)

Choose a nonzero relation with coefficient vector c and a supported index i
of maximal |c_i|. For every prime p its valuation equation gives

    v_p(a_i) <= sum_(j!=i) v_p(a_j).

Truncating each term on the right at v_p(a_i) preserves the inequality:
either one term reaches that value or all terms remain unchanged. This
proves (1). Consequently strict numerical dominance over the gcd product
at EVERY index suffices for independence. The conclusion singles out one
maximal-coefficient index, not every supported index.

For an actual positive odd Collatz orbit, write

    m_0=n, m_(j+1)=(3m_j+1)/2^k_(j+1), S_j=sum_(r<=j) k_r,
    C_(i,j)=sum_(r=i)^(j-1) 3^(j-1-r)2^(S_r-S_i), i<j.

The exact identity

    2^(S_j-S_i)m_j=3^(j-i)m_i+C_(i,j)

implies gcd(m_i,m_j) divides C_(i,j). Put
Q_i=product_(j!=i) C_(min(i,j),max(i,j)). Thus

    m_i>3^(N-2)Q_i for every i                                (2)

implies independence of 3m_0,...,3m_(N-1), and hence of the pairs
(-3m_i,3m_i+1). For a fixed valuation word, the explicit sufficient cutoff

    n > max_i 3^(N-2) Q_i 2^S_i/3^i                          (3)

follows from m_i>=3^i n/2^S_i. This supplies elementary eventual
independence on each fixed word, without an S-unit theorem.

## Explicit pointwise no-dip theorem

Let N>=2 and define the exact rational bound

    H_N = 2^(N-1)(N-1)!3^(N(N-2))/2^((N-1)(N-2)/2).

If n>=N, n>H_N, and the first N odd nodes m_0,...,m_(N-1) are all
at least n, then their N first slots are multiplicatively independent.
Here H_2=2, H_3=108, H_4=39366, and H_5=86093442.

**Proof.** The product identity and m_i>=n give, for ell<N,

    2^S_ell <= 3^ell product_(i<ell)(1+1/(3m_i))
             <=3^ell(1+1/(3n))^ell <2*3^ell.

Hence C_(i,j)<2(j-i)3^(j-1)/2^i. Multiplying the incident bounds yields

    Q_i < R_i
      =2^(N-1)i!(N-1-i)!3^[(N-1)(N-2)/2+i(i-1)/2]
       /2^[i(N-1)-i(i+1)/2].

Its successive ratio is

    R_(i+1)/R_i = (i+1)/(N-1-i) * 3^i * 2^(i-N+2).

These ratios increase, so the maximum occurs at an endpoint. The last
endpoint is largest because R_(N-1)/R_0=(3/2)^((N-1)(N-2)/2).
Thus 3^(N-2)Q_i<H_N<n<=m_i for every i, establishing (2).

Stirling's formula gives

    log_2 H_N=(log_2 3-1/2)N^2+N log_2 N+O(N).

Therefore for each fixed c<1/sqrt(log_2 3-1/2)=0.960047311978290...,
every sufficiently large no-dip source has full rank through
N=floor(c sqrt(log_2 n)) odd nodes. The endpoint constant is not claimed.

If no descent is known only through floor(log_2 n) shortcut steps, put
ell=ceil(N log_2 3)+1 and additionally require n>=max(N,2^ell).
Were S_N>ell, the first ell steps would have at most N odd actions and
all intermediate values at least n. Thus

    T^ell(n)/n <= (3^N/2^ell)(1+1/(3n))^N
                <=exp(1/3)/2<1,

contradicting the horizon assumption. This puts the needed odd nodes
within that known no-descent horizon.

## Orbit-level consequences and boundaries

For every odd n>max(N,H_N), either one of the first N odd nodes is below
n, or their N first slots are independent. This is a finite descent-or-rank
alternative, with no bound imposed on the valuation word.

Any infinite injective positive Collatz orbit has arbitrarily long
consecutive full-rank odd blocks and infinite rank of its whole first-slot
family. Indeed its odd nodes tend to infinity, since an injective sequence
of positive integers eventually leaves every finite set. Future tail
minima also tend to infinity and are attained. For any N, choose such a
minimum above max(N,H_N); the next N nodes satisfy the theorem. This does
not imply that every prefix is independent and does not exclude divergence.

The [full proof and controls](../../05-knowledge/results/crossroads233_20260926_carry.md)
also establishes:

- The depth-233, 153-odd-step collision lifts to arbitrarily high sources.
  Both colliding trajectories have 153 independent first slots, certified
  by private prime factors without factoring the large nodes.
- For positive coprime integers A,B and digits D distinct modulo A,
  concatenating blocks (Ax+C_0+Ad)/B gives at most |D| integer sources
  at a fixed endpoint: all digits except the first are determined.
  For the two 233-step blocks, capacity stays two at every block depth.
- (6,10,15) is independent without strict gcd dominance or private primes;
  (2,3,6) has one nontrivial relation; the positive-definite log-gcd
  kernel of (2,4) does not imply multiplicative independence.
- Arbitrarily large (4^a-1)/3 descend immediately to 1 and then repeat.
  Height alone, with variable words and no no-dip condition, is insufficient.

Theorem [THM-4490](THM-4490-collatz-affine-word-sunit-specialization.md)
has a different advantage: a longer growing prefix for almost all starts,
uniformly over intervals. The present result is pointwise on high no-dip
sources. No publication-priority claim is made for eventual independence
under translation; the full note links the relevant primary literature.
