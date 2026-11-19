# A three-parameter +2 family with cubic logarithmic count growth

**Proved:** the following literal two-source construction produces
exact neighboring odd counts at all three lengths. It gives
Theta((log X)^3) distinct odd targets below X with a witnessed
neighbor, not a positive density theorem.

## Physical sources

Use outgoing bitmasks
A=[14,92,8,48,101,71,13], B=[88,17,91,98,40,71,18].
On marks 0,1,2 in that order substitute TT_m, TT_n, TT_k.
The marks induce a transitive tournament in each source; hence
no run between consecutive outside vertices uses a mark twice.
With four outsiders, each mark has at most five visits.
The inherited ordered-cover substitution rule gives

    H = sum c_rst T(m,r) T(n,s) T(k,t),

where T(m,r)=r!S(m,r). The full finite compressed-word coefficient
tables agree except at (1,1,1), which is 89 versus 91. This proves
a difference of two at every choice of positive m,n,k. The
standalone replay independently enumerates the coefficient tables
and checks the literal 3×3×3 grid in both variants. For substitutions
of general nontransitive modules the difference scales with the
product of their one-path counts and need not be two.

The exact counts are C±1 with C=8+V+2^k W and

    V = -7·2^n + 2·3^n - 9·2^m + 11·2^(m+n)
        - 4·2^m·3^n + 2·3^m - 4·3^m·2^n + 2·3^(m+n),
    W = 24 - 28(2^m+2^n) + 8(3^m+3^n) + 66·2^(m+n)
        - 40(2^m·3^n+3^m·2^n) + 44·3^(m+n)
        + 6(2^m·4^n+4^m·2^n) - 12(3^m·4^n+4^m·3^n)
        + 6·4^(m+n).

This is the expansion of the exact coefficient table by the
inclusion–exclusion formula for T; the replay also checks its
integer exponential coefficients. From the base89/91, increasing
only one length from 1 to 2 gives 267/269, 271/273, or 177/179.

## The three lengths really produce cubically many counts

Set S=2(m+n)+k. Divide C by 6·2^S. Its dominant term is
W/(6·4^(m+n)), which equals

    alpha_m + beta_m (3/4)^n + gamma_m (1/2)^n + delta_m (1/4)^n,

where alpha_m=1−2(3/4)^m+(1/2)^m, beta_m tends to −2,
and gamma_m,delta_m are bounded independently of m.
The normalized contribution from 8+V is O((3/4)^(m+n)),
uniformly for k>=1. These facts follow term by term from the formula.

Restrict to m>=M, n>=m+L for sufficiently large fixed constants
M,L. Then C/(6·2^S) lies in [0.9,1.1], so disjoint dyadic bands
recover S exactly. Its deficit from 1 equals

    2(3/4)^m (1+O((2/3)^m + (3/4)^(n−m))).

Choose M,L so that the relative error is at most 0.1. The
resulting geometric bands are disjoint (1.1·3/4 < 0.9), so this
recovers m. Given m, its deficit from alpha_m similarly equals

    −beta_m (3/4)^n (1+O((2/3)^n + (3/4)^m)),

again with disjoint relative bands after increasing M. This
recovers n, then k=S−2(m+n). Thus C is injective on this wedge.
Taking all three lengths below a sufficiently small multiple of
log X produces Omega((log X)^3) eligible triples below X.
For the upper bound, use the positive c_211, c_121 and c_112 terms:
each individual length has a contribution at least 2^length−2.
Thus only O(log X) lengths per coordinate are possible. No
count-density or universal neighbor conclusion follows.

Run python3 triple_port_neighbors_replay.py.
