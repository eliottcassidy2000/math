# Positive-density obstruction for the gap-one family

**PROVED in the stated scope; independently audited in `audit/REPORT.md`.**
Derived in coordination with the one-return coverage investigation. This is
a separate audit from that of `../canonical-deficit-cores/`.

Let t range over odd integers at least 3, a,k over nonnegative integers, and

    R=2^(a+1), M=2^(k+1)-1,
    A=M(R-1)+2, B=MR+1.

The family is `n=A g(t)+B t-1`. The physical unit-deletion source at odd cover
t has E=d=1 and, by the positive-endpoint core bound, `g(t)>=(t+3)/2`.
The weaker bound `g(t)>=t/2+1` is enough. Hence

    n >= (A/2+B)t + A-1.

Since A>=3, any output n<=N has `t<=2N/(A+2B)`. The number of odd t>=3
in that interval is at most `N/(A+2B)`; no additive counting error is needed.
Different t or exponent choices can collide, which only decreases the number
of distinct targets. Thus for every N>=1,

    |family intersect [1,N]| <= c N,
    c=sum_(a,k>=0) 1/[(2^(k+1)-1)(3*2^(a+1)-1)+4].       (1)

The sum converges geometrically. Here is a short exact upper bound. The
four terms a,k<2 sum to

    S=1/9+1/15+1/19+1/37.

For a>=2, all k, discard +4 in the denominator and separate the sums:

    sum_(k>=0)1/(2^(k+1)-1) <= 1+2/3=5/3,
    sum_(a>=2)1/(3*2^(a+1)-1) <= 2/23.

The corresponding tail is at most 10/69. For the remaining tail k>=2 and
a=0,1, the same argument gives

    sum_(k>=2)1/(2^(k+1)-1) <= 2/7,
    tail <= (2/7)(1/5+1/11)=32/385.

The elementary geometric comparisons used here are
`2^(k+1)-1>=3*2^(k-1)` for k>=1,
`3*2^(a+1)-1>=23*2^(a-2)` for a>=2, and
`2^(k+1)-1>=7*2^(k-2)` for k>=2.
Consequently

    c <= S+10/69+32/385
       =27199346/56025585 < 17/35 < 1/2,

with exact positive gap `17/35-upper=13081/56025585`.

It follows that this *entire infinite-core family* misses at least
`N/70-O(1)` odd targets through N, or asymptotically at least 1/35 of all
odd targets. This is stronger than ruling out a fixed finite core bank:
t, a, and k were all free. It does not show that the union with other physical
core families is sparse, and does not obstruct the unrestricted C1/C2/C3
spectrum.

The verifier also checks a separate coarse geometric bound from a 6-by-6
initial square. The omitted tails satisfy denominator >=5*2^(a+k), giving
a combined bound 1/40. That alternative exact-rational computation yields
an upper coefficient below 0.471, but the four-term rational proof above is
sufficient and is the displayed conclusion.
