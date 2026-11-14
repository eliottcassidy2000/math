# Two independently movable dyadic wings, with a different physical core

For all independently chosen positive integers s,t, put u=2^s, v=2^t.
There is a strong tournament for each of

    336uv + 141u + 228v + 101,
    336uv + 141u + 228v + 103.                           (1)

The two cores have different path counts (9 and 33), different
start and end pair counts, and different orders. This is not obtained
by using the same fixed full profile at two offsets.

Use [the bilateral readout](../chanpak_movable_exact_windows_20261114/MOVABLE_TRIPLES.md)
with h=core path count, f=ordered 2-path-cover count allowing empty,
P=one path starts in the contact, Q=one path ends outside the contact,
K=count after adding the apex to the core. Its two transitive wing
addresses x,y contribute fxy+Qx+Py+h(g(x)+g(y)-2)+K plus the two
endpoint parity terms. When x=u and y=v are powers of two at least 2,
the latter terms vanish, and g(u)=u+1, g(v)=v+1. Thus the coefficients
are exactly (f,Q+h,P+h,K). This triangular projection allows h to
change while preserving all three leading coefficients.

    core outgoing rows        contact  h    f    P    Q    K
    (0,5,9,19,39,15,63)           81    9  336  219  132  101
    (8,5,33,22,7,27)              45   33  336  195  108  103

Contact bit i means the apex beats old vertex i. The left wing goes
forward with contact word 1 then s zero bits; the right goes forward
with t one bits and a terminal zero bit. The apex reaches the first
vertex and is reachable from the last, so the whole tournament is
strong. The first pair is 2183,2185. No ratio constraint on s,t is used.

## Exact inverse and a location bound

For N and each d=101,103, enumerate powers u>=2 with 141u<N. Compute

    v=(N-d-141u)/(336u+228).

Accept precisely if the numerator is divisible by the denominator
and v is a power of two >=2. This uses O(log N) arithmetic operations
on O(log N)-bit integers; a witness needs O(log N) vertices. An
equivalent factor certificate is
(336u+228)(336v+141)=336(N-d)+228*141, with both quotients dyadic.
This is an exact absolute selector within this chart, not a rule
asserting existence of factors for arbitrary N.

There are Theta((log X)^2) distinct realizable counts <=X in (1).
The upper bound is the number of parameter pairs. For either d,
restrict s<=t; then v_2(N-d)=s, since 141 is odd, 228 is divisible
by 4, and the product coefficient is even. The valuation recovers
s and the linear quotient recovers t. Every pair with s<=t and
s+t<=floor(log_2(X))-10 contributes for sufficiently large X,
giving the lower bound. The added axis increases the number of
centers but does not give a covering family.

Literal proof replay: python3 two_axis_replay.py verifies both core
profiles and eighteen enlarged tournaments by subset DP. Running
python3 two_axis_replay.py 2185 emits actual adjacency rows.
