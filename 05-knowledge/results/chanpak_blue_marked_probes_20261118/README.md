# What the blue calculator can read beyond the total Hamiltonian count

Continuation of the [self-converse core packet](https://github.com/eliottcassidy2000/math/tree/84409f9a3d631d94eacb17b33167307a9c95bc68/05-knowledge/results/chanpak_self_converse_core_20261118).
The triangular extraction and the order-six unmarked collision came
from the follow-up research pass. This note independently replays
them and gives two literal marked probes that separate the collision.

## Recover the entire slot polynomial with linear-order blue probes

For a nonempty source X of order n let q_j be its number of ordered
j-slot covers with every slot nonempty. Put
P_X(z)=sum_j binom(z,j) q_j, its ordered-cover polynomial allowing
empty slots. The [Grinberg–Stanley tournament expansion](https://arxiv.org/html/2307.05569)
specialized to z variables equal to 1 identifies it as a monic
polynomial with possible degrees n,n-2,... and no constant:
a packing of directed odd cycles contributes 2 per nontrivial
cycle and one power of z per cycle or uncovered singleton.
The empty packing gives z^n, and each odd cycle of length l
reduces the exponent by the even number l-1.

Consequently for d=floor(n/2), d evaluations P_X(1),...,P_X(d)
suffice: subtract z^n, divide by z when n is odd or z² when n is even,
and interpolate a polynomial of degree <=d-1 in z². It thus
suffices to extract q_1,...,q_d, since the conversion is triangular.
The earlier two strong blue probes give q_1=H(X) and q_2=P_X(2)-2q_1.

For m>=2 suppose q_1,...,q_m are known. Compute the Hamilton count C_m
of the strong self-converse cyclic blowup C3[X,X^op,TT_m]. Its
Hamilton order contracts to a cyclic succession of nonempty runs
in the three blocks. A transitive center has exactly m! ordered
covers into m nonempty slots. The only run patterns using q_(m+1)
on either arm have counts (m+1,m,m),(m,m+1,m),(m+1,m+1,m), each
with a unique starting module. All other terms use only known
q's (the center may be split into at most m runs). Let R_m be the
exact count formula using q_(m+1)=0. Then

    C_m - R_m = m! (q_(m+1)^2 + 2 q_m q_(m+1)).

The positive branch q_(m+1) = sqrt(q_m²+(C_m-R_m)/m!)−q_m
is exact. For n>=4, the two bootstrap probes and d-2 cyclic probes
use d Hamilton reads, all at order <=max(2n+3,2n+d-1). The small
orders n<=3 can use the two bootstrap probes directly. The
conditioned *unmarked* profile is the whole information that goes
through a cyclic center, no matter how the center is changed.

## Port coordinates that unmarked centers lose

For a marked vertex r of X, let e_r be the number of Hamilton
paths of X ending there. Let q_r^start denote the number of ordered
two-slot covers whose second slot starts at r. The first slot may
be empty, the second may not. Their sum over r is P_X(2)-H(X),
so the unmarked polynomial only sees their sum at this arity.

Take X -> TT_s -> X^op and reverse the corresponding arc r -> r^op.
This is self-converse under the same exchange of arms. If X is
strong of order >=3, the wrapper is strong, since the return
r^op -> r is reachable from every vertex of the right arm, and r
reaches every vertex of the left arm. The remaining cross-arcs
and optional center complete the reachability argument. At s>=1
an arbitrary path in X, the center and its reflected path furnish
an explicit blue frame. At s=0 use a path of X ending at a vertex
different from r, available from a directed Hamilton cycle of X.

If a Hamilton path uses the reversed arc, the two segments are
ordinary ordinal runs respecting the old block order. Restricting
to X gives two ordered slots with the second starting at r. On
X^op reversal gives the same read on the two slots. The middle
transitive block assigns each vertex freely to either slot, giving
2^s. If the path does NOT use the return, the original ordinal count
h² survives unless the middle is empty: there the loss of the
original forward arc removes exactly e_r² concatenations. Thus

    H_0(r)=h²+(q_r^start)²-e_r²,
    H_1(r)=h²+2(q_r^start)².

For general positive s the second coefficient 2 becomes 2^s.
Two same-source blue reads per root thus return both marked
coordinates by positive square roots. A useful negative boundary:
this simple one-return edit alone cannot give a +2 repair from a
strong source. At s=0 a difference of squares cannot be 2, and at
s>=1, q_r^start>=2: its second slot may be {r} or {r,v} for any
outneighbor v, and the first slot has a Hamilton path by Redei.

## Literal separation at order six

Use outgoing masks

    A=(58,4,1,6,14,30)
    B=(38,8,2,5,15,30).

Both have P_X(z)=z^6+8z^4+8z² and h=17. The last coefficient combines
four directed 5-cycles versus two plus one packing of two disjoint
triangles. The constant center and every larger cyclic center
give identical reads. In contrast, the multisets of marked start
coordinates are distinct:

    A: 31,18,17,28,44,69
    B: 33,17,20,21,61,55.

Even the unlabelled collections of H_1(r) distinguish the sources;
for example their maxima are 9811 and 7731. These are literal counts
in strong order-13 self-converse tournaments. This is not complete
tournament reconstruction: it is an explicit marked coordinate
that recovers information every unmarked cyclic center misses.
The next useful constructive problem is to combine two marked
contacts so that their square responses cancel to a small signed
difference at a prescribed absolute target. Chan–Pak remains open.

Run `python verify.py` for the full subset-DP witness checks,
triangle and 5-cycle enumeration, cyclic coefficient tests and
polynomial interpolation. Normal and optimized transcripts agree.
