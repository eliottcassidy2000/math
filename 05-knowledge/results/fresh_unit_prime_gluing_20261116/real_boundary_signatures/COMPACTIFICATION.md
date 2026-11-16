# When a ready continuation compresses to finitely many dyadic phases

The proof below gives an all-order obstruction to selecting a prescribed absolute
Hamilton count with a normalized root, small deletion count, and bounded
relative defect. It includes both strong and decomposable tournaments.
It neither excludes a large deletion budget nor resolves Chan–Pak.

For a rooted tournament write h=H(T), g=H(T−r), and E,O for the weighted
even and odd rooted-path sums (weight H of the omitted subtournament).
Call it normalized when O−E=g. Define the low intercept and defect as
```
I=h−2O+2g,    Δ=I+g,
```
so Δ=h−2E+g on normalized roots, and Δ=0 is the ready equation used by
the synchronized splice. Logarithms below are natural unless marked.

**Theorem.** Fix C≥0 and G≥1. All counts of normalized roots with
`g≤G, |Δ|≤C g` can be enumerated as a finite exceptional
set and dyadic rays (I_0+Z_0 2^m, m≥1), with the original specified
(g,Delta) held fixed along each ray. They admit at most
```
M_C(G) ≤ exp( κ (ln G)^2 + O_C(ln(G+2)) ),
κ = 9 ln(2) / (2 ln(5)^2),
```
core presentations. Thus the number of distinct such counts at most X is
`O_C(M_C(min(G,X)) log(X+2) + min(G,X))`; the separate
sink term can also be incorporated into M. In particular, for any
```
c < sqrt(2) ln(5) / (3 sqrt(ln(2))) ≈ 0.911,
```
the number of odd H≤X admitting such a root with g≤exp(c sqrt(ln X))
is o(X). The same is true if the condition is (g≤exp(c sqrt(ln H))).

## The exact discrete cutoff

Take the ordered SCCs of the deletion D=T−r and let Q be the root's set
of outgoing contacts. For a path in D its weight is H of its complement
in D. Let K sum the weights of all even-order paths whose first and last
vertices both lie in Q. The odd-cycle identity applied at the root gives
h−g = twice the corresponding sum whose first endpoint lies in Q and
last endpoint outside Q. The odd weighted rooted sum O has the singleton
term g plus the entire even-order sum starting in Q. Therefore
```
I=g−2K,    Δ=2(g−K).
```
This calculation holds even before imposing normalization. It gives
(0≤K≤(1+C/2)G). When Q is empty h=g and the root is a sink; count
those exceptional scalar values directly.

For Q nonempty take the deletion component containing the first contact
and the component containing the last contact. Everything strictly before
the first contact component is an ordinal prefix of the whole rooted
tournament. Its singleton factors can be discarded, since their path
contribution is 1 and a path starting at the root never enters this prefix.
If the endpoint components are different and j≥1 singleton SCCs occur
strictly between them, use any even subset of those j singleton vertices
as intermediate vertices between a chosen contact in the first and one
in the last component. This yields (2^{j−1}) distinct even paths with
weight at least 1. Consequently, over the active component interval,
there are at most (3+\lfloor\log_2 K\rfloor) singleton SCCs for K>0;
for K=0 there are at most two endpoint singletons, since if there were
any intervening singleton then the two endpoint contacts alone would
make K>0. The same bound is immediate if the endpoint components coincide.

Finally the suffix of deletion components strictly after the last contact
is entirely incontact, and every remaining nonroot vertex beats it.
Thus it is exactly an opaque near-sink block at the root. In this
operation its effects depend only on
((h_U,q_U)=(H(U),P_U(2))), where (P(2)) counts ordered two-slot covers
with possible empty slots. Namely it reads (h_U I+q_U Z) and lands with
(E=q_U Z/2, O=E+h_Ug, g'=h_Ug), where (Z=E+O−g) before appending.
This is the archived rooted block identity (the earlier variable-reader block transfer).
Both h and q multiply along an ordinal sum: each one- or two-slot path
splits in component order. Consequently rearranging the suffix's singleton
components to the far right does not change a single output profile
coordinate. This is a re-encoding of the boundary, not a relabeling of arcs
in general. Delete those far-right singleton components to obtain a
bounded core, keeping their number m as its address. For m≥1 the
singleton near-sink laws give ((H,E,O,g)=(I_0+2^mZ_0,2^{m−1}Z_0,
2^{m−1}Z_0+g_0,g_0)). The m=0 unpadded core output is read separately.
The core may itself be unnormalized; this causes no difficulty because
one near-sink append automatically normalizes.

## Capacity and its subexponential deletion threshold

Every nonsingleton deletion SCC has intrinsic Hamilton count g_i≥3.
Their product is g. Busch's bound gives (ln g_i≥(n_i−1)\ln(5)/3);
hence all nonsingleton vertices together have total bounded by
((3/\ln5+1/\ln3)\ln G). Combined with the active singleton bound,
the reduced rooted core has order at most
((3/\ln5+1/\ln3+1/\ln2)\ln G+O_C(1)).
More importantly, all arcs between deletion SCCs are prescribed by
their order. Writing a=3/ln 5, the number of binary internal arc slots is
```
Σ_i binom(n_i,2) ≤ (a²/2) Σ_i (ln g_i)^2 +(a/2) Σ_i ln g_i
                  ≤ (a²/2) (ln G)^2 + (a/2) ln G.
```
This upper bound applies to every ordered SCC size composition actually
possible under the cutoff. There are at most exponential in the core
order many ordered size compositions and root contact subsets. One can
label consecutively within each component and freely overcount internal
orientations, so no factorial for rearranging the labels is needed.
The displayed exponent in the theorem follows.

On a nonempty-contact ray (Z_0≥1), and
(|I_0|=|\Delta-g|≤(C+1)G). Since g≤H by induced monotonicity, replace
G by min(G,X) when counting outputs below X. There are (O_C(\log(X+2)))
addresses per core, and at most one count at m=0. This proves the bound.
Putting G=exp(c sqrt(ln X)) gives at most
(X^{\kappa c^2+o(1)}) counts and hence density zero whenever c is below
the stated threshold (with c≥0). An effective target-first decoder for
the capped layer is simply to enumerate this finite pool and test
whether ((N-I_0)/Z_0) is a power of 2 with exponent m≥1, together with
the direct core possibilities. An exact profile requirement for g or
Δ can be checked as well. It is an exact all-height decoder inside
the cap, not a proof that some core works.

## What projective compactness forgets

For fixed G the reduced family has only finitely many logarithmic
phases (\log_2 Z_0\bmod 1). Along all ready rays the normalized state
((H,E,g)/H) tends to the same point ((1,1/2,0)), but the underlying
absolute counts cannot interpolate between these finite phases. In
every sufficiently high octave [X,2X], each of the M core rays contributes
only a bounded number of points, since its shift is bounded and the
address doubles. Thus some interval of length
(Omega_C(X/(M_C(G)+1))) in each octave has no such reader, provided
G=o(X); for fixed G this is a missing interval of positive relative
length. The hypothesis G=o(X) also removes the sink exceptions from
the octave. A continuum/projective limit can close and still leave
large absolute holes. The missing coordinate is the discrete phase
of the stripped core, not the normalized limiting location.

Run **python3 check_compactification.py** for independent literal checks of the normal form, including a nontrivial suffix component. The capacity bound itself is proved above rather than inferred from the finite checks.
