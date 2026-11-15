# A literal three-move crossing, and why a fixed mirror cannot prolong it

## A unit route from 721 to 727

Every row below is an actual strong tournament with one apex and an
11-site transitive spine. A set bit means the apex beats that site.
Starting from the prior bottom witness 1158, follow this route to the
prior top witness 1146:

| Mark word as binary integer | H | XOR with next | Number of marks reversed |
|---:|---:|---:|---:|
| 1158 | 721 | 80 | 2 |
| 1238 | 723 | 140 | 3 |
| 1114 | 725 | 32 | 1 |
| 1146 | 727 | | |

This is a monotone physical route of increments +2, in the same
residue modulo 2. The step is an atomic bundle of one, two or
three legitimate apex-arc reversals; there is no claim about counts
while a bundle is partially reversed. The second word, 1238, is not
a corner of the earlier additive rectangle. This is how the route
bypasses two independent controls without introducing more branches.
Each row starts with bit 1 and ends in 0, so is strong via its literal
spanning cycle. Direct permutation-subset DP verifies the counts.

There is a scoped minimality claim: on the 11-site apex-word cube,
allow arbitrary moves changing at most two marks and having |ΔH|≤2.
The component of 1158 is exactly {1158, 1190, 1270, 1238}, so cannot
contain 1146. Replay is exhaustive within that finite, precisely
declared architecture, not a universal lower bound for all routes.
The new route reaches both interior odd targets; it is not a
construction of endpoint witnesses for every odd N.

## A fixed-width firewall even with arbitrary interior moves

For an m-bit core w with ordinary address u and reverse-complement
address v, let σ(w)=u+v. In the fixed-width mirror family use
exterior addresses a=b=t (padding representations does not matter).
By the exact two-cut seam,

    H_m(t,w) = 2^m t² + t σ(w) + 2 H(t) + H(w) − 2 + e(t,w),

where 0≤e≤2 is the sum of the two boundary parity products.
The four words in the displayed physical route have σ respectively
2428, 2428, 2440, 2440. Thus mirror extension preserves its first
and last steps, but the middle step becomes 2+12t. This is an
exterior response proof, not a numerical experiment.

In fact for fixed m even **arbitrary** interior replacements and
changes of t do not yield a bounded-step route at sufficiently
large positive t. For fixed t, different σ give differences of
magnitude at least t−(2^m−1)−2, since 1≤H(w)≤2^m and e is bounded.
Within each layer t, a move of size ≤M must preserve σ once
t>2^m+1+M.

Also adjacent layers have disjoint absolute bands. Both u and v
are between 0 and 2^m−1, but cannot simultaneously attain the
same endpoint. Thus σ_min≥1 and σ_max−σ_min≤2^(m+1)−4.
Using H(t)≤t+1 and H(t+1)≥1 gives

    min_w H_m(t+1,w) − max_w H_m(t,w)
      ≥ 2^m(2t+1) + σ_min − (σ_max−σ_min)t
        −2t −(2^m−1)−2
      ≥ 2t.

For any fixed allowed step M, large layers therefore cannot be
connected by such steps even if the exterior is changed arbitrarily
and all fixed-width core words are made available. Along with
the different-σ obstruction this explains exactly what the mirror
preserves and what it cannot transport. It does not rule out
growing core width, asymmetric exteriors, or a change of reader.
