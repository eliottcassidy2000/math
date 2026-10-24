# Balanced admission removes the bulk boundary congestion

Status: **PROVED refinement** of `PROOF.md`, with exact finite controls in
`balanced.py`. This retains all accepted states; it is not an interior-only
restriction. The incoming full conditional-law and shared-frame requirements
are unchanged.

The deterministic parent in the first bound is deliberately conservative:
many boundary states can load the same representative. Instead distribute
each boundary state's mass over all legal nearest-band admissions. They all
have the same universal clock bound. This is a genuine source-weighted flow,
not a count of paths without demand weights.

## 1. The overlapping-parent variance lemma

Let C be the complete interior band, with its chain H_0 of gap gamma.
For every outside state x choose a probability K_x(c) on legal interior
admissions; put K_c=delta_c for c in C. Use the ACTUAL visible masses M and
suppose

    h(c) = sum_x M(x) K_x(c) / M(c) <= H.                 (1)

Define a reversible chain on the entire visible space by

    P(x,c)=K_x(c)/(4H),                                  (x outside C)
    P(c,x)=M(x)K_x(c)/(4H M(c)),
    P(c,d)=H_0(c,d)/4,                                  (c!=d in C)

and fill the diagonal. Rows are lazy since (1) bounds admission load.
Exactly as for the deterministic parent, compare f(x) with f(c), now average
over K_x(c), and put `tilde_pi(c)=sum_x pi(x)K_x(c)`. Then

    Var_pi f <= 2 sum_{x,c} pi(x)K_x(c)(f(x)-f(c))^2
                  + 2 Var_tilde_pi(f restricted to C)
             <= (8H/gamma) E_P(f,f).

Thus **gap(P)>=gamma/(8H)**. Full conditional reservoir refresh preserves
this bound. Overlapping admission fibres are allowed: they are edges of a
coupling, not disjoint charts incorrectly treated as a partition.

## 2. An explicit balanced choice and its source-weighted load

Let reservoir size S, complete band a<b, and total ordinary occupancy k.
The base case is(S,a,b)=(20,5,15), the pinned cases(20,7,13). For low boundary
r=a-d, choose uniformly among all deletions of d selected outside slots.
For high boundary r=b+d, choose uniformly among all additions of d unselected
outside slots. Every destination is authentic by complete conditional charge
coverage. In the twelve-pin case keep z fixed and recompute the two1532
controls; in the single-pin case recompute176 from368.

For a low-side parent c with t_a=|c|=k-a the incoming load is bounded by

    1 + sum_{d=1}^a (A_{a-d}/B_a)
                         e_d(w on E minus c)/binom(t_a+d,d).     (2)

The denominator is the number of legal parents of its particular leaf,
whose outside size is t_a+d. For an upper-side parent with t_b=k-b it is

    1 + sum_{d=1}^{S-b} (A_{b+d}/B_b)
                         e_d(w^{-1} on c)/binom(m-t_b+d,d).     (3)

The A and B are the same actual coefficient maxima/minima as in the
un-pinned or pinned theorem, including control weights where required.
The residual charge is NOT averaged away: bounding each coefficient is
why these formulas hold for every specified charge q.

An entirely explicit source-sensitive bound replaces the restricted e_d
in(2)/(3) by e_d over all outside slots. No enumeration of visible states
is needed. If K^{-1}<=w_i<=K, a sharper closed bound is

    H_bal = 1 + max(L_k,U_k),
    L_k = sum_{d=1}^a (A_{a-d}/B_a) K^d
                      binom(m-t_a,d)/binom(t_a+d,d),
    U_k = sum_{d=1}^{S-b} (A_{b+d}/B_b) K^d
                      binom(t_b,d)/binom(m-t_b+d,d).       (4)

Omit a side when its parent size is outside[0,m]; it has no leaves. A binomial
with d exceeding its nonnegative upper index is zero. Low and high parents
belong to distinct grades, which licenses the maximum rather than a sum.

## 3. Uniform inverse-quadratic control in the bulk

If `t_a>=delta*m` and `m-t_b>=delta*m`, with fixed delta>0, each low ratio
in(4) is at most `((1-delta)/delta)^d`, as is each high ratio. This follows
factor by factor from the falling numerator and rising denominator. If a
side is absent its term is zero. The bank, control weights and K are fixed,
so H_bal is independent of m. Combining with the reference exchange bound
gives **gap >= c/(m+b-a)^2** on the ENTIRE accepted fibre in this occupancy
range. It retains boundary sources, unlike conditioning on the interior.

At extreme total occupancies the denominators can be one; the earlier
degree-seven (or degree-nine pinned) bound remains a valid all-grade bound.
This proof does not claim that those exponents are sharp. The cyclic heavy
boundary hostile still prevents a full bound independent of relative odds.

For the actual81-source binary branch with all source odds one, m=61,k=40,
the exact H_bal and full-fibre gap bound are recorded in `balanced.normal.json`.
They improve the deterministic-parent bound while keeping every original
source and all events within22/54. The weights of this illustrative uniform
law are stated explicitly; an actual nonuniform law uses its own coefficients.

## 4. Remaining target

This moves the bottleneck from interior weight mixing to extreme-occupancy
admission and native supply. A useful next problem is an exchange/circuit
bound within the sparse boundary layers that avoids evacuating through the
bank. Such a theorem must retain actual labels and pins: the cyclic-charge
hostile shows that arbitrary alphabets can require large direct exchanges.
It is also still necessary to compare any prescribed ROOT transport with the
displayed source-preserving rates. No deterministic Collatz entry follows
from the auxiliary reversible chain alone.
