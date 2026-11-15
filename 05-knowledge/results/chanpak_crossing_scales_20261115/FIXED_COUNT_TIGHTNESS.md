# Tightness and a finite stopping rule at fixed Hamilton count

**Proved for all tournaments.** The equal-H switch of the [Heisenberg session](../chanpak_heisenberg_target_20261115/PHYSICAL_PROFILE_SWITCH.md) moves Q_2 with the same absolute H. The theorem below says what any continuation of that idea can and cannot do. It is a concrete tournament analogue of requiring tightness *before* promoting a limit argument. It makes no additional claim about Yang–Mills itself.

Let H=t be fixed. Every tournament is an ordinal sum of its strongly connected components (SCCs), and Q_j multiplies across such sums by splitting each path at the boundary. Strip all singleton SCCs, with H=1 and Q_j=j. Let n_0 be the remaining order. The triangle has H=3; for strong components at order m≥4, [Busch's lower bound](https://doi.org/10.37236/1141) gives H≥5^((m−1)/3)≥3^(m/3). Thus 3^n_0≤t^3. If t is prime, there is only one nonsingleton SCC and m≤1+3 log_5 t by Busch, also valid at t=3. The reduced bank is finite. If empty it has H=1,Q_j=1 by convention.

There are two useful exact consequences.

1. The simultaneous slot readout is a finite union of rays (q_2·2^s,q_3·3^s,q_4·4^s,…), with the *same* s≥0 in each coordinate. In particular, Q_2 by itself has finitely many odd parts. A uniform bound for them is B_t=floor(2^n_* (t+3)/4), where n_* is the reduced order cutoff. In the odd-cycle sum Q_2/2^n_0, the empty packing contributes 1, and the remaining terms are each at most one quarter of the corresponding H term, since their exponent deficit is positive and even. Hence q_2≤B_t before singleton inflation.
2. For a proposed q=2^v o with o odd, take b so 2^b>max(o,B_t). If a physical core has H=t **exactly** and Q_2≡q mod 2^(v+b), it must have Q_2=q as an integer: the valuation is v, and after dividing by 2^v the two odd parts below 2^b agree modulo 2^b. The promise H=t cannot be replaced by a congruence on H. This upgrades a fixed-H residue construction to an exact construction at an explicit stopping precision.

In particular in the bounded progression N=LQ_2+(M−L)t, for each specified t recover q=(N−(M−L)t)/L, reject nonpositive or nonintegral q, then apply the finite stopping condition or enumerate the bounded reduced cores and add the dictated number of singleton SCCs. There is no unbounded search at fixed t. This is not a guarantee that the progression has a suitable t or that the resulting finite bank contains its q.

## No indefinite exact additive actuator at fixed H

Fix H and suppose Q_2 repeatedly increases strictly via an affine integer update q↦a q+b with a,b integers. If a≤0 there is no strictly increasing infinite orbit (a=0 reaches a constant, and a<0 reverses order). If a is odd, the update permutes residues modulo every power of 2. Returning infinitely often to the initial class modulo 2^(v+1), where v is the initial 2-adic valuation, would yield infinitely many increasing q with valuation v. The finite odd-part bank forbids that. If a is even and b≠0, the orbit tends in Z_2 to b/(1−a)≠0, so valuations are eventually fixed, again impossible. If a is even and b=0, the odd part of a must be 1 or it grows without bound. Thus only multiplication by a positive power of 2 survives, physically realized by singleton padding. This is about infinite orbits of a fixed affine scalar update; switching h or using a target-dependent non-affine move is not excluded.

If r odd parts are available, an increasing arithmetic progression of Q_2 values at fixed H has at most 2r terms: at least half the terms lie in (M/2,M] for the largest term M, and each dyadic ray can supply at most one in that interval.

## Exact example H=9

The nontrivial SCC factors are either a single H=9 strong core or two triangles. In the first case Busch gives order at most five, and the literal census of labelled strong cores through five leaves only Q_2=84 at H=9. In the second case Q_2=12²=144. The full fiber is therefore precisely the rays 84·2^s and 144·2^s, for s≥0. In particular odd parts 21 and 9 are the entire bank. The 21 in the Q_2 channel is not an H=21 witness. See [replay](replay_fixed_count.py).

This result shows the difference between proving compactness and proving exact placement. The [November 10 Yang–Mills assessment](https://github.com/yang-mills-group/YM) uses compactness to obtain infinite-volume subsequences and does not assert uniqueness or a mass gap. In this finite enumerative setting the strong-component bound is extra information: it turns a would-be subsequence problem into finite identification at fixed H. There is no claim that the same bound or the same approach transfers back to gauge theory.
