# Arrival-frontier dissipation cones for a uniformly populated native relation bank

Session30, 2026-10-22. **PROVED, scoped:** an exact arrival-boundary identity, a negative-definite completed-relation cone, its sharp real Schur bound, and a source-preserving sufficient gate retaining exterior interactions. **FINITE-EXACT:** the declared literal source/clock/register controls. **OPEN:** a uniform exterior-plus-boundary estimate, native 153N for all scales and horizons, and arbitrary-source entry.

## 1. Inheritance and provenance

The two growing parent banks, exact inverse edges and counts are inherited from `05-knowledge/results/scale_credit_supply_20261021a/README.md`, sections4–5. The source carrier, frontier and energy are inherited from `session18_native_annulus/README.md`. Synchronous blocks and exact destination-conflict clustering are inherited from `future_tail_blocks_20261021b/README.md` and `tail_geometry_20261021c/README.md`. The new calculation below uses the **actual one-depth lag and linked clock signs**, which independent future block choices discard; it does not claim discovery of those banks or of quadratic completion.

Closest proved mechanism: the literal parent-child seam. Hostile: the child-alone arrival state gives positive increment1 or3 but is dynamically impossible. Corrected near miss: negative completed-pair energies do not authorize erasing exterior interactions. Least-used sidecar: which parents arrived **exactly at the present cut**, rather than eventual parent coverage.

The incoming values520738,1039055,3848613 remain user-reported in this packet, not independently recovered or audited. Targeted physically present canon/results text and JSON checks, followed by a read-only recovery scan of699 Collatz/tail/session20–29 Markdown paths in HEAD, did not recover those exact values. The accessible linked-tail source instead records551844,1118881,4142214 at native cuts105,110,130. This is a provenance gap, not a refutation of the incoming report. Neither set is used below. Composition of certified inequalities and supply of source-owned relations are different obligations.

## 2. Two source-injective relation families at every scale

Let Q(v)=v1²−v1v2+v2², β its polarization, and δ(w,h)=Q(h)+β(w,h). Use the original primitive sources and stopped map P, the same627-point frontier V (max63665), and original old/new charges.

For c=ceil(3N/2), retain:

- A: c≤u<2N, u mod72 in{9,27,63}; s=(16u+3)/3; P(s)=u, full edge clock4; old charge L0=(1,2), new charge U1=(−1,−1).
- B: N≤u<c, u mod72 in{17,35,71}; s=(64u+19)/3; P(s)=u, full edge clock6; old charge L2=(1,−1), new charge −L2=(−1,1).

Each family is injective in u. Their old ranges are disjoint, their new ranges lie respectively in(8N,16N) and(16N,32N), and old and new ranges are disjoint. Therefore **all original endpoints are distinct**, stronger than merely distinct pair addresses. The exact counts are inherited floor sums over the displayed residue classes. If144 divides N, each family has N/48 pairs. This is a genuinely linear original-source supply, not exponent density in a lacunary family.

For every N≥7959 all children exceed maxV. At smaller scales, delete pairs whose child lies in V before applying the arrival theorem. No such deletion is needed for N≥7959; at most627 deletions suffice at any smaller scale. A parent may itself lie in V: then it has actual arrival depth0 and clock0, and its child has depth1 and clock4 or6.

For every retained pair, if its parent never arrives, neither does the child. If the parent has actual first-arrival terminal t, depth d, and full clock k, then the child has terminal t, depth d+1, and clock k+4 or k+6. The sign is consequently the same. This theorem uses no assertion of eventual arrival. The mod15 and mod35 registers are always read from the literal u,s, not independent choices.

## 3. Exact completed-cone and one-layer boundary theorem

Fix any finite horizon d. At each terminal t define signed integer counts:

- x_t,y_t: A/B parents first arriving by d−1 (both endpoints now arrived), weighted by(−1)^parent_clock;
- a_t,b_t: A/B parents first arriving exactly at d (old parent present, new child still absent), with the same weights.

All later or nonarriving parents are absent. The bank's old and new vectors are

    W_t=(x_t+a_t)L0+(y_t+b_t)L2,
    H_t=x_t U1−y_t L2.

Thus its exact internal increment is

    F_t=−2x_t²+3x_t y_t−3y_t²−3a_t x_t+3a_t y_t−6b_t y_t.       (1)

This is a direct expansion in the original energy, not a new credit account.
For a_t=b_t=0,

    F_t=−2(x_t−3y_t/4)²−15y_t²/8≤0,                         (2)

with equality only x_t=y_t=0. Hence arbitrarily many completed pairs, at arbitrary actual terminals and with arbitrary actual clock signs, remain in a negative-definite two-dimensional cone. The previous individual β=−3/−6 observations did not alone imply this closure under collisions.

Only the **current one-arrival-layer boundary** can produce positive internal energy. Completing the square gives

    F_t≤(6/5)(a_t²−a_t b_t+4b_t²).                          (3)

Proof: let M=[[2,−3/2],[−3/2,3]] and l=(−3a,3a−6b). Then F=−zᵀMz+lᵀz. M is positive definite, determinant15/4; its real maximum is lᵀM⁻¹l/4=(6/5)(a²−ab+4b²). This is sharp over real z; it is not asserted to be the exact integer maximum or to have a native maximizer.

A useful scalar corollary: if m_t counts unsigned currently arriving A parents and n_t counts unsigned currently arriving B parents, then

    F_bank(d)≤(6/5) Σ_t(m_t²+m_t n_t+4n_t²).

Thus all past completed mass is controlled by the instantaneous relation boundary. This is an energy reduction, not a bound on that boundary's actual size.

## 4. Exterior interactions cannot be omitted: exact Schur gate

At a terminal t let O_t,H_t^ext be the actual arrived old/new vectors from **every original source outside the bank**. Put B_t=a_tL0+b_tL2, C_t=(O_t+B_t,H_t^ext), and define polarization on two-channel blocks

    κ((w,h),(v,j))=β(h,j)+β(w,j)+β(v,h).

Write A=(L0,U1), B=(L2,−L2), l1=κ(C_t,A), l2=κ(C_t,B). The complete original increment at this terminal is exactly

    δ(O_t,H_t^ext)+β(B_t,H_t^ext)
      −2x_t²+3x_ty_t−3y_t²+l1 x_t+l2 y_t.                  (4)

Consequently a sufficient pointwise gate is

    Σ_t[δ(O_t,H_t^ext)+β(B_t,H_t^ext)
          +(3l1²+3l1l2+2l2²)/15] ≤153N.                  (5)

Equation(5) retains every exterior diagonal, bank/exterior interaction and boundary contribution. It only relaxes the feasible integer completed counts to real variables. Therefore it introduces no new source norm, no repeated negative pair, and no assumed arrival. It is not guaranteed sharper than every earlier stopping gate; valid alternatives can be minimized, not added as credits.

For the existing ledger Δ=A_diag+R′+J′ and A_diag+M=153N, any valid bound G from(5) yields R′≤G−A_diag−J′. Passing G≤153N gives the **compensated** condition R′≤M−J′, not the stronger R′≤M unless separately proved.

Hostile: an exterior old vector O=(−10,0), no exterior new vector, and one completed B pair have total increment27, although the exterior alone contributes0 and the pair alone contributes−3. This algebraic hostile proves why negative-cone membership is not a free deletion rule; no native realization of this arbitrary exterior is claimed.

## 5. General source-owned relation profiles

A wider criterion is available without changing the carrier. Partition pending original-source blocks into disjoint authenticated relation clusters. In one cluster every current node x has a certified P-path to a common root r, of length l_x and full clock c_x. Require that none of the path nodes before r is in V. For root r outsideV, let T be its future first-arrival time, possibly infinite. The block from x arrives at T+l_x, at the same terminal, with sign(−1)^c_x times the root continuation sign.

Therefore its only possible future cluster vectors are the **depth-prefix profiles**

    p_j=Σ_(l_x≤j)(−1)^c_x b_x,

or zero, with a single common terminal and root sign. If r is already inV, its arrival is known instead. This is strictly more constrained than allowing each pending node an independent sign, absence and terminal. A finite cluster envelope may maximize over those profiles; cross-cluster terms must still be retained exactly once. The A/B theorem is the explicit linear-density one-edge instance, and its common negative-definite cone is what permits completed counts to be eliminated without enumerating them.

## 6. Precise remaining uniform obligation

A uniform relation supply is now explicit for A/B at all sufficiently large scales. What remains is **not** another arithmetic existence/count theorem for these two families: their exact density is already known. A sufficient route is to prove, for every N≥7959 and every finite horizon d, the exterior-plus-current-boundary bound(5), with separate treatment of the finite scale range4096≤N<7959; alternatively use(5) only until an independently authenticated all-future stopping gate succeeds.

More informative future work would enlarge the relation alphabet while preserving a negative-semidefinite completed Gram matrix and bounding the new boundary width. The first missing lemma is a quantitative bound on **outside channels plus current arrival boundary**, not eventual convergence of every selected parent. Pair supply alone cannot pay all cross terms.

The Session29 odd-exponent family u_k=(3^k+1)/4 remains a hostile control: density-one finite descent in the exponent parameter supplies at most two original sources in[4N,256N). That family cannot replace these linear-size banks in source accounting. Conversely, the present linear-size banks do not prove153N merely because they are numerous.

## 7. Reproduction and scope

Run `python verify.py --repo-root /path/to/math --out verification.json` and the same with `python -O`.

The checker uses explicit exceptions (not removable assertions), exact integer/rational arithmetic, independent literal B-chain charges and primitive edges. It checks14641 arrival-boundary cases,6561 exterior identities,29595 literal family rows, and5385 parent/child first-arrival pairs at scales63666 and65536. Every one of those finite pairs arrives within500 steps; that is a finite observation only. Their full clocks and literal mod15/mod35 registers enter a reproducible digest. It also checks exact counts on20 unbounded-ray sample scales and one50-digit off-ray scale; the proof of the floor formula, not those samples, is the all-N argument.

The output records544533 checks. Optimized replay must agree byte for byte. These tests do not evaluate the whole153N ledger at the tested scales and do not certify(5) uniformly. No new theorem identifier or literature priority claim is assigned.
