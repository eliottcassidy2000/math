# Quantitative escape from the transitive multiport chart

2026-11-04. **PROVED:** an effective bound on the forced affine digit, and impossibility of an infinite bounded-digit lineage inside the transitive-prefix/transitive-suffix chart, even with varying integer clock shifts. This is a boundary for a specified operation language, not a global Hamilton-spectrum obstruction or a two-digit construction.

## Sources and increment

The exact chart and physical joint-port formula are inherited from `multiport_defect_grammar_20261103/TRANSITIVE_CHART.md`, read at pinned main `3f2fa002801868c280fafed0c05a1d9060df6c86`, and S182 `induction/PROOF.md`. The former already proves at most one uniform digit per source chart, shift invariance, and finite absence of +2/uniqueness of +4 in a bounded scan. We do not claim those as new. The new ingredient here is a quantitative reversed-bit separation estimate, which turns bounded digit questions into a finite certified problem and excludes every infinite bounded-alphabet lineage. The companion S183 `other_targets` packet independently closes the +2/+4 classification by an exact binary induction. The estimate here provides a separate finite-cutoff proof of the same classification and a general bounded-alphabet theorem.

## 1. Actual chart and inherited identities

For a positive mask f with m bits, take an actual named suffix TT_m with its physical port set W given by the bits of f. The highest vertex is a port. Attach a TT_t prefix and a root r as in S182: prefix beats suffix, ports beat r, and r beats prefix and all nonports. Every such attachment is strong by the terminal-port criterion. Write J(f) for the spanning ordered two-slot count with first endpoint and second start both in W. Then

    H_t(f) = 1 - J(f) + 2^t f,
    J(f) = sum_{i<k; b_i=b_k=1} (2^(k-i-1) + [k=i+1]),
    J(2^a f) = J(f).

The directed joint contributions are individually physical: for i<k in W, the contribution whose first endpoint is k and second start is i is 2^(k-i-1); the opposite directed contribution is [k=i+1]. Equal endpoints contribute zero. Thus J is not a phase or an independently adjustable scalar.

An all-clock relation H_(t+k)(f') = 3 H_t(f)+d forces

    2^k f'=3f,    d+2=3J(f)-J(3f).                 (1)

Here k may be any integer for which both clocks are legal, and (1) is forced already by equality at two legal input clocks. These facts are inherited, not inferred from the finite controls below.

## 2. Reversed-bit estimate

For odd u>0 with m bits let

    R(u)=sum_i b_i 2^(-i),
    w=sum_i b_i,
    a=sum_i b_i b_(i+1),
    B=sum_{i<k} b_i b_k 2^(i-k).

The exact product expansion and the physical joint formula give

    J(u) - u R(u)/2 = a - (w+B)/2.                 (2)

Indeed the terms of uR(u) with larger numerator index contribute twice the non-adjacency part of J, its diagonal contributes w, and its remaining terms contribute B. For each fixed k the sum over i<k is less than one. Hence 0<=B<w<=m, while 0<=a<=w-1. In particular

    |J(u)/u - R(u)/2| <= m/u.                      (3)

This deliberately loose error bound is sufficient and holds for all positive u, not only odd ones.

For odd u, the first three low bits place R(u)/2 in these intervals:

    u mod 8 = 1: [1/2,5/8),    3: [3/4,7/8),
    u mod 8 = 5: [5/8,3/4),    7: [7/8,1).

The tail after the first three bits is strictly less than 1/8. Multiplication by 3 modulo 8 swaps 1 with 3 and 5 with 7. Therefore

    |R(3u)/2 - R(u)/2| >= 1/8.

Apply (3) to u and 3u, and let m'=bitlength(3u)<=m+2. The reverse triangle inequality yields the effective escape estimate

    |J(3u)-3J(u)| >= 3u/8 - 3m - m'
                         >= 3u/8 - 4m - 2.       (4)

All bounds hold without any asymptotic or finite-search assumption.

## 3. Effective bounded-alphabet exclusion

For any proposed digit d, reduce f to its odd part u. Equations (1) and (4) imply

    |d+2| >= 3u/8 - 4 bitlength(u) - 2.            (5)

Consequently, for any fixed finite digit alphabet D, the possible odd source masks are bounded effectively. One can choose M so that

    3*2^(M-1)/8 - 4M - 2 > max_{d in D}|d+2|

and the dyadic lower bound remains increasing thereafter. This gives the certified bound u<2^(M-1). No unbounded census is required. The dyadic increment is 3*2^(M-1)/8-4, positive already for M>=5.

For D={2,4}, M=8 works: at u>=128 the lower bound is at least 48-32-2=14>6, and increases at all later dyadic starts. Thus only the 64 odd masks below128 can possibly occur. The exact checker records those finite masks (only u=7 has d=4; none has d=2). Combined with shift invariance this is an independent proof of all-mask +4 uniqueness and absence of +2 within this chart. The companion classification packet supplies a different exact binary-induction proof. This packet's central result is the general estimate and arbitrary bounded alphabet consequence.

## 4. No infinite bounded-digit lineage, with varying clock shifts

Suppose actual charts f_n, arbitrary integer shifts k_n, and digits d_n in a fixed bounded set satisfy, at each n,

    H_(t+k_n)(f_(n+1)) = 3 H_t(f_n)+d_n

at two legal clocks (and therefore all legal clocks). Equation (1) implies that the odd parts obey

    u_(n+1)=3u_n, so u_n=3^n u_0.

This contradicts the effective boundedness established in Section 3. Thus no such infinite lineage exists, even if the shifts, port sets, suffix orders, physical realizations and bounded digits vary at every step.

An equivalent asymptotic proof is informative: (1) makes J(u_n)/u_n a convergent sequence, because its successive differences equal -(d_n+2)/u_(n+1), an absolutely summable sequence. But the reversed-bit cylinders alternate between intervals separated by1/8, and the errors in (3) tend to zero. Therefore J(u_n)/u_n cannot converge.

## 5. Limits, physical ports, and fields

The theorem is restricted to transitive suffixes and all-clock affine returns. It excludes neither a one-clock physical return nor any transformation with a nontransitive output suffix. In particular the current-main `transport/TWO_DIGIT_ENTRY.md` correctly realizes +2 from the all-marked TT3 by replacing its terminal singleton with a marked directed triangle. That output leaves this chart, as it must; the +4 output TT5 remains in the chart for one step. No claim of general two-digit closure follows from either construction or this obstruction.

Changing finite-field observer labels cannot change the physical joint count or the necessary count relation (1). No vertices are identified because their phases agree, and no ancestor phase is assigned to a child. This theorem constructs no new child, so there is no new phase-transport claim. Native guards are also outside its scope.

## 6. Reproducibility

`verify.py` independently counts every directed physical joint by literal ordered subset covers for all masks through order8, checks the exact bit-product identity with rational arithmetic, and audits (4) on a larger finite sample. It verifies the bounded cutoff, all64 odd masks below128, arbitrary legal clock shifts over small masks, and an explicit large-mask orbit. The all-order results rely on Sections2–4, not these tests. All guards are explicit exceptions and survive `python -O`; both modes produce identical JSON.
