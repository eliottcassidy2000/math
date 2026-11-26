# Weighted expanding passages and the exact first-crossing carry

2026-11-26. Extends the [Collatz passage-flags and Erdős 1135 packet](../collatz_passage_flags_and_farey_gate_20261126/README.md), not the unrestricted conjecture. There are two independent certificates: quantitative exclusion of positive cycles in a specified compressed grammar, and a first-descent implication for the shortcut map. The stronger cyclic period lower bound cannot be imported into an arbitrary first-descent argument.

## 1. Put both clocks on the same height

Use the old native guard for an expanding packet: a repeated valuation word of total valuation b and odd length a, with multiplier 3^a/2^b>1. Keep the actual packet start x and Y=log2(x+C)+L. For this packet set β=a/b, g=158496251/100000000 > log2(3), and γ=gβ. The native guard gives at most βY odd steps in the packet and Y_next ≤ γY. Every individual valuation at least two lowers this *same* shifted height by more than μ=83/200 if x≥N=2^71: the rational check is

    (4(N+C))^200 > 2^83 (3N+1+4C)^200.

The upper bound for (next+C)/(x+C) decreases with x, so checking N suffices. The cyclic logarithms now telescope in the same coordinates. For r expanding packets put P_0=1 and P_i=γ_1···γ_i. If H bounds the starting height, the period satisfies

    k ≤ H [ Σ_i β_i P_(i−1) + (P_r−1)/μ ].

The contraction term telescopes as Σ(γ_i−1)P_(i−1)=P_r−1; there is no extra charge per packet for changing height coordinates. The bound increases in each β. For fixed β multiset the largest bound uses descending β: an adjacent swap changes the local two-term expression by its common prefix product times β_i−β_(i+1), because γ_i is proportional to β_i. Thus a few aggressive motifs need not be charged as if every motif were aggressive.

### Exact constants and period cutoff

For the finite motif library b≤32, with any stated β cap, take the maxima over the feasible integer pairs (a,b). The largest additive numerator is

    B_max = 3^(a−1) + 2^(b−a) [3^a−2^a−3^(a−1)].

For a fixed pair each coefficient is positive and depends on a prefix valuation; putting all b−a excess at the first valuation simultaneously maximizes every nonempty prefix. The resulting anchor bound C=max ceil(B_max/(3^a−2^b)) is 17,920 in all three libraries. The loose upper bound L=max ceil(log2(3^a−2^b)) is 51, 33, or 27 for the three caps. An internal cut at the cyclic minimum occurs in the first repetition by monotonicity at a fixed phase. A proper prefix of i≥1 letters has inverse slope at most 2^(b−a)(2/3)^i, giving a uniform cut charge of 11 bits. The reciprocal-sum constant M=max[a·3^a/(3^a−2^b)] is at most 6377292/7153.

As in the predecessor, the logarithmic gap estimate of [Simons–de Weger, Lemma 12](https://deweger.net/papers/%5B35%5DSidW-3n%2B1-ActaArith%5B2005%5D.pdf) gives n_min < 2^8 min(k^14.3,r(M+5)k^13.3). This uses the bound on reciprocals in each expanding packet or contracting stretch; it is independent of the number of individual contracting moves. Consequently one can take

    H = 20 + L + min(14.3 log2 k, 13.3 log2 k + ceil(log2[r(M+5)])).

[Hercher, Corollary 29, with linked corrigendum](https://cs.uwaterloo.ca/journals/JIS/VOL26/Hercher/hercher5.html) gives k≥137,500,000,000 once verification covers 3·2^69. [Bařina's published 2^71](https://link.springer.com/article/10.1007/s11227-025-07337-0) suffices. Dividing the upper bound by k gives a decreasing function beyond that barrier; compare there with outward-rounded BigInt interval arithmetic. This excludes the following parsings of any nontrivial positive cycle:

| β cap on every expanding motif with b≤32 | Maximum excluded packet count | Rounded period upper bound at the barrier |
|---|---:|---:|
| 1 | 38 | 95,123,825,321 |
| 2/3 | 302 | 134,648,219,857 |
| 7/11 | 1,747 | 137,158,439,782 |

One extra packet in each row exceeds the barrier in this test, not in general. If every other motif is capped at 7/11, permit at most 1, 2, 5 or 10 arbitrary-slope motifs: the excluded **total** packet counts are 1,689, 1,636, 1,478 and 1,215 respectively. Individual contracting steps remain unrestricted. This neither bounds arbitrary cycle descriptions nor proves convergence.

## 2. Exactly which parity word maximizes carry at a first drop

Use the Terras shortcut map. Write α=log2(3); let p be the first length at which the trajectory's cumulative coefficient drops below 1, with j≥1 odd steps at zero-indexed times i_m, 1≤m≤j. This drop is an even step, so 2^(p−1)≤3^j<2^p and p=ceil(αj). Before odd step m the coefficient has not dropped, so i_m≤floor(α(m−1)). All these upper positions can occur together: at each next odd time and at p−1 the current coefficient is at least one, which also checks the earlier even steps between them. The resulting binary parity word is native on its unique Terras cylinder. At the end, the normalized affine carry is Σ_m 2^i_m/3^m. Every summand is strictly increasing in its own odd time, so the maximizing word is unique, with carry

    (1/3) Σ_(t=0)^(j−1) 2^(−{αt}).

This is a candidate self-contained lemma for the [published Erdős 1135 parity/affine API](https://conjectures.io/problems/erdos1135-erdos-1135). No Lean check or external submission is claimed.

The predecessor's first-drop bound already forces j≥72,057,431,991 if actual descent fails from n≥N. For that range the average in the displayed sum is less than 0.7214. Indeed, use q=190537 and the certified |α−301994/q|<1/q². Every block of q rotated fractional parts is within 1/q on the circle of a shifted equally spaced grid. Consequently the empirical cumulative distribution at t is at most t+3/q. Integrating against the decreasing function 2^(−t), this raises the block average above its integral 1/(2 ln 2) by at most 3/(2q), and the incomplete block contributes at most q/j. The exact upper sum is less than 7214/10000.

If actual descent fails at the first drop from n≥cN, it follows that

    0 < p/j − α ≤ (7214/10000) / (3 c N ln 2).

Use the same rational Farey certificate as the predecessor. A bound on the denominator gives the bound on p by p>αj. Applying it at various source heights gives:

| Source lower bound | First coefficient drop forces descent through shortcut step |
|---|---:|
| N | 114,208,327,603 |
| 1.33163 N | 217,976,794,616 |
| 2 N, 4 N, or 8 N | 217,976,794,616 |
| 16 N = 2^75 | 1,411,629,234,714 |

The higher gates do not apply uniformly to the least counterexample unless one independently raises its source lower bound accordingly. They condition on a coefficient drop, not its existence. The exact arithmetic and brute small-word controls are in **node weighted.js**. All logarithmic rounding is directed, not floating-point.
