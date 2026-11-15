# Arbitrary-prefix carry locks: exact classification for two fixed arms

Work in the [binary-cut two-core glue](../chanpak_heisenberg_target_20261115/ANCHOR_DECODER.md), without requiring separated high and low parts. A fixed arm specifies two physical cores (U,W) and a binary suffix (v), and can be applied after **any** nonempty prefix (p) beginning 1. This note characterizes when two fixed arms differ by a constant absolute count for all prefixes. It gives reproducible non-separated locks, not a universal +2 move.

## 1. The common response forces comparable literal suffixes

Write (h_X=H(X)), (Q_X=Q_2(X)), (P_X=Q_X-h_X). The core coordinates are
```
d=P_U Q_W,    t=h_U h_W,    a=P_U h_W,    G=d-a+2t.
```
Use (P(w)=A(w)-1) for the excess apex count of a binary word, and (val(w)) for its binary value including any zero padding. Let (rev_k) reverse exactly (k) bits. For an empty suffix these values are zero. With (p) of value (x), (v) of width (k), value (b), and reverse (r), the exact count is
```
N(pv)=G+d b+tP(v)+L x+tP(p)+t e(v) last(p)
L=(d+t)2^k-t(1+r)
e(v) = 1 if v is nonempty and starts 0, else 0.
```
This follows either from the gap polynomial or the seam identity in the previous host decoder: weighted cross inversions between (p) and (v) give (x(2^k-1-r)), plus the adjacent extra exactly when (p) ends 1 and (v) starts 0.

If the difference of two arms is constant on all (p), testing (p=1^m) forces (L=L') because (P(p)=0) and (last(p)=1). Testing (p=10^m, m≥1), then forces (t=t') because (P(p)=2^m=x) and the last bit is zero. Conversely (L=L') and (t=t') give equal adjacency terms too: (L mod  2=e(v)) since (d) is even, (t) is odd, and (r mod 2) is the first suffix bit when nonempty. Thus ((L,t)) is a complete response for constant differences, not merely a sufficient invariant.

Having (t=t') and (L=L') implies, modulo (2^{min(k,k')}), that the reversed suffix values have identical initial low bits (since (t) is odd). Thus the literal suffixes are prefix-comparable. Relabel the arms with (k≤k'). Remove their shared suffix prefix (v); the rest is one empty arm and one extension by a word (z) of width (ell=k'-k). The response condition becomes
```
d = 2^ell d' + t rho,
rho = rev_ell(val(complement(z))).
```
If (z) is nonempty, its first bit must be 1: parity on the displayed equation forces (ρ) even. So the extra word has **no** seam adjacency correction, whatever preceding prefix ends it receives. In particular removing the shared suffix prefix is an exact physical response comparison, not an invalid cancellation of unknown counts. The difference (longer-suffix arm minus empty arm) is
```
Delta = a-a' - d' val(complement(z)) - t D(z)
D(z) = rho-P(z)
     = sum over two zero positions i<j of [2^(j-i-1)+[j=i+1]].
```
The last identity subtracts the actual preceding 1 contributions at each zero position (j) from their value if all earlier bits were 1. Thus (D(z)) is a completely explicit nonnegative penalty of having multiple zero marks.

This also gives a complete search over **all suffix widths** for any chosen finite catalogue of core pairs, without enumerating suffixes. Given (d≥d') and common (t), use (ell=⌊log_2(d/d')⌋), then check (ρ=(d-2^ell d')/t) is an even integer in ([0,2^ell)). The width is forced because (d'≥2t) and hence (2^ell≤d/d'<(3/2)2^ell) for nonzero width (the zero width case is (d=d')). Recover (z) by complementing the (ell)-bit reversal of (ρ); then calculate (Delta) as above. Suffixes of all ones have (ρ=D(z)=0) and can be absorbed into a leading transitive block of (W'), physically moving the boundary without changing arcs.

## 2. True switches for every prefix; carries are not being suppressed

Two empty-suffix choices give ((d,t)=(252,9)):
```
U rows (10,12,9,0), W rows (2,4,1): G=207.
U rows (2,0), W rows (14,12,8,16,7): G=243.
```
For every nonempty binary prefix of value (x) beginning 1, both (252x+9P(p)+207) and (252x+9P(p)+243) are Hamilton counts of actual glues. The difference +36 does not require a dyadic prefix or a separated window.

There is also a genuinely different-width comparison. The empty-suffix arm has transitive (U=TT_4) and (W) rows ((62,124,56,112,32,64,21)), with ((h_W,Q_W)=(31,644)). Thus ((d,t,G)=(9660,31,9257)). The other arm has (U'=TT_1), (W') rows ((76,109,8,32,47,5,60)) with ((h'_W,Q'_W)=(31,596)), and **extra physical suffix** (z=1101). Its raw coordinates are ((d',t,G')=(596,31,627)). Here (ρ=P(z)=4), (L'=596·16+31·4=9660), and its constant is (627+596·13+31·4=8499). This gives an exact −758 switch over every prefix. The zero in (z) prevents identifying this switch with a slide of transitive 1s.

Either switch therefore holds when one intentionally lets the low quotient carry: for example, set (p=10v) with (v) arbitrarily wide and large in binary value, and split at the first gap of length one. No bounding of the subsequent term below 2 is used. The two shifts having gcd 2 does **not** show they can be chained: their response states differ.

## 3. Two all-order obstructions to a universal +2 in this mechanism

**Equal original slopes, including extra all-one suffixes.** Absorb all-one extras into (W'), so both suffixes are empty. At common (d) and (t), the gap is (a-a'). The odd-cycle packing expansion gives (Q_X≡(-1)^{|X|}h_X mod 3): each packing of (j) odd cycles on (s) vertices contributes (2^j) to (h) and (2^{|X|-s+2j}) to (Q_X), and (s≡ j mod 2). Thus (a≡ epsilon_U t mod 3), with (epsilon_U=0) for even order and 1 for odd order. If (3∤t), the condition (d≡0 mod 3) is ≡alent to (epsilon_U=0), since (Q_W) is a unit modulo 3. If (3|t), both (a,a') vanish modulo 3 anyway. Both are odd, so their gap is a multiple of 6.

Moreover (Q_X≡0 mod 4) for nonsingleton (X): each packing exponent (|X|-s+2j) is at least 2. When both (U,U') are nonsingleton, (a≡ a'≡-t mod 4), and the gap is a multiple of 12. When both are singleton it is zero. When just (U) is singleton, (a=t) and a gap of magnitude 6 would force (h'_W(Q'_{U}-2h'_{U})=6). Here (h'_W) is 1 or 3. The cut floor (Q'_U≥ (|U'|+1)h'_U) bounds the possible order and forces (U'=TT_2,TT_3,) or (C_3); indeed orders at least four imply (h'_U=1) and (Q'_U=2^{|U'|}), which is too large. Path count one means transitive: a directed triangle contributes to the positive odd-packing formula. Count three means precisely one directed triangle as a nonsingleton SCC, by the same formula and [Busch's strong bound](https://doi.org/10.37236/1141), so its cover counts are (12·2^s). In the three cases (TT_2,TT_3,C_3), equality of (d) would equate a power of two or (12·2^s) to one with an additional odd factor 3 or 7, or (for (C_3)) compare odd parts 3 and 9. Each is impossible. Therefore every nonzero gap here is at least 12.

**A singleton (W) in the shorter-suffix arm.** Suppose (z≠∅) and the shorter-suffix arm has (W=TT_1). Then (d=2a≡2 mod 4), so (ρ≡2 mod 4) and (z) begins 10 (an extra suffix of width 1 would give (ρ=0)). The first 1 alone contributes (ρ/2+1) to (P(z)), including its adjacency bonus. Hence
```
Delta = d'(2^(ell-1)-val(complement(z))) + t(P(z)-rho/2) - a'
        >= d'/2+t.
```
Here (a'≤ d'/2) and the first parenthesis is at least 1. Equality with 2 would require (d'=2,t=1), hence both longer-arm cores singleton. Then the shorter cores are transitive too and (d=2(2^m-1)) for some (m≥1). This cannot lie between (2^{ell+1}+2) and (3·2^ell-2), the necessary range from (d=2^{ell+1}+ρ). So this channel never gives magnitude 2 either. Neither obstruction covers every possible matching width with both middle cores nonsingleton.

Replay: `python replay_equal_response.py` verifies the literal core counts and selected complete glues by independent subset-endpoint DP, the cancellation/penalty formulas on finite padded words, and the width decoder. The infinite assertions follow from the formulas.
