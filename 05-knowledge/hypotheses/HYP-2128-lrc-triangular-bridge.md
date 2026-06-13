---
id: HYP-2128
status: SYNTHESIS + exact identity 8*C(n,2)+1=(2n-1)^2 + creative hypothesis H
source: opus-2026-06-03-S586
related: [HYP-2134, HYP-2133, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2127, HYP-2126, HYP-2097, HYP-2116, HYP-2117, THM-401]
---
See `07-reflections/lrc-triangular-bridge-add-mult-odd-even-s586.md`. Triangular = add/mult bridge (T_k=1+..+k=k(k+1)/2, odd*even/2). Arcs C(n,2)=T_{n-1}. +2 ladder adds odd gnomon 2n+1 (->squares, T_k+T_{k-1}=k^2). x2 doubling geometric/2-adic. KEY: 8*C(n,2)+1=(2n-1)^2 => 2n-1=sqrt(8*pairs+1) = odd-square root (8=2^3); n=14: 27^2=729=8*91+1, 2n-1=27 = the modulus of the 64 self-converse worry-classes (S570). The worry-modulus Z/(2n-1) is the odd-square face of the triangular pair-count. 7*T_k=7,21,42,70 (21=C(7,2), 42=heptagon triangulations). H: additive-face (2n-1 modular transversal) survives while multiplicative/doubling face fragments at 2q (S585).

Power-tower extension: the same `(n+1,n)` block geometry persists for all power
weights `k^p` if one lets the anchor move. In midpoint coordinates
`c=a+n, u=n(n+1)`, the balancing anchor satisfies
`a_p(n)=p*n^2+(p-1)*n+alpha_p+beta_p/u+O(u^-2)` with
`alpha_p=(p-1)(p-2)/(12p)`, so the additive tower (`p=1`) and square tower
(`p=2`) are the first two exact integer members of a wider asymptotic family.
See `03-artifacts/drafts/triangular-power-anchor-asymptotics.md`.
