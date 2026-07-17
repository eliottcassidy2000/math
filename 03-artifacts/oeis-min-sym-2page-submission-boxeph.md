# OEIS submission draft (boxeph): minimal crossings of rotation-invariant 2-page drawings of K_n
Sequence (n odd, 5..17): 5, 21, 72, 176, 377, 705, 1224.
Definition: min over partitions p of the gap set {1..(n-1)/2} into two pages of
#{quadruples a<b<c<d in Z_n : p(gap(a,c)) = p(gap(b,d))} — the crossing count of the
best circular 2-page book drawing of K_n whose page rule is rotation-invariant.
Comments: ratios to C(n,4) tend to ~1/2 (the random-page law); Guy's conjectured
cr(K_n) = Z(n) scales at 3/8: the two densities cross uniquely at n = 9 where the
term 72 = 2·Z(9) — the triple-tangency point (cf. A000241, A001110). Program: repo
04-computation/lrc14_* (boxeph S38/S39/S52 referees).
