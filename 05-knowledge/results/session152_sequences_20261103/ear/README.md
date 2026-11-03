# Order-five ear image: exact finite certificate

The number set is precisely the image of the nonconstant one-vertex cuts of the strong order-five tournament code 8 in THM-4118; it is not evidence of an arithmetic or cross-domain correspondence. H counts directed Hamiltonian paths. Label 0,...,4 and direct edges (1,0),(2,0),(3,0),(0,4),(2,1),(3,1),(4,1),(3,2),(4,2),(4,3). Its H is 9. Let x -> i iff i is in S. Literal subset path DP yields

P(z,t) = sum_{S subset V} z^{|S|} t^{H(T+x_S)}
= (1+z^5)t^9 + (z+z^4) A(t) + (z^2+z^3) B(t),
A=t^15+t^17+2t^23+t^33,
B=t^17+t^19+t^23+t^25+t^27+t^29+t^33+2t^37+t^41.

The equality of cardinality layers is specific to this parent; F(S) need not equal F(V\S), e.g. F({0})=17, F(V\{0})=23. The univariate nonconstant-cut polynomial is 2(A+B), giving frequencies 2,4,2,6,2,2,2,4,4,2 at 15,17,19,23,25,27,29,33,37,41 respectively. This polynomial counts labeled cuts, not distinct children up to isomorphism.

For a general tournament, the exact stateful recurrence is D({v},v)=1 and D(U,v)=sum_{u in U\{v}} D(U\{v},u) 1_{u->v}; H(T[U])=sum_v D(U,v). Equivalently with F(S)=H(T+x_S), H=H(T), A_i=F({i})-H, K_ij=F({i})+F({j})-F({i,j})-H,
F(S)=H+sum_{i in S} A_i-sum_{i<j in S}K_ij.
This is the exact exposed-slot quadratic (THM-4118); its quadratic nature follows by splitting any child Hamiltonian path around x, so the insertion depends on zero, one or two incident indicators. In the above labels A=(8,6,14,24,14) and K_01,02,03,04,12,13,14,23,24,34=(6,4,4,8,10,10,4,10,4,6). This permits a compact exhaustive check independent of enumerating permutations.

Inductive reduction, with the theorem's cited Moon strong-subtournament result: all strong tournaments of order >=4 admit a strong one-vertex deletion, so in principle their next-order image is the union of the nonconstant-cut images of every strong parent. However a recurrence retaining H alone is not well-defined: code 8 and code 759 both have H=9 but different images, the latter containing 31,43 and missing 25,27,41. Thus an induction aimed at the global odds-except-{7,21} conjecture must retain enough adjacency/boundary information, prove coverage/overlaps, and separately establish attainment at unbounded values. The n=5 ear image does not supply that induction. Existing all-order omissions of 7,21 and global completeness open are stated in THM-1370 and THM-4094; this computation does not prove completeness.

Supplementary finite check: enumerating all 1024 labeled order-5 tournaments yields 544 strong parents in six image families, with H multiplicities 9:240, 11:120, 13:120, 15:64. Their ear-image union is {15,17,19,23,25,27,29,31,33,37,41,43,45}; notably 35 and 39 also fail at this fixed order. These are finite results only. The accompanying `test.py` literally constructs children and runs subset path DP; further output contains full cut and family witnesses.
