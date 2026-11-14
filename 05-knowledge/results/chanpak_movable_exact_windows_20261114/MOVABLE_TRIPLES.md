# Two independently movable exact three-odd windows

Let R be an old tournament, z an apex, U its set of outneighbors in R, and K the Hamiltonian-path count with the apex. Define h as the count for R and f as the number of ordered full two-path covers of R, allowing an empty part (with one empty path). Let P count those covers with the second path nonempty and starting in U; let Q count those with the first path nonempty and ending outside U. Let S count full paths of R ending in U and T full paths starting outside U.

Join transitive wings as L → R → B, all arcs L → B. The apex's outneighbor bits on L are the usual binary bits of X in transitive order; its bits on B are the reverse complement of the usual binary bits of Y. For nonnegative X,Y set p(X)=X mod 2, g(0)=1, g(2X)=g(X)+X+p(X), g(2X+1)=g(X). Then literal whole Hamiltonian count is
[
K+h(g(X)+g(Y)-2)+QX+PY+fXY+Tp(X)+Sp(Y).                 (1)
]
Widths can include leading address zeros. Minimal widths at positive X,Y ensure strong connectivity: the apex beats the first left vertex, and the last right vertex beats the apex; within the ordinary join every vertex can be reached from the first left vertex and can reach the last right vertex.

For a direct proof of (1), split an apex-using path into ordinary subpaths numbered 1 and 2. If L2 is nonempty, its minimum must be an outneighbor, with X possible assignments of L to the two subpaths. Dually a nonempty B1 has Y possible assignments. Both nonempty give fXY. When only L2 is nonempty and the first core portion is nonempty, this gives QX; when only B1 is nonempty and the second core portion is nonempty, this gives PY. If the applicable core portion is instead empty, the wing must itself support both apex conditions, allowing an empty subpath. The left number D(X) now requires both a nonempty L2 with its minimum colored 1 and L1 empty or its maximum colored 0. Appending bit b sends D(z) to D(z)+(1-b)z+b(1-p(z)), so D(X)=g(X)-1+p(X); reversal gives the right counterpart. These two cases contribute h(D(X)+D(Y)). With both L2 and B1 empty, the original apex paths with an empty first core portion lose (h-T)p(X), and those with an empty second core portion lose (h-S)p(Y); the remaining apex paths are unchanged. This proves (1).

## Exact finite certificates, all-size transfer

The outgoing-bitmask rows below include the apex in the last row. A direct subset DP checks each tuple (h,f,P,Q,S,T,K); it is implemented in replay_independent.py.

| Whole rows | (h,f,P,Q,S,T,K) |
|---|---|
| (22,108,24,113,98,69,5) | (45,432,129,258,15,30,135) |
| (14,12,120,112,99,67,3) | (45,432,129,258,15,30,137) |
| (14,76,120,48,99,67,9) | (45,432,129,258,15,30,139) |
| (1022,508,952,1008,480,960,900,256,0,402) | (9,1344,427,245,9,9,89) |
| (1022,988,504,496,992,450,896,256,0,428) | (9,1344,427,245,9,9,91) |
| (1022,508,1016,368,992,448,384,776,0,362) | (9,1344,427,245,9,9,93) |

Thus choosing the puncture changes only the Hamiltonian path count, not any of the six coordinates subsequently used by arbitrary independent left/right wings. With K centered at 137 or 91, (1) gives an exact window {B-2,B,B+2} at every such center B, including one representative of each odd residue modulo 6. For h45 the first puncture core is prime, while the second and third are the same tournament C3[TT2,TT2,TT2]. Interface equivalence does not require identical old modular decomposition.

## Exact inverse and remaining absolute obligation

For a target N, the smaller positive address x satisfies f x² ≤ N. Orient this as the left wing, or swap Q↔P and S↔T to test the other orientation. For each x and each of the three initial K values, put λ=f x+P and D=K+h(g(x)-1)+Qx+Tp(x). The unknown y satisfies y≤Y=floor((N-D)/λ). Pick dyadic M≥2 such that λM>h(floor(Y/M)+1). On every remainder class y=Mz+v, unrolling the recurrence gives g(Mz+v)=g(z)+g(v)-1+z(M-1-rev(v))+p(z)(1-msb(v)), where rev reverses exactly log₂ M bits. The elementary bound g(z+1)-g(z)≥-z follows directly from the recurrence. Hence the target expression increases strictly in z, and each class is decided by binary search. The sum of the tested residue counts is O(√N log N) for fixed profile since M=O(1+√(hN)/λ). Including per-state arithmetic evaluation gives O(√N log³ N) basic integer operations. This is a complete inverse **within these two charts**, not a general tournament nonexistence procedure.

The center itself can be placed within O(√N) of any sufficiently large prescribed height: choose x=floor(√(N/f)) and y=floor(N/(fx)); then fxy ≤ N < fxy+fx, and 1≤g(u)≤u+1 bounds the other addends by O(x+y+1). The triple chooses the phase modulo 6 without disrupting any outer-wing coordinates. This does not shrink the absolute error to <6 and does not place an exact window at every target. The next target is a profile-equivalent family with K-window radius growing proportionally to the outer-wing mesh, or an independent adjustment of the center at fixed external response.

## Joining two profile-equivalent families at a shared apex

Write (h,f,P,Q,S,T,K) for one puncture and (j,g,p,q,s,t,L) for another. Put the two old cores in ordinal order, and use a single apex with the original two masks on their respective cores. The composed six-coordinate interface is
[
(hj,fg,Pg+hp,fq+Qj,hs,Tj).
]
The count is
[
K'=j(K-(h-S)) + h(L-(j-t)) + Pq + (h-S)(j-t).
]
Indeed, partition full paths at the common apex and split into four cases according to whether the second subpath meets the first core and whether the first subpath meets the second core. If both do, the conditions separate into P and q, giving Pq. If just the first condition holds, apex paths on the first puncture with nonempty second subpath can extend by all full paths of the second core, giving j(K-(h-S)); the other singleton case is dual. If neither condition holds, each core is entirely on its own side of the apex, giving (h-S)(j-t).

This gives a calculus for amplifying exact windows at identical interfaces: changing the two input counts by d and e changes the output by jd+he. For r copies of the h9 triple, choose input changes in {-2,0,2}; the output has an exactly certified ladder of 2r+1 successive points of spacing 2·9^(r-1), retaining the same six coordinates at every point. It enlarges the window in a sublattice only. In particular, stacking cannot be mistaken for radius growth at step 2; a supply of different compatible count-preserving interfaces or a second instrument that fills the enlarged gaps would still be needed.
