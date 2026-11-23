# Smallest scalable templates and pure-curvature sevenfold squares

**FINITE-EXACT** for the census; **PROVED** for the all-filler identities and modular consequences. The template keeps the edge \(b\to a\), places a uniform nonempty tournament module \(M\) in one gap of a transitive spine with \(k\) fixed vertices, and permits two distinct apex-to-fixed-vertex reversals only. It requires \(H_{11}=7H_{00}\) with both edit orders strictly increasing for every filler. This is not a minimum for arbitrary two-edit tournaments.

## Exact search and all-filler certificate

For each \(k\), enumerate \((k+1)4^{k+1}\binom{2k}{2}\) descriptions: tag each fixed vertex and the module by a mask in \(\{0,1,2,3\}\), choose the gap, and choose the two contact bits excluding the module. Bit 1 of the mask value means \(a\) beats the site, bit 2 means \(b\) beats it. The [checker](template-scan.cpp) caches the exact subset/path-end DP at width one. Retain only exact ratio and both strict inequalities there; these are necessary. Test the remaining candidates at transitive-module widths two and three.

To see why those three probes certify the ratio for **every** filler, put \(D=H(M)\), \(Q=Q_2(M)\), \(R=Q_3(M)\). Fix the order of apices, assign each fixed site to one of the three spine runs, and choose which runs contain a nonempty module part. For one, two and three specified occupied runs the multiplicities are \(D,\ Q-2D,\ R-3Q+3D\). Boundary contacts depend only on the tag, so each corner is an integer form \(aD+bQ+cR\). For the three transitive probes the coordinates are \((1,2^t,3^t)\); the evaluation matrix is invertible. Thus exact sevenfold equality on these probes is an identity of readings. For strict monotonicity check all four increase vectors in the nonempty-slot basis: the coefficient of \(D\) must be strictly positive and the other two nonnegative. This sufficient test certified **every** retained candidate through \(k=7\); none was unresolved.

| Fixed sites \(k\) | Marked descriptions | Lower reads changing with transitive width (number of descriptions) |
|---|---:|---|
| 0–4 | see replay | none, even without the variable-width requirement |
| 5 | 1,105,920 | none; the only reads are \(5D,11D,19D\) |
| 6 | 7,569,408 | \(5(Q+D)\) (6), \(6Q+5D\) (2), \(15Q+D\) (4) |
| 7 | 47,710,208 | \(48Q-D\) (4), \(15Q+D\) (8), \(30Q+D\) (8), \(5(Q+D)\) (12), \(6Q+5D\) (4), \(5(2Q+D)\) (12), \(12Q+5D\) (4), \(9(Q+D)\) (2), \(11Q+9D\) (4), \(15(Q+D)\) (4), \(29Q+15D\) (2) |

At \(k=6\), six is the minimum for variable width, and A/B are exactly the prime-capable choices. The third frame uses tags \(2,M_0,2,2,3,2,2\), toggling both apex contacts to the second fixed site after \(M\). The corners are \((5,11,15,35)(Q+D)\), always composite below, and the half-mixed derivative is \(7(Q+D)\), singular mod 7. The \(k=7\) forms \(12Q+5D,\ 30Q+D,\ 5(2Q+D)\) are obtained from earlier forms by multiplying \(Q\) by 2 via one singleton added ordinally to the module. None of the width-dependent reads through seven involves \(R\).

## Three genuinely new lower forms at seven sites

In the table, tags include a single placeholder at zero-based index \(g\). Replacing that site by a uniform module gives the desired frame. Contact bit \(e\) toggles the apex \(a\) for even \(e\), \(b\) for odd \(e\), at site \(\lfloor e/2\rfloor\). The half-mixed derivative is \(d=(H_{11}+H_{00}-H_{10}-H_{01})/2\).

| Tags, \(g\), marks | \(H_{00}, H_{10}, H_{01}, H_{11}\) | \(d\) |
|---|---|---:|
| \(\mathtt{00020131},6,(1,8)\) | \(11Q+9D,\ 75Q+53D,\ 13Q+11D,\ 77Q+63D\) | \(4D\) |
| \(\mathtt{00203132},6,(0,9)\) | \(29Q+15D,\ 197Q+95D,\ 35Q+19D,\ 203Q+105D\) | \(3D\) |
| \(\mathtt{01023113},4,(1,15)\) | \(48Q-D,\ 144Q+13D,\ 112Q-21D,\ 336Q-7D\) | \(64Q\) |

Each new form gives a prime **outside both earlier A/B frames even with arbitrary fillers**: respectively \(97\to679,\ 73\to511,\ 191\to1337\). At filler order one the old starts are 17 and 31. At order \(n\ge2\), \(Q\ge3D\) by the cut bound and \(4\mid Q\) by the two-slot odd-cycle expansion. In A a candidate \(p\) would require \(D\le p/23\) and \(D\equiv5p\bmod24\); none of \(p=73,97,191\) passes. In B the conditions \(D\le p/46\) and \(D\equiv p\bmod15\) already exclude each. These are obstructions to A/B, not to all older unrestricted witnesses. The second frame gives \(73\to511\) at \(Q=2,D=1\); the third gives the prime seams \(191\to1337\) and \(383\to2681\) at \(Q=4,8,\ D=1\). No assertion of infinitely many prime members or minimal unrestricted edit distance is needed. Here the interaction is a **pure filler reading**: it sees either \(D\) alone or \(Q\) alone, unlike the earlier two mixed forms.

In particular the last frame maps \((Q,D)\) to \((h,d)=(48Q-D,64Q)\) by an integer map of determinant 64. The [physical slot controller](../chanpak_fork_curvature_20261123/ALL-PRIME-WORDS.md#complete-two-reading-image-including-powers-of-three) therefore prescribes both coordinates modulo **every odd modulus**, subject only to \(d\equiv-h\bmod3\) for even filler order or \(d\equiv h\bmod3\) for odd filler order when 3 is present. Unit \(h\bmod3\) can now coexist with unit \(d\bmod3\), so the gradient chart is nonsingular mod 3, unlike A/B. This is an exact-ratio construction, not admission at an arbitrary integer \(h\).

For the same frame, a simpler wheel works at any prescribed unit cofactor phase \(r\) modulo a squarefree odd wheel. Set \(D=48Q-r\); the two middle readings then have residues \(768Q-13r\) and \(7(3r-128Q)\). At each prime \(\ell\ge5\), exclude the roots of \(D,\ Q,\ 768Q-13r,\ 3r-128Q\): at most four forbidden points. At three choose even filler parity and \(D=Q=-r\), so the same normalized readings are units. The simultaneous selector gives an exact square with rough cofactor and two rough middle readings after stripping the fixed 7 from the second, and with \(d\) a unit at every odd wheel prime. Dividing the arithmetic reading by 7 does not itself construct the divided height.

These three examples separate curvature scales too. On \(n\) filler vertices, \(Q\ge(n+1)D\) by all cuts of each Hamilton order; different orders and cuts give distinct ordered covers. The first two frames therefore have \(d/h\to0\) uniformly as \(n\to\infty\), while the third has \(d/h\to4/3\), with errors \(O(1/n)\). More generally in A, \(k_A=Q+D\) obeys \(0<k_A-h/6\le h/[6(6n+11)]\); in B, \(0<4h/5-d\le34h/[5(15n+16)]\). In A the actual half-mixed derivative is \(-3k_A\). These are bounds at the same physical height and same marks, not limits across independent endpoint witnesses.

Run `g++ -O3 -std=c++17 template-scan.cpp -o template-scan && ./template-scan 7` to reproduce all counts and representative corner tables. A search for \(Q_3\) dependence must now allow more fixed sites or change the mark/filler layout.
