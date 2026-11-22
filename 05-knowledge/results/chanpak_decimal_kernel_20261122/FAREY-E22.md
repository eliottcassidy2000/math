# A 37-vertex Farey flip at the twenty-second Euclid center

**FINITE-EXACT.** Two strongly connected tournaments on 37 vertices, differing only at the arc joining vertices 5 and 6, have Hamilton path counts
\[
                  8,796,090,925,055,\qquad 8,796,090,925,057.
\]
Their midpoint is \(E_{22}=2^{21}(2^{22}-1)\). Together with the preceding certificates, every integer exponent from 3 through **23** now has an exact single-arc crossing. Exponent 22 is composite. At 23 the Mersenne factor is composite, and an independent marked gate already supplies the upper and lower pair. [Earlier consecutive band](../chanpak_e61_fibonacci_target_20261122/ADDITIONAL-SEAMS.md) · [Exponent-23 marked gate](../chanpak_e31_physical_20261122/README.md)

Order three tagged blocks \(A,D,V\) by ordinal sum, then introduce a common apex. As in the earlier star construction, outgoing row masks use zero-based vertex indices; the apex beats the vertices in the tag mask. The exterior left and right addresses are \(x=0\) and \(y=185591\). In rank order the transitive right flank uses the complemented binary digits of \(y\) read from low to high.

| Block | Outgoing rows | Tags |
|---|---|---:|
| \(A\) | 2, 4, 1 | 1 |
| \(D\) | 4, 17, 34, 7, 13, 27 | 22 |
| \(V\) | 310, 364, 56, 465, 482, 456, 389, 7, 132 | 420 |

Reverse \((2,3)\) inside \(D\) to obtain the upper configuration. Its first four chart entries \(h,q,c,t\) agree across the flip, and \((\Delta k,\Delta e,\Delta r)=(-2,2,2)\). The left triangle has three full paths, one ending at a tag. The right module has 2053 full paths, 1369 starting at a nontag. The normalized change is the determinant
\[
                   3\cdot1369-2\cdot2053=1 .
\]
Consequently the signed change stays two at *every* pair of exterior addresses. This is a second absolute placement with the Farey method, with a different pair of rational proportions from those at exponent 15. [Marked join and selector](E15-FAREY.md)

In the order \((h,q,c,t,e,r,k)\), the lower profile is
\[
 (178611,132280320,47383164,61153374,99876,119074,22266611).
\]
The upper one increases only the last entry by two. The integer \(y\) has 18 binary digits, yielding \(3+6+9+18+1=37\) vertices. Substitution in \(k+cy+e[y\ {\rm odd}]+hf(y)\) gives the announced lower. The replay also counts the two full tournaments independently of the joined profile, splitting at the apex and enumerating paths on the occupied subsets of each block. It checks both strong-connectivity predicates and the literal single edit. [Exact certificate](clock_obstruction.py)

## The same fixed edge can match every Euclid center at every finite modulus

**PROVED.** For each \(M\ge1\) and \(p\ge2\), some pair of exterior addresses on these same tagged blocks gives the two counts \(E_p-1,E_p+1\) simultaneously modulo \(M\). More generally the permitted lower phases are all odd targets \(N\), with the additional condition \(N\equiv0,2\bmod3\) only if \(3\mid M\).

In fact \(\gcd(h,q,c,t)=3\), \(v_3(h)=1\), and the lower height satisfies \(L\equiv2+[x\ {\rm odd}]\bmod3\) identically. For a prescribed allowed target choose the parity of \(x\) accordingly and take \(y\) even. At all other prime powers dividing \(\gcd(h,M)\) at least one of \(q,c,t\) is a unit. If \(q\) is a unit, fix \(x\) so \(c+qx\) is a unit and solve for \(y\) after subtracting the now fixed parity terms. If \(q\) is a nonunit, one of the two linear coefficients is a unit; set the other variable to zero modulo this prime power and solve. Combine the residues by CRT, preserving the selected parities. The finite reset automaton then controls \(f(y)\) at the unchanged address residue to solve the rest of the congruence. All \(E_p-1\) pass since their residues modulo three are zero at odd \(p\) and two at even \(p\). [Reset proof and contrasting incompatible clocks](UNSYNCHRONIZED-CLOCKS.md)

For example, take the product \(M=1,618,964,990,108,856,390\) of the 15 absolute supersingular primes. Start with \(x=1,y=36\) and append 3,332,926,730,265,235 copies of \(1^{18489240}0^{18489240}\) to \(y\). The lower seed height is 6,557,140,035. Each copy adds \(h\bmod(M/2)\) at fixed address residue, so the constructed lower and upper counts match \(E_{61}-1,E_{61}+1\bmod M\). The absolute exponent-61 problem in this chart remains open. Its subfamily in which each exterior word has at most one run of ones, with arbitrary trailing zeros, has no exact hit at any order.

The concurrent quadratic-length word controller reduces this certificate to **one 7,562-bit right address** with the same one-bit left address. The replay constructs it and evaluates the exact integer count before reducing modulo the 15-prime product. [Short finite-state inverse](LOCAL-WORDS.md) · [Physical replay](clock_obstruction.py)
