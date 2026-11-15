# Complete dyadic classification for cyclic diagonal covers, in every grade

This extends the rank-three Gram theorem (CAYLEY.md). Write
\(F_k(A)=q_k(C(A,A))/k\) for \(k\ge1\), with \(|A|\ge2\). For every
\(L\ge4\), the literal tournaments attain **precisely** these residues:

| Grade \(k\) | Residues of \(F_k\) modulo \(2^L\) |
|---|---|
| odd | all residues \(3\bmod8\) |
| \(k\equiv2\bmod4\) | all residues \(2\bmod8\) |
| \(k\equiv0\bmod4\) | all residues \(2\bmod16\) |

These are physical local results, not absolute coverage. In particular the change
from 8 to 16 is predicted by the exact reciprocity of slot polynomials, not by
the number of coordinates alone.

## Reciprocity supplies the extra frozen bit

The slot function \(q_j(A)\) is an integer-coefficient polynomial in \(j\),
because \(P_i/i!\) is the number of *unordered* covers by \(i\) nonempty
paths and \(q_j=\sum_i(P_i/i!)j(j-1)\cdots(j-i+1)\).
It satisfies \(q_{-j}=(-1)^{|A|}q_j\). Indeed, concatenate the ordered slots
to a permutation and let \(d\) be the number of its backward adjacencies.
Slot labels form a weakly increasing sequence with a strict rise at every
backward adjacency, contributing \(\binom{j+n-1-d}{n}\) for \(n=|A|\).
At \(-j\) the polynomial value is \((-1)^n\binom{j+d}{n}\), which is the
original contribution of the reverse permutation: reversing turns all good
adjacencies bad and vice versa because \(A\) is a tournament. Summing proves it. This is also the principal specialization of the converse/duality identity of [Grinberg–Stanley, Theorem 1.39](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge.pdf), as discussed in BOTT.md.

For odd k, the Gram form forces 3 mod 8 as previously shown. For even k,
the middle Gram coordinate \(q_k\) is divisible by 4, and its two neighboring
grades are odd. Expanding the squares gives
\[
 F_k(A)\equiv q_{k-1}(A)^2+q_{k+1}(A)^2\pmod{16};
\]
the cross term contains \(2q_k(q_{k-1}+q_{k+1})\), divisible by 16.
If 4 divides k, then \(k+1\equiv -(k-1)\pmod8\), so reciprocity gives
the same odd square modulo 16 from both grades. Their sum is 2 mod 16.

## Saturation by literal ordinal sums

\(q_j\) multiplies under ordinal sums. A singleton contributes \(j\);
the directed triangle contributes \(j^3+2j\) (its nonempty counts are
\(3,6,6\)). An ordinal run of \(m\ge L\) singletons kills any even-grade
\(q_j\) modulo \(2^L\), and \(q_0=0\) already.

**Odd k:** choose the auxiliary block \(D=TT_1\) if \(k\equiv\pm3\bmod8\),
and \(D=C_3\) otherwise. Then \(q_k(D)\equiv\pm3\bmod8\) and its square
generates all of \(1+8\mathbb Z\pmod{2^L}\): beginning at valuation 3, repeated
squaring increases valuation by one at every step. Take the ordinal sum of
\(t\) copies of \(D\) and \(m\ge L\) singletons. The resulting normalized norm is
\(3 k^{2m} q_k(D)^{2t}\bmod2^L\), which fills the entire required sheet.

**\(k\equiv2\bmod4\):** write \(j_-=k-1,j_+=k+1\). Let \(G=1+8\mathbb Z\bmod2^L\).
The two elements of \(G^2\) contributed to the *squared* grades by a
singleton and a triangle are respectively
\[
 (j_-^2,j_+^2),\qquad
 ((j_-^3+2j_-)^2,(j_+^3+2j_+)^2).
\]
Via base-9 logarithms modulo the order of \(G\), these two columns form a
matrix with odd determinant: modulo 2 the singleton has exactly one active
coordinate (exactly one of \(j_-,j_+\) is \(\pm3\bmod8\)), and the triangle
toggles both (its odd-grade profile is \(3j\bmod8\)). Therefore the two
columns generate all of \(G^2\); finite positive repetitions realize any
target pair. Start with \(m\ge L\) singletons so the even grade is zero; this
only multiplies the odd-grade squared pair by a fixed element of \(G^2\).
Now \(F_k\) reduces to the sum of any two elements of \(G\), exactly the
desired class.

**\(k\equiv0\bmod4\):** use \(m\ge L\) singletons and \(t\) triangles.
Writing \(a=j_-^{2m}, b=j_+^{2m}, r=(j_-^3+2j_-)^2,
s=(j_+^3+2j_+)^2\), the output modulo \(2^L\) is
\(f(t)=a r^t+b s^t\). Here \(a,b,r,s\equiv1\bmod8\). For
\(k=0,4,8,12\bmod16\), the pairs \((r,s)\bmod32\) are respectively
\((9,9),(1,17),(25,25),(17,1)\), so \(r+s-2\equiv16\bmod32\).
By squaring induction for any \(z\equiv1\bmod8\),
\(z^{2^\ell}-1\equiv 2^\ell(z-1)\pmod{2^{\ell+5}}\). Hence
\[
f(t+2^\ell)-f(t)\equiv2^{\ell+4}\pmod{2^{\ell+5}}.
\]
(The factors \(a r^t,b s^t\equiv1\bmod8\) have no effect at this precision.)
Choose the binary bits of \(t\) from the lowest upward. The \(\ell\)-th bit
changes exactly the next desired residue bit without disturbing earlier
ones. This fills every residue \(2\bmod16\), using nonnegative \(t\) throughout.

This is a second instance of the same design principle as the complementary
binary clock: isolate a move whose first nonzero bit is in the right
position, and prove that every such move is a literal composition. A new
frozen bit is compatible with complete local control *on the smaller sheet*.
No claim is made for arbitrary integer hitting, high-density physical norms,
or reduction of Chan–Pak to a finite congruence problem.
