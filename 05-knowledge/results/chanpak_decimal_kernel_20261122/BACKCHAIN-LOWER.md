# An exact exponent-61 lower witness with a one-arc exit from divisibility by three

**PROVED with finite-exact contact data.** There is a strong 133-vertex tournament with exactly \(E_{61}-1\) Hamilton paths and a single reversal to a strong tournament of height \(E_{61}+1\bmod 69{,}654\), where \(69{,}654=2\cdot3\cdot13\cdot19\cdot47\). Thus the lower endpoint is exact and this edge escapes the mod-three obstruction found in the compact 93-vertex families. Another order-133 lower has a strong neighbor of upper phase modulo \(108{,}914=2\cdot13\cdot59\cdot71\). The jumps are larger than two in both cases. A full neighborhood calculation excludes an exact upper exit in all 10,080 reorderings and reversals of these seven factor words. [All single-arc neighbors](ACTIVE-ENDPOINTS.md) Each construction extends to infinitely many formal Euclid exponents with the same local information and an exact lower count. [Earlier lower grammars and their marked obstructions](../chanpak_e61_fibonacci_target_20261122/NEW-LOWER.md) · [The mod-three plateau](../chanpak_e61_fibonacci_target_20261122/PROFILE-FIBERS.md)

The companion pressure theorem bounds the full transitive-gap backchain image through \(X\) by \(O(X^{23/24})\). Exact hits on the lower Euclid sequence can therefore coexist with a sparse integer image. Its terminal inverse also allows zero runs inside a factor word. [Pressure bound and terminal inverses](../chanpak_terminal_profiles_20261122/README.md)

### Combining positive resets at an actual lower height

Use the concurrent [three-register backedge chain](../chanpak_absolute_descent_20261122/BACKCHAIN.md). A width word lists the numbers of ordinary vertices between consecutive boundaries. Reverse the arc between each consecutive pair of boundaries of an otherwise transitive order. At a positive gap of width \(m\) following a prefix of height \(F\) that ends with a positive gap, the register state resets to \((2^m F,0,F)\). This is exactly \(F\) times the first-gap state for an independent word. Thus concatenating positive-start, positive-end words multiplies their heights **inside a single strong tournament**. For words with at least two gaps, the first boundary reaches the last, which reaches every earlier boundary along the backwards chain. Every ordinary vertex lies on a forward route from the first boundary to the last. This also proves strongness after the two particular contact reversals below, which leave the boundary chain and those routes present.

The following words have heights given by the positive transfer. The last column is the order increase when appending the word to an existing boundary.

| Word | Height | Order increase |
| --- | ---: | ---: |
| \((2)\) | 5 | 3 |
| \((1,0,0,1)\) | 31 | 6 |
| \((1,0,0,0,1,0,0,4,0,0,0,2)\) | 119,743 | 20 |
| \((1,0,0,0,0,1,0,0,0,2,0,0,0,1,0,1)\) | 379,127 | 22 |

Use the first word twice, then the second, then the fourth, then the **reverse** of the third. Reversing a word takes the chain to its converse up to relabelling, preserving height. The resulting prefix has \(F=35{,}183{,}298{,}379{,}775\), and its height just before the last gap is \(G=10{,}504{,}481{,}267{,}175\). Append the two positive widths 15 and 61 in either order. The height and order are exactly
\[
 F(2^{15}+1)(2^{61}+1)=(2^{60}-1)(2^{61}+1)=E_{61}-1,
 \qquad 1+6+6+20+22+16+62=133 .
\]

### A target-selecting boundary register

Let \(b\) be the boundary just before those last two gaps, \(A\) the vertices before it, and \(B\) the vertices after it. All arcs run from \(A\) to \(B\). Read \(b\) as the common apex. Its only tag in \(A\) is the preceding boundary, followed in \(A\) by a nonempty string of ordinal singletons. Hence \(h_A=G,e_A=0,k_A=F\). We also have \(c_A=F-G\). Indeed append just one positive test gap of width \(s\): its deletion block has height one, its apex height is \(2^s+1\), and its terminal unsigned count is \(2^s\). The join law gives \(F+2^s(G+c_A)=(1+2^s)F\). For a contact toggle between \(b\) and one vertex in \(B\), \(h_B\) stays fixed and the exact jump is therefore
\[
                         \Delta H=G\,\Delta k_B+(F-G)\,\Delta t_B.
\]
Whenever \(\Delta k_B\ne\Delta t_B\), the contact selects the *only* possible register \(G=(2-F\Delta t_B)/(\Delta k_B-\Delta t_B)\) for a +2 exit from this lower. This argument holds for every positive-ending prefix with that height. [Shared-apex join](E15-FAREY.md)

Index vertices of the suffix starting at zero at \(b\). The following physical changes give the claims above.

| Widths | Contact reversed | Exact increase from \(E_{61}-1\) | Modulus at which the increase is 2 |
| --- | --- | ---: | ---: |
| 15, 61 | \((0,53)\) | 5,806,327,764,816,848,693,124,291,406,459,391,128,371,200 | 69,654 |
| 61, 15 | \((0,56)\) | 41,538,374,868,278,609,977,132,081,004,954,400 | 108,914 |

To evaluate \(k_B\), order the three boundary vertices on the path and distribute the other vertices, in their transitive order, among the four gaps. Each gap remembers whether it is empty, or whether its current last vertex can exit to the next boundary. To evaluate \(t_B\) replace \(b\) by a separator with two-way arcs and require the path to finish at a nontag in \(B\). Splitting there gives the required ordered two-slot cover after exchanging the two slots. The replay compares both evaluations to subset path counts on small actual row matrices. Run `python3 -B backchain_lower.py`. [Physical construction and contact replay](backchain_lower.py)

### Every two-positive-gap ending at this exact target

**FINITE-EXACT.** Allow arbitrary positive terminal widths \(s,t\) and *any* positive-ending prefix attaining the required factor \(F\). No contact from the shared boundary into the suffix gives an exact +2 at \(E_{61}-1\). Divisibility by \((1+2^s)(1+2^t)\) leaves widths drawn from 1, 2, 3, 5, 6, 10, 15, 30, 61. There are 60 ordered pairs whose product divides the lower target. Over all 1,906 resulting contacts, the displayed register inverse has **no integer solution even in the relaxed interval** \(0<G\le F\). When the two differences agree, the condition reduces to \(F\Delta t_B=2\), which also fails throughout. Run `python3 -B backchain_lower.py --all-widths`.

The same coalescence formula for a general arc edit supported in \(B\cup\{b\}\) reads \(\Delta H=G\Delta k_B+(F-G)\Delta(h_B+t_B)\). Therefore every such jump is divisible by \(\gcd(F,G)\), which is 293,823,425 in the two displayed 133-vertex hosts. Edits supported entirely in \(A\cup\{b\}\) also miss +2 in those two hosts: their coefficients in the join formula are \(h_B,k_B-h_B,t_B,r_B\), all zero modulo three. Both odd widths can be reduced to one modulo this prime by the width law below. The suffix profile then has \((h_B,k_B,t_B,r_B)=(3,9,6,0)\). Only a strictly crossing contact between \(A\) and \(B\) could still supply a +2 exit before evaluating the full flip neighborhood.

### Transporting the exact lower and its marked upper phase

Put \(p_j=1+60\cdot2^j\). Prepend the single-gap factors with widths \(60,120,\ldots,60\cdot2^{j-1}\) ahead of the original prefix, and replace its terminal width 61 by \(p_j\) in the same position. The prefix registers \(F,G\) are multiplied by
\[
       P_j=\prod_{i<j}(1+2^{60\cdot2^i})=(2^{p_j-1}-1)/(2^{60}-1).
\]
The complete lower height becomes \((2^{p_j-1}-1)(2^{p_j}+1)=E_{p_j}-1\), exactly for every \(j\ge0\).

For each prime \(q\) in the modulus of a row of the table, the marked suffix readings agree modulo \(q\) whenever \(p_j-61\) is a multiple of \(q-1\). Split the marked growing transitive gap at its fixed contact rank if needed. The only changing segment is a positive homogeneous tail. A skeleton visiting it \(r\) times has filling count \(r!{w\brace r}\), periodic in the positive width \(w\bmod(q-1)\) modulo \(q\) by the alternating power formula. This applies to the separator reading as well. [Transitive-width residue law](../chanpak_e61_fibonacci_target_20261122/SEED-LIFTS.md)

For the second row, take any multiple \(j\) of 84. For \(q=13\), \(P_j\equiv2^j\equiv1\pmod{13}\). For 59 and 71, the order of \(2^{60}\) divides 29 and 7 respectively, and \(2^j=1\) modulo both. If \(2^{60}\ne1\bmod q\), the geometric product gives \(P_j=1\bmod q\), and direct evaluation gives \(2^{60}\ne1\) at both primes. Also \(p_j=61\bmod(q-1)\) at all three. Parity is automatic. The first row works for every multiple \(j\) of 132: the same argument uses \(2^{60}=1\bmod 3,13\), and orders of \(2^{60}\) dividing 3 and 23 at 19 and 47, where it is nontrivial. Hence both the suffix differences and the prefix multipliers return, so the same physical contacts realize the announced upper phases with exact lower endpoints. No primality assumption is made on \(2^{p_j}-1\).
