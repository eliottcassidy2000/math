# The five-corridor clock passes every finite congruence test simultaneously

Use only gap words with a positive initial width, then any positive
number of complete repetitions of the ordered zero corridors
\[
                            (4,2,5,3,7),
\]
with positive widths between them. No factorizing corridor of length
0 or 13 is allowed. **Every odd residue modulo every integer \(m\)
has a strong witness in this single family.** More precisely its
closure among the profinite integers is the integers odd at 2.

This is a stronger local/global separation than separate surjectivity
for powers of two and for odd moduli: the same physical tournament
must satisfy all of a given finite list of congruences at once.
It is not an absolute-target construction. The whole first terrace
has \(O(X^{0.952})\) distinct outputs by the published pressure
bound; this particular subfamily inherits that obstruction.

## 1. What we import from the odd clock

Conventions are column-state: a step of length \(q\) and width
\(y=2^r\) is \(\operatorname{diag}(1,y)C_q\).
Read count by \((1,1)\) on the left, and start from \((1,y_0)^t\).
Let \(A=C_7C_3C_5C_2C_4\), i.e. all widths equal 1
as a *modular* choice, physically realizable by large positive
exponents. The three standard variants are \(A\), doubling its
last width, and doubling its penultimate width.

Write the modulus as \(2^k b g\) with \(b=3^a7^d\),
\(\gcd(g,42)=1\). Enlarge it if necessary so that \(k\ge1,a\ge2,d\ge1\).
On the good part \(g\) the standard variants generate a group
containing \(\mathrm{SL}_2(\mathbb Z/g)\), by the odd-modulus proof
in this packet. On \(b\) a suitable power of \(A\) is the unique
stable rank-one idempotent \(E_b\); each sandwiched standard
variant acts as a unit on its image. Its reduction modulo 63 is
\[
                         E_{63}=\begin{pmatrix}12&51\\53&52\end{pmatrix}.
\]
That proof's other input is a *terminal* surjectivity statement:
given \(E_b\), the entries of the seed \((1,y_0)^t\) and the five
widths of one terminal clock can be chosen so that the scalar
readout of that terminal clock after \(E_b\) is any residue modulo
\(b\). The exhaustive statement is only modulo 63; the affine
unit derivatives and the width lifts \(y\mapsto64^t y\) prove
the remaining prime powers. All six exponents can afterward be
made at least \(k\), without changing any chosen odd residue,
by adding a multiple of the order of 2 modulo \(bg\).

## 2. Build the binary digits without breaking the odd repair

Choose the common baseline width exponent \(L\ge k\) a multiple
of the order of 2 modulo \(bg\). Thus \(y=2^L\) is 1 on the
whole odd part and 0 modulo \(2^k\). Use this same physical
baseline and the two doubled-width variants throughout.
Put \(\gamma=e_4e_2e_5e_3e_7\), an odd integer.
Modulo \(2^k\) the baseline matrix has the form
\[
 A_2=\gamma\begin{pmatrix}1&f_4/e_4\\0&0\end{pmatrix}.
\]
A sufficiently chosen positive power of the *physical* baseline
word now gives a sandwich \(E\) that acts as \(E_b\) on \(b\),
as the identity on \(g\), and as the displayed rank-one matrix
without \(\gamma\) on \(2^k\). This is just simultaneous periodic
powers of the three finite actions. In particular \(E_2\) projects
onto the first coordinate axis and acts there as the identity.

One binary adjustment is a sandwich \(E A(y) E\), with \(A(y)\)
a full clock whose *first* width is \(y\), and whose four other
widths are the large baseline. On the first coordinate axis
modulo \(2^k\) this multiplies by
\[
                           \gamma(1+(11/9)y).
\]
Indeed only the path through the first two old transfer states
remains uncollapsed for one step:
\(e_2e_4+f_2f_4y=27+33y\); every subsequent step contributes
its \(e_q\). The denominator 9 is a unit modulo \(2^k\).

For bit \(r=1,\ldots,k-1\) choose \(y=2^r\) or the large baseline.
After stripping a factor \(\gamma\), these act respectively by
\(1+(11/9)2^r\) and by 1. The former toggles bit \(r\) of an odd
number while retaining all lower bits. Use exactly one sandwich
per bit, even when not toggling. When both the physical seed
exponent and the five terminal exponents are at least \(k\),
the initial seed projects to \((1,0)^t\) in this part and the
terminal clock contributes one more \(\gamma\). Thus the binary
target can be programmed by aiming at the target divided by
\(\gamma^k\) in the bit toggles.

These choices do not force a zero or nonunit on the bad odd line.
For all six residues \(y=2^r\bmod21\) (including the large baseline
choice \(y=1\)), the compressed matrix \(E_b A(y)E_b\) is nonzero
modulo both 3 and 7. The exact replay checks the six matrices.
Therefore every chosen sandwich multiplies the image of \(E_b\)
by a known unit. Their aggregate \(\lambda\) can be compensated
by choosing the terminal bad target to be the original one
divided by \(\lambda\). That choice changes no binary bit because
all the terminal and seed exponents are then shifted large.

## 3. Good correction at the end

The three standard sandwich generators lie in the local matrix
group with identity \(E\): their restrictions to each projected
bad line are units, and the good part is invertible. They generate
a finite group by taking positive powers (the subgroup identity
is realized by \(E\)). Its bad projection is abelian, being a
product of the scalar unit groups on the three bad lines. Its
good projection is the group generated by the original standard
variants, as \(E_g=I\). Hence its commutator subgroup fixes both
the binary and bad odd lines while surjecting onto a group
containing \(\mathrm{SL}_2(\mathbb Z/g)\), which is perfect for
the prime powers in question. Choose a physical word realizing
the needed commutator image before the terminal clock.

Explicitly, after the chosen seed and binary sandwiches the good
column is primitive: the seed is primitive and all these matrices
are invertible on \(g\). The left readout row after the terminal
clock is primitive there too. An \(\mathrm{SL}_2\) action can send
the good column to any primitive column, and a primitive row has
any desired scalar value on some primitive column. Both claims
follow by completing the row or column to a determinant-one
matrix over \(\mathbb Z/g\). The good correction thus sets the
remaining scalar residue without touching either of the other
parts. Every operation just concatenates whole copies of the
same five-corridor clock. No subtractive or rational weights become
physical arcs.

This computation exhibits why local saturation is inadequate
even when imposed *simultaneously at every finite modulus*.
For any individual target \(N\), there are sequences in this
sparse physical family converging to \(N\) modulo every modulus.
The witnesses may escape in ordinary size and count. What is
missing remains an exact absolute-height descent, not a more
refined congruence completion.
