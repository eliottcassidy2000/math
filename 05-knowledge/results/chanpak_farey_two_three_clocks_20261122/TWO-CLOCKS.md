# Farey modular pruning and its absolute limit

For the canonical lower Farey profile write

\[
  H(x,y)=k+t x+c y+qxy+r[x\text{ odd}]+e[y\text{ odd}]
          +h(f(x)+f(y)),
\]
where \(f(2v+1)=f(v)\), \(f(2v)=f(v)+v+[v\text{ odd}]\), and
\(f(0)=0\). In particular \(f\) is even-valued. The upper orientation
has height \(H(x,y)+2\) at the same addresses.

## Any modular system has a solution on every relevant single axis

The receiver recurrence gives, for each \(n=4,8,\ldots,112\),
\(\gcd(h,t)=1\). The small replay below is an exact check using the
existing finite-run receiver machine. These are all potentially relevant
spines: every core has \(h=145(2^n+1)\), and its ordinary apex height
\(k\geq h\) by the odd block-gluing formula (or extend each core
Hamilton path through an apex insertion, as guaranteed in a tournament).
For \(n\geq116\), this already exceeds \(E_{61}-1\).

In fact for **every odd** \(N\), every \(M\geq1\), and any one of these
spines, there is an even positive \(x\) such that
\(H(x,0)\equiv N\pmod M\). Replace \(M\) by \(2^sD\) with \(D\) odd
and \(s\geq1\), doubling if necessary. Put \(g=\gcd(h,D)\). Choose
positive \(x_0\equiv0\pmod{2^s}\) and
\(t x_0\equiv N-k\pmod g\), possible by the gcd check. It follows
that
\[
 h z\equiv N-k-tx_0\pmod M
\]
has a solution with \(z\) even: the congruence modulo \(2^s\) has an
even solution because \(N,k\) are odd, and the odd congruence is
solvable by construction.

Append a block \(1^a0^b\) to the binary address, where \(a,b\geq s\)
are positive multiples of \(\operatorname{ord}_D(2)\) (use 1 if \(D=1\)).
The new address is \(x'=2^b(2^a(x+1)-1)\), so \(x'\equiv x\pmod M\).
Writing \(w=2^a(x+1)-1\), the increase in the read is
\[
 \Delta f=(2^b-1)w+1,
\]
because the ones preserve the read and the first following zero is the
only one with an odd parent. Modulo \(2^s\) this is 2, and modulo \(D\)
it is 1. Its additive order is therefore \(2^{s-1}D\), and its
iterations step through the entire even subgroup mod \(M\). Pick the
iteration whose read is \(z\).

This proves simultaneous modular solvability by a single address for
every finite system, but is not a hit at the exact target. In
particular it explains why unrestricted modular exclusion cannot prove
the E61 target missing from any canonical spine, even with the right
address pinned to zero. Size restrictions can still make modular
pruning effective.

## Exact joint mod 3 and 2-adic pruning on an affine line

The read cannot be reduced by dropping the high address bits. Instead,
write its binary digits \(u_j\) in little-endian order. Rearranging
its defining sum shows that when reading bit \(u_j\), its exact
contribution is
\[
 u_j(W_j+L_j),\qquad
 W_j=\sum_{i<j:u_i=0}2^{j-1-i},\quad
 L_j=[j>0,\ u_{j-1}=0].
\]
Advance by \((W,L)\mapsto(2W+1-u_j,1-u_j)\), reducing \(W\)
and the partial read modulo the desired modulus. This gives an
LSB-first weighted automaton; modulo four it is also immediate that
\(f/2\) mod 2 counts the overlapping binary triples 010 and 100.

For an affine line \((x,y)=(A+Bj,D+Ej)\), read the bits \(b\) of
\(j\) in the same order and use signed carries. Initialize the
two carries at \((A,D)\); divide \((A_{\rm carry}+Bb)\) and
\((D_{\rm carry}+Eb)\) by 2 with nonnegative remainders, obtaining
the next carries and physical address bits. The filter tracks both
addresses modulo \(m\), their two weighted automata, and the combined
read. Track interval endpoints with signed carries initialized at
\((-L,U)\), updating \(z\mapsto\lfloor(z+b)/2\rfloor\) and
\(w\mapsto\lfloor(w-b)/2\rfloor\), respectively. After enough
bits these endpoint carries are zero precisely when \(L\leq j\leq U\).
Exact arithmetic at terminal states gives an exact count modulo
\(m\), potentially clearing an entire line slab when zero. The
state space is bounded for fixed slopes and modulus independently of
interval length: at bit depth \(i\), an address carry is
\(\lfloor(A+B p)/2^i\rfloor\) with \(0\le p<2^i\), hence has at
most \(|B|+2\) possible values (similarly for \(E\)), and the
endpoint carries each have at most two values.

The standalone script compares 150 signed-affine cases with brute
force. On a 4,096-position example, the survivor counts are
1,399 mod 3; 1,988 mod 4; 675 mod 12; 366 mod 24; and 205 mod 48.
These are modular pruning counts, not a new absolute E61 exclusion.

## Closed form for the left slope

For every \(n=4m\), the tag word is \((1001)^{m-1}1010\):
its binary integer is \((3\cdot16^m+2)/5\), and the difference
at consecutive \(m\) is precisely the new leading block \(9\cdot16^m\).
The left slope of the full lower chart is
\(t=10240 T_R+4592(2^n+1)\), where \(T_R\) is the receiver's
untagged-end two-slot total. The exact identity is
\[
 1170T_R=1000\cdot 3^n-441\cdot 2^n-1404.
\]
Indeed for each of the six boundary-marker orders, the receiver machine
for \(T_R\) has \(3^4=81\) states. The prefix block \(1001\) supplies
one repeated linear operator; the final block \(1010\) and terminal test
supply a fixed output functional. The difference from the displayed sum
is a scalar read of the direct sum of the 486-dimensional operator and
\(\operatorname{diag}(81,16,1)\). By Cayley–Hamilton, checking the
first 489 terms proves the identity for all positive \(m\).
The axis-unit replay performs all checks exactly.

## The first failed axis and a coupled 2–3 lift

**PROVED, with exact finite first-face calculation.** Put \(S=2^n+1\)
for \(n>0\) divisible by four. The closed form gives the exact identity
\[
 \gcd(t,h)=\gcd\bigl(S,1000\cdot3^n-963\bigr).
\]
Indeed \(S\) is coprime to \(2,3,5,13\): modulo \(13\), \(2^4\) has
order three and cannot produce \(-1\). Thus the denominators and the
factor 10240 are units modulo \(S\); substitute \(2^n=-1\) in the
closed form. The additional factors of \(h=145S\) make no extra
contribution: \(t=4\bmod5\), and for \(n/4=1,2,\ldots,7\bmod7\) its
residues modulo 29 are \(8,18,15,17,10,16,13\). The last check uses the
orders of \(2^4\) and \(3^4\), which both divide seven modulo 29.

Consequently an odd prime \(p\) dividing the chart's axis gcd is
exactly a simultaneous root of
\[
 F_2(n)=2^n+1,\qquad F_3(n)=1000\cdot3^n-963\quad\bmod p.
\]
The first failing spine is \(n=8200\), with gcd \(19681\) (prime).
This concerns the whole axis, not just one fixed parity: when the gcd
has an odd prime factor the two parity branches occupy at most two
classes modulo that prime; they cannot cover every odd target.
At \(p=19681\), \(2\) and \(3\) have orders 3280 and 9840,
respectively. Both equations hold exactly at \(n=8200\bmod9840\).

Take a shared residue class \(n=n_0+Lj\), with
\(L=\operatorname{lcm}(4,\operatorname{ord}_p(2),
\operatorname{ord}_p(3))\). In the first lift from \(p\) to \(p^2\), the
two normalized readings are affine in \(j\bmod p\):
\[
 \frac{F_i(n_0+Lj)}p\equiv A_i+B_i j,\qquad
 A_i=F_i(n_0)/p,\quad
 B_i=d_i\,i^{n_0}(i^L-1)/p,
 \qquad d_2=1,\ d_3=1000.
\]
The congruences are taken after the indicated integral division.
When at least one slope is a unit, simultaneous lifting is equivalent
to \(A_2B_3-A_3B_2=0\) in \(\mathbb F_p\). If both slopes vanish, it
instead requires both offsets to vanish. This is two conditions on
one exponent, so two separately liftable units need not lift jointly.
Here \((A_2,B_2)=(1745,2094)\) and \((A_3,B_3)=(6270,7167)\). The
individual choices of \(j\) are \(16400\) and \(1828\), respectively;
**19681 never divides any axis gcd to depth two**. These are exact
calculations, not a prediction about other primes.

The slopes are scaled Fermat quotients:
\((a^L-1)/p\equiv L(p-1)^{-1}q_p(a)\bmod p\).
This is the actual Farey instance of the [Mirimanoff lift
warning](../chanpak_mirimanoff_unit_lifts_20261122/MIRROR-UNIT.md):
a base-3 slope is not the whole obstruction. Here the binary and
ternary linear equations share their exponent and must agree. The
first failure is far above the potentially relevant exponent-61 spine
lengths, and does not eliminate a two-word chart.

The first axis failure does **not** force a failure of the full chart:
at \(n=8200\), \(q=16538\), \(c=18124\), and \(r=0\) modulo
19681. The left axis lands in one residue class even allowing both
parities, whereas the right coordinate restores the missing phase. In
fact \(\gcd(h,t,c,2q)=1\), so the [complete local image theorem](../chanpak_decimal_kernel_20261122/LOCAL-WORDS.md) gives all odd
classes modulo any modulus even in the even-even two-coordinate
branch. This is a concrete instance of a second port repairing a
local obstruction without promising an exact absolute target.
