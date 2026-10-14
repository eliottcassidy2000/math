# Prime shapes, the predecessor of 561, and marked Pythagorean fibres

2026-10-14. **PROVED:** the domain comparisons and explicit maps below.
The divisor and primitive-triangle classifications are **INHERITED**,
with no priority claim. **FINITE-EXACT:** the declared source checks and
exported triangle receipts. No prime-density, universal Collatz, or BSD
conclusion follows from these comparisons.

## 1. Recover the definitions before changing the exponent

For $N\ge2$, the repository's $F=U+S$ means

\[
F(N)=\#\{d\mid N:1<d<N\},\quad
S(N)=\#\{d\mid N:1<d<N,\ d\text{ squarefree}\},
\]
\[
U(N)=\#\{p\mid N:p\text{ prime},\ p<N\}.
\]

This is the exact convention in
[arithmetic braids, DB1--DB3](arithmetic_braids_20260917_divisors.md).
$U$ is not the number of factors counted with multiplicity. Its underlying
set is contained in the squarefree-divisor set, so the equality of counts
is not a disjoint-set partition. At a prime all three counts are zero;
$N=1$ is excluded.

The closest proved mechanism is the exponent box of divisors and its
Boolean squarefree subcube. The hostile is $p^2$, where deleting the
repeated-factor direction would hide the failure $1\ne1+1$.
The corrected near miss is the universal conservation reading already
rejected in the [odd-square synthesis](odd_square_20260921_synthesis.md)
and its 2026-09-21 correction. The least-used marker is the **freeness
threshold**: squarefree and fourth-power-free are different predicates.
Our board is exponent / freeness threshold / predecessor / source role /
Gaussian content / tree address.

For distinct primes $p_i$ and $N=\prod_{i=1}^r p_i^{a_i}$,

\[
F=\prod_i(a_i+1)-2,\quad
S=2^r-1-[N\text{ squarefree}],\quad
U=r-[N\text{ prime}].                                    \tag{1}
\]

The inherited complete classification is
$F=S+U$ iff $N=p$, $p^3$, or $p^2qr$, with distinct labels in the last
shape. This packet does not rediscover that classification.

For the proposed marked family $N_a=p^aqr$, $a\ge2$, (1) gives

\[
(F,S,U)=(4a+2,7,3),\qquad F-S-U=4(a-2).                 \tag{2}
\]

Therefore $p^4qr$ has counts $(18,7,3)$ and defect eight, not zero.
The least example is $240=2^4\cdot3\cdot5$; $60=2^2\cdot3\cdot5$
is the positive comparison. Multiplication by $p^2$ preserves prime
support and changes the defect by eight on this marked family. It does
not preserve divisor balance.

There is a valid changed-domain statement. Let $S_k$ count proper
nontrivial divisors free of every prime $k$th power. Then

\[
S_k=\prod_i\min(a_i+1,k)-1-[N\text{ is }k\text{-free}].
\]

At $N_a$, $S_a=4a-1$, so $F=S_a+U$. In particular
$18=15+3$ at $p^4qr$ uses **fourth-power-free** divisors. This is already
[divisor balance family, Theorem 3, §3.2](collatz_mod6_20260917_divisor_balance_family.md):
for $k\ge3$ the complete shapes are $p,p^{k+1},p^kq^2,p^kqr$.
Keeping $S=S_2$ makes the proposed extension false; changing it explicitly
produces a true inherited extension. The rank layers themselves remain
the inherited $3,4,\ldots,4,3$ profile from
[divisor-layer geometry, §3](divisor_layers_pair_geometry_20261004.md).

## 2. The actual source of the primes 2, 5 and 7

The exact predecessor calculation is

\[
561=3\cdot11\cdot17,\qquad
560=561-1=2^4\cdot5\cdot7,
\]
\[
L(561)=\operatorname{lcm}(3-1,11-1,17-1)
       =\operatorname{lcm}(2,10,16)=80,\qquad 560=7L(561). \tag{3}
\]

Thus $\{2,5,7\}$ is the prime support of the predecessor **560**.
Only $2,5$ occur in the prime-predecessor lcm; the additional factor 7
is in its integer quotient. These are distinct objects.
For general $N$, the code calls $\operatorname{lcm}_{p\mid N}(p-1)$
`prime_predecessor_lcm`: it is not the Carmichael function when repeated
prime factors are allowed. On a squarefree $N$ it is that group exponent.

Prime support alone does not certify divisibility. For $85=5\cdot17$,
$L=16$ and $84=2^2\cdot3\cdot7$. The radical of 16 divides the radical
of 84, but $16\nmid84$. The missing coordinate is the valuation at 2.
The retained witness is the complete factorization and the exact
quotient $(N-1)/L(N)$, not just a set of primes.

For the three explicit squarefree Carmichael inputs the quotient is

| $N$ | prime factors | $L(N)$ | $(N-1)/L(N)$ |
|---:|---|---:|---:|
| 561 | $3,11,17$ | 80 | 7 |
| 1105 | $5,13,17$ | 48 | 23 |
| 1729 | $7,13,19$ | 36 | 48 |

Korselt's criterion is the classical interface; see Kellner--Sondow,
[*On Carmichael and polygonal numbers, Bernoulli polynomials, and sums
of base-p digits*, Theorem 1.1, p. 2](https://math.colgate.edu/~integers/v52/v52.pdf).
The companion verifies
the divisibilities and independently checks $a^N\equiv a\pmod N$ for
every residue at these three finite inputs. These checks are not an
infinitude argument. Since each input is squarefree with three factors,
all three have $(F,S,U)=(6,6,3)$, so none satisfies $F=S+U$.
More generally, every squarefree composite has $F=S$ and $U>0$.
This separates the Carmichael and divisor-balance predicates immediately.

The taxicab identity $1729=1^3+12^3=9^3+10^3$ is also checked, but no
implication between this identity and Korselt's criterion is asserted.

## 3. An explicit prime-partition map into primitive triangles

Use the ordered primitive triple $(A,B,C)$ with odd leg $A$, even leg
$B$ and hypotenuse $C$. The classical Euclid map and the repository's
Gaussian convention are

\[
\Phi(m,n)=(m^2-n^2,2mn,m^2+n^2),\quad
m>n>0,\quad\gcd(m,n)=1,\quad m-n\text{ odd}.
\]

See [THM-3333, Gaussian-square carrier](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md)
and [THM-3756, odd-square ordinal chart](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md).
The new use here is a source-role receipt for the supplied number, not
a new classification of primitive triples.

**Odd-leg fibre.** For every odd $N>1$, all primitive triples with
marked odd leg $N$ are exactly

\[
\left(N,\frac{v^2-u^2}{2},\frac{v^2+u^2}{2}\right),
\quad uv=N,\quad 1\le u<v,\quad\gcd(u,v)=1.             \tag{4}
\]

Indeed, set $m=(v+u)/2$, $n=(v-u)/2$ for the forward map. They are
coprime and of opposite parity. Conversely Euclid's formula gives
$u=m-n$, $v=m+n$; these are coprime odd integers. The inverse can also
be read from the triangle as $u=\sqrt{C-B}$ and $v=\sqrt{C+B}$.
Every full prime-power block of $N$ must be assigned wholly to $u$ or
$v$. Complementing all assignments exchanges them, giving exactly
$2^{\omega(N)-1}$ unordered coprime splits.

For 561, its three prime blocks give precisely four triangles:

| $(u,v)$ | $(A,B,C)$ | Euclid $(m,n)$ | forward branch word |
|---|---|---|---|
| $(1,561)$ | $(561,157360,157361)$ | $(281,280)$ | $U^{279}$ |
| $(3,187)$ | $(561,17480,17489)$ | $(95,92)$ | $AU^{30}$ |
| $(11,51)$ | $(561,1240,1361)$ | $(31,20)$ | $ADAU$ |
| $(17,33)$ | $(561,400,689)$ | $(25,8)$ | $U^7D$ |

The branch gauge is explicitly
$U(m,n)=(2m-n,m)$, $A(m,n)=(2m+n,m)$,
$D(m,n)=(m+2n,n)$, starting at $(2,1)$ and applying letters left to
right. It agrees with
[THM-3334, Berggren spine](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md).
The inverse chooses respectively $m<2n$, $2n<m<3n$, or $m>3n$,
and strictly decreases $m+n$. The ratio-two primitive case is the root;
ratio three is excluded by parity. This authenticates ancestry without
calling 561 itself the root or claiming Carmichael status survives an
edge of the tree. The root remains $(3,4,5)$.

**Even-leg fibre.** A primitive even leg must be divisible by four.
For $4\mid B$, split $B/2=mn$ into coprime factors, one even and one
odd, and apply $\Phi$. Hence there are $2^{\omega(B/2)-1}$ such triples.
At 560 there are four, with odd legs $78399,3111,1551,1161$.
Passing to $\operatorname{rad}(560)=70$ preserves its three prime labels
but destroys this predicate: 70 is twice an odd number, so cannot be
a primitive even leg. It is even, so cannot be a primitive odd leg or
hypotenuse either. This is an exact obstruction to replacing the source
integer by its prime support.

## 4. Fermat's two-square predicate concerns the hypotenuse role

An odd $C>1$ is the hypotenuse of a primitive integer right triangle iff
every prime divisor of $C$ is $1\pmod4$. This classical Gaussian
factorization result is already used in THM-3334's fixed-hypotenuse
fibre; the current primary teaching reference is Ben Lynn's
[sum-of-two-squares proof](https://crypto.stanford.edu/pbc/notes/numberfield/sumsquares.html).

The obstruction is elementary: if a prime $p\equiv3\pmod4$ divides
$x^2+y^2$, then it divides both $x,y$. Thus it cannot divide a primitive
Gaussian norm, and an odd exponent of such a prime forbids any norm
representation. For sufficiency, choose one Gaussian prime above each
$p\equiv1\pmod4$ and place its full exponent on that choice, avoiding
the conjugate; multiplication supplies a primitive norm representation.
This retains the Gaussian factor choice, not just $\omega(C)$.

The contrast is exact:

\[
1105=23^2+24^2,
\qquad\Phi(24,23)=(47,1104,1105).
\]

Its four primitive hypotenuse triangles are already the first affine
four-point fibre in THM-3334. Both 561 and 1729 have primes $3\pmod4$
to odd exponent, so neither is even a sum of two integer squares, and
neither is a primitive hypotenuse. They still have primitive **odd-leg**
triangles by (4). These predicates must not be interchanged.

Allowing nonprimitive triangles is another change of domain:

\[
264^2+495^2=561^2,\qquad(264,495,561)=33(8,15,17).
\]

The exact norm enumeration of $561^2$ has only $(0,561)$ and $(264,495)$
up to signs and exchange. Primes 3 and 11 force the common content 33.

Finally, multiplying $N$ by $p^2$ preserves whether it is a Gaussian
norm, by the sum-of-two-squares exponent criterion; scaling an explicit
representation by $p$ realizes one direction. That scaling destroys
primitivity. For $p\equiv3\pmod4$, increasing its exponent from two to
four forces every representing pair's common content to increase from
a multiple of $p$ to a multiple of $p^2$. The cases $3^2\cdot5\cdot13$
and $3^4\cdot5\cdot13$ are checked. Thus the proposed exponent change
can preserve one predicate while failing divisor balance and changing
primitive content.

## 5. Receipt scope and reproduction

`Triangle` retains the ordered legs, hypotenuse, Euclid pair and branch
word. `audit_triangle` additionally consumes a source integer and an
explicit role. A valid odd-leg receipt for 561 is rejected as a
hypotenuse receipt for 561. It verifies primitivity, actual coordinates
and ancestry; it proves no property of a different integer obtained by
projection or rescaling. Exact prime labels and distinctness are checked
before constructing $p^aqr$.

Run [the standalone script](../../04-computation/experiments/predecessor_prime_shapes_20261014b.py)
with `python -B -X utf8` and `python -B -O -X utf8`.
The finite universe is direct divisor comparison for $2\le N\le1500$;
1,554 ordered exponent profiles of lengths one through four with entries
one through six; four marked prime triples at exponents two through
twelve; all odd-leg inputs $3,5,\ldots,801$; every primitive Euclid pair
with $m\le70$ checked independently against all three role fibres; the
named norm and Carmichael controls; and 13 malformed or wrong-role
receipts. All checks remain active under optimization. The
[saved JSON](predecessor_prime_shapes_20261014b.json) records 17,426 checks.
Normal and optimized outputs agree; their LF-normalized SHA256 is
`390bd070c51bf9d1e3e6defe43c3065a82aac85bc74cf89b3d5e2db4e022431d`.

This package recovers the prior closed classification and supplies an
exact source-role diagnostic. A subsequent transition theorem must say
which integer, freeness threshold, primitive content and side of the
triangle it retains. Matching prime counts or matching four-element
fibres alone supplies no such transition.
