# Odd-square brackets choose an absolute Farey search

**FINITE-EXACT.** For the lower target $N=2^{121}-2^{60}-1$, the same primitive normal $(a,b)=(7,10)$ minimizes the number of admissible values of $S=ax+by\bmod 2h$ on both unresolved canonical spines. The minima are 336 of 4,930 classes on spine four and 5,264 of 74,530 classes on spine eight. Here admissibility means compatibility with the parity-corrected bilinear equation modulo $2h$. A separate absolute inversion excludes every pair for which $7x$ and $10y$ lie within 512 odd-square brackets of each other, as well as seven other directions at both spines. Neither whole spine is exhausted.

## The inherited triple and two different signs

Set $B_n=((2n-1)^2,(2n+1)^2]$ for $n\ge1$. Exactly the primes $p=2,3,11$ occupy the same bracket as an integer multiple $kp, k\ge2$, since $2p$ lies between them. The inherited proof compares the width to doubling and leaves only the first two brackets. For the midpoint $2p$, exactly one of the three proposals has both neighbors connected by a physical one-arc flip. At $p=2$, a triangle followed by one singleton gives $3\to5$ when a cross-arc is reversed. At $p=3$ the upper neighbor is $7$, and at $p=11$ the lower neighbor is $21$. Their absence uses the tournament theorem. The square-bracket condition alone supplies the three centers. [Odd-square prime brackets](../procgen_numerology_20260923_odd_square_brackets.md) · [THM-1370, H-spectrum omits 7 and 21](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) · [One-arc component grammar](../chanpak_arc_fiber_grammar_20261122/README.md) · [Chan and Pak's published conjecture](https://doi.org/10.1016/j.tcs.2024.114776)

For the accelerated maps $C_\epsilon(p)=(3p+\epsilon)/2^{v_2(3p+\epsilon)}$, every positive odd $p$ satisfies
$$
 C_{-\epsilon}(2p+\epsilon)=C_\epsilon(p)
$$
because the two numerators differ by a factor of two. In particular $C_-(7)=C_+(3)=5$ and $C_+(21)=C_-(11)=1$. An image of either odd accelerated map is prime to three, giving 21 no odd predecessor. The shortened minus map uses a single division on odd inputs and halves even inputs. Its cycle $5,7,10$ supplies the direction $(7,10)$. The known negative cycles of ordinary Collatz have absolute odd minima $1,5,17$, which supply the other six ordered directions below. The argument uses these three explicit cycles and allows additional cycles. [THM-4591, integer cycles of the Terras map to period 301993](../../../01-canon/theorems/THM-4591-integer-cycles-of-the-terras-map-to-period-301993-are-the-five-pythagorean-cycles.md)

## A bracket bounds the number of weighted sums uniformly in the target

**PROVED.** Assume a canonical Farey height
$$
 H=k+qxy+tx+cy+r[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y))
$$
with the inherited $0\le f(x)\le x$, even-valued $f$, and nonnegative coefficients with $q,h>0$. Let $a,b$ be coprime positive integers. Adjoin $B_0=[0,1]$ to allow zero and one. If $H(x,y)=N\ge k$ and the bracket indices of $U=ax,V=by$ differ by at most $d$, put $M=\max(b(t+h),a(c+h))$. Then $S=U+V$ is in the integer interval
$$
 \min\{s\in\mathbb Z_{\ge0}:q s^2+4Ms\ge4ab(N-k-r-e)\}
 \ \le S\le
 \left\lfloor 2\sqrt{ab(N-k)/q}\right\rfloor+4(d+1)^2.
$$

Indeed $\lvert\sqrt U-\sqrt V\rvert\le 2(d+1)$, as the bracket boundaries are spaced by two in the square root coordinate. The upper bound follows from $S=2\sqrt{UV}+(\sqrt U-\sqrt V)^2$ and $qxy\le N-k$. For the lower bound, $4UV\le S^2$ and $b(t+h)U+a(c+h)V\le MS$. The number of sums is at most $4(d+1)^2+O_{a,b,H}(1)$ independently of $N$. In fact $2\sqrt{ab(N-k)/q}\le S+2M/q+2\sqrt{ab(r+e)/q}$. This bounds the number of sums, while the work per sum still depends on the address sizes.

For a fixed sum we apply the inherited quadratic envelope to $(U,V)$, target $abN$, and coefficients $(h\max(a,b),q,ac,bt,ab e,ab r,ab k)$ in the order $(h,q,c,t,e,r,k)$. This bounds the digital and parity allowances even when multiplication changes the parity. Exact reconstruction uses the original physical parities. [Absolute antidiagonal inverse](../chanpak_farey_complementary_exact_20261122/CUBE-WINDOW.md)

## The optimal projected phase, including singular directions

**PROVED.** Form the phase relation $\mathcal R\subset(\mathbb Z/2h)^2$ by requiring
$$
 k+qxy+tx+cy+r(x\bmod2)+e(y\bmod2)\equiv N\pmod{2h}.
$$
The digital contribution vanishes modulo $2h$ since $f$ is even-valued. For any coprime weights, first select all pairs $(u,v)\in\mathcal R$ with $au+bv\equiv S\bmod 2h$. Then all lifts on this line are obtained without inverting either weight modulo $h$. Put $x=u+2h j$ and solve
$$
 j\equiv a^{-1}(S-au-bv)/(2h)\pmod b,\qquad y=(S-ax)/b.
$$
Intersect the progression with the absolute envelope and evaluate the actual digital reads. In particular both 5 and 17 divide $h$ at spine four, and the directions $5:17,17:5$ are covered by the same formula.

For an odd prime $p\mid h$, fix the physical parities $i,j$. If $q\ne0\bmod p$, put $u=qx+c, v=qy+t$. The local equation reads $uv=D=q(N-k-ri-ej)+ct$, and $au+bv=qS+ac+bt$. When $Dab\ne0\bmod p$, the projection takes exactly
$$
 \frac{p+\chi_p(abD)}{2}
$$
values. Two inputs $u,u'\ne0$ give the same projection exactly when $u'=u$ or $u'=bD/(au)$. This involution has $1+\chi_p(abD)$ fixed points. If exactly one weight vanishes and $D\ne0$, the image has $p-1$ values. The replay handles $q=0$ or $D=0$ directly.

**FINITE-EXACT optimum.** For these squarefree $2h$, each primitive normal determines one projective point at every prime $p\mid 2h$. Multiplying the normal by any local unit scales the set of acceptable $S\bmod p$ by that same unit, preserving all intersection sizes. For each of the four physical parity pairs compute its projected set at each prime, including $p=2$. The allowed sums modulo $2h$ are the union of their four CRT products. Inclusion–exclusion over the 15 nonempty subsets of the four pairs gives the exact cardinality, multiplying the local intersection sizes. The code takes the minimum over all $\prod_{p\mid 2h}(p+1)$ projective tuples and produces an actual positive coprime normal attaining it.

| Spine | Projective tuples | Tuples achieving minimum | Minimum projected classes | Full modulus | Smallest max weight achieving minimum |
|---:|---:|---:|---:|---:|---:|
| 4 | 9,720 | 112 | 336 | 4,930 | 10, with $(7,10)$ |
| 8 | 139,320 | 1,792 | 5,264 | 74,530 | 10, with $(7,10)$ |

For $(7,10)$ the sum reveals left parity modulo 2 and the constrained left address modulo 5. At 29 the parity pairs share the same 14 allowed values. At the remaining prime the two right-parity images have union size 12 at 17 and 188 at 257. Thus the products are $2\cdot1\cdot14\cdot12=336$ and $2\cdot1\cdot14\cdot188=5264$. Optimality is for this specific projection and phase. Relative runtime also depends on slab width and the distribution of roots along a line.

**PROVED periodic reduction, FINITE-EXACT classification.** For every odd exponent $e\ge3$, use the lower target $N_e=2^{e-1}(2^e-1)-1$, whether or not $2^e-1$ is prime. The orders of 2 at $5,17,29,257$ are $4,8,28,16$, so all four odd-prime phase tables depend only on $e\bmod112$ and the parity table is fixed. The two E61 tables agree with the corresponding $N_e$ tables simultaneously precisely for $e\equiv61\bmod112$ when agreement means equality of the target residues. Indeed the individual residue conditions on $e$ are $1,2\bmod4$ at 5, $5\bmod8$ at 17, $5,15\bmod28$ at 29 and $13\bmod16$ at 257, obtained from $X(X-1)=X_0(X_0-1)$ with $X=2^e,X_0=2^{61}$. Their intersection gives the assertion.

The optimal direction has a wider range than equality of those targets. Of the 56 odd classes, $(7,10)$ attains the minimum on 32 classes at spine four and on 22 at spine eight. It attains both minima exactly for
$$
 e\bmod112\in\{1,5,21,29,33,43,49,61,71,75,77,85,103\}.
$$
Thus the same two consecutive cycle values give the strongest phase projection for 13 of the 56 odd exponent classes. The classification concerns the projection and places no absolute height by itself.

**PROVED, compatibility with tangent approximation.** For any positive irrational slope $\theta$, there are infinitely many distinct reduced positive normals $(a,b)$ projectively equal to $(7,10)$ modulo every prime of either spine and satisfying
$$
 |\theta-a/b|<4m^2/b^2,\qquad m=\mathrm{lcm}(4930,74530)=1267010.
$$
Take two successive convergent vectors $(p_k,q_k),(p_{k-1},q_{k-1})$. Their determinant is a unit, so a combination with coefficients $1\le u,v\le m$ reaches $(7,10)\bmod m$. For its numerator $A$ and denominator $B\le2m q_k$, $|B\theta-A|<2m/q_k$ by the adjacent convergent bounds. Any common factor of $A,B$ is a unit modulo $m$ since $(7,10)$ is primitive, so reduction preserves every projective class. The approximants tend to the irrational slope, which forces infinitely many distinct fractions. For a single spine one may use its own smaller modulus. The constant does not guarantee a useful short normal in every finite search slab.

## Absolute exponent-61 search

**FINITE-EXACT.** Each entry lists all weighted sums in the proved interval at $d=512$, followed by the number of original height evaluations. The replay tests the full interval, a superset of the bracket-conditioned pairs. It found zero hits in all 16 runs. The digital addresses otherwise vary independently and may have arbitrary bit patterns.

| $(a,b)$ | Sums, spine 4 | Reads, spine 4 | Sums, spine 8 | Reads, spine 8 |
|---|---:|---:|---:|---:|
| $(1,1)$ | 1,052,677 | 53,331 | 1,052,677 | 303 |
| $(1,5)$ | 1,052,679 | 79,972 | 1,052,680 | 480 |
| $(1,17)$ | 1,052,686 | 107,874 | 1,052,690 | 645 |
| $(5,1)$ | 1,052,682 | 80,189 | 1,052,682 | 484 |
| $(5,17)$ | 1,052,686 | 32,711 | 1,052,689 | 200 |
| $(17,1)$ | 1,052,696 | 108,707 | 1,052,696 | 665 |
| $(17,5)$ | 1,052,696 | 32,591 | 1,052,696 | 190 |
| $(7,10)$ | 1,052,684 | 22,196 | 1,052,684 | 133 |

The bracket index, dynamics, congruence phase and literal tournament height retain different data. The new step toward an absolute inverse is to use the optimal local phase along a geometrically appropriate normal, keeping the exact line and both address reads. Run from this directory with the Python standard library.

```sh
python3 -B square-directions.py
python3 -B square-directions.py --exponents
python3 -B square-shells.py --check
python3 -B square-shells.py --distance 512
```

## Exact cost of using a single direction on both spines

**FINITE-EXACT across the full exponent period.** Let \(m_4,m_8\) denote the minimum projected support sizes attained separately on spines four and eight at \(N_e\). In 40 of the 56 odd classes modulo 112 one direction attains \((m_4,m_8)\). In the other 16 classes, exactly those with \(e\equiv7,11,13,23\pmod{28}\), the complete Pareto frontier consists of
\[
 (m_4,15m_8/14),\qquad (15m_4/14,m_8).
\]
Thus optimizing either sieve forces at least \(1/14\) more residues at the other. Both tradeoffs are attained. At exponent 61 there is no such conflict. The underlying opposition between the two preferred quadratic characters modulo 29 is identified in the [all-exponent Pell classification](PELL-EXPONENTS.md). For a nondegenerate hyperbola the two quadratic-character projections at 29 have 14 and 15 values, respectively.

To calculate the Pareto frontier, fix the normal projectively at the shared primes 2, 5 and 29. This gives 540 assignments. With them fixed, minimize each support independently over the remaining direction at 17 or 257. The two chosen private directions give a joint primitive positive normal. Indeed, use CRT to prescribe a nonzero residue pair representing each local projective point. Choose a positive lift \(b\) of the second coordinate. For each prime of \(b\) outside the modulus, additionally prescribe \(a=1\) there. The primes of \(b\) inside the modulus already force a nonzero first coordinate. CRT now gives a positive \(a\) coprime to \(b\). Every other pair with this shared assignment is weakly worse in each coordinate, so removing dominated pairs across the 540 assignments is exact. Run `python3 -B projection_tradeoff.py` for all exponent classes and the five distinct two-point frontiers.
