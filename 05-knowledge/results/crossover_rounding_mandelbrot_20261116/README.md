# The crossover and the carry: how close becomes exact

**Proved:** a converse-ratio cutoff, an exact finite-precision rounding law for binary tournament continuations, a no-go for fixed dyadic rails tracking a fixed decimal offset, and literal realizability of Mandelbrot hyperbolic-component counts at prime periods and at twice an odd prime. **Open:** target-dependent construction of a physical source with a simultaneous error and residue certificate for every admissible prime.

## 1. Order five is the first and largest peak once chiral classes exist

Let $a_n$ count tournament isomorphism classes, $b_n$ the self-converse classes, and $c_n=(a_n-b_n)/2$ the chiral pairs. Direct census gives
$$
(n,a_n,b_n,c_n)=(4,4,2,1),(5,12,8,2),(6,56,12,22),(7,456,88,184).
$$
In particular $b_5/c_5=4$. For $n\le3$, $b_n/a_n=1$ and $b_n/c_n$ is undefined. Among $n\ge4$ the two ratios order the sizes the same way. Every self-converse class has a reversing involution: take the odd power of a reverser yielding its order-two part, since the tournament automorphism group has odd order. Such an involution is a product of $m=\lfloor n/2\rfloor$ transpositions with at most one fixed vertex. Standardize its pairs and orient the arc inside each pair in a prescribed direction by relabeling. The remaining freely oriented orbits count
$$
\frac{\binom n2-m}{2}=\left\lfloor\frac{(n-1)^2}{4}\right\rfloor.
$$
This gives a surjection from the resulting binary frames to self-converse classes. Each unrestricted class accounts for at most (n!) labeled tournaments, so $b_n/a_n\le n!/2^{\lfloor n^2/4\rfloor}$. At $n=8$ this is (315/512<2/3), and its next-to-current quotient is $(n+1)/2^{\lfloor(n+1)/2\rfloor}<1$ thereafter. Thus order five is the unique maximum on the nontrivial domain. On the converse-merged quotient, the fraction of nodes fixed is $b_n/(b_n+c_n)$, also maximizing at five (there (8/10)); this is a different denominator.

The [archived anti-Burnside analysis](../session29_framed_profiles_20261022/session29-selfconverse-asymptotics/REPORT.md) uses an odd partition $\lambda$ of $m$ with $\ell$ parts and weights $z_\lambda=\prod d^{a_d} a_d!$. Its relative defect corrections shrink exponentially, but a decisive *absolute* bound goes the other way: each omitted type costs at least 1 in the unnormalized count. Indeed,
$$
z_\lambda\le 2^{m-\ell}\ell!\le 2^{m-\ell+\ell(\ell-1)}
\le 2^{m-\ell+2\sum_{i<j}\gcd(d_i,d_j)},
$$
the last number being the even-order numerator; the odd-order numerator is larger. The inequality is strict for $m\ge2$. Thus even very sharp **relative** truncated defect expansions cannot simply be rounded to the exact integer. They must aggregate the omitted types or bring additional arithmetic information.

The [order-five switching packet](../planar_shadows_petersen_switching_20261116/README.md) is a third object: the labeled switching space has dimension 6 (63 nonzero vectors), while the even sheet quotiented by converse is dimension 4 (15 nonzero vectors). These are not the 12 tournament isomorphism classes. There is not even a surjective map from length-six binary words onto length-four words intertwining the two cyclic shifts: its image would be fixed by shift to the sixth power, hence by the square of the target shift. The four possible images then have period dividing 2. Matching sizes does not transport the dynamics.

## 2. MRB: an actual crossover and an absolute tail correction

For the [Burns alternating constant](https://www.marvinrayburns.com/), set $a(x)=\exp((\log x)/x)-1$. Its maximum is at $e$. The sign of (a''$x$) is that of
$$
(1-\log x)^2+x(2\log x-3).
$$
This is negative on ([1,e]); after $e$ it strictly increases and crosses zero just once, at $x\approx4.36777$. Thus the macroscopic curve is not just falling; its tail is convex, so its positive successive drops decrease. With $S_{m-1}$ the alternating partial sum and $m\ge5$, summing first differences gives
$$
0\le (-1)^m(C_{\rm MRB}-S_{m-1})-a(m)/2
\le (a(m)-a(m+1))/2=O(\log m/m^2).
$$
This is a one-sided correction, one order better than the uncorrected alternating envelope. It is the **bound on the remainder**, not the resemblance of the crossover to 5, that transports to the tournament problem. The [previous MRB audit](../session31_mirror_mrb_20261022/mrb/MRB_CONTOUR_AUDIT.md) proves that the complete inverse-power tail, every negative-axis branch jump and any finite prefix *together* do not determine the sum in its analytic class. The [archived branch-residue synthesis](../signature_axis_cat_map_rr_mrb_20261022.md) likewise distinguishes equal-looking small signatures from an exact transfer identity.

## 3. A real-plus-dyadic rounding theorem for tournament counts

Use the normalized physical binary address of the [prime inverse](../prime_packing_inverse_20261116/README.md):
$$
H(z)=A(z+1)+g f(z)+J,
\qquad f(0)=1,quad f(2u+1)=f(u),\quad f(2u)=f(u)+u+(u\bmod2)
$$
(the last case for positive $u$). A zero run $[a,b)$ means bits $a$ through (b-1) are zero and bit $b$ equals one; do not include leading zeros. Define $R=R(z)$, $w_{a,b}=1-2^{a-b}$, and $\alpha_a=(z\bmod 2^a)/2^a$. Expanding the recurrence run-by-run gives
$$
f(z)=1+z\sum_{[a,b)}2^{-a}w_{a,b}
   +\sum_{[a,b)}(1-\alpha_a w_{a,b}).
$$
For the proof, every quotient term $\lfloor z/2^{i+1}\rfloor$ within the run uses the same discarded low part $z\bmod 2^a$, whose contribution is geometrically summed. Exactly one of the next bits in the run is 1, so all of the parity terms total 1. The first sum reads the reversed zero pattern, while the second records the interaction with preceding digits.

Replace each $\alpha_a$ by its first $\ell$ binary fractional digits, rounded down, and call the resulting real upper approximation $U_\ell(z)$ after applying $A(z+1)+g(\cdot)+J$. Then if (R>0),
$$
0\le U_\ell(z)-H(z)<gR2^{-\ell}.
$$
If $R=0$ it is exact. Independently, obtaining $r=H(z)\bmod 2^m$ from the recurrence uses only the $m$ bits immediately above each zero bit, including the next-bit term for $m\ge1$; long-distance interactions vanish in this modulus. Provided
$$
2^{m+\ell}\ge gR,
\quad H(z)=r+2^m\left\lfloor (U_\ell(z)-r)/2^m\right\rfloor.
$$
This is an **absolute integer reconstruction** law at logarithmic total local precision in the number of zero runs. At $m=1$ the residue is known odd, so for the 30,011, 30,013, 1,000,000,007, 1,000,000,009 literal witnesses only 2, 3, 4, 5 preceding bits per run suffice, respectively, together with the exact reversed zero pattern. The billion-scale upper overshoots are under 0.37, even though the main term is about a billion. The formula certifies a *given* address; it is not a polylogarithmic algorithm for selecting one.

**Inverse target theorem.** A prescribed odd $P$ lands precisely when some physical address has $H(z)\equiv P\pmod{2^m}$ and $P\le U_\ell(z)<P+2^m$, at the stated precision. This puts the target into an explicit one-sided digital tube. The conditions must hold on the *same* source and word. Odd composites give a counterexample to any argument treating them separately: they hit every permitted residue for each modulus $M$, via ((1+2M)b) for odd $b\ge3$ in the desired class, and they come within distance 3 of every sufficiently large integer (use positive odd multiples of 3). Still none is prime. The next nontrivial task is to select the source and the residue in the real tube jointly, not to increase residue reachability separately.

## 4. Baselines, envelopes, and physical corrections

| Inputs | Exact decomposition | What transports |
|---|---|---|
| 10007,10009 | $10^4+8\pm1$ | The archived two-edit mirror lineage uses a variable exterior; the squared decimal anchors are not a fixed-size-step tournament operation. |
| 30021,30027 | $3(10^4+8\pm1)=30024\pm3$ | An ordinal directed triangle preserves the literal seam but triples the digit step. |
| 30011,30013 | $30012\pm1$ | Different sources and contacts; actual midpoint shift is $-12$ relative to the tripled window. |
| 1000000007,1000000009 | $10^9+8\pm1$ | Same decimal offsets, but different normalized sources and different runs. |
| 6015,6017,6019 | $756\cdot2^3-33+\{0,2,4\}$ | $756\cdot2^s-33=3(63\cdot2^{s+2}-11)$ explains the first entry; the other two require different placements, not the ordinal triangle. |

The [756 row](../first25_spread_target_returns_20261116/TARGETS.md) has an exact decoder: the digit $j\in\{0,1,2\}$ is read modulo 756 from $N+33=756\cdot2^s+2j$; then the remaining quotient must be a power of two. At $s=3$ the baseline is $3\cdot2005$, and the other entries are not multiples of 3. The equal growth coefficient alone does not match all endpoint data; the original triple even has three different deletion counts. This is precisely why a scalar equality should not be promoted to a stable continuation.

**Fixed-rail obstruction.** In any finite bank of fixed dyadic rails $A2^s+B$ with (A>0), only finitely many exact $D10^k+d$ can occur for fixed positive integer $D$ and fixed integer $d$. If $B\ne d$, once $k+v_2(D)>v_2(B-d)$, the valuation of $A2^s=D10^k+d-B$ is exactly $v_2(B-d)$; this fixes $s$. If $B=d$, valuation at 5 fixes $k=v_5(A)-v_5(D)$. Thus the recurring decimal offsets cannot be an infinite fixed-bank law. Families with moving exteriors, varying dyadic leading coefficient, or a new exact carry are genuinely needed.

## 5. Mandelbrot counts become literal tournament counts on two infinite period families

[Schleicher's center count](https://www.numdam.org/article/AST_2000__261__405_0.pdf) gives $M_n=\sum_{d\mid n}\mu(n/d)2^{d-1}$: the degree counts all centers of periods dividing $n$, and exact periods come from subtraction. At $n>1$, $M_n\equiv\mu(n)\pmod2$, so nonsquarefree periods give even counts (and therefore never tournament Hamiltonian counts). At a prime period $p$, $M_p=2^{p-1}-1$; at $2p$ with $p$ an odd prime,
$$
M_{2p}=2^{2p-1}-2^{p-1}-1=(2^p+1)(2^{p-1}-1).
$$
The [archived Fibonacci fan, THM-4748](../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md) realizes every $2^k-1$ for $k\ge4$ in one strong tournament. A triangle handles $k=2$ and a singleton $k=1$; the excluded exponent 3 never equals (p-1) for prime $p$. A single longest reversed arc realizes $2^p+1$. Ordinal join then realizes $M_{2p}$ with *literal factors*. Thus **every prime period and every twice-odd-prime period has its Mandelbrot component count realized**, including 15, 63, 495, 8127. This is a construction for infinite families, not a bijection of components with paths.

More generally for distinct primes (p,q), with $X=2^q$,
$$
M_{pq}=M_q\big(X^{p-1}+2X^{p-2}+\dots+2^{p-1}-1\big).
$$
This follows by dividing $X^p-X-2^p+2$ by (X-2); by symmetry $M_p$ also divides $M_{pq}$. The next concrete term is the $p=3$ quotient $2^{2q}+2^{q+1}+3=(2^q+1)^2+2$. A uniform physical (+2) seam near those square products would extend the theorem to periods $3q$; no such seam is asserted here.

The $15,63$ fixed grids for angle doubling at exponents 4 and 6 are **not** themselves the period-5 and period-7 centers; primitive angles at exponents 4 and 6 count 12 and 54. At exponent 6 the fixed grids of sizes 3 and 7 meet at 0 (9 total old angles removed), while angle (1/21) already has exact period 6. Independently $2^6-1=9\cdot7=3\cdot21$ is the old-prime exception in the Mersenne sequence. The fan explains how to physically attain 63 without realizing either forbidden factor, not why those factors are forbidden. This keeps the [archived Collatz/Mandelbrot angle dictionary](../../../01-canon/theorems/THM-4714-collatz-parity-numbers-are-mandelbrot-external-angles-and-minkowski-farey-addresses.md) distinct from the tournament count theorem.

Run `python replay.py` for exact rational run bounds and the displayed arithmetic identities. Small exhaustive controls test the rounding algorithm, not the unproved joint source-selection condition.
