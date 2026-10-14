# Order groups, distinct-prime spacing, and the exponent's divisor budget

**PROVED:** source-specific reciprocal-support bounds; no Giuga Mersenne with
an odd prime-power exponent; necessary conditions `tau(E)>=1920` and
`omega(E)>=8` for a jointly Giuga--Carmichael Mersenne with odd exponent.
The divisor-count statement uses a finite exact certificate for every exponent
size, not a finite search over Mersenne integers.
**FINITE-EXACT:** the declared rational certificates and independent controls.
No joint source is constructed; Collatz ROOT remains **OPEN**.

Artifacts: [checker](../../04-computation/experiments/mersenne_order_budget_20261014h.py)
and [saved JSON](mersenne_order_budget_20261014h.json).

## 1. Inheritance and the retained objects

Put $M_E=2^E-1$ for odd $E\ge3$ and
$S(M_E)=\sum_{p\mid M_E}1/p$, over distinct prime divisors.
The [preceding source-owned packet](mersenne_lattice_child_guards_20261014g.md)
proved

\[
S(M_E)<\frac{\tau(E)-1}{5}.
\]

It also proved that a hypothetical joint Giuga--Carmichael source has all
$p\mid M_E$ congruent to seven modulo eight and Giuga index $\kappa\ge4$,
with $4\mid\kappa$. Thus $S=\kappa+1/M_E>4$. The resulting bound
$\tau(E)\ge22$ was a necessary threshold obtained by assigning every order
group the same coarse capacity $1/5$.

The live concepts are unchanged exponent / exact multiplicative order /
distinct-prime spacing / divisor multiplicities / positive weighted sums /
source-specific terminal. The least-used sidecar is distinctness *inside*
an order group. Its use changes the bound substantially. The hostile boundary
is that a capacity estimate above four asserts no actual prime support and
no simultaneous Korselt/Giuga solution.

The current Collatz source is still $E_*=99\,708\,991\,741$.
Its earlier arithmetic exclusion and the legal but unpaid later descent are
retained. Stronger arithmetic restrictions do not authenticate a new Collatz
edge at this source.

## 2. Exact budgets for each marked order

For $p\mid M_E$, let $d=\operatorname{ord}_p(2)$. Then $d>1$ is an odd divisor
of $E$. Every such prime is $1\pmod{2d}$. For a joint source also require
$p\equiv7\pmod8$. The exact allowed integer progressions are

\[
p_i\ge L_d+s_d(i-1),\qquad
(L_d,s_d)=
\begin{cases}
(2d+1,2d),&\text{without the joint premise},\\
(2d+1,8d),&\text{joint and }d\equiv3\pmod4,\\
(6d+1,8d),&\text{joint and }d\equiv1\pmod4.
\end{cases}\tag{1}
\]

Here the $p_i$ are the distinct order-$d$ primes in increasing order.
The joint rows follow by solving $p=kd+1\equiv7\pmod8$: the least positive
even $k$ is two or six, respectively. Composite progression entries may be
included in an upper bound, but are not prime candidates certified by it.

Their product divides $2^d-1$. If there are $r_d$ primes, then
$L_d^{r_d}<2^d$, so

\[
S_d<\frac{d}{L_d\log_2 L_d}.
\tag{2}
\]

Keep the strict integer packing capacity as well:

\[
R_d=\max\{r\ge0:L_d^r<2^d\}.
\]

The checker computes $R_d$ with certified rational logarithm intervals,
without forming either enormous power. Since $L_d>1$ is odd,
$d/\log_2L_d$ is not an integer, so successive interval refinement terminates.
We obtain two further, non-strict upper bounds:

\[
S_d\le\frac{R_d}{L_d},\qquad
S_d\le\frac1{L_d}+
 \frac1{s_d}\log\frac{L_d+s_d(R_d-1)}{L_d}\quad(R_d\ge1).
\tag{3}
\]

For $R_d=0$ the group is empty. The second inequality uses the decreasing
reciprocal function: its sum over the first $R_d$ progression entries is at
most its first value plus the integral from one to $R_d$. Equality is possible
when $R_d=1$; these are not advertised as strict bounds. The exact native order
and primality tests are stronger than the relaxed progression and product
conditions, so (2)--(3) are upper certificates only.

### The unchanged giant source

The full, independently checked exponent factorization is
$E_*=47779\cdot2086879$. Its three nontrivial divisors give:

| $d$ | Unconditional $L_d$ | Exact $R_d$ |
|---:|---:|---:|
|47779|95559|2887|
|2086879|4173759|94888|
|99708991741|199417983483|2656285255|

The rational upper bounds in the JSON prove

\[
\boxed{S(M_{E_*})<10^{-4}.}\tag{4}
\]

The computed upper is approximately $0.00009683380371624$; the approximation
is display only. Since $1/M_{E_*}<1/E_*$ and the exact upper plus $1/E_*$
is still below $10^{-4}$, both signed indices lie strictly between zero and
$10^{-4}$. Positivity of the minus index uses that $E_*$, hence $M_{E_*}$,
is composite. The joint-only progression would reduce the upper to about
$0.00003594589634409$, but **that is an additional hypothetical bound**, not
an unconditional assertion that the actual prime factors are all seven mod8.

`source_budget` authenticates the full factorization of the same exponent.
Its explicit divisor-enumeration cap causes rejection rather than silently
dropping groups. No prime factorization or expansion of $M_{E_*}$ is used.

## 3. A closed divisor sum from distinct-prime spacing

Use the unconditional progression in (1). Since $r_d<d$, we have

\[
S_d\le\sum_{i=1}^{r_d}\frac1{2di+1}
 <\frac{H_{r_d}}{2d}
 \le\frac{1+\log r_d}{2d}
 <h(d):=\frac{1+\log d}{2d}.
\tag{5}
\]

Empty groups also satisfy $0<h(d)$. Thus

\[
S(M_E)<\mathcal H(E):=\sum_{d\mid E,\ d>1}h(d).
\tag{6}
\]

Let $E=\prod p^{a_p}$ and define the exact geometric sums

\[
g_p=\sum_{j=0}^{a_p}p^{-j},\qquad
j_p=\sum_{j=0}^{a_p}j p^{-j},\qquad R=\prod_p g_p.
\]

Expanding the divisor product and differentiating each finite geometric sum,
or simply collecting each prime's contribution to $\log d$, gives

\[
\boxed{\mathcal H(E)=\frac12\left[
R-1+\sum_p j_p\log p\prod_{q\ne p}g_q\right].}\tag{7}
\]

This is an exact product-of-choices representation of the divisor weights.
Its coordinates are exponent-prime choices, with multiplicity retained.
It is not an identification with a Euclidean square packing or an Apéry
coefficient formula. The useful capacity mechanism is literal: distinct
primes occupy separated arithmetic slots, and their product has a fixed
logarithmic budget.

## 4. Two all-height consequences

### Prime-power exponents cannot produce Giuga Mersennes

The function $h(x)$ is decreasing for $x>1$, since
$h'(x)=-\log x/(2x^2)$. For $E=p^a$, replacing $p$ by three and allowing
all exponents $j\ge0$ only increases (6). Therefore

\[
\mathcal H(p^a)
 <\frac1{2(3-1)}+\frac{3\log3}{2(3-1)^2}
 =\frac14+\frac{3\log3}{8}<\frac23<1.
\tag{8}
\]

Thus **no odd prime-power exponent, of any multiplicity, yields a Giuga
Mersenne number**. The final inequality is certified with rational bounds for
$\log3$ (for example $\log3<10/9$ suffices).

### Joint exponents need at least eight distinct prime factors

For fixed $k=\omega(E)$, replace its prime bases by the first $k$ odd primes
and let every multiplicity tend to infinity. The same monotonicity gives

\[
\mathcal H(E)<\frac12\left[
\prod_{i=1}^k\frac{p_i}{p_i-1}-1+
\left(\prod_{i=1}^k\frac{p_i}{p_i-1}\right)
\sum_{i=1}^k\frac{\log p_i}{p_i-1}\right].\tag{9}
\]

For $k=7$, with $p_i=3,5,7,11,13,17,19$, the exact certificate bounds this
by $3.99$ (displayed value about $3.987929563322$). Since a joint source has
$S>4$, it follows that

\[
\boxed{\omega(E)\ge8.}\tag{10}
\]

This holds for arbitrary multiplicities. The eight-base infinite envelope
exceeds four; that is a limitation of this bound, not an example of a joint
source.

## 5. A finite profile certificate improves 22 to 1920

There are only finitely many nonincreasing positive tuples
$a=(a_1,\ldots,a_k)$ with $\prod_i(a_i+1)\le B$.
They cover every odd exponent with $\tau(E)\le B$, regardless of its size.
The following two reductions justify using each tuple on the smallest odd
prime bases.

1. Replacing a prime base by a smaller distinct base decreases each matched
   divisor, and hence increases its $h$ weight.
2. For $p<q$ and multiplicities $a\ge b$, moving the larger multiplicity to
   $p$ increases the sum. The common exponent square is unchanged; pair each
   remaining term $p^i q^j$, $b<i\le a$, $0\le j\le b$, with
   $p^j q^i$. The first is smaller. Other prime factors may be fixed throughout.

This changes the exponent only to compute an upper envelope; it does not
substitute another source into a Collatz certificate.

The exact checker enumerates **all 16,894 profiles with divisor count at most
1919**, including the harmless empty profile. For speed, a multiplicity above
24 is replaced by its infinite geometric sums in (7). This remains an upper
bound: before simplifying, the right side is a positive polynomial in the
$g_p,j_p$, apart from the constant minus one. Logarithms have rational upper
bounds, never floating-point acceptance tests.

Every resulting rational envelope is below $3999/1000<4$. The largest occurs
at profile $(3,2,1,1,1,1,1,1,1)$, whose divisor count is 1536; its displayed
value is $3.998100372920$. Thus

\[
\boxed{M_E\text{ jointly Giuga--Carmichael},\ E\ge3\text{ odd}
\quad\Longrightarrow\quad\tau(E)\ge1920.}\tag{11}
\]

A lower rational bound for the relaxed envelope of profile
$(4,2,1,1,1,1,1,1,1)$, with divisor count 1920, exceeds four.
This explains the stated stopping boundary of the present certificate.
It does not establish that all other profiles at or beyond 1920 are feasible.
For squarefree $E$, (11) in particular requires at least eleven prime factors.
For general $E$, (10) and the full profile bound are useful additional guards.

The numerical threshold 22 was therefore a coarse capacity artifact, not an
identified modular level. The current improvement is a proof about the
arithmetic divisor set; it is not a Langlands correspondence.

## 6. A precise Paley representation, and what it forgets

There is a direct, elementary connection to the inherited circulant discussion
in [THM-4644, Fermat primes, rigid circulants and Wieferich orbits](../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md).
For a prime $p\equiv3\pmod4$, orient $x\to y$ in the Paley tournament when
$y-x$ is a nonzero square. Multiplication by two preserves every edge iff two
is a square modulo $p$. Thus:

- at a hypothetical joint factor $p\equiv7\pmod8$, doubling is a tournament
  automorphism; zero is fixed, and every nonzero orbit has the retained odd
  length $d=\operatorname{ord}_p(2)\mid E$;
- at $p\equiv3\pmod8$, doubling reverses every edge instead. The case $p=11$
  is an exact small hostile to forgetting this sign.

The code checks complete pairwise relations for $p=7,23,31,47,71$ and the
reversing controls $3,11,19$, together with their full nonzero orbit partitions.
This represents the *order label*. It does not authenticate the product
$M_E$, the Giuga index, a Hamiltonian count, or a Collatz transition.

The earlier [level-22 oldform packet](collatz_level22_oldform_bridge_20261005.md)
uses a different specified interface: Hecke prime two, reduction characteristic
three, conductor eleven, and a retained two-dimensional operator. The role of
11 there is not the order eleven of doubling modulo 23 here. Shared numbers
do not identify these representations.
`THM-4652` and `HYP-9274` were **RESERVED** at intake; incoming commit
`294e754ecd` records subsequent Paley results. They are not dependencies here.
No Apéry irrationality or Paley-count theorem is inferred
from the geometric products in (7).

## 7. Exact arithmetic, reproduction, and scope

For $1\le x\le2$ put $z=(x-1)/(x+1)\le1/3$ and use

\[
2\sum_{j=0}^{n-1}\frac{z^{2j+1}}{2j+1}
\le\log x\le
2\sum_{j=0}^{n-1}\frac{z^{2j+1}}{2j+1}
+\frac{2z^{2n+1}}{(2n+1)(1-z^2)}.
\]

Range reduction by an exact power of two handles every larger rational input.
All comparisons, including the profile cutoff, use fractions. The JSON retains
the exact envelopes and a SHA256 of the complete deterministic profile stream.

```text
python -B 04-computation/experiments/mersenne_order_budget_20261014h.py
python -B -O 04-computation/experiments/mersenne_order_budget_20261014h.py
```

Independent controls use literal defining powers for every odd order 3 through
129, full prime support of $M_E$ only for odd $3\le E\le31$, direct divisor
enumeration versus (7), an independently generated small profile set, and the
complete finite Paley graphs just listed. Malformed types, partial source
factorizations and exhausted explicit enumeration caps are rejected. No
production API discovers a ROOT suffix or claims that a failed capacity
exclusion is an actual arithmetic or Collatz completion.

Normal and optimized runs reproduce **9,701 checks** and identical saved JSON.
LF-normalized SHA256:
`fc225d1228889f561f71b6c14bf0538b1d70bdc476e85b9dedb6b95fe27b173a`.
The separately enumerated profile certificate covers 16,894 profiles.
Independent proof/API and optimized-versus-saved audit: **PASS**, including an
independent multiplicative-partition count of all 16,894 profiles.
