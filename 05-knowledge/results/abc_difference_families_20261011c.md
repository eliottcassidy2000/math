# Difference families, golden norms, and an arithmetic geometry boundary

**Status:** PROVED elementary identities, factorizations, finite-sieve obstruction,
and rational-coordinate colouring; CONDITIONAL radical limit assuming ABC;
FINITE-EXACT controls described below. No ABC proof, prime-density theorem,
twin-prime infinitude, or six/seven-colour plane result is claimed.

The useful connections here retain different data: an additive packet retains
its primitive scale; a golden norm retains its sign; a unit-distance construction
retains exact coordinates and a component origin. Forgetting those data does not
create a prime predictor or a geometric colouring obstruction.

## 1. Inheritance and the live concepts

The five objects being compared are gap-two prime pairs, primitive ABC packets,
Fibonacci/Lucas norms, finite residue sieves, and exact unit-distance carriers.
The closest current mechanisms are:

- [THM-2413, prime-index affine drift and twin-centre weld](../../01-canon/theorems/THM-2413-prime-index-affine-drift-and-twin-center-weld.md): the unique overlapping twin-prime pair is already inherited.
- [THM-3833, ABC-conditional cube radicals](../../01-canon/theorems/THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness.md) and [THM-4108, reciprocal-power radicals](../../01-canon/theorems/THM-4108-abc-conditional-reciprocal-power-radicals-and-lrc-gauge-obstruction.md): radical growth requires the stated ABC premise and primitive height.
- [Golden prime clocks](golden_prime_clocks_20261003.md): exact signed norm and coefficient arithmetic precede any primality interpretation.
- [Difference diagonals and golden phase recursion](difference_diagonals_golden_phase_recursion_20261004.md): a difference quotient loses position along its diagonal.
- [Signed boundary parity](plane_six_reframes_parity_20261011b.md): a retained arithmetic colouring can be checked by the actual graph, without replacing unit geometry by a graph analogy.

The [ABC/IUT reference ledger](../reference/CORE-PAPERS-ABC-IUT.md) keeps ABC
OPEN and does not permit a disputed implication to become a proved dependency.
The conditional statement used here is the usual one: for every
$\epsilon>0$ there exists $K_\epsilon>0$ such that positive coprime
$a+b=c$ implies

$$c\le K_\epsilon\operatorname{rad}(abc)^{1+\epsilon}.$$

This is also the explicitly conjectural formulation in
[Darmon--Granville, Section 5.2](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf).
No effective numerical value of $K_\epsilon$ is assumed.

Historical boundary, now repaired in the producer and saved-output prose: the original
[Fibonacci/Lucas census](fib_lucas_primes_s90bo.out) printed the wrong sign in
its first general norm formula, although its following odd-index identity was
correct. It also incorrectly made $\gcd(F_n,L_n)$ depend on parity. The correction
banners preserve the old finite numerical data. The correct
formulas are (4) and (7) below. Coprimality alone supplies no probabilistic
independence of the two primality events. Historical finite percentages are not
asymptotic prime densities. The present packet does not rerun that old census.

## 2. Twin centres: an exact local classification and a sieve obstruction

“Only 2 and 3 are consecutive primes” is correct when *consecutive* means
integers differing by one. It does not mean consecutive entries in the ordered
list of primes. Apart from 2, primes are odd, proving the difference-one claim.

If $p,p+2,p+4$ are prime, then $p=3$. The case $p=2$ fails; otherwise one of
the three odd numbers is divisible by 3 and must be 3 itself. Thus the graph
with prime integers as vertices and difference-two edges has exactly one
three-vertex component, $3-5-7$. Every other nontrivial component is one edge.
This classifies possible overlap, not which further edges exist.

The centre $m$ of an odd twin pair is even. Except for $m=4$, it satisfies
$m\equiv0\pmod6$. For a finite set $S$ of distinct odd primes, put

$$M=2\prod_{p\in S}p,\qquad
\mathcal W_S=\{m\bmod M:m\text{ even},\ m\not\equiv\pm1\pmod p\ (p\in S)\}.$$

The Chinese remainder theorem gives exactly

$$|\mathcal W_S|=\prod_{p\in S}(p-2). \tag{1}$$

This is a necessary twin-prime condition once both endpoints exceed all primes
in $S$. At the finite boundary an endpoint equal to the sieving prime is not
composite and must not be discarded for that reason.

**An infinite failure family for every finite wheel.** Choose distinct odd
primes $q,r\notin S$. Solve

$$m\equiv0\pmod M,\qquad m\equiv1\pmod q,\qquad m\equiv-1\pmod r.$$

Every sufficiently large positive member of this arithmetic progression lies
in $\mathcal W_S$, while $q\mid m-1$ and $r\mid m+1$ are proper factors.
For example, with $S=\{3,5,7\}$ the compiler gives
$m=9450+30030t$, $t\ge0$, with proper factors 11 and 13. No fixed finite wheel
turns survival into primality, even for one of the two neighbours.

## 3. The additive maps and their radical information

For even $m\ge2$, both

$$ (m-1)+2=m+1,\qquad (m^2-1)+1=m^2 \tag{2}$$

are primitive positive ABC packets. Their radicals are, respectively,

$$R_1=\operatorname{rad}(2(m^2-1)),\qquad
R_2=\operatorname{rad}(m)\operatorname{rad}(m^2-1). \tag{3}$$

If $m-1,m+1$ are prime, then $R_1=2(m^2-1)>m+1$ and
$R_2=\operatorname{rad}(m)(m^2-1)>m^2$. Thus these twin-prime packets have
ABC quality $\log c/\log\operatorname{rad}(abc)<1$; they are not examples
of exceptionally small radicals. Conversely the same packet families include
$675+1=676$ with radical 390 and $243+2=245$ with radical 210.

Under ABC, (2) gives lower radical bounds
$R_1\ge((m+1)/K_\epsilon)^{1/(1+\epsilon)}$ and
$R_2\ge(m^2/K_\epsilon)^{1/(1+\epsilon)}$. Neither inequality asserts that
any term is prime. The map to a radical discards valuations, allocation of
prime factors to terms, and height: $1+2=3$ and $1+8=9$ both have radical 6.

## 4. Golden addition gives a global Fibonacci twin exclusion

Let $F_0=0,F_1=1,L_0=2,L_1=1$, with the usual common recurrence. The signed
golden norm is

$$L_n^2-5F_n^2=4(-1)^n. \tag{4}$$

The addition identity, with $F_{-j}=(-1)^{j+1}F_j$, is

$$F_aL_b=F_{a+b}+(-1)^bF_{a-b}. \tag{5}$$

For odd $n\ge9$, take $a=(n+3)/2$ and $b=(n-3)/2$. Then

$$F_aL_b=F_n+2(-1)^b,\qquad
F_bL_a=F_n-2(-1)^b. \tag{6}$$

Every displayed factor exceeds 1 because $a\ge6$ and $b\ge3$. Therefore
both $F_n-2$ and $F_n+2$ are composite for every odd $n\ge9$.

For completeness, $F_d\mid F_n$ when $d\mid n$, by the Fibonacci addition
formula and induction on $n/d$. If $n$ is composite and $n\ne4$, a proper
divisor $d\ge3$ exists, except the irrelevant indices below 4. Then
$1<F_d<F_n$, so a prime Fibonacci value has prime index or index 4. Checking
indices 3, 4, 5, and 7 and using (6) proves the exact, classical conclusion:

**The prime Fibonacci values that have a twin-prime neighbour are precisely
$3,5,13$.**

This is an elementary factorization result, not a new priority claim. The
surviving prime-index condition is not sufficient for primality:
$F_{19}=37\cdot113$ and $L_{23}=139\cdot461$.

## 5. What ABC would imply for the golden norm

Since $L_n=F_n+2F_{n-1}$ and consecutive Fibonacci numbers are coprime,

$$g_n=\gcd(F_n,L_n)=
\begin{cases}2,&3\mid n,\\1,&3\nmid n.\end{cases} \tag{7}$$

Here the divisibility by 2 follows from the exact period-three recurrence
modulo 2. In particular, prime indices $p>3$ give coprime $F_p,L_p$;
$p=3$ gives $(2,4)$, not a coprime pair.

For $n\ge1$, divide (4) by $g_n^2$ and place the smaller square first:

$$a_n=\frac{\min(L_n^2,5F_n^2)}{g_n^2},\quad
b_n=\frac4{g_n^2},\quad
c_n=\frac{\max(L_n^2,5F_n^2)}{g_n^2}.$$

This is a positive primitive ABC packet. Indeed the two square terms have
gcd $g_n^2$ (5 never divides $L_n$, by (4) modulo 5), and their difference
is 4. Put $H_n=F_nL_n$. The identity for $L_n$ gives
$1\le L_n/F_n\le3$, hence

$$c_n\ge\frac54F_n^2\ge\frac5{12}H_n,
\qquad \operatorname{rad}(a_nb_nc_n)\le10\operatorname{rad}(H_n).$$

Consequently ABC implies, for every $\epsilon>0$,

$$\operatorname{rad}(H_n)\ge\frac1{10}
\left(\frac{5H_n}{12K_\epsilon}\right)^{1/(1+\epsilon)}.$$

As $H_n\to\infty$ and $\operatorname{rad}(H_n)\le H_n$, this proves the
**conditional** limit

$$\lim_{n\to\infty}
\frac{\log\operatorname{rad}(F_nL_n)}{\log(F_nL_n)}=1. \tag{8}$$

The fixed factors 5 and $b_n\in\{1,4\}$ affect only constants. This is a
Pell-family instance of the inherited radical-growth mechanism. It gives
neither a count of Fibonacci/Lucas primes nor a probability of simultaneous
primality. The checker validates the algebraic premises, not ABC itself.

## 6. A genuine geometric transfer, with a two-colour cap

Difference of squares gives the exact rational unit direction

$$v(m,d)=\left(\frac{m^2-d^2}{m^2+d^2},
                    \frac{2md}{m^2+d^2}\right),\qquad 0<d<m.$$

Prime restrictions on $m,d$ leave this inside rational-coordinate geometry.
That carrier cannot produce a six- or seven-chromatic unit-distance witness.
Here is a constructive proof independent of any ABC premise.

Write any rational unit vector as $(a/c,b/c)$ with integer $c>0$ and
$\gcd(a,b,c)=1$. The identity $a^2+b^2=c^2$ modulo 4 forces $c$ odd and
exactly one of $a,b$ odd. Thus each coordinate is in
$\mathbb Z_{(2)}=\{u/v:v\text{ odd}\}$ and the sum of the two coordinates
reduces to 1 modulo 2.

For a finite unit-distance graph on rational points, retain one named point
$o$ per connected component. Every difference $p-o$ is a sum of rational
unit steps, so its coordinates lie in $\mathbb Z_{(2)}$. Define

$$\chi(p)=(p_x-o_x)+(p_y-o_y)\pmod2.$$

Every unit edge changes $\chi$. This proves two-colourability and supplies
an exact arithmetic colouring algorithm; traversal only chooses component
origins, not a bipartite-colouring oracle. Arbitrary rational translations
are allowed. The origin sidecar matters: points $(1/2,0)$ and $(3/2,0)$ have
no direct odd-denominator coordinate reduction, but their relative colours
are 0 and 1.

The API either constructs all exact unit edges or reauthenticates each edge
of a supplied graph with the complete original vertex set. Its output can
be checked independently by the signed-boundary compiler. The cap does not
apply to arbitrary real coordinates, which already admit an equilateral
unit triangle. This is a carrier restriction, not a plane-colouring bound.

## 7. Connection contract and the surviving target

| Source and map | Preserved predicate | Information lost or required |
|---|---|---|
| Twin centre $m$ to (2) | exact primitive additive identity | primality is not implied by its radical |
| Golden coefficients to (4) | signed norm and exact integer values | a norm alone does not determine factors |
| Centre to residues modulo $M$ | finite divisibility exclusions | unseen prime factors; the CRT hostile supplies both |
| Rational parameters to $v(m,d)$ | exact unit length | component origin is needed for arithmetic parity |
| Prime-gap graph | actual difference-two relation | it has nonedges, so is not a tournament |

The finite sieve and golden factorization are different operations: the
former supplies local exclusions, while the latter supplies a global exact
product. Neither is a count of labelled tournaments. The useful next plane
target is an exact irrational coordinate field together with its norm and
valuation sidecars; adding primality filters to the rational carrier cannot
overcome its proved two-colour cap.

## 8. Reproduction and finite scope

From the repository root:

```text
python -B 04-computation/experiments/abc_difference_families_20261011c.py
python -B -O 04-computation/experiments/abc_difference_families_20261011c.py
```

The declared controls are even centres 2 through 2000; Fibonacci/Lucas
identities through index 500 (trial-factor radical premises only through 26);
five complete finite wheels and 64 members of each composite-survivor ray;
and 91 rational rotated squares with 364 exact unit edges. Exact type,
nonprimitive-packet, translated-origin, and forged-unit-edge hostiles are
included. These finite checks supplement the all-parameter proofs above.

Normal and optimized runs both pass **5,685 checks** and reproduce the
[saved output](abc_difference_families_20261011c.out). LF-normalized SHA-256:
`0aa19a9fc3c6a24fc290b16a1870688a7cd34d76ba3ceec1e76a63290d8b8348`.
