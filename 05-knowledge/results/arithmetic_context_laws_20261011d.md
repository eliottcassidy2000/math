# Arithmetic laws need a sampling context

**Status:** PROVED elementary finite identities, two context-transfer lemmas,
and their stated limits; CITED classical average orders and prime theorems;
CONJECTURAL binary prime-representation asymptotics; FINITE-EXACT checks.
The new instrument preserves local divisibility data and its exact weights.
It does not solve a prime-representation conjecture.

## 1. Inheritance and the three prime triples

The live board is **sampling measure / multiplicity / local congruences /
global prime correlation / retained arithmetic roles**. The nearest proved
mechanism is CRT factorization of a finite residue space. Its hostile is
the [infinite both-composite family surviving every fixed finite prime
wheel](abc_difference_families_20261011c.md). The corrected near miss is
inferring primality independence from coprimality; the historical
Fibonacci/Lucas census was repaired in the preceding session. The less-used
sidecar is the multiplicity of each boundary/divisibility state.

[THM-2447, twin-gap local prime harmonics](../../01-canon/theorems/THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum.md)
already separates an exact local residue count from a fitted global spectrum.
[Arithmetic braids, SW3](arithmetic_braids_20260917_divisors.md) already keeps
Chen's prime-or-semiprime union distinct from its individual outcomes.
The [analytic-number-theory source sidecar](../reference/CORE-PAPERS-ANALYTIC-NUMBER-THEORY.md)
is inherited; none of its unrelated preprint claims is needed here.

These three sets encode different predicates:

| Numbers | Exact inherited role | What the set alone does not retain |
|---|---|---|
| $\{2,3,11\}$ | Hecke prime 2, reduction characteristic 3, conductor 11 in the specified oldform example | the curve, basis, lattice and operator |
| $\{3,5,13\}$ | exactly the prime Fibonacci **values** with a twin-prime neighbour | their indices $4,5,7$ and the factorization proving exclusion of all later odd indices |
| $\{3,5,7\}$ | the unique overlapping pair of twin-prime edges | the difference-two relation and its mod-3 obstruction |

The first row is the [level-22 oldform bridge](collatz_level22_oldform_bridge_20261005.md):
the conductor remains 11, while the two-dimensional oldspace has ambient
level 22. In its specified basis,

$$A=\begin{pmatrix}-2&1\\-2&0\end{pmatrix},\qquad
A^4=-4I,\qquad
A\bmod3=\begin{pmatrix}1&1\\1&0\end{pmatrix}.$$

The last matrix has order 8 and golden characteristic polynomial
$X^2-X-1$ over $\mathbb F_3$. This is an actual operator reduction, not
an identification of three prime triples. The other rows have elementary
all-index proofs in [the difference-family packet](abc_difference_families_20261011c.md).
Their shared entries do not transfer their predicates. The golden
discriminant 5 and the Fibonacci value $13=F_7$ are also different roles.

A second, precisely specified comparison uses the same ring for all six
primes: $R_p=\mathbb F_p[t]/(t^2-t-1)$ and its unit $t$. The inherited
[golden prime clocks](golden_prime_clocks_20261003.md) give the quotient
mechanism; exact coordinate multiplication supplies the following controls.

| $p$ | Ring type | Order of $t$ | Scalar-root orders, when present |
|---:|---|---:|---|
| 2 | inert, $\mathbb F_4$ | 3 | none |
| 3 | inert, $\mathbb F_9$ | 8 | none |
| 5 | ramified, $\mathbb F_5[\eta]/(\eta^2)$ | 20 | root 3 has order 4 |
| 7 | inert, $\mathbb F_{49}$ | 16 | none |
| 11 | split, $\mathbb F_{11}\times\mathbb F_{11}$ | 10 | roots 4 and 8 have orders 5 and 10 |
| 13 | inert, $\mathbb F_{169}$ | 28 | none |

Here irreducible/repeated/split factorization decides the ring type, and
successive exact powers decide the listed finite orders. Thus the two sets
$\{3,5,7\}$ and $\{3,5,13\}$ have the same splitting-type pattern but
different clocks. At 5, forgetting the nonzero nilpotent $t-3$, whose square
is zero, loses the factor 5 in the order. No such quotient determines
whether a Fibonacci value is prime. The table is classical finite-ring
arithmetic, not a new number-pattern classification.

## 2. Four familiar scales, correctly sampled

Let $U_N$ be uniform on $\{1,\ldots,N\}$ and let all logarithms be natural.
Write $\tau(n)$ for the number of positive divisors, $\omega(n)$ for distinct
prime factors, and $\Omega(n)$ for prime factors with multiplicity. The exact
finite definitions and standard asymptotic meanings are:

| Statistic of $U_N$ | Exact expression | Asymptotic statement |
|---|---|---|
| $\mathbb E\tau(U_N)$ | $N^{-1}\sum_{d\le N}\lfloor N/d\rfloor$ | $\log N+2\gamma-1+O(N^{-1/2})$ |
| $\Pr(U_N\text{ prime})$ | $\pi(N)/N$ | $\sim1/\log N$ |
| $\mathbb E\omega(U_N)$ | $N^{-1}\sum_{p\le N}\lfloor N/p\rfloor$ | $\log\log N+O(1)$ |
| $\mathbb E\Omega(U_N)$ | $N^{-1}\sum_{p^a\le N}\lfloor N/p^a\rfloor$ | $\log\log N+O(1)$ |
| $\Pr(U_N\text{ squarefree})$ | $N^{-1}\sum_{d\le\sqrt N}\mu(d)\lfloor N/d^2\rfloor$ | $6/\pi^2+O(N^{-1/2})$ |

The divisor estimate follows from the hyperbola identity; see
[Edinburgh number-theory notes, Proposition 9.13](https://webhomes.maths.ed.ac.uk/~chris/NTh/Number_Theory_Notes2012.pdf#page=42).
The prime row is the prime number theorem, with its meaning stated in the
[Oxford analytic-number-theory course description](https://courses.maths.ox.ac.uk/course/section.php?id=13623).
The two factor means and the squarefree limit are standard;
[Tao's elementary multiplicative-number-theory notes, Sections 3--4](https://terrytao.wordpress.com/2014/11/23/254a-notes-1-elementary-multiplicative-number-theory/)
also distinguish these statistics from typical values.

The last identity comes from $1_{\rm sf}(n)=\sum_{d^2\mid n}\mu(d)$.
Replacing floors by $N/d^2$ contributes $O(\sqrt N)$ to the count; the
discarded tail of $\sum\mu(d)/d^2$ has the same bound. Its Euler product is
$\prod_p(1-p^{-2})=1/\zeta(2)=6/\pi^2$.

An expectation is not the value at a fixed integer. In particular, for
almost all integers up to $N$, the classical normal-order statement is
$\tau(n)=(\log N)^{\log2+o(1)}$, rather than $\log N$; this follows from
$2^{\omega(n)}\le\tau(n)\le2^{\Omega(n)}$ and the Hardy--Ramanujan theorem.
The larger mean is sensitive to the tail. Under the different experiment
“choose uniformly from the primes up to $N$”, the four statistics become
$\tau=2$, $\omega=\Omega=1$, and squarefree indicator 1 exactly.

## 3. A useful conditioning law: $1/p$ becomes $1/(p+1)$

Fix a finite set $S$ of primes, and choose a residue uniformly modulo
$M_S=\prod_{p\in S}p^2$. Condition on no $p^2$ dividing it for $p\in S$.
For each $p$, there are $p^2-1$ allowed residues, of which $p-1$ are
divisible by $p$. Therefore

$$\Pr(p\mid n\mid p^2\nmid n)=\frac1{p+1}. \tag{1}$$

CRT makes these flags independent **in this finite conditional sample**.
More precisely, for bits $e_p\in\{0,1\}$,

$$\Pr((1_{p\mid n})_{p\in S}=(e_p))
=\prod_{p\in S}\frac{1\text{ if }e_p=1\text{, else }p}{p+1}. \tag{2}$$

The same joint limit holds for each fixed $S$ when $U_N$ is conditioned
to be globally squarefree. To justify the passage, first impose square
exclusions for all primes through a fixed $Y\ge\max S$ and apply CRT.
The upper density of integers with any omitted prime square divisor is
at most $\sum_{p>Y}p^{-2}$, which tends to zero; the squarefree density
stays positive. Let $N\to\infty$ first and then $Y\to\infty$.
No independence claim for a set of primes growing with $N$ is needed.

This is an exact example of why admissible states alone are insufficient:
the divisibility flag has two values, but their weights are $1:p$, not
$1:1$. Products are legitimate before normalization because CRT separates
the coordinates. Conditioning changes the weights without changing the
flag alphabet. This is the arithmetic counterpart of retaining boundary
multiplicity when gluing finite objects; arbitrary correlated objects need
their own factorization proof.

## 4. Additive prime statements: count, probability, theorem or conjecture

Define the **ordered** Goldbach count

$$R(N)=\#\{(p,q):p,q\text{ prime},\ p+q=N\}.$$

If $X$ is uniform on $\{1,\ldots,N-1\}$, the probability that both
$X,N-X$ are prime is exactly $R(N)/(N-1)$. A representation count and
this probability differ by a factor of order $N$.

- **CITED theorem:** Vinogradov's sufficiently-large odd three-prime theorem
  concerns $p+q+r$. Helfgott proves every odd integer greater than 5 is a
  sum of three primes; see [the primary paper](https://arxiv.org/abs/1312.7748).
- **CONJECTURE:** every odd integer greater than 5 is $p+2q$ with both
  $p,q$ prime. This is Lemoine's conjecture, not the preceding theorem;
  [Sun, Remark 1.4](https://arxiv.org/html/0803.3737v7) records the formulation.
- **CONJECTURE:** every even $N\ge4$ has $R(N)>0$ (binary Goldbach).
- **CONJECTURAL asymptotic:** with
  $C_2=\prod_{p>2}(1-(p-1)^{-2})$,

  $$R(N)\sim2C_2\prod_{\substack{p\mid N\\p>2}}
       \frac{p-1}{p-2}\;\frac{N}{(\log N)^2}
       \qquad(N\to\infty\text{ through even integers}). \tag{3}$$

  The Hardy--Littlewood count includes the target's local divisor factor.
  Dropping that factor does not give a uniform formula over even targets.
  An unordered convention has half the leading term.
- **CITED theorem:** Chen gives a prime plus a number with at most two
  prime factors counted with multiplicity, for every sufficiently large
  even target. It does not determine which prime/semiprime outcome occurs.
  [Li's primary research paper, introduction and (1)--(3)](https://arxiv.org/html/2405.05727v4)
  records Chen's classical lower bound for this count and the local factor
  used in (3). We use that classical statement, not the paper's new numerical
  constant.

The expression $C_2(\log\log N+1)^2/(\log N)^2$ is not Chen's theorem.
Landau's classical fixed-factor-count law gives
$\#\{n\le N:\Omega(n)=2\}\sim N\log\log N/\log N$;
see [the primary almost-prime study](https://arxiv.org/abs/1401.2694).
Multiplying two such marginal scales suggests a squared logarithm for two
almost-prime events only as an **independence heuristic**. One prime and
one almost-prime marginal instead suggest one $\log\log N$ factor. Neither
product proves the required additive correlation, and the additive $+1$
does not acquire a controlled lower-order meaning from leading asymptotics.

## 5. An exact local transfer between two different prime problems

For integer affine forms $L_i(x)=a_ix+b_i$, let $\nu_p$ count residues
where at least one form vanishes modulo $p$. On a squarefree wheel
$W=\prod_{p\in S}p$, CRT gives exactly

$$\#\{x\bmod W:p\nmid\prod_iL_i(x)\ (p\in S)\}
=\prod_{p\in S}(p-\nu_p). \tag{4}$$

Relative to the unit density of $k$ generic independent integers, the
normalized factor is $(1-\nu_p/p)/(1-1/p)^k$. This definition does not
assert independence of the actual forms.

Take an odd positive target $N$. Compare Lemoine's two forms
$q,N-2q$ with Goldbach's forms $x,2N-x$. For odd $p$, the bad sets are
$\{0,N/2\}$ and $\{0,2N\}$ respectively. Multiplication by 4 maps the
first to the second, including the coincident-root case $p\mid N$.
At $p=2$, both contexts allow exactly the odd input residue. Hence

$$\nu_p=\begin{cases}1,&p=2\text{ or }p\mid N,\\2,&p>2,\ p\nmid N,
\end{cases} \tag{5}$$

and all their normalized local factors agree.

**Explicit finite-wheel compiler.** Suppose $W=2Q$, with $Q$ odd and
squarefree. The map

$$x=(Q+4)q\pmod{2Q} \tag{6}$$

is a bijection of the entire residue space because $Q+4$ is a unit modulo
$2Q$. It equals $4q$ modulo each odd prime divisor and $q$ modulo 2.
Consequently it maps the complete Lemoine survivor set onto the complete
Goldbach survivor set. The inverse is multiplication by $(Q+4)^{-1}$.
Refinement gives compatible residue maps, but their least nonnegative
integer representatives can change: $q=7$ maps to 13 modulo 30 and 133
modulo 210. There is no single integer multiplier satisfying the prescribed
rules at every prime, and no primality-preserving integer map is claimed.

The archimedean intervals also differ: Goldbach at target $2N$ samples
$1\le x\le2N-1$, while Lemoine samples
$1\le q\le\lfloor(N-2)/2\rfloor$. Their lengths have ratio tending to 4.
Already $N=13$ has two actual Lemoine representations and five ordered
Goldbach representations of 26. Identical local factors do not remove the
interval or global prime-correlation problem.

There are two additional sharp boundaries. On wheel 30, the two coordinates
$x,N-x$ each have the same unit marginal $8/30$ for every $N$, but the joint
survivor counts are 8 for $N=30$ and 3 for $N=32$. Marginals lose the target
correlation. Also sieving out a prime divisor rejects the exceptional value
equal to that prime: the genuine representation $7=3+2\cdot2$ is absent
from the parity wheel's surviving $q$ residues. Exact finite counts must
restore these bounded prime-equals-modulus exceptions.

## 6. What can be reused

The compiler first binds the affine forms, target, tested primes and sample
space, then returns exact residue counts and an authenticated transport.
It can transfer a local sieve computation between two contexts. To transfer
an actual prime theorem one would additionally need interval control and
the joint prime correlation; these are not present in the packet.

The productive connection is therefore a division of responsibilities:
local congruences determine an exact finite weight, whereas a global theorem
must justify how actual integer values occupy that weighted space. Shared
numbers such as 2, 3, 5, 11 and 13 do not supply the missing map by themselves.

## 7. Reproduction

[Script](../../04-computation/experiments/arithmetic_context_laws_20261011d.py)
and [saved output](arithmetic_context_laws_20261011d.out):

```text
python -B -X utf8 04-computation/experiments/arithmetic_context_laws_20261011d.py
python -B -O -X utf8 04-computation/experiments/arithmetic_context_laws_20261011d.py
```

The checks compare independent summatory formulas for every $N\le256$ and
at $100,1000,10000$; exhaust four squarefree residue wheels through modulus
44100; and verify all residues of five pair wheels through modulus 2310
for every odd $N\le65$ (165 contexts, 84,414 inverse controls).
Eight small actual prime-representation counts are receipt checks, not
asymptotic evidence. Exact type, nonprime modulus, wrong parity, repeated
prime, zero-form and prime-equals-modulus controls are retained.

Normal and optimized runs reproduce **86,338 checks** and the saved output.
LF-normalized SHA-256: `a5f9f1408ead7c8cf2fd6dda32d627e695130021d99702751e789e4031ca3e93`.
