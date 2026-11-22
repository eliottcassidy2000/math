# Affine sieve on a positive tournament grammar

**PROVED using established expander and sieve theorems.** A fixed three-letter
physical port grammar gives infinitely many distinct, strong tournament
Hamiltonian-path counts with a bounded number of prime factors, with
multiplicity counted. This also holds in every fixed reduced residue
class of modulus coprime to 6, with one uniform bound on the number of
prime factors. No RH or GRH input is needed for this theorem. It does
not give prime counts or an exact selector for arbitrary targets.

The new ingredient is a spectral gap for the **positive** grammar,
not just transitivity with modulus-dependent positive inverse powers.
The closest input is [the rooted migration transfer](../chanpak_bank_fusion_20261117/migration/PROOF.md);
[odd-tag steering](../chanpak_branching_inverse_20261122/ODD-TAGS.md)
already gives joint residues without the uniform length or distribution control.

## 1. A fixed physical alphabet with one certified strong finish

Use the marked column state $(h,R)$ of the previous packet. Here $h$ is
the ordinary Hamilton count; $R$ counts a nonempty path ending at the
marked terminal together with its complementary path, possibly empty.
The transitive pair marked at its last vertex starts at $(1,2)$, with
exactly one full path ending at the root and none starting there. An
inherited leading-zero migration preserves these last two guards. For
the following widths and binary addresses its transfer matrices are
$$
 A=M(3,2)=\begin{pmatrix}3&3\\6&8\end{pmatrix},\qquad
 B=M(5,8)=\begin{pmatrix}9&9\\30&32\end{pmatrix},\qquad
 C=M(4,2)=DA,\quad D=\operatorname{diag}(1,2).
$$
The first three matrices are physical, with determinants $6,18,12$.
$D$ by itself is used only in the analysis. The leading-one finish at
width one and address one gives $N=h+2R$, as an ordinary count of a
strong tournament. There are at most $6n+4$ vertices after $n$ letters
and this finish. The largest column sum of the alphabet is 41, and
each step at least triples $h$. Hence $3^n\le N\le6\cdot41^n$.
The finish count is odd by Rédei. Modulo 3, after the first letter
and thereafter, $h=0$ and $R$ is nonzero: each second row reduces to
$(0,1)$ or $(0,2)$. It follows that every finished count at positive
depth is coprime to 6. This deterministic handling of a bad prime is
essential; it cannot simply be discarded from the sieve.

## 2. Uniform expansion from the relative two-letter motion

Let $q>1$ be coprime to 6 and let $\Omega_q$ be the primitive vectors in
$(\mathbb Z/q)^2$. The letters act by permutations on this set. Let
$P$ be the average of the three permutation operators in counting
measure. We show there are constants $\theta<1$ and $K>0$, independent
of $q$, such that on the mean-zero subspace
$$
 \|P^n\|\le\theta^{\lfloor n/2\rfloor}.
$$

Write $T=P^2$ and let $u,v$ range over the nine chronological two-letter
products (with multiplicities). For the unitary permutation action
$\rho$ one has the exact identity
$$
 \|f\|^2-\|T^* f\|^2
  =\frac{1}{162}\sum_{u,v}\|f-\rho(uv^{-1})f\|^2.
$$
This is just the pairwise squared difference identity for an average
of equal-length vectors. Every relative move is bounded by this
energy. If a fixed group word $g$ has length $\ell$ in relative moves,
the triangle inequality and Cauchy–Schwarz bound
$\|f-\rho(g)f\|^2$ by $162\ell^2$ times the same energy.
Inverses appear in this comparison, not in the growing tournaments.

Comparing $BA,CA,AB,AC$ with $A^2$ gives respectively
$G=BA^{-1},D,AGA^{-1},ADA^{-1}$. Direct multiplication gives
$$
 G=\begin{pmatrix}3&0\\8&1\end{pmatrix},\qquad
 J=DGD^{-1}G^{-1}=E_{21}(8/3).
$$
Thus both $J$ and $AJA^{-1}$ are bounded words in relative moves. In
the rational basis $Q=\begin{pmatrix}3&0\\8&1\end{pmatrix}$ these become
$$
 Q^{-1}JQ=E_{21}(8),\qquad
 Q^{-1}AJA^{-1}Q=E_{12}(-4/3).
$$
Taking a sixth power of the second and allowing inverses gives the
conjugates by $Q$ of the two integral shears $E_{12}(8),E_{21}(8)$ and
their inverses. Their subgroup is Zariski dense in $\mathrm{SL}_2$:
the powers of either shear are dense in its one-parameter unipotent
subgroup, and those two subgroups generate the algebraic group.
Modulo $q$ they generate $\mathrm{SL}_2(\mathbb Z/q)$, since 8 is a
unit, and arbitrary elementary shears generate that group (reduce
over each prime-power ring, then combine with CRT). Its action on
$\Omega_q$ is transitive: a primitive vector can be completed to a
determinant-one matrix. The basis $Q$ is invertible modulo $q$.

[Bourgain–Varjú, Theorem 1](https://arxiv.org/abs/1006.3365) gives a
uniform gap for the symmetric integral shears on all their congruence
quotients. It gives the same lower bound for the Dirichlet energy
on this transitive permutation representation, by identifying
functions on a coset space with a subspace of functions on the group.
Comparison with the displayed exact energy for $T^*$ gives
$\|T^*\|_{L^2_0}\le\theta<1$. The norm of $T$ on the invariant
mean-zero space is the same, and the asserted bound for $P^n$ follows.
There is no use of a mixing conjecture for this positive language.

In particular, starting from any primitive residue state and using
uniform independent physical letters,
$$
 \left|\Pr[x_n\in U]-\frac{|U|}{|\Omega_q|}\right|
 \le \sqrt{|U|}\,\theta^{\lfloor n/2\rfloor}
$$
for every $U\subseteq\Omega_q$. Apply Cauchy–Schwarz to the evolved
point mass minus the uniform distribution. Taking $U$ a singleton
and $|\Omega_q|\le q^2$ shows that every primitive vector is reached
by a positive word of any given length $n\ge K\log q$ with a
fixed alphabet, before the finish. This also shows that any prescribed
ordinary residue $a\bmod q$ can be obtained by a strong tournament on
$O(\log q)$ vertices, since $(a-2,1)$ is primitive and the finish is
$h+2R$. This short realization does not prescribe the actual value.

## 3. Sieve the actual integer outputs, not an ambient progression

Count words, including multiplicity if counts collide. Put $X_n=3^n$.
For squarefree $d$ coprime to 6, measure $d\mid N$ using the same walk
modulo $d$. At a prime $p$, there are $p-1$ nonzero vectors in the
kernel of $h+2R$ among $p^2-1$ primitive vectors. Thus the exact local
density is $\beta(p)=1/(p+1)$, and CRT gives
$\beta(d)=\prod_{p\mid d}1/(p+1)$. The mixing estimate yields
$$
 \#\{w:\ |w|=n,\ d\mid N(w)\}
 =X_n\beta(d)+O(X_n d e^{-\eta n})
$$
for some $\eta>0$; use $|U|\le d^2$. Choosing
$D=\exp(\delta n)$ with $0<\delta<\min(\eta/4,\log3)$ makes the sum
of errors for $d\le D$ power-saving in $X_n$. The sieve dimension is
one since $\beta(p)=1/p+O(p^{-2})$, and the standard Mertens estimate
for the sum of $\log p/p$ gives the dimension-one condition.
Apply the [Bourgain–Gamburd–Sarnak combinatorial sieve, Theorem 3.1](https://www.ias.edu/sites/default/files/math/Bourgain_Gamburd_Sarnak_Affine.pdf)
to this nonnegative sequence with primes 2 and 3 excluded, and take
$z=D^{1/10}$. It gives $\gg 3^n/n$ words whose $N$ has no prime factor
at most $z$ (2 and 3 were already excluded). The bound $N\le 6\cdot41^n$
gives an absolute bound on $\Omega(N)$, the number of prime factors
with multiplicity. The lower bound $N\ge3^n$ proves infinitely many
distinct outputs as $n$ varies. The word estimate is not an estimate
for the number of distinct counts at one depth.

For the stated progression refinement fix $q$ coprime to 6 and any
primitive vector $b\bmod q$ for which $b_1+2b_2$ is a unit. Keep only
words ending at $b$ modulo $q$. With squarefree $d$ coprime to $6q$,
mix modulo $qd$: the main term is
$3^n\beta(d)/|\Omega_q|$ with error $O(3^n qd e^{-\eta n})$.
For fixed $q$ the same exponent $\delta$ works, and the excluded
prime set is the primes dividing $6q$; none divides the selected
counts. This proves the same absolute upper bound on $\Omega(N)$
in every such phase. For a desired unit class $a\bmod q$, take
$b=(a-2,1)$.

The argument also permits prescribed nonzero valuations at finitely
many primes $p\ge5$: encode each exact valuation in $q$ using one
extra digit, and sieve only at primes not dividing $6q$. The same
bound holds for the number of **remaining** prime factors. For example
the phase $N\equiv7\bmod49$ supplies infinitely many distinct strong
counts $7M$ with $7\nmid M$ and with a uniform bound on $\Omega(M)$.
It does not yet show $M$ prime, much less realize all required $7p$.

## 4. RH clarifies why Type I is not enough

[The Lean and Mathlib joint statement](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/)
certifies the formal classical RH proof; [the independent MIT verification](https://math.mit.edu/~drew/grh/)
also identifies the Dirichlet-GRH statements checked. This distinction
does not assert RH for arbitrary other $L$-functions. The link from the
[new proof's announcement](https://alpo.ge/remark.pdf) to this theorem
is a concrete proof design: turn signed comparisons into a positive
quadratic expression, prove a gap, then use enough iterations to beat
an arithmetic remainder. Here the operator is an exact finite
two-step variance, not a theta lift or a zeta calculation.

There is also a direct consequence of classical RH that explains
the missing step. With $\lambda$ the completely multiplicative
Liouville sign, compare the nonnegative arithmetic weights
$1+\lambda(n)$ and $1-\lambda(n)$ on integers coprime to 6.
One has zero weight on all primes; the other has weight two on
each prime greater than three. Yet for $(d,6)=1$ both give
$$
 \sum_{\substack{n\le Y\\ (n,6)=1\\d\mid n}}(1\pm\lambda(n))
 =\frac{Y}{3d}+O_\epsilon((Y/d)^{1/2+\epsilon}),\qquad d\le Y.
$$
Indeed the classical Littlewood estimate
$\sum_{n\le t}\mu(n)=O_\epsilon(t^{1/2+\epsilon})$ follows from RH
([Soundararajan, introduction](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/135214/blms0213.pdf?sequence=1)).
Since $\lambda(n)=\sum_{r^2\mid n}\mu(n/r^2)$, the same bound holds
for its summatory function. Removing multiples of 2 and 3 by
inclusion–exclusion, then using complete multiplicativity to remove
$d$, gives the formula. Summing the errors over $d\le Y^{1-\tau}$
still gives a power saving for any fixed $\tau>0$ (take
$0<\epsilon<1/2$). These are arithmetic examples, not tournament
outputs. They show why even very strong marginal divisibility
estimates do not detect primes by themselves.

**Next obligation:** obtain an additional, genuinely orbit-sensitive
estimate (for example, parity-sensitive bilinear cancellation), or
find a deterministic exact selector with improving height. Even an
infinitude of prime-valued words would not by itself hit every
required prime in the Chan–Pak lane.

Finite controls of the new physical transfers and algebraic identities: run `python3 check_affine_sieve.py`.
