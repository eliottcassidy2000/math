# Parity of the small-prime factors on physical tournament words

**PROVED.** For the fixed positive tournament alphabet and strong finish
of the [affine-sieve packet](../chanpak_affine_sieve_metagraph_20261122/AFFINE-SIEVE.md),
the parity of the total multiplicity of prime factors up to a fixed
positive multiple of word length is asymptotically balanced. There are
also two simultaneous physical ordinary counts from the same word whose
small-factor parities become jointly uniform, even with the parity of
their scaled gap included. This is not cancellation of the full
Liouville function and does not show that any of these outputs is prime.

## 1. Two actual counts and one gap

Let $(h,R)$ be the positive column state after $n$ independent uniform
physical letters $M(3,2),M(5,8),M(4,2)$ of the inherited fixed grammar.
At positive depth $h$ is odd, $R$ is even, $h=0\bmod3$, and
$R\ne0\bmod3$. The parity assertions follow from the seed $(1,2)$ and
the three integer matrices; the mod-3 last-row guard is in the
affine-sieve proof. Each state has an actual marked terminal with a
unique full path ending there and none starting there.

Apply one extra leading-zero migration $M(4,6)$ with matrix
$$
 D=\begin{pmatrix}5&7\\10&16\end{pmatrix}.
$$
It preserves those physical terminal guards. Use that rooted tournament
as the right flank of either of the two $(q,t)=(2,1)$ middle gadgets in
[the unit-switch packet](../chanpak_affine_sieve_metagraph_20261122/UNIT-SWITCHES.md),
and use the marked transitive pair as the left flank (after taking the
converse as prescribed). The four-channel graft formula is valid for
**every** rooted flank with the signed endpoint sums zero; its
specialization to invariant lines was only needed for a constant
difference. The new right flank is not assumed to be on that line.
The left ordinary count and cover count are $(1,2)$, so the middle
signatures give the rows $(1313,2447)$ and $(1319,2445)$ on the right
flank after $D$. Consequently the two literal strong tournament counts
and their scaled difference are
$$
 L_0=31035h+48343R,\qquad L_1=31045h+48353R,\qquad
 G=(L_1-L_0)/10=h+R.
$$
Both middle gadgets admit the stipulated internal path from the second
mark to the first. Hence the strongness part of the old graft proof
applies with the same marked guards. All three positive integers are
coprime to 6: for 2 use the parity of $h,R$; modulo 3 the coefficients
of $R$ in these three rows are all 1 or 2. The third linear form is
an integer obtained from the two counts, **not** asserted to be itself
a tournament count. The gap is generally unbounded.

## 2. Amplification of congruence information to truncated Liouville

For $z\ge5$ set
$$
 \lambda_{\le z}(N)=(-1)^{\sum_{p\le z}v_p(N)}.
$$
There is an absolute $c>0$ such that, with $z=cn$, for every nonempty
subset $I$ of $\{L_0,L_1,G\}$,
$$
 \mathbb E\prod_{L\in I}\lambda_{\le z}(L(h,R))
 \sim \frac{C_I}{(\log z)^{2|I|}},\qquad C_I>0.
$$
The expectation is over the $3^n$ input words with multiplicity, not
over distinct tournaments or distinct count values. The same result for
the original strong finish $h+2R$ holds with $|I|=1$.

**Proof.** The inherited affine-sieve gap gives, on primitive vectors
modulo $q$ coprime to 6, total-variation distance at most
$O(qe^{-\delta n})$ from uniform, uniformly in $q$ and from the
primitive seed $(1,2)$. This is the event bound in that packet with
$|U|\le q^2$. For $5\le p\le z$, set
$k_p=\lceil 10\log n/\log p\rceil$ and $q=\prod p^{k_p}$.
Chebyshev's estimates give $\log q\ll z$ for $z=cn$, since
$\log q\le 10\pi(z)\log n+\sum_{p\le z}\log p$. Choose $c$ small
enough that the total-variation error is exponentially small.
Each row is primitive at every $p$ by its coefficients.
For a uniform primitive vector over $\mathbb Z/p^{k_p}$, a primitive
linear form $\ell$ has
$$
 \Pr[p^j\mid\ell]=\frac{p^{1-j}}{p+1}\quad(1\le j\le k_p):
$$
complete the row to an invertible matrix; then its first coordinate
is divisible by $p^j$ and its second is a unit. CRT gives independence
between the local factors. A valuation at least $k_p$ might not be
read exactly, but its probability is at most $p^{-k_p}\le n^{-10}$.
A union bound over all rows and primes, plus total variation, shows
that the expected bounded sign product differs by $O(n^{-9})$ from
the corresponding product of means for Haar-uniform primitive
$p$-adic pairs (the inverse limit of the finite uniform laws).

The pairwise determinants of the three rows have absolute values
$173080=2^3\cdot5\cdot4327$ and $17308=2^2\cdot4327$
(the latter twice). The prime $4327$ can be checked by trial division.
For all $p\ge5$ except $5,4327$ the lines on which the selected
forms vanish modulo $p$ are distinct. A primitive pair cannot
lie on two of them. For one form the probability of an odd
valuation is
$$
 \sum_{\substack{j\ge1\\j\ {\rm odd}}}
 \frac{p-1}{p+1}p^{-j}=\frac p{(p+1)^2},
$$
by subtracting successive tails in the preceding formula.
Thus at an ordinary prime the mean of the local sign product
is exactly $1-2|I|p/(p+1)^2$. At $p=5$ the first two lines
coincide and the gap line is distinct; at $4327$ all coincide.
Regardless of higher-power coincidences, at these primes the
probability that **any** selected form has positive valuation is
at most $2/(p+1)$ and $1/(p+1)$ respectively. Both resulting sign
means are strictly positive. At 2 and 3 all the signs are +1.
The ordinary local factors are positive and equal to
$1-2|I|/p+O(p^{-2})$; Mertens' product theorem gives the assertion.
For the original single read $h+2R$, there are no exceptional
lines and its row is primitive at every prime. ∎

By sign-character inversion, the three signs of the pair and its scaled
gap are jointly uniform up to $O((\log n)^{-2})$ across the eight patterns.
This is sensitive to prime-factor **multiplicity**, not merely whether
small primes divide. The cutoff is linear in word length, which is only
an order of logarithm of the output magnitude. Prime factors above the
cutoff can change each sign. There is no representation of those
unseen signs in this finite-modulus argument.

## 3. An unconditional simultaneous almost-prime consequence

The three forms also admit a sieve to exponentially large small-prime
cutoff, where one no longer tracks parity. At primes $p\ge5$
the local density for their product to vanish is
$3/(p+1)$ generically, $2/(5+1)$ at 5, and $1/(4327+1)$ at 4327.
These follow from the same line count; CRT multiplies them
over squarefree moduli. The mixing remainder for $d\mid L_0L_1G$
is $O(3^n d e^{-\delta n})$. Summing to $d\le e^{\delta' n}$ for
small fixed $\delta'$ gives a power-saving error, and the finite
exceptions leave sieve dimension three. The
[Bourgain–Gamburd–Sarnak combinatorial sieve, Theorem 3.1](https://www.ias.edu/sites/default/files/math/Bourgain_Gamburd_Sarnak_Affine.pdf)
then gives $\gg 3^n/n^3$ words for which the entire product
has no factor below $\exp(\delta' n/30)$ (2 and 3 are absent).
All three positive forms are bounded above by a fixed exponential
in $n$ and grow at least exponentially. Thus there are infinitely many
distinct physical pairs of ordinary counts for which both counts and
the gap divided by 10 have bounded total prime-factor multiplicity.
This is not a fixed-gap prime theorem.

## 4. What the linked preprint contributes

[Allahyar–Yianni, *Infinitely many small prime gaps of 180 or less*,
Theorems 1.1–1.2 and §§2, 4–6](https://github.com/sahibzada-allahyar/math/blob/main/preprint/prime-gaps-180/prime_gaps_at_most_180.pdf)
gets below the tuple diameter by subtracting the weighted case where
only the two extreme shifts are prime. The status matters: the claimed
180 is conditional on the stated analytic and numerical hypotheses;
in particular its A4 Type-III correlation is not established by the
Lean implication. The two useful transfers here are to keep a labelled
pair of actual outputs and to isolate correlations at bad primes rather
than sacrifice the entire pair estimate. The integer determinants
make the exceptional local overlaps explicit. Our finite-prime parity
calculation and simultaneous sieve do not use A4 or an ambient
prime-gap theorem.

The missing input for **full** parity is a correlation estimate on the
part of each integer not described by this small modulus. The preprint's
smooth-modulus prime estimates concern a different ambient sum. They
cannot simply be substituted for a von Mangoldt or full-Liouville
estimate along the literal word outputs. That is an explicit
difference of sums, not a verdict that the full-parity problem is
impossible.

Exact finite checks: `python3 check_small_prime_parity.py`.
