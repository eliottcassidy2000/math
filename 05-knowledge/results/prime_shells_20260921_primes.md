# Prime shells: identifying the two harmonic sequences

**Status: PROVED** elementary harmonic/binomial identities and the affine-period template below. **CITED** Bernoulli congruence, historical search bounds, and the complete irregular-index lists for the larger named primes. **FINITE-EXACT** all stated local computations and bounded searches. **OPEN** new exceptional-prime discovery outside the cited searches and any implication for Collatz convergence. No new-prime or literature-priority claim. Date: 2026-09-21.

The two supplied sequences are the same prime-dependent residue in two coordinates. For a prime $p\ge5$, put

\[
 \eta_p=\frac{H_{p-1}}{p^2}\pmod p,
 \qquad
 \beta_p=\frac{\binom{2p-1}{p-1}-1}{p^3}\pmod p.
 \tag{PSP1}
\]

Then $\beta_p=2\eta_p\pmod p$. In particular:

| $p$ | 5 | 7 | 11 | 13 | 17 | 19 | 23 |
|---|---:|---:|---:|---:|---:|---:|---:|
| $\eta_p$ | 3 | 6 | 6 | 7 | 10 | 14 | 18 |
| $\beta_p$ | 1 | 5 | 1 | 1 | 3 | 9 | 13 |

Here one divides in the rational numbers **before** reducing modulo $p$; division by $p^2$ inside the field $\mathbb F_p$ is undefined. The proof guarantees that the resulting rational denominator is prime to $p$.

## 1. Inheritance, board, and the literal data

The closest inherited mechanism is the order-lifting valuation in [row-braid typing, §§2–3](collatz_mod6_20260917_row_braid_typing.md) and the LTE reading of square factors in [wild typing, C3](collatz_mod6_20260917_wild_typing.md). Those notes already use the exceptional prime 1093; this lane does not claim to discover that obstruction. The canonical hostile is $1093^2\mid(4^{182}-1)/3$, defeating a squarefreeness rule that retained the index but discarded the lifting valuation. The corrected near miss is the distinction between base-2 and base-4 pseudoprimality at $85$. The least-used sidecar here is the **index** of a vanishing Bernoulli number: “irregular” means some allowed index, whereas “Wolstenholme” singles out $p-3$.

[THM-307, Burnside-QR decomposition](../../01-canon/theorems/THM-307-burnside-qr-decomposition.md) already keeps Euler and Wilson quotients as different inputs. Its special-prime examples are inherited context, not evidence that their zero sets agree.

The live concept board is: reciprocal reflection; binomial products; Bernoulli indices; prime-power lifting; odd-square shells. The anchor is identifying the user's residues; the niche is a bounded irregular-index search; the wildcard is transporting the extra harmonic digit to an affine residue permutation. The board's useful update is that reflection explains the factor 2 between the supplied sequences, while it does **not** identify their zeros with the half-harmonic zeros controlling the base-2 shell.

The elementary arithmetic is:

\[
 H_4=1+\tfrac12+\tfrac13+\tfrac14=\frac{25}{12},
 \qquad H_6=\frac{49}{20}.
\]

These are sums of the first four or six reciprocals. If “the sum of the first four harmonic numbers” is read literally, then

\[
 \sum_{j=1}^nH_j=(n+1)H_n-n,
 \qquad \sum_{j=1}^4H_j=\frac{77}{12},\quad
 \sum_{j=1}^6H_j=\frac{223}{20}.
\]

Also $\binom92=36\ne126=5^3+1$. The matching binomial coefficients are $\binom94=\binom95=126$. This repair is especially relevant: $\binom94-1=5^3$, so it is precisely the $p=5$ instance of $\beta_p$.

In fact **5 is the unique positive integer** $n$ such that $\binom{2n-1}{n-1}=n^3+1$. Direct calculation covers $1\le n\le5$. For $n\ge5$, the binomial's next-term ratio is $4-2/(n+1)\ge11/3$, whereas the ratio for $n^3+1$ is strictly less than $((n+1)/n)^3\le216/125<11/3$. Induction gives a strict inequality thereafter. This elementary growth argument explains the exact isolated equality, without treating 126 as a perfect power or invoking Catalan's theorem.

An exact sieve gives

\[
 \pi(16843)=1944,\quad \pi(2124679)=157504,
 \quad \#\{p:16843<p<2124679\}=155559.
\]

Both endpoints must be excluded. Including both gives 155561; including exactly one gives 155560. A separately implemented segmented sieve agrees on the strict interval.

## 2. Why the two sequences agree after multiplication by 2

All congruences in this section are in the rationals whose denominators are prime to $p$. Let

\[
 S_j=\sum_{k=1}^{p-1}k^{-j}.
\]

For $1\le j<p-1$, $S_j=0\pmod p$. Indeed inversion permutes the nonzero residues, and there is a nonzero $a$ with $a^j\ne1$: the polynomial $X^j-1$ has fewer than $p-1$ roots. Multiplication by $a$ then gives $S_j=a^{-j}S_j$.

Reflection $k\mapsto p-k$, expanded modulo $p^3$, gives

\[
 2S_1\equiv-pS_2-p^2S_3\pmod {p^3}.
 \tag{PSP2}
\]

For $p\ge5$, $S_2=0\pmod p$ and $S_3=0\pmod p$. Reducing (PSP2) modulo $p^2$ proves **Wolstenholme's harmonic congruence** $S_1=0\pmod {p^2}$, and therefore makes $\eta_p$ well-defined.

Now expand the finite product

\[
 Q_p:=\binom{2p-1}{p-1}
 =\prod_{k=1}^{p-1}(1+p/k).
\]

The first three elementary symmetric sums give, modulo $p^4$,

\[
 Q_p\equiv
 1+pS_1+\frac{p^2}{2}(S_1^2-S_2)
 +\frac{p^3}{6}(S_1^3-3S_1S_2+2S_3)
 \equiv1+pS_1-\frac{p^2}{2}S_2.
\]

The denominators 2 and 6 are units for these primes. Equation (PSP2) reduces to $pS_2\equiv-2S_1\pmod {p^3}$, hence

\[
 \boxed{Q_p\equiv1+2pH_{p-1}\pmod {p^4}},\qquad
 \boxed{\beta_p=2\eta_p\pmod p}.
 \tag{PSP3}
\]

This proves both directions:

\[
 p\mid\operatorname{num}(H_{p-1}/p^2)
 \iff p^3\mid\operatorname{num}(H_{p-1})
 \iff Q_p\equiv1\pmod {p^4}.
 \tag{PSP4}
\]

The $p=3$ boundary is excluded: $H_2=3/2$ is not divisible by $3^2$. Nor does the $p^2$ factor imply that the complete harmonic numerator is a square: $H_{10}=7381/2520=11^2\cdot61/2520$.

## 3. Half-harmonic zeros detect a different exceptional prime

For an odd prime $p$, define the Fermat quotient

\[
 q_p(2)=\frac{2^{p-1}-1}{p}\pmod p.
\]

From $2^p-2=\sum_{k=1}^{p-1}\binom pk$ and

\[
 \frac1p\binom pk=\frac1k\binom{p-1}{k-1}
 \equiv\frac{(-1)^{k-1}}k\pmod p
\]

we obtain

\[
 2q_p(2)\equiv
 \sum_{k=1}^{p-1}\frac{(-1)^{k-1}}k
 =H_{p-1}-H_{(p-1)/2}
 \equiv-H_{(p-1)/2}\pmod p.
 \tag{PSP5}
\]

Thus $p\mid\operatorname{num}(H_{(p-1)/2})$ iff $p$ is base-2 Wieferich. This classical Eisenstein congruence is also recorded in [Gy, *Extended Congruences for Harmonic Numbers*, (1.3), 2019 preprint](https://arxiv.org/pdf/1902.05258); the elementary proof above is sufficient here.

The two harmonic zero tests use different cutoffs and different division depths:

| Test | Harmonic coordinate | Exceptional condition |
|---|---|---|
| Base-2 Wieferich | $H_{(p-1)/2}\pmod p$ | $2^{p-1}=1\pmod {p^2}$ |
| Wolstenholme | $H_{p-1}/p^2\pmod p$ | $Q_p=1\pmod {p^4}$ |

The [shell-lifting lane](prime_shells_20260921_shells.md) supplies the precise bridge from the first test to doubling on prime-power odd-square fibres. Neither harmonic test alone is a Collatz descent criterion.

## 4. Bernoulli indices, the named primes, and the small search

Use $t/(e^t-1)=\sum_{j\ge0}B_jt^j/j!$, so $B_1=-1/2$. An odd prime is **irregular** if some even index $2\le j\le p-3$ has $B_j=0\pmod p$. McIntosh's Theorem 2 gives, for $p\ge11$,

\[
 \eta_p\equiv-\frac{B_{p-3}}3
 \equiv-\frac1{63}\sum_{p/6<k<p/4}k^{-3}\pmod p.
 \tag{PSP6, CITED}
\]

This follows by combining (PSP3) with the theorem's binomial congruence. The Bernoulli equality also holds at 5 and 7 by direct evaluation; the short-sum formula is used only at $p\ge11$. Source: [R. J. McIntosh, *On the converse of Wolstenholme's Theorem*, Acta Arith. 71 (1995), pp. 386–387](https://matwbn.icm.edu.pl/ksiazki/aa/aa71/aa7144.pdf), original PDF retrieved and read. Consequently Wolstenholme primes are irregular at the particular index $p-3$; an arbitrary irregular index is insufficient.

| $p$ | $\eta_p$ | $\beta_p$ | $q_p(2)$ | Complete irregular-index list |
|---:|---:|---:|---:|---|
| 1093 | 659 | 225 | 0 | none |
| 2851 | 3 | 6 | 671 | none |
| 3511 | 3321 | 3131 | 0 | 1416, 1724 |
| 12101 | 2 | 4 | 3686 | 7718 |
| 16843 | 0 | 0 | 4352 | 16840 |
| 2124679 | 0 | 0 | 1913833 | 701898, 2124676 |

The complete lists are the exact rows of [Buhler's author dataset](https://people.reed.edu/~jpb/bernoulli/) ([compressed table](https://people.reed.edu/~jpb/bernoulli/To163.gz)). Our script independently recomputes **every** allowed index for 1093 and 2851, and verifies each positive listed index for the other four primes. Completeness of the larger lists remains CITED, not a claim about what this script exhaustively searched.

The positive witnesses use the finite-field moment identity in [Hart–Harvey–Ong, *Irregular primes to two billion*, 2016 preprint, §2, Eq. (1)](https://arxiv.org/pdf/1605.02398). For even $2\le r\le p-3$, choose $1<c<p$ with $c^r\ne1\pmod p$, and set

\[
 f_c(x)=\left\lfloor\frac{c\,[x/c]_p}{p}\right\rfloor-\frac{c-1}{2}.
 \qquad
 \frac{c^r-1}{r}B_r=\sum_{x=1}^{p-1}x^{r-1}f_c(x)\pmod p.
\]

This is a cited formula used as an independent exact computation, not a new derivation. The paper reports a complete search below $2^{31}$; the present search is a much smaller recovery of already known examples.

Our exhaustive search of primes $p\le300$, with two independent Bernoulli recurrences, gives the following irregular pairs, grouped by prime:

\[
\begin{array}{c|l@{\qquad}c|l}
37&32&59&44\\
67&58&101&68\\
103&24&131&22\\
149&130&157&62,110\\
233&84&257&164\\
263&100&271&84\\
283&20&293&156
\end{array}
\]

All other odd primes at most 300 are regular. In particular $37\mid\operatorname{num}(B_{32})$, but $\eta_{37}=24\ne0$: irregular does not imply Wolstenholme. Conversely 1093 is a Wieferich prime and is regular, so Wieferich does not imply irregular.

The bounded small-residue search over **all primes $5\le p\le20000$** gives:

| $\eta_p$ | Primes in this bounded universe |
|---:|---|
| 0 | 16843 |
| 1 | 103 |
| 2 | 12101 |
| 3 | 5, 2851 |

This supplies a concrete reading of the appearances of 12101 and 2851. Small **nonzero** representatives do not satisfy the exceptional congruence, and smallness depends on the chosen residue coordinate; multiplication by 2 preserves zero but need not preserve a chosen notion of nearness to zero. No distribution law follows from this small table.

For clarity, **Wilson** means $((p-1)!+1)/p=0\pmod p$. **Fibonacci-Wieferich**, also called **Wall–Sun–Sun**, means $F_{p-(5/p)}/p=0\pmod p$, for primes $p\ne2,5$. The recurrence $F_0=0,F_1=1,F_{n+1}=F_n+F_{n-1}$ fixes the indexing. Our complete tests to 10000 give Wieferich ({1093,3511}), Wilson ({5,13,563}), and no Fibonacci-Wieferich prime. None of the six named primes in the table is Wilson or Fibonacci-Wieferich; their nonzero quotients are saved in JSON. These statements do not classify all primes.

The definitions and classical examples are in [Crandall–Dilcher–Pomerance, *A search for Wieferich and Wilson primes*, Math. Comp. 66 (1997), pp. 433–449](https://math.dartmouth.edu/~carlp/PDF/paper111.pdf). For historical scope, [McIntosh–Roettger, *A search for Fibonacci-Wieferich and Wolstenholme primes*, Math. Comp. 76 (2007), pp. 2087–2094](https://doi.org/10.1090/S0025-5718-07-01955-2) reports no Fibonacci-Wieferich primes below $2\cdot10^{14}$ and no additional Wolstenholme primes below $10^9$. Its original article was read through [this full-article mirror](https://manuals.plus/m/a8e6070239b460b7f0495935e115355357ab8a074e0348f4b6e6fe04a38e61db); these are dated publication results, not an assertion about the latest search frontier.

## 5. A genuine shared lifting template, with the lost predicate explicit

Let $p$ be odd, $Q>1$ even with $Q\equiv1\pmod p$, and write

\[
 r=v_p(Q-1),\quad c=(Q-1)/p,\quad R(n)=Qn+c.
\]

On the $p^s$ odd residue classes modulo $2p^s$, $R$ is a permutation and **every** cycle has length

\[
 \boxed{p^{\max(0,s-r+1)}};
 \qquad\#\text{cycles}=p^{\min(s,r-1)}.
 \tag{PSP7}
\]

Proof: $c$ is odd and, on writing $n=2z+1$, the induced slope on $z\pmod {p^s}$ is the unit $Q$. Moreover

\[
 pR^t(n)+1=Q^t(pn+1),\qquad
 v_p(R^t(n)-n)=v_p(Q^t-1)-1=r+v_p(t)-1.
\]

The last step is the odd-prime LTE identity; $pn+1$ is a unit. Both iterates are odd, so the factor 2 in the modulus is automatic. The displayed valuation gives the least positive return time and the cycle count. In particular every point is fixed for $s\le r-1$, and the period first becomes $p$ at $s=r$. This is the inherited affine-braid valuation mechanism, applied to a different multiplier.

Take now $Q=Q_p=\binom{2p-1}{p-1}$, for a prime $p\ge5$. Legendre's factorial valuation gives

\[
 v_2(Q_p)=s_2(p)-1\ge1,
\]

where $s_2(p)$ is the binary digit sum: $v_2\binom{2p}{p}=s_2(p)$, and $\binom{2p}{p}=2Q_p$. Thus $Q_p$ satisfies the hypotheses. Ordinary Wolstenholme divisibility gives $r\ge3$, while the exceptional prime condition is $r\ge4$. The extra harmonic vanishing therefore gives an **extra fixed residue level** in this deliberately constructed affine family. Higher $r$, if present, gives more fixed levels; no universal bound on $r$ is asserted.

Exact controls give $r=3$ at $p=5,7,11,13,17$, and $r=4$ at 16843. The program completely decomposes the odd residue permutations for $p=5,7,11$ and $1\le s\le4$, and independently checks closed-iterate returns and minimality through $s=5$ for all six selected controls.

The source/target map is $n\mapsto pn+1$, which intertwines $R$ with multiplication by $Q$ on the class $1\pmod p$. It preserves the exact local return valuation, with a one-level shift from division by $p$. It discards order on the positive integers and all forward Collatz trajectory data. Its needed sidecar is $r$.

In particular the binomial multiplier generally has a nontrivial odd part. At $p=5$, $Q=126=2\cdot63$, $c=25$, and $R(1)=151$. The odd core of $5\cdot1+1$ is 3, but that of $5\cdot151+1=756$ is 189. Thus this affine family does **not** preserve the odd core that makes a pure-power-of-two inverse fibre a Collatz object. A common LTE proof is a real arithmetic connection; it is not an isomorphism of those dynamical systems.

## 6. Reproduction and stopping boundary

Run from the worktree root:

```text
python 04-computation/experiments/prime_shells_20260921_primes.py
python -O 04-computation/experiments/prime_shells_20260921_primes.py --output PATH_OUTSIDE_REPOSITORY.json
```

The [script](../../04-computation/experiments/prime_shells_20260921_primes.py) uses standard-library exact integers and fractions, with explicit exceptions instead of `assert`; the [JSON](../../04-computation/experiments/prime_shells_20260921_primes.json) records the source hash, universes, tables, and 19605 active checks. Positive controls include both named Wolstenholme primes, both named Wieferich primes, all three small Wilson primes, and the irregular-index witnesses. Hostiles include $p=37$, regular Wieferich 1093, the nonsquare harmonic numerator at $p=11$, the literal binomial mismatch, the two distinct harmonic summations, and the odd-core loss at $Q=126$.

Independent paths are: full versus segmented sieve; rational harmonic sums versus modular sums and exact binomial coefficients for every prime $5\le p\le300$; full reciprocal sum modulo $p^3$ versus the cited short cubic sum for each named prime; rational versus modular Bernoulli recurrences for the full small census; finite-field moment witnesses for the larger listed indices; and direct residue cycles versus the valuation formula. Completeness of the larger irregular-index lists is explicitly inherited from the author data.

The finite work identifies the data and yields an exact harmonic-to-affine lifting bridge. It discovers no previously unknown exceptional prime and does not estimate the frequency of future exceptions. The useful next search would require a stated new range or a structural restriction on the Bernoulli index. Treating the two residue sequences as independent evidence, or treating an affine modular clock as a proof of Collatz termination, would discard exactly the coordinates exposed here.
