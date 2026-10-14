# Wilson defects, inverse pairs, and a reversible ternary lift tree

**Status:** PROVED elementary classifications and section-lift compiler; FINITE-EXACT declared controls. The Wilson examples $5,13,563$ are certified inputs, not a claim that no others exist. The polynomial Wilson theorem in §6 is CITED in its own coefficient domain. No universal Collatz consequence is claimed.

Companions: [exact checker](../../04-computation/experiments/wilson_defect_trees_20261014c.py) and [saved JSON](wilson_defect_trees_20261014c.json).

## 1. Inheritance and the object being multiplied

The closest repository mechanism is [THM-307, Burnside–QR decomposition](../../01-canon/theorems/THM-307-burnside-qr-decomposition.md), which retains the Wilson quotient as a separate modular correction. We recover its definition, not its unrelated tournament proof. The [Brocard/doubling note, §3](collatz_doubling_tower_brocard_20260927.md) distinguishes $4!+1=5^2$ from $5!+1=11^2$ and $7!+1=71^2$. In particular 11 and 71 in those equations are not further Wilson examples by that argument.

The canonical hostile is the composite 4: it shares the factorial remainder 2 with the prime 3. The corrected near miss is replacing the fixed ordinary residue representatives by adjustable lifts and then claiming to have proved a property of the original factorial. The least-used sidecar is the labelled section and its next defect digit.

Board: fixed source / full ring / unit group / inverse involution / defect digit / free lift address. The anchor is the exact factorial remainder; the niche is a complete reversible finite-depth lift compiler; the wildcard is the genuinely different function-field Wilson criterion.

## 2. Complete factorial-remainder classification

For every integer $n\ge2$,

\[
(n-1)!\bmod n=
\begin{cases}
n-1,&n\text{ prime},\\
2,&n=4,\\
0,&n\text{ composite},\ n\ne4.
\end{cases}
\tag{1}
\]

For a prime, pair the nonzero field elements with their inverses; the only unpaired elements are $1,-1$ (coincident at 2). This proves the first line. For a composite with unequal proper factors $n=ab$, both $a$ and $b$ occur in $(n-1)!$, so $n$ divides it. If the chosen factor pair has $a=b\ge3$, the distinct entries $a,2a<n=a^2$ already contribute a multiple of $a^2$. The only excluded square is 4, where $3!=6$. This proves the classification, including the prime 2 boundary.

It follows immediately that

\[
\boxed{(n-1)!\bmod n=2\quad\Longleftrightarrow\quad n\in\{3,4\}.}
\tag{2}
\]

This is an all-integer theorem, not the result of a finite search. The companion independently multiplies residues for every $2\le n\le2000$.

For an odd prime $p$, define the integer Wilson quotient

\[
W_p=\frac{(p-1)!+1}{p}.
\]

Wilson primality is the stronger predicate $W_p\equiv0\pmod p$. The finite exact census through 2000 finds $5,13,563$. Its next quotient digits are

| $p$ | $W_p\bmod p$ | $(W_p/p)\bmod p$ |
|---:|---:|---:|
| 5 | 0 | 1 |
| 13 | 0 | 4 |
| 563 | 0 | 290 |

Thus each named factorial has exact divisibility $p^2\mid(p-1)!+1$ but $p^3\nmid(p-1)!+1$. The checker compares modular products against literal full factorials; it makes no search-exhaustion assertion above 2000. A receipt retains the prime and precision, and its consumer can require the independently supplied prime. The identical verdict zero at 5 and 13 is not an identifier for either source.

## 3. The exact lift tree and its inverse

Fix an odd prime $p$. A depth-$k$ **labelled section** is

\[
r=(r_1,\ldots,r_{p-1}),\quad
1\le r_a<p^k,\quad r_a\equiv a\pmod p,
\quad\prod_{a=1}^{p-1}r_a\equiv-1\pmod{p^k}.
\tag{3}
\]

Labels are fixed; a permutation is not silently identified with the same section. At depth one the unique section is $(1,\ldots,p-1)$, by Wilson's theorem. Define its next defect

\[
\delta_k(r)=\frac{\prod r_a+1}{p^k}\pmod p.
\]

For digits $t_a\in\{0,\ldots,p-1\}$, expand the lifted product. Terms using two new digits vanish modulo $p^{k+1}$ because $2k\ge k+1$. Since $\prod r_a\equiv-1\pmod p$,

\[
\frac{\prod(r_a+p^kt_a)+1}{p^k}
\equiv\delta_k(r)-\sum_{a=1}^{p-1}a^{-1}t_a\pmod p.
\tag{4}
\]

Therefore the new leaf guard is exactly one nonzero linear equation:

\[
\sum_a a^{-1}t_a=\delta_k(r)\quad\text{in }\mathbb F_p.
\tag{5}
\]

Choose the digits of labels $1,\ldots,p-2$ freely. The last digit is forced:

\[
t_{p-1}=(p-1)\left(\delta_k(r)-\sum_{a=1}^{p-2}a^{-1}t_a\right)\pmod p.
\tag{6}
\]

Every section has exactly $p^{p-2}$ children. Conversely reduce a valid depth-$k+1$ section modulo $p^k$ and read the free upper digits. Equation (5) recovers the pivot uniquely. Repeating this operation reaches the unique depth-one root. This proves a genuine explicit inverse pair

\[
\left(\mathbb F_p^{p-2}\right)^{k-1}
\longleftrightarrow
\{\text{labelled sections at depth }k\},
\qquad \#=p^{(p-2)(k-1)}.
\tag{7}
\]

The map preserves reduction of prefixes, the residue labels, and the guarded product. Forgetting the free address loses the chosen section; retaining it gives exact reconstruction. A fixed pivot is an encoding choice, not a canonical order-independent arithmetic object.

### The genuinely ternary case

For $p=3$ there is one free digit at every depth. At the root $(1,2)$, $\delta=1$, and the three children modulo 9 are exactly

\[
(1,8),\qquad(4,2),\qquad(7,5).
\tag{8}
\]

Each has product $-1\pmod9$. The unchanged section $(1,2)$ does not. Thus a completely branching ternary tree exists even though 3 is not a Wilson prime. The tree derives from one free coordinate over $\mathbb F_3$; it has no asserted bijection with the three named integer Wilson examples.

More generally, the ordinary section $(1,\ldots,p-1)$ is an admissible depth-$k$ node iff $p^k\mid(p-1)!+1$. Its zero-digit continuation is a marked path, distinct from freely correcting representatives. For each fixed prime this unchanged path stops at a finite depth, since the positive integer $(p-1)!+1$ is fixed. The lift tree itself never stops. These quantifiers are the main boundary of the construction.

## 4. A distinguished branch and the Fermat-defect sum

There is a natural branch different from the ordinary integer section. For each $a\in\mathbb F_p^*$, lift the simple root $a$ of $X^{p-1}-1$ uniquely to every $p^k$. The one-digit correction at a current root $r$ is

\[
t=-\frac{r^{p-1}-1}{p^k}
\big((p-1)r^{p-2}\big)^{-1}\pmod p.
\tag{9}
\]

The derivative is a unit, so this is an exact elementary Hensel step. All $p-1$ lifted roots are closed under inversion. Their only self-inverse roots are $1,-1$: if $(x-1)(x+1)$ is divisible by $p^k$, the two factors are coprime at odd $p$. Pairing inverses therefore gives product $-1\pmod{p^k}$. These Teichmüller sections form one coherent branch of (7).

At the first step write $q_p(a)=(a^{p-1}-1)/p\pmod p$. Then (9) becomes $t_a=a q_p(a)$. Equation (5) gives the classical Lerch congruence

\[
\sum_{a=1}^{p-1}q_p(a)\equiv W_p\pmod p.
\tag{10}
\]

For completeness, multiply $a^{p-1}\equiv1+p q_p(a)$ modulo $p^2$. On the other side, $((p-1)!)^{p-1}=(-1+pW_p)^{p-1}\equiv1+pW_p\pmod{p^2}$. Comparing coefficients proves (10). See also Sondow's primary discussion of [Lerch quotients](https://arxiv.org/abs/1110.3113).

At $p=5$ the first lifted section is $(1,7,18,24)$. Individual corrections are nonzero although their weighted sum vanishes. This is an exact example of a zero aggregate defect with retained internal data, not a statement that each factor is unchanged.

## 5. Inverse pairing, the full ring, and Korselt are different tests

In any finite abelian group, inverse pairs cancel in the product of all elements. The product is therefore the product of the subgroup $T$ of self-inverse elements. This is an elementary vector space over $\mathbb F_2$. Its product is the unique nonidentity involution if $|T|=2$, and the identity if $|T|=1$ or $|T|\ge4$: in the latter case each coordinate occurs $2^{\dim T-1}$ times.

Applied to $(\mathbb Z/N\mathbb Z)^*$, this computes the product of **units**. It does not compute $(N-1)!$ when nonunits occur. Exact controls are

| Source | Unit product modulo source | Full factorial remainder | Separate power test |
|---:|---:|---:|---|
| 9 | $-1$ | 0 | not prime |
| 15 | 1 | 0 | $2^{15}\not\equiv2\pmod{15}$ |
| 561 | 1 | 0 | $a^{561}\equiv a\pmod{561}$ for every residue |

Thus the inverse-pair product does not distinguish 15 from the Carmichael input 561. For odd squarefree $N$ with at least two prime factors, CRT supplies at least four involutions, so the unit product is always 1, independently of the extra Korselt divisibilities. A finite field, its higher residue ring, and a squarefree composite ring require different retained markers.

The numerical shift requested in this session also has a precise finite reading: $5-2=3$, $13-2=11$, and $563-2=561$ all satisfy the absolute Fermat power condition. Only three inputs have been established here; no general implication from Wilson primality to a property of $p-2$ is claimed. The parent's successor source has the separate exact factorization

\[
684106401=3\cdot11\cdot17\cdot401\cdot3041,
\quad \lambda=15200,\quad\lambda+1=27\cdot563.
\]

The checker verifies these prime factors and their predecessor lcm. This is a source-specific arithmetic link: the Wilson test uses modulus $563^2$, whereas the successor's divisibilities use its five prime factors. Neither modulus may be discarded.

## 6. What characteristic three changes

Thakur's [*Differential characterization of Wilson primes for $\mathbb F_q[t]$*, Theorems 2.9 and 2.10, pp. 1845–1846](https://www.sas.rochester.edu/mth/sites/dinesh-thakur/wilson2pub.pdf) uses the product of all nonzero polynomial representatives of degree below that of an irreducible polynomial $f$. In odd characteristic, the Wilson property is equivalent to $f''=0$, and there are infinitely many such polynomial primes. These are cited function-field results, not integer Wilson statements. The different coefficient field and its characteristic are essential.

The companion performs a small independent control over $\mathbb F_3$: it enumerates all eight monic irreducible cubics, proves irreducibility by absence of roots, and multiplies all 26 nonzero degree-below-three representatives modulo $f^2$. Exactly two pass the Wilson test, agreeing with the vanishing-second-derivative condition. For example $t^3-t-1$ passes and $t^3+t^2+2t+1$ fails. This supplies a literal positive and hostile example of the domain change; it is not an alternative proof of the cited infinitude theorem.

## 7. Transfer contract and verification

The reusable operation is finite smooth lifting with a marked defect: a rank-one linear constraint on $p-1$ new digits has exactly $p-2$ free digits, with a reconstructible pivot. It can compile all permitted lifts without confusing existence of some lift with acceptance of the supplied unchanged source.

Compare the actual Collatz [terminal-phase compiler](collatz_successor_holonomy_20261014b.md), which solves a scalar principal-unit congruence by ternary digit lifting. Its nonsingular one-variable equation has one permitted next digit. Our two-variable section at $p=3$ has three, because one variable remains free. The common operation is exact guarded refinement; the source objects, predicates and free dimensions differ. A modified Wilson section does not certify a fixed factorial, just as a newly constructed Collatz source does not prove a route for a previously supplied Mersenne integer. No payment or ROOT guard is produced by the Wilson tree.

Run `python -B -X utf8 04-computation/experiments/wilson_defect_trees_20261014c.py` and repeat with `-O`. The declared controls are: every factorial input through 2000; all odd primes through 2000 with direct full-factorial comparisons; Lerch/first-lift controls through 199; every ternary node through depth 8; independent full tuple enumeration at $p=3$ through depth 5 and $p=5$ at depth 2; canonical branches through depth 8 at seven primes; every unit-product input through 500; the eight irreducible cubics; and 20 malformed or wrong-source interfaces. Public receipts validate exact types and labels before equality.

During development a proposed hostile changed a zero Wilson receipt from 5 to 13. That remains a valid receipt for 13: the failure was in the test's assumed source binding. The consumer now accepts an explicit expected prime, and the regression distinguishes internal validity from validity for the caller's unchanged prime.

Normal and optimized runs agree at **13,754 exact checks**. Saved LF-normalized JSON SHA256: `53470a3b3e32e10f1928b937e0a58bf80c69ec94a719a176c8d419c60491b156`.
