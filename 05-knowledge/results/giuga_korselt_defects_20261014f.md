# Signed Giuga defects, unique prime completion, and retained old-prime guards

2026-10-14. **PROVED:** the scoped elementary interfaces below.
**FINITE-EXACT:** the declared integer census and independent compiler controls.
**CITED:** the classical Giuga criterion and primary-pseudoperfect constructions.
The classical Giuga conjecture, twin-prime conjecture, and binary Goldbach conjecture
remain **OPEN**. No historical novelty is claimed for the arithmetic identities.

The executable gain is a complete one-prime completion consumer: for any supplied
squarefree core $N>1$, at most one candidate prime can complete each of two signed
local conditions. A separate exact CRT consumer retains all old-prime Giuga and
Korselt obligations, while explicitly leaving new-prime obligations unpaid.

## 1. Inheritance, names, and the missing coordinates

The closest proved mechanism is the classical criterion already reconstructed in
[the Giuga circuit audit, §§1–4](collatz_mod6_20260921_grand_circuit_typing.md).
Its corrected near miss is the expression $p^2(p-1)\mid N-p$: this is exactly
the classical joint criterion, not a stronger condition or a Collatz theorem.
[The Bernoulli boundary, §7](bernoulli_boundary_20260925.md) retains the
canonical hostile pair: 30 is Giuga but not Carmichael; 561 is Carmichael but
not Giuga. [The Carmichael block compiler](carmichael_successor_blocks_20261014b.md)
already separates successful simultaneous extensions from failed intermediate
products. The least-used sidecar here is the **signed cofactor sum modulo the
unchanged core**, whose inverse determines the only possible new prime.

The board is: signed reciprocal index; source factorization; predecessor period;
prime insertion; neighbouring primes; old versus newly introduced obligations.
The anchor is exact arithmetic completion, the niche is its signed inverse
grammar, and the wildcard is the twin-prime interpretation. The research cards
used are “search the statement before the method” and “retain local incidence
through counting.”

**Guy D19 is a different route.** Guy's book labels D19 as rational distances
from square corners; its final paragraphs also discuss integer-sided triangle
base/height ratios. The latter identification is explicit in MacLeod's primary
paper [introduction, p9](https://ftp.gwdg.de/pub/misc/EMIS/journals/SWJPAM/Vol2_2003/3.pdf).
The [book's table of contents](https://books.google.com/books?id=EbLzBwAAQBAJ)
confirms the section heading. These are geometric rational-distance questions;
D19 is not the name of the Giuga criterion. No arithmetic-to-distance map is
asserted here.

## 2. Two signed local predicates and the full classical criterion

For squarefree $N>1$, retain its distinct prime factors and put

$$
A(N)=\sum_{p\mid N}\frac Np,\qquad
\lambda_N=\mathop{\rm lcm}_{p\mid N}(p-1),\qquad
\Delta_\epsilon(N)=A(N)-\epsilon,\quad \epsilon\in\{-1,1\}.
$$

Write $S_\epsilon(N)$ when $N\mid\Delta_\epsilon(N)$, and retain the integer

$$
\kappa_\epsilon(N)=\frac{A(N)-\epsilon}{N}
=\sum_{p\mid N}\frac1p-\frac\epsilon N.                 \tag{1}
$$

Modulo an old prime $p$, all summands of $A(N)$ but $N/p$ vanish. Thus

$$S_\epsilon(N)\quad\Longleftrightarrow\quad
N/p\equiv\epsilon\pmod p\quad\hbox{for every }p\mid N.   \tag{2}$$

The plus sign includes every prime with index zero; its composite members are
Giuga numbers, with positive integer index. The minus sign is the weak primary
pseudoperfect condition; **index one** defines primary pseudoperfect numbers.
We never replace an arbitrary integral index by one. The parameter convention
in the $\mu$-Sondow literature is $\mu=-\epsilon$; see
[Grau–Oller-Marcén–Sadornil](https://arxiv.org/abs/2111.14211).

For arbitrary $N$, the multiplicity-aware arithmetic derivative is
$N'=\sum_{p^e\Vert N}eN/p$. If $p^2\mid N$, then $p\mid N'$, so
$N'\equiv\pm1\pmod N$ is impossible. On squarefree inputs $N'=A(N)$;
hence the derivative and signed local tests are equivalent, with squarefreeness
forced. The plus characterization is classical; see
[Grau–Oller-Marcén, JIS 15 (2012)](https://cs.uwaterloo.ca/journals/JIS/VOL15/Oller/oller5.html).
This is an arithmetic derivative, not a differential or Jacobian of a map.

Let $P(N)=1+\sum_{a=1}^{N-1}a^{N-1}$. Reduction modulo a prime $p\mid N$
gives the exact two-channel gate

$$P(N)\equiv
\begin{cases}1-N/p,&p-1\mid N-1,\\1,&p-1\nmid N-1\end{cases}\pmod p. \tag{3}$$

Indeed each residue is repeated $N/p$ times, and the finite-field power sum
is $-1$ or 0 according to whether its exponent is divisible by $p-1$.
If $p^2\mid N$, (3) is 1, excluding the congruence. For squarefree $N$, CRT
therefore proves

$$P(N)\equiv0\pmod N\quad\Longleftrightarrow\quad
S_+(N)\ \hbox{and}\ \lambda_N\mid N-1.                 \tag{4}$$

For composites this is precisely “Giuga and Carmichael.” The equivalence is
classical, not the open conjecture; the latter says no composite passes both.
Primary source: [Borwein–Maitland–Skerritt, Theorems 3–4](https://www.carmamaths.org/resources/jon/giuga2013.pdf),
following Giuga and Borwein–Borwein–Borwein–Girgensohn. We rederive (3) so the
implementation consumes the full source and both gates.

Since $N-1\equiv N/p-1\pmod{p-1}$ and $\gcd(p,p-1)=1$, (4) is also exactly
$p^2(p-1)\mid N-p$ for every $p\mid N$. A prime input satisfies it, so a
composite guard is required when discussing a counterexample.

## 3. A complete one-prime compiler for an arbitrary squarefree core

Fix squarefree $N>1$, not necessarily signed. Since
$A(N)\equiv N/p\ne0\pmod p$, $A(N)$ is a unit modulo $N$.
For a prime $q\nmid N$,

$$S_\epsilon(Nq)\quad\Longleftrightarrow\quad
q\equiv\epsilon A(N)^{-1}\pmod N,\qquad q\mid N-\epsilon. \tag{5}$$

The first condition is (2) at all old primes; the second is (2) at the new
prime. They are independent obligations and together are sufficient.

Let $r\in\{1,\ldots,N-1\}$ be the displayed residue. For $\epsilon=1$,
$q\le N-1$, so the only possible prime is $q=r$. For $\epsilon=-1$,
$q\le N+1$; if $r>1$ the only candidate is again $r$, while $r=1$ leaves
the single prime candidate $N+1$ because 1 is not prime. Thus **each sign
has at most one candidate**, followed by exact primality and divisor checks.
This is complete without a search cutoff on $q$. The excluded source $N=1$
would have infinitely many prime completions and is a necessary domain guard.

For example, $N=15$ is neither signed, but $A(15)=8$ and $8^{-1}=2\pmod{15}$:
the plus compiler returns $q=2$, producing 30. The source need not already
pass the target predicate.

## 4. Complete signed transition table and its small boundary

Now assume $S_\alpha(N)$. Equation (5) becomes
$q\equiv\alpha\epsilon\pmod N$, $q\mid N-\epsilon$.

| Source sign $\alpha$ | Target sign $\epsilon$ | All possible new primes |
|---|---|---|
| $+1$ | $+1$ | none |
| $-1$ | $-1$ | $q=N+1$, if prime |
| $-1$ | $+1$ | $q=N-1$, if prime |
| $+1$ | $-1$ | only $(N,q)=(2,3)$ or $(3,2)$ |

For the last row, $q\equiv-1\pmod N$ and $q\mid N+1$. If $N\ge3$, the
only representative below $N+1$ is $N-1$, and $N-1\mid2$ forces $N=3$.
For $N=2$, the next representative is 3 and is valid. In particular every
composite Giuga source has **no one-prime extension on either sheet**.

The exact index law is

$$\kappa_\epsilon(Nq)-\kappa_\alpha(N)
=\frac{N+\alpha q-\epsilon}{Nq}.                         \tag{6}$$

It is zero for both minus-source rows. It is one for the two exceptional
plus-source transitions to 6. The number 2 has **both** signed presentations:
$A(2)=1=0\cdot2+1=1\cdot2-1$. Forgetting this marker caused the initial
hand table to miss $(2,3)$; the independent divisor-path test caught it before
publication. Both presentations and the resulting index laws are retained.

At index one, these are the classical successor/predecessor constructions
$N(N+1)$ and $N(N-1)$; see [Sondow–Tsukerman, Theorem 8(i–ii)](https://arxiv.org/pdf/1401.0322).
The complete table and arbitrary-core compiler are supplied here as exact
consumers, without a claim that the constructions are historically new.

The verified forks are

$$6\longmapsto(30,42),\quad42\longmapsto(1722,1806),\quad
47058\longmapsto(2214408306,2214502422),$$

where each pair lists plus and minus outputs. For these inputs the two
neighbours are twin primes. The implications are conditional on **each**
neighbour passing its own prime test; one test has not been shown to imply
the other. At 1806 both branches stop because
$1805=5\cdot19^2$ and $1807=13\cdot139$. This says nothing about infinitely
many twin primes or representations of an arbitrary even integer as a prime
sum. It records a specified pair $N-1,N+1$ whose sum is $2N$.
The source-bound `primary_fork` returns both authenticated completions and
thereby supplies the literal witnesses

| Primary source $N$ | Binary witness for $2N$ | Ternary witness for $2N+3$ |
|---:|---|---|
| 6 | $5+7=12$ | $5+7+3=15$ |
| 42 | $41+43=84$ | $41+43+3=87$ |
| 47058 | $47057+47059=94116$ | $47057+47059+3=94119$ |

Failure of this fork is not failure of Goldbach: it concerns only these
prescribed two neighbours, with their primary-pseudoperfect source retained.

Removing every prime factor and retesting both signs is a finite complete
inverse to this one-prime grammar. The Giuga number $858=2\cdot3\cdot11\cdot13$
has no signed parent: the grammar does not cover all Giuga numbers. In
particular “every successful object has a successful one-factor predecessor”
is false. Root's separate block compiler addresses simultaneous repairs.

## 5. Block defect law and the limited product closure

For coprime squarefree $N,B>1$, the exact carrier is

$$A(NB)=B A(N)+N A(B).$$

More explicitly, for any signs $\alpha,\beta,\epsilon$,

$$\Delta_\epsilon(NB)=B\Delta_\alpha(N)+N\Delta_\beta(B)
+\alpha B+\beta N-\epsilon.                              \tag{7}$$

For a signed source $S_\alpha(N)$ and a squarefree block of new primes $B$,

$$S_\epsilon(NB)\Longleftrightarrow
B\equiv\alpha\epsilon\pmod N\quad\hbox{and}\quad
q\mid NB/q-\epsilon\text{ for every }q\mid B.             \tag{8}$$

The source index updates by
$\kappa_\epsilon(NB)=\kappa_\alpha(N)+\sum_{q\mid B}1/q
+\alpha/N-\epsilon/(NB)$.
The block API verifies (8) and the full output independently; it does not
assume the block itself is signed or scan all possible blocks.

There is a sharp description when both factors are already signed. Suppose
$2<N<B$, $S_\alpha(N)$ and $S_\beta(B)$. Then their product has sign $\epsilon$
**iff** $B=N+1$, $\beta=-\alpha$, and $\epsilon=\alpha$.
Indeed (8) at both factors forces $N\equiv\epsilon\beta\pmod B$, hence
$N=B-1$ and $\epsilon\beta=-1$; then $B\equiv1\pmod N$ forces
$\epsilon\alpha=1$. The converse is immediate. For smaller factor 2, the
sole extra presentation is $2\cdot3=6$ with both factors on the plus sheet
and output on the minus sheet. Equal-sign product nonclosure in its original
form is classical [Sondow–Tsukerman, Theorem 10](https://arxiv.org/pdf/1401.0322).

## 6. A complete old-prime joint guard, not a completed number

Keep the squarefree core $N$, $A=A(N)$ and $\lambda=\lambda_N$. For a positive
integer multiplier $B$, require the two target conditions only at **old** primes:

$$B(N/p)\equiv1\pmod p,\qquad p-1\mid NB-1\qquad(p\mid N).$$

This is exactly

$$B\equiv A^{-1}\pmod N,\qquad NB\equiv1\pmod\lambda.   \tag{9}$$

It has a solution iff $\gcd(N,\lambda)=1$. Necessity follows from the second
congruence. Under this condition, $B\equiv N^{-1}\pmod\lambda$ and CRT gives
**one class modulo $N\lambda$**. The normal-support condition says no prime
factor of $N$ divides a different prime factor minus one. It is the classical
normal-set restriction, not sufficient for Korselt; $N=15$ is a hostile.

For 561 the exact old-prime guard is

$$B\equiv24241\pmod{44880}.$$

Its least positive member does not complete the full Giuga congruence. Every
new prime still owes both $q\mid NB/q-1$ and $q-1\mid NB-1$, and repeated
prime factors are forbidden. Conversely $\gcd(30,\lambda_{30})=2$, so no
multiple retaining all factors of 30 can ever be Carmichael. This obstruction
concerns the actual retained source, rather than the number of passing tests.

The connection contract is: source = factored squarefree core; target = signed
or joint completion; map = multiply by a guarded new factor/block; preserved
predicate = the explicitly checked local conditions; lost data in a scalar
index = prime labels, new obligations, and period residues; restoration =
factor tuple, sign, source and native completion congruence. The same formal
CRT operation occurs in Collatz inverse guards, but no orbit or convergence
claim follows from these number-theoretic predicates.

## 7. Verification and boundaries

The companion [script](../../04-computation/experiments/giuga_korselt_defects_20261014f.py)
uses exact integers, rational comparisons and trial-division primality.
Private factor caching is reached only after exact public type checks.
No `assert` statement carries a verification gate. `Core`, `Signed`,
`Completion`, and `OldJointGuard` consumers reauthenticate the source and all
stored fields before equality, including bool/float aliases.

The declared universe is all $2\le N\le10000$ for signed tests and complete
one-prime completion versus an independent prime-divisor path; literal power
sums through 400; full-ring power identities through 200; all derivative
products with factors 1 through 100; 18,200 old-guard comparisons on squarefree
cores through 300 and multipliers through 100; 104 specified block controls;
the listed successor/predecessor forks and source/sign/type hostiles. The
all-prime completeness claim comes from (5), not from the finite census.

Reproduce from the repository root:

```text
python -B -X utf8 04-computation/experiments/giuga_korselt_defects_20261014f.py
python -B -O -X utf8 04-computation/experiments/giuga_korselt_defects_20261014f.py
```

The saved [JSON](giuga_korselt_defects_20261014f.json) contains all finite
counts and selected receipts. Normal and optimized runs agree exactly after
LF normalization: **101,343 checks**, including 25 malformed-input controls.
SHA256: `502447be42afde3c5623331710c958760376b794748f36f969bdfbaed66f5873`.
The exact local criterion is not a resolution of the simultaneous composite
exclusion in (4).
