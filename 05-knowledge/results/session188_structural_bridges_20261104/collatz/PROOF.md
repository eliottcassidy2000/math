# A four-letter reset compiler and the ternary address of a periodic Collatz template

**Status: proved elementary construction and sharp bounds, with finite exact tests. Collatz remains open.**

This note extends the repository's native inverse-phase construction from initial-one runs to arbitrary ordered words. Its main additions are the complete reset-phase spectrum with the sharp four-letter threshold, and an exact bridge from a template's rational cycle center to its compatible ternary reset addresses. No claim of historical novelty is made.

## 1. The question and the quantifiers

Write

\[
U(n)=\frac{3n+1}{2^{v_2(3n+1)}}
\]

for positive odd integers. Given a nonempty valuation word \(w=(a_0,\ldots,a_{L-1})\), can one construct positive sources that follow **exactly those letters**, then reach 1 in one further odd step? Can the required final reset be compiled without constructing the potentially astronomical source? Can every permissible reset phase be encoded by a uniformly bounded alphabet?

The answers are yes, yes, and yes with the sharp alphabet \(\{1,2,3,4\}\).

The source is constructed *after* the word is prescribed. This is not a procedure that takes an arbitrary supplied source and proves that it reaches 1. In particular, density in the 2-adic topology below is not natural density or universal convergence.

## 2. Exact completion and its full exponent phase

Let

\[
A=\sum_{i=0}^{L-1}a_i,\qquad
B=\sum_{i=0}^{L-1}3^{L-1-i}2^{a_0+\cdots+a_{i-1}}.
\]

The formal forward composition is

\[
F_w(n)=\frac{3^Ln+B}{2^A}.
\]

### Theorem 1. One-reset completion

There is a unique even residue \(\beta(w)\) in \([0,2\cdot3^L)\) satisfying

\[
2^{\beta(w)}\equiv 1+3B2^{-A}\pmod{3^{L+1}}. \tag{1}
\]

It always satisfies \(\beta(w)\equiv2\) or \(4\pmod6\). The positive integer exponents \(b\) for which an integral source follows \(w\) and then a valuation-\(b\) step to 1 are exactly

\[
b\equiv\beta(w)\pmod{P_L},\qquad P_L=2\cdot3^L. \tag{2}
\]

For each such positive \(b\), the source is

\[
n(w,b)=\frac{2^{A+b}-2^A-3B}{3^{L+1}}. \tag{3}
\]

It is positive and odd, and every prescribed valuation is exact. The route reaches 1 for the **first time** at its final step if and only if \(b\ge4\). Thus its least first-hit reset exponent is

\[
b_*(w)=
\begin{cases}
\beta(w),&\beta(w)\ge4,\\
P_L+2,&\beta(w)=2.
\end{cases} \tag{4}
\]

**Proof.** The endpoint immediately before the reset is \(m_b=(2^b-1)/3\). Solving \(F_w(n)=m_b\) gives (3), whose integrality is precisely (1) with \(b\) in place of \(\beta\).

For completeness, \(v_3(4^{3^j}-1)=j+1\): cubing \(1+3^ku\), with \(3\nmid u\), increases the valuation by exactly one. Hence 4 generates the group of units equal to 1 modulo 3, modulo \(3^{L+1}\), of order \(3^L\). The right side of (1) lies in that group. There is one even exponent modulo \(P_L\). Moreover

\[
B2^{-A}\equiv2^{-a_{L-1}}\not\equiv0\pmod3,
\]

so the right side is not 1 modulo 9. Thus \(3\nmid\beta\), proving the asserted phase restriction. In particular every positive allowed \(b\) is even and at least 2.

Integrality of (3) makes every intermediate inverse step integral. To see this without assuming actual forward behavior, undo the letters from the terminal positive odd integer \(m_b\), using \(E_a(x)=(2^ax-1)/3\). The final numerator congruence modulo \(3^L\), reduced modulo 3, first implies legality of the last forward letter; divide that congruence by 3 and continue. Each legal inverse step from a positive odd integer has a positive odd output: its numerator is positive and odd. Its forward numerator is exactly \(2^a x\), with odd \(x\), so its valuation is exactly \(a\).

If \(b=2\), then \(m_b=1\), so the final step is a repeated root step. If \(b\ge4\), then \(m_b>1\). No earlier state can equal 1 because its subsequent actual trajectory would stay at 1, contradicting \(m_b>1\). This proves the first-hit assertion and (4). QED.

### A compressed algorithm

Compute \(B\) modulo \(3^{L+1}\), using the recurrence \(B\leftarrow3B+2^A\), \(A\leftarrow A+a\). Let the right side of (1) be \(g\). Write \(\beta=2j\). Starting at \(j=0\pmod1\), lift

\[
4^j\equiv g\pmod{3^{k+1}}
\]

to the next level by testing \(j,j+3^k,j+2\cdot3^k\). Exactly one succeeds. This takes \(L\) ternary lifting stages on moduli with at most \(L+1\) ternary digits. It does **not** scan \(3^L\) candidate exponents.

These are stage and modular-precision bounds, not a runtime bound depending only on \(L\). If the valuation letters are binary encoded, the compiler must also read their total input bit length \(N=\sum_i\lceil\log_2(a_i+1)\rceil\). Computing the cumulative cost and modular powers depends on those input bits; a power exponent can be reduced modulo the known order \(2\cdot3^L\). Ordinary integer arithmetic gives a bit complexity polynomial in \(L+N\) for the **phase** and symbolic certificate. The word and cumulative-cost sidecars are not claimed to fit in only \(O(L)\) bits. None of these statements bounds expanded-witness construction by that compact input size.

The default `completion(expand=False)` path retains the word, (3), and the carry formula symbolically, computing only the carry residue modulo \(3^{L+1}\). It does not construct the full carry \(B\). The `expand=True` path deliberately constructs the carry and source. Materializing \(2^{A+b}\) still requires \(A+b+1\) bits; this cost is not hidden by the compressed certificate.

## 3. A sharp finite alphabet for all reset phases

### Theorem 2. Complete spectrum and sharp threshold

For every \(L\ge1\),

\[
\{\beta(w):w\in\{1,2,3,4\}^L\}
=
\{\beta:0\le\beta<2\cdot3^L,\ \beta\equiv2\text{ or }4\pmod6\}. \tag{5}
\]

The spectrum has exactly \(2\cdot3^{L-1}\) elements. There is an explicit compiler from an allowed phase to a word, using only \(O(L)\)-digit ternary residues and \(L\) inverse steps.

The bound 4 is sharp among initial alphabets \(\{1,\ldots,M\}\). At length 2, the alphabet \(\{1,2,3\}\) misses \(\beta=14\). At every length \(L\ge2\), it misses every phase \(\beta\equiv14\pmod{18}\).

**Proof of surjectivity.** Fix an allowed \(\beta\), and consider the residue of \(m=(2^\beta-1)/3\) modulo \(3^L\). It is a unit. From any odd positive unit \(x\), the two candidate inverse exponents are

- \(a=1,3\) if \(x\equiv2\pmod3\);
- \(a=2,4\) if \(x\equiv1\pmod3\).

Both predecessors are integral, and

\[
E_{a+2}(x)=4E_a(x)+1. \tag{6}
\]

Consequently they cannot both be divisible by 3. Select a unit predecessor until the last of the \(L\) inverse steps; at the last step no unit requirement is necessary. This gives a length-\(L\) word in \(\{1,2,3,4\}\), read in reverse order, followed by the reset exponent \(\beta\). Theorem 1 then identifies its phase as \(\beta\).

The actual huge integer \(m\) is unnecessary. Compute

\[
x=\frac{(2^\beta\bmod 3^{L+1})-1}{3}\pmod{3^L}.
\]

At a stage with \(r\) ternary digits, compute the two inverse residues modulo \(3^{r-1}\), discard a 3-divisible option if another step remains, and consume one digit. This supplies exactly the same legal residue choices. Positivity belongs to the actual integer lift supplied by Theorem 1, not to the arbitrarily chosen residue representatives.

**Proof of sharpness.** If \(\beta\equiv14\pmod{18}\), then \(2^\beta\equiv22\pmod{27}\), so

\[
m=(2^\beta-1)/3\equiv7\pmod9.
\]

Because \(m\equiv1\pmod3\), the only legal exponent in \(\{1,2,3\}\) is 2. But

\[
E_2(m)=(4m-1)/3\equiv0\pmod3,
\]

and a multiple of 3 has no integral odd predecessor under any inverse exponent. Thus no prefix of length at least 2 over that alphabet is possible. Exponent 4 supplies the missing alternative through (6). QED.

This is a terminating finite-branch controller that consumes ternary precision. It does not assert that a fixed finite state stores an unbounded residue. The control alphabet is constant; the residue tape has length \(L\), as in the carry-strip interface of S178.

## 4. The reset cost is genuinely exponential

### Theorem 3. Sharp extremal phases

For every \(L\ge1\),

\[
\begin{aligned}
\beta(1^L)&=3^L+1,\\
\beta(1^{L-1},3)&=2\cdot3^L-2,\\
\beta(2^L)&=2,\qquad b_*(2^L)=2\cdot3^L+2.
\end{aligned} \tag{7}
\]

Hence the largest least representative among all allowed reset phases is attained by a word over \(\{1,3\}\). The uniform first-hit bound

\[
b_*(w)\le2\cdot3^L+2 \tag{8}
\]

is sharp, including when every prefix letter is 2. In particular no bound polynomial in \(L\) can replace (8) for this **one-reset completion** problem.

**Proof.** For \(1^L\), \(A=L\), \(B=3^L-2^L\), so the target of (1) is \(-2\) modulo \(3^{L+1}\). Since \(2^{3^L}\equiv-1\pmod{3^{L+1}}\), its phase is \(3^L+1\).

For \((1^{L-1},3)\), the carry is unchanged but \(A=L+2\), and the target becomes \(1/4\). Its phase is \(-2\) modulo \(P_L\), hence \(P_L-2\). For \(2^L\), \(A=2L\), \(B=4^L-3^L\), and the target is 4, giving phase 2. The first-hit condition raises that exponent to \(P_L+2\). QED.

For example, the least first-hit completions of \((1)\), \((1,1)\), \((1,3)\), and \((2,2)\) have sources 3, 151, 38835, and 621377, respectively. Their reset exponents are 4, 10, 16, and 20. These are exact constructions, not a trajectory search used as production input.

The last pair illustrates the cost precisely: a very simple local word can have a huge direct-to-root reset. Allowing several additional odd steps is a different optimization problem; this theorem supplies no lower bound for that more flexible completion.

## 5. The ternary address behind a periodic template's rational center

S179 assigns a word the rational cycle center

\[
D=2^A-3^L,\qquad\rho=B/D.
\]

It uses \(v_2(Dn-B)\) to predict the exact length of a repeated-template shadow. The root completion compiler supplies a complementary **ternary** address derived from the same center.

### Theorem 4. Center-to-reset address

For every \(r\ge1\), the phase of \(w^r\) is the unique even residue modulo \(2\cdot3^{rL}\) satisfying

\[
2^{\beta(w^r)}\equiv1+3\rho
=\frac{D+3B}{D}\pmod{3^{rL+1}}. \tag{9}
\]

The denominator \(D\) is a unit modulo 3. These residues are compatible:

\[
\beta(w^{r+1})\equiv\beta(w^r)\pmod{2\cdot3^{rL}}. \tag{10}
\]

Equivalently, \(\beta(w^r)/2\) is the first \(rL\) ternary digits of the uniquely defined 3-adic exponent \(j\) with

\[
4^j=1+3\rho.
\]

The compatible finite congruences (9) are sufficient to define this address; no numerical 3-adic logarithm approximation is required.

**Proof.** The repeated carry is

\[
B_r=B\frac{2^{rA}-3^{rL}}{D}.
\]

Therefore

\[
1+3B_r2^{-rA}
=1+3\rho-\frac{3^{rL+1}B}{D2^{rA}}. \tag{11}
\]

The error has exactly 3-adic valuation \(rL+1\), since \(B,D,2\) are units modulo 3. Substitution into (1) proves (9). Uniqueness of lifting proves (10) and defines the infinite compatible address. QED.

### Theorem 5. Rational reset addresses occur only at the two fixed templates

The 3-adic address \(j\) in Theorem 4 is rational over \(\mathbb Q\) if and only if the word is entirely 1s or entirely 2s. The corresponding addresses are

\[
j(1^L)=\tfrac12,\qquad j(2^L)=1. \tag{12}
\]

Every other template has an irrational 3-adic address. Its ternary digit sequence is therefore not eventually periodic.

**Proof.** Put \(C=1+3\rho\). The carry identity gives

\[
C=\frac{2^{a_0}B_1}{D},\qquad v_2(C)=a_0\ge1,
\]

because \(B_1,D\) are odd. Suppose \(j=p/q\in\mathbb Q\cap\mathbb Z_3\), in lowest terms, with \(q>0\). Raising the compatible finite congruences to the power \(q\), or using the continuous exponent law in \(1+3\mathbb Z_3\), yields the equality of rational numbers

\[
C^q=4^p.
\]

Taking ordinary prime valuations, every prime other than 2 has valuation zero in \(C\). Thus \(C=\pm2^k\), where \(k=v_2(C)=a_0>0\). The condition \(C\equiv1\pmod3\) fixes the sign and gives \(C=(-2)^k\). Therefore

\[
\rho=\frac{(-2)^k-1}{3}
\]

is an odd integer. The cyclic carry identities show that this center realizes the periodic word under the signed odd Collatz map. But its very next value is \((-1)^k\). A periodic point whose forward orbit reaches a fixed point must be that fixed point itself. Hence \(\rho=-1,k=1\), or \(\rho=1,k=2\). Their actual valuation streams are respectively all 1s and all 2s, forcing exactly the claimed templates.

Conversely \(\rho=-1\) gives \(C=-2\), whose square is 4 and whose residue is 1 modulo 3, so its unique address is \(j=1/2\). The root center gives \(C=4\), hence \(j=1\).

Finally, an eventually periodic ternary expansion sums to a rational number by a geometric-series formula; conversely the finite-remainder long-division algorithm for a rational 3-adic integer is eventually periodic. This proves the digit statement. QED.

For a concrete signed control, \(w=(1,2)\) has center \(\rho=-5\) and the actual signed cycle \(-5\to-7\to-5\), while its address satisfies \(4^j=-14\) and is irrational. Thus this irrationality cannot itself exclude a cycle.

This is a nonperiodicity theorem for the **reset-address digits**, not for the actual Collatz orbit: the template itself is periodic by definition. It does not imply equidistribution of digits, a lower bound on each individual digit, a transcendence result, or an obstruction to more general finite automata reading a non-unary input.

### Corollary 6. Exponential reset cost recurs for every fixed template

Fix any template \(w\) of length \(L\ge1\), and let \(b_r=b_*(w^r)\) be the **least first-hit** one-reset exponent, as in (4). Then

\[
\limsup_{r\to\infty}\frac{\log_3 b_r}{rL}=1.
\]

For the all-1 and all-2 templates the full limit exists and equals 1. For every other template, infinitely many repetition counts satisfy

\[
b_r\ge 2\cdot3^{(r-1)L}.
\]

**Proof.** Uniformly, (8) gives \(b_r\le2\cdot3^{rL}+2\), so the limsup is at most 1. Write \(\beta_r=\beta(w^r)\) and \(P_r=2\cdot3^{rL}\). Compatibility and the least-nonnegative-representative convention give

\[
\beta_r=\beta_{r-1}+k_rP_{r-1},\qquad k_r\in\{0,\ldots,3^L-1\}.
\]

For a template other than all 1s or all 2s, Theorem 5 says the limiting address is irrational. Hence \(\beta_r\) cannot eventually be constant: if it were constantly \(B\), all finite truncations would give the rational address \(j=B/2\). There are therefore infinitely many \(r\) with \(k_r\ge1\). At each of them,

\[
b_r\ge\beta_r\ge P_{r-1}=2\cdot3^{(r-1)L}.
\]

Dividing logarithms by \(rL\) along this infinite subsequence proves the matching lower limsup.

If \(w=1^L\), formula (7) gives \(b_r=3^{rL}+1\), already at least 4, so no first-hit adjustment is needed. If \(w=2^L\), the least positive exponent 2 gives the padded route from 1 to itself; the least **first-hit** exponent is instead \(b_r=P_r+2\). Both explicit formulas have the claimed full limit. QED.

This concerns the cost of appending exactly one reset to the prescribed repeated prefix. It does not bound the cost of a more flexible suffix, and it gives neither a full limit nor an effective bound on the gaps between expensive reset phases for the other templates.

### Exact interface with the S179 mismatch parser

Let \(n=n(w^r,b)\) be a first-hit root completion, and let \(a_0\) be the first template letter. A direct cancellation gives

\[
Dn-B=\frac{2^{rA}\bigl(D2^b-D-3B\bigr)}{3^{rL+1}}. \tag{13}
\]

The cyclic-carry identity gives \(D+3B=2^{a_0}B_1\), where \(B_1\) is odd. Thus if \(b\ne a_0\),

\[
v_2(Dn-B)=rA+\min(b,a_0). \tag{14}
\]

The S179 parser consequently predicts **exactly \(rL\) matching letters**, followed by the constructed reset: if \(b<a_0\) the exit valuation is smaller; if \(b>a_0\) it is larger. When \(b=a_0\), cancellation supplies extra precision and the template may match further; no exact-boundary claim is made in that case.

This gives actual first-hit-rooted witnesses for arbitrarily long shadows of **every** periodic valuation template, expanding or contracting. The binary shadow center and the ternary reset address are linked by the rational affine map \(\rho\mapsto1+3\rho\). That link does not produce a global rank or make template changes descend.

## 6. What is new here, what was inherited, and what remains open

Targeted sources read:

- [S178: exact local detector and carry-strip interface](https://github.com/eliottcassidy2000/math/blob/06316c266ad2fa82cb2632f1b7093ef4f4d4bfa8/05-knowledge/results/session178_prime_flow_honeycomb_20261104/detector/detector.md). The packet was also read from its verified local publication snapshot. It supplies the constant-alphabet versus unbounded-tape distinction; it is not a convergence argument.
- [S179: periodic-template parser and rank](https://github.com/eliottcassidy2000/math/blob/4043ce7be6880f7ef0234aace2fe530e6e4003bb/05-knowledge/results/session179_niche_methods_20261104/collatz/RESULT.md). Its exact mismatch criterion, rational centers, and repeated-word redundancy are inherited.
- [Unbounded run blocks, October 5](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/collatz_run_block_family_20261005.md). The ternary digit lifting for \((1^r,a)\), two unit choices among inverse phases, and honest output-bit cost are inherited. Theorem 1 generalizes this mechanism to arbitrary ordered prefixes.
- [Guarded affine lifts, October 4](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/collatz_affine_guarded_lifts_20261004.md). It already constructs legal ancestors with finite affine prescriptions using power-of-two phases; changing the constructed source is not supplied-source coverage.
- [Branch-toll rank, October 4](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/collatz_branch_toll_rank_20261004.md). Rational-anchor controllers consume precision; that general principle is inherited.
- [Weak-reset family, October 5](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/collatz_weak_reset_family_20261005.md). Positive rational centers already characterize word descents; those thresholds are not newly claimed here.

The concrete additions are (i) the **full** arbitrary-prefix one-reset exponent class with first-hit normalization; (ii) its complete inverse spectrum and sharp alphabet-four threshold; (iii) sharp extremal reset exponents within tiny alphabets; and (iv) the compatible ternary reset address determined by the S179 rational template center, including the exact binary mismatch interface; (v) the classification of rational reset addresses, which proves nonperiodicity of the reset digits for every other template; and (vi) recurring exponential one-reset cost for every fixed template.

Every finite valuation cylinder contains explicitly constructed first-hit-rooted sources. The finitely many words and their certificates are not an all-source cover: an input in the same cylinder can take an entirely different reset. No natural-density improvement, global selector, inter-template rank decrease, exclusion of nontrivial cycles, or exclusion of divergent trajectories is claimed.

A genuinely stronger next target would be to complete a prescribed source, rather than a prescribed word, or to bound the cost of a multi-reset completion while preserving independently supplied source data. The present compiler solves neither problem.

## 7. Reproduction

Run, without third-party packages:

```
python reset_phase.py > RESULTS.json
python -O reset_phase.py > RESULTS.optimized.json
cmp RESULTS.json RESULTS.optimized.json
```

The tests exhaust all 5,460 words over \(\{1,2,3,4\}\) of lengths 1 through 6, independently replay two first-hit root routes for each word, and compare the phase with a brute finite logarithm table. The converse controller is checked for every permitted phase at lengths 1 through 10. Further controls cover the alphabet-three obstruction, closed-form extremes through length 60, periodic centers, exact mismatch precision, and rejected malformed inputs. All production checks raise exceptions explicitly and remain active under optimized Python.
