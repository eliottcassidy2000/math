# Global ballot quadratics, bronze uniqueness, and two correction controls

Status: **PROOF CANDIDATE, awaiting independent audit**. The proofs below are elementary exact identities; the finite controls corroborate rather than establish their universal quantifiers. No assertion about LRC, an untouched Laurent coefficient row, or general Newton no-return follows.

## Inheritance and scope

The source is [THM-3010, ballot-column Newton ratios and metallic alternation](../../01-canon/theorems/THM-3010-ballot-column-newton-ratios-and-metallic-alternation.md). Its four displayed central-binomial/Catalan/ballot identities survive. Two nearby claims need correction: integer-order Fuss ratios are rational also for orders at least three; scaled-reciprocal coefficient polynomials can have maximally alternating circuits. The old finite ballot rectangle did not establish a global metallic classification. This note supplies that classification, with the divisor reduction independently suggested by the parent synthesis agent.

The map from a positive sequence to its Newton ratio is the literal map (h_k\mapsto h_k^2/(h_{k-1}h_{k+1})). There is no identification between recurrence characteristic roots and roots of a polynomial whose coefficients happen to satisfy that recurrence.

## 1. The all-integer-shift quadratic

Let (a,b\in\mathbb Z), (h_k=\binom{2k+a}{k+b}), and restrict to

\[
k\ge\max(1,1-b,1-a+b).
\]

All three adjacent binomial coefficients are then defined and positive. Put (\delta=a-2b). Factorial cancellation gives

\[
R_k=\frac{(2k+a)(2k+a-1)(k+b+1)(k+a-b+1)}
{(k+b)(k+a-b)(2k+a+2)(2k+a+1)},
\]

and, with the positive denominator (D_k=(k+b)(k+a-b)(2k+a+2)(2k+a+1)),

\[
1-R_k=\frac{2Q_{a,b}(k)}{D_k},\qquad
Q_{a,b}(k)=k^2+(a+1-\delta^2)k+
\frac{a(a+2)-(2a+1)\delta^2}{4}.
\tag{1}
\]

This identity proves the degree-two statement for every integer shift before cancellation. A reduced numerator can have smaller degree. Whenever it still has degree two, its monic normalization is the same (Q_{a,b}).

**Global metallic classification.** If the monic numerator is (k^2-\nu k-1) for an integer (\nu\ge1), then precisely

\[
(a,b)=(0,-1),(0,1),\qquad \nu=3.
\]

Indeed the constant coefficient condition is (a(a+2)+4=(2a+1)\delta^2). Setting the nonzero odd integer (D=2a+1) gives

\[
4D\delta^2=D^2+2D+13,
\qquad D\mid13,
\qquad \delta^2=(D+2+13/D)/4.
\]

For (D=-1,-13) the last expression is (-3). For (D=1,13) it is (4). Thus (a=0,b=\pm1) gives (k^2-3k-1), whereas (a=6,b=2,4) gives (k^2+3k-1), outside the positive metallic parameter scope. Both quadratics have nonsquare discriminant (13), so neither cancels against a rational linear denominator factor. The classification therefore also holds if one starts from the reduced monic degree-two numerator.

The displayed short circuit words in THM-3010 concern only their displayed indices. For example (a=1,b=-1) has numerator (k^2-7k-6), negative at (7) and positive at (8); it also changes concavity beyond that short bank.

## 2. The integer-order Fuss extension is rational

For any integer (p\ge2), let (H_p(k)=\binom{pk}{k}), (k\ge0). Its adjacent ratio for (k\ge1) is

\[
A_p(k)=\frac{H_p(k)}{H_p(k-1)}=
\frac{\prod_{i=1}^{p}[p(k-1)+i]}
{k\prod_{i=1}^{p-1}[(p-1)(k-1)+i]}.
\tag{2}
\]

For (F_p(k)=H_p(k)/((p-1)k+1)), the adjacent ratio is

\[
B_p(k)=A_p(k)\frac{(p-1)(k-1)+1}{(p-1)k+1},
\qquad R^{F_p}_k=B_p(k)/B_p(k+1).
\tag{3}
\]

These are rational functions at every integer order. For example (B_3(k)=3(3k-1)(3k-2)/(2k(2k+1))). Rationality does not preserve the quadratic numerator phenomenon: at (p=4),

\[
1-R^{F_4}_k=
\frac{432k^4+404k^3-129k^2-101k+24}
{k(2k+1)(3k-1)(3k+1)(4k+1)(4k+3)}.
\]

## 3. Reciprocal symmetry does not prevent maximal alternation

For (N(n)=(n+1)^2(n+3)^2(n+9)^2), write (e_k) for the elementary symmetric root parameters and (h_k=e_k/\binom6k). The five normalized ratios (h_k^2/(h_{k-1}h_{k+1})) are exactly

\[
\left(\frac{65}{57},\frac{4693}{4005},\frac{71289}{61009},
\frac{4693}{4005},\frac{65}{57}\right).
\]

Hence the signs of the four successive log-ratio differences are (+,-,+,-): antipalindromic **and** maximally alternating. The multiset is fixed by (r\mapsto9/r), so this is an exact hostile to a general norm-(+1) implies nonalternation assertion. The disjointness of the particular two-element positive reciprocal pair and a metallic characteristic pair (\{\lambda,-1/\lambda\}) remains valid. It is not a classification of coefficient-polynomial circuit words.

## Reproduction and boundary

The [standalone source](../../04-computation/continuing8_20260906_ballot_repair.py) checks all index-valid rows with (-2\le a\le6), (-3\le b\le4), (k\le11), the symbolic identity (1), all four divisors in its complete global proof, literal binomial/Fuss ratios for (2\le p\le8,1\le k\le9), and the exact reciprocal hostile. There are **2,196 always-active gates**. The [certificate](continuing8_20260906_ballot_repair_certificate.json) retains the exact formulas and hostile ratios; [normal output](continuing8_20260906_ballot_repair.out) and [optimized output](continuing8_20260906_ballot_repair_optimized.out) agree as raw LF bytes.

```text
python continuing8_20260906_ballot_repair.py
python -O continuing8_20260906_ballot_repair.py
```

SHA-256 source `3d537f9a32a4b6f18c152d1986231a2fe4f5def78de0e072f558ad78b434902c`; output `42f65a43f6e9a1ff3c074a2cd28cffbb759603eed16cc34a978a97dd7a043adc`; certificate `151ddbd2240c210ccaf87f5f40f101439af2d58bf4fcebd765a6d6c06cc379d2`.

The successful operation is exact adjacent-ratio cancellation plus a global integer divisibility obstruction. It does not transfer the global classification to arbitrary hypergeometric columns or higher-order Fuss numerator shapes.
