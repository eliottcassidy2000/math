# Optimizing the inherited strip bound: a logarithmic spectrum improvement

Status: PROVED deduction from THM-4887, subject to independent parent audit. No claim of improved power exponent, positive density, linear energy, or Chan–Pak completeness. No new theorem ID reserved.

## Result

Let g(0)=1, g(2v+1)=g(v), and g(2v)=g(v)+v+(v mod 2). Let E_n count ordered equal-g pairs in [2^(n-1),2^n). For every integer n>=12,

    E_n <= 7.3 * 2^(4n/3) * n^(2/3).

Consequently, if A(X) counts distinct tournament Hamiltonian-path counts <=X, then for every real X>=4096,

    A(X) >= X^(2/3) / (60 (ln X)^(2/3)).

One can take the maximum of this bound and THM-4887's earlier X^(2/3)/(13 ln X+102). The new inequality improves the asymptotic logarithmic loss from power 1 to power 2/3. It is an optimization of that theorem's proved general bound, not a new arithmetic collision mechanism.

## Proof

Put N=2^n. Choose m=floor((n-log_2 n)/3)+1 and C=2^m. Write

    C=c (N/n)^(1/3),  1<c<=2,
    P=N/(2C), B=P/C.

The endpoint c=1 can harmlessly be included in all bounds. For n>=12, B>=8: indeed log_2 B >= n/3+(2/3)log_2 n-3, which is increasing in n and exceeds 3 already at n=12. THM-4887 (D1), its Delta=0 version, gives

    E_n <= PC + 8P^2 + PC^2(33.74+12 ln B) + 3PC(1+ln B).

Divide by Q=N^(4/3)n^(2/3). Since

    ln B=(n ln2)/3+(2 ln n)/3-ln2-2 ln c,

we get the exact upper-bound expression

    2/c^2 + 2c ln2
    + c(33.74+8 ln n-12 ln2-24 ln c)/(2n)
    + N^(-1/3)n^(-2/3)(2+(3/2)ln B).

For 1<=c<=2, the first line is a convex function of c. Its maximum occurs at an endpoint and is at most 2+2 ln2 < 3.386295.

The second line is at most (33.74-12 ln2+8 ln n)/n. This decreases for n>=12 and at n=12 is less than 3.775458.

Finally ln B <= n ln2/2 for n>=12. To see this, discard -2 ln c and use log_2 n <= n/4+3/2, true at 12 and preserved thereafter because the right derivative is 1/4 and 1/(n ln2)<1/4. In fact the needed equivalent inequality is log_2 n<=n/4+3/2. Thus the last line is at most

    2^(-n/3)n^(-2/3)(2+(3/4)n ln2).

This is decreasing for n>=12: expand it as a positive multiple of 2^(-n/3)n^(-2/3) and of 2^(-n/3)n^(1/3), both decreasing there. At n=12 it is less than 0.098235. The total is less than 7.2600<7.3, proving the energy inequality.

Cauchy–Schwarz gives at least N^(2/3)/(29.2 n^(2/3)) different g-values on this interval. All are tournament counts and <=N, by the inherited one-apex realization and range g(v)<=v+1. For n=floor(log_2 X), N>X/2 and n<=ln X/ln2, whence

    A(X) >= X^(2/3)/[29.2 (2/ln2)^(2/3) (ln X)^(2/3)].

The numeric constant 29.2(2/ln2)^(2/3)<60 finishes the proof.

## Why split optimization alone stops here

This statement concerns the numerical right-hand side of THM-4887 (D1), not a lower bound on actual energy. Across all admissible C with P=N/(2C)>=8C, that right-hand side is at least

    (3 ln2/2) N^(4/3)n^(2/3).

Indeed let C0=(N/n)^(1/3). If C<=C0, the term 8P^2 alone is at least 2Q. If C>=C0, retain 12PC^2 ln(P/C)=6NC ln(N/(2C)). The function C ln(N/(2C)) increases on the admissible interval because its derivative is ln B-1>=ln8-1>0. Its value at C0 is at least C0(n ln2)/4 for n>=12, since ln(N/(2C0))=(n ln2)/3+(2 ln n)/3-ln2 >=(n ln2)/4. Therefore this term alone is at least (3 ln2/2)Q. Hence no further optimization of this same displayed bound removes its (ln X)^(2/3) loss or improves its 2/3 exponent. Actual arithmetic information remains necessary.

## Provenance and scope audit

Inherited: THM-4887-one-apex-energy-at-scale-two-to-the-n-over-three.md, (D1), zero-width constants, one-apex realization, range and Cauchy–Schwarz; THM-4890-one-apex-collision-balance.md gives exact symmetric block identity used by the checks. Repository snapshot supplied by parent: 0f3dad8c.

Closest proved mechanism: thin-strip lattice counting. Hostile: its smoothing floor at scale N^(1/3) cannot be read backward as linear exact energy. Corrected near miss: THM-4890 random micro model is not a bound on true micro values. Least-used sidecar here: the freedom to choose every admissible split in (D1), including below floor(n/3).

Explored but unproved: (i) dyadic recurrences of g(v)+tv; (ii) general small separable perturbations hc+f(h)+f(c); (iii) exact gcd counting in the balanced-window equation. None yields a stronger theorem here. In particular no random-model estimate has been substituted for actual energy. The existing multiplicative join closure does not by itself bootstrap a power-law counting lower bound to positive density.

## Reproduction

Run `python3 check.py > check.json` from this directory. The finite universe is all words in [2^(n-1),2^n), n=12..22. The script computes g by its defining recursion, the exact collision energy and support size, verifies the optimized inherited bound numerically, and independently exhausts the symmetric block identity on every word at each selected split. These checks support implementation, while the all-n proof is above. Decimal constants are deliberately padded; none of the computation is needed for the proof.
