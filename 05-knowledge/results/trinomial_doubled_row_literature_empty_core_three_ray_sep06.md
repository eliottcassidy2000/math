# Doubled factorial rows: bounded literature recovery and transport obstruction

Status: **CITED SOURCE STATEMENTS; ANALYTICALLY VERIFIED MAP OBSTRUCTION.**
The general carry-normalized cross-return sign claim remains **OPEN** here.
Read-only research subtask for the overnight session, 2026-09-06; no theorem ID.

Inheritance: THM-4436, `01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md`, proves simple negative roots for each complete factorial row. Its Laguerre polynomials are finite diagonal-operator symbols; it does not identify successive return rows as one orthogonal sequence. Incoming carry correction from `certificate_audit` is retained below.

## 1. Exact source-to-target map

Write the THM-4436 row as

\[
P(t)=\sum_{j=0}^h\frac{m!t^j}{(x+(B-A)j)!\,(Bh+r-Bj)!\,(z+Aj)!},
\quad m=x+Bh+r+z.
\]

Set \(d=B-A\), \(u_k=-h+(k-r)/B\), \(0\le k<B\), and form the lower-parameter multiset

\[
L=\{(x+k)/d:1\le k\le d\}\cup\{(z+k)/A:1\le k\le A\}.
\]

Remove one occurrence of \(1\) from \(L\), possible because \(0\le z<A\), to obtain \(L'\). Direct factorial grouping gives the terminating identity

\[
P(t)/P(0)={}_BF_{B-1}\!\left(\begin{matrix}u_0,\ldots,u_{B-1}\\L'\end{matrix};
\frac{(-1)^BB^B}{d^dA^A}t\right).
\]

One upper parameter is \(-h\); every upper parameter is negative when \(h\ge1\). All lower parameters are positive. This preserves the entire coefficient row and zero set up to an explicit nonzero scaling of the variable and polynomial. Merely retaining the label “hypergeometric” destroys these parameter signs and the relation between returns.

[Wolfs, *Applications of multiple orthogonal polynomials with hypergeometric moment generating functions*, arXiv:2401.08312v2](https://arxiv.org/pdf/2401.08312), Theorems 2.13 and 3.11, gives type-II rows with upper parameters \(-N,\mathbf b+\mathbf n+1\), lower parameters \(\mathbf a+1\), and \(\mathbf a,\mathbf b>-1\). Thus the nonterminating upper parameters there are positive. Our direct map fails that hypothesis: it has additional negative upper parameters. Theorem 3.13 also requires a step-line index and explicit cross-parameter inequalities; it proves individual simple positive roots, not doubled-row separation. A transformation producing a common admissible measure system would be an additional theorem, not an automatic consequence of this identity.

## 2. What neighboring-zero theorems actually supply

[Kozhan, *Zeros of multiple orthogonal polynomials: location and interlacing*, Theorem 5.1](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/blms.70281), assumes normal neighboring indices \(\mathbf n,\mathbf n+\mathbf e_j\) for the same measure system. With the first polynomial real-rooted, their zeros interlace precisely when \(\mathbf n\) remains normal after every real double Christoffel transform \((x-z_0)^2\boldsymbol\mu\). Neither a common measure system nor neighboring indices has been supplied for our returns. Iterating neighboring interlacing also loses the parity of the number of second-row zeros between first-row zeros, which is necessary for a pointwise sign statement.

[Arvesú–Driver–Littlejohn, *Interlacing of zeros of Laguerre polynomials of equal and consecutive degree*, Theorem 2.1](https://arxiv.org/pdf/2009.10206), compares \(L_n^{(\alpha)}\) and \(L_{n+1}^{(\alpha+1)}\), with \(\alpha>-1\) and no common zeros. Even here an extra point \(n+1\) can replace a polynomial zero in the interlacing pattern. It supplies neither arbitrary parameter shifts nor degree doubling.

The bounded primary-source search recovered no theorem whose checked hypotheses imply our proposed cross-return sign. This is a retrieval stopping reason, not a claim that no such theorem exists.

## 3. Two hostile controls and retained sidecars

The incoming support \((-15,1,9)\) at return mass \(8\) has first tuple \((2,3,2,0,1,1)\) and doubled tuple \((2,3,5,0,0,1)\). Direct factorial reconstruction gives

\[
P=56(1+10t+t^2),\quad
Q=16+10920t+400400t^2+1681680t^3+720720t^4+8008t^5.
\]

Here \(Q\bmod(t^2+10t+1)=-47087024-466126752t\). Its values at \(-5\pm\sqrt{24}\) are both positive: the smaller is \(2283546736-466126752\sqrt{24}>0\), as follows by squaring the positive terms. Thus canonical \(Q<0\) fails; the candidate uses \(t^{-\epsilon_z}Q\), where \(\epsilon_z=\lfloor2z/A\rfloor=1\). Preserve both carries:

\[
h_2=2h+\epsilon_r+\epsilon_z,\quad r_2=2r-B\epsilon_r,\quad
z_2=2z-A\epsilon_z,\quad x_2=2x-(B-A)\epsilon_z,
\]

with \(\epsilon_r=\lfloor2r/B\rfloor\). For an actual positive Laurent first row, \(b=Ax-(B-A)z>0\), so \(x_2>0\); an arbitrary abstract factorial row requires its own doubled-tuple admissibility check. Losing the monomial changes the sign on the negative real axis.

Incoming `e5a45df05f`, [the independent complete-row path proof and sign bank](overnight3_20260906_moments_root_geometry.md), Section 6, independently retains this same factor and certifies the corrected sign at all 1,015 roots in its declared 221-support bank. That finite evidence does not supply the general cross-row theorem. Incoming THM-4440 is a **RESERVED / UNPROVED EMPTY STUB** for a signed-duplication direction and is not a dependency here.

Even genuine ordinary orthogonality is insufficient by itself. Take the monic recurrence \(p_{n+1}=(t-a_n)p_n-p_{n-1}\), \(p_0=1,p_1=t\), with \(a_3=3\) and all other \(a_n=0\). Positive recurrence coefficients give a positive orthogonality measure by [Favard's theorem, DLMF §18.2(viii)](https://dlmf.nist.gov/18.2.viii). Yet \(p_2=t^2-1\) and \(p_4=(t-3)(t^3-2t)-(t^2-1)\) satisfy \(p_4(1)=2\), \(p_4(-1)=-4\). So even an orthogonality realization would need an additional special structure for the desired same-sign conclusion.

The cheapest next decisive test is an exact remainder or signed-subresultant identity for the carry-normalized pair, retaining the tuple and both carry bits. The successful no-carry width-four family from the algebra siblings is a positive control; the support above is the normalization hostile. No further support census was run here. The hypergeometric identity was independently checked coefficient by coefficient on these two hostile rows and the width-four tuples \((1,2,2,0,0,1)\), \((1,2,4,0,0,2)\), totaling 17 exact rational identities.
