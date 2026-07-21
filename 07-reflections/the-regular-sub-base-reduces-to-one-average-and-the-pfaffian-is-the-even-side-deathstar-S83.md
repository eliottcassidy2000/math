# The regular sub-base of H ≥ disc reduces to one clean average; the Pfaffian injection is the even/odd duality

**death-star-2026-07-21-S83** (HYP-8698). Owner: work the Pfaffian injection and the regular sub-base (S82's
two handles toward klein THM-1950's open strong base $H(C)\ge\max(1,s(C))\mathrm{disc}(C)$ for HYP-8636). The
regular sub-base now reduces to a **single, clean, well-evidenced claim**; the Pfaffian injection resolves to
the even/odd (cycle-space vs OCF) duality.

## The regular sub-base, PROVED for n ≥ 7 modulo one average
For a regular tournament ($n$ odd, all out-degrees $(n-1)/2$, so $s=n$), the base is
$H(\mathrm{reg})\ge n\cdot\mathrm{disc}(\mathrm{reg})$. It follows from the chain
$$H(\mathrm{reg})\ \overset{(i)}{\ge}\ \frac{n!}{2^{n-1}}\ \overset{(ii)}{\ge}\ \frac{n\,(n+1)^{(n-1)/2}}{2^{n-1}}\ \overset{(iii)}{\ge}\ n\cdot\mathrm{disc}(\mathrm{reg}),$$
of which **(ii) and (iii) are proved and (i) is the single crux:**

- **(iii) — PROVED (AM–GM).** $\mathrm{disc}(\mathrm{reg})=\prod_{j}(1+\mu_j^2)/2^{n-1}$ where $\pm i\mu_j$ are
  the eigenvalues of $K=A-A^\top$, subject to the fixed energy $\sum\mu_j^2=\operatorname{tr}(-K^2)/2=\binom n2$
  over the $(n-1)/2$ pairs. By AM–GM (concavity of $\log(1+x)$), $\prod(1+\mu_j^2)$ is maximized at equal
  $\mu_j^2=\binom n2/\tfrac{n-1}2=n$, so $\mathrm{disc}(\mathrm{reg})\le(n+1)^{(n-1)/2}/2^{n-1}$, with equality
  iff **doubly-regular** ($\mu_j^2=n$ ∀$j$, i.e. Paley — all $K$-eigenvalues $\pm i\sqrt n$). Verified tight at
  Paley-3,7,11.
- **(ii) — PROVED (elementary).** $(n-1)!\ge(n+1)^{(n-1)/2}$ **fails** at $n=3,5$ but holds at $n=7$
  ($720\ge512$), and the ratio $(n-1)!/(n+1)^{(n-1)/2}$ is increasing for $n\ge7$ ($1.4,4,\dots$), so it holds
  for all $n\ge7$.
- **(i) — THE CRUX (conjecture, strongly evidenced): every regular tournament has at least the Szele average
  number of Hamiltonian paths, $H(\mathrm{reg})\ge n!/2^{n-1}$.** Verified **exhaustively** $n=3$ ($H=3\ge1.5$)
  and $n=5$ (all 24 regular, $H=15\ge7.5$), and on Paley/rotational/random-regular samples at $n=7$
  ($\min H=171\ge78.8$) and $n=9$ ($\min H=3243\ge1417.5$) — always with large margin.

So: **for $n\ge7$ the regular sub-base is proved iff (i) holds; $n=3,5$ hold directly** (exhaustive: ratio
$H/(n\,\mathrm{disc})=1$ at $C_3$, $3$ at $n=5$). The doubly-regular (Paley) tournaments are the tightest — max
disc, ratio $\to3.375$ at Paley-7 — exactly the "maximally-symmetric configuration is the wall" pattern
(S75/S76). **This reduces the whole regular crux to a single named statement about Hamiltonian-path counts** —
"regular tournaments beat the average" — plausibly known (Moon/Alon/Busch circle) and far more tractable than
the eigenvalue-product original.

## The Pfaffian injection is the even/odd (cycle-space vs OCF) duality
The base (at $s\le1$) is $2^{n-1}H\ge\det(I+K)=\sum_{S\ \mathrm{even}}\mathrm{Pf}(K[S])^2$ (S82). Confirmed with
room: $\min(2^{n-1}H-\det(I+K))=112$ ($n=5$), $416$ ($n=6$) over all strong tournaments. The structural reading:
- $\mathrm{Pf}(K[S])$ is the signed count of **oriented perfect matchings** of $S$; $\det(I+K)=\sum_S\mathrm{Pf}^2$
  is a signed count of disjoint **even** cycle covers — an **even/cycle-space** object.
- $H=I(\Omega,2)$ (Grinberg–Stanley) is built from the **odd**-cycle collection $\Omega$ — an **odd/OCF** object.
- So $H\ge\mathrm{disc}$ is literally **"the odd (OCF) count dominates the even (Pfaffian) count"** — the same
  even/odd, cut/cycle, $E_n$-vs-$\Omega$ duality as my S80 (even-graph metagraph $=$ cycle space) and the
  cut$\oplus$cycle GF(2) split. The per-subset injection $\mathrm{Pf}(K[S])^2\le(\text{Ham paths compatible with }S)$
  remains open, but the object it must exploit is now named: an odd-structure/even-structure comparison, not an
  eigenvalue estimate. (Tested $\mathrm{Pf}(K[S])^2\le H(T[S])H(T[V\setminus S])$ — not clean; the right
  compatibility is subtler.)

## Status and honest scope
Genuine progress on HYP-8636's open base: (a) the regular sub-base is **reduced to one clean average (i)**, with
(ii),(iii) proved — so the "regular is the wall" crux is now a single tractable Ham-path statement; (b) the
Pfaffian injection is re-identified as the even/odd duality (aggregate confirmed, per-subset open). No new
theorem beyond the reduction; (i) is the residual, and if it is known/proved the regular crux closes for $n\ge7$.
Cross-links: klein THM-1950 (the reduction), S82 (Pfaffian-mean, regular crux), S80 (even/odd $E_n$ duality),
THM-002/209 ($H=I(\Omega,2)$), THM-474 (skew-det), S75/S76 (Paley-is-the-wall). Script
`h_ge_disc_regular_subbase_deathstar_S83.py` (+out). HYP-8698.
