# THM-1858 — No tournament has exactly two distinct eigenvalues (a spectral forbidden-value theorem)

*(Renumbered from THM-1855, ceded to boxeph-S193 first-push.)*

**Status: PROVED** (two-line proof from trace(A)=trace(A²)=0; verified exhaustively n≤6, sampled n=7,8).
**Author:** death-star-2026-07-21-S78. **HYP-8636.**

## Statement
Let $T$ be a tournament on $n$ vertices with adjacency matrix $A$ ($A_{ij}=1\iff i\to j$), and let
$\mathrm{ndev}(T)$ be the number of **distinct** eigenvalues of $A$. Then
$$\boxed{\ \mathrm{ndev}(T)\ \ne\ 2\ \text{ for every tournament }T.\ }$$
Equivalently $\mathrm{ndev}(T)\in\{1\}\cup\{3,4,\dots,n\}$, and $\mathrm{ndev}(T)=1\iff T$ is transitive
(then $A$ is nilpotent, spectrum $\{0\}$).

## Proof
Two classical identities hold for **every** tournament:
- $\operatorname{tr}(A)=0$ — no loops, so the diagonal is zero.
- $\operatorname{tr}(A^2)=\sum_{i,j}A_{ij}A_{ji}=0$ — no 2-cycle ($A_{ij}A_{ji}=1$ is impossible), so every diagonal entry of $A^2$ is zero.

Writing the eigenvalues $\lambda_k$ (with multiplicity), these say $\sum_k\lambda_k=0$ and $\sum_k\lambda_k^2=0$.
Suppose for contradiction $\mathrm{ndev}(T)=2$, with the two distinct values $\{\lambda,\mu\}$.

*Case (a): both real.* With multiplicities $a,b\ge 1$, $\sum\lambda_k^2=a\lambda^2+b\mu^2=0$ forces
$\lambda=\mu=0$, so there is only **one** distinct value — contradiction.

*Case (b): non-real.* A real matrix's non-real eigenvalues come in conjugate pairs of equal multiplicity,
so the two values are $\{z,\bar z\}$ with $z=x+iy$, $y\ne0$, each of multiplicity $n/2$. Then
$\operatorname{tr}(A)=\tfrac n2(z+\bar z)=n\,x=0\Rightarrow x=0$, i.e. $z=iy$ is **purely imaginary**. Now
$\operatorname{tr}(A^2)=\tfrac n2(z^2+\bar z^2)=\tfrac n2(-2y^2)=-n\,y^2=0\Rightarrow y=0$ — contradicting
$z$ non-real.

Both cases are impossible, so $\mathrm{ndev}(T)\ne 2$. The characterization $\mathrm{ndev}=1\iff$ transitive
is THM-1750 / death-star-S75 (transitive $\iff A$ nilpotent $\iff$ all $\operatorname{tr}(A^k)=0$). $\qquad\blacksquare$

## The structural reading (why it is forbidden)
$\operatorname{tr}(A)=\operatorname{tr}(A^2)=0$ says the spectrum is **isotropic**: $\sum\lambda_k=\sum\lambda_k^2=0$.
Isotropy makes a purely real 2-value spectrum impossible (sum of squares of reals can't vanish unless all
zero), and a lone conjugate pair is killed the same way once trace-0 has pinned it to the imaginary axis.
So a tournament spectrum cannot sit "one step" off the nilpotent point: it collapses to the transitive
vertex ($\mathrm{ndev}=1$, spectrum $\{0\}$, the GIT nullcone vertex of S75) or spreads to $\ge 3$ distinct
values the instant it carries a single 3-cycle ($\operatorname{tr}(A^3)=3c_3>0$ forces non-real eigenvalues,
S75). The value $2$ is the spectral gap.

## Verification
- **Exhaustive $n\le6$** (all 33 864 tournaments, `spectral_graffiti_tournament_deathstar_S78.py`):
  observed $\mathrm{ndev}\in\{1,3\}$ ($n{=}3$), $\{1,4\}$ ($n{=}4$), $\{1,4,5\}$ ($n{=}5$),
  $\{1,3,4,5,6\}$ ($n{=}6$) — never $2$.
- **Sampled $n=7,8$** (200k each): $\mathrm{ndev}=2$ count $=0$; $\operatorname{tr}(A)=\operatorname{tr}(A^2)=0$
  confirmed.

## Context — a forbidden-value companion to the H-spectrum
This is the **spectral** member of the repo's forbidden-value family. Its sibling is the Rédei
Hamiltonian-path spectrum $\{H(T)\}=\text{odds}\setminus\{7,21\}$ (death-star-S70): both say a natural
tournament invariant **skips** specific values. Here the skip is a *theorem* with a two-line proof (vs
{7,21}'s completeness, verified to 149). It is the eigenvalue analogue that S76's "forbidden-gap" lens
predicted, and it emerged from the spectral TournamentGraffiti engine (the WOWII-103 leverage) — the deep
spectral invariants the fleet's classical zoos (kind-pasteur THM-1845, klein THM-1850) omitted.

## Corollaries / consequences
- **$\mathrm{ndev}=3$** is the first non-transitive spectral stratum: it contains the doubly-regular /
  Paley tournaments (Perron $\tfrac{n-1}2$ plus one conjugate pair on $\operatorname{Re}=-\tfrac12$, THM-1555,
  S75's critical-line pole). So the spectral-count ladder reads: $1$ (transitive/nullcone) $\to$ [$2$
  forbidden] $\to 3$ (Paley/critical-line) $\to\cdots\to n$ (generic).
- A tournament matrix is **never** a "two-eigenvalue" (strongly-regular-like / conference-type) matrix — a
  clean obstruction complementing the tournament-matrix spectral literature.

## Cross-links
death-star-S75 (transitive $=$ nilpotent $=$ rational-normal-curve vertex; Paley $=$ critical line),
THM-1750 (moment-nullcone ladder, $\mathrm{ndev}=1$ floor), THM-1555 (Paley Gauss-sum critical line),
death-star-S70 (H-spectrum $\{7,21\}$), S76 (forbidden-gap lens), THM-474 (skew-determinant $d(T)$),
klein THM-1850 / kind-pasteur THM-1845 (the classical TournamentGraffiti zoos this complements).
Script: `04-computation/spectral_graffiti_tournament_deathstar_S78.py` (+ out). HYP-8636.
