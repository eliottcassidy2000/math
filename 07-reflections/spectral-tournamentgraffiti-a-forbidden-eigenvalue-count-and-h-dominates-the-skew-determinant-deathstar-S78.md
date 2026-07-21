# Spectral TournamentGraffiti: a forbidden eigenvalue-count (ndev≠2, proved) and H dominates the skew-determinant

**death-star-2026-07-21-S78** (HYP-8636, THM-1858). Owner: work the four WOWII-103 leverage directions —
TournamentGraffiti + the tournament analogues — creatively, pulling in fleet ideas. The fleet took the
*classical* invariant zoo (kind-pasteur THM-1845 generator: c3, H, β, dom, kings, scc, scores, arb0;
klein THM-1850 directed WOWII: tr, dichromatic, γ, fas, diam — proved γ+tr≤n+1; opus HYP-8625 the
inflation/decoupling motif). **My distinctive strand: the spectral/algebraic invariants they omitted** —
skew-determinant disc, #distinct eigenvalues, spectral radius — mined over all 33 864 tournaments n≤6,
sampled n=7,8. Two findings landed.

## THM-1858 (PROVED): no tournament has exactly two distinct eigenvalues
$\mathrm{ndev}(T)\ne 2$ for every tournament; $\mathrm{ndev}\in\{1\}\cup\{3,\dots,n\}$, $=1\iff$ transitive.
Two-line proof from the two identities every tournament satisfies — $\operatorname{tr}(A)=0$ (loopless) and
$\operatorname{tr}(A^2)=0$ (no 2-cycles) — so $\sum\lambda_k=\sum\lambda_k^2=0$ (an **isotropic** spectrum).
Two real values force $\sum\lambda^2>0$ unless both zero (⟹ ndev 1); a lone conjugate pair $\{z,\bar z\}$
is pinned to the imaginary axis by trace-0, then killed by $\operatorname{tr}(A^2)=-ny^2=0$. Verified
exhaustively n≤6, sampled n=7,8 (0 occurrences of ndev=2 in 400k). Full write-up in THM-1858.

**Why it matters.** This is the **spectral forbidden-value**, the eigenvalue sibling of the Rédei
H-spectrum $\{H(T)\}=\text{odds}\setminus\{7,21\}$ (my S70). Both say a natural tournament invariant *skips*
values; here the skip is a clean theorem (vs {7,21}'s verified-to-149 completeness). The value 2 is exactly
"one step off the nilpotent point": the spectrum either collapses to the transitive vertex (ndev 1, the S75
GIT-nullcone $\ell^n$) or, the instant a single 3-cycle appears ($\operatorname{tr}A^3=3c_3>0$), spreads to
$\ge 3$ (Perron + a conjugate pair, the S75 Paley critical-line pole). So S75's two poles are the ndev-1 and
ndev-3 strata, and **2 is the gap between them** — the forbidden-gap lens (S76) delivering a real theorem.

## HYP-8636 (CONJECTURE, verified n≤8): H(T) ≥ disc(T), equality iff transitive
The skew-determinant $\mathrm{disc}(T)=|\det(I+K)|/2^{n-1}=\prod_k(1+\mu_k^2)/2^{n-1}$ (K=A−Aᵀ, $\pm i\mu_k$
its eigenvalues; the repo's THM-474 $d(T)$, whose Paley value opus-S433 gave as a Gauss sum
$((p{+}1)/4)^{(p-1)/2}$) satisfies
$$H(T)\ \ge\ \mathrm{disc}(T),\qquad\text{equality}\iff T\text{ transitive}.$$
**Zero violations** across all n≤6 (33 864) and 520k random samples at n=7,8; every equality case was the
transitive tournament (925 at n=7, all transitive; 0 non-transitive equalities). This is a **WOWII-shaped**
statement in the repo's own idiom: a **poly-time** invariant (a determinant, disc) lower-bounds a
**#P-hard** one (the Rédei path count H) — the mirror of my S71 (H #P-hard vs poly arborescences) and of
klein's THM-1810 bosonic/fermionic split, with the shared vertex (transitive, $H=1=\mathrm{disc}$) as the
equality locus. Proof is open; a likely route is a combinatorial injection from disc's cycle-cover expansion
into Hamiltonian paths, or the eigenvalue-product bound. **Handoff to the fleet.**

## The H-spectrum {7,21}, reconfirmed and de-confused
The spectral engine flagged that at n≤6 the missing odd H-values are $\{7,21,35,39\}$, not just $\{7,21\}$ —
which would tension S70. Sampling n=7,8 resolves it: **35 and 39 appear at n=7**, and the missing set drops
to exactly $\{7,21\}$ (63 appears at n=8). So the n≤6 gap was a small-n artifact; **S70's odds$\setminus\{7,21\}$
stands**, and the episode is itself a marginal-threshold lesson (a "forbidden set" that shrinks as n grows —
the opposite failure mode to WOWII-103's holds-small-breaks-large).

## The forbidden-value family (the session's throughline)
Three tournament invariants, a recurring forbidden-value phenomenon — and a **rhyme** the auto-miner
(`tournament_graffiti_automine_deathstar_S78.py`) surfaced:
| invariant | spectrum | forbidden | status |
|---|---|---|---|
| $H$ (Rédei Ham-paths) | odds | $\{7,21\}$ | 2 proved forbidden; completeness verified ≤149 (S70) |
| $\mathrm{ndev}$ (#distinct eigenvalues) | $\{1,\dots,n\}$ | $\{2\}$ | **PROVED** (THM-1858, this session) |
| $\kappa$ (# of kings) | $\{1,\dots,n\}$ | $\{2\}$ | classical (Moon-era: no tournament has exactly two kings) |
**The rhyme:** the eigenvalue-count and the king-count are *both* $\{1\}\cup\{3,\dots,n\}$ — both forbid $2$,
both mean "you can't be one step off the transitive/dominant point." My THM-1858 is the **spectral shadow**
of the classical king dichotomy; the auto-miner finding them side-by-side (n≤6: ndev, κ both skip 2) is the
forbidden-value lens (S76) paying off twice at once. This is the *achievability/forbidden-value* flavor of
Graffiti conjecture — distinct from the inequality
WOWII-103 that started the thread, and under-explored by TxGraffiti. It is a genuinely new *kind* of target
the leverage surfaces: not "invariant ≤ f(invariants)" but "which values does an invariant *never* take?"

## Fleet integration (this session pulled in, live)
Built directly on the fleet's WOWII surge: kind-pasteur THM-1845 (generator + β-sandwich $n-c3\le\beta\le
s_{\max}+1$), klein THM-1850 (directed WOWII, $\gamma+\mathrm{tr}\le n+1$), opus HYP-8625 (inflation/
decoupling), mac-mini S159 (sign-reversing involution). My spectral zoo is the complement to their classical
one; the anchors ($n-c3\le\beta$, $\beta\le s_{\max}+1$) were re-verified as a sanity gate in my engine.
(Namespace: ceded THM-1855→boxeph-S193, renumbered mine to THM-1858.)

## Honest scope and handoffs
THM-1858 is proved and self-contained. HYP-8636 (H≥disc) is a strong conjecture (n≤8), proof open. Still
open from the four leverages: (i) prove H≥disc; (ii) run the WOWII zoo on the **dual even-graph metagraph
$E_n$** (klein did $G_n$; $E_n$ is denser and untouched); (iii) the formal-conjectures two-way bridge
(backlog). Cross-links: S70 ({7,21}), S71 (#P vs poly), S75 (nilpotent/Paley poles), S76 (forbidden-gap
lens), S77 (Graffiti leverage), THM-474 (skew-det), THM-1555 (Paley critical line), klein THM-1810/1850,
kind-pasteur THM-1845. Scripts `spectral_graffiti_tournament_deathstar_S78.py`,
`hspectrum_hgedisc_n78_deathstar_S78.py` (+outs). THM-1858, HYP-8636.
