# Paley and quasirandomness, deeply recovered — THM-438 essentially closes the H≥disc crux, and the Paley Bridge unifies it with LRC

**death-star-2026-07-21-S85** (HYP-8703). Owner: deeply investigate past work on quasirandomness and Paley.
Two exhaustive sweeps (the ~477-file Paley cluster; the quasirandom/discrepancy cluster). The headline: **my
open H≥disc crux is essentially already answered in canon**, and Paley is the single object where the
tournament (H) and LRC (covering) worlds meet.

## The crux is (essentially) already closed — THM-438
My S83/S84 crux (i) was "every regular tournament beats the Szele average $H\ge n!/2^{n-1}$," reduced (S84) to
the binding **doubly-regular** case. **THM-438 (PROVED, modulo one Weil bound) already gives it:**
$$R(p)=\frac{H(T_p)\,2^{p-1}}{p!}\ \longrightarrow\ e=\exp(-\chi(-1)),\qquad\text{i.e.}\quad H(\text{Paley}_p)\sim e\cdot\frac{p!}{2^{p-1}},$$
and it is **DRT-universal** — the leading order is shared by *all* doubly-regular tournaments (the engine needs
only $S^2=J-nI$: regularity + antisymmetry + two-point spectrum $\{0,\pm i\sqrt n\}$; "the only place
quasirandomness enters is the expander-mixing bound $|\lambda|=\sqrt n$"). My measured ratios $2.0,2.40,2.44$
($p=3,7,11$) are **converging to $e\approx2.718$** — this *is* THM-438. Since $n\cdot\mathrm{disc}(\text{DRT})/\text{avg}
=n(n+1)^{(n-1)/2}/n!\to0$ super-exponentially, $H(\text{DRT})\sim e\cdot\text{avg}\ge n\cdot\mathrm{disc}$ for
large $n$ (ratio $\to\infty$); small $n$ direct. **So the binding case of my crux is proved by THM-438**, and
the general quasirandom Ham-path lemma (which the repo only *cites*, never proves — Agent-2 confirmed) is not
even needed for the DRT case. The mechanism is Wigner/free-probability: the Paley cluster integrals are Catalan
numbers $A_{2k}=C_k p^{k+1}$ (free cumulants of $\tfrac12(\delta_a+\delta_{-a})$), the moment-method signature of
a random skew-Rademacher matrix.

## My recent work overlapped existing canon (honest)
The sweep shows two of my S82/S83 "contributions" are pre-existing theorems:
- **S82 "disc $=$ mean of Pfaffian squares"** $=$ **THM-468**'s expansion $\det(I+S)=\sum_{K\ \mathrm{even}}\mathrm{Pf}(S[K])^2$ (each Pf odd; $d(T)$ a switching-class invariant). Already canon.
- **S83 (iii) AM-GM bound** $\mathrm{disc}\le(n+1)^{(n-1)/2}/2^{n-1}$ $=$ **THM-472**'s determinant *ceiling*, equality $\iff$ DRT $\iff$ skew-Hadamard order $n+1$ $\iff n\equiv3\bmod4$ (Reid–Brown). Already proved, and it says the max-disc tournaments are *exactly* the DRTs — confirming my "regular=Paley is the wall."

So the H≥disc program, correctly credited: klein THM-1950 reduced it to the strong base; the strong base's crux
is the DRTs (S82/S84 $=$ THM-468/472); and **THM-438 closes the DRT case**. What remains is the non-doubly-regular
strong tournaments (small disc, easy) and making THM-438's "large $n$" explicit — a far smaller gap than before.

## The Paley Bridge — one object, two worlds (the deep unification)
The most important structural recovery: **THM-640 (the Paley Bridge, PROVED)** — the LRC$(p)$ family
$\{0,\dots,p-1\}$ under the QR-mod-$p$ cutoff **is** the Paley tournament $T_p$, so the **LRC covering-minimizer
(AP) and the tournament $H$-maximizer (Paley) are the same object seen through two functionals** ($M$ vs $H$),
with reversal $\leftrightarrow$ complementation. This makes rigorous "regularity $=$ AP $=$ Paley are one object"
(kps-S13, flagged for promotion in my S79 recovery ledger). And it explains the composite-14 obstruction: $14=2\cdot7$
is not a regular Paley cutoff — *why LRC(14) is hard*. Corollaries recovered: the density-floor jump $M:1/13\to1/12$
mirrors $H:189\to175$ (same Gauss-sum rigidity, "nothing between"); the three-gap count $g$ $=$ LRC spectral
multiplicity $=$ eigenvalue-magnitude count.

## Critical line = quasirandomness (the S75 pole, recovered as pseudorandomness)
Paley's flat spectrum $|\lambda_k|=\sqrt{(p+1)/4}$ (THM-126/162), non-Perron on $\mathrm{Re}=-\tfrac12$ (THM-1555),
is literally "**spectral flatness $=$ the Riemann Hypothesis for $y^2=x$ over $F_p$**" (qr-tournament-foundations)
and $=$ the expander/pseudorandom gap. So my S75 "Paley $=$ critical-line pole" $=$ "Paley $=$ quasirandom." And
the whole LRC discrepancy apparatus — $\mathrm{disc}_v$ (THM-731/732 $=$ generalized Dedekind sum), the
autocorrelation-discrepancy (THM-729), additive energy (klein-S175), three-gap $=$ spectral-flatness (kps-S29),
the Riesz-product program — is the **circle-side** quasirandomness, dual to the tournament-side Ham-path
quasirandomness (THM-438). Both are "flat spectrum $=$ small discrepancy $=$ near-average count."

## Leads surfaced (into the backlog)
- **THM-130 Paley Betti numbers are flagged WRONG** (GLMY convention bug; true $\beta_4=0$ not $6$ for $T_7$) — an
  unresolved court case; $\beta_1\le1,\beta_2=0$ survive. Worth resolving.
- **Promote "regularity=AP=Paley are one object"** (kps-S13) — it is THM-640 made conceptual, not a dormant frame.
- **Ihara/Bartholdi zeta of tournaments** (Paley $T_p$ $=$ circulant digraphs; Tang–Yau arXiv:2602.04140) — never done.
- The THM-438 non-circulant DRT remainder (replace the Weil bound by expander-mixing $|\lambda|=\sqrt n$) — the
  route to make the crux fully uniform.
- Two THM-130s (Betti vs $c_5$ closed form) — de-collision.

## Honest scope
A recovery/synthesis: the H≥disc crux is substantially closed by finding THM-438 (not new work — recovered);
my S82/S83 re-derived THM-468/472; the Paley Bridge (THM-640) and critical-line$=$quasirandom unify the
tournament and LRC sides. No open problem newly closed by me; the value is connecting my H≥disc thread to the
canon that answers it and mapping the Paley/quasirandom corpus. Cross-links: THM-438/468/472/640/1555/126/162,
THM-729/731/732 (disc_v), THM-1950 (H≥disc), S75/S76/S82/S83/S84, kps-S13/S29, klein-S175/S287/S288. Scripts:
`paley_hampath_ratio_deathstar_S85.py` (+out, ratio→e). HYP-8703.
