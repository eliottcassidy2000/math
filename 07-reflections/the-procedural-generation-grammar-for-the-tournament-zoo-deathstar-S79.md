# A procedural generation grammar for the tournament zoo — objects × lenses × invariants × operations, and the gaps it exposes

**death-star-2026-07-21-S79** (part of the zoo-atlas sweep; companion to the census/lost-threads atlas).
Owner: keep adding to the zoo, recover every past thread, **procedurally generate** new frames / methods /
angles / things to compute, and find the gaps between and around them. This file is the *generative engine*;
the companion atlas is the *census*. The engine says: the research surface is the product

$$\textbf{OBJECTS}\ \times\ \textbf{LENSES}\ \times\ \textbf{INVARIANTS}\ \times\ \textbf{OPERATIONS},$$

and a **gap** is any cell of that product not yet filled. Below is each factor, the generation rules that
combine them, and a concrete batch of targets each rule emits (each is backlog-able).

## The four factors (the alphabet)

**OBJECTS** $O$: tournaments $G_n$; even graphs $E_n$ (dual); merged metagraph $G_n/\mathbb Z_2$;
tilings/staircase $\delta_{n-2}$; LRC speed sets; GMC polynomials; the odd-cycle collection $\Omega(T)$;
Paley/QR tournaments; the tournament pencil $M(t,u,v)$; the Cayley–Dickson tower.

**LENSES** $L$ (the 12 thinking-ways, S76): reify · nullcone · two-poles · complexity-as-depth ·
self-duality/apolarity · composition-mode · tangent-at-origin · Reynolds/charge-killing · critical-line ·
forbidden-gap · big-stabilizer-boundary · substrate-recursion.

**INVARIANTS** $I$: $H$ (Rédei), the cycle-trace spectrum $p_k=\operatorname{tr}A^k$, $\#$Ham-cycles,
$\beta$ (largest transitive sub), $\gamma$ (domination), kings, scc, fas/Slater, dichromatic, $\sum a_r$
(arborescences), $\mathrm{ndev}$, spectral radius $\rho$, $\rho_2$, skew-det $\mathrm{disc}=d(T)$,
$\det(I+A)$, $\operatorname{per}(A)$, $I(\Omega,x)$ (independence poly), $\chi(E_n)$, width, Betti.

**OPERATIONS** $\mathrm{Op}$: complement $T^{\mathrm{op}}$ ($A\mapsto A^{\!\top}$); ordinal sum $\oplus$;
vertex-deletion (OCR, $n\!\to\!n\!-\!1$); the dual $G_n\!\leftrightarrow\!E_n$; Mode A ($n\!\to\!n\!-\!1$,
hypotenuse) / Mode B ($n\!\to\!n\!-\!2$, both legs); transpose/relabel; waggly/wiggly tile-flips;
tensor / Clebsch–Gordan.

## The five generation rules (the grammar), with emitted targets

### Rule 1 — $I \times L=\text{forbidden-gap}$: "which values does invariant $I$ never take?"
The single most productive rule this month. Confirmed cells: $H\to\text{odds}\setminus\{7,21\}$ (S70);
$\mathrm{ndev}\to\{1\}\cup\{3..n\}$, forbids $2$ (THM-1858); kings $\to\{1\}\cup\{3..n\}$, forbids $2$
(classical). **Emitted (open) cells — one research question each:** the achievable spectrum of $\gamma$,
dichromatic, fas, $\sum a_r$, $\det(I+A)$, $\operatorname{per}(A)$, and each $p_k$. *Procedure:* run the
full-zoo battery, read off gaps, then seek the proof/mechanism (as ndev≠2 got its 2-line isotropy proof).

### Rule 2 — $I \times L=\text{critical-line}$: "does $I$ have a Paley/Gauss-sum extremal?"
Confirmed: $\mathrm{ndev}=3$ is the Paley stratum; $\mathrm{disc}(\text{Paley}_p)=((p{+}1)/4)^{(p-1)/2}$ is a
Gauss sum (opus-S433); non-Perron eigenvalues on $\operatorname{Re}=-\tfrac12$ (THM-1555). **Emitted:** is
there a critical-line/Gauss-sum closed form for $\sum a_r(\text{Paley})$? for $I(\Omega,x)$ at Paley? for
$\det(I+A)$? Paley maximizes/extremizes *which* invariants (the "Paley is the wall" program)?

### Rule 3 — $I \times \mathrm{Op}=\oplus$: "the composition algebra of invariants."
Under ordinal sum $A_{T\oplus S}=\begin{smallmatrix}[A_T&J]\\[0&A_S]\end{smallmatrix}$ (block-triangular):
$H$ is **multiplicative** (Rédei monoid, S70); eigenvalues are the **union** so
$\mathrm{ndev}(T\oplus S)=|\mathrm{spec}(T)\cup\mathrm{spec}(S)|$ and $p_k$ **add**. **Emitted:** the
$\oplus$-law for every remaining invariant — is $\mathrm{disc}$ multiplicative (its $K$-block is *not*
triangular, so maybe not)? $\sum a_r$? $\gamma$, kings (additive-ish)? Each $\oplus$-law is a mini-theorem;
together they make the invariant vector a **ring homomorphism from the $\oplus$-monoid** — the algebraic
skeleton the pencil $M(t,u,v)$ (THM-1760) only glimpsed for one invariant.

### Rule 4 — $I \times \mathrm{Op}=\text{complement}$: "the functional equation of $I$."
$M(t,u,v)$ has a complement functional equation (THM-1760). Spectral invariants are transpose-invariant so
$\mathrm{ndev},\rho,\mathrm{disc},H$ are complement-**invariant** (quick: $A^{\mathrm{op}}=A^{\!\top}$).
**Emitted:** the *non*-invariant ones carry data — $\gamma(T)$ vs $\gamma(T^{\mathrm{op}})$ (in- vs
out-domination, klein flagged), fas, dichromatic, $\sum a_r$ under complement. The complement-**odd** part of
each invariant is the SC/NS-sensitive content (spine/ribs/sea).

### Rule 5 — $O \times \mathrm{Op}=\text{dual}$: "run every $G_n$ result on $E_n$ (and vice-versa)."
The biggest single gap. Klein ran the WOWII zoo on $G_n$; **nobody ran it on the dual even-graph metagraph
$E_n$** ($V(E_n)=2,3,7,16,54$; denser; $\chi$ grows faster). **Emitted:** the entire invariant battery on
$E_n$; the forbidden-value program on $E_n$; the $G_n\leftrightarrow E_n$ bridge matrix's spectrum; whether
$E_n$ has its own $\{7,21\}$-type gap. Also Mode A vs Mode B: every $n\to n-1$ result has an untested
$n\to n-2$ shadow.

### Rule 6 — $O \times L$ (the S76 grid, extended): the blank cells.
S76 tabulated 10 objects × lenses with many blanks. **Emitted:** the least-developed cells were
self-duality/apolarity and the metagraph-spectrum; procedurally, *every blank in that grid is a prompt* —
e.g., "composition-mode of the Cayley–Dickson tower," "tangent-at-origin of the LRC speed form," "Reynolds
lens on arborescences."

## The bosonic/fermionic doubling (a cross-cutting generator)
klein THM-1810: determinant/Vandermonde (fermionic, cancels, transitivity survives) vs permanent/Gaussian
(bosonic, no cancellation). This **doubles every determinant invariant**: $\det(I+A)\leftrightarrow
\operatorname{per}(I+A)$; $\mathrm{disc}=$ skew-**det** $\leftrightarrow$ skew-**per** (hafnian); Matrix-Tree
$\sum a_r$ (a determinant) $\leftrightarrow$ its permanent analogue. **Emitted:** compute every permanent-twin
and ask Rule-1 (forbidden values) and Rule-3 ($\oplus$-law) of it. The permanent side is where $H$ (the
#P-hard flagship) lives, so the permanent-twins are the ones most likely to rhyme with $H$.

## What this session already emits and computes
The full-zoo battery (`tournament_zoo_full_deathstar_S79.py`, ~20 invariants incl. the new
$p_k$-spectrum, dichromatic, fas, $\sum a_r$, ham-cycles, scc, $\gamma$) executes Rule 1 (forbidden values)
and the WOWII miner (Rule "$I\times I$") over all n≤6 at once, saving a **reusable dataset**. The candidate
$H\ge\mathrm{disc}$ (HYP-8636) is the $I\times I$ rule's tightest current output.

## The meta-gap
The repo has **objects and invariants in abundance but the OPERATION rows are nearly empty**: only $H$ has a
full $\oplus$-law, only $M(t,u,v)$ has a complement functional equation, only $G_n$ has a run zoo. The
highest-leverage move is not a new invariant but **filling the operation columns** — the $\oplus$-algebra, the
complement functional-equations, and the $E_n$ dual — because each fills a whole row of cells at once.

## Cross-links
S76 (the 12 lenses + the object grid), S70 ($\{7,21\}$), S75 (nilpotent/Paley poles), S78/THM-1858
(ndev≠2), HYP-8636 (H≥disc), klein THM-1760 (pencil complement eq), klein THM-1810 (bosonic/fermionic),
opus-S433 (Paley Gauss-sum), THM-1750/1775 (moment-nullcone ladder), even-graphs-as-first-class ($E_n$).
Companion: the invariant-atlas census (this sweep). Script `tournament_zoo_full_deathstar_S79.py`.
