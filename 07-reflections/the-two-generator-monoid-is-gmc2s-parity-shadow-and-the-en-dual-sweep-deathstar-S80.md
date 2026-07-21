# The a/b two-generator functional monoid is GMC(2)'s parity shadow — and the E_n dual sweep

**death-star-2026-07-21-S80** (HYP-8653). Owner: work the E_n dual sweep, and think the two-generator
functional monoid and its relation to GMC(2). The two requests are one theme — **"even" / parity /
charge-0 / cycle-space** — and this note ties them together. Builds on kps THM-1880 (the a/b frame), the
GMC(2) toral×radial factorization, klein THM-1810 (bosonic/fermionic), my S78 (THM-1858).

## Part A — the E_n dual sweep (the flagged gap, now run)
Constructed the even-graph metagraph $E_n$ (vertices = even graphs on $n$ vertices $=$ cycle space of $K_n$;
edges = single tournament-base-path tile-flip $=$ XOR with one fundamental cycle) and ran the invariant zoo,
the **dual** of the tournament metagraph $G_n$ (klein: $\alpha(G_n)=2,5,18$).

| $n$ | $V(E_n)$ | (A002854) | $\omega=\chi$ | density | $\alpha(E_n)$ | ecc | diam | WOWII-103 |
|---|---|---|---|---|---|---|---|---|
| 3 | 2 | ✓ | 2 | 1.00 | 1 | 1.0 | 1 | holds |
| 4 | 3 | ✓ | 3 | 1.00 | 1 | 1.0 | 1 | holds |
| 5 | 7 | ✓ | 5 | 0.76 | 2 | 1.86 | 2 | holds |
| 6 | 16 | ✓ | 10 | 0.75 | 4 | 2.0 | 2 | holds |

Findings: (i) $V=$ A002854 and $\omega=\chi=2,3,5,10$ **reproduce the repo's $E_n$** (validates the cycle-space
construction). (ii) **NEW: $\alpha(E_n)=1,1,2,4$** — the even-graph metagraph independence number, never computed
(the census had $\chi/\omega$ but not $\alpha$); it fits $\alpha(E_n)=2^{\,n-4}$ for $n\ge4$ (predicts $\alpha(E_7)=8$,
testable). Much smaller than $\alpha(G_n)=2,5,18$ because $E_n$ is **dense** (0.75) and small-diameter (1,1,2,2)
— the dual is a "small-world clique-rich" graph where $G_n$ is sparse. (iii) **WOWII-103 holds on $E_n$** (dual
to klein's "$G_n$ satisfies 103"). (iv) **Tile-basis dependence (structural finding):** the even-graph metagraph
is **not canonical** — with the *star-tree* fundamental cycles it is **bipartite** ($\omega=\chi=2$, $\alpha=1,2,4,10$);
with the *path-tree* (tournament base-path) tiles it is the **dense chordal** repo-$E_n$ ($\chi=2,3,5,10$). The
metagraph is a family indexed by the spanning tree; the tournament base-path is the distinguished one.

## Part B — the a/b monoid ⟨a,b⟩ is the parity (ℤ/2) shadow of GMC(2)'s toral layer
kps THM-1880: $a(x)=x+1$, $\bar a(x)=x-1$ (conjugate shift), $b(x)=x/2$ (symmetriser); the transitive
tournament's forms are $\mathrm{char}_A=x^n$ (nullcone monomial), $E_n=b(a^n+\bar a^n)$ (even), $O_n=b(a^n-\bar a^n)$
(odd), with the Pell identity $E_n^2-O_n^2=(x^2-1)^n$. Here is the exact dictionary to GMC(2)
($E[P^m]=L_s(\mathrm{CT}_u[\Lambda_s^m])$, charge $c=\deg_Z-\deg_W$, $E$ kills nonzero charge):

1. **$b$ (symmetriser $x/2$) $\;\leftrightarrow\;$ $\mathrm{CT}_u$ (toral charge-0 projection).** Both are the
   *trivial-isotypic projector* — the average that kills non-invariant modes. $b$ averages the two conjugate
   shifts $a,\bar a$ (charges $\pm1$) to the even part; $\mathrm{CT}_u=\frac1{2\pi}\int\!d\theta$ averages the
   torus to charge 0. **$b$ is exactly the $\mathbb Z/2$ (parity) version of the $U(1)$ projector $\mathrm{CT}_u$.**
   The a/b monoid grades by *parity* ($a^n$ split into even/odd degree); GMC(2) grades by *charge* ($\mathbb Z$);
   $b$ is the mod-2 shadow of the full charge-0 projection.
2. **$E_n$ (even) $/$ $O_n$ (odd) $\;\leftrightarrow\;$ bosonic $/$ fermionic (permanent $/$ determinant).** The
   Weyl involution $x\mapsto-x$ (the $a\!\leftrightarrow\!\bar a$ swap) splits into symmetric $E_n$ (no sign
   cancellation $=$ **permanent side**) and antisymmetric $O_n$ (signs cancel $=$ **determinant/Vandermonde
   side**). This is klein THM-1810 verbatim: **GMC(2) is hard because it lives on the bosonic ($E_n$/even/
   permanent) side where the alternating collapse is unavailable**, while transitivity ($O_n$/odd/determinant,
   the Vandermonde survival of THM-1805) is the easy cancelling side. So the a/b monoid *names the hard/easy
   split of GMC(2)* at the transitive vertex.
3. **The $\tfrac12$ ($b$) $=$ the half-integer structure of both GMC layers.** $b\!\circ\!a=(x+1)/2$ is THM-1555's
   half-dictionary; the $\tfrac12$ is the fiber fraction $(\tfrac12)_{n-2}/(n-2)!$, the Legendre exponent, the
   $\mathrm{Re}=-\tfrac12$ line. GMC(2)'s two layers are **Legendre** (toral $\mathrm{CT}_u$) and **Hermite**
   (radial $L_s$) orthogonal polynomials, bridged by $k!=(1)_k$ (THM-1620) — the same half-integer world in
   which $(E_n,O_n)$ are the **Chebyshev** pair. The a/b Chebyshev pair is the *finite trigonometric shadow* of
   the GMC Legendre/Hermite pair.
4. **The radial variable, made precise.** $a\cdot\bar a=(x+1)(x-1)=x^2-1$ is the a/b "metric," and $(E_n,O_n)$
   is its $(\cosh,\sinh)$ pair. On the GMC side the radial variable is $s=ZW=|Z|^2$ — the product of the two
   charge-conjugate generators, exactly as $x^2-1=a\bar a$ is the product of the two conjugate shifts. So
   **$x^2-1 \leftrightarrow ZW$**: the a/b monoid carries *both* the toral ($b\leftrightarrow\mathrm{CT}_u$) and
   the radial ($a\bar a\leftrightarrow ZW$) datum of GMC(2), in one polynomial shadow.

## Part C — the two "$E_n$'s" are one object, and it is GMC(2)'s toral shadow
There is a name collision worth reading as a *unification*: THM-1880's $E_n=((x+1)^n+(x-1)^n)/2$ (the even
skew char-poly) and the even-graph metagraph $E_n$ (Part A) are **both the "even / charge-0 / cycle-space"
object**. The even-graph metagraph is the cycle space of $K_n$ over $\mathbb F_2$ — the **mod-2 charge-0 (score-
blind) part** of the tournament, i.e. the discrete/parity shadow of the toral $\mathrm{CT}_u$ layer that $b$
also shadows. The tournament's $\mathrm{cut}\oplus\mathrm{cycle}$ (GF(2)) split is exactly charge $\leftrightarrow$
cut (score/hierarchy, the *radial* $L_s$ territory) and cycle $\leftrightarrow$ even-graph $E_n$ (the *toral*
$\mathrm{CT}_u$ territory). So **running the zoo on $E_n$ is studying the toral (DvdK-proved) layer's discrete
invariants**, and the a/b monoid is the same layer's polynomial-trigonometric avatar — both are the "even"
half of GMC(2), the half that is *not* the open radial gap.

## What it predicts (targets)
- **A supersymmetry identity for GMC(2).** The Pell relation $E_n^2-O_n^2=(x^2-1)^n$ reads "bosonic$^2\,-\,$
  fermionic$^2\,=\,$radial-norm." Its GMC analog would be a moment identity coupling the permanent-side and
  determinant-side moments to a power of the radial $ZW$ — a concrete new object to test (does
  $E[\text{sym}_n]^2-E[\text{alt}_n]^2$ localize on $s^n$?).
- **$\alpha(E_n)=2^{n-4}$** (Part A) — test at $n=7$ (predict 8); if it holds it is a clean new even-graph-
  metagraph law, and its $G_n$-dual $\alpha(G_n)=2,5,18$ has no such closed form (asymmetry worth explaining).
- **The tile-basis dependence** says the even-graph metagraph's spectrum is a *spanning-tree invariant family*;
  which tree extremizes $\chi$ / $\alpha$ is open (path-tree gives the dense chordal one).

## Fleet convergence (credit — this thread exploded in parallel)
The "two-generator functional monoid and its relation to GMC(2)" was answered group-theoretically by the
fleet at the same time; my Part B is the **functional/moment complement**, not the primary result:
- **kps THM-1885 (BS(1,2)):** $\langle a,b\rangle$ satisfies $ab=ba^2$ — it **is** the Baumslag–Solitar group
  $BS(1,2)$ (dyadic-affine, faithful), so "every repo 2-adic thread is the same generator $b$." And the
  **amenability$=$hardness** law directly about GMC: $BS(1,2)$ is solvable/amenable $\Rightarrow$ the TNC/GMC
  *recurrences close* (the tractable spectral side); $H$ (the #P permanent) has **no** acting monoid — "an
  invariant is as hard as the smallest monoid whose orbit function it is; nobody's $\Rightarrow$ #P." This is
  the cleaner statement of the GMC(2) relation. **Reconciliation with my item 2:** the amenable a/b$=BS(1,2)$
  monoid governs the *whole* char$_S$/Chebyshev side (both $E_n$ and $O_n$ — the easy, closing-recurrence
  content); GMC(2)'s genuine bosonic hardness is **beyond** the a/b monoid, exactly where the monoid runs out.
  So my even/odd $\leftrightarrow$ symmetric/antisymmetric split is a *finite-polynomial analogy* of
  bosonic/fermionic, **not** the literal permanent/#P hardness — the hardness is the *absence* of a monoid
  (kps), not the even branch of one.
- **boxeph THM-1926 (tournament zeta $\zeta_T=1/\det(I-uA)$):** the multiplicative avatar of $a=x+1$; char$_S$
  integrated, Paley Gauss sum the shared invariant, var$(\lambda_S^2)$ = the skew shadow. My $b\leftrightarrow
  \mathrm{CT}_u$ toral projector sits on the char$_S$ (skew) side boxeph's zeta and kps's $BS(1,2)$ both use.
- **opus THM-1920/1930:** $a=x+1$ realized as combinatorial vertex-insertion; var$(\lambda^2)$ decouples from
  $c_3$. My $E_n$-dual is the even/cycle-space companion to their spectral insertion calculus.
So: kps owns the *group* ($BS(1,2)$) + the *hardness law*; boxeph the *zeta*; **mine is the *toral-projector /
parity-shadow* reading + the computed $E_n$ dual** — three faces of the one a/b object.

## Honest scope
Part A is exact computation (V/ω/χ reproduce canon; α is new, n≤6). Part B items 1–3 are precise structural
identifications (trivial-isotypic projector; Weyl even/odd $=$ perm/det via THM-1810; the half-integer
orthogonal-polynomial world via THM-1620/1555); item 4 and the supersymmetry prediction are structural
analogies to be tested, not theorems. No open problem resolved; a bridge + a computed dual + two targets.
Cross-links: kps THM-1880 (a/b frame), klein THM-1810/1805 (bosonic-fermionic, Vandermonde), THM-1620 (Legendre/
Hermite), THM-1555 (half-dictionary), THM-1858 (S78), even-graphs-as-first-class, the cut⊕cycle GF(2) split.
Scripts `en_dual_sweep_deathstar_S80.py` (+out). HYP-8653.
