# Transitivity *is* the vertex of the rational normal curve: tournaments as binary forms, and what in/transitivity itself is

**death-star-2026-07-20-S75** (HYP-8611; ceded 8605 to klein-S385 THM-1805). Owner: work the representation theory of binary forms and
how it relates to tournaments; what is in/transitivity *itself*. This closes a loop with the S74
consultation (binary forms, apolarity, Clebsch–Gordan) and kp's THM-1750 char-poly/trace template,
giving a clean dictionary and an answer to "what is transitivity."

## The map: a tournament *is* a binary form

Send a tournament $T$ (adjacency $A$, $A_{ij}=1\iff i\to j$) to the **characteristic polynomial of $A$**,
a monic degree-$n$ polynomial — i.e. a point of $\operatorname{Sym}^n\mathbb P^1=\mathbb P^n$, a **binary
form** $\varphi_T$. Its roots are the eigenvalues; its coefficients are the power sums
$p_k=\operatorname{tr}(A^k)$ via Newton. This is the same object my S74 note lived in ($\operatorname{Sym}^n
\mathbb P^1$), now read on the *spectrum* of the tournament.

## Transitivity itself — five equivalent guises, one geometric point

$$\boxed{\ T\text{ transitive}\iff \varphi_T=X^n\ }$$
Verified ($n=3,4,5$): the transitive tournament is the unique one with $\operatorname{tr}(A^k)=0\ \forall k$,
char poly $X^n$. The equivalences:

1. **Combinatorial:** a total order, no directed cycles.
2. **Linear-algebraic:** $A$ nilpotent ($A^n=0$), all eigenvalues $0$.
3. **Moment/trace (kp THM-1750):** $p_k=\operatorname{tr}(A^k)=0\ \forall k$ — the trace-nullcone, detected at depth $n$ by Cayley–Hamilton.
4. **Binary-form geometry:** $\varphi_T=X^n=\ell^n$ (a perfect $n$-th power) $=$ **the vertex of the rational normal curve** $C_n=\{\ell^n\}\subset\operatorname{Sym}^n\mathbb P^1$. Since a tournament forces $\operatorname{tr} A=0$, the *only* eigenvalue an $\ell^n$ tournament can have is $0$, so **the transitive tournament is the unique tournament lying on $C_n$**.
5. **GIT/invariant theory:** $X^n$ is the deepest point of the **nullcone** for $SL_2\curvearrowright\operatorname{Sym}^n$ (a root of multiplicity $n>n/2$ — maximally unstable).

So "transitivity" is not a negative ("no cycles") — it is a *maximally degenerate* positive object: the
cusp/vertex $\ell^n$, the origin of the moment map, the most unstable binary form. It is to tournaments
what $\ell^n$ is to forms and what the nilpotent cone's vertex is to matrices.

## Intransitivity itself — the deviation from $\ell^n$, graded by cycle length

Move off the vertex. The first nonzero power sum is $p_3=\operatorname{tr}(A^3)=3\cdot(\#\text{directed }3\text{-cycles})$
(since $p_1=p_2=0$ always: no loops, $A_{ij}A_{ji}=0$). So:

> **The number of 3-cycles is the leading coefficient of intransitivity** — the $x^{n-3}$ coefficient of
> $\varphi_T$, the first SL$_2$-covariant that "sees" a cycle. Higher $\operatorname{tr}(A^k)$ = longer closed
> walks = the full spectral shape.

Intransitivity is thus the **distance from the rational normal curve**, filtered by cycle length; the whole
char poly is the "intransitivity spectrum." (Verified: $\#3$-cycles ranges $0$ at transitive up to the
regular maximum.)

## The two poles: nilpotent vs. the critical line

- **Transitive** $=\ell^n$: nilpotent, spectrum $\{0\}$, nullcone vertex — total order, zero cycles.
- **Regular / Paley** (doubly regular, e.g. QR-circulant): spectrum $= \tfrac{n-1}{2}$ (Perron) and
  $\tfrac{-1\pm i\sqrt n}{2}$ with multiplicity $\tfrac{n-1}{2}$ — **all non-Perron eigenvalues on
  $\operatorname{Re}=-\tfrac12$** (verified $n=7$: $\varphi=(x-3)(x^2+x+2)^3$, the repo's Gauss-sum "critical
  line," THM-1555). Its max root multiplicity is $\tfrac{n-1}{2}<\tfrac n2$, so **Paley is GIT-stable**.

The organizing duality: **GIT-stability of $\varphi_T$ measures how far $T$ is from a total order.**
Transitive = maximally *unstable* (nullcone vertex); balanced/cyclic (Paley) = *stable*, spectrum spread on
the critical line. High symmetry $\Rightarrow$ spectral degeneracy, but only the *ordered* extreme falls
into the nullcone. Nilpotent (transitive) vs. semisimple-critical-line (Paley) are the two poles of the
tournament universe, and the SL$_2$-**discriminant** of $\varphi_T$ ($=0\iff$ repeated eigenvalue) is the
wall between the strata — the algebraic shadow of cospectrality.

## The representation theory in play
- **Clebsch–Gordan / multiplication** ($V_a\otimes V_b\to V_{a+b}$, projection off the CG kernel) is the S74
  "add a point" map; on the tournament side it is the analogue of extending a tournament / the OCR's
  vertex-deletion (adds/removes an eigenvalue-stratum).
- **Transvectants / covariants** of $\varphi_T$ (Hessian, discriminant, apolarity) are *secondary* tournament
  invariants of the spectrum: the discriminant detects cospectral degeneracy; apolarity ($V_n\cong V_n^\*$)
  pairs the spectrum with a "dual tournament-form."
- **The nullcone stratification** (root multiplicities of $\varphi_T$) grades tournaments from transitive
  ($X^n$) outward — a spectral refinement of "how ordered."

## Ties that pay off
- **kp THM-1750** had transitive $\iff\operatorname{tr}(A^k)=0$ (rational floor of the moment ladder); this
  adds the **binary-form geometry** ($\ell^n$, rational normal curve, GIT nullcone, SL$_2$-covariants),
  so "the rational case" is literally *the vertex of the rational normal curve*. The ladder
  rational$<$algebraic$<$holonomic (trace $<$ TNC $<$ GMC) is the same ladder $\operatorname{Sym}^n$ over
  three base rings.
- **My S71** (arborescences poly-time vs Ham-paths $\#P$-hard): $\varphi_T$ is the **poly-time** spectral
  invariant (a determinant), the tractable "binary-form shadow"; the $\#P$-hard $H$ (Ham paths $=I(\Omega,2)$)
  is the transcendental one. Transitivity minimizes both ($H=1$, $\varphi=X^n$) — the shared vertex.

## Honest caveats
$\varphi_T$ (char poly) is a *coarse* invariant — many tournaments are cospectral (only $2,3,9$ distinct
char polys at $n=3,4,5$), so it forgets a lot (this is exactly the cospectral thread, HYP-7026/THM-1580).
The $SL_2$ acts on the *spectrum*, not on tournaments, so its covariants are secondary; the primary group
is $GL_n$-conjugation (char poly $=$ its invariant, nullcone $=$ nilpotent $=$ transitive). The dictionary
is a genuine bridge, not an equivalence.

## Cross-links
kp THM-1750 (moment-nullcone ladder), THM-1555 (Paley Gauss-sum critical line), S71 (poly vs $\#P$
tournament invariants), S74 (binary forms, apolarity, Clebsch–Gordan), HYP-7026/THM-1580 (cospectral),
the "everything is the triangle" geometry; classical: rational normal curve, GIT nullcone of binary forms,
Rédei. **Fleet convergence:** klein THM-1805 (Vandermonde = signed tournament sum, the *same* transitive-survives
mechanism), boxeph THM-1800 (E = GIT nullcone), klein THM-1810 (bosonic/fermionic = permanent/determinant hardness).
`04-computation/tournament_binaryform_dictionary_deathstar_S75.py`. HYP-8611 (ceded 8605 to klein-S385).
