# "$X\cong\mathbb A^3$" is equivalent to disproving JC(3): the symmetric-product étale witness

**death-star-2026-07-20-S74** (HYP-8600). Owner (consultation): is there an easy non-computational
way to see $X\cong\mathbb A^3$ for a specific construction, and achieve complete structural
understanding? The construction turns out to be a **candidate Jacobian-Conjecture(3) counterexample**,
and the "easy proof" the owner hoped for is exactly the thing that would disprove JC(3). Recording
because it sits squarely on the repo's GMC⇒JC thread.

## The construction
$\pi\colon \mathbb P^1\times\operatorname{Sym}^2\mathbb P^1\to\operatorname{Sym}^3\mathbb P^1$,
$(p,\{q,r\})\mapsto\{p,q,r\}$ (= multiplication of binary forms, degree 3). $R$ = ramification
(a marked root becomes double). $H\subset\operatorname{Sym}^3\mathbb P^1=\mathbb P^3$ a hyperplane
**tangent but not osculating** to the twisted cubic $C$ (triple points). $Y=\operatorname{Sym}^3\setminus
H\cong\mathbb A^3$ (obvious). $X=(\mathbb P^1\times\operatorname{Sym}^2)\setminus(R\cup\pi^{-1}H)$,
claimed $\cong\mathbb A^3$; $\pi|_X\colon X\to Y$ the "counterexample."

## The key fact (reorganizes everything)
$\pi|_X$ is **étale** (deleted the ramification $R$), generically **3:1** (non-injective), and
**non-surjective**: over a triple point $\{p^3\}\in Y$ the only root $p$ is not *simple*, so
$(p,\{p,p\})\in R$ is deleted — no preimage. Hence $\operatorname{im}(\pi|_X)=Y\setminus(C\cap Y)$,
missing $C\setminus H\cong\mathbb G_m$ (a **codim-2** curve). Therefore:

> **If $X\cong\mathbb A^3$, then $\pi|_X$ is a polynomial endomorphism of $\mathbb A^3$ with
> $\det J\in\mathbb C^\*$ (étale) that is not an automorphism (degree 3) — a counterexample to the
> Jacobian Conjecture in dimension 3.**

So "$X\cong\mathbb A^3$" is **not** a cheap lemma: it is (essentially) equivalent to disproving JC(3).
The naive Serre-vanishing argument ("$\mathbb A^1$-bundle over affine $\mathbb A^2$ is trivial") silently
assumes the fibration is a **torsor** (Zariski-locally trivial); JC $\Rightarrow$ that must fail. The
$\mathbb A^1$-**fibration** is real; the $\mathbb A^1$-**bundle** upgrade is the whole ballgame.

## What $X$ almost surely is
An **exotic contractible affine 3-fold** ("fake $\mathbb A^3$"), à la Koras–Russell / Zariski
cancellation — $\mathbb A^1$-fibered over $\mathbb A^2$ by a *non-trivial* (Danielewski-type) fibration.
Model case, degree 1: $\mathbb A^2\setminus\{0\}\hookrightarrow\mathbb A^2$ is étale, misses a codim-2
point, and its source is not $\mathbb A^2$. $X$ is the affine degree-3 analogue.

## Complete structural scaffolding (representation theory of binary forms)
- $\pi$ = Segre 3-fold $\mathbb P^1\times\mathbb P^2\subset\mathbb P^5$ projected from the
  **Clebsch–Gordan kernel line** $\mathbb P(V_1)$ (since $V_1\otimes V_2\cong V_3\oplus V_1$).
- Hyperplanes in $\mathbb P(V_3)$ **are binary cubics $h$** via apolarity ($V_3\cong V_3^\*$), and
  **contact order of $H$ with $C$ at $[\ell_p^3]$ = multiplicity of $p$ as a root of $h$**. So:
  secant $\leftrightarrow h$ type $(1,1,1)$; **tangent-not-osculating $\leftrightarrow h$ type $(2,1)$**;
  osculating $\leftrightarrow h=\ell^3$. Residual symmetry: $\operatorname{Stab}_{PGL_2}(h)=\mathbb G_m$.
- Clean model: $X=\{(p,f):p$ a **simple** root of $f$, $f\not\perp h\}$; $\pi$ = forget the marked root.
- Tangent-not-osculating is Goldilocks: the second polar $\ell_p(\partial)^2h$ (which controls the
  natural vertical $\mathbb G_a$) vanishes at exactly **one** point $p=p_0$ for type $(2,1)$ (double
  degeneration if osculating; stray $\mathbb G_m$ if secant), keeping $X$ smooth & contractible.

## Bearing on the repo's GMC⇒JC thread
GMC realizes the Jacobian obstruction as a Mathieu–Zhao (kernel-of-a-functional) condition; its
failure at $n\ge3$ is a **codimension phenomenon** ("one extra degree of freedom breaks affine-space
rigidity") — the *same shape* as tangent-vs-osculating and the codim-2 missed curve here. Concrete
follow-up to *decide* fake-vs-real $\mathbb A^3$: compute the **Makar–Limanov invariant** or $H^\*(X)$
for the type-$(2,1)$ choice; a nontrivial ML invariant would prove $X\not\cong\mathbb A^3$ (and dissolve
the JC scare); a trivial one keeps the tension live.

## Honest status
A consultation write-up, not a repo theorem: the structural picture (apolarity/contact, Segre-CG
projection, universal-simple-root) is solid; "$X\cong\mathbb A^3$" is the open crux, equivalent to a
JC(3) disproof, so the responsible prior is $X$ = exotic $\mathbb A^3$. If an independent proof
$X\cong\mathbb A^3$ exists, the two things to stress-test are (i) $R$ = the *entire* ramification and
(ii) the source iso being algebraic — any leak that saves JC lives in (i)/local-triviality.

## Cross-links
GMC/JC thread (GMC(n)⇒JC(n)), THM-1300 (Alpoge JC counterexample), Koras–Russell / Zariski cancellation,
Makar-Limanov invariant, apolarity of binary forms; classical: twisted cubic, its tangent developable
(the discriminant), Clebsch–Gordan. HYP-8600.
