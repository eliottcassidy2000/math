# The NC2 tied-core weights are the central trinomial — the free-probability bridge, completed

> **CORRECTION (MISTAKE-214).** The central-trinomial coefficient identity is
> correct. The regular/Paley identification inherited from S89 and the claimed
> equivalences `NC2 noncancellation <=> Laguerre--Polya <=> Paley critical-line
> spectrum` are unsupported and withdrawn. Equal Vandermonde nodes are not
> equal tournament-score exponents. Keep the free-probability language as an
> analogy, not a reduction. THM-2022 gives the actual NC2 certificate `Q^p`.

**death-star-2026-07-21-S90** (HYP-8790). Owner: keep finding tournament↔NC2 connections. This completes the
free-probability/semicircle bridge that boxeph-S203 (THM-2033) flagged as "next" and that my S88 free-prob lens
predicted — with a concrete identification.

## The identification
At the NC2 resonance central offset (the fully-confluent Vandermonde = the fully-regular/Paley tournament, my
S89), every channel has the same radial degree, so the channel *weights* are the pure multinomials
$\binom{m}{i,i,m-2i}$, and their sum is
$$W(m)=\sum_{i}\frac{m!}{i!\,i!\,(m-2i)!}\;=\;1,3,7,19,51,141,393,1107,3139,8953,\dots\;=\;\textbf{A002426, the central trinomial coefficients},$$
$W(m)=[x^0](1+x+x^{-1})^m$ (verified). The ratio $W(m)/W(m-1)\to3$. These are the **moments of the free
convolution / spectral measure of three atoms** — a Wigner/free-probability object, growth $\sim 3^m/\sqrt{\pi m}$.

## Why it closes the bridge
This is the concrete form of the S88 free-probability lens and the THM-438 recovery: THM-438 proved the Paley/DRT
**cluster integrals are Catalan = free cumulants of the two-point spectrum $\tfrac12(\delta_a+\delta_{-a})$**, giving
$H(\text{Paley})\sim e\cdot\text{avg}$. Here the NC2 **tied-core moment weights are the central trinomial**, the
free-moment sibling. So the whole "wall" object — regular/Paley tournament (S89) — carries a **free-probability
moment/cumulant structure on both its H-count (THM-438) and its NC2-channel weights (this note)**: one
free-probabilistic law, two appearances.

The noncancellation reading becomes sharp: NC2 fails on the wall iff the central-trinomial-weighted signed
channel sum vanishes for all $m$ — i.e. iff the associated free-cumulant/generating series has a real positive
zero. That is exactly a **Laguerre-Pólya failure** (boxeph-S202: the boundary $\Phi(x)=\sum x^k/((q_0k)!(p_0k)!)$
is Laguerre-Pólya, all zeros real-negative), which is the **reality of the Paley spectrum** (char$_S=\prod(x^2+p)$
on $\mathrm{Re}=-\tfrac12$, THM-1555/213). So: **NC2 noncancellation on the wall $\iff$ the free-cumulant series
has no real positive zero $\iff$ the Paley spectrum stays on the critical line** — three faces (combinatorial /
analytic / spectral) of one tournament fact.

## The completed tournament↔NC2 picture (S88–S90 + boxeph THM-2033)
- **S88 / boxeph THM-2033:** NC2's channels form a tournament; noncancellation = its transitivity =
  Vandermonde(radial degrees) = signed tournament sum (klein THM-1805); det/fermionic rigidity lent to the
  bosonic permanent (THM-1815).
- **S89:** the wall (confluent Vandermonde) = the regular/Paley tournament (equal scores); NC2 wall = H≥disc wall
  = LRC wall = one object.
- **S90 (this):** the wall's channel weights = central trinomial = a free-probability moment (THM-438 sibling);
  noncancellation $\iff$ Laguerre-Pólya $\iff$ Paley spectrum on the critical line.
Transitive channels = easy (distinct scores, S75 nullcone vertex); regular/Paley channels = the wall (equal
scores, free-probability semicircle core, S75 critical-line pole). The two S75 poles govern NC2 exactly as they
govern H≥disc and LRC.

## Honest scope
A concrete identification + synthesis, not a proof (NC2's wall stays open, being the regular case). Value: the
free-prob bridge is now explicit (central trinomial A002426 = the tied-core weights), and the three faces
(combinatorial channel sum / analytic Laguerre-Pólya / spectral Paley reality) are named as one target. Cross-links:
boxeph THM-2033/S202, codex THM-2017/2023, S88/S89, THM-438 (Wigner/Catalan), THM-1555/213 (Paley spectrum),
THM-1815/1805 (det bridge), S85 (Paley recovery). A002426 (OEIS central trinomial). Script inline (central-trinomial
verification, in the S89/S90 results). HYP-8790.
