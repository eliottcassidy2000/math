# Rédei through the sign-reversing engine: the descent it DOES give, the determinant collapse it does NOT, and the hypotheses it spins off

*mac-mini-2026-07-21-S159 (follow-up). The owner: do a fresh proof of the founding theorem
(Rédei — every tournament has an ODD number of Hamiltonian paths) with the sign-reversing
tournament-involution engine, and see how the same idea creates new proofs and hypotheses. This
is the honest report: the engine gives Rédei cleanly by DESCENT (verified), it provably does NOT
give it by a determinant COLLAPSE (verified — a real boundary of the technique), and it spins off
one clean new object and one instructive refuted hypothesis.*

---

## What the engine is, and the crucial distinction

The S159 engine is: a **signed sum over tournaments with a sign-reversing involution** (reverse a
3-cycle) that cancels the cyclic/intransitive terms and leaves the transitive fixed points. It
proves things that are **sums over ALL tournaments** (Vandermonde → transitive; Burnside
even-cycle → 0). **Rédei is different in kind:** it is a sum over the Hamiltonian paths of **one**
tournament, $h(T)=\sum_{\pi}\prod_i A_{\pi_i\pi_{i+1}}$ — an *unsigned, permanental* count. So the
engine cannot be dropped in verbatim; the interesting content is exactly *how much* of it transfers.

## 1. What the engine DOES give: Rédei by descent (verified)

The engine's slogan "cyclic content cancels, the transitive core survives" becomes, for one
tournament, a **parity descent**:

> **$h(T)\equiv h(T-v)\pmod 2$** for every vertex $v$ — verified for **all** tournaments and all
> $v$ at $n=4,5$ (256 and 5120 checks, zero failures).

Mechanism (the engine's shadow): build the Ham paths of $T$ by inserting $v$ into linear orders of
$T-v$. Deleting $v$ and splicing its neighbors yields a Ham path of $T-v$ **unless** $v$ sat inside
a 3-cycle $a\to v\to b\to a$ (the one *intransitive* configuration), which leaves a single defect.
Those cyclic insertions are the terms that cancel in pairs; the transitive skeleton descends. By
induction $h(T)\equiv h(\text{single vertex})=1$: **odd**. This is Rédei's classical induction
*read through the engine* — the same "intransitive cancels, transitive survives," now as a descent
rather than a collapse. (Honest: the theorem and induction are classical; the engine reading and
the exhaustive verification are what this file contributes.)

## 2. What the engine does NOT give: no determinant collapse (verified — the real finding)

For the sum-over-all-tournaments results the involution *collapses* the signed sum to a determinant
(the Vandermonde). One would hope Rédei's parity equals some $\det M(T)\bmod 2$. **It does not** —
tested exhaustively:

| candidate mod 2 | $n=3$ | $n=4$ | $n=5$ | always odd? |
|---|---|---|---|---|
| $\det(A)$ | 2/8 | 24/64 | 400/1024 | no |
| $\det(A+I)$ | 6/8 | 24/64 | 304/1024 | no |
| $\det(A-A^{\mathsf T})$ | 0/8 | **64/64** | 0/1024 | no |
| $\det(A+A^{\mathsf T}+I)$ | 0 | 0 | 0 | no |

None is universally odd, so **none can equal $h(T)\bmod2\equiv1$.** The near-miss is illuminating:
over $\mathbb F_2$, $A-A^{\mathsf T}=A+A^{\mathsf T}=J-I$ is **tournament-independent**, so
$$\det(A-A^{\mathsf T})\bmod2=(n-1)\bmod2=[\,n\text{ even}\,],$$
which is exactly **THM-1440's forced-zero parity = the blue parity** (odd $n$ → the skew-adjacency
is singular). It is a clean odd/even invariant — but *because* it is tournament-independent it
**cannot** be $h(T)$. So:

> **Rédei's parity is genuinely permanental**, not linear-algebraic. The engine reaches it only by
> descent, never by a determinant collapse — and the one determinant that "wants" to be it
> (the skew-adjacency) is the blue-parity invariant instead. This is a precise boundary of the
> S159 technique: *collapse* is for signed sums over all structures; *descent* is what survives
> for permanental counts over one structure.

## 3. Hypotheses the idea spins off

**(a) NEW object — the signed Ham-path count.** $R(T)=\sum_{\pi\ \text{Ham path}}\operatorname{sgn}(\pi)$
(the engine-aligned *signed* Rédei) is **always odd** (verified $n\le5$), so "R(T) odd" is a signed
restatement of Rédei. Its distribution is symmetric and **gapped**: at $n=5$,
$|R|\in\{1,3,5,7,11,15\}$ — **9 and 13 never occur.** Characterizing the achievable $|R(T)|$ (and
the gaps) is a clean new question; $R(T)$ is not a plain determinant, so it is genuinely new
structure, not a repackaging.

**(b) REFUTED — a mod-4 Rédei by cycle count.** The natural engine guess "the deviation of $h$ from
1 is the cyclic content" predicts
$$h(T)\equiv 1+2\,c_3(T)\pmod4,\qquad c_3=\#\text{3-cycles}.$$
It holds at $n=3,4$ (all tournaments) and **fails at $n=5$** (624/1024) — the classic $n=5$ wall.
No simple odd-strongly-connected-subset count mod 2 rescues it either (best, $c_3+c_5$: 880/1024 at
$n=5$, 25280/32768 at $n=6$). **Verdict:** the mod-4 digit of $h(T)$ is *higher-order* — not a
linear function of small-cycle counts. (Dead-end recorded so no one re-derives it: the deviation
from Rédei's $1\bmod2$ is not carried by 3-cycles once $n\ge5$, exactly because strongly-connected
blocks larger than a triangle contribute.)

## The transferable lesson (how the idea makes new proofs)

The engine is a **diagnostic**: attack a repo quantity with the sign-reversing involution when it
is a *signed sum over all structures* (then it collapses to the transitive core — Vandermonde,
Burnside, the GMC discriminant); expect only a *descent/parity recursion* when it is a *permanental
count over one structure* (Rédei); and expect **no help at all** when the weights are not clean
$\pm1$ (LRC's sinc covering — the S157 barrier). Rédei sits in the middle box, and this pins why:
its parity is permanental, its one candidate determinant is the blue-parity invariant, and its
finer (mod-4) structure is genuinely nonlinear.

## Honest scope

- $h(T)\equiv h(T-v)\pmod2$ and "no determinant candidate is universally odd" are **verified
  exhaustively** ($n\le5$; the descent also $n=5$). The Rédei theorem/induction is **classical**;
  the engine reading, the determinant-limit finding, and the $\det(A-A^{\mathsf T})=$ blue-parity
  cross-link are the contributions.
- $R(T)$ odd is verified $n\le5$; the gap structure ($9,13$ absent at $n=5$) is observed, not
  explained.
- The mod-4 hypothesis is **REFUTED** at $n=5$ and logged as such.

---

*Cross-links: the S159 engine reflection (the-sign-reversing-tournament-involution-as-a-repo-wide-engine),
THM-1805/1815 (Vandermonde/discriminant = transitive survive), THM-1440 (skew-adjacency
forced-zero = blue parity), THM-1840-C (blue parity = $\mathbb Z/2$ cyclotomic character), and the
project's Rédei/OCF core. Artifacts:
`04-computation/redei_signed_involution_freshproof_macmini_S159.py`,
`04-computation/redei_mod4_refinement_macmini_S159.py` (+outs).*
