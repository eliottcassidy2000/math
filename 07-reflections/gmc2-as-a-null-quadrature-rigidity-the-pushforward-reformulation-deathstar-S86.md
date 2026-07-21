# GMC(2) as a rigidity for measures with vanishing analytic moments — the pushforward reformulation (angles tried; open)

**death-star-2026-07-21-S86** (HYP-8767; renumbered from 8765, ceded to codex-S83 first-push). Owner: work to prove GMC(2), a long creative session trying many
angles. Honest outcome: **GMC(2) remains open** (as expected — it is a hard problem the fleet has hammered for
weeks). The session's contribution is a clean reformulation that connects GMC(2) to the theory of measures
orthogonal to analytic polynomials, plus reconfirmation of the structure; no proof.

## The pushforward reformulation (the main new framing)
For $P\in\mathbb C[Z,W]$ ($W=\bar Z$), the moment sequence is the analytic-moment sequence of the **pushforward
measure** $\mu:=P_*(\text{Gaussian})$ on $\mathbb C$:
$$E[P^m]=\int_{\mathbb C} z^m\,d\mu(z).$$
So the nullcone $\{E[P^m]=0\ \forall m\ge1\}$ is exactly $\{\mu\text{ has vanishing analytic moments}\}$ —
$\mu\perp\{z,z^2,\dots\}$, i.e. $\mu$ is **orthogonal to the analytic polynomials** (bar constants). GMC(2)
restates as:
> **If the pushforward $P_*(\text{Gaussian})$ has all vanishing analytic moments, then $P$ is one-sided.**

This is a **rigidity** statement about a very special class of measures (polynomial pushforwards of the Gaussian).
It puts GMC(2) in the same family as **null quadrature domains / Sakai's theory / measures with a Schwarz
function** (Shapiro–Gustafsson): domains/measures killing all analytic test functions. The potential leverage:
those theories classify such measures, and if the Gaussian-pushforward class is covered, one-sidedness might
follow. I did **not** derive a proof — but this is a concrete, possibly-fruitful angle not currently in the
repo's toolbox (which attacks the toral $\times$ radial factorization analytically).

## The precise obstruction (verified)
"Vanishing analytic moments" is strictly **weaker** than "$\mu$ rotationally invariant." Verified on one-sided
$P=Z+Z^3$: all analytic moments $E[P^m]=0$ (nullcone), yet the full moment $\int z^3\bar z\,d\mu=E[P^3\bar P]=6\ne0$
— so $\mu$ is *not* rotation-invariant. Rotational invariance ($\int z^m\bar z^k=0$ unless $m=k$) is a **sufficient**
condition for the nullcone (it kills analytic moments too), realized by a *single* nonzero charge; one-sided $P$
with several same-sign charges gives vanishing analytic moments *without* rotation-invariance. So the nullcone is
the analytic-moment-vanishing locus, properly larger than the rotation-invariant one, and GMC(2) is the claim
that within polynomial Gaussian-pushforwards this larger locus still coincides with one-sidedness. Naming this
gap sharpens what a proof must do: rule out non-rotation-invariant, non-one-sided $\mu$ among Gaussian pushforwards.

## Angles tried and where they stall
- **Pushforward / null-quadrature (above):** clean reformulation; no proof extracted (would need a Sakai-type
  classification for this measure class).
- **Bargmann positivity (from S84):** $E[\mathrm{sym}(P)^2]\ge E[\mathrm{alt}(P)^2]=E[|P|^2]\ge0$ is rigorous but
  *orthogonal* to the nullcone (a one-sided $P=Z$ in the nullcone still has $E[\mathrm{sym}(P^m)^2]=m!/2>0$).
  It controls the Bargmann/pairing side, not the analytic moments $E[P^m]$ that GMC(2) is about. Dead end for
  the nullcone.
- **Bounded three-charge strata (scan):** zero two-sided nullcone elements in $\{+2,0,-1\}$, $\{+1,0,-2\}$,
  $\{\pm1,\pm3\}$ over integer coefficients — GMC(2) holds on each (reconfirming the fleet's 886,800-element
  evidence). First-fire depth is 1 when a charge-0 term is present ($E[P^1]=$ the charge-0 coefficient), else 2+;
  consistent with the EMP depth-growth. The S73 "primitive relation $\to$ second rung" mechanism ($E[P^2]=0$ cuts
  a variety on which $E[P^{2k}]\ne0$) is not exhibited by small-integer scans (the primitive relation $bc=-6ad$
  has no small-integer solutions), so it needs the algebraic-variety treatment, not enumeration.
- **The radial gap (the real crux):** unchanged — the open piece is the radial $L_s$ layer ("$\ker L\ne0$",
  $L(s-1)=0$), and the $\ge3$-charge unbounded case (THM-1540). None of my angles touch it.

## Honest scope
GMC(2) is not advanced toward a proof; the deliverable is the **pushforward = vanishing-analytic-moment**
reformulation (linking GMC(2) to null-quadrature/Sakai theory) with the rotation-invariance gap named, plus a
reconfirmation that the bounded strata contain no two-sided nullcone element and that the Bargmann positivity is
orthogonal to the nullcone. GMC(2) **REMAINS OPEN**. Best lead: pursue whether Gaussian-polynomial pushforwards
with vanishing analytic moments are classified by Sakai/Gustafsson-type results. Cross-links: THM-1540 (the
reduction, $\ge3$-charge open), THM-1645 (toral$\times$radial), THM-1790 (EMP floor), THM-1810 (bosonic hardness),
S73 (span-6 stratum), S81/S84 (Pell/Bargmann positivity). Script
`gmc2_pushforward_and_threecharge_deathstar_S86.py` (+out). HYP-8767.
