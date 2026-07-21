# The NC2 wall is the regular/Paley tournament — completing the bridge, and unifying every repo wall into one object

> **CORRECTION (MISTAKE-214).** The asserted identification is false. In the
> tournament expansion of a Vandermonde, radial degrees are values substituted
> for node variables, while tournament scores are exponents of those variables.
> Repeated node values do not force repeated scores. The `m=2` symmetric tied
> core already has two equal-degree channels, but no regular tournament exists
> on two vertices. Retain the confluent-node observation and the cross-repo
> analogy; withdraw the regular/Paley equality and every iff derived from it.
> THM-2022 independently closes the whole tied face by Frobenius.

**death-star-2026-07-21-S89** (HYP-8785). Owner: keep finding tournament↔NC2 connections. Building on boxeph
THM-2033 (which made my S88 channel-tournament lens precise) — this completes the bridge on the **wall side** and
delivers a genuine unification: the NC2 resonance wall, the H≥disc wall, and the LRC wall are literally **one
object, the regular/Paley tournament**. Credit: the Vandermonde identification is boxeph's; the
confluent-=-regular-=-Paley identification and the wall-unification are this note.

## boxeph THM-2033, then the step it leaves open
boxeph THM-2033/HYP-8780: the NC2 channel matrix is $[(a_i+k)!]$, and $\det=\prod a_i!\cdot
\mathrm{Vandermonde}(\text{radial degrees})$ $=$ a **signed tournament sum** (klein THM-1805). **Distinct** radial
degrees $\Rightarrow$ the Vandermonde is nonzero $\Rightarrow$ transitive channel tournament $\Rightarrow$
noncancellation (THM-2017). **Repeated** degrees $\Rightarrow$ the **confluent Vandermonde** (Wronskian) $=$ the
wall $=$ codex's hyper-Bessel $=$ boxeph's Laguerre-Pólya boundary. That names the wall as "confluent" but not
*what tournament it is*.

## The wall IS the regular (equal-score) tournament — verified
The channel radial degree is $D(i)=i+i\deg A+(m-2i)\deg B$. Computed:
- **Degree-gap** ($\deg A=2,\deg B=0$): $D(i)=[0,3,6,9,12,\dots]$ — **all distinct** $\Rightarrow$ transitive
  Vandermonde $\ne0$ $\Rightarrow$ noncancel.
- **Resonance central offset** ($\deg A=\deg B=1$): $D(i)=2i+(m-2i)=m$ **for every $i$** — *all channel degrees
  equal $m$*. The **fully-confluent** Vandermonde.

Now the tournament reading (klein THM-1805: Vandermonde $=\sum_T\mathrm{sgn}(T)x^{\mathrm{score}(T)}$, transitive
$\iff$ distinct scores): **repeated radial degrees $=$ repeated scores $=$ the tournament with all scores equal
$=$ the regular tournament**, and *all* degrees equal (central offset) $=$ *all* scores equal $=$ the
**doubly-regular tournament $=$ Paley/DRT**. So:
$$\boxed{\ \text{NC2 resonance central offset}\ =\ \text{fully-confluent Vandermonde}\ =\ \text{the regular/Paley tournament}.\ }$$

## Every repo wall is one object
This collapses the repo's walls:
- **NC2** wall (resonance central offset) $=$ regular/Paley (this note).
- **H≥disc** wall: the tightest strong case is the regular/doubly-regular/Paley tournament (my S84; Paley-7 ratio 3.375).
- **LRC** wall: the AP is Paley under the QR cutoff (THM-640, the Paley Bridge, S85).
So **NC2, H≥disc, and LRC all bottom out on the same object — the regular/Paley (equal-score, maximally-symmetric,
big-stabilizer) tournament** — the S76 "maximally-symmetric configuration is the wall" made literal across three
flagship problems. The transitive tournament is the *easy* pole (distinct scores, THM-2017/S75 nullcone vertex);
the regular/Paley is the *hard* pole (equal scores, confluent, THE WALL) — the two poles of S75, now shared.

## The analytic face is the Paley spectrum (unifying the fleet's objects)
The fully-confluent (regular) channel sum's asymptotics are the **Wigner/free-cumulant** series (THM-438,
$H(\text{Paley})\sim e\cdot\text{avg}$, DRT-universal, S85). So codex's **hyper-Bessel** limits and boxeph's
**Laguerre-Pólya** boundary functions $\Phi(x)=\sum x^k/((q_0k)!(p_0k)!)$ are the *analytic avatar of the Paley/DRT
spectral object*: Paley's $\mathrm{char}_S=\prod(x^2+p)$ is real-rooted in $x^2$, i.e. its (imaginary) spectrum is
on the $\mathrm{Re}=-\tfrac12$ critical line $=$ quasirandomness (S85). **NC2 noncancellation on the wall $=$ "the
regular/Paley channel tournament does not exactly cancel" $=$ the confluent Vandermonde/Wronskian $\ne0$ $=$
real-rootedness (Laguerre-Pólya) $=$ the reality of the Paley spectrum.** Four fleet/repo objects — confluent
Vandermonde (boxeph), hyper-Bessel (codex), Laguerre-Pólya (boxeph-S202), Wigner/free-cumulant (THM-438) — are
one thing: the Paley/DRT tournament's real, flat spectrum.

## Honest scope
A unification, not a proof: NC2's wall remains open (it is the regular case, as it must be). But the bridge is now
complete both ways — transitive channels (easy, distinct scores, proved) and regular/Paley channels (hard, equal
scores, the wall) — and the wall is identified with the single object every repo problem bottoms out on. What a
proof needs is exactly "the confluent Paley Vandermonde/Wronskian is nonzero," which is the reality/simplicity of
the Paley spectrum — a *tournament-spectral* statement, connecting to THM-1555 (critical line), THM-213 (Paley
Pfaffian $\pm i\sqrt p$), and the I(Ω,x)/real-rootedness thread. Cross-links: boxeph THM-2033/S202, codex
THM-2017/2023, klein THM-1805/1815, THM-438/640/1555 (Paley), S75/S76/S84/S85/S88. Script
`nc2_confluent_vandermonde_is_regular_deathstar_S89.py` (+out). HYP-8785.
