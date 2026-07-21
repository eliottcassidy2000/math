# Dividing by the common factorial (pA₀)! turns NC2 into one tournament-discriminant condition

**death-star-2026-07-21-S91** (HYP-8795). Owner: think about dividing E[P^m] by the full common factorial
(pA₀)!. Developed — it is the canonical normalization that peels off the trivial factorial growth and exposes
the tournament/Vandermonde discriminant as the *entire* NC2 condition, unifying the S88–S90 arc and boxeph
THM-2033 into a single clean statement.

## What (pA₀)! is, and what dividing by it does
For the general three-weight $P=Z^p a(s)+b(s)+\bar Z^q c(s)$, charge-0 forces $pj=qk$; the channel-$t$ term uses
$j$ copies of $Z^p a$, and its Laplace value carries the factorial $(pj)!$. The **dominant common factorial** is
$(p\,j_{\max})!$ — that is, $\boxed{(p A_0)!}$ with $A_0=j_{\max}$ the **maximum multiplicity of the top charge
atom** $Z^p a$ (the endpoint/source channel). Verified:
- $p{=}q{=}1$: dividing by $m!$ turns $E[P^m]$ (factorial-growing) into a **bounded rational sequence** — the
  factorial growth is entirely the common $(pA_0)!=m!$, and the residual is the tournament/free-prob core (S90).
- $p{=}2,q{=}1$: the dominant factorial is $(2\lfloor m/3\rfloor)!=(pA_0)!$, and $E[P^m]/(pA_0)!$ is O(1)-scaled
  per channel — the **endpoint channel $\to$ a nonzero constant** (codex's "one endpoint ratio one").

## Why it collapses NC2 to one condition
After dividing by $(pA_0)!$:
1. **Endpoint = a single term.** The top channel's normalized value is the leading radial coefficient
   $a^{j_{\max}}c^{k_{\max}}$ evaluated — *one* nonzero term, which cannot self-cancel. So if the endpoint is a
   **strict source** (its degree strictly exceeds all others $=$ degree-gap $=$ transitive channel tournament,
   S88), $E[P^m]/(pA_0)!\to$ that term $\ne0$: **noncancellation, immediately** (this is codex THM-2017 and the
   S88 transitivity, made a one-liner by the normalization).
2. **The normalized sum is the peeled Vandermonde.** boxeph THM-2033: the channel determinant is $\prod a_i!\cdot
   \mathrm{Vandermonde}(\text{radial degrees})$. Dividing by the common factorial $(pA_0)!$ is exactly *removing
   $\prod a_i!$* — leaving the pure Vandermonde $=$ the signed tournament sum (klein THM-1805). So:
   $$\boxed{\ \text{NC2}\iff \text{the (possibly confluent) Vandermonde of the channel radial degrees is nonzero.}\ }$$
3. **The wall is the confluent case.** When degrees tie (resonance $=$ regular/Paley, S89), the ordinary
   Vandermonde vanishes and the normalized sum is the **confluent Vandermonde/Wronskian $=$ central trinomial
   (S90) $+$ hyper-Bessel** (codex); noncancellation $=$ that confluent discriminant $\ne0$ $=$ Laguerre-Pólya
   $=$ Paley spectrum on the critical line (S90).

## The payoff: NC2 stripped to its algebraic core
"Divide by $(pA_0)!$" is the move that **separates the (trivial) factorial growth from the (hard) sign
structure**: the factorial is the common $(pA_0)!$; the sign structure is the Vandermonde discriminant of the
channel degrees. Two regimes, one condition:
- **distinct degrees** (transitive channel tournament, endpoint source): Vandermonde $\ne0$ by inspection $=$
  the degree-gap proved region (THM-2017) $=$ a *single surviving normalized term*;
- **repeated degrees** (regular/Paley wall): confluent Vandermonde, the open residual $=$ central-trinomial /
  hyper-Bessel / Laguerre-Pólya / Paley-spectrum-reality — the four faces (S90).

This is the cleanest statement of the tournament↔NC2 bridge yet: **after dividing out the common factorial
$(pA_0)!$, NC2 is precisely "the channel-degree Vandermonde (transitive tournament sign-sum) is nonzero,"** and
the whole difficulty is concentrated in its confluent (regular/Paley) limit. It also explains the factorial
structure of the EMP floor: the depth grows because the common factorial $(pA_0)!$ grows with the radial
degree/charge, and only after peeling it does the finite discriminant condition appear.

## Honest scope
A development of the owner's normalization into the clean condition "NC2 $\iff$ channel-degree Vandermonde $\ne0$
(confluent at the wall)"; verified $p{=}q{=}1$ and $p{=}2,q{=}1$. Not a proof (the confluent/wall case is the
open residual, as always). Value: it names $(pA_0)!$ as the dominant channel factorial, shows dividing by it
peels boxeph's $\prod a_i!$ to expose the pure tournament discriminant, and reduces NC2 to one algebraic
statement with the difficulty isolated at the regular/Paley confluent limit. Cross-links: boxeph THM-2033/S202,
codex THM-2017/2023, S88 (channel tournament), S89 (wall=regular/Paley), S90 (central trinomial/free-prob), klein
THM-1805/1815, THM-1790 (EMP factorial growth). Script `nc2_divide_common_factorial_deathstar_S91.py` (+out).
HYP-8795.
