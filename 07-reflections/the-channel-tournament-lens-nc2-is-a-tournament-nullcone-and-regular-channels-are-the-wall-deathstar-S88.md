# The channel-tournament lens: NC2 is a tournament-nullcone on its radial channels, and the regular-channel core is the wall

> **CORRECTION (MISTAKE-212 / THM-2021).** The claimed equivalence
> "NC2 noncancellation iff channel transitivity/source dominance" is false.
> A strict source is a sufficient domination certificate, not a necessary one.
> On the fully tied central slice `h=kappa*b^r`, THM-2021 factors the sum as a
> toral return polynomial times `L(b^m)` and proves noncancellation generically.
> Explicitly, `b=s(s-2), a=s-2, c=s(s-2)` has `h=b^2`, `E[P]=0`, and
> `E[P^2]=24` despite identical factorial degree in every channel. Exact ties
> also require an extrinsic tie-breaker before the relation is a tournament.
> Retain this reflection as a useful regime-classification lens and an
> explanation of why domination stops at resonance; withdraw its iff and
> "moment-nullcone is a tournament-nullcone" assertions.

**death-star-2026-07-21-S88** (HYP-8772). Owner: work on new lenses and unifying theorems for NC2. Here is a
new lens that connects NC2's open residual (the radial-channel noncancellation, HYP-8765/8766/8770) to the
repo's tournament framework and to the "regular/Paley is the wall" meta-principle — and, as a bonus, explains
*why the domination route was refuted* (MISTAKE-202).

## The lens: NC2's channels form a tournament
For the single-straddle $P=Z\,A(s)+B(s)+\bar Z\,c_0$, $E[P^m]=\sum_i \binom{m}{i,i,m-2i}L(s^iA^ic_0^iB^{m-2i})$
is a sum over primitive-return **channels** $i$. codex already noted these channels carry a *tournament*
structure (the "channel-degree tournament", observable $D(k)-D(l)$, sign gauge, one Hamiltonian path). Make it
the lens: **orient $i\to j$ iff channel $i$ dominates channel $j$ as $m\to\infty$.** NC2 (noncancellation) becomes a
**transitivity statement** on this tournament. Verified (`nc2_channel_tournament_lens_deathstar_S88.py`):

- **Degree-gap regime $\Rightarrow$ TRANSITIVE channel tournament.** One channel outgrows the rest by factors of
  $10^9$ and up ($P=Z s^2+1+\bar Z s^2$, $m=6,8,10$): a **clear source**, a total order, the single Hamiltonian
  path codex saw. A transitive tournament is the **S75 nullcone vertex ($X^n$)**: the source channel's sign
  survives, so $E[P^m]\ne0$ — **noncancellation**. This *is* codex's THM-2017 ("one endpoint ratio one").
- **Resonance band $\Rightarrow$ TIED (REGULAR) channels.** At the central offset ($P=Z(1{+}s)+(1{+}s)+\bar Z$),
  the top/second channel ratio $\to1$ ($1.70,1.33,1.09$ as $m$ grows): the domination order **degenerates to a
  balanced/regular core**. That is the **wall** — the same regular/doubly-regular/Paley configuration that is the
  tightest case for H≥disc (S84) and the LRC AP/Paley pole (S75). Only here can channels cancel, and only if their
  **signs and phases conspire** (S87).

## The unifying theorem (statement)
> **NC2's radial-channel residual is a tournament-nullcone: $E[P^m]$ fails to vanish for all $m$ iff the
> channel-degree tournament is transitive-dominated (has a strict source). The proved cases (degree-gap,
> THM-2017) are exactly the transitive channel tournaments (S75 vertex); the open case (resonance band) is
> exactly the regular/tied channel tournament (the S84/S75 wall).**

So the **holonomic rung of the moment-nullcone ladder (NC2) reflects into a tournament-nullcone on its own
channels** — and more generally every rung (kp THM-1775: rational/trace, algebraic/CT, holonomic/E) carries a
channel tournament whose transitivity is the nullcone-detection. The **open case is always the regular-channel
configuration**, i.e. the repo's universal "the maximally-symmetric object is the wall" (S76 big-stabilizer
lens) now realized *inside* the moment functional itself.

## Why the domination route was refuted (MISTAKE-202), explained by the lens
The refuted "top-term dominates the $r$-average" step (klein-S351, refuted by kp-S128c120, my MISTAKE-202
concession) is precisely **"assume the channel tournament has a source."** It is *true* in the degree-gap
(transitive) regime — which is exactly where THM-2017 succeeds by a domination estimate — and *false* in the
resonance band, because there the channels are **regular** (no source, ties). So the tournament lens gives a
one-line reason the whole domination program split the way it did: **domination = channel-source = transitivity;
it works iff the channels are not regular.** The resonance band is the regular locus where it must fail, so the
correct weapon there is not domination but a **positivity/SOS** certificate (S87: cancellation needs sign-
intransitivity among the tied channels; an SOS/Hankel-PD form forbids it), or the analytic no-common-zero
(codex hyper-Bessel = my S62/S64 Sheffer).

## Two more lenses (shorter)
- **Free-probability lens (via THM-438).** The channel weights are free cumulants of the two-point spectrum; the
  regular (tied) core is a **free-convolution semicircle regime**, and the central-offset "entropy saddle"
  (codex) is a free-probability rate function. Prediction: the resonance-band asymptotics have a closed form as
  a free-cumulant generating series — the same Wigner/Catalan machinery (THM-438) that gave $H(\text{Paley})\sim
  e\cdot\text{avg}$. (Unifies the NC2 residual with the Paley/quasirandom recovery, S85.)
- **Rigidity lens (via S86).** NC2 $=$ "$P_*(\text{Gaussian})$ with vanishing analytic moments is one-sided"; the
  channel tournament is the *combinatorial skeleton* of that pushforward measure, and its transitivity is the
  rigidity. Null-quadrature-domain classification (Sakai) would be the analytic counterpart of "the channel
  tournament is transitive."

## Honest scope
Lenses + a unifying statement, not a proof: the channel-tournament transitivity is *verified* to track NC2's
proved/open split (degree-gap transitive, resonance regular) and to *explain* the domination refutation, but the
resonance (regular-channel) case is still open — as it must be, being the wall. Value: it places NC2's residual
inside the repo's tournament + "regular-is-the-wall" framework, unifies four attack lines (domination = source,
positivity = anti-sign-intransitivity, hyper-Bessel = no-common-zero, free-prob = tied-core asymptotics) as
facets of one channel tournament, and points the resonance band at SOS + free-probability. Cross-links: codex
THM-2017/HYP-8766, S87 (Sheffer/positivity), S86 (null-quadrature), S84 (regular=wall for H≥disc), S75 (nullcone
vertex vs Paley pole), S76 (big-stabilizer lens), THM-438 (Wigner/free cumulants), THM-1775 (ladder), THM-1810
(bosonic), MISTAKE-202. Script `nc2_channel_tournament_lens_deathstar_S88.py` (+out). HYP-8772.
