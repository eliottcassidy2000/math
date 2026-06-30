# Rédei is Brouwer-odd, LRC is Borsuk–Ulam-even: a Tucker / Ky-Fan / Ham-sandwich reading, and why the lonely count is provably even

*opus-2026-06-29. Owner: think Ky-Fan alternating count, Tucker, Ham sandwich. These are the
Borsuk–Ulam combinatorial family, and they pin down exactly why LRC needs the *signed/alternating*
certificate (the R-odd eigenspace, the i√p odd degree) where Rédei needs only a parity.*

## The clean elementary result (proved + verified)
> **For every covering 13-set `S` (it contains an even speed) and every modulus `D`, the lonely
> count `N(S,D) = #{a∈Z/D : ‖sa/D‖≥1/14 ∀s}` is EVEN.**

*Proof.* `a` lonely ⟺ `−a` lonely (`‖s(−a)/D‖=‖sa/D‖`). The only `a` with `a=−a` are `a=0`
(`‖0‖=0`, never lonely) and `a=D/2` (D even; `‖s_0·D/2/D‖=‖s_0/2‖=0` for the even speed `s_0`, so not
lonely). No fixed lonely point ⇒ lonely points pair up ⇒ `N` even. ∎ (Verified across D=27..89: always
even, always mirror-paired, no fixed lonely point.)

**Consequence.** The LRC witness for a covering set is genuinely an **antipodal pair `{t,−t}`** with
no symmetric representative — so the certificate is **Borsuk–Ulam (even index), not Brouwer (odd
fixed point)**. This is the rigorous basis for kps's "Borsuk–Ulam not Brouwer" (S31av) and refines my
prime-density criterion: `N(S,p)≥1 ⟹ N(S,p)≥2` (a pair).

## The two-index dichotomy (= mac-mini THM-582), now grounded
| | parity of the count | certificate | why |
|---|---|---|---|
| **Tournament (Rédei)** | `H(T)` **odd** | Brouwer / parity: `H≡1 (mod 2) ⇒ H≥1` ⇒ a Ham path exists | a Ham path's reverse is a *different* path; the diagonal `t=t` is unobstructed |
| **LRC (lonely)** | `N(S,D)` **even** | Borsuk–Ulam: need a *nonzero signed/alternating count*, not a parity | witnesses are antipodal pairs; `N` even could be `0` — parity is silent |

So **Rédei is the Brouwer/odd shadow; LRC is the Borsuk–Ulam/even shadow** — and the even parity is
*exactly why LRC is hard*: a mod-2 argument (which proves Rédei in one line) is vacuous, because the
lonely count being even is consistent with `0`.

## Where Ky Fan, Tucker, Ham sandwich each enter
- **Tucker** (combinatorial Borsuk–Ulam): label `a∈Z/D` by the *signed tightest runner*
  `λ(a)=±s`, with `λ(−a)=−λ(a)`. If the danger sets *cover* (no lonely point), `λ` is a total
  `Z₂`-map and Tucker forces *complementary edges* (the danger-arc centers). But in dimension 1 this
  is trivial (arcs always have centers) — verified: comp-edge count `≈ D/14`, nonzero even when
  covered (D=45). **The honest content is higher-dimensional:** LRC lives on the `(k−1)`-torus (the
  view-obstruction diagonal), and the real certificate is the Ky-Fan count on *that* sphere.
- **Ky Fan** (the *alternating count*): the signed number of alternating simplices is **odd** — a
  robust topological degree that survives where a plain parity dies. This is the certificate the EVEN
  case needs; it is the discrete incarnation of the **i√p "odd degree"** (kps) and of the **R-odd
  eigenspace's nonzero signed residual** (mac-mini HYP-3538). The tournament analogue — the *signed*
  HP count `S(T)=Σ_P sgn(P)` — is always **odd** (`S≡H≡1 mod 2`, Rédei), i.e. the alternating count is
  automatically nonzero on the Brouwer side; on the LRC side it is the open question.
- **Ham sandwich / necklace splitting** (the measure side): the 13 danger sets have total measure
  `13·(1/7)=13/7 > 1`, so they *must* overlap (the union bound is vacuous) — the same `13/7>1`
  obstruction I hit in the exp-sum work. Fair-division/necklace bounds (Borsuk–Ulam corollaries) count
  the cuts; here the "cuts" are the `Σ_s s` arc boundaries, and the lonely point is a gap the
  overlap-forced cover must still leave.

## Net
The three tools say one thing: **LRC is the even/antipodal/Borsuk–Ulam member of a family whose
odd/Brouwer member is Rédei.** Rédei closes by parity (`H` odd); LRC cannot (`N` even) and instead
demands the *signed alternating count* (Ky Fan) be nonzero — which is precisely the R-odd eigenspace /
i√p odd-degree / disjoint-family residual that every route this arc has identified as the irreducible
core. The Tucker/Ky-Fan frame doesn't close it, but it **names the certificate type exactly** and
proves the parity (`N` even) that forces it.

## Status
- **Proved/verified:** `N(S,D)` even for covering sets (mirror-pair argument); Rédei-odd vs lonely-even
  dichotomy; the 1-dim Tucker triviality; the `13/7>1` overlap.
- **Open (the core):** the higher-dimensional Ky-Fan alternating count / R-odd signed certificate —
  nonzero ⇒ LRC(14). Same wall, now correctly typed as a *signed degree*, not a parity.

Related: mac-mini THM-582 (two-index), HYP-3537 (cap = measure-valued Claim A), HYP-3538 (R-odd),
kps S31av (Borsuk–Ulam), THM-001 (Rédei), THM-088 (signed F), LEM-003, OPEN-Q-108.
