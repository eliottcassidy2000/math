# THM-430 — The antipodal involution σ:x↦−x unifies the shell-partner floor (mod q) and the torsion leak (mod n); the signed floor is always a genuine σ-orbit, never the half-turn

**Status:** PROVED (σ-orbit recast, the never-self-partner law, the CRT fiber split, the
odd-shell/even-clock face split) + VERIFIED (n=5,6,7 minimizers + the 5 published HYP-2296
minimizers + the face table n=5..15). Answers **poke Steering Task 1** (connect the binding
shell-partner `q=a+b` to the prime-torsion / fiber picture).
**Source:** opus-2026-06-06-S702. Builds on THM-425 (synchronization), THM-426 (cuts), THM-429 /
HYP-2296 (signed floor = max-cut LRC, binding pair `a+b=q=denom t*`), THM-421/427 (torsion leak),
THM-428 (two-tower clock⊗shell). Convention as THM-429.

## The unifying object

> The norm `‖·‖` is invariant under the **antipodal involution** `σ_N : x ↦ −x` on `ℤ/N`
> (`‖−x‖=‖x‖`). Every loneliness modulus in the LRC apparatus — clock `n`, shell `2n−1`, and the
> *actual optimal witness denominator* `q = a+b` (HYP-2296) — carries this same `σ`. The fixed-point
> set of `σ_N` is the **2-torsion** `Fix(σ_N) = {x : 2x≡0} = T_2^{(N)} = {0, N/2}` (the second
> element only when `N` even). This single involution organizes both the shell-partner floor and the
> torsion leak.

## The theorem

> **(A) Synchronization = σ-invariance (recast of THM-425 L0).** `a+b ≡ 0 (mod q) ⟺ b ≡ σ_q(a)`, so
> a `q`-shell-partner pair `{a,b}` is precisely a **`σ_q`-orbit**, and `‖a·k/q‖ = ‖b·k/q‖` for all
> `k` is exactly the `σ`-invariance of `‖·‖`. Shell-partners *are* the antipodal orbits.
>
> **(B) Self-partners = the 2-torsion (the half-turn).** A pair is a **self-partner** (`a≡b`, i.e. a
> single relative speed equal to its own partner) iff `a ∈ Fix(σ_q) = {0, q/2}`. A nonzero
> self-partner exists **iff `q` is even**, and it is the half-turn `q/2`. (n=14 clock: `q=n=14`,
> half-turn `7 = 14/2` — poke's leak; verified `q=20→{0,10}`, `42→{0,21}`, `8→{0,4}`.)
>
> **(C) The signed floor is NEVER a self-partner (PROVED).** The binding pair of `Gstar(S)` is always
> a *genuine* `σ_q`-2-orbit (`a≢b mod q`). *Proof.* A half-turn relative speed `w=q/2` gives
> `‖w·k/q‖ = ‖k/2‖ ∈ {0, ½}`: it is `0` when `k` even (would kill loneliness, so the maximin forces
> `k` odd) and `½` when `k` odd (the global maximum of `‖·‖`). Either way `w=q/2` is never the
> *minimising* (binding) speed of a positive floor `M=k/q < ½`. ∎ So the floor lives on the genuine
> orbits, and the half-turn 2-torsion — the very leak poke names — is structurally *excluded* from
> setting the signed floor. **Verified:** all 12 minimizers (searched + published) have
> `self_partner = False`.
>
> **(D) CRT fiber split — σ is trivial on the 2-fiber (PROVED + the answer to "mod 2 vs mod p").**
> `σ` commutes with CRT: `σ_q ≅ ∏_p σ_{p^{e_p}}`, so a shell-partner decomposes into per-prime
> `σ`-orbits. On the **2-fiber `σ_2 = id`** (`−1≡1 mod 2`), so the pair is *always* a (trivial)
> self-partner mod 2; the **genuine antipodal content lives in the odd-prime fibers.** Verified: for
> `n=7`, `{19,23}`, `q=42=2·3·7` — `(a,b) mod 2 = (1,1)` self/fixed, but `mod 3 = (1,2)` and
> `mod 7 = (5,2)` are genuine orbits. So poke's "mod 2 vs mod 7" question resolves as: *mod 2 is
> always the trivial fiber; the antipodal (shell-partner) structure is an odd-prime phenomenon.*
>
> **(E) Face split via parity (PROVED + VERIFIED n=5..15).** The shell modulus `2n−1` and the
> consecutive-block witness `q = 4n−5` are **always odd**, so `Fix(σ) = {0}` only — the shell face is
> **antipodally free** (every nonzero residue in a genuine 2-orbit; `(N−1)/2` orbits). The clock `n`
> is even for half the cases, and then carries the half-turn 2-torsion `n/2`. So the prime `2` acts
> oppositely on the two faces: it can be a *clock*-tower prime creating the half-turn self-leak
> (n even), but **never appears on the shell** as an antipodal fixed point.

## Why it matters — the unified answer to Task 1

- **One involution, three moduli.** The shell-partner floor (mod `q`, THM-425/429/HYP-2296) and the
  torsion leak (mod `n`, THM-421/427) are the *same* antipodal structure `σ` read on different
  moduli. Synchronization, the half-turn, and the 2-torsion are all `σ`-statements.
- **The half-turn is a red herring for the floor (C).** Poke's n=14 half-turn `r=7` is the clock's
  `σ`-fixed 2-torsion — maximal-leak in the *count* sense (THM-427 C2), but it can **never set the
  signed floor**, which is always a genuine odd-fiber orbit. The two roles of the 2-torsion (maximal
  cell-leak vs. excluded-from-floor) are reconciled by `σ`: a fixed point leaks maximally but is
  never a binding minimiser.
- **Antipodal structure is odd-prime (D,E).** Because `σ_2 = id`, all genuine shell-partner /
  antipodal content sits at odd primes — which is exactly *why* the shell modulus is forced to be the
  odd `2n−1`, and why (THM-428) n=14's hardness is the **odd** prime-cube `3³` shell tower, not the
  `2` in its clock. The antipodal involution explains the parity asymmetry of the two-tower group.

## Scope / honesty

- (A)(C)(D)(E) are proved; (B) is a definition + the parity fact. The minimizers at `n=6,7` are
  **`B`-limited** (`B≤11`,`≤10`): `3/20`, `1/8` are small-`B` floors, not the true `inf` (HYP-2296
  Part C shows it keeps dropping) — but the `σ`/torsion structure of the *binding pair* is
  `B`-independent and holds for every minimizer found.
- An **observed-not-proved** correlation: in the sample the optimal `q` tends to share a prime with
  the *clock* `n` (`gcd(q,n)>1` common: q=25,20,5,20,27,42) while `gcd(q,2n−1)=1` throughout —
  suggesting the witness denominator aligns with the clock face. Recorded as part of **HYP-2297**;
  needs a larger census before any claim. This does **not** resolve any open LRC case.

**Tournament reading (directive):** `σ_N` is the order-2 automorphism (`−1`) of the cyclic
group/round tournament `C_N`; shell-partners are its 2-orbits, self-partners its fixed points. The
odd-shell / even-clock split is the statement that the round tournament `C_{2n−1}` is
**self-converse with no central involution** (odd order), whereas `C_n` (n even) has the central
half-turn — the same parity that makes `2n−1` the self-converse Farey modulus (THM-401).

**Artifacts:** `04-computation/lrc_antipodal_shellpartner_torsion_s702.py` (+`.out`). Reflection
`07-reflections/the-antipodal-involution-unifies-shell-and-leak-s702.md`. New: **HYP-2297**. Builds
on THM-425/426/429/HYP-2296, THM-421/427/428, THM-401/403.
