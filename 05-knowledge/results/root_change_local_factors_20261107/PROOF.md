# THM-4916: local factors of the two-root-change multiplicity, product families, and the cutoff barrier

kind-pasteur-2026-11-07-S15. Status by section:
- **PROVED**: identities in sections 1 and 3.
- **FINITE-EXACT and pre-registered**: section 2; R1 is in section 4.
- **CONJECTURE**: HYP-9506.

Pre-registration [PREREGISTRATION.md](PREREGISTRATION.md) was pushed in 0cda244632. Companion notes:
[ATLAS.md](ATLAS.md) compares small-number coincidences; [README.md](README.md).

**Concurrent same-prompt work.** Session223
(`05-knowledge/results/session223_root_changes_and_local_factors_20261107/`) has the two-run product rectangle and the
square-sum boundary on 2^a + 1. It lists a "height-weighted 3-adic factor" and the cutoff as under investigation.

**Method credited.**
- Seed confinement mod 3 is the mechanism opus-S82 used for its bank (`local_factor_mod3_prediction.py`), from
  Session193 modular_sync and Session206 weighted capacities.
- The K-form is Session199's.
- The D + 2^k U family is Session200's Theorem 2.

Notation is that of THM-4904. Words over {O,Z,J} act on t = (D,E,U,d) = (1,0,2,1). F(s) = D - E + U, and M_2(N)
counts words with at most two J of value N.

## 1. Exact identities (PROVED)

**(1a) Confinement.** A seed s = (h,R,e,d) has tail values

    H_s(m) = h + R m - d Q(m) - e (m mod 2) = h - e (m mod 2)   (mod gcd(R,d)).

*Proof.* R m - d Q(m) is divisible by gcd(R, d). QED.

For the root-change seed s = J(p) of a state p = (D,E,U,d), we have (h,R,e,d) = (D-E+U, D+U, D-E, D), so
gcd(R, d) = gcd(D, U). Odd tails land in U (mod gcd(D,U)), and even tails land in D - E + U.

**(1b) The tail of a root change in pre-change coordinates.** For every state p = (D,E,U,d) and m >= 0, the value of
the word p J (tail m) is

    m odd:   N = (m - Q(m)) D + (m + 1) U,
    m even:  N = (m + 1 - Q(m)) D + (m + 1) U - E.

The odd case does not involve E or d. Hence gcd(m - Q(m), m + 1) = gcd(Q(m) + 1, m + 1) divides N for every odd tail,
whatever the state. We call this tail confinement.

*Proof.* Substitute (h,R,e,d) = (D-E+U, D+U, D-E, D) into the tail formula. QED.

**(1c) Fermat-type multipliers.** For j >= 2, Q(2^j + 1) = 2^(j-1). So for every state p,

    value(p J Z O^(j-1) Z) = (2^(j-1) + 1) (D + 2U) = (2^(j-1) + 1) G(p).

The multipliers are 3, 5, 9, 17, 33, 65, 129, 257, 513, ... = 2^i + 1, among them the Fermat primes 3, 5, 17, 257.
Together with (1b) at m = 2^j - 1 (Session200's D + 2^j U), every state p with at most one J yields the families
(2^i + 1) G(p) and D + 2^j U inside F_<=2.

Checked on 889 words with exact arithmetic.

**(1d) The defect.** Put a = D - E - U + d. Then:
- Z fixes a;
- J resets a to 0;
- O sends a to a + E - U.

Session199's K-form K(v,w) = v^T Q w, with Q its matrix (2), satisfies

    K(v,w) = U_v U_w + E_v E_w - d_v d_w + a_v d_w + a_w d_v.

Session199 shows K(t.alpha, t.beta) = F(t. beta alpha^dagger). So F_<=2 is the set of values K(v,w) with
#J(v) + #J(w) <= 2.

**(1e) Product family.** Let (x, eps, y, 1) and (x', eps', y', 1) be J-free states, and p, q >= 1. Put
v = (x,eps,y,1) J Z^p and w = (x',eps',y',1) J Z^q. Both have a = 0, so

    (x + 2^p y)(x' + 2^q y') = K(v,w)   lies in F_<=2.

So F_<=2 contains V V, where V = {D + 2^p U : J-free states, p >= 1}. Checked on 3,000 random pairs with exact
arithmetic.

## 2. The local factors (derivation, FINITE-EXACT, pre-registered tests)

**Mechanism.** By (1a)-(1c), N mod q is the image of the word under a finite automaton on (D,E,U,d) mod q. The
automaton is not ergodic. Two kinds of confinement make the classes N = 0 (mod p) heavy:
- **seed confinement**: d = U = 0 (mod p) freezes D mod p under O and Z;
- **tail confinement**: p | gcd(m - Q(m), m + 1) for odd tails.

The exact decomposition at [2^23, 2^24) (`decomp_2p24.txt`) shows both parts:
- seeds with 3 | gcd(R,d) put all their odd-tail mass on N = 0 (mod 3^j), exactly as (1a) says;
- seeds with 3 ∤ gcd still favour N = 0 (mod 3), with a max/min ratio of 1.29 over the classes mod 729. This ratio
  was 1.87 at 2^18 and decreases with scale. It comes from tail confinement (1b).

**Residue automaton (`resdp.c`).** Count all prefixes with n letters O/Z and at most two J, with the M_2 weights, by
class mod q. This is a first-principles model: no census data. Its local factors L_n(r) = q * share(r):

- **mod 27**, L by v3(N) = 0, 1, 2, >= 3:
  - n = 20: 0.800, 1.244, 1.630, 1.877;
  - n = 400: 0.866, 1.205, 1.369, 1.451.
- **mod 3**, Richardson limit from n = 4000, 8000, 16000: L(0) = 1.2592, L(1) = L(2) = 0.8703. The conjectured exact
  values are 34/27 = 1.259259 and 47/54 = 0.870370.
- **asymptotic L(0)** at n = 1000, with the nonzero classes uniform:

  | p | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | L(0) | 1.263 | 1.2125 | 1.1687 | 1.1219 | 1.1023 | 1.0877 | 1.0754 | 1.0644 | 1.0511 | 1.0479 |

  (L(0) - 1) p tends to about 1.5. This is a divisor-type local factor, as expected from the product structure (1e).
- **Finite n (about 20).** The primes that divide a cheap multiplier 2^i + 1 (1c) are boosted strongly: 3, 5, 17 (the
  Fermat primes) and 11 (which divides 33). This effect decays like 1/n.

**Pre-registered census tests** (odd N in [2^23, 2^24), exact 16-bit census, `rc_lf.c`; outcomes in
[PREREG_OUTCOMES_R2_R5.txt](PREREG_OUTCOMES_R2_R5.txt)):

| Test | Outcome |
| --- | --- |
| R2: rho_p = mean(p \| N) / mean(p ∤ N) in the DP band (n = 18..23) +- 3% | **FAILED**. Inside for 7, 23, 31. Above the band for 5 (1.614 vs <= 1.567), 11 (1.374 vs <= 1.357), 13 (1.1710 vs <= 1.1707) and 17 (1.534 vs <= 1.486) |
| R3: rho_5 > rho_17 > rho_11 > max(rho_7, rho_13, rho_23, rho_31) | **CONFIRMED** (1.614 > 1.534 > 1.374 > 1.240) |
| R4: the nonzero classes are uniform, max/min <= 1.10 | **CONFIRMED** (worst 1.045, at p = 7) |
| R5: at least 25 of 27 classes mod 27 in the DP band +- 5% | **CONFIRMED** (27 of 27) |

**Reading.**
- The automaton explains the class pattern: which classes are heavy (R5), uniformity of the rest (R4), and the
  ordering of the primes (R3).
- It underestimates the zero-class boost at census scale for the multiplier primes 5, 11, 13 and 17 (R2).
- The census weights words by value, not by letter count. Words in the multiplier families (1c) have systematically
  different heights for the same letter count, so a height-weighted automaton is needed for quantitative agreement.
  This matches Session223's description of a "height-weighted" factor.

The observed 3-adic profile is 205.8, 332.0, 457.3, 510.7, 556.8, 614.2, 618.7 by v3(N) = 0..6 (S14). It is this
automaton's finite-n pattern, with an asymptotic mod-3 limit near (34/27, 47/54, 47/54).

**Near-misses are rough.** In both lists, the factor 3 never occurs.
- The 15 misses of F_<=2 above 10^5 include 5 primes, and 14 of the 15 are 1 (mod 3).
- The 12 values with M_2 <= 6 in (2^22, 2^28) include 3 primes. The only prime factor below 20 among them is 13, in
  three of them.
- 277,438,727 (M_2 = 8, the auditor's list beyond 2^28) is prime.

The lowest-intensity N are those with no small prime factor, the low points of the divisor-type local factor.

## 3. The cutoff: what is proved, and the barrier

**Goal.** An effective N_0 with every odd N > N_0 in F_<=2. Combined with the census, this would prove that every odd
N > 2647 is a uniquely Hamiltonian count, and hence Chan-Pak. **No proof is given here.** What is proved:

1. **Explicit families inside F_<=2** (1b, 1c, 1e):
   - D + 2^j U for every state with at most one J;
   - (2^i + 1) G(p);
   - the product set V V;
   - more generally the bilinear values K(v,w).
2. **Barrier lemma (elementary).** Every family that is a product (1e) or a multiple (1c) covers only composite N. A
   proof of the cutoff must therefore represent every large prime N through the non-product terms
   E_v E_w - d_v d_w + a_v d_w + a_w d_v of (1d).

   Coverage by one seed, by any finite bank, or by any fixed linear form in one free state is impossible:
   - each single tail image misses a positive proportion (HYP-9480 (A), (C));
   - Session217 shows a finite Q/Z bank has only O((log N)^2) outputs up to N.

   So every proof needs the joint freedom of two states.
3. **Why the Gerbicz template does not transfer directly.** The square-sum theorem has a finite base (41..2032) and a
   self-similar 49-fold blow-up n -> 49n + res, by Gerbicz in 2018 (`ATLAS.md`). An analogue here needs an affine
   self-map N -> aN + b of F_<=2 with a >= 2. The value N is not a sufficient statistic for the state (D,E,U,d), and
   every suffix acts linearly on the state, not on N. So a blow-up must act on state cones: a set C of states with
   K(C, C) closed under a fixed pair of word insertions that multiply K by a constant. Session216's digit probes
   (H -> 2H +- 1 macros) found none in the two-apex grammar; that is the closest prior negative evidence.
4. **What would suffice.** One of the following:
   - a state cone C, closed under a word map that scales K, together with a finite base verified by the census;
   - an analytic lower bound M_2(N) >= 1, which needs a second-moment (collision) estimate for the bilinear form over
     two independent states.

   With either, the census to 2^32 (section 4) would supply the finite part.

## 4. Census to 2^32 (R1)

[filled from `census_2p32.txt`: the auditor's independent segmented enumerator, segments of 2^28 from 805,306,368]
