# THM-4916: local factors of the two-root-change multiplicity, product families, and the cutoff barrier

kind-pasteur-2026-11-07-S15. Status by section:
- **PROVED**: section 1 (identities, mostly KNOWN and credited) and section 2.1 (the closed form for every odd prime,
  in the letter-count model; new for general p).
- **FINITE-EXACT and pre-registered**: section 2; R1 is in section 4.
- **CONJECTURE**: HYP-9506.

The independent audit is PASS WITH CORRECTIONS, applied ([audit/AUDIT.md](audit/AUDIT.md); MISTAKE-734). The
pre-registration [PREREGISTRATION.md](PREREGISTRATION.md) was pushed in 0cda244632. Companion notes:
[ATLAS.md](ATLAS.md) compares small-number coincidences; [README.md](README.md).

**Concurrent same-prompt work.** Session223 (`05-knowledge/results/session223_root_changes_and_local_factors_20261107/`)
has the two-run product rectangle and the square-sum boundary on 2^a + 1. Its checkpoint 2 (1c524e45a9, 13:31:03,
pushed 32 s before this packet's checkpoint) PROVES the fixed-length 3-adic law at every precision, by the
Cesaro-projector method of section 2.1:
- mod 3: 34/27, 47/54, 47/54;
- with one J: 41/36, 8/9, 35/36.

It also gives the exact two-J membership equation, whose eq. (2) is (1b). The height-weighted factor and the cutoff
remain open there.

**Method credited.**
- Seed confinement mod 3 is the mechanism opus-S82 used for its bank (`local_factor_mod3_prediction.py`), from
  Session193 modular_sync and Session206 weighted capacities.
- The K-form is Session199's.
- The D + 2^k U family is Session200's Theorem 2.

Notation is that of THM-4904. Words over {O,Z,J} act on t = (D,E,U,d) = (1,0,2,1). F(s) = D - E + U, and M_2(N)
counts words with at most two J of value N.

## 1. Identities (PROVED; priority as stated)

The exact checks are reproducible with `audit/identities_audit.py` (3,000 to 1,740,800 cases each).

**(1a) Seed confinement (KNOWN: Session206 weighted-parameter-capacity and exact-modular-image; THM-4888 (C)).** A seed
s = (h,R,e,d) has tail values

    H_s(m) = h + R m - d Q(m) - e (m mod 2) = h - e (m mod 2)   (mod gcd(R,d)).

For the root-change seed s = J(p) of a state p = (D,E,U,d), (h,R,e,d) = (D-E+U, D+U, D-E, D), so gcd(R,d) = gcd(D,U).
Odd tails are divisible by gcd(D,U), since h - e = U. Even tails are congruent to D - E + U, that is to -E, mod
gcd(D,U).

**(1b) The root-change tail in pre-change coordinates (concurrent with Session223 cutoff/REPORT.md eq. (2)).** For every
state p and m >= 0, the value of the word p J (tail m) is

    m odd:   N = (m - Q(m)) D + (m + 1) U,
    m even:  N = (m + 1 - Q(m)) D + (m + 1) U - E.

New here: this m, Q(m) form, and the corollary that gcd(Q(m) + 1, m + 1) divides every odd-tail value (tail
confinement).

**(1c) Fermat-type multipliers (a one-line corollary of s.JZ = s.OJ and of the mechanism of Session223's rectangle).**
For j >= 2 (it fails for j = 1) and every state p, Q(2^j + 1) = 2^(j-1), so

    value(p J Z O^(j-1) Z) = (2^(j-1) + 1) (D + 2U).

**(1d) The defect (Session199's coordinate, a = -f).** a = D - E - U + d is fixed by Z, reset to 0 by J, and sent to
a + E - U by O. Session199's form satisfies

    K(v,w) = U_v U_w + E_v E_w - d_v d_w + a_v d_w + a_w d_v

and K(t.alpha, t.beta) = F(t. beta alpha^dagger). So

    F_<=2 = {1} ∪ {K(t.alpha, t.beta) : #J(alpha) + #J(beta) <= 2}.

**(1e) Product family (Session199 physical-adjoint section 4 on balanced pairs).** For J-free states (x,eps,y,1) and
(x',eps',y',1) and p, q >= 1 (it fails for p = 0),

    (x + 2^p y)(x' + 2^q y') = K((x,eps,y,1) J Z^p, (x',eps',y',1) J Z^q)   lies in F_<=2.

The realizing word has exactly two J. New as stated: the J-budget form V V in F_<=2, with
V = {D + 2^p U : J-free, p >= 1}.

## 2. Local factors: the residue automaton and the census (FINITE-EXACT, pre-registered)

**Mechanism.** N mod q is the image of the word under a finite automaton on (D,E,U,d) mod q. Seed confinement makes the
classes N = 0 (mod p) heavy; tail confinement contributes only at finite length.

**Exact decomposition at [2^23, 2^24)** (`decomp_2p24.txt`):
- Seeds with 3 | gcd(R,d) put all their odd-tail mass on N = 0 (mod 3^j), exactly as (1a) says.
- Seeds with 3 ∤ gcd still favour N = 0 (mod 3):
  - their max/min over the classes mod 729 is 1.286 at 2^24, against a sampling-noise floor of 1.04. At 2^18 it is
    1.870 against a noise floor of 1.50, so that larger value is mostly noise;
  - their class-0-mod-3 ratio falls from 1.126 (even tails) and 1.256 (odd tails) at 2^18 to 1.088 and 1.166 at 2^24.
- Tail confinement alone has no net effect for uniform residues. At census scale, the odd-tail bias comes from an
  over-representation of tail-confined tails (15.3% against 1/9). The even tails are biased by a finite-size effect
  not identified here.
- In the letter-count limit, the whole zero-class excess is seed confinement at the last root change (section 2.1).

**Residue automaton (`resdp.c`).** The first-principles model counts all prefixes with n letters O/Z and at most two J,
with the M_2 weights, by class mod q. Its local factors are L_n(r) = q times the share of class r.
- At n = 20, by v3(N) = 0, 1, 2, >= 3 (mod 27): 0.800, 1.244, 1.630, 1.877.
- Its n -> infinity limit is the closed form of section 2.1: L_p(0) = 1 + (3p-2)(p-1)/(2p^3). The values are 34/27,
  151/125, 400/343, 1486/1331, 2419/2197, 5305/4913, 7354/6859 and 12904/12167 for p = 3, 5, 7, 11, 13, 17, 19, 23.
- (L(0) - 1)p = (3p-2)(p-1)/(2p^2), which tends to 3/2. It is 1.420 at p = 31.
- The values in `dp_primes_n1000.txt` are finite-n values with an O(1/n) bias, not the limits.
- For p = 3, Session223 proved these limits at every 3-adic precision; the letter-count limit coincides with its
  fixed-length law.
- The excess is seed confinement at the last root change; product words are a vanishing fraction of all words.

**Pre-registered census tests** (odd N in [2^23, 2^24), exact 16-bit census, `rc_lf.c`; outcomes in
[PREREG_OUTCOMES_R2_R5.txt](PREREG_OUTCOMES_R2_R5.txt); reproduced by the audit):

| Test | Outcome |
| --- | --- |
| R2: rho_p = mean(p \| N) / mean(p ∤ N) in the DP band (n = 18..23) +- 3% | **FAILED**. Inside for 7, 23, 31. Above the band for 5 (1.614 vs <= 1.567), 11 (1.374 vs <= 1.357), 13 (1.1710 vs <= 1.1707) and 17 (1.534 vs <= 1.486) |
| R3: rho_5 > rho_17 > rho_11 > max(rho_7, rho_13, rho_23, rho_31) | **CONFIRMED** (1.614 > 1.534 > 1.374 > 1.240) |
| R4: the nonzero classes are uniform, max/min <= 1.10 | **CONFIRMED** (worst 1.045, at p = 7) |
| R5: at least 25 of 27 classes mod 27 in the DP band +- 5% | **CONFIRMED** (27 of 27). Only 9 of 27 lie inside the unwidened band: at every class the census is more contrasted than the automaton, in the same direction as R2's failure |

**Reading.**
- Qualitatively, the census shows the automaton's finite-n pattern.
- It is more contrasted than the letter-count model at every class.
- The census weights words by value, not by letter count. Whether the height-weighted census factor tends to the
  closed form is OPEN. Session223 section 2 shows that residues correlate with real-state history. At [2^23, 2^24) the
  census factor at N = 0 (mod 3) is 1.444.

**Where the zero-class excess comes from at census scale** (audit M3, `audit/census_mechanism.c`, [2^23, 2^24)):

| p | 3 | 5 | 17 | 11 | 13 | 7 |
| --- | --- | --- | --- | --- | --- | --- |
| root-rectangle words Z^a O^b J, divisible by 2^a + 1 | 24% | 21% | 19% | 9% | 11% | 0 |
| (1c) multiplier words | 24% | 23% | 22% | 12% | 12% | 5% |
| seed confinement at the last root change | 54% | 56% | 48% | 70% | 97% | 101% |

The R3 ordering persists with both Fermat-type families removed: rho = 1.363, 1.339, 1.324, 1.247 for p = 5, 17, 11, 7.

**Near-misses.**
- The 130 values with M_2 <= 12 in (2^22, 2^28) avoid 3 and are depleted in 5, 7, 11, 17 and 19. Only 21 of 130 have a
  prime factor in 5..19, against 49% of random odd numbers prime to 3, and none is divisible by 5 or 17. They are not
  depleted in 13: 11 of 130 are divisible by 13, against 10.0 expected.
- Of the 15 misses of F_<=2 above 10^5, 14 are 1 (mod 3) and none is 0. Five are prime, and six have a prime factor
  of at most 19, for example 116,215 = 5 * 11 * 2113.

## 2.1 The closed form for every odd prime (PROVED, letter-count model)

**Lemma (M), block mixing.** Proved by the audit for every odd prime p (audit section 5, by an SCC/doubly-stochastic
argument and characteristic functions). For a long O/Z block with i.i.d. uniform letters, in the Cesaro limit:
- **(M1)** If d != 0 (mod p), the final (D,U) is uniform on (Z/p)^2 and independent of E, which is 0 or d with
  probability 1/2 each. So F is uniform.
- **(M2)** If d = 0 and U != 0 at the block start, then E = 0 after the first letter, U stays nonzero, D is uniform,
  and F is uniform.
- **(M3)** If d = U = 0 at the start, right after a J (so D = E), the value class freezes after the first letter: F = 0
  if the first letter is Z, and F = E_J if it is O.

**Reduction (Session223 section 5; justified in the audit for every odd prime).**
- Prefixes with at most one J weigh O(n 2^n), against C(n+2,2) 2^n for exactly two J.
- Exactly-two-J prefixes are uniform over the block lengths.
- The first block starts in the aperiodic d = 1 class.

So L_n = readout(Pi J Pi J Pi delta_t) + O(1/n).

**Cascade.**
1. **Block 1** (d = 1): (D1,U1) is uniform, and E1 is in {0,1}.
2. **J1** gives d2 = D1 and U = D1 + U1. Three cases:
   - **A** (D1 != 0, probability (p-1)/p): (D2,U2) is uniform, and E2 is in {0, D1}.
   - **B** (D1 = 0, U1 != 0, probability (p-1)/p^2): D2 is uniform, U2 != 0, E2 = 0.
   - **C** (D1 = U1 = 0, probability 1/p^2): by (M3) with E_J = -E1, D2 = 0 with probability 3/4, U2 = 0, E2 = 0.
3. **Confinement at J2:**

       c = P(D2 = U2 = 0) = (p-1)/p^3 + (3/4)/p^2 = (7p-4)/(4p^3),
       e = P(E2 = 0 and confined) = (5p-2)/(4p^3).

4. **Block 3**, by (M1)-(M3): P(F = 0) = (1/p)(1 - c) + c/2 + e/2. Hence

       L(0) = 1 + c(p/2 - 1) + (p/2) e = 1 + (3p-2)(p-1)/(2p^3).

   The confined value -E2 = -D1 is uniform over the nonzero classes, so L(r) = 1 - (3p-2)/(2p^3) for r != 0.

The same computation with one J gives L(0) = 1 + (3p-4)/(4p^2) and L(-1) = 1 - 1/p^2 + 1/(4p); the other classes have
1 - 1/p^2.

**Verified.**
- Exact rational arithmetic at p = 3, 5, 7 (`audit/cesaro_exact.py`).
- Within 1.2e-12 for every prime p <= 31, by two projector methods (`exact_limit.py`, `audit/cesaro_lazy.c`).
- The mod-15 limit is not the product of the mod-3 and mod-5 factors. The ratios are 1.0088, 0.9972, 0.9937 and
  1.0020, so the local factors are not exactly multiplicative across primes.

## 3. The cutoff: what is proved, and the barrier

**Goal.** An effective N_0 with every odd N > N_0 in F_<=2. Combined with the census, this would prove that every odd
N > 2647 is a uniquely Hamiltonian count, and hence Chan-Pak. **No proof is given here.**

1. **Explicit families inside F_<=2** (1b, 1c, 1e): D + 2^j U for every state with at most one J; (2^i + 1) G(p); the
   product set V V; and in general the bilinear values K(v,w).
2. **Barrier lemma (elementary).** Every product family (1e) and every multiplier family (1c) consists of composite
   numbers. So for a prime N every representation N = K(v,w) has a nonzero non-product part
   E_vE_w - d_vd_w + a_vd_w + a_wd_v. If that part vanished, N = U_vU_w with U_v, U_w >= 2.

   Related facts, with their status:
   - a single O/Z seed misses a positive proportion of odd integers (PROVED: Session206 parameter-capacity and
     root-residue-capacity);
   - whether every finite O/Z bank does is OPEN (HYP-9480 (C));
   - Session217's O((log N)^2) bound is for banks of the commuting Q/Z construction, not for O/Z tail banks.

   One free state is not excluded: K(t,s) = F(s), so F_<=2 \ {1} is the image of a single free state with at most two
   J. Preferring two-state arguments is a heuristic choice.
3. **Why the Gerbicz template does not transfer directly.** The square-sum theorem has a finite base (41..2032) and a
   self-similar blow-up n -> 49n + res (Gerbicz 2018; ATLAS.md). An analogue here needs an affine self-map N -> aN + b
   of F_<=2 with a >= 2. The value is not a sufficient statistic for the state (D,E,U,d), and every suffix acts linearly
   on the state, not on N. So a blow-up must act on state cones: a set C of states with K(C, C) closed under a fixed
   pair of word insertions that multiply K by a constant. Session216's digit probes (H -> 2H +- 1 macros) found none in
   the two-apex grammar; that is the closest prior negative evidence.
4. **What would suffice.** One of the following:
   - a state cone with such a scaling word map, together with a finite base verified by the census;
   - an analytic lower bound M_2(N) >= 1, which needs a second-moment (collision) estimate for the bilinear form.

   With either, the census (section 4) would supply the finite part.

## 4. Census to 2^32 (R1)

Pending at the time of the audit. The result is recorded in [R1_OUTCOME.md](R1_OUTCOME.md) when the run (the S14
auditor's independent segmented enumerator, segments of 2^28 from 805,306,368) completes.
