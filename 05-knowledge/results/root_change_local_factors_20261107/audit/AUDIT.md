# Independent audit: THM-4916 and HYP-9506 (packet `root_change_local_factors_20261107`)

Independent adversarial audit, 2026-11-07, of `README.md`, `PROOF.md`, `ATLAS.md`, `PREREGISTRATION.md`,
`PREREG_OUTCOMES_R2_R5.txt`, the C/Python files, `01-canon/theorems/THM-4916-...md`, `05-knowledge/hypotheses/HYP-9506-...md`,
and the two coordinator addenda (`exact_limit.py`, `exact_limit.out`, `exact_limit_q15.out`, `closed_form_check.out`, and the
three-block derivation with lemma (M)). Every number below was recomputed with the auditor's own code in this folder. The
packet's code was read only to compare definitions.

## Verdict: PASS WITH CORRECTIONS

- **Mathematics: correct.** The five identities (1a)-(1e) were re-derived by hand and confirmed by exact integer checks
  (3,000 to 1,740,800 cases each).
- **Registered numbers: all reproduced.** Every registered census number matches to print precision. R2 FAILED and R3,
  R4, R5 CONFIRMED, exactly as reported. The pre-registration precedes the outcomes and was not changed afterwards.
- **Corrections are needed in three areas:**
  - **Priority.** Most of THM-4916 is KNOWN or follows in one line from known identities. Session223 proved the mod-3
    limits (34/27, 47/54) that the packet calls conjectures.
  - **The cutoff section.** It contains an "impossibility" paragraph that is not proved and is partly false.
  - **The Fermat boost.** The mechanism the packet asserts for it ("this is why 3, 5 and 17 dominate") is contradicted by
    the auditor's decomposition.
- **Coordinator's closed form: correct and provable.** The closed form L_p(0) = 1 + (3p-2)(p-1)/(2p^3) is correct for
  every odd prime p, with the stated reduction. Lemma (M) is provable for every odd prime (proof in section 5), so the
  addendum can be labelled PROVED. The method and the p = 3 case are Session223's.

## 1. Recomputed numbers

| Item | Packet | Auditor (own code) | Match |
| --- | --- | --- | --- |
| Tail formula H_s(m), all 4,463 states with <= 2 J and length <= 9, plus 400 random integer tuples, m < 2^9 | used | 1,740,800 literal evaluations equal | yes |
| (1a) confinement mod gcd(R,d); gcd(R,d) = gcd(D,U) for J(p); odd tails divisible by gcd(D,U) | proof | 870,400 / 3,400 / 384,000 cases | yes |
| (1b) both parities; gcd(Q(m)+1, m+1) divides odd-tail values | proof | 1,740,800 cases each | yes |
| (1c) value(pJZO^(j-1)Z) = (2^(j-1)+1)(D+2U) | 889 words | 63,219 cases: every state with <= 2 J reachable in <= 9 letters plus 400 random integer tuples, all j = 2..14. It fails for j = 1 on 4,859 of 4,863 states (j >= 2 is needed) | yes |
| (1d) a-dynamics; K expansion; K(t.alpha, t.beta) = F(t.(beta alpha^dagger)) | proof | 4,863 / 32,000 / 43,720 cases (all alpha up to 6 letters) | yes |
| (1e) built literally as the word beta alpha^dagger Z and evaluated letter by letter | 3,000 pairs | 3,000 words, each with value (x+2^p y)(x'+2^q y') and exactly two J. p = 0 fails on 252 of 500 pairs | yes |
| Automaton against literal words | none | weighted class counts identical for q = 27, 5, 7, 11 and n = 0, 1, 2, 5, 8, 9 (24/24; JJ, J first and J last all included) | yes |
| `dp_bands.txt`, `dp27_classes.txt`, `dp27_long.txt`, `dp_primes_n1000.txt`, `q3_long.txt` | files | 512 printed numbers compared. All L(0) and per-class values are within print rounding. The L(nz)avg and ratio columns of `dp_bands.txt` come from 4-decimal per-class values (42/42 rows reproduced that way; at most 8.1e-5 from unrounded). Registered R2 bands move by at most 1e-4 | yes |
| q = 3 at n = 4000 / 8000 / 16000 | 1.2603 / 1.2598 / 1.2595 | 1.2602511 / 1.2597556 / 1.2595075. Successive differences have ratio 1.998 (1/n). Richardson gives 1.2592595 (34/27 = 1.2592593). n and n+1 differ by about 2e-7: no parity oscillation | yes |
| Census at odd N in [2^23, 2^24): class means mod 27, 5, 7, 11, 13, 17, 23, 31 | `lf_2p24.txt` | largest difference 5.0e-5 (print rounding). Overall mean 264.618519, max M_2 = 3792, 2,044,971,045 counted words with value <= 2^24 | yes |
| Mean M_2 by v3(N) = 0..6 | 205.8, 332.0, 457.3, 510.7, 556.8, 614.2, 618.7 | 205.8461, 332.0336, 457.3275, 510.7218, 556.7526, 614.1891, 618.6854 | yes |
| R2 | FAILED | FAILED. rho = 1.61397, 1.23951, 1.37393, 1.17096, 1.53378, 1.12400, 1.10195 for p = 5, 7, 11, 13, 17, 23, 31. p = 13 is outside by only 0.00026; p = 5, 11 and 17 by 0.047, 0.017 and 0.048 | yes |
| R3, R4 | CONFIRMED | CONFIRMED (worst nonzero max/min 1.04518, at p = 7) | yes |
| R5 | 27/27 | 27/27 inside the 5%-widened band; only 9/27 inside the unwidened band | yes (m4) |
| `decomp_2p24.txt` | 14 j3/b lines, TOTAL, p = 5 and 7 lines | all identical. Confined odd-tail mass outside N = 0 (mod 3^a): 0 | yes |
| Unconfined max/min over classes mod 729 | 1.28574 (2^24); 1.87 (2^18) | 1.28574; 1.86983. Permutation-noise null: 1.042 (2^24) and 1.498 (2^18) | numbers yes; reading no (m3) |
| Misses of F_<=2 below 2^24 | 140, largest 3,780,967 | same; the list above 10^5 is identical to `rc_1e9.txt` | yes |
| M_2 <= 12 in (2^22, 2^24] | 90 entries of `rc_low_2p28.txt` | identical values and counts | yes |
| Near-miss arithmetic: 5 of 15 misses prime, 14/15 = 1 (mod 3), none = 0; 3 of 12 (M_2 <= 6) prime, only 13 below 20; 277,438,727 prime | as stated | all exact | yes |
| Pre-registration | 0cda244632 before 05b4c174db | 0cda244632 is an ancestor and on origin/main. `PREREGISTRATION.md`, `resdp.c` and `dp_bands.txt` are unchanged in 05b4c174db; `eval_s15.py`'s bands equal the registered table | yes |
| Square-sum counts (ATLAS) | paths 10, 12, 35, 52, 19, 20, 349, 392 (n = 25..32); cycles 1, 1, 11, 57, 31, 20 (n = 32..37) | identical (own DFS) | yes |
| Mean M_2 at the last miss (ATLAS) | "about 180" | 226.94 on [3.5e6, 4.0e6) (exact 32-bit counts) | no (m8) |
| Closed form (addendum), two J and one J | 1 + (3p-2)(p-1)/(2p^3); 1 + (3p-4)/(4p^2) | exact rational equality at p = 3, 5, 7. Within 1.2e-12 for every prime up to 31 (lazy-chain projector, a different method) | yes |

## 2. MAJOR findings

### M1. Priority: THM-4916 is not "new as stated", and Session223 already proved the 3-adic limits

**Prior sources for each identity.** These are checked exactly in `priority_check.py`.

- **(1a) is KNOWN.** Session206 `weighted-parameter-capacity/WEIGHTED_CAPACITY.md` says "g divides R and d, so
  H_s(m)=h-eb modulo g". See also `exact-modular-image`.
  - The repository already restates it as THM-4888 (C), whose status says "KNOWN (Session 206)".
  - The root-change form is THM-4888's consequence 2, stated mod 3 in `finite_oz_bank_poisson_root_changes_20261106/PROOF.md`.
- **(1b) is Session223 `cutoff/REPORT.md` eq. (2):** F(s.JC) = alpha D + beta E + gamma U.
  - It was committed in 1c524e45a9 at 13:31:03, 32 s before this packet's checkpoint. The two are concurrent and
    independent.
  - Session223's coefficients equal (1b)'s exactly, with h = m and r = c_C + q_C + e_C = Q(m). Checked on all 2,046
    blocks of length 1..10.
  - New in (1b): the m, Q(m) form and the corollary that gcd(Q(m)+1, m+1) divides every odd-tail value.
- **(1c) follows in one line from two known facts:**
  - s.JZ = s.OJ = (D+2U, D, D+2U, D): Session205, THM-4888 (D)1, and Session199
    `hidden-structure/dependencies/BALANCED_SOURCES.md` Theorem 2 (balanced sources are exactly W.OJ);
  - F((h,e,h,e).O^k) = (2^k+1)h: the mechanism of Session223's rectangle F(t.Z^a O^b) = (2^a+1)(2^b+1), committed in
    8865bb666b at 13:20:51, before this packet's pre-registration.
- **(1d) is Session199's coordinate up to sign.**
  - a = D-E-U+d = -f, where f = U-D+E-d is Session199's positive coordinate, and Session199 already gives the O/Z/J laws
    of f.
  - The same expression is Session201 `digital-inverse`'s delta.
  - The K expansion rewrites Session199's (2).
- **(1e) combines two known results.**
  - Session199 `physical-adjoint` section 4: K((h,e,h,e),(h',e',h',e')) = hh'.
  - Applied to the balanced states w J Z^p = w O^p J (Session223 `bridges` section 1).
  - New as stated: only the J-budget form V V ⊆ F_<=2.

**Session223 is described as out of date.** Session223 checkpoint 2 (1c524e45a9) proves the fixed-length 3-adic law at
every precision: 34/81, 47/162, 47/162 mod 3, the one-J law (41/108, 8/27, 35/108), and the valuation factors. It uses
the same Cesaro three-block projector the coordinator's addendum uses, and handles the period-two component. The packet
still describes Session223's 3-adic factor and cutoff as "under investigation" and calls 34/27 and 47/54 "conjectured".

**Replacement text.**

THM-4916 `status:` replaces everything before "Session223 (same owner prompt, concurrent) has...":

> PROVED by direct substitution (audited). Priority:
> - (1a) is KNOWN: Session206 `weighted-parameter-capacity` (g = gcd(R,d,p) divides R and d, so H_s(m) = h - eb mod g)
>   and `exact-modular-image`, restated as THM-4888 (C). Its root-change form gcd(R,d) = gcd(D_p,U_p) generalizes
>   THM-4888's consequence 2 (mod 3).
> - (1b) is Session223 `cutoff/REPORT.md` eq. (2), from commit 1c524e45a9, pushed 32 s before this packet's checkpoint
>   (concurrent and independent), written with h = m and r = Q(m). New here: the m, Q(m) form and the tail-confinement
>   corollary gcd(Q(m)+1, m+1) | N.
> - (1c) is a one-line corollary of two facts: s.JZ = s.OJ = (D+2U, D, D+2U, D) (Session205; THM-4888 (D)1; Session199
>   `hidden-structure` Theorem 2), and F((h,e,h,e).O^k) = (2^k+1)h, the mechanism of Session223's rectangle (commit
>   8865bb666b).
> - (1d): a = -f for Session199's positive coordinate f = U-D+E-d, whose O/Z/J laws are Session199's. The K expansion
>   rewrites Session199's (2).
> - (1e) is Session199 `physical-adjoint` section 4 (K = hh' on balanced pairs), applied to the balanced states
>   w J Z^p = w O^p J. New as stated: the J-budget form V V ⊆ F_<=2.
>
> Credited as before: Session199's K-form and Session200's Theorem 2.

README "Concurrent same-prompt work" and the PROOF preamble: replace "a height-weighted 3-adic factor and the cutoff,
both under investigation" (README) and the sentence "It lists a "height-weighted 3-adic factor" and the cutoff as under
investigation." (PROOF) with:

> Its checkpoint 2 (1c524e45a9, 13:31:03, pushed 32 s before this packet's checkpoint) PROVES the fixed-length 3-adic law
> at every precision, by the Cesaro-projector method of section 2:
> - mod 3: 34/27, 47/54, 47/54;
> - with one J: 41/36, 8/9, 35/36.
>
> It also gives the exact two-J membership equation, whose eq. (2) is (1b). The height-weighted factor and the cutoff
> remain open there.

PROOF section 2 (mod 3 bullet) and HYP-9506 (A): replace "The conjectured exact values are 34/27 = 1.259259 and 47/54 =
0.870370." with:

> These are the values PROVED by Session223 (`threeadic/PROOF.md`). The letter-count limit coincides with its
> fixed-length law.

### M2. The cutoff section's "impossibility" paragraph is not proved, and is partly false

PROOF section 3, item 2 is listed under "What is proved". It says that coverage "by one seed, by any finite bank, or by
any fixed linear form in one free state is impossible".

- **Single seeds: proved, but by the wrong citation.** The fact is PROVED, by Session206 `parameter-capacity` (for every
  seed with log(R/(R-d))/d < 1/2, which is every physical seed except the root) and `root-residue-capacity` (the root).
  It does not follow from HYP-9480 (A) or (C), which are a conjecture and an OPEN item.
- **Finite O/Z banks: OPEN.** This is HYP-9480 (C).
- **Session217 does not apply.** Its O((log N)^2) bound (`paired/audit/approved/COMMUTING.md`) concerns banks of the
  commuting Q/Z two-counter construction. A single O/Z seed already has about cX outputs up to X.
- **"Any fixed linear form in one free state" is false as stated.** F itself is such a form: K(t,s) = F(s) and
  F({s : <= 2 J}) = F_<=2 \ {1}. So HYP-9506 (E), if true, would contradict this sentence.
- **"Every proof needs the joint freedom of two states" is heuristic.**
- **The barrier lemma itself is correct but trivial.** If the non-product part vanishes, then N = U_v U_w with both
  factors at least 2.

**Replacement for PROOF section 3, item 2 (whole item):**

> **Barrier lemma (elementary).** Every product family (1e) and every multiplier family (1c) consists of composite
> numbers. So for a prime N every representation N = K(v,w) has a nonzero non-product part
> E_vE_w - d_vd_w + a_vd_w + a_wd_v; if that part vanished, N = U_vU_w with U_v, U_w >= 2.
>
> Related facts, with their status:
> - a single O/Z seed misses a positive proportion of odd integers (PROVED: Session206 `parameter-capacity` and
>   `root-residue-capacity`);
> - whether every finite O/Z bank does is OPEN (HYP-9480 (C));
> - Session217's O((log N)^2) bound is for banks of the commuting Q/Z construction, not for O/Z tail banks.
>
> One free state is not excluded: K(t,s) = F(s), so F_<=2 \ {1} is the image of a single free state with at most two J.
> Preferring two-state arguments is a heuristic choice.

**HYP-9506 (E):** replace "Single seeds and finite banks fail." with:

> A single O/Z seed misses a positive proportion (Session206). For finite O/Z banks this is OPEN (HYP-9480 (C)).

### M3. The asserted mechanism of the Fermat boost is not supported

**The claims.**
- ATLAS section 3: "the Fermat primes 3, 5, 17, 257 are exactly the primes p for which a single cheap tail makes every
  state's balanced readout a multiple of p. This is why 3, 5 and 17 dominate the local factors at census scale (R3
  CONFIRMED). This is a mechanism, not a coincidence".
- README item 4 says the same, and HYP-9506 (B) repeats it.

**Why "exactly" fails.** Every prime with even ord_p(2) divides some multiplier 2^i + 1: 11 | 33, 13 | 65, 19 | 513,
41 | 1025. The Fermat primes (including 65537) are exactly the primes that are themselves multipliers.

**What the decomposition shows** (`census_mechanism.c`, window [2^23, 2^24)). The zero-class excess
mean(N = 0 mod p) - mean(N != 0 mod p) splits into four categories, first match wins:

- **PREF:** the word begins with the root rectangle Z^aO^bJ (a, b >= 1) and p | 2^a + 1.
  - Then t.Z^aO^bJ = (h, g, h, g), with g = 2^a + 1 dividing h, so every continuation is divisible by 2^a + 1.
  - All 174,482,210 such counted words up to 2^24 are divisible by 2^a + 1; none fails.
  - This Fermat-type family is not mentioned in the packet.
- **MULT:** the prefix ends in a (1c) multiplier block after its last J.
- **SEED:** seed confinement at the last root change, p | gcd(D,U).
- **REST:** everything else.

| p | 3 | 5 | 17 | 11 | 13 | 7 | 23 | 31 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rho_p | 1.857 | 1.614 | 1.534 | 1.374 | 1.171 | 1.240 | 1.124 | 1.102 |
| share PREF | 24.0% | 21.4% | 19.3% | 8.7% | 10.5% | 0 | 0 | 0 |
| share MULT (1c) | 23.7% | 23.2% | 22.0% | 11.5% | 11.7% | 5.4% | 1.4% | -0.2% |
| share SEED | 54.1% | 55.5% | 48.1% | 70.1% | 96.9% | 101.4% | 121.6% | 128.5% |
| rho without PREF and MULT | 1.463 | 1.363 | 1.339 | 1.324 | 1.145 | 1.247 | 1.133 | 1.111 |

**Conclusions.**
- The (1c) multiplier words carry only about a quarter of the boost at 3, 5 and 17.
- Seed confinement carries about half at 3, 5 and 17, and more than all of it at 7, 23 and 31. There the unconfined
  remainder is depleted at 0, as m3 predicts.
- R3's ordering 5 > 17 > 11 > the rest survives the removal of both Fermat-type families. So R3's confirmation does not
  demonstrate the multiplier mechanism.

**Replacement for the second paragraph of ATLAS section 3** ("So the Fermat primes ... measured consequence."):

> The Fermat primes 3, 5, 17, 257, 65537 are exactly the primes that are themselves multipliers 2^i + 1. Every prime
> with even ord_p(2) divides some multiplier (11 | 33, 13 | 65, 19 | 513).
>
> Fermat-type numbers enter the census through two exact families:
> - the (1c) multiplier words;
> - the words that begin with the root rectangle Z^aO^bJ, all of whose values are divisible by 2^a + 1, because
>   t.Z^aO^bJ = (h, g, h, g) with g = 2^a + 1 dividing h.
>
> At [2^23, 2^24) these carry about 23% and 21% of the zero-class excess at p = 5, and 22% and 19% at p = 17. Seed
> confinement at the last root change carries the rest, about half. The R3 ordering persists with both families removed:
> rho = 1.363, 1.339, 1.324, 1.247 for p = 5, 17, 11, 7. So R3 confirms the ordering but not this mechanism.

**README item 4:** replace "Fermat primes enter the root-change grammar as exactly the cheap multipliers of THM-4916
(1c)." with:

> The Fermat primes are the prime multipliers 2^i + 1 of THM-4916 (1c). At census scale, the multiplier words and the
> root-rectangle words together carry under half of the boost at 3, 5 and 17 (audit, M3).

**HYP-9506 (B):** append to the first sentence:

> (they carry about a quarter of the zero-class excess at 3, 5 and 17, and an eighth at 11 and 13; the rest is
> root-rectangle words and seed confinement, and the R3 ordering survives the removal of both Fermat-type families)

**HYP-9506 (B) status:** see m5.

## 3. MINOR findings

1. **m1. (1d)'s set identity misses N = 1.** PROOF (1d) and THM-4916 item 4 say "F_<=2 is the set of values K(v,w) with
   #J(v) + #J(w) <= 2". But K >= U_vU_w/2 + d_vd_w >= 3, while 1 lies in F_<=2 (the empty word).

   Replace with:

   > F_<=2 = {1} ∪ {K(t.alpha, t.beta) : #J(alpha) + #J(beta) <= 2}.

2. **m2. (1a) wording.** "Odd tails land in U (mod gcd(D,U))" hides the fact that U = 0 mod gcd(D,U).

   Replace with:

   > Odd tails are divisible by gcd(D,U) (h - e = U), and even tails are congruent to D - E + U, that is to -E
   > (mod gcd(D,U)).

3. **m3. The decomposition paragraph (PROOF section 2) misreads its own data.**
   - **The 2^18 ratio is mostly noise.** The ratio 1.870 at 2^18 is largely sampling noise: shuffling the window's
     unconfined counts gives 1.498 on average, and 1.573 at worst over 20 runs. At 2^24 the noise floor is 1.042.
   - **A noise-robust statistic** is the class-0-mod-3 ratio of the unconfined part. It moves from 1.126 (even tails)
     and 1.256 (odd tails) at 2^18 to 1.088 and 1.166 at 2^24.
   - **"It comes from tail confinement (1b)" is not right.**
     - Tail confinement alone gives no net bias. If (m - Q(m), m + 1) and (D,U) are uniform, the confined odd tails
       (all at 0) are exactly compensated by the other odd tails. Those are a nonzero linear form that vanishes on
       p - 1 of the p^2 - 1 nonzero pairs.
     - Measured: the not-tail-confined odd tails have class-0 ratio 0.683.
     - The bias comes from an excess of tail-confined tails: 15.3% of the unconfined odd-tail mass at 2^24, against 1/9.
     - Even tails have no tail confinement at all, yet are biased at 1.088.

   Replacement for the second decomposition bullet:

   > Seeds with 3 ∤ gcd still favour N = 0 (mod 3). Their max/min over the classes mod 729 is 1.286 at 2^24, against a
   > sampling-noise floor of 1.04 (1.870 at 2^18 against 1.50, mostly noise). Their class-0-mod-3 ratio falls from 1.126
   > (even tails) and 1.256 (odd tails) at 2^18 to 1.088 and 1.166 at 2^24. Tail confinement alone has no net effect
   > for uniform residues. At census scale the odd-tail bias comes from an over-representation of tail-confined tails
   > (15.3% against 1/9), and the even tails are biased by an unidentified finite-size effect. In the letter-count limit
   > the whole zero-class excess is seed confinement at the last root change (closed form, section 5 of the audit).

   In the "Mechanism" paragraph, replace "Two kinds of confinement make the classes N = 0 (mod p) heavy" with:

   > Seed confinement makes the classes N = 0 (mod p) heavy. Tail confinement adds to it only at finite length.

4. **m4. R5 and R1 bookkeeping.**
   - **R5.** R5 passes 27/27 only because of the registered 5% widening; 9/27 classes lie inside the unwidened band. At
     every class the census is more contrasted than the automaton, in the same direction as R2's failure. Add this
     sentence to the R5 row of PROOF section 2.
   - **R1.** `census_2p32.txt` is listed in the README's Files table but is absent, and PROOF section 4 is a
     placeholder. Mark both "(pending)".

5. **m5. HYP-9506 (B) status.** "(B) CONFIRMED by pre-registered R3" overstates what was tested: R3 tests the ordering at
   one scale only. The "1/log N" decay is untested; rho_5 is 1.667 on [2^19, 2^20) and 1.614 on [2^23, 2^24).

   Replace with:

   > (B) ordering CONFIRMED at one scale (R3); the decay rate is untested; the attribution to (1c) is partial (audit
   > M3).

6. **m6. The n = 1000 table is not asymptotic.** PROOF section 2's "asymptotic L(0) at n = 1000" values carry an O(1/n)
   bias: n(L_n - L_inf) tends to 4.48, 2.60, 5.46 and 7.94 for p = 5, 7, 11 and 17.

   - Replace the table with the exact limits L_p(0) = 1 + (3p-2)(p-1)/(2p^3): 34/27, 151/125, 400/343, 1486/1331,
     2419/2197, 5305/4913, 7354/6859, 12904/12167, ...
   - Replace "(L(0) - 1) p tends to about 1.5" with:

     > (L(0) - 1)p = (3p-2)(p-1)/(2p^2), which tends to 3/2 (it is 1.420 at p = 31, not 1.485).

   - Replace "as expected from the product structure (1e)" with:

     > The excess is seed confinement at the last root change (the c-term). Product words are a vanishing fraction of
     > all words.

   PROOF section 2 also says "It is this automaton's finite-n pattern, with an asymptotic mod-3 limit near (34/27, 47/54,
   47/54)". The census factor at N = 0 (mod 3) is 1.444 at [2^23, 2^24). Whether the height-weighted factor tends to
   34/27 is OPEN; Session223 section 2 shows that residues correlate with real-state history.

   Replace with:

   > Qualitatively this is the automaton's finite-n pattern. The limit of the height-weighted census factor is OPEN.

7. **m7. Near-misses.** The PROOF's literal counts are correct, but "rough" is uneven.
   - **The 130 values with M_2 <= 12 in (2^22, 2^28) are strongly rough.** 21 of 130 (16%) have a prime factor in
     5..19, against 49% for random odd numbers prime to 3. None is divisible by 5 or 17.
   - **13 is not depleted:** 11 of the 130 are divisible by 13, against 10.0 expected.
   - **The 15 misses above 10^5 are only weakly rough.** 6 of them have a prime factor of at most 19, for example
     116,215 = 5·11·2113 and 797,083 = 7^2·16267.

   Replacements:
   - HYP-9506 (D): replace "near-misses ... are rough (no factor 3, few small factors)" with:

     > near-misses avoid 3 and are depleted in 5, 7, 11, 17 and 19, but not in 13

   - ATLAS section 2: "Every large miss is rough and 1 (mod 3)" is false; 136,607 is 2 (mod 3). Replace with:

     > 14 of the 15 misses above 10^5 are 1 (mod 3), none is 0, and 6 have a prime factor below 20.

8. **m8. ATLAS growth rates.** "M_2 mean about 0.5 ln^2 N" (two rows) is wrong by a factor of about 2. HYP-9480 (D)'s
   0.500 is the coefficient of the cumulative mass 0.500 X ln^2 X, so the mean per odd N is about
   ln^2 N - 2.06 ln N + 4.2. The exact census gives 222.0 on [2^21, 2^22) and 264.6 on [2^23, 2^24).
   - Replace "M_2 mean about 0.5 ln^2 N" with:

     > mean M_2 about ln^2 N

   - Replace "mu(N) is about 180 at the last miss" with:

     > mu(N) is about 227 at the last miss (exact counts on [3.5e6, 4e6); the 180.8 in `rc_1e9.txt` is an 8-bit
     > saturated lower bound)

   - Change e^-144 to e^-182. The conclusion is unchanged.

9. **m9. ATLAS Chan-Pak row.** "21 = 3 * 7 or falls in (15, 25)" omits the essential exhaustive step: strong tournaments
   on 5 or 6 vertices are not excluded by the floor (f(5) = 9 and f(6) = 15 are both below 21).

   Replace with:

   > 21 needs a strong component with H = 21, which would need n <= 6 because f(7) = 25. THM-1370's exhaustive spectra
   > show that none exists, and 21 = 3·7 would need a strong 7.

   Further corrections:
   - Cite the file `THM-1370-h-spectrum-omits-7-21-all-n.md`, because two canon files carry the ID THM-1370.
   - The 2^32 - 1 coverage is S78's THM-4875 (D), recertified by Session203. Change "(Session203)" to:

     > (S78, THM-4875 (D); recertified by Session203)

10. **m10. ATLAS section 3, the Session219 totient law.** Replace "a + 2 for a <= 31, then exactly 32" with:

    > a + 2 for 0 <= a <= 31, and 32 for 32 <= a <= 63 (Session219; beyond 63 it depends on further Fermat numbers)

11. **m11. ATLAS row for HYP-9451.** Replace "'eventually always' fails for a long time" with:

    > arrangements need Chebyshev's mod-4 bias in {0, 1, 2} and six mod-36 conditions, which below 5·10^6 hold on only
    > two windows; whether infinitely many bounds admit one is OPEN

12. **m12. ATLAS section 4, Session211.** "Session211's port bridge sends square-sum path counts to tournaments by
    H = 2p + 1" inverts the result. Replace with:

    > Session211 proves that no map from square-sum graphs to tournaments can have H = 2p(Q_n) + 1: p(Q_23) = 3 and
    > p(Q_25) = 10 would force H = 7 and 21.

13. **m13. Reproducibility.**
    - The "889 words" check of (1c) and the "3,000 random pairs" check of (1e) have no script in the packet. Cite
      `audit/identities_audit.py` or add one.
    - The L(nz)avg and ratio columns of `dp_bands.txt` come from 4-decimal per-class values. Say so. It changes no
      outcome.
    - `dp27_long.txt` and `dp_primes_n1000.txt` are truncated excerpts, not raw outputs. Say so.

14. **m14. The new files `exact_limit.py` and `exact_limit.out`.**
    - **The docstring is wrong about periods.** "averaging over 12 consecutive steps, which kills periods up to 12" is
      false: averaging kills only periods that divide 12. The d = 0 classes have period ord_p(2), which is 10, 8, 18
      and 11 at p = 11, 17, 19 and 23.
    - **The results are nevertheless correct.** The U-phase affects neither the F-law nor the confinement event, and
      the auditor's lazy-chain projector lim((I+P)/2)^k agrees to 1.2e-12.
    - **Some displayed fractions are wrong.** The "~ a/b" fractions come from `limit_denominator(10000)` and are wrong
      whenever the denominator exceeds 10^4:
      - p = 19, nonzero classes: the true value is 13663/13718, not 7701/7732;
      - p = 23: the true values are 12904/12167 and 24267/24334, not 6951/6554 and 9417/9443.
    - **Fix:** use the closed form, or raise the bound to 10^6.

## 4. Task-specific answers

1. **(1a)-(1e) are correct.** The proofs are one-line substitutions and were re-derived independently. The exact checks
   are in the table above.
   - (1c) holds for every state, including states with one or two J, but only for j >= 2.
   - (1e) needs p, q >= 1. Its realizing word has exactly two J.
2. **The automaton is correct.**
   - **State space.** (Z/q)^4 is correct.
   - **J-closure.** The J-closure (A_1 += J A_0, then A_2 += J A_1 at each stage) allows a first J, a last J and JJ.
     Literal word enumeration agrees exactly (24/24).
   - **Reproduction.** All of `dp_bands.txt`, `dp27_*.txt`, `dp_primes_n1000.txt` and `q3_long.txt` reproduce.
   - **Limits.** (L_p(0) - 1)p tends to exactly 3/2, and the mod-3 limit is exactly (34/27, 47/54, 47/54) (section 5).
3. **The registered census reproduces.** It matches exactly, and R2-R5 come out as the packet says when the registration
   is applied as written: the registered band table for R2, and the multiplicative 5% widening for R5, the convention the
   R2 table demonstrably uses.
   - **Ordering.** The pre-registration commit precedes the outcome commit and was pushed to origin/main.
   - **Not verifiable.** That no registered census value was computed before 13:22 cannot be checked from the files.
   - **R5's per-class band.** It was not in the pre-registration commit, but it is fixed by the committed `resdp.c`; the
     auditor's automaton reproduces `dp27_classes.txt`.
4. **The decomposition reproduces.** All confined odd-tail mass lies on N = 0 (mod 3^j), and both ratios reproduce. The
   reading needs m3.
5. **The near-miss statements reproduce.** Every literal statement is correct. "Rough" needs m7. 277,438,727 is prime.
6. **The cutoff section and ATLAS.** The barrier lemma is correctly stated but trivial. The impossibility paragraph
   overclaims (M2).
   - **Correct ATLAS facts:** THM-1370's floors 3, 5, 9, 15, 25, 45; the coverage to 2^32 - 1; the square-sum facts
     (Gerbicz base 41..2032 with the step 49n + res, A071983 and A071984 counts, terminal holes at 33 and 39 with 50
     cycles at 39); Boklan-Conway at under 10^-9 (about 2^-30); HYP-9410 at 6.6·10^-9; THM-4847; THM-4884; THM-4805;
     Session216's {47, 63}; Session219's five-run escape.
   - **ATLAS misstatements:** m7-m12 and M3.
7. **Priority:** see M1. KNOWN: (1a). Concurrent: (1b), Session223 eq. (2). Corollaries of known identities: (1c), (1d)
   and (1e). New: the tail-confinement corollary of (1b), the general-p identity (1c) as a formula, and the J-budget form
   of (1e).

## 5. Coordinator addenda: the n -> infinity limit, lemma (M) and the closed form

**Verdict: CONFIRMED, and upgradable to PROVED for every odd prime p** (letter-count model). The algebra of the
three-block derivation is correct line by line:

- c = (7p-4)/(4p^3);
- e = (5p-2)/(4p^3);
- L(0) = 1 + c(p/2 - 1) + (p/2)e = 1 + (3p-2)(p-1)/(2p^3);
- L(r) = 1 - (3p-2)/(2p^3) for r != 0.

With one J, L(0) = 1 + (3p-4)/(4p^2), L(-1) = 1 - 1/p^2 + 1/(4p), and the other classes have 1 - 1/p^2.

**Checks with the auditor's own code.**

- **Exact rational arithmetic** (`cesaro_exact.py`): a generic SCC decomposition with exact stationary laws and
  absorption probabilities, no structural assumption.
  - The two-J and one-J limits equal the closed forms exactly at p = 3, 5, 7: 34/27, 151/125, 400/343; 41/36, 111/100,
    213/196.
  - (M1)-(M3) hold exactly on the class structure. For d != 0, the 2p^2 states with E in {0, d} form one class, which
    has a self-loop and a uniform stationary law. For d = 0, the classes are D × (doubling orbit of U) with a uniform
    law, plus fixed points (x, 0, 0, 0).
- **A different numerical projector** (`cesaro_lazy.c`), lim((I+P)/2)^k: every prime p <= 31 agrees with the closed forms
  within 1.2e-12, for both two J and one J. The values at q = 9 and 27 equal Session223's f(v) and C_k.
- **Extrapolation of the auditor's automaton.**
  - At p = 5, 7, 11, 17, n(L_n - L_inf) converges, to 4.48, 2.60, 5.46 and 7.94, so the convergence is O(1/n).
  - Richardson extrapolation reaches the closed form within 1e-6.
  - At p = 17, where ord_17(2) = 8, nine consecutive n show no oscillation.
- **q = 15.** The limit is not the product of the mod-3 and mod-5 factors. The ratios are 1.00877 (class 0), 0.99721
  (0 mod 3 only), 0.99366 (0 mod 5 only) and 1.00202 (units), which confirms "0.994 to 1.009".

**Proof of (M)** (p an odd prime; (M1) holds for any odd modulus with d a unit).

- **(M1)**
  - After one letter, a state with d != 0 lies in S_d = {E in {0, d}}, which is closed.
  - Each state of S_d has exactly two incoming letter-edges:
    - O-preimages of (D, 0, U): (D, 0, U/2) and (D, d, U/2);
    - Z-preimages of (D, d, U): U_0 = (U+d)/2 and D_0 = D + E_0 - U_0, for each E_0 in {0, d}.
  - So the chain is doubly stochastic, the uniform law is stationary, and every state is recurrent.
  - Every state reaches Phi = {(x, d, d, d)}. Steer U to d with 2^k >= q letters (2^kU - d·sum z_i 2^(k-i) runs over all
    residues), then apply one Z if the last letter was O.
  - Z fixes every point of Phi (a self-loop, hence aperiodicity), and O^r Z with r = ord_q(2) maps (x, d, d, d) to
    (x+d, d, d, d).
  - So S_d is a single aperiodic recurrent class with uniform law: (D, U) is uniform and independent of E in {0, d}.
- **(M2)**
  - After the first letter E = 0 and U_k = 2^k U_1 != 0.
  - Then D_(k+1) = D_k + z_k U_k.
  - For xi != 0, each step multiplies the characteristic function by (1 + e(xi U_k/p))/2, of modulus at most
    cos(pi/p) < 1.
  - So D_k becomes uniform, independently of the U-phase, and F = D + U is uniform.
- **(M3)**
  - O sends (D, E, 0, 0) to (D, 0, 0, 0) and Z sends it to (D-E, 0, 0, 0); both letters fix (x, 0, 0, 0).
  - Right after a J, a = 0 gives D = E.

**The reduction is justified.** This is Session223 section 5's argument, valid for every odd q.

- **Weights.** Prefixes with at most one J have total weight O(n·2^n), against C(n+2,2)·2^n for exactly two J.
- **Compositions.** Exactly-two-J prefixes are uniform over the compositions (a, b, c) of n.
- **Spectral split.** P^k = Pi + sum over |lambda| = 1, lambda != 1 of lambda^k Pi_lambda + O(rho^k), with semisimple
  peripheral eigenvalues. The first block starts at t, in the d = 1 class, which is irreducible and aperiodic by (M1),
  so P^a delta_t = pi_1 + O(rho^a).
- **Averaging.** Averaging lambda^b mu^c with (lambda, mu) != (1, 1) over b + c <= n gives O(n) against about n^2/2.
- **Conclusion.** L_n = (readout of Pi J Pi J Pi delta_t) + O(1/n).
- **Periodic components.** These are the d = 0 classes, with period ord_p(2). They are harmless because the first block
  is aperiodic, and because their phase affects neither the F-law nor the confinement event.

**Suggested label:**

> PROVED in the letter-count model for every odd prime p. The method and the case p = 3 are Session223's
> (`threeadic/PROOF.md`, `ALL_K.md`, commit 1c524e45a9). New here: the closed forms for general p, two J and one J.

The census (height-weighted) factor remains OPEN (R2 FAILED). The fix list for the new files is m14.

## 6. Files (this folder)

| File | Purpose |
| --- | --- |
| `identities_audit.py`, `identities_audit.out` | exact checks of (1a)-(1e), the tail formula and the incremental Q rule |
| `priority_check.py`, `priority_check.out` | (1b) equals Session223 eq. (2); (1c), (1e) from JZ = OJ, balanced products, K = hh'; a = -f |
| `resdp_audit.c`, `resdp_brute.py`, `resdp_brute_compare.txt` | own letter-count automaton, checked against literal words |
| `dp_bands_audit.txt`, `dp27_long_audit.txt`, `dp_n1000_audit.txt`, `q3_long_audit.txt`, `dp_convergence_audit.txt` | raw automaton outputs |
| `dp_compare.py`, `dp_compare.out` | comparison with the packet's automaton files |
| `census_lf_audit.c`, `census_2p24.txt`, `census_2p18.txt`, `census_local_3p78e6.txt`, `decomp_compare.out` | own enumerator: class means, misses, low values, decomposition, noise null, tail-confinement split |
| `eval_audit.py`, `eval_audit.out` | R2-R5 from the auditor's census, applying the registration as written |
| `census_mechanism.c`, `census_mechanism_2p24.txt`, `census_mechanism_2p20.txt` | M3 attribution (PREF / MULT / SEED / REST) |
| `cesaro_exact.py`, `cesaro_exact.out`, `cesaro_lazy.c`, `cesaro_lazy.out` | exact and lazy-chain limits; (M1)-(M3) |
| `nearmiss_audit.py`, `nearmiss_audit.out` | near-miss factorizations and small-prime statistics |
| `squaresum_counts.c`, `squaresum_counts.out` | ATLAS square-sum counts |

**Reproduce:**

    gcc -O2 -fopenmp -o census_lf_audit census_lf_audit.c && ./census_lf_audit 16777216 8388608 16777216 4 12
    gcc -O2 -fopenmp -o census_mechanism census_mechanism.c && ./census_mechanism 16777216 8388608 16777216 4
    gcc -O2 -o resdp_audit resdp_audit.c && ./resdp_audit 27 23 18,19,20,21,22,23
    gcc -O2 -o cesaro_lazy cesaro_lazy.c -lm && ./cesaro_lazy 31 2000
    python3 cesaro_exact.py 3 5 7
    python3 identities_audit.py
    python3 priority_check.py

The 2^24 census takes about 72 s on 4 threads.
