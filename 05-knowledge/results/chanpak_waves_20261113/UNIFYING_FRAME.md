# What our tournament mathematics represents: Chan–Pak as a local–global principle for a positive thin semigroup

Session mac-mini-2026-11-15-frame (wave 39 of this packet). **No open problem is solved here. Chan–Pak remains OPEN.**

The owner asked what our new mathematics fundamentally represents, as a deeper unifying theorem. The suggested
connections were:
- Hasse and Artin on the order of 2 mod p;
- Ramanujan's form x² + y² + 10z²;
- Tunnell's theorem and BSD;
- Friedlander–Iwaniec;
- the Hadamard conjecture;
- class numbers;
- Furstenberg's ×2×3 problem.

This file proposes one frame that places all of them, says which parts are proved, and computes the frame's
parameters for span-3 ribbons. Labels: PROVED, KNOWN (cited), FINITE-EXACT, NUMERICAL, CONJECTURE, ANALOGY.

## 0. The answer in one paragraph

Every count we construct is a **matrix coefficient of a semigroup of nonnegative integer matrices**: H(W) = h·A_W·v₀.
- For ribbons, A_W is the product of the letter transfers (THM-4880).
- For general compositions, the semigroup is replaced by multilinear operators on port vectors (THM-5242 (A)).

Chan–Pak says that the value set of this **positive thin semigroup** contains every odd number except 7 and 21.
- That is a **local–global principle with a finite exceptional set**.
- It has the same shape as Zaremba's conjecture (denominators of bounded continued fractions, also a positive
  semigroup; Bourgain–Kontorovich proved density one) and the Apollonian local–global problem.
- In the classical quadratic-form world, the same shape appears in Ramanujan's ternary form and in Tunnell's theorem.

Our results of waves 28–38 are, one by one, the standard components of that machinery:
- local obstructions;
- a multiplicative "torus" part that must be escaped;
- a growth gap for escaping it;
- bilinear structure;
- equidistribution;
- archimedean placement;
- an exact last step.

The seven suggested connections are the classical models of those components.

## 1. The frame (definitions)

- **Semigroup and values.** Fix the span w. The letters are ℓ ⊆ {1, …, w}, and A_ℓ ∈ M_d(ℤ_{≥0}) is the profile
  transfer (d = 119 raw states for w = 3; the minimal representation after the first vertex has dimension 61). The
  semigroup is Γ_w = ⟨A_ℓ⟩, and the values are V_w = {h A_W v₀ : W a word}.
- **Admissibility.** N is admissible if, for every modulus q, it lies in the image of V_w mod q (equivalently, in the
  closure of V_w in the profinite integers).
- **Local–global principle.** Every sufficiently large admissible N lies in V_w, and the exceptions form an explicit
  finite set. For w = 3, HYP-9473 conjectures the exceptions are exactly {7, 21}, which would prove Chan–Pak.
- **Density-one form (CONJECTURE, new here).** #{odd N ≤ X : N ∉ V_3} = O(X^(1−η)) for some η > 0. This is the form
  Bourgain–Kontorovich proved for Zaremba. The ingredients it needs are listed in §4.

## 2. Our results as components of the frame

| component of a local–global proof | what we have | status |
|---|---|---|
| local conditions at 2 | H odd (Rédei); Q ≡ 0 (mod 4/8) and the 12/24 fibre steps (pair-line census §4); the index-2^16 coset of span-2 states (THM-4881) | PROVED |
| local conditions at 3 (and 5) | Q ≡ (−1)^n H (mod 3) (THM-4756); major arcs at 3^i and 5^j in the span-2 local limit law (HYP-9474) | PROVED / CONJECTURE |
| a multiplicative "torus" part | rank-one bottlenecks factor counts exactly, H = L·R (THM-5202). The cheap ones carry the floor multipliers 3 and 5. Rank-one words in general carry many eigenvalues (11, 13, …, 43 at length 6, mostly from cuts), so the torus is not ⟨3, 5⟩ beyond the cheap part | PROVED / FINITE-EXACT (wave 39) |
| a growth gap for escaping the torus | bottleneck-free words grow at least 1.7699 per vertex over the 13 floor units (THM-5212 (E)), against the floor 5^(1/3) = 1.7100; literal span 3 open (wave 39 §3 of the README) | PROVED for units; OPEN for letters |
| bilinear (Type II) structure | concatenation is a Hankel pairing, H(W₁W₂) = ⟨h A_(W₂), A_(W₁) v₀⟩; cyclic joins Re(β_A β̄_B) in ℤ[√−2]; the triangle law (THM-5242) | PROVED |
| equidistribution (spectral gap) | H mod ℓ over all span-3 words of length n: TV 0.114 → 0.020 (ℓ = 5) from n = 6 to 9, roughly halving per letter; the two-apex family at TV 0.005–0.013 (HYP-9792 (ii)); minor-arc decay for span 2 (HYP-9474 (2)) | NUMERICAL / CONJECTURE |
| counting exponent (more words than values) | **CORRECTED 2026-11-15** (audit: chanpak_target_ribbon_audit_20261116/DENSITY_AUDIT.md): 3.32260 is only a lower bound for the joint spectral radius, and a lower bound for α needs an upper bound. The certified bound there is H ≤ 5·(39/10)^n, so α ≥ log 8 / log 3.9 = 1.528; the value 1.732 is what one gets IF the joint spectral radius equals 3.3226 (not proved). A word count is also not a count of distinct values. Original text: span 3: 8^n words, joint spectral radius 3.32260 (the all-back-arcs letter), so α ≥ log 8 / log 3.3226 = 1.732. The typical growth is 2.2449 per letter (log 8 / log 2.2449 = 2.57) | NUMERICAL (exact JSR lower bound) |
| archimedean placement | two-apex Bernoulli rule; placement within O(N^(3/8)) (target-address packet) | PROVED |
| the exact last step | thin readers and fibre swaps (waves 37–38); exact repairs of residuals; seams and switches (concurrent packets) | partial |
| the exceptional set | 7 and 21, explained by odd-cycle crowding below Busch's threshold (THM-5192) | PROVED that they are skipped |

## 3. The seven connections, each as one component

1. **Hasse and Artin: local images of cyclic and torus sub-semigroups.**
   - Strong approximation fails for tori. The image of ⟨a⟩ mod p is a proper subgroup for a positive proportion of
     primes, and Artin's conjecture asks how often it is everything.
   - In our frame the torus part is the bottleneck semigroup. Its cheap multipliers 3 and 5 make ⟨3, 5⟩ mod p the local
     image of near-floor counts. This is why the ribbon local limit law has major arcs exactly at denominators 3^i and
     5^j (HYP-9474).
   - The same phenomenon at the prime 2 is the Collatz-table alternation: −1 ∈ ⟨2⟩ mod p, with density 17/24 by Hasse
     (THM-5242 (C)).
   - *Prediction.* The major-arc masses at 3^i 5^j are governed by the orders of 3 and 5 mod those prime powers, as
     Artin-type local densities.
2. **Furstenberg ×2×3: equidistribution of multiplicative orbits.**
   - ⟨3, 5⟩ is non-lacunary (log 3/log 5 is irrational), so its orbits are dense on the circle (Furstenberg). Mod p they
     are large unless p is special (Bourgain–Lindenstrauss–Michel–Venkatesh).
   - The bulk of the semigroup mixes much faster. That is what the measured geometric TV decay of H mod ℓ says.
   - So "×3×5 inside the floor" plus "expansion outside it" is the tournament version of the ×a×b picture.
3. **Friedlander–Iwaniec: the bilinear form.**
   - Primes in thin sets need Type II sums. Concatenation gives them for free: H(W₁W₂) is a bilinear pairing of a
     left-state vector and a right-state vector, and cyclic joins are Re(β_A β̄_B).
   - The norm (diagonal) route is density-zero (HYP-9792 (i)), so primes must come from genuinely bilinear pairs, as in
     X² + Y⁴ and x³ + 2y³.
4. **Ramanujan's x² + y² + 10z²: a local–global principle with finitely many exceptions.**
   - Ono and Soundararajan proved, under GRH, that this form represents every odd number except 18 listed. GRH enters
     only to exclude Siegel zeros, that is, exceptional characters.
   - Chan–Pak has the same shape. The tournament analogue of an exceptional character is a major arc: the torus at 3
     and 5.
   - That 7 and 21 are on both lists is a coincidence (no claim). The mechanism differs: odd-cycle crowding below
     Busch's threshold.
5. **Tunnell and BSD: main term versus cuspidal error.**
   - Tunnell decides congruent numbers by comparing representation counts of x² + 2y² + 8z² and x² + 2y² + 32z². The
     comparison is a statement about the cusp-form part of theta series; BSD supplies the converse.
   - In the frame, the representation number r_n(N) = #{|W| = n : H(W) = N} splits the same way:
     - main term: archimedean density × local densities (HYP-9474 (1));
     - error term: minor arcs, the "cuspidal" part. Bounding it is exactly the equidistribution component.
   - The binary core x² + 2y² of Tunnell's forms is our diagonal-join norm (THM-5172 (C), THM-5222).
6. **Class numbers: the mass formula.**
   - Siegel's mass formula says Σ over the classes of a genus of 1/|Aut L| = Π_p (local densities) × (archimedean
     factor). Class numbers measure the gap between "locally represented" and "represented by this class".
   - For tournaments the same identity is exact: #{labelled n-tournaments with H = N}/n! = Σ_(classes, H = N) 1/|Aut T|.
     The local limit law of HYP-9474 is then a **Siegel–Weil formula for the mass of the fibre H = N**: archimedean
     log-normal density × local factors at 2, 3, 5.
   - THM-4666 already showed class numbers inside tournaments: the ranking excess of a Paley tournament equals h(−p).
7. **Hadamard: the top of the spectrum and optimal pseudorandomness.**
   - Paley and doubly regular tournaments (skew Hadamard matrices) maximize H at small orders and have flat spectra
     (THM-162).
   - In the frame:
     - the floor 5^(n/3) is the joint spectral subradius of the torus part;
     - the top is the joint spectral radius (3.3226 per letter in span 3, from the all-back-arcs letter, a local
       Paley-like block);
     - Chan–Pak says everything in between is filled.
   - Spectral flatness is the extremal form of the expansion that the equidistribution component needs.

## 4. What a proof in this frame needs, and where each piece stands

1. **(Local) the admissible set.**
   - Conditions at 2 and 3: PROVED.
   - For ℓ ≥ 5, admissibility of every residue: FINITE-EXACT here. Among all words of length 9, no residue mod 5, 7, 11
     or 13 is missed, since each TV is below 1/(2ℓ). A proof needs the image of Γ₃ mod ℓ, a strong-approximation
     statement for the semisimple part of the bottleneck-free semigroup.
2. **(Growth gap) the torus is expensive.**
   - PROVED for floor units at 1.7699.
   - For literal letters the certificate is open, with two obstructions now identified:
     - slow bottlenecks need long or rank-adaptive memory (concurrent LITERAL_CONES.md; this wave's targeted dictionary
       handles periods ≤ 6);
     - letter-level cones need a bounded-growth comparison lemma (wave 39).
3. **(Expansion) minor-arc decay** for the bottleneck-free part, uniformly in the modulus. This is the
   Bourgain–Gamburd–Sarnak input. CONJECTURE (HYP-9474 (2)); the measured TV decay is geometric.
4. **(Bilinear) a Type II estimate** for the Hankel pairing over two independent word sets: the Bourgain–Kontorovich
   circle-method step. Together with 1–3 and α > 1 it should give density one. CONJECTURE.
5. **(Exact) from density one to every N** above an explicit bound: the repair and placement programme (target-address
   packet, waves 37–38), then a finite computation below the bound. OPEN.

## 5. Why this is the right frame (and not only an analogy)

- **Positivity explains every no-go.** The semigroup preserves a cone.
  - There are no cancellations, so no Pell supply (wave 35).
  - There is no unipotent counter or linear tape (HYP-9762).
  - Readers are monotone, and module swaps inside hosts are coarse (wave 38).
  - Zaremba's semigroup is positive too, and its local–global theory works because of positivity: proximality and
    contraction on the projective cone.
- **The bottleneck theorem is the torus/semisimple decomposition.** Rank-one elements are the "parabolic/torus"
  boundary of the projective action. THM-5202 says prime counts must avoid them. THM-5212 says avoiding them costs a
  definite growth premium. A Zaremba-type proof works in the bulk and treats the torus by its local densities.
- **The numbers fit.** Bourgain–Kontorovich-type density-one theorems need the counting exponent to exceed 1 with room.
  For Zaremba, Huang's threshold is a dimension δ > 5/6, i.e. counting exponent 2δ > 5/3 = 1.667. Span-3 ribbons have
  α ≥ 1.528 rigorously (CORRECTED 2026-11-15, see the table above; the earlier figure 1.732 assumed that the joint
  spectral radius is 3.3226). So the comparison with Huang's threshold is not established.

## 6. Receipts (wave 39)

| quantity | value | script |
|---|---|---|
| joint spectral radius lower bound (words ≤ 6 letters) | 3.32260 at the letter {1,2,3} | `scripts/w39_frame_params.py` |
| Lyapunov growth of H, uniform random words | 2.2449 per letter (sd of log H / n: 0.0198 at n = 400) | same |
| H mod ℓ, all words of length 6, 7, 8, 9 (TV to uniform) | ℓ = 5: 0.114, 0.047, 0.034, 0.020; ℓ = 7: 0.264, 0.085, 0.043, 0.023; ℓ = 11: 0.203, 0.128, 0.068, 0.034; ℓ = 13: 0.208, 0.082, 0.039, 0.028 | same |
| eigenvalues of minimal rank-one words (length ≤ 6) | 1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45 (4,928 minimal words) | `scripts/w39_torus_census.py` |
| Ramanujan exceptions | 3, 7, 21, 31, 33, 43, 67, 79, 87, 133, 217, 219, 223, 253, 307, 391, 679, 2719 | `scripts/w38_bilinear_cover.py` |
