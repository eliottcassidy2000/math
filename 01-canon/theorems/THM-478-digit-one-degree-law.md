# THM-478 — The Digit-1 Degree Law: H mod 4 has tile-ANF degree ≤ 2⌊(n-1)/2⌋, and H mod 4 is annihilated by (D+1)-flats

**Status:** Upper bound PROVED and ADVERSARIALLY VERIFIED (2026-06-11: pairing
identity symbolic at L=3,5; reversal shown fixed-point-free on directed cycles;
flat-annihilation equivalence verified exhaustively over all 2^16 functions at
m=4, D=1; all digit ladders reproduced with independent conventions; pointwise
digit_1 = c_odd mod 2 at n=4..6; Grinberg-Stanley arXiv:2307.05569 Thm 7.1 cited
verbatim — NOTE it is a v1 draft, so the congruence itself was ALSO re-verified
exhaustively over all tournaments n=4..6, removing the dependence. See
05-knowledge/results/thm478_adversarial_recheck_macmini_0611.out, 67 PASS 0 FAIL.)
Equality for all n is HYP-2406.
**Provenance:** mac-mini-2026-06-11-S1 (Reed–Muller dispatch). Builds on:
Grinberg–Stanley 2023 (H ≡ 1 + 2·c_odd mod 4, c_odd = # nontrivial odd directed
cycles), THM-163 (band limit D = 2⌊(n-1)/2⌋ for H in the arc model), THM-466
(kind-pasteur: 2-adic digits of H are odd-cycle collection counts), T779.

## Statement

View H as a function on the tiling cube F_2^m (m = C(n-1,2); fixed base path,
tile bits = non-base arcs). Write digit_k for the k-th binary digit of H, a
Boolean function on F_2^m, and let D = 2⌊(n-1)/2⌋.

1. **(Rédei = repetition.)** digit_0 ≡ 1: the constant function, RM(0,m) — the
   repetition-code end of the Reed–Muller ladder.
2. **(Degree law.)** digit_1 (= c_odd mod 2, by Grinberg–Stanley) has GF(2) ANF
   degree ≤ D, i.e. digit_1 ∈ RM(D, m). Equality holds at n = 4, 5, 6, 7
   (degrees 2, 4, 4, 6) — conjecturally all n (HYP-2406). Notably this equals the
   REAL Walsh degree of the integer function H itself (verified 2, 4, 4, 6 in the
   tiling model at n = 4..7): the mod-4 shadow is as high-degree as the full
   function, no higher.
3. **(Dual/flat annihilation — the RM(r)^⊥ = RM(m-r-1) statement.)** For every
   (D+1)-dimensional affine flat F ⊆ F_2^m:  ⊕_{x ∈ F} digit_1(x) = 0.
   (Indicators of (D+1)-flats are the minimum-weight codewords of RM(m-D-1,m) =
   RM(D,m)^⊥, and they generate it, so this family of XOR-recurrences is
   EQUIVALENT to the degree bound.) Equivalently: all (D+1)-fold finite
   differences of H mod 4 vanish — H mod 4 on any subcube of dimension D+1 is
   determined by any 2^{D+1} - 1 of its corners.
4. **(Saturation above.)** The structured zone is shallow: digit_2 has degrees
   3, 6, 7, 11 at n = 4..7 and digits k ≥ 3 saturate to ≥ m-1 (full or near-full
   degree). The 2-adic depth filtration of H enters the Reed–Muller hierarchy at
   the repetition end, climbs through RM(D,m) at depth 1, and exits to
   pseudo-randomness by depth ≈ 3. All digits are σ-invariant (grid transpose).

## Proof of the upper bound (2)

By Grinberg–Stanley, digit_1 = c_odd mod 2 = Σ_{odd ℓ ≥ 3} c_ℓ mod 2. Fix an odd
ℓ-cycle C (a vertex ℓ-subset with a cyclic order, up to rotation). Its two
directed orientations C→ and C← have indicator product Π_{e ∈ C} z_e and
Π_{e ∈ C} (1 + z_e), where z_e ∈ {0,1} are the arc variables (affine in the tile
variables: non-base arcs ARE tile bits, base arcs are constants). Mod 2,

   Π z_e + Π (1 + z_e)  =  Σ_{S ⊊ C} Π_{e ∈ S} z_e   (the two top terms cancel),

a polynomial of degree ≤ ℓ - 1. Summing over all cyclic orders of all odd
subsets: c_odd mod 2 has degree ≤ max odd ℓ ≤ n minus 1 = 2⌊(n-1)/2⌋ = D
(ℓ_max = n for odd n, n-1 for even n; substituting base-arc constants only
lowers degree). ∎

(3) is the standard equivalence: deg f ≤ D ⟺ f ⊥ RM(m-D-1,m) ⟺ XOR over every
(D+1)-flat vanishes (minimum-weight codewords of RM(m-D-1,m) are exactly the
(D+1)-flat indicators and generate the code). (1) is Rédei. (4) is computation.

## Verification

04-computation/rm_digit_ladder_macmini_0611s1.py Part B (exhaustive over all
2^m tilings at n = 4..7; ANF via Möbius transform); output
05-knowledge/results/rm_digit_ladder_macmini_0611s1.out. Also: c_3 mod 2 has
degree exactly 2 at n = 4..7 (the triangle case of the pairing argument:
abc + (1+a)(1+b)(1+c) = 1 + a+b+c + ab+ac+bc).

## Notes

- The degree-2 triangle identity is the seed: a single cyclic-triple's parity is
  QUADRATIC in the arcs — this is the same quadratic that makes the oriented
  two-graph g(x,y,z) = f(x,y)f(y,z)f(z,x) a 2-cocycle (Babai–Cameron), and the
  c_odd-sum keeps the quadratic ceiling only through the band limit D.
- Reading kind-pasteur's THM-466 through this lens: the 2-adic digit tower of H
  consists of odd-cycle COLLECTION counts; this theorem says the depth-1 layer is
  a low-degree (RM(D,m)) function while depth ≥ 3 layers are full-degree — the
  tower's complexity is graded by Reed–Muller level.
- Open problem leverage (logged HYP-2408): Ax 1964 (precise form: # ZEROS of a
  degree-d polynomial on F_2^m is divisible by 2^(ceil(m/d)-1); the ones-count
  version follows by complement) constrains level-set sizes of digit patterns; the
  forbidden values H ∈ {7, 21} are joint digit-level-set emptiness statements —
  a possible second proof route / strengthening of the forbidden-H phenomena.
- Literature check (2026-06-11): the ANF degree of 2-adic digits of combinatorial
  counting functions (permanents, Hamiltonian counts) as functions of input bits
  appears UNPUBLISHED. The model case is classical: digit i of the Hamming weight
  is e_{2^i} mod 2, degree exactly 2^i (Canteaut-Videau 2005). Note the contrast:
  H's ladder is 0, D, ~m — flat, jump, saturate — NOT the 2^i geometric climb of
  the weight function. Computational adjacents: permanent mod 2^k in P (Valiant;
  Bjorklund-Husfeldt-Lyckberg IPL 2017), parity of Hamiltonian cycles in
  O(1.619^n) (Bjorklund-Husfeldt FOCS 2013).
