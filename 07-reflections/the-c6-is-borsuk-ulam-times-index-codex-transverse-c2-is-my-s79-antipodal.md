# The C₆ is Borsuk-Ulam × index: codex's transverse C₂ = Q(√−7) is exactly my S79 antipodal, and the magnitude layer is the doubling-shadow 2U = 2·U

*mac-mini-2026-06-28-S85. The owner: continue the S84 CRT factorization, lose no work, extend it, connect with
the incoming concurrent ideas. codex extended my S84 reflection in real time — it named the field `ℚ(ζ₇)` with
`Gal = C₆ = C₂ × C₃` and flagged a coordinate I had missed: the **C₂ = ℚ(√−7) quadratic character is transverse**
to my C₃ binding-pair skeleton. This note NAMES that C₂ — it is exactly the **Borsuk-Ulam antipodal** of my S79
index-theorem — and links the two layers by the doubling map. Builds on
[[the-census-factors-via-crt-7-adic-residue-c3-skeleton-times-2-adic-magnitude-doubling-hinge]] (S84 + codex C₆ guardrail, HYP-3310/3311),
[[the-hidden-c3-the-witness-space-is-one-galois-orbit-leverage-the-proved-binding-pair-plus-c3]] (S83 C₃, HYP-3257),
[[the-index-theorem-frame-lrc-is-the-nonvanishing-of-an-index-analytic-equals-topological]] (S79), kps census split (HYP-3258).*

## What codex added, and the coordinate it left unnamed
codex's C₆ guardrail (its extension of my S84 note) is the correction the slogan needed: the residue field is
`ℚ(ζ₇)`, with `Gal(ℚ(ζ₇)/ℚ) = C₆ = C₂ × C₃`, and the proof packet has **three** legal coordinates, not one:
```
  C₃  = the binding-pair orbit  {1,13}→{3,11}→{5,9}     (my S83; the real cubic ℚ(cos 2π/7))
  C₂  = the quadratic character χ₇, TRANSVERSE to C₃     (codex; the imaginary quadratic ℚ(√−7))
  2-adic = the height/flex ledger on 2U ∪ {7}            (kps's hinge; the blind magnitude data)
```
codex proved C₂ is transverse: **every binding pair `{a,−a}` carries one QR and one NQR mod 7** (verified —
`{1,6}: (+,−)`, `{3,4}: (−,+)`, `{5,2}: (−,+)`). But it left C₂ unnamed. Here is its name.

## C₂ is the Borsuk-Ulam antipodal (the S79 structure), because −1 is a non-residue mod 7
The transversality has one cause: `χ₇(−1) = −1` (since `7 ≡ 3 mod 4`), so `a` and `−a` always have **opposite**
quadratic character. That is *precisely* the fact my S79 index-theorem runs on:
> `7 ≡ 3 (mod 4)` ⟹ `−1` is a non-residue ⟹ complement `×(−1)` is an **anti-automorphism** ⟹ a free `ℤ₂`
> action ⟹ **Borsuk-Ulam** (not Brouwer), witness = an antipodal pair, certified by an **odd degree**.
So **`C₂ = {±1} = χ₇ = the antipodal map = the Borsuk-Ulam `ℤ₂`.`** codex's transverse quadratic coordinate and
my S79 antipodal structure are the same object, seen from the number field vs. the topology. Therefore:
```
  C₆  =  C₂           ×        C₃
      =  (BU antipodal, ℚ(√−7))  ×  (index degree (p−1)/2 = 3, ℚ(cos 2π/7))
```
- **C₃** says *which* binding pair (the index/degree, S83); fixed field `ℚ(cos 2π/7)` (the real cubic).
- **C₂** says *which element of the pair* (the QR/NQR sign = the antipodal orientation, S79); fixed field
  `ℚ(√−7)` (the imaginary quadratic — the Heegner field, where the S77 norm-gem lived).
The index-theorem (S79) is exactly the statement that the **residue** layer is `C₆`: it uses the C₂ antipodal to
force an odd degree equal to the C₃ index. **What it cannot see is the 2-adic height** — codex's "blind data," the
hard residual (kps S256: "the frame describes; the floor proves").

## The two layers are linked by ×2: the magnitude is the doubling-shadow of the skeleton
The residue and magnitude layers are not independent — the magnitude evens are the **doubling-image** of the unit
skeleton:
```
  2·U mod 14 = {2,4,6,8,10,12} = the even covering layer (the magnitude shadow);   nonzero(14) = U ⊔ 2U ⊔ {7}
```
So `nonzero(ℤ/14) = U` (skeleton) `⊔ 2U` (doubling-shadow = magnitude) `⊔ {7}` (ramified leftover). The mod-14
**residue** of every magnitude runner is forced (it is `2·`a unit); only its **integer height** (the lift `12` vs
`24` vs `26`) is free. That is the exact split: *residue determined by ×2 from the skeleton; height = the blind,
free coordinate where the one Jacobsthal hinge `12→24` lives.* The doubling map is the bridge between the C₆
residue structure and the 2-adic height ledger.

## Connecting to kps S257: why the topological degree closes only the residue half
kps S257 (same day) independently fused my C₃ into the full Galois picture and added the decisive fact:
- **The 13 runners = two regular C₆-orbits + the fixed apex:** units `{1,3,5,9,11,13}` (binding skeleton) ⊕ evens
  `{2,4,6,8,10,12}` (covering) ⊕ `{7}` (fixed, `7u≡7`); `6+6+1=13`. Two regular reps + the trivial.
- **The two subfields carry the two proof halves:** `ℚ(cos2π/7)` (C₂-fixed, real cubic) carries the cap = C₃-trace
  + equioscillation; `ℚ(√−7)` (C₃-fixed, imaginary quadratic) carries the floor / Gauss sum `i√7` / the ±
  complement. kps's C₂ = "complex conjugation = complement `a↦−a`" is **exactly my Borsuk-Ulam antipodal** — the
  field-theoretic and the topological name of the same involution (THM-568 = kps's name for the proved binding
  pair HYP-2909).
- **The dichotomy (verified by kps):** the residue half is C₆-symmetric (`×u` permutes the AP, `M=1/14` preserved)
  — CLOSED. The magnitude half is C₆-BROKEN — only `12` of the evens-orbit admits a tight doubling (`12→24`=GW);
  `2,4,6,8,10` admit none.

**The synthesis my C₂ = BU naming forces:** the Borsuk-Ulam / topological-degree argument is **q-uniform** — it
runs identically for every LRC(2p) (odd degree `(p−1)/2`, no dependence on the specific `p` beyond parity). But
kps S257 shows the magnitude break is **q-specific**: a tight doubling exists for `q=4,7` yet NOT for `q=5,6`
(`n=10,12` have census `{AP}` alone). **A q-uniform invariant cannot detect a q-specific phenomenon.** So the whole
topological/Galois/charge apparatus (S79 BU degree, the C₆-symmetry, the conserved charges) can only ever prove the
q-UNIFORM part — the residue half — and is *structurally blind* to the magnitude break, which is gated by the apex
prime's arithmetic (the gcd-window / Jacobsthal switch, kps). This is the precise mechanism behind "the frame
describes, the floor proves" (kps S256): the conserved (topological) charges are q-uniform = the residue seed; the
un-conserved magnitude coordinate is q-specific = the hard core. The BU degree and the q-switch sit on opposite
sides of the conservation line.

## Connecting to the conservation language (codex HYP-3400 / HYP-3301)
codex frames the three coordinates as **shadow-charge reservoirs** that a scalar proof-shadow must preserve,
transfer, or pay off as named debt. Naming C₂ makes one charge concrete: **the C₂ charge is the Borsuk-Ulam
degree — odd, a topological invariant, hence automatically conserved.** A proof may *never* silently forget C₂:
dropping it would change a homotopy degree from odd to even, which no continuous deformation can do. So of the
three reservoirs, two are "rigid" conserved charges (C₃ index, C₂ BU-degree — the residue layer, the seed) and the
third (2-adic height) is the *un-conserved* coordinate where the proof must spend real work (the deformation-flex
rigidity, kps's pointer). codex's HYP-3301 first-obstruction cocycle, in this packet, is therefore the obstruction
living in the **2-adic height** factor — the only coordinate not protected by a Galois/degree invariant.

## Honest status
- **EXTENSION (structural, verified):** codex's transverse C₂ = ℚ(√−7) character = my S79 Borsuk-Ulam antipodal
  (both are `χ₇(−1)=−1`, `7≡3 mod 4`). So `C₆ = C₂(BU) × C₃(index)`; the index-theorem (S79) IS the C₆ residue
  structure. The magnitude layer is the doubling-shadow `2U = 2·U`; residues forced, heights free.
- **CONNECTION:** unifies S79 (BU), S83 (C₃), codex's C₆ guardrail (C₂-transverse), kps's hinge (2-adic), and
  codex's conservation language (C₂ charge = BU degree, conserved). One field `ℚ(ζ₇)`, two subfields
  `ℚ(√−7)`/`ℚ(cos2π/7)`, one doubling map, one blind height ledger.
- **NOT a proof.** This is the structural accounting: the residue layer (`C₆` = two conserved charges) is the seed;
  the 2-adic height (one un-conserved charge, the hinge) is the hard open residual — exactly where codex's blind
  data, kps's deformation flex, and HYP-3301's first obstruction all live. LRC(14) open.

Related: HYP-3312 (this), HYP-3311/HYP-3310 (codex C₆ synthesis + sidecar), HYP-3257 (my S83 C₃), HYP-3258 (kps
census split), HYP-3400 (codex shadow-charge), HYP-3301 (codex first-obstruction), HYP-2909 (proved binding pair),
S79 index-theorem, S77 (ℚ(√−7) norm-gem), OPEN-Q-108.
