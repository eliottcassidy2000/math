# Chan–Pak: reverse-complement collisions, an independent thirteen-unit certificate, and certified hosts for prescribed 16-digit primes

Session klein-2026-11-15-rc. **Chan–Pak remains OPEN.**

Owner prompt (2026-11-15):
- choose a host that hits the prescribed absolute count;
- build a cone certificate for the 13 units;
- attack the twisted-collision lemma, aiming at new methods rather than a longer verified range.

A concurrent packet, [chanpak_cone_collision_placement_20261115](../chanpak_cone_collision_placement_20261115/README.md),
finished the cone certificate first and gave the two-ended lattice-line equation. This packet credits both and
adds the following.

1. **[Reverse-complement collisions](RC_COLLISIONS.md) — NEW, PROVED ([THM-5232](../../../01-canon/theorems/THM-5232-reverse-complement-collisions-every-odd-core-with-its-reverse-complement-collides-on-an-anti-diagonal-so-the-row-collision-bound-is-sharp-and-one-explicit-family-carries-linear-energy.md)).**
   - *The lemma.* The checksum p is invariant under reverse-complement rc_a(x) = 2^a − 1 − rev_a(x), and
     x − rc x = rev x − rev rc x.
   - *The family.* Hence every odd core with top bit 1 collides with its reverse-complement along a whole
     anti-diagonal of (prefix, reversed suffix) pairs, whenever k | c 2^(j−1).
   - *Sharpness.* The family gives N(l, w) ≥ (c/2k) 2^l/w on a positive-log-proportion set of prefixes in every
     dyadic block. So the conjectured bound of [HYP-9752](../../hypotheses/HYP-9752-row-collision-lemma-a-carrier-row-cu-plus-kp-u-has-linear-energy-which-would-give-the-chan-pak-spectrum-positive-lower-density.md)
     is sharp in shape, with C ≥ 3 for the 3-cycle row.
   - *Energy.* This single explicit family has linear size. It carries 44% of the 3-cycle row's collision energy
     at L = 18, and explains the linear growth of the per-level totals.
   - *Dichotomy.* Rows with odd(k) ∤ c lose the family and have 50–250 times less energy.
   - *What a proof of the upper bound needs.* Subtract the explicit lattice-ray families, then prove an
     equidistribution of (p(x) − p(y))·d⁻¹ mod δ for the generic single-point lines.
2. **[Thirteen-unit cone, independent replication](CONE13_REPLICATION.md) — PROVED (computer-assisted).**
   - Every bottleneck-free word over all 13 span-3 floor units has H ≥ 0.013402·1.7699^|W|.
   - Cones are images of the nonnegative orthant of the ribbon profile transfer. Uniform depth 4 replaces the
     adaptive schedule, and the constant improves to 0.0134 from about 1.5·10⁻⁷.
   - The verifier is exact and independent; uniform depth 3 stops near 1.7675.
3. **[Certified hosts for prescribed primes](HOSTS.md) — DEMONSTRATION.**
   - The two-tape carrier closed form turns a prescribed N into an arithmetic search.
   - Explicit tournaments with exactly N Hamiltonian paths are certified for the largest primes below 10^12,
     10^13, 10^14, 10^15 and 10^16 (42–54 vertices). Each count is checked twice: by the closed form and by an
     independent cut transfer.
   - Small instances are brute-forced end to end; for example 65537 is H of an 18-vertex tournament.
   - The cost is about √N per hit, matching the concurrent inverse.

## Concept board after this session

- *Collisions.* The collision energy splits into explicit lattice-ray families plus a generic remainder. The
  families are reverse-complement cores and the fixed-core rays. The remainder is visibly small when the family
  is absent (189U + 11p). The open estimate is now a first statement about generic lines.
- *Converse symmetry.* Converse symmetry (H(T) = H(T^op)) is the source of the largest collision family. A
  "partial converse" reading of other near-collision families is a lead, not a theorem.
- *Hosts.* The pointwise first moment (every N is hit), not the second moment, is what Chan–Pak needs. Second-moment
  progress (collision energy) gives density at best.

## Replay

| command | time | checks |
|---|---:|---|
| `python3 verify_rc_collisions.py` | ~10 s | Lemma 1, family members, closed count, sharpness |
| `python3 verify_cone13.py ../../../04-computation/experiments/ribbon_spectrum_20261104 ../chanpak_waves_20261113/scripts/w30_hankel_13units.pkl cert13_depth4.pkl` | ~50 s | the 13-unit theorem |
| `cc -O2 -o hostsolve hostsolve.c -lm && echo 999999999989 \| ./hostsolve cores35.txt 1 64` | < 1 s | one 12-digit host |
| `cc -O2 -o hp2 hp2.c && python3 bf_check.py small_hosts.txt ./hp2` | ~1 s | brute-force end-to-end controls |

`cores35.txt` comes from `coredata.c`, a main() appended to the shared functions of
`../chanpak_waves_20261113/scripts/unioncar.c`; it needs `gentourng`.
