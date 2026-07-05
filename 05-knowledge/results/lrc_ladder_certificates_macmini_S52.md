# The 14r-ladder witness certificates (kps-S2 rational_point_margin atom rows)

mac-mini-2026-07-05-S52, HYP-4103.  Machine-checked exact (fractions); consumers:
kps LRCHarmonicGate.rational_point_margin (each row is one atom instance), klein
S134/S135 assembly (the single-lift stratum's floor witnesses), opus LRCLiftRigidityRows.

## The law

For r = 7..12, the sieve-surviving single lift at residue r is k = r, i.e. the family
W_r = ({1..12} \ {r}) u {14r}.  Its exact covering-min and witness:

    M(W_r) = 14/(13(r+1)),   t*_r = a_r/(13(r+1)),   (13-r)*a_r == -14  (mod 13(r+1))

(killer 14r == (r+1) - 14 mod 13(r+1), so the dilation puts the killer at residue
EXACTLY +14 while the base binder 13-r sits at -14: a TWO-POINT EQUIOSCILLATION
(THM-618 species), binding pair (13-r, 14r) on opposite sides of 0.  The law is a
CONGRUENCE, not an inverse: g = gcd(13-r, 13(r+1)) in {1,2} and solvability is
g | 14 -- exactly the condition for the killer to land at +14.  g solutions, base
clearance selects; mirrors a <-> q-a give the same margin.)

## The rows (all verified exact; margin numerator uniformly 14)

| r  | W_r                              | q=13(r+1) | a_r (mirror)  | margin mu/q      | binders (at -+14)  |
|----|----------------------------------|-----------|---------------|------------------|--------------------|
| 7  | {1..6, 8..12, 98}                | 104       | 15 (89)       | 14/104 = 7/52    | 6 = 13-r, 98 = 14r |
| 8  | {1..7, 9..12, 112}               | 117       | 44 (73)       | 14/117           | 5, 112             |
| 9  | {1..8, 10..12, 126}              | 130       | 29 (101)      | 14/130 = 7/65    | 4, 126             |
| 10 | {1..9, 11, 12, 140}              | 143       | 43 (100)      | 14/143           | 3, 140             |
| 11 | {1..10, 12, 154}                 | 156       | 71 (85)       | 14/156 = 7/78    | 2, 154             |
| 12 | {1..11, 168}                     | 169       | 14 (155)      | 14/169           | 1, 168 (deep well) |

Atom shape per row: for every v in W_r: 14 <= (v * a_r) % q <= q - 14   (pure integer).

## Floor statements these feed

- SINGLE-LIFT FLOOR (closed this session at floor level): every single lift
  r -> r+13k (k >= 1) has margin >= 14/169; the structural cutoff for the floor
  question is killer <= 2016 = 14*144 (swept k <= 155; beyond, the beta*-level
  window closes: bad arc 2*beta*/w shorter than the Lipschitz window 2*(1/12-beta*)/12).
  Floor attained uniquely at W_12 = {1..11,168}.
- MULTI-LIFT FLOOR = 2/25, ATTAINED: M({1..12}\{4,6} u {17,19}) = 2/25 EXACTLY
  (witness t = 6/25, binders runner 8 at -2, runner 17 at +2 mod 25); the unique
  structured slice below 14/169.  Height-1 hypercube (4094 primitive C): zero
  scan failures at 2/25.  Full l=2 structural domain (600,756 sets, w_b <= 258,
  w_a <= 24 w_b, killers to ~6200): ZERO sets below 2/25, zero exact escalations
  -- with the fee (both >= 259) and window (w_a > 24 w_b) closures, EVERY double
  lift has margin >= 2/25.  The spectral gap (1/13, 2/25) is EMPTY on all swept
  lift ground.  See lrc_multilift_floor_macmini_S52.out.
- Modulus-bound data for kps-S3 (certificate completeness): the ladder's needed
  modulus is 13(r+1) = killer + 1 + (13 - r) <= killer + 7 for r >= 7; block-lift
  certs use q in {25, 50, 13s}.  Supports a linear-in-w_max modulus bound.
