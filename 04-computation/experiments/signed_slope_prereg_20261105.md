# Pre-registered predictions: signed one-apex responses a n - b g(n) by parity and slope (opus-2026-11-05-S80)

Written before any computation. Setting: g(1) = 1, g(2n) = g(n), g(2n-1) = g(n) + n - (n mod 2) (Session198's
arithmetic image of the one-apex O/Z orbit: g(n) is the Hamiltonian-path count of the one-apex tournament whose
ordered two-cover count at the apex is n, i.e. whose chain word has binary value n - 1). S(a, b) = { a n - b g(n) :
n >= 3 odd }. Session198 PROVED: for coprime a even, b odd, S(a, b) = all odd integers iff b/2 < a < b, and gave the
3/4 control (positive values of 3n - 4g(n) are 1 mod 4). The probe computes, for n < 2^24, which targets in [-999, 999]
of the right parity are attained ("small-witness coverage": a theorem's witnesses may need larger n, so a miss here
refutes nothing proved; a hit is exact).

- P1 (Session198's class, a even, b odd, 1/2 < a/b < 1: (2,3), (4,5), (6,7), (8,9), (4,7), (6,11), (8,11), (8,13),
  (10,11), (10,13)): every odd target in [-999, 999] is attained with n < 2^24 except possibly a few of the largest.
- P2 (both odd, 1/2 < a/b < 1: (3,5), (5,7), (5,9), (7,9), (7,11), (9,11), (7,13), (9,13)): the image consists of even
  integers, and every even target in [-998, 998] is attained with n < 2^24 except possibly a few of the largest
  (prediction: Session198's clock extends with a boundary term).
- P3 (a odd, b even, 1/2 < a/b <= 3/4: (3,4), (5,8), (7,10), (7,12), (9,14), (11,16), (9,16)): every positive value
  is congruent to -(a + b) modulo 4; hence at least half of the positive odd targets are missed (this is the claim to
  be PROVED: the m-odd branch is negative for a/b <= 3/4 by the strong cone, and the m-even branch is a fixed class
  modulo 4).
- P4 (a odd, b even, 3/4 < a/b < 1: (7,8), (5,6), (9,10), (11,12), (13,16), (11,14), (13,14), (15,16)): both odd
  classes modulo 4 occur among positive values; prediction (low confidence): every odd target in [-999, 999] is
  attained with n < 2^24 except possibly a few of the largest, i.e. the parity hypothesis of Session198 is not
  necessary above 3/4.
- P5 (outside the interval: (1,3), (2,7), (4,3), (5,4)): the image misses one sign entirely, as Session198's sharpness
  argument says.
