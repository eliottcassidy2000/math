# Pre-registration: finite O/Z tail banks - Poisson law with local factors, root-change trichotomy

kind-pasteur-2026-11-06-S13. Committed and pushed **before** the confirmatory runs. Reservations: THM-4887, HYP-9480.

## Objects (inherited)

Tail formula (Session 205, audited): a physical unit-start seed s = (h, R, e, d) = (D, U, E, d) followed by an O/Z
word with parameter m has H_s(m) = h + R m - d Q(m) - e (m mod 2), Q(0) = 0, Q(2m) = Q(m),
Q(2m+1) = Q(m) + m + (m mod 2). Capacity (Session 206, audited): #{m : H_s(m) <= X} = c_s X + O(sqrt X),
c_s = ln(R/(R - d))/d. Root (rooted edge t = (1, 2, 0, 1)): converse involution T with H(T(m)) = H(m) (Session 206).
Root changes: J(D, E, U, d) = (D - E + U, D - E, D + U, D). F_k = H-values of t . w_1 (J w_2) ... (J w_(k+1)).

## Claims to be proved (THM-4887; proofs in PROOF.md, not decided by computation)

- A (mod-3 law). Q(m) = m^2 - r(m) (mod 3), r(m) the number of maximal runs of 1-bits of m. Hence
  H_s(m) = h + R m - d m^2 + d r(m) - e (m mod 2) (mod 3) for every seed.
- B (involution rigidity). For every seed, H_s(T(m)) - H_s(m) = (R - 2d)(T(m) - m) - e((T(m) mod 2) - (m mod 2)).
  Among physical seeds, H_s(T(m)) = H_s(m) for all m holds iff R = 2d and e = 0 iff s = t.

## Local-factor Poisson prediction (the model under test)

Split the parameters of s into the classes of m mod 2. If an odd prime p <= 13 divides gcd(d, R), every value with
m of given parity lies in one residue class mod p (h, resp. h - e). Let P be the product of such primes. A piece j
with intensity lambda_j (parameters per odd integer) landing in a set of odd residues mod 2P of relative size f_j
is predicted to cover the fraction f_j (1 - exp(-lambda_j / f_j)) of the odd integers; pieces landing in the same
classes are merged first. With no such prime, the prediction is 1 - exp(-2 c_s). For the root it is
1 - exp(-ln 2) = 1/2 (the involution halves the intensity).

## Confirmatory predictions (new code; normal and -O receipts must agree)

- **P1 (single seeds).** For every physical J-born seed (last letter J) with 64 < U <= 160, d > 1, R != 2d and
  R != h, the hit density among odd N in [2^24, 2^25) differs from the local-factor prediction by at most 0.004 for
  at least 95% of the seeds and by at most 0.010 for all of them.
- **P2 (root).** The root image's hit density among odd N in [2^27, 2^28) lies in [0.4978, 0.5005].
- **P3 (fixed bank).** Bank A (the root and the 84 physical J-born seeds with d > 1 and U <= 64) misses a fraction
  in [0.0205, 0.0235] of the odd N in each of [2^26, 2^27) and [2^27, 2^28).
- **P4 (two root changes).** For F_(<=2), computed exactly to 10^8: (a) at most 6 odd N in (4 x 10^6, 10^8] are
  missed; (b) of all misses above 10^5, at most two are not 1 (mod 3). For F_(<=1): (c) the miss counts per dyadic
  window keep growing, with every ratio between consecutive windows from [2^16, 2^17) on in [1.45, 1.75].
- **P5 (THM checks, not evidence for the proofs).** A on all m < 2^26; B's identity on every physical seed with
  U <= 64 and m < 2^16; among all physical states with U <= 4096, only t has R = 2d and e = 0.

A failed prediction is recorded as REFUTED in HYP-9480, not repaired after the fact.
