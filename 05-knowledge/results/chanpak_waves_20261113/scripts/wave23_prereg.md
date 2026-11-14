# Wave 23 pre-registrations (written before any wave-23 computation was run)
timestamp: 2026-11-14T21:15:54Z

P46 (rank ladder, exact algebra). Transfer algebras (associative algebra spanned by the letter-matrix words):
  (a) one-apex O/Z (4-dim state) and one-hub words: triangularizable (semisimple quotient commutative);
  (b) O/Z/J root-change words and two-hub words: NOT triangularizable (semisimple quotient contains a
      matrix block of size >= 2).
  Prediction: (a) and (b) both hold.

P47 (no Jacobi-symbol invariant in the GL3-dense orbit). On the orbit of the root-change block semigroup
  (states (H,E,lam) reached by block words), for every pair of primitive integer linear forms (a.v, b.v) with
  coefficients in [-2,2] and at least 2000 usable samples (b.v odd positive, gcd 1), the Jacobi symbol
  (a.v / b.v) takes each value +1 and -1 with frequency in [0.40, 0.60], except pairs whose symbol is
  forced by a congruence visible mod 8 (to be listed).

P48 (Apollonian test for F1). F1 = g(h1)g(h2) + 2^s(2h1+1)(2h2+1) - [s=0] on odd N in [2^22, 2^27]:
  the miss rate of F1 on odd squares m^2 divided by the miss rate on locally matched controls
  (same residue modulo 8*3*5*7*11*13) lies in [0.5, 2.0], and likewise for t*m^2, t in {3,5,7,11,15,21}.
  (No reciprocity obstruction of square-class type in the rank-one family F1.)

P49 (registered 2026-11-14T21:26:54Z, after the in-sample F1 square-class run on [2^22,2^27), before any
  run on [2^27,2^29)).  In-sample, with strata (w0,w1 over primes 11..97, N mod 16, N mod 105), the classes 7m^2
  (7 not dividing m) and 3m^2 had observed/expected F1 misses 17/5.6 (z = +4.8) and 10/3.8 (z = +3.1); squares 22/22.5.
  Prediction (the effect is real): on the disjoint window [2^27,2^29) with the same strata model, the 7m^2 class has
  observed/expected >= 1.8 and the 3m^2 class >= 1.5, while squares stay in [0.6, 1.5].
