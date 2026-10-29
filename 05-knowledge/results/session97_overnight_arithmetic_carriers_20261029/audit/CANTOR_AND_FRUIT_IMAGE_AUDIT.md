# Independent audit: exceptional-address recurrence and fruit image density

## Cheap-controller Cantor set

**PASS** `controllers/EXCEPTIONAL_CANTOR_RECURRENCE.md`, including its
self-contained Hausdorff-dimension proof. Inverse similarities have ratios
2^-1,2^-6,2^-11,2^-11,2^-11 and disjoint full source cylinders. The mass
bound uses the preceding cylinder before a precision crossing; its depth is
at least L-11, so the stated Frostman constant is valid. The lower occupied-
cell bound must not be read as a positive-integer existence result.

`check_cantor_independent.py` takes a different computational route: start
with *all* odd residues at each finite precision, repeatedly apply the inverse
system, and stop at the greatest fixed set. Its occupied counts agree through
precision17 with the initial data and all applicable recurrences. Exact
Fraction checks certify the lambda bracket, exponent interval(0.441,0.442)
and upper constant17/8. Normal/-O outputs agree.

For a real bound X>=1, L=floor(log2 X) gives modulus2^(L+1)>X. At most one
positive integer up to X lies in each occupied cell, giving the claimed
(17/8)X^s bound. This counts only integers whose cheap first returns all
expand forever. It is not a bound for all divergent Collatz integers and
neither establishes existence nor emptiness of the exceptional set.

## Full fruit isogeny image

**PASS conditional on the explicitly CITED saturated source group**
E(Q)=ZG+C6. The target torsion argument uses distinct good primes11 and17,
with counts12 and18, and a complete exclusion of3-torsion via the division
polynomial. The dual-isogeny relation forces H primitive, rather than
assuming target saturation from a bounded search. The target quotient is C3.
A mod2 obstruction for H then injects this global quotient into the local
quotient, proving rational-image membership iff Q2-image membership.

`check_fruit_full_image_independent.py` independently replays the torsion
counts, division polynomial and root obstruction. At the good prime3 it
finds six target points, with H reducing to a generator and3H to U. It
verifies162 signed-index/torsion targets and the equivalent cheap test for a
rational cover coordinate q=a/d: 3 divides d or a+d. A rational W satisfying
the cover equation, actual point-at-infinity handling, and the torsion
fruit-denominator boundary remain mandatory. The test is not for arbitrary
rational q without its cover point.

## One-third positive-window law

**PASS** `fractions/DENSITY.md`. H lies in the bounded real component, U in
the identity component, and positive real roots occur only on the bounded
component. For fixed e, the three odd index classes modulo6 are each an
irrational rotation by6H on that same circle component. They hit the open
positive window with the same Haar frequency mu. Only index3 modulo6 is
arithmetically liftable, giving the exact limiting conditional share1/3.
Canonical-height ordering is equivalent to index ordering via n²*hhat(H);
no claim under naive numerator/denominator ordering follows.

The stated integral for mu correctly includes both W branches, whose common
factor cancels. Substitution q=q_minus+Delta*sin²(theta) introduces the
common factor1/7, also cancelling. Endpoint singularities are integrable.
The general clock rule gcd(c,d)/d follows from compatible residue classes
and irrational rotation on each progression; incompatible classes give0.
There is no finite-cutoff independence or effective convergence-rate claim.
