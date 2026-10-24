# Dyadic support of the weight-one root modulo four

Status: PROVED CONSEQUENCE of inherited TheoremB at2 and its separately audited exceptional Frobenius lift. This is an elementary extraction from those inputs, not an independent reproof of their modular hypotheses or a novelty claim. Earlier frozen packets are unchanged.

Let G(t)=sum_{n>=0}c(n)t^n be the integral normalized formal square root of the s7 generating series, c(0)=1. Then in (Z/4Z)[[t]],

G(t)=1+2 sum_{j>=0}t^(2^j).

Equivalently, for every positive integer n, c(n)=2 modulo4 if n is a power of2, and c(n)=0 modulo4 otherwise. At n=0 the coefficient is1, not0 or2. In particular 4 divides c(n) on a set of natural density1.

## Exact input pins

The included inherited-proof.md is the repository source cooper_s7_weight_one_root_frobenius_20261018.md from commit7818149bf0a1eca778549fdfd1dde34f479421fd (source SHA256 fd7935e86a41854f67817d4f1b3b16d4fcd52134a72757debd2a9d1489ea66b4).

- Section3, TheoremB at split prime2: G=(1+2t)G(t^sigma) modulo4.
- Section5.6: t^sigma=t²(1+2t) modulo4. Its original proof clears denominators into a weight12 level14 modular form and checks past Sturm bound24. This supplement inherits that identity and does not claim a new modular certificate.
- Section1: G is the integral normalized series under the stated modular parametrization. Normalization and integrality are not inferred from an arbitrary square root of an integer series.

## Proof with precision audit

Reducing the two inputs modulo2 gives G(t)=G(t²). If there were a nonzero positive coefficient modulo2, choose its least index m. If m is odd, the coefficient at t^m on the right is zero. If m is even, that coefficient equals c(m/2), zero by minimality. Both are contradictions. Thus G=1 modulo2, meaning every positive c(n) is even; the constant is1.

For the formal Euler derivative D=t d/dt, the binomial identity gives

G(t²(1+2t))=G(t²)+2t(DG)(t²) modulo4.

This formula uses integer coefficientwise expansion and no division by an index or factorial. Because every positive coefficient of G is even, every coefficient of DG is even. Its constant coefficient is zero. Therefore the correction is zero modulo4, and the excellent lift may be replaced here by the naive one:

G(t)=(1+2t)G(t²) modulo4.

Iterating r times yields

G(t)=product_{j=0}^{r-1}(1+2t^(2^j)) G(t^(2^r)) modulo4.

For any fixed coefficient index N, choose 2^r>N. The tail factor then equals1 through degree N, and every product involving two nonconstant factors is divisible by4. Thus the product is exactly 1+2 sum_{j>=0}t^(2^j) in the coefficientwise limit. This proves all indices without analytic convergence or finite extrapolation.

The positive exceptional indices below X are exactly the powers of2, at most 1+log_2 X. Their proportion tends to zero. The n=0 exception does not affect natural density. This is a coefficient-divisibility statement, not an orbit-coverage statement.

## Sharp scope and countermodels

- The parity argument is needed before discarding the lift correction. For the integral series A=1+t, A(t²(1+2t))-A(t²)=2t³ modulo4, which is nonzero. Hence a naive lift replacement is not a formal identity for arbitrary normalized integral series.
- The constant is essential: A=0 also solves A(t)=A(t²) modulo2, but is not normalized. The actual result uses G(0)=1 explicitly.
- The lacunary expression does not lift unchanged modulo8: c(2)=22 is6 modulo8, not2. Also c(4)-c(1)=6006-2=6004 is4 modulo8. Thus even the two-step fixed-point congruence C4G=G cannot be asserted modulo8.
- No claim is made here about the separate observed identity C4G-G=4t modulo8. It is outside this proved supplement, not an inferred theorem.
- All positive coefficients are even, but not all are divisible by4: c(1)=2 is the smallest hostile.

## Targeted prior extraction search

A bounded case-insensitive search for powers of two, powers of2, modulo4/mod4 and 2^j/2^k in the inherited root note and the related materialized notes beukers_remark114_power_structure_and_inert_cm_20261018.md and pr93_modular_bridge_20261018a.md found the inherited section5.6 one-step congruence, but no explicit power-of-two support extraction. The three candidate notes were selected by a bounded A279619/weight-one-root/Cooper-s7 text search. This is not an exhaustive repository or literature search and establishes no originality.

## Reproduction

check.py computes c through10000 by exact recurrence with checked integral division, independently verifies the defining binomial-sum square through240, tests the full mod4 support formula and finite infinite-product truncation, and checks all countermodels above. Coefficients are never computed using the conjectured residue law. Normal/-O outputs match. These finite checks support the separately proved all-index consequence; they are not its proof.

Run python3 replay.py from any working directory. It checks recursive exact file membership, all hashes, and normal/-O fixture equality. Run hostile.py normally and with -O to verify rejection of corruption, missing members, and unexpected nested files. No external package, network, or repository access is needed for replay.
