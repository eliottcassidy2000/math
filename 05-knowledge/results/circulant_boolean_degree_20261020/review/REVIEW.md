# Independent adversarial review: Boolean-degree boundary

2026-10-20. **ACCEPTED: PROVED in the stated scope, with reproduced
FINITE-EXACT controls.** No mathematical correction is required.

This review covers `session13-nonlinear-degree/BOOLEAN_DEGREE_BOUNDARY.md`,
especially §§3–6, together with the explicit witnesses, verifier, quadratic
census, and inherited prime-circulant dictionary. The reviewed main-note SHA-256
is:

    cbdf5e4566ddec0be6e337302df683d18bf7524db7eff8d5f518a8866b300b55

All source hashes are recorded in `SOURCE_SHA256SUMS`. The five sources frozen
at review startup remained byte-identical at close. No producer file was
changed, and nothing was committed or published by this reviewer.

## 1. Verdict and exact scope

Let \(m=e q\), \(q\) odd, \(e=2^a\ge4\). Every map from the antiperiodic
selector torsor \(X_m\) to the primitive-even binary-word set \(B_m\), equivariant
for twisted source shift and ordinary target shift and bijective on cyclic
orbit sets, has Boolean algebraic degree at least

\[
 \left\lceil\frac{e-1}{a}\right\rceil .
\]

The proof is an all-order argument, not an extrapolation from quadratic
enumeration. The necessary degree is unbounded along actual primes
\(p=2m+1\). The minima at \(m=4,6,8\) are respectively \(2,2,3\).

The prime statement concerns labelled multiplier-equivariant maps inducing
bijections on isomorphism types. It is not a labelled bijection, an obstruction
to unrestricted nonlinear maps, a composite/prime-power graph-isomorphism
extension, or a theorem about first-order definability. The stated
multiplicative-depth consequence is correct only for constants, XOR, and binary
AND gates. The note keeps these distinctions explicit.

The mathematical dependency needed for the prime interpretation is the
prime-circulant dictionary in
`math/05-knowledge/results/circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md`,
§2. Its elementary Sylow/normalizer and sign-character reasoning was reread.
The general degree bound itself needs only the word objects and the dyadic
count proved within the reviewed note, not that predecessor's routing
construction or any finite computation.

## 2. Audit of the complete scalar parametrization

From \(F(Tx)=RF(x)\), the coordinate relation is
\(F_i(x)=f(T^{-i}x)\), where \(f=F_0\). Since \(T^m x=x+\mathbf1\) but
\(R^m=1\), \(f\) must be complement-invariant. Conversely, this invariance
supplies the wraparound coordinate relation and hence full equivariance.

The linear coordinate transformation

\[
 (x_0,\ldots,x_{m-1})
 \longmapsto
 (x_0+x_{m-1},\ldots,x_{m-2}+x_{m-1},x_{m-1})
\]

is invertible. Complement invariance means that the resulting Boolean function
is independent of its last coordinate. Thus its algebraic degree is exactly
its degree in the \(m-1\) quotient variables, and all other coordinates have the
same degree because \(T\) is invertible affine. Choosing another antipodal
transversal only permutes and complements input coordinates, which also
preserves degree.

Consequently the quadratic search really has dimension \(1+7+21=29\) at
\(m=8\), including the constant coefficient. It is not merely a local or
translation-invariant subfamily of equivariant rules.

## 3. Dyadic source orbits and target count

For dyadic \(m\), the order of \(T\) divides \(2m\), a power of two. Every
proper divisor of \(2m\) divides \(m\). Therefore a shortened orbit would be
fixed by \(T^m\), which is fixed-point-free complementation. Every orbit has
length exactly \(2m\), so there are \(2^m/(2m)\) source types.

For \(d=2^b>1\), every nonprimitive length-\(d\) word repeats a length-\(d/2\)
word and consequently has even total weight. The \(2^{d-1}\) even words
therefore contain exactly \(2^{d-1}-2^{d/2}\) primitive even words. Dividing
each primitive count by its actual orbit length \(d\), then adding the
period-one zero type, gives the displayed telescoping count
\(A(m)=2^{m-1}/m\).

The boundary cases are harmless: \(m=1,2\) each have one source and target
type, and \(m=2\) permits the degree-zero map. The theorem deliberately starts
at \(m\ge4\), where \(A(m)\ge2\).

## 4. Complement pairing and the exact Walsh normalization

For a nonzero good target type at dyadic length, its primitive period \(d\)
is even. Complementation preserves primitive period, and its primitive weight
changes from \(w\) to \(d-w\), preserving evenness. Thus nonzero good types
are closed under complementation. Zero is the one excluded-complement
exception, because the all-one word has primitive root \(1\).

Weight is invariant under rotation. Two distinct complement-paired types
have representative weights summing to \(m\). A type fixed by complementation
has weight \(m/2\), since its complement is a rotation of itself. Therefore

\[
 S:=\sum_{\text{good target types }[y]}\operatorname{wt}(y)
   =\frac m2(A(m)-1).
\]

A source orbit has \(2m\) inputs. Along that orbit, the zeroth output
coordinate samples each of the \(m\) output coordinates twice. This is true
even if the target's primitive period is a proper divisor of \(m\). Its
contribution to the full scalar truth-table weight is \(2\operatorname{wt}(y)\).
Type bijectivity then gives full scalar weight \(2S\).

Each quotient point has exactly two source representatives. Thus the
quotient scalar weight is \(S\), and the two distinct Walsh values are

\[
 W_{\mathrm{full}}=2^m-4S=2m,\qquad
 W_{\mathrm{quotient}}=2^{m-1}-2S=m.
\]

No orbit-size cancellation, self-complementary-type exception, or factor-of-two
error was found. This normalization is the decisive part of the lower bound.
As an additional hostile control, the valid \(m=6\) witness has Walsh value
zero in both spaces; one must not extend the dyadic forced-Walsh formula
directly to arbitrary \(m\).

## 5. Elementary Boolean divisibility and the degree bound

The integer expansion

\[
 (-1)^f=(-1)^c\prod_j(1-2M_j)
\]

is valid on Boolean inputs. A selected set of \(k\) monomials with support
union of size \(u\) contributes, after summing, a signed
\(2^{k+r-u}\). Since \(u\le\min(r,dk)\),

\[
 k+r-u\ge\max(k,r-(d-1)k)\ge\lceil r/d\rceil .
\]

The final inequality holds for every integer \(k\ge0\): if \(k\) is already
at least the ceiling it is immediate; otherwise \(k<r/d\), so the second
quantity is at least the required integer ceiling. The empty selected set
also satisfies the bound. Cancellation can only increase divisibility, not
invalidate it. Thus \(2^{\lceil r/d\rceil}\mid W(f)\).

The hypotheses \(1\le d\le r\) are met for a dyadic type bijection with
\(m\ge4\), because there are multiple types and hence the scalar is
nonconstant. Substituting \(r=m-1\) and \(W=2^a\) gives
\(\lceil(m-1)/d\rceil\le a\), equivalent to
\(d\ge\lceil(m-1)/a\rceil\).

This proves the all-degree claim without an external divisibility theorem,
quadratic-form classification, or enumeration. The finite all-functions
checks below are only controls.

## 6. Odd-repeat fixed-subtorsor restriction

Work on full antiperiodic words of length \(2m\). Imposing period \(2e\)
reduces the antiperiodic relation at shift \(m=eq\) to the relation at shift
\(e\), because \(q\) is odd. Thus
\(\operatorname{Fix}(T^{2e})\) is precisely the affine embedded copy of
\(X_e\): on the first \(m\) coordinates, blocks of length \(e\) alternate
between the block and its complement.

On this copy, the restricted shift is exactly the smaller twisted shift.
Each large input coordinate is either a small coordinate or that coordinate
plus one, so restriction cannot increase Boolean degree.

Equivariance gives \(R^{2e}F(x)=F(x)\), while the ambient target has
\(R^m=1\). Since \(\gcd(2e,m)=e\), the output has period dividing \(e\).
It is therefore an ordinary repetition of a length-\(e\) block. Repetition
preserves the primitive word, so target membership is exactly membership of
that block in \(B_e\).

Both source and target embeddings intertwine the respective shifts.
They also reflect orbit equivalence: if two embedded words are related by a
large shift \(k\), taking the first small block identifies them by the
corresponding small shift \(k\). Hence global injectivity on types yields
injectivity between the small source and target type sets.

When \(e\) is the dyadic part, Lemma 2 gives equality of their finite
cardinalities without importing an all-length existence theorem. Therefore
the restricted type injection is a surjection. Selecting the first \(e\)
output coordinates also cannot raise degree. This establishes the stated
general dyadic-part bound.

Oddness is indispensable: for the even-repeat hostile \(m=2e\),
\(T^{2e}=T^m\) is complementation and its fixed set is empty.

## 7. Fermat divisors and circuit scope

For any prime divisor \(p\) of \(2^{2^n}+1\), \(p\) is odd and the order of
\(2\bmod p\) divides \(2^{n+1}\) but does not divide \(2^n\). As all divisors
of \(2^{n+1}\) are powers of two, the order is exactly \(2^{n+1}\).
Consequently \(2^{n+1}\mid p-1\), giving
\(v_2((p-1)/2)\ge n\). The function
\((2^a-1)/a\) is strictly increasing for positive integer \(a\), since its
successive-difference numerator is \((a-1)2^a+1>0\).

Therefore the lower bound is at least \(\lceil(2^n-1)/n\rceil\), which tends
to infinity. Prime divisors exist for every Fermat number; no assertion of
infinitely many Fermat primes occurs. The differing exact orders also prevent
a fixed prime from supplying arbitrarily many indices.

For circuits over constants, XOR, and binary AND, degree is nonincreasing
under XOR relative to the maximum input degree and at most the sum under AND.
Induction on multiplicative depth therefore gives degree at most \(2^h\).
The displayed lower bound for \(h\) is correct. It does not apply unchanged
to arbitrary lookup gates, unbounded-fan-in AND gates, field primitives, or
arbitrary first-order descriptions, and the note does not claim otherwise.

The listed prime examples and integer lower bounds were independently
recomputed: \(17\mapsto3\), \(97\mapsto4\), \(193\mapsto7\),
\(641\mapsto11\), \(257\mapsto19\), \(65537\mapsto2185\).
These are not asserted exact minima beyond \(p=17\).

## 8. Independent exact checks

`independent_checks.py` imports no producer module. It uses explicit tuples
for shifts, the primitive-root definition via divisors, full forward source
orbits to evaluate inverse shifts, and the direct subset-parity formula for
ANF coefficients rather than a Möbius-transform routine.

Normal and optimized runs were byte-identical. They executed **249,374
explicit checks**, including:

1. Every source word for the \(m=2,4,6,8\) witnesses: target landing,
   equivariance, type bijection, full labelled surjectivity, exact degree,
   displayed images, and fibre counts.
2. Complete source and target orbit partitions at \(m=1,2,4,8,16\), including
   freeness, equal type counts, complement closure and target weight sums.
   At \(m=16\), both type counts are \(2048\) and the target weight sum is
   \(16376\).
3. Exact full fixed-subset equality at \((e,q)=(1,3),(2,3),(3,3),(4,3)\);
   both embedding intertwinings, target predicate preservation, and orbit
   reflection. Even-repeat cases are separately hostile.
4. Every Boolean function on one through four variables, checking Walsh
   divisibility at its actual ANF degree.
5. Every equivariant polynomial of degree at most one or two at
   \(m=2,3,4,6\), reproducing the complete small census, including the
   \(432\) quadratic bijections at \(m=6\).
6. A byte-identical independent reconstruction of all 29 packed basis
   columns and all 256 target classification entries in `input8.txt`.
7. Exact primality and degree-bound checks for the displayed table, plus
   order calculations for selected prime divisors of Fermat numbers through
   index six, including composite Fermat-number cases.

The sharp \(m=8\) witness is independently confirmed to have:

- seven ANF monomials and exact degree three;
- full scalar weight \(120\), quotient weight \(60\);
- full Walsh value \(16\), quotient Walsh value \(8\);
- \(117\) labelled outputs;
- \(112\) output fibres of size two, four of size four, and one of size
  sixteen.

At \(m=6\), the degree-two witness has \(28\) labelled outputs, with \(24\)
fibres of size two and four of size four. Its degree-one impossibility follows
both from the complete finite affine census and the inherited proved affine
obstruction. At \(m=4\), the displayed quadratic attains the lower bound.
These finite verifications prove each explicit witness but do not prove
generality beyond their tested universes.

### Independent full quadratic enumeration

`independent_exhaust8.cpp` is a separately written enumerator. It precomputes
the contributions of 15 lower and 14 upper coefficient bits, then visits all
pairs in ordinary binary coefficient order. It therefore does not reuse the
producer's Gray-code transition.

All **536,870,912** maps were checked against the independently reconstructed
input. The complete coverage histogram matched `exhaust8.json` byte-for-byte:

- target-valued maps: **22,136,841**;
- type bijections: **0**;
- maximum coverage: **15 of 16** target types;
- number attaining that maximum: **24,576**.

The producer C++ program was also recompiled and rerun, independently of this
second enumerator. The producer Python verifier completed its stated
**268,536** explicit checks, matched the stored receipt, and produced identical
normal and optimized receipts.

No cubic census or monomial-count optimality was checked or inferred. The
explicit cubic plus the proved lower bound establishes exact minimum degree
three, which is the claim actually made.

## 9. Publication and presentation notes

There is no mathematical blocker. The producer's `VALIDATION.md` still says
that a separate adversarial reviewer is pending. After attaching this review,
the parent may update that sentence to record this completed audit.

The main note's all-order proof is appropriately marked PROVED, while the
enumerations and implementation checks must remain FINITE-EXACT. The absence
of a historical priority claim is appropriate; this review does not establish
novelty or perform a new literature survey.

Review startup found the shared `math` checkout clean. The read-only
coordination peek could not run because `.machine-id` is absent; no
registration or identity file was invented. No mathematical check depended on
that optional coordination status. The shared checkout was still clean at
review close. Publication and any fetch/push reconciliation remain with the
parent, as requested.

## 10. Reproduction

From `/workspace/scratch/07d5532f2fce`:

```bash
python3 session13-degree-review/independent_checks.py
python3 -O session13-degree-review/independent_checks.py \
  --output session13-degree-review/independent_checks.optimized.json
cmp session13-degree-review/independent_checks.json \
    session13-degree-review/independent_checks.optimized.json

g++ -O3 -std=c++17 session13-degree-review/independent_exhaust8.cpp \
  -o /tmp/session13-independent-exhaust8
/tmp/session13-independent-exhaust8 \
  session13-degree-review/input8.independent.txt \
  > session13-degree-review/exhaust8.independent.json
cmp session13-nonlinear-degree/exhaust8.json \
    session13-degree-review/exhaust8.independent.json

g++ -O3 -std=c++17 session13-nonlinear-degree/exhaust8.cpp \
  -o /tmp/session13-review-exhaust8
/tmp/session13-review-exhaust8 session13-nonlinear-degree/input8.txt \
  > session13-degree-review/exhaust8.rerun.json
cmp session13-nonlinear-degree/exhaust8.json \
    session13-degree-review/exhaust8.rerun.json

python3 session13-nonlinear-degree/verify_degree_boundary.py \
  --output session13-degree-review/verification.rerun.json
python3 -O session13-nonlinear-degree/verify_degree_boundary.py \
  --output session13-degree-review/verification.optimized.rerun.json
cmp session13-nonlinear-degree/verification.json \
    session13-degree-review/verification.rerun.json
cmp session13-degree-review/verification.rerun.json \
    session13-degree-review/verification.optimized.rerun.json

sha256sum -c session13-degree-review/SOURCE_SHA256SUMS
```
