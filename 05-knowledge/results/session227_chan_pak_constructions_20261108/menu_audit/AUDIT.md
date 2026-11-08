# Independent audit: fixed finite graft menus and bounded intercepts

**Verdict: PASS.** The all-order finite-menu extension is correct. The sharp
ratio constants, finite-core argument, and inherited common-residue-gap
obstruction have been independently checked. No mathematical correction was
needed. One requested wording clarification makes the integer line-intercept
quantifier explicit. Chan–Pak all-odd completeness remains OPEN.

This audit was produced independently of the producer's code. It makes no
changes to the producer packet or shared repository. The independent checker
has no repository or third-party imports and uses explicit exceptions, so its
checks remain active under Python optimization.

## 1. Exact scope audited

The construction fixes all three parameters before counting:

- a finite bank of actual nonempty rooted tournament seeds;
- a finite menu of fixed nonempty rooted tournament grafts, with O and Z also
  available;
- an integer bound C on the final endpoint intercept e.

Every finite word is allowed, including adaptive choices, repetitions, the
empty word, and arbitrarily long words. Closing means O followed by J followed
by one additional J. The theorem restricts e measured at the OJ endpoint; it
does not bound all scalar counts from this grammar.

The finite-line consequence additionally requires a common odd residue gap
for the retained strict-slope lines. It does not say that every finite affine
line bank is incomplete. No graph inverse, physical realizability criterion,
strong-connectivity theorem, scalar image classification, or unrestricted
Hamilton-count coverage statement follows.

## 2. All-order argument review

### The physical cone and the exact ratio

For x=(D,E,U,d), the four coefficients

    E, d−E, D−d, U−D−d+E

are nonnegative on every physical state. Their combination of the four rays
in the producer proof reproduces x exactly. These are the previously audited
physical cone inequalities, not a new enumeration-based assertion.

Applying the graft's first and third rows independently gives

    D on rays: H, H+Q, Q, Q−H;
    U on rays: 2alpha, 2alpha+2T−s, 2T, 2T−alpha.

All four denominators are strictly positive because H≥1 and Q≥2H. Thus U'/D'
is a convex combination of the four ray ratios, with weights proportional
to cone coefficients times denominator values. The maximum is an upper bound
on every physical input.

Optimality has a separate justification. The first two rays are physical,
and the other two are approached by actual physical families from the
already audited exact-cone theorem. The quotient is homogeneous and continuous
at each of these rays because its denominator is positive there. Consequently
every ray quotient is attained or approached, and the maximum is the supremum
over physical inputs. This does not require every point of the cone to be
physical. The proof correctly avoids confusing a hull with realizability.

The singleton graft's exceptional fourth row is retained. The first and third
rows used in the ratio argument remain valid, yielding kappa_J=2. Direct ray
evaluation gives kappa_Z=2. O has a zero denominator and positive numerator
on the fourth ray, and actual transitive source-rooted inputs show that no
finite bound is available. The named values kappa_L=4, kappa_C3=2, and
kappa_B4=2 are correct. A transitive source-rooted gadget of order n applied
to a singleton gives D'=1 and U'=2^n, so no menu-independent constant is
silently being used.

### Finite cores, uniformity, and exact closing counts

The OJ endpoint has state

    (D+2U, D, D+2U, D),

and the extra closing J has count N=D+4U. Hence h=D+2U, g=2U, and e=D.
In particular h>g>0, g is even, and g/h≥2/3.

A word consisting solely of O's yields its seed pair (D_s,U_s) followed by
a dyadic tail. Otherwise split after its last non-O operation F. The output
there is physical, D_1 equals the final intercept, and

    1≤U_1≤kappa D_1≤kappa C.

Trailing O's preserve D_1 and double U_1. Thus all counts lie among

    D+4u·2^k, k≥0,

with 1≤D≤C and 1≤u≤floor(kappa C), together with the eligible seed pairs.
The displayed number of possible pairs

    B_C=C floor(kappa C)+|S_C|

is valid even with repeated pairs, since repetition only enlarges the bound.
Each pair contributes at most 1+floor(log_2 X) counts through any real X≥1.
The more precise count with (X−D)/(4u), including its zero-contribution case,
is also correct.

This is uniform in word length. It is not a finite-depth argument, and it
needs no bound on the number of non-O letters. Monotonicity of D follows
from insertion into an induced extension and justifies prefix pruning, but
is not needed to replace the last-operation bound by a word-length bound.
Singleton seeds, arbitrary physical seed orders, and sink roots are included.
For a sink root x=(D,D,D,D), Zx=x; arbitrarily long stationary histories are
therefore harmless. All-O seed tails are essential: physical seed states with
D_s=1 and U_s larger than kappa C exist.

### Transfer of the residue obstruction

For a line a h−b g=c with coprime positive integers a<b, substituting N=h+g
gives

    g=(aN−c)/(a+b), h=(bN+c)/(a+b).

Integer h and even integer g are equivalent to aN≡c modulo 2(a+b). This is
exact only as an arithmetic candidate test, not a physical reachability test.
For a/b<2/3, the endpoint cone forces

    h≤floor(−3c/(2b−3a)),

so the finite ceiling N≤2 floor(−3c/(2b−3a))−1 is valid when there are positive
candidates. The producer's conservative bound on the number of excluded
counts is safe.

The lcm modulus is even. A common odd residue avoiding every retained line
therefore supplies X/m+O(1) odd integers through X before boundary and finite
exceptions are removed. The bounded-intercept theorem removes at most
O(log(X+2)) further counts. The claimed lower bound follows with constants
depending only on the fixed bank, menu, and lines. The explicit representative
formula for progression counting is also correct, including X below the
first positive representative.

The five lines 2h−3g=c, c in {0,−2,−4,−6,−8}, are a valid hostile: exactly one
line matches each odd residue, and every odd N≥17 has a positive chamber
candidate with even g on its matching line. This genuinely prevents dropping
the common-gap hypothesis. Neither the producer nor this audit claims that
those candidates are reachable.

## 3. Independent finite controls

The standalone `independent_check.py` was written without importing or
copying a producer checker. It uses direct permutation/two-bin counts for
source states, direct permutation/bin-cut counts for gadget statistics,
and a separately implemented memoized deletion recursion for literal graphs.
Its finite universe is explicit:

- Every labelled tournament through order 5, with every marked root:
  5,405 presentations and 46 distinct profiles. There are 361 sink-root
  presentations, with their equality case checked directly.
- Every rooted gadget through order 4: 285 presentations, checked against
  all 46 source profiles, giving 13,110 ratio controls.
- Every rooted source through order 3 and every rooted gadget through order
  4: 8,265 literal graft pairs, with no source, gadget, or root filtering.
- 87 literal O, Z, and final OJJ controls.
- 14 literal extreme-ray family profiles and 3,420 exact rational limiting
  quotient controls, using orders and scales recorded in the checker.
- A different eight-operation menu from the producer's larger-ratio menu:
  O, Z, singleton graft, both roots of an edge, a directed triangle, and
  B4 at roots 0 and 1. Its common ratio constant is 4; the second B4 root
  has constant 12/5, so a nonintegral individual constant is exercised.
- Six physical seed states, including a nontrivial sink and the source-rooted
  transitive order-9 state (1,0,256,1). Every word through length 5 is checked,
  totaling 224,694 presentations. At C=25, there are 2,883 bounded-intercept
  word presentations and 182 distinct counts in this finite word universe.
  The pair envelope bound is 2,506.
- Separate controls include 10,000 stationary Z steps, exact O tails through
  exponent 1,000, and all-O seed tails that exceed the non-O core cap.
- 260,145 exact even-g congruence checks, 372 positive discarded-line
  candidate checks, and 1,000 odd-count controls for the five-line hostile.
  A concrete common-gap bank and exact progression counts are checked too.

These finite controls supplement, and do not prove, the all-order statements.

## 4. Producer replay and reproducibility

The producer checker was reviewed and independently replayed. It tests all
5,405 rooted gadgets through order 5; 841 literal source/gadget pairs with
all roots through order 3; 374,490 mixed-menu words; and its additional
boundary, stationary-run, ray-family, and arithmetic controls. Its larger
mixed menu includes the source-rooted transitive order-5 gadget and has
common ratio constant 32. This is intentionally a different finite control
menu from the independent check above.

Both checkers passed in normal and optimized modes. Each pair of resulting
receipts is byte-identical. Producer replay receipts also equal the producer's
checked-in receipts. The audited source and receipt hashes are recorded in
`audited_inputs.sha256`; this packet's own files are covered by `SHA256SUMS`.

Reproduce the independent controls in this audit directory:

    python3 -I -B independent_check.py > replay.json
    cmp replay.json independent.json
    python3 -I -B -O independent_check.py > replay.optimized.json
    cmp replay.optimized.json independent.optimized.json

`producer.json` and `producer.optimized.json` preserve the separate producer
replays. No finite success is promoted to an unrestricted completeness or
impossibility result.

After assembly as sibling `menu_boundary` and `menu_audit` directories,
replay the audited producer from this directory with:

    python3 -I -B ../menu_boundary/verify.py > producer.replay.json
    cmp producer.replay.json producer.json
    python3 -I -B -O ../menu_boundary/verify.py > producer.replay.optimized.json
    cmp producer.replay.optimized.json producer.optimized.json
    sha256sum -c audited_inputs.sha256
    sha256sum -c SHA256SUMS
