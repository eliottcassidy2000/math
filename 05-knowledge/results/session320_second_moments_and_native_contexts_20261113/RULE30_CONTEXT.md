# Rule 30 has exact three-way continuation; homogeneous tournament dominance does not

## Native setting and provenance

This concerns one **spatial** Rule 30 output word on a periodic ring. It says nothing about the temporal center of an isolated seed. The inherited [left-permutive trace compiler](../../../01-canon/theorems/THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary.md) and [inverse de Bruijn reset theorem](../../../01-canon/theorems/THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md) supply four boundary states in order `00,01,10,11`. The two inverse matrices are column-deterministic: their unique nonzero row in each column is respectively

    g_0=(0,2,3,3),     g_1=(2,0,1,1).

Products use transformation composition, and the number of periodic predecessors of output word w is the trace of its product M_w. By trace cyclicity, arbitrary left and right contexts reduce to a single suffix probe. A pair is **context-equivalent** if every suffix probe has equal trace. This is equality of the full integer readout, not just the inherited 13-state test for “exactly one predecessor.”

## Same-depth sign balance for every left-permutive binary rule

For Rule 30 put D=A_0+A_1. It is the adjacency matrix of the full order-two de Bruijn graph; D²=J, the 4×4 all-ones matrix, and D^m=2^(m-2)J for m≥2. Each M_u has column sums one, so

    sum_{|s|=m} tr(M_u M_s) = tr(M_u D^m) = 2^m,  m≥2.

Thus any trace discrepancy at a fixed suffix depth m≥2 must have a discrepancy of the opposite sign at another suffix of exactly that depth. More generally, for a binary left-permutive rule with boundary memory d, use the full de Bruijn graph of order d, with D^d=J of order 2^d, to get the same conclusion at m≥d. This proof does not require a synchronizing word.

Rule 30 supplies an exact equality probe too. By the inherited cyclic reset criterion, any spatial cyclic word containing 010 has exactly one predecessor. Thus the suffix 010 followed by arbitrary padding gives equality at every depth m≥3, independent of u and v. The replay verifies that 010 is the unique universal equality probe at length three and that there are none shorter.

## The complete count interface and exact-depth separation

The 41-state inverse monoid has exactly 23 context classes, or 22 when restricted to nonempty words. The four suffix probes

    empty, 0, 1, 01

already distinguish all classes, and no choice of three probes of any lengths can do so, even on the nonempty states alone. The bound on maximum probe length is sharp: 01 and 10 agree on every probe of length ≤1, but disagree at 01.

More strongly, **at every fixed depth m≥3**, probes of length exactly m separate all 22 nonempty classes. Consequently every context-distinct nonempty pair has negative, zero, and positive discrepancies at every depth m≥3. This combines exact-depth separation, same-depth balance, and the 010 equality probe; none alone implies the three-way conclusion.

The finite certificate is fully reproducible. Close the two generators by right composition starting at the identity; the 41 products have ranks distributed 4,30,6,1. Tabulate the traces of pairwise products. Test the four-probe partition against the full table; test the at most binomial(41,3) possible triples of probes. If S_m is the set of exact-length m products, its sizes for m=0,…,11 are

    1,2,4,8,16,24,30,34,36,36,36,36.

Check the set identity S_8=S_11, which propagates under the common one-letter update. The exact-depth partitions on nonempty states agree with the full table at m=3,…,10, exhausting the initial range and one full future period. These are finite exact certificates for all suffix depths, not fits to growing word lengths. The empty state is intentionally excluded from the exact-depth assertion.

## A sharp pair with exact information-decay rate

For u=01 and v=10 every discrepancy belongs to {-1,0,1}. For each nonzero sign, the number a_m of length-m suffix probes of that sign has generating function

    sum_m a_m z^m = z² / [(1-z)(1-z-2z³)].

In particular a_0=a_1=0, a_2=1, a_3=2, a_4=3. Both signs have the same counts by same-depth balance (their magnitudes are one), while their proportion decays exponentially. The dominant growth constant is the positive root rho of x³=x²+2, the same constant appearing in THM-4204's distinct unique-fibre language. The transfer certificate evaluates the sign indicator f on the monoid, with T f(a)=f(a g_0)+f(a g_1). For each sign form r=(T⁴-2T³+T²-2T+2I)f. The replay verifies (T^k r)(identity)=0 for 0≤k<41. Cayley-Hamilton makes this valid for every k, proving the displayed rational series from the four initial coefficients. Neither equality of the growth constant nor the transfer recurrence says anything about the isolated seed's temporal center.

## Precise tournament comparison

In the [one-arc reservoir switch](../session319_arc_usage_and_signed_cycles_20261113/RESERVOIR_POLYNOMIAL.md), equal Hamilton counts coexist with a strict one-sided order in the ordered nonempty path-cover coordinates. Every common homogeneous continuation reads these coordinates with nonnegative itinerary coefficients. A strong quotient exposes the two-run difference; after that, equality cannot return within this common homogeneous grammar.

Rule 30 has the **opposite geometry of probes**: equally long contexts average to a constant independently of the incoming boundary transformation, so there is no strictly positive contextual order between distinct persistent count classes. A separator can distinguish in either direction, and a physical reset provides equality at exactly the same depth. This contrast points to a concrete tournament design target rather than a numerical analogy: a common continuation cannot recover equality from componentwise dominance, so a reusable signed repair must first create *incomparable cover profiles*, or change the context with the target. The Rule 30 reset is not itself a tournament construction.
