# A reusable exchange interface with a third-cover response

**Status: PROVED by the finite compressed-run certificate below, with literal
controls.** The predecessor's count kernel is INHERITED. Arbitrary target
coverage and a post-birth return to this source class remain OPEN. No new
numerical-spectrum or external-priority claim is made.

## 1. The finite example becomes an all-module operation

Use the inherited quotient

    Q=(198,196,136,211,167,207,20,64),  ports 5 and 7.

Rows are outgoing bit masks. Substitute arbitrary nonempty tournaments A,B
at ports 5,7 and retain the exterior vertex 0 as a singleton mark. Write

    h=H(A), c=C₂(A), f=C₃(A),
    g=H(B), d=C₂(B), e=C₃(B),

where Cᵣ counts ordered covers by r nonempty directed paths. The inherited
[bilinear kernel](../../session261_creative_chanpak_20261111/constructors/BILINEAR_KERNEL.md)
is

    K(A,B)=149hg+72(hd+cg)+37cd.                         (1)

Thus exchanging A,B preserves the actual source Hamilton count. The new
rooted cover response, with first path allowed empty and second starting
at the retained vertex 0, is

    U₀(A,B)=430hg+403hd+120he+202cg+196cd+62ce.          (2)

Every coefficient is nonnegative. The third-cover coordinate of the SECOND
module has become visible to the root. The formula is valid for arbitrary
module orders, not only chains, and includes the initial cases with C₂ or
C₃ zero.

Retain also E₀, the number of source Hamilton paths ENDING at0:

    E₀(A,B)=14hg+7cg=7Q₂(A)H(B).                      (3)

Adjoin a new vertex z with z→0 and every other source vertex →z. Splitting
a Hamilton path at z gives exactly

    H(Q[A,B]+z)=K(A,B)+U₀(A,B)−E₀(A,B).                (4)

The z-last paths contribute K−E₀: an old path ending at0 cannot append z.
Every other path is uniquely P,z,R, where R is nonempty and starts at0;
P can be empty. Every such cover is legal, since the last vertex of a
nonempty P cannot also be0. This proves (4) directly.

For t≥1, put t exterior unselected sinks between the source and z,
ordered transitively, retaining z→0 and all other vertices→z. Now the
inherited selected-tail law gives

    N(A,B,t)=K(A,B)+2^t U₀(A,B),   t≥1.                (5)

The empty tail is NOT obtained by substituting t=0 into (5); its correction
is (4). This initial-case distinction was detected by the new literal
constructor control before publication. For positive tails, increasing t
reuses the same source and gives N→2N−K. The zero-to-one transition instead
has increment U₀+E₀. The base offset is never freely selectable.

## 2. Why the certificate proves all module sizes

Compress every maximal run of A or B along a path to its port symbol.
All six exterior symbols occur once; equal port symbols cannot be adjacent.
A skeleton with r A-runs and s B-runs expands in Cᵣ(A)Cₛ(B) ways: order the
local paths by their occurrence in the skeleton. This is a bijection, so
its coefficients are nonnegative integers.

Between successive occurrences of a fixed port there must be an exterior
vertex. Otherwise the intervening string uses only the other port and
would require both directions of the A/B quotient edge. With k exterior
vertices this proves r,s≤k+1. The ordinary source has k=6; the born quotient
has k=7. Hence a finite enumeration of all skeleton words is complete at
all module sizes. The supplied code uses the harmless larger cutoff k+2.

`response_kernel.py` enumerates these words by the finite state

    (exterior subset, A-run count, B-run count, last symbol).

A transition appends an unused exterior symbol or another port run when
the quotient arc permits it. Starts, termination and all arc tests are
explicit. The complete nonzero born coefficients are

    (r,s)       (1,1) (1,2) (1,3) (2,1) (2,2) (2,3)
    Hborn         565   475   120   267   233    62.

Subtracting the inherited (and independently enumerated) four coefficients
of K and adding the endpoint coefficients (14 at(1,1),7 at(2,1)) gives(2).
The endpoint coefficients use the same exhaustive skeleton DP, retaining
only skeletons whose last symbol is0. The subtraction has the positive
cover interpretation in(4). All other coefficients are zero because ALL allowed
run numbers have been enumerated, not because a few module lengths fit a
formula. Sixteen literal substituted tournaments test the expansion,
including nontransitive C₃ modules.

## 3. Exact balance for the positive ray and its finite head

Let Δ=U₀(A,B)−U₀(B,A). From (2),

    Δ=201(hd−cg)+120(he−fg)+62(ce−fd).                 (6)

Define L(M)=201H(M)−62C₃(M) and R(M)=201C₂(M)+120C₃(M). Then

    201Δ=L(A)R(B)−R(A)L(B).                           (7)

Thus an exchange preserves the positive-tail ray exactly when the actual
profiles have the same projective [L:R] coordinate. This includes the
singleton boundary R=0 without division. The skew pairing has null vector
(62,−120,201) in (H,C₂,C₃) coordinates. Formal points on that null line do
not manufacture physical modules.

There is a sharper criterion when the bare-apex initial case must ALSO
survive. Since K is already symmetric, preserving both the t=0 output and
any positive-tail output is equivalent to preserving U₀ and E₀. Equation(3)
first gives hd=cg. Substituting in(6), and dividing by hg>0, leaves

    (120+62c/h)(e/g−f/h)=0.

Its first factor is positive. Therefore the complete pointed-tail response
including its finite head is exchange-invariant if and only if

    C₂(A)/H(A)=C₂(B)/H(B),
    C₃(A)/H(A)=C₃(B)/H(B).                             (8)

This is an exact necessary-and-sufficient response guard on actual modules,
not a scalar count coincidence. The inherited six-vertex full-profile pair
in the constructor lane supplies nonidentical physical modules meeting(8),
with different future rooted menus; no new pair search is needed to inhabit
the criterion. After arbitrary inflation those menus may cease to match,
so such inflation is not asserted to return to the guard.

The universal kernel has a genuine third-cover coefficient. This is not a
claim that C₃ is functionally independent of (H,C₂) on every restricted
family. An exploratory collision search through labelled order six found
no same-(H,C₂), different-C₃ example; that observation is not used in any
proof or all-order exclusion.

## 4. All-size chain supply, with exact initial objects

For TTₐ, C₁=1, C₂=2^a−2 and C₃=3^a−3·2^a+3. Substitution gives, for a,b≥1,

    K(a,b)=37·2^(a+b)−2^(a+1)−2^(b+1)+9,
    U₀(a,b)=62·2^a3^b+10·2^(a+b)−4·2^a
             +23·2^b−4·3^b−8.                        (7)

Together with (4) and (5), this is a literal three-parameter constructor. Its
initial source Q[1,1] has (K,U₀,E₀)=(149,430,14), with bare-apex count565
and first positive-tail count1009.
The first unequal pair is

    Q[TT₁,TT₂]: (K,U₀,E₀)=(293,1236,14), bare output1515;
    Q[TT₂,TT₁]: (K,U₀,E₀)=(293,834,28),  bare output1099.

No negative-length module or absent-module extrapolation is used. Their
positive-tail rays are 293+1236·2^t and 293+834·2^t, t≥1; their
first values are2765 and1961. Their separate bare outputs are as displayed.

**Strict response ordering.** If b>a≥1 then U₀(a,b)>U₀(b,a). Indeed

    Δ=62(2^a3^b−2^b3^a)−4(3^b−3^a)+27(2^b−2^a).

Putting d=b−a≥1 gives

    2^a(3^d−2^d)≥3^d−1.

It suffices to take 2^a≥2 and use 3^d−2^(d+1)+1≥0, true at d=1 and
increasing thereafter. Thus Δ≥58(3^b−3^a)+27(2^b−2^a)>0. So the finite293
phenomenon persists at EVERY unequal pair, while the diagonal is exactly
balanced within this chain grammar.

The generators are actual graph operations:

* L: add one transitive twin inside A;
* R: add one transitive twin inside B;
* E: exchange the two module graphs while retaining the quotient labels;
* P: add one exterior tail vertex.

L and R commute; E² is the identity and E L E=R. Source exchange preserves
K at every state. P preserves the source and its marked interface. A
state (2^a,3^a,2^b,3^b,t) determines every displayed observable and updates
exactly; the 3^a and 3^b channels cannot be dropped merely because K itself
uses only powers of two. After a born output, this note does not supply an
operation returning the WHOLE output to a Q-source.

The quotient is strong (a direct reachability check on Q), and replacing
vertices of a strong quotient by nonempty modules preserves strong
connectivity. The output is also strong: z reaches the base through0,
the strong base reaches every tail vertex and z, and each tail reaches z.
All guards therefore survive arbitrary positive a,b and nonnegative t, with the zero-tail branch treated separately.

## 5. What the universal defect inequality adds

At every new apex, the inherited inequality gives the exact constraint

    H(output)=H(base)+d(G_z)−1+σ,  σ≥0.

For the t=0 birth here, H(base)=K and the increment is the ACTUAL U₀−E₀.
Hence d(G_z)−1≤U₀−E₀. Formula (2) provides the missing interface information
that the inequality alone cannot select. The constructive task is not to
make domination large; it is to route a target into a legal tuple (A,B,t)
with its coupled K and U₀, or find a new physical return that changes this
response while paying the base count. The exchange changes that response
without changing K, uniformly and with known initial objects. It gives
supply and a balance equation, not arbitrary independent increments.

## 6. Recorded historical sources and actual transfers

The four paths and reproducible random seed are in `random_sources.json`.
They are idea inputs at their stated status; no historical claim was
re-audited or promoted here.

* **Mills shell decoder:** its lost fractional remainder suggested retaining
  a missing response coordinate. Map: source module → (H,C₂,C₃). Preserved:
  exact rooted continuation. Lost by (H,C₂): the third-run coefficient of
  this universal polynomial. Decisive transfer: equations(2)–(5), not a claim
  that prime-generating constants model tournaments.
* **Finite O/Z involution rigidity, THM-4888:** a count symmetry need not
  preserve a tail observer. Map: Q module exchange → alternating pairing
  (6). Preserved: K; sidecar: projective response [L:R]. Decisive test:
  initial293 exchange and the all-a,b strict ordering theorem. Its O/Z
  formulas are not imported into this different graph grammar.
* **THM-1370's omitted7,21 boundary:** finite initial cases must be part of
  the constructor, not supplied by an ambient odd-number ring. Map:
  abstract profile domain → actual positive module lengths. Preserved:
  literal tournaments and exact counts; excluded: empty-module extension.
  The initial table149,293 is a complete initial specification for this
  grammar, not a new spectrum census.
* **THM-102's filling/cancellation discussion:** its global-rank versus
  local-projection distinction suggested retaining an entire skew form,
  rather than one response difference. Actual transfer is the rank-two
  determinant identity(7). No path-homology claim is a proof dependency.

## Reproduction

    python3 -B core/response_kernel.py
    python3 -B core/construct_exchange.py --a 1 --b 2 --tail 0

Run from this packet, or use each repository-relative script path. The
coefficient certificate and sixteen literal controls are recorded in
`response_kernel_receipt.json`; these do not prove general Chan–Pak.
