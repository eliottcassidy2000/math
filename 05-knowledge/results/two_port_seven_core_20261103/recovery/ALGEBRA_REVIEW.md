# Independent review of B1-B5

Status: **ACCEPTED — PROVED**, with finite controls separately scoped.
Independent recovery/audit lane,2026-11-03. Reviewed
[algebra/PROOF.md](../algebra/PROOF.md), its complete checker, and the
inherited parity, uniform-run and multi-in-port identities. No blocking
defect found. One independent replay matches the frozen receipt exactly;
the normal and optimized producer receipts also agree.

## B1: coefficient divisibility

Every Hamilton word uses at most two old-port incidences per port, and
cannot use the same incidence twice. Its product of attachment bits and
their complements is therefore already Boolean multilinear, of total
degree at most2m. This degree bound concerns the complete response,
including its three-path terms.

For a selected monomial, the pointwise stabilizer of its mentioned old
labels fixes that monomial while merely permuting variables outside it.
Its coefficient in each word's product is consequently constant on
orbits, including the zero coefficient. The action on fully named words
is free: a fixed word fixes every old label. Signed coefficients therefore
sum in multiples of the stabilizer order. There is no assumption that
all coefficients are positive. Taking the gcd of the orders over all
sets of at most2m labels indeed gives a common divisor of all values.

## B2: exact triangle content, including zero boundaries

Independent triangle rotations leave at least s-2m whole components
untouched by any monomial. This establishes the3^max(s-2m,0) lower bound.
The all-source attachment has value3^s, excluding every other prime from
the content. The local first difference2 and mixed difference-10 are
verified by literal Hamilton permutations of the four- and seven-vertex
extensions, with the specified mask convention.

The block construction is a genuine ordered join: earlier triangles
point into a block's port, and its port points into later triangles;
ports are ordered by their blocks, with unused ports universal sinks.
The number of assigned ports is ceil(min(s,2m)/2), at most m. Thus the
restricted mixed coefficient factors as

    3^(s-t)(-10)^floor(t/2)2^(t mod2),  t=min(s,2m).

It is an integer alternating sum of actual query values, so the content
divides it. Together with3^s its gcd is precisely the stated power of three.
For a restricted polynomial, coefficients can collect contributions from
variables frozen to one, but the argument needs only their expression as
integer query combinations and is unaffected. In the highlighted m=2,t=4
case, total degree four also identifies the full quartic coefficient.

The cases s=0 and m=0 are correctly covered: an empty old tournament with
transitively ordered ports has one Hamilton word; zero ports give exactly
3^s. A common factor is not a collision or a pair of proportional response
tables. The text correctly rejects only the asserted universal primitive
content shortcut.

## B3: fixed-slot endpoint descent

For k>=3, permuting the other k-1 slots freely makes each fixed-slot
endpoint count divisible by(k-1)!. For k=2, the squarefree F2 identity
e=Us, P=U-1 and U^2=1 gives e_uP=s_uP. This is the required evenness
of E_k(u)-B_k(u), with both paths nonempty and disjoint; ordinary endpoint
parity alone would not justify it.

The ambient skeleton coefficients are integers. In the endpoint channel,
the distinguished M-run is the last one; in the cross-derivative channel,
it is the unique run immediately before the chosen exterior vertex w.
Move that distinguished run to the first fixed slot and keep the remaining
runs in their ambient order. This gives a bijection without division by a
factorial, and its coefficient is independent of the named endpoint u
because all module cross incidences are uniform.

For the cross channel, differentiating u->w removes exactly that factor.
It need not be an existing arc; this is why the one-run coefficient is
the full quotient derivative C_zw(Q). All other required quotient arcs
remain in the exterior skeleton. At least one usable leading coefficient
is odd because

    E_z(Q)+sum_(actual z->w) C_zw(Q)=H(Q).

Every Hamilton word either ends at z or uses one of its actual outgoing
arcs. This identity also covers a one-vertex quotient, where E_z=1 and
there are no cross channels. Subtracting the converse expansion kills all
higher-run terms modulo two and yields the claimed intrinsic endpoint
parity condition. The module and named-arc hypotheses are essential and
explicit; no arbitrary-subset converse claim is made.

## B4: ordered-cover skew descent modulo four

The unordered block partition formula has the correct factor k! for
ordered slots. A directed-arc derivative contributes only in the unique
block containing both of its endpoints. The uv-minus-vu derivative of
that block is even by the previously proved parity symmetry of C on
every induced tournament. Hence the full difference is divisible by2k!,
and in particular by four for k>=2.

The one-run multiplier is H(Q), odd. Internal converse transposes the
ambient gradient for all module sizes, not just odd sizes. Equality of
the two ambient C matrices makes that gradient symmetric and forces the
intrinsic module skew to vanish modulo four. Multiplication by the odd
quotient count is legitimate in this residue ring. The claimed order-eight
hostile has C_02-C_20=150, nonzero modulo four, so its endpoint parity
balance cannot hide its converse inside any uniform ambient quotient.

Both congruences are correctly stated as necessary conditions. The proof
does not claim that jointly satisfying them suffices for a collision, or
that all larger even nonselfconverse cores violate at least one.

## B5: physical ports and exact zero-readiness

The imported Session181 defect is H(A)J_W, with
J_W=sum_(u,v in W) C_uv(B). Diagonal entries are zero. Each unordered pair
contributes C_uv+C_vu=2w_uv, where the inherited cut theorem supplies
positive integer w_uv. Thus

    J_W >= |W|(|W|-1),

with equality exactly when every pair weight in W is one. The nonempty
prefix has positive Hamilton count. For the stipulated nonempty W,
zero defect is therefore equivalent to a singleton W. No strongness,
clock truncation or bounded suffix census is used in this implication.
The result concerns that particular multi-in-port construction, as stated.

## Computation and remaining scope

The checker uses arbitrary-precision integers and explicit exceptions.
The ambient controls contain889 pairs,9,564 exact fixed-slot channel
comparisons and14,612 internal skew entries, with up to four old module
runs. Complete response tables comprise4,750 actual queries in the
eleven specified(s,m) cases; the additional factorization witnesses test
the larger sharpness coefficients directly. The physical-port identity
is checked in126 suffix/prefix cases. These are controls of the all-order
arguments above, not their extrapolation.

The independent replay, source match, source stability and agreement of
normal/optimized expectations are recorded in
[algebra_review_receipt.json](algebra_review_receipt.json). The complete
three-path response is retained throughout. Universal two-port
reconstruction remains open, and none of B1-B5 supplies a six-colour
unit-distance obstruction.
