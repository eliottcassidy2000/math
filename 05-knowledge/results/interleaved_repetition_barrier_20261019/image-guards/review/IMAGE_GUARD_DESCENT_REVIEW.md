# Independent review: image guards and periodic-tail endpoint descent

Date: 19 October 2026.

## Verdict

**PASS, within the stated expansion-first first-crossing bank language and
conditional on the independently audited matched-block height theorem and its
cited logarithmic-form input.** No mathematical obstruction was found in the
frozen packet. Its new conclusion is actual strict endpoint descent, not merely
a conditional cycle obstruction, for the rational-bracket family in Section 5
and consequently an unbounded family indexed by upper continued-fraction
convergents. The exact enormous example belongs to that excluded family.

The audit does not claim arbitrary bank entry, exclusion of all Collatz cycles,
general Collatz convergence, or a complete description of any realized S_d
boundary. The frozen packet's separate eventual typed-factor question is not
used in any of its proved statements. This review concerns the assigned frozen
packet, not subsequent extensions of that question.

The principal new arithmetic fact is particularly robust: a lawful V output
is 82 modulo 243, whereas a lawful VF output is 1 modulo 243. The selected
matching checkpoints are therefore unequal. The already-audited matched-block
theorem turns that inequality into strict descent at the actual endpoint.

## 1. Frozen material, dependencies and reproduction

The primary assigned hash is verified:

    IMAGE_GUARD_DESCENT.md
    8e32d46d23eceefee6a45b2bc7c00548d5409cc14dda99171298d0d56dce3c84

Every file size and SHA-256 digest listed in the packet's MANIFEST.json matches.
Its SHA256SUMS also verifies. The hostile companion has hash
4608853b0348341174dd1e98b8ba1e5fbf665c8add2601d2d165c4f0aadacd48.

Both ordinary and optimized runs of the submitted verifier reproduce their
receipts byte for byte: **30,113 explicit checks**. No Python assertion is being
silently removed by optimization; the submitted checks use explicit exceptions.

The independent program in this review imports no submitted code. Both its
ordinary and optimized runs agree byte for byte: **238,861 explicit checks**.
It uses a different logarithm basis and series truncation, a reversed closed-sum
affine-coefficient implementation, actual parity/orbit replay, rational model
controls, and independent finite-word continuation tests. The universal
quantifiers are justified below, not by these finite check counts.

The local primary INTERLEAVED_REPETITION_BARRIER.md is byte-identical to the
published reviewed version. The greedy dependency is also byte-identical.
The matched-block dependency differs only in its introductory review-status
sentence: the published version replaces the staging statement with
"Independent audit: PASS." Its mathematical text is unchanged, and its source
and published hashes match the published review/provenance record. Thus the
packet uses the theorem that was actually audited, including its arbitrary-bit
start extension.

This review inherits rather than re-proves the primary logarithmic-form theorem
application. It checks that all its hypotheses are retained in the new use.

## 2. Exact bank output types

For a chronological Terras word of length n and h odd bits, write its branch
as (3^h x+b)/2^n. An independent closed-sum reconstruction gives:

| Bank letter | h | b | M=2^n | v_3(b-M) |
|---|---:|---:|---:|---:|
| K | 19 | 4912087879 | 1073741824 | 3 |
| V | 77 | 83829391903560839836419813590128166431 | 5316911983139663491615228241121378304 | 4 |
| H | 4 | 101 | 512 | 1 |
| F | 1 | 1 | 4 | 1 |

For K, V and H, the numerator of y-1 is 3^h x+(b-M). Since h exceeds the
listed finite valuation, the first summand has strictly larger valuation for
every integer x. Cancellation at the lower valuation is impossible. Division
by M does not change a 3-adic valuation. This proves the exact type uniformly
on the lawful integer cylinder, including inputs divisible by arbitrary powers
of three. Sampled outputs are not the proof.

F requires the separate identity 4(F(x)-1)=3(x-1). Its type increases by one
when finite; input 1 remains the fixed point 1 and has infinite type. The
packet explicitly retains this exception. It does not incorrectly apply the
strict-valuation argument when h=v_3(b-M).

The resulting reset-and-increment rule is a typed counter with potentially
unbounded values, not a finite arithmetic quotient. The text preserves that
distinction.

## 3. Past images, future cylinders and the CRT gap

For every lawful past suffix A, the exact identity

    2^n y=3^h z+b

implies y=b 2^(-n) modulo 3^h. The unknown state z may itself follow an
arbitrary earlier context: the suffix identity is enough, and no earlier
context is discarded.

For two actual checkpoints admitting the same next r parity bits, the
odd-coefficient cylinder theorem gives y_i-y_j=2^r k for an integer k. If
their certified past residues differ modulo 3^s, then k is nonzero and is
congruent to

    a=(rho_i-rho_j)2^(-r) modulo 3^s.

For 0<a<3^s, the least absolute integer in that residue class is exactly
min(a,3^s-a). Negative differences are included. Consequently

    |y_i-y_j| >= 2^r min(a,3^s-a).

There is no distribution assumption here. The independent controls explicitly
exercise negative as well as positive representatives.

For unequal finite types t_i,t_j, factoring the lower power of three out of
(y_i-1)-(y_j-1) proves v_3(y_i-y_j)=min(t_i,t_j); this supplies a simple
special case of incompatible past-image residues.

Under the non-descent supposition, the audited first-crossing height estimate
strictly bounds each positive checkpoint value and therefore its nonzero
difference from another positive value. Combining that strict upper bound
with the CRT lower bound produces the packet's strict logarithmic inequality
(5), with +122 log 2 at arbitrary bit starts. At macroblock starts the audited
bound has no such factor and the term may be removed. Replacing full L by a
prefix expansion length means using a prefix long enough to dominate both
checkpoints, for example the later checkpoint's prefix.

For the convenient condition 100r>=L, the audited short/long overlap makes
non-descent with unequal checkpoint values impossible. The packet removes
the exact-return alternative by an actual integer image contradiction, not
by comparing formal multipliers alone.

## 4. The V/VF separation

The independent arithmetic reproduces

    b_V=103 mod 243, M_V=22 mod 243, 22^(-1)=232 mod 243,
    rho_V=82, rho_VF=1.

The VF formula is computed by exact affine composition and also follows from
(3*82+1)/4 modulo 243. Its finite output type is 5, while V's is 4.

Their difference 81 times a unit remains either 81 or 162 modulo 243. Indeed
only that unit's residue modulo 3 matters, so the statement holds for every
nonnegative integer r, not just for the finite r values used in a control.
The resulting gap is at least 81*2^r. In particular, these checkpoints cannot
be an exact return even if their symbolic future continuations never diverge.

## 5. Rational-floor construction and the final reset

Assume the displayed hypotheses, write delta=p/q-theta, and put

    H_n=floor(np/q)-1_{q divides n}.

The inequalities m>2q and delta<1/((m-q)q) imply
m delta<2/q<=1. For 1<=n<=m, classify k=np mod q:

- If k>=2, the displacement n delta is smaller than k/q, so no integer
  boundary is crossed.
- If k=0, a positive displacement less than 1 lowers the floor exactly once,
  as recorded by H_n.
- If k=1, invertibility of p modulo q and mp=1 modulo q show that n=m modulo
  q. Every earlier such n is at most m-q, where n delta<1/q. At n=m the lower
  bracket gives m delta>1/q, producing precisely one additional loss.

This proves the floor identities for every n in their full stated range.
No sampled-floors extrapolation is involved.

Because 0<p<q, the increments g_n=H_n-H_(n-1) are binary and q-periodic with
p ones in a period. H_0=-1 gives g_1=1. At the period boundary, the ordinary
rational mechanical sequence's 10 pair becomes 01. Thus the first actual
greedy indicator f_1 is 0, all middle indicators equal g_n, and at the final
index f_m=g_m-1. Both indicators are binary, so f_m=0 and g_m=1. The one
extra final F restores exactly that missing last 1.

Consequently the entire finite literal word, including its final forced
reset, is

    W_m=V B^k R,
    B=encode(g_2,...,g_q,g_1), 0 -> V, 1 -> VF.

R is the appropriate prefix of B. B ends in VF and has P=122q+2p bits and
77q+p odd bits. Its multiplier logarithm is

    log rho=log(4/3)(q theta-p)<0.

The independent program additionally checks the complete floor construction
for 1,385 reduced rational models with q<=30 and five admissible inverse
congruence lifts each. These are controls of the general residue proof.

## 6. Equal futures fit, and they force the claimed endpoint descent

Take r=ceil(122m/100). The literal periodic tail gives equal factors starting
at a=122 and a+P. Since r<=2m, P<124q, m>2q and q>=2,

    122+P+r < 122+124q+2m < 122m <= N.

For the middle strict inequality, 120m>240q and 116q>122 suffice. Thus both
future blocks lie wholly in the actual finite word. This conclusion does not
require the much stronger hostile condition r>=122+P.

The first checkpoint is a V output, of type 4. The second is the endpoint
of B after the initial V; B ends in VF, so it has type 5. If the actual
endpoint did not descend, the matched-block theorem would force these two
checkpoints to be equal. The type contradiction proves strict endpoint
descent for every positive lawful source.

Every finite parity cylinder has infinitely many positive integer lifts,
so the result is not vacuous. On the other hand, the existence of those lifts
is not substituted for a non-descending endpoint witness. The proof rules out
all d>=0 endpoint displacements for the specified family.

For an independent whole-word control, the first two convergent members have
m=19 and m=128, with N=2322 and N=15642. Three lawful source lifts of each
were replayed bit by bit. Their endpoint descent, selected types 4 and 5,
literal equal future factors, and exact 81*2^r difference lower bound all
check. No enormous orbit or full enormous word was materialized.

## 7. Continued fractions and unboundedness

The prime-factor contradiction establishes irrationality of theta. For an
upper simple continued-fraction convergent p/q and its predecessor p'/q',
the determinant is pq'-p'q=1. With next complete quotient t>1,

    theta=(pt+p')/(qt+q'),
    p/q-theta=1/[q(qt+q')].

Let a=floor(t). Irrationality makes a<t<a+1 strict, and set
m=(a+1)q+q'. The determinant gives mp=1 modulo q. Moreover q'>0 for the
upper convergents at issue and a>=1, so m>2q even when a=1. Thus

    m-q=aq+q' < qt+q' < m

is exactly the required two-sided approximation bracket after inversion.

The independently certified interval 1/10<theta<1/9 makes the first upper
convergent 1/9; all relevant upper convergents have 0<p<q and q>=2. The
standard denominator recurrence and alternating determinant prove infinitely
many reduced upper convergents with denominators tending to infinity.
Equivalently m=q_next+q, so these constructed m are unbounded. Since theta>0,
the explicit-reset count floor(m theta)+1 is unbounded as well.

No hypothesis that partial quotients are large, unbounded, equidistributed,
or typical is used. Even a hypothetical uniformly bounded partial-quotient
sequence would support this descent-family argument. The prior pure-expansion
repetition avoidance applies eventually because these m tend to infinity,
including the audited m>915 sufficient threshold under K^4F/V reparsing.

The independent logarithm enclosure recovers all 39 upper-convergent controls,
including the small and enormous displayed certificates. The infinite-family
claim relies on the determinant/denominator proof, not this finite list.

## 8. Exact enormous certificate and the stronger formal hostile

The independent logarithm enclosure uses 192 terms for log 2 and log(3/2),
rather than the submitted 512-term log 2/log 3 calculation. Writing

    theta=(77 log(3/2)-45 log 2)/(log 2-log(3/2))

and applying rational interval bounds directly recovers both strict brackets.
It independently verifies the enormous final floor

    1522163965844627930660358762978517428452,

and length

    N=1833237826965784108483587883633758044204176.

The exact large dimensions, 187 complete periods and final partial-period
length agree with the receipt. N exceeds 2^52, and r>122+P. The smaller
control has N=136401326 and is below the old short-range cutoff. Thus the
large arithmetic exclusion genuinely goes beyond that cutoff.

For any prefix of A B-infinity, equal r-factors with r>=|A|+|B| have equal
entire infinite continuations. After at most |A| aligned positions, both
streams are in the periodic tail; the original common prefix then contains
a complete |B|-period of aligned agreement. All later agreement follows by
periodicity. This is valid for arbitrary factor starts, including preperiod
starts, and does not assume that B is primitive.

If starts i<j have identical infinite continuations, their difference d=j-i
is a period of the stream from i onward. Its odd density equals the eventual
B-tail density h_B/P. Therefore the intervening multiplier has logarithm
(d/P) log rho<0. This proves the claimed contracting-separator property.
The independent controls check 930 short binary preperiod/period examples and
11,060 matched pairs, including nonprimitive periods, exact separator-density
equalities and a below-threshold mismatch countercontrol.

Importantly, the infinitude of the descent family does not prove infinitude
of this stronger formal hostile property. That property additionally needs
r>=122+P, imposing a substantially larger m/q ratio. Only two of the first
39 controlled convergent members satisfy it. The hostile note explicitly
does not claim an asymptotic abundance or infinitude of those stronger
hostiles. There is no hidden large-partial-quotient assumption.

## 9. Cycle-equation equivalence and its exclusion

For W_m=V B^k R, non-descent and the chosen long future match force y=T^122(x)
to return after B. Since B contracts, its necessary fixed point is

    y=b_B/(M_B-a_B).

The companion note's three-test equivalence is correctly scoped to the full
endpoint condition: positive integral y; positive integral V preimage x; and
the required final-phase equality T_R(y)=x+d.

Sufficiency does not hide an extra parity hypothesis. For a parity word with
odd affine coefficient, endpoint integrality is exactly its unique source
congruence modulo 2^n, which authenticates the entire chronological cylinder.
Thus an integral B fixed point authenticates B, an integral V preimage
authenticates V, and the prefix R is lawfully followed. Necessity uses the
matched-block theorem and the actual initial/final phases, as stated.

Those conditions are incompatible here: the same y would simultaneously be
a V output of type 4 and a B output whose terminal VF gives type 5. The
image argument supplies a concrete arithmetic contradiction to the cycle
side condition; it is not a claim that every contracting formal branch lacks
positive integral fixed points.

## 10. Audit boundaries and final disposition

No correction is required for the frozen theorem. Its published README identifies this completed audit; this is a presentation change rather than a mathematical repair.

The important retained boundaries are:

1. The matched-block height theorem and its cited explicit logarithmic-form
   input remain dependencies.
2. The particular upper-convergent family has actual endpoint descent for
   every positive lawful source, with unbounded episode count.
3. The much stronger formal hostile property is proved for the two exact
   controls and conditionally for any member meeting its extra length ratio;
   no unproved asymptotic partial-quotient claim is used.
4. Arbitrary bank entry, the general realized S_d boundary and Collatz
   convergence are not concluded.

Reproduction:

    python3 review/independent_image_guard_controls.py
    python3 -O review/independent_image_guard_controls.py

The outputs should equal independent-receipt.json and
independent-optimized.json. The submitted reruns are preserved separately.
