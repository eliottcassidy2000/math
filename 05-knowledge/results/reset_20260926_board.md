# Precision resets: sharp prices, sufficient precision, and faithful small codes

**PROVED scoped results / FINITE-EXACT certificates / independently audited
as specified in the audit note / Collatz OPEN, 2026-09-26.** No novelty claim.
This continues the [carry-defect forest](nextforest_20260926_board.md).

## 1. The inequality now has an exact object and clock

Let U(n)=oddpart(3n+1). Retain the full state y=q2^r+u, with q,u positive
odd and r>=1, and the split product K=q*u*2^r. Let m be the ACTUAL number
of binary divisions in the next U step. Away from a collision,

    K_next/K=(9+3/u)/4^m.

This is unchanged by exchanging coefficient/core roles. In particular
m>=2 contracts K by at least the factor3/4. The newly large precision
at a swap is accounted for by the other register, rather than treated as
free resource. Valuation-one steps still expand it.

At a collision, call the split balanced if y/4<=u<=3y/4. The
[product note](reset_20260926_inequality.md) and
[clock note](reset_20260926_clock.md) prove:

| Exact operation | Uniform product bound | Scope |
|---|---|---|
| Collision, then split its immediate endpoint | K_new/K <3/4 | sharp over balanced admissible triples |
| Collision, one EXTRA U step, then split (inherited P5) | K_new/K <27/16 | sharp; can exceed1 on an actual canonical trajectory |
| Same P5 operation, extra step divides by at least4 | K_new/K <27/64 | sharp over balanced admissible triples |

All constants have explicit equality-limit families and integer proofs.
For the middle row, an actual canonical trajectory is
177->133->25->19->29. The split at25 has K=126; the fresh split at29
has K=208. A clock that forgets the extra arithmetic step would falsely
claim contraction here.

Outside balance, the immediate-reset estimate retains the explicit debt

    K_new/K <= [(3y+5)/(4(y+3))] * [3(y-1)(y+3)/(16K)].

That last factor can be unbounded. Also, K alone does not measure integer
descent: writing B=4K/y^2, actual descent is equivalent to
K_end/K_start < B_end/B_start. Endpoint balance must remain in the account.

## 2. A complete collision phase is still not a descending return

There is an exact hostile starting from canonical decompositions:

    n_H=2^H-1, H even,
    first collision endpoint=(3^H-1)/2^(2+v2(H)).

The phase consists of H-2 consumes, a swap, and a collision. For H=4k+2
the endpoint/source ratio tends to infinity. Examples include63->91
and1023->7381, where the arrow denotes the WHOLE first-collision phase.
The final collision is almost perfectly balanced and locally contracts;
its refund does not pay the accumulated growth of the phase.

The ternary consume-run clock reaches H-2 before the swap. Its depth is
linear in H, while the extra collision cancellation is only2+v2(H).
Thus the ancestral debt cannot be erased at a reset merely because the
new representation is canonical. This refutes a specific return selector,
not Collatz, nor every variable-depth selection of several phases.

## 3. Sufficiently large recreated precision can instead certify descent

The productive complement is the
[large-swap lifting certificate](reset_20260926_swaplift.md). At a swap
with r>=2, the actual integer already decreases. The difficult case is
r=1. Write R=a-1 for the new precision and v=U(u); then

    y=(2^(R+1)*v+6q-1)/3 -> v*2^R+3q.

This exposes the exact smaller computational problem: replay a finite
prefix of the core3q until its value b is below2q. If that prefix has J
steps and total division budget A, then R>=A is sufficient for a smaller
actual iterate within J+1 steps, uniformly over EVERY admissible positive
odd coefficient v. The endpoint equality R=A is a final collision whose
extra divisions only help. The companion note sharpens the threshold by
also treating a partially consumed final division.

The optimized finite bank proves a concrete sufficient reset inequality:

    q odd, 1<=q<=341, R>=64
       => a smaller actual iterate within41 U steps.

The worst core in this bank is3q=27, at q=9. Per-coefficient thresholds
are retained in the certificate and are usually much smaller.

This is a witness-lifting theorem: a checked finite core prefix is expanded
into an all-height family of actual source-preserving descent certificates.
It covers arithmetic-progression tails of positive density, not merely
isolated values close to powers of two. Its precise finite bank, thresholds,
and hostile boundary are recorded with a reproducible certificate.

At the opposite boundary, fixed q=1,R=1 admits the family

    n_k=2^(3k+2)-5=27,251,2043,16379,..., k>=1,
    n_(k+1)=8n_k+35,
    U^(2j)(n_k)=4*9^j*8^(k-j)-5, 0<=j<=k.

Every one of its first2k iterates exceeds its starting value. Thus even
a fixed small coefficient and a fixed small precision permit arbitrarily
long first excursions; the remaining coefficient v retains essential
information for any UNIFORM bounded-horizon statement. Each source is
different, so this does not construct a divergent positive orbit.
Unlike the earlier constant-defect family27,91,347,..., this family has
a PROVED unbounded lower bound on its first-descent time. For X>=27 its
count up to X is floor(log_8((X+5)/4)), hence it has density zero. Its
three-bit recursion reflects the exact valuation word(1,2), consuming
three binary divisions for two multiplications by3 and multiplying n+5
by9/8. This is a concrete instance of the rational-shadow mechanism in
[THM-4507 / finite-valuation obstruction](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md),
with negative cycle -5->-7->-5; no new general shadowing theorem is claimed.

The universal coverage question remains: arbitrary q may require an
unproved core crossing, and a given swap may have too little precision.
No assumption that every remaining core converges is used in the certified
bank. The hard residual is the coupling between that core's required
division budget and the precision supplied by the actual orbit.

## 4. What the three colours really encode

[The marked Zeckendorf construction](reset_20260926_colours.md) makes the
two extra units precise. Use ordered weights1,1,1,2,3,5,..., with the first
three units black, blue, red, and forbid consecutive selected indices.
Every positive integer has exactly2 or3 representations. Exactly2 preserve
the inherited full charge: its ordinary Zeckendorf word and one uniquely
determined marked alternative. A single sheet bit distinguishes those two.

The supplied35-colour sequence is a legal isolated row under this grammar,
including its final blue. But all legal next rows of length36 force red at
position35. Thus the extra units resolve the atom-splitting problem; they
do not by themselves define a prefix-coherent continuation of the supplied
sequence. Another context or ordering rule remains possible, but must be
specified rather than inferred by changing a supplied symbol.

There is also an intrinsic ternary state in the arithmetic controller.
After every noncollision, exactly one of q,u is divisible by3. Which one
is divisible records the last branch. The nonzero3-adic valuation records
a consume-run length, which can be arbitrarily large. A finite role marker
and an unbounded depth register are different resources.

At fixed y the full split fibre has(y-1)/2 states; for odd y coprime to3,
floor(y/3) remain after the ternary role restriction. Hence a couple of bits
cannot preserve every split over the same y. This does not obstruct computing
U from y alone, choosing one canonical split, or finding a successful
quotient specifically for descent.

## 5. The PA/ZFC analogy has a concrete, useful reading

[The logic note](reset_20260926_logic.md) gives an exact pure-membership
presentation of finite set theory equivalent to PA: Extensionality,
Adjunction, and one structural-adjunction induction schema. The empty set
is derivable. Thus the proposed small syntactic package is realizable;
the schema carries unbounded mathematical strength.

For ZFC the note gives an explicit reflection presentation, and shows how
it too can be packaged into two fixed sentences and one schema. Choice
is retained explicitly. The substantive distinction is reflection of
formulas and their existential witnesses, not the sentence count.

One tempting distinction is false: finite-set arithmetic already supports
transitive containers containing every subset of every member. In Ackermann
coding, the container with member codes0,...,n does exactly this. That is
different from requiring the POWER-SET OBJECT itself to belong to the
container. An unbounded reflection schema supplies still more information.

The mathematical bridge to Collatz is precise. Encoding finite prefixes is
computable. Checking a supplied descent certificate is decidable. The open
part is proving that a certificate will be produced for every source. A
smaller representation must preserve actual witnesses, and its constructor
must terminate. The swap-lift theorem above is an explicit example of such
arithmetic witness transport on a proved domain; no set-theoretic axiom is
being imported as a Collatz descent estimate.

## 6. Inheritance, live board, and incoming work

Closest proved mechanism: exact swap/collision transport and the bounded
carry-defect first-descent bank. Canonical hostile:27 returns as the core
inside91. Corrected near miss: core convergence without enough precision
does not preserve the source shadow. Least-used sidecars: the two labelled
summands, final-step collision, ternary run length, and existential witness.

Anchor: precision-reset control. Niche: the three-unit representation lift.
Wildcard: finite-set coding and reflection. The board after the calculations:

| Live concept | Preserved object | Missing coordinate / next decisive test |
|---|---|---|
| Split product | exact noncollision multiplier | accumulated valuation-one growth and endpoint balance |
| Reset clock | number of actual U steps | distinguish immediate and P5 policies |
| Ternary role | which register inherited3q | retain unbounded run length and core value |
| Core lifting | exact finite descent witness | characterize insufficient-precision residuals |
| Marked Zeckendorf words | value, full charge, chosen sheet | state a rule extending35-blue to the next row |
| Finite-set reflection | formulas and transported witnesses | construct descent witnesses for every arithmetic source |

Incoming work at8287bf6b6 repaired and strengthened
[THM-4511 / Gilbreath size-four wall](../../01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md).
The useful comparison is an invariant safe boundary: in its stated
0/2/4 universe a2-wall contains arbitrarily many size4 defects. Our bounded
coefficient/large-precision domain is not proved invariant under Collatz;
q can grow and the next precision can fail its threshold. Therefore the
Gilbreath barrier is a prompt to prove residual-domain closure, not a
transported Collatz theorem. The incoming tournament symmetry corrections
do not alter any arithmetic certificate here.

## 7. Strongest surviving targets

1. Extend the verified core bank by structural families whose crossing
   below2q has an explicit division budget, rather than only more samples.
2. Characterize or charge the insufficient-precision residual using the
   full q/u history. An inequality must survive both the all-ones phase
   and unbalanced reset families.
3. Seek a quotient that preserves DESCENT WITNESSES, without requiring it
   to preserve every decomposition. The colour fibre obstruction leaves
   this weaker and potentially useful objective open.

Proofs, exact universes, controls and commands are in the five linked
notes. The [audit](reset_20260926_audit.md) records independent checks and
the one scope wording repaired before publication. No global Lyapunov
function, universal core crossing, or Collatz proof is claimed.
