# Charge receipts: polynomial carry loss and a source-preserving reroute

2026-10-11. **PROVED, scoped and self-audited:** the inequalities, precision
invoice, charged grammar and two all-height cylinders below. **FINITE-EXACT:**
the declared census and the complete smaller-source obstruction at39.
**OPEN:** independent budgets for every source and universal Collatz. No
independent peer audit or formalization is claimed.

## Inheritance and concept board

The nearest proved mechanism is [C1--C5, fixed-source shadow charge](collatz_shadow_charge_20261011.md):
a bounded staircase gives a finite deadline, and a genuine ROOT shadow at
sufficient binary precision transfers that bound to the *requested* integer.
The exact supremum `E_d(n)=2^-o_d(n)` is already available. We do not replace
that source by its shadow or infer boundedness from the existence of shadows.
The stronger prior [terminal basis](collatz_terminal_basis_20261007.md)
already grounds every odd source through32767 by retained actual edges. The
finite compiler census here is a comparison of restricted grammars, not new
finite Collatz coverage.

The canonical hostile is the reset-two inverse shield in
[joins and lifts, section3](collatz_join_shields_and_lifts_20261004.md).
The corrected near miss is equating a greater inverse depth with a new
source exit. The least-used sidecar is the *signed difference of odd charges*
on two legal routes with a common future. The banked source, the actual
child, and the guard all remain in the receipt.

Board: **marked source / distinct states / charge debt / first exposed join /
quintic port attachment / complete boundary relation**. Anchor: independent
charge budgets on recognizable source families. Niche: the
[eleventh-root port obstruction](eleven_port_layers_20261011.md). Wildcard:
compare a rigid port layer with an inverse controller whose target cannot
move. The comparison is about preserved interfaces, not a dynamical
conjugacy. Incoming main was fetched at startup and at research checkpoints;
no newer commits than b4eb117d0b were present at those checks. That checkpoint
already integrates the independently audited Stechkin and multiport work.

Use `U(n)=oddpart(3n+1)` on positive odd integers, stopping at the first1.
A ROOT word has letters `a_i=v2(3n_i+1)` and odd charge `o=sum(a_i mod2)`.
The inherited score is

    Z(n)=(3n+1)f(n)/4
        =2^-o product_(i=1)^(tau-1) 3n_i/(3n_i+1),

with `Z(1)=f(1)=1`. Scores on hypothetical unrooted positives are zero.

## R1. Distinct states improve exponential carry loss to polynomial loss

Let `s_j` enumerate the positive odd integers coprime to3 starting at5:

    s_(2r-1)=6r-1,   s_(2r)=6r+1,
    C_b=product_(j=1)^b 3(s_j-1)/(3s_j+1),   C_0=1.

Every rooted positive odd n satisfies

    Z(n) > 2^-o C_(o+1) > 2^(1-o)/(o+3).                 (1)

No sharpness of these constants is asserted. This strengthens the previous
`Z>(3/4)(3/8)^o`; odd charge still has to be supplied independently.

**Proof.** A strict ROOT trajectory has no repeated states: a repetition
would produce a cycle before its first1. All internal states are odd,
coprime to3, and at least5. Partition them into maximal blocks connected
by even valuation letters. There are at most `o+1` blocks.
For `z=n-1` put `H(z)=3z/(3z+4)`. An even letter contracts z by at least3/4,
and its state's carry factor is `H(z)/H(4z/3)`. As in C1, the product in
one block is strictly greater than `H(z_last)`.

The block endpoints are distinct. If their number is b, ordering them
increasingly puts the jth endpoint at least `s_j`. Thus their product is
at least `C_b>=C_(o+1)`. The empty internal list has product1 and also
satisfies the stated inequality. Finally `s_j>=2j+3`, so

    3(s_j-1)/(3s_j+1) >= (3j+3)/(3j+5) > (j+1)/(j+2).

Multiplication gives `C_b>2/(b+2)` for b>0. This proves(1), including n=1.
The proof never assumes distinctness for an unproved infinite orbit. QED.

An independent budget `o<=B` therefore pays the original source's floors

    f(n) > 4*2^-B C_(B+1)/(3n+1)
         > 2^(3-B)/[(3n+1)(B+3)],
    W(n) >= (9/16) f(n)^2.                               (2)

The existing native refinement in C5 now uses either improved floor in(2).
Its alias bill tends to zero while the source stays fixed. This transports
positivity already paid by B; it does not supply B.

## R2. The finite shadow receipt needs far fewer source bits

Keep the previous deadline: for n>1, let t be the least nonnegative integer
with `(n-1)2^B 3^t < 2*4^t`, and put `T=B+t`. By time T the actual source
has reached1 or spent its `(B+1)`st odd letter. For n=1 put T=0.

For n>1 it is enough to request binary precision

    D=1+floor(log2((3n+1)*(16/5)^(T-1))).                 (3)

At n=1 take D=1. The floor in the checker uses integer comparison of a
rational number; it is not a floating logarithm.

**Proof.** For an actual strict prefix of length j ending at y>=1 and with
total valuation A,

    2^A = [(3n+1)/y] product_(i=1)^(j-1) (3+1/n_i)
         <= (3n+1)(16/5)^(j-1).

Every intermediate nonroot state after the initial source is at least5.
At the stopping event j<=T, hence A<D. The full event, including each exact
valuation and its final oddness bit, is therefore determined by n modulo
`2^D`.

Supply a strict ROOT certificate for any `m>=n` with `m=n mod2^D` and total
odd charge at most B. The same event word occurs for m. A `(B+1)`st odd
letter would contradict its certificate, so the requested n must reach1
and have charge at most B. The shadow may use a changing terminal controller;
its source congruence and total charge remain mandatory. QED.

| n | supplied B | T | previous sufficient precision | new D |
|---:|---:|---:|---:|---:|
| 3 | 1 | 4 | 21 | 9 |
| 7 | 3 | 15 | 121 | 28 |
| 27 | 28 | 105 | 3676 | 181 |
| 223 | 18 | 78 | 2185 | 139 |
| 4591 | 41 | 167 | 9186 | 293 |

These table budgets come from separately checked ROOT words, so the table
validates the invoice rather than proving independent all-source bounds.
The theorem is conditional on *any* independently supplied budget/receipt.

## R3. A signed charge grammar has independent source-only bounds

If two actual strict prefixes u from n and v from h meet, transporting a
ROOT certificate from h gives

    o(n)=o(h)+delta,      delta=o(u)-o(v).                (4)

The sign matters. [THM-4555, trailing-ones collision switches](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md)
already proves the switches; the new consumer is their exact charge cost.
In its notation the cost is `D+o(u)-o(u')`. The reset switch costs+1,
whereas `(2,6,c)~(4,1,1,c+2)` and
`(2,2,10,c)~(6,3,2,1,c+2)` both cost-1 when deleting one trailing one.
ROOT is stopped at its first occurrence, with no padding used as evidence.

Put `E(n)=floor(log2(n+1))`. A finite grammar with a checked seed bank and
strictly smaller children h has a source budget whenever every rule pays

    delta <= r*(E(n)-E(h)).                              (5)

If `C=max_seed[o(seed)-r E(seed)]`, strong induction gives

    o(n) <= r E(n)+C                                     (6)

on the grammar's accepted set. This is a source-only formula; its verifier
recursively checks guards and smaller children, not the target's completed
forward orbit. A failed guard or uncovered leaf rejects this grammar only.
Termination of the recognizer follows from h<n, independently of Collatz.

The inherited ordinary one-step descents, proper sibling reduction and
uncapped reset all pay(5) at r=1. For a one-step odd valuation>=3, the child
obeys `h+1<=(n+1)/2`; an even valuation costs0. A proper sibling costs0.
The reset halves n+1 exactly. The two sporadic switches cost-1.

The [inverse-word compiler](collatz_uncovered_join_routes_20261007.md)
has child `h=(C_w n+D_w)/3^L`, with

    Q=2^sum(w), C_w=Q*4^k<3^L,
    D_w=Q*(4^k-1)/3-B_w,    U_w(h)=S^k(n), S(n)=4n+1.

Its exact integrality, oddness, positivity and h<n guards remain in place.
The final sibling join changes its valuation by an even number, so the
charge cost is exactly `-o(w)`. Every such rule pays(5) at r=1.

The experiment enumerates *all* positive valuation words of length<=L
with `2^sum(w)4^k<3^length`, and all allowable k>=0. Thus no arbitrary
letter cap is hidden in the following inverse-template counts. These are
complete inside this contraction compiler, not all possible controllers.
The checked bank is exactly the15 odd seeds1..29, giving r=1,C=24.

| inverse length limit | templates | accepted among 50,001 odd n<=100001 |
|---:|---:|---:|
| 2 | 4 | 3648 |
| 4 | 24 | 3652 |
| 6 | 137 | 3923 |
| 8 | 1045 | 4074 |
| 10 | 6577 | 4074 |
| 12 | 77441 | 4074 |

All eligible smaller-child rules are considered, rather than stopping at a
first unproved child. The census stores a proof receipt hash; final-policy
receipts and seed words are saved. A separate forward replay validates
charges after recognition. The budget `o(n)<=E(n)+24` is **not** universal:
the least counterexample to this inequality is327, whose actual charge is35
and whose E is8. Thus a policy enforcing this potential cannot cover327,
regardless of how its free inverse bank is enlarged.

## R4. Source39 has an all-depth inverse obstruction and a paid new exit

The complete ROOT trajectories of the19 positive odd sources below39 have
55 distinct odd states, with maximum3077. None is any of

    S^k(39)=39,157,629,2517;   the next value is10069>3077.

Therefore **no** smaller positive odd h, no valuation-word length, and no
k>=0 satisfies `U_w(h)=S^k(39)`. This excludes all smaller inverse-sibling
controllers at this source, including those outside the finite contraction
bank. It uses a complete finite union with an explicit height cutoff; it is
not an inference from a timed-out inverse search. The union is stored in the
reproduction data. It also shows that the first source39 endpoint occurring
on any smaller-source trajectory is19, at odd step5.

The reroute is the entire native cylinder, with the requested source kept:

    n=39+1024t --(1,1,2,1,4)--> h=19+486t,   t>=0.

The carrier is `(243n+251)/512`, the charge cost is3, and
`2(h+1)<=n+1`. Consequently this rule pays(5) at r=3. It transports any
supplied child certificate; it does not assert all these children are rooted.
The native word guard proves the statement at every t, not just at t=0.

A second controller is

    n=43+1024t --(1,2,2,4)--> h=7+162t,

with carrier `(81n+101)/512`, charge1, and `4(h+1)<=n+1`.
It already pays r=1. Adding it to the length12 grammar increases accepted
sources from4074 to4356. Adding the39 controller gives4769, with r=3,C=16
for the same seed bank. These are strict relative gains in the declared
compiler, not claims beyond the inherited finite ROOT atlas.

This is the reusable move: the failed inverse interface specifies the next
forward event to expose. Moving the controller's terminal is useful only
when the resulting child is actually smaller and its charge invoice is paid.
The next pending source in the final policy is75. Universal coverage would
need a proved entry mechanism for all such residuals, not a longer version
of the inverse search just excluded at39.

## Validation and transfer ledger

Reproduction:

    python3 04-computation/experiments/charge_ports_eleven_20261011.py
    python3 -O 04-computation/experiments/charge_ports_eleven_20261011.py

[Saved output](charge_ports_eleven_20261011.json) includes the finite universes,
seed words, charged proof receipts, the complete shield39 union, exact
precision invoices and the companion geometric colour lift. The checker
uses always-active exceptions, integer arithmetic and rational products.
The infinite statements above have explicit proofs; the census is not one.

| Transfer | Preserved predicate | Lost information / retained sidecar |
|---|---|---|
| ROOT word -> odd charge | Exact event count and score floor | Order and source are kept in the certificate |
| Even blocks -> distinct endpoint product | Lower carry bound | Legal order lost; used only as a lower bound |
| Shadow -> requested n | Stop event at precision D | Equality of integers is not claimed; source congruence and charge retained |
| Common future -> signed delta | Exact ROOT charge difference | Legality, actual child and strict smaller rank retained |
| Finite inverse failure -> reroute | All-depth exclusion at39 | No universal failure; changed terminal19 supplies the positive exit |

META-PATTERNS applied: search the statement before the method; correct the
object before sharpening the technique; retain selected guards; re-evaluate
a certificate after a fibre-changing operation. No new card is promoted.

Final replay completed with **6,419,872 always-active exact checks**. Normal and
optimized outputs are byte-identical. Source SHA256:

    676cac9e363d681c5fef3e4f82d2030124fc161b407ef375106dd83475ea228e

Saved-output SHA256:

    ce32d86c7f64bc4c6b291225f802072066d127a7943514c224fa44284c1c855a

The repository documentation checker, whitespace check and tracked-tree link
check pass. These are executable and self-audit checks, not an independent
mathematical peer review.
