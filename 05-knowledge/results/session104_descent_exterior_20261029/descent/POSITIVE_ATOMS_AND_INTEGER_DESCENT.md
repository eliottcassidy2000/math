# Positive target atoms and literal natural-number descent

**PROVED (elementary, scoped):** the extraction and rank conversion below.
**INHERITED:** atomic coverage, rank/credit conversion, affine carries, and the (3,3) potential.
**NEW SPECIALIZATION/ASSEMBLY:** an exact integer rank for each active (3,3) return, a source-marked carry debt, an explicit aggregate-budget-to-bounded-rank interface, and a bounded-reset compiler.
**OPEN:** a source-complete Collatz rank or the necessary atomic tail/budget estimate. No priority claim.

## 1. What positive weight actually buys

Let X be countable, T deterministic and total outside a genuine terminal set H. Stop upon first entry into H. Let tau(x) be that time (infinity if none), and fix one probability law w with w(x)>0 for every target x. All clocks below count actual T steps, not controller actions unless explicitly declared.

**Atomic extraction.** Write A_K=sum_x w(x) 1[tau(x)>K]. If A_K<w(x), then tau(x)<=K. Consequently A_K -> 0 is equivalent to termination at every target. Proof: a surviving target alone contributes w(x); conversely, if all terminate, split the probability into a finite head and a tail of arbitrarily small mass, then wait until every head target has stopped. This is the already inherited atomic coverage mechanism, now stated for stopped native clocks. Fixed full-support weights prevent a single exception from disappearing; they do not prove A_K -> 0.

**Budget extraction.** Suppose for every integer K>=1,

    sum_x w(x) min(tau(x), K) <= M < infinity.                 (1)

Then tau(x)<=floor(M/w(x)). Indeed w(x)min(tau(x),K)<=M for every K; an infinite tau or a larger integer tau contradicts this. Conversely, a finite weighted expected termination time implies (1). Universal termination alone need not imply this stronger finite-mean property for a specified w: countdown targets of length 2^j with weights 2^(-j) (j>=1) all stop, but have infinite mean.

One sufficient way to obtain (1) is a NONNEGATIVE stopped potential V_t on the same source probability space with finite E V_0<=M and

    E V_(t+1) <= E V_t - A_t                                (2)

for every t. Summing gives sum_(t<K) A_t = E min(tau,K)<=M. This is not an arbitrary negative-average heuristic: the exact survival indicator is charged on every active native step. Conditional drift is a stronger familiar sufficient condition; the displayed aggregate condition itself suffices for deterministic targets. For a genuinely randomized transition process, almost-sure termination does not exclude zero-probability infinite realizations, and it does not automatically rank every permitted edge.

If w(x), M are effective rationals and T,H are decidable, the bound B(x)=floor(M/w(x)) allows bounded simulation to compute r(x)=tau(x). Then r(Tx)=r(x)-1 before H. On active states r is a positive integer; r=0 at terminals. This is a constructive compiler only **after independently proving** (1). Computing tau until termination without such a proof is not the missing proof.

More generally a computable certified tail modulus K(x) with A_K(x)<w(x) supplies a bounded rank search. The bare limit A_K->0 does not supply an effective modulus by itself. A primitive-recursive T,H plus a primitive-recursive certified time bound produces a primitive-recursive tau by bounded search. Requiring a primitive-recursive rank is stronger than merely requiring termination: for a deterministic machine computing a total non-primitive-recursive function, a primitive-recursive bound on its primitive-recursive step iteration would make its output primitive recursive, a contradiction. We make no assertion that Collatz requires a non-primitive-recursive rank.

Weights must stay attached to the original source. Resampling a fresh prior at each endpoint, conditioning on survival without its normalizing mass, or replacing atomic source weights by Haar cylinder weights invalidates the lower bound w(x). For a fixed event G containing x, the correctly conditioned atom is w(x)/w(G); it is not an unchanged prior. An expectation of a signed defect can cancel; (1) uses nonnegative costs.

## 2. Existence of a rank is not a recipe

For deterministic T with a specified terminal H, the following are equivalent: every source reaches H; a map r:X->N with r=0 on H and r(Tx)<r(x) off H exists. Forward choose first hitting time; reverse no infinitely decreasing nonnegative integer sequence exists. Using positive integers everywhere simply replaces r by r+1 and terminal zero by terminal one. Stop transitions at H: even the 3x+1 accelerated fixed point 1 cannot strictly decrease if one continues its self-loop.

This equivalence alone restates the problem. A useful certificate supplies r and proves every native transition inequality without assuming termination, or supplies a simpler well-founded dependency reduction with authenticated endpoints. There is no claim that every nondeterministic terminating relation has an N-valued rank: unbounded branching can give arbitrarily long finite paths from one root, requiring more general well-founded orders.

## 3. Literal positive-integer rank for the inherited auxiliary (3,3) controller

For positive odd n, use T_3(n)=(3n+1)/2^a with exact a=v_2(3n+1). Active queue Q>=1 changes by +2 if 3 divides a, otherwise -1; stop at Q=0. Inherited Session100 proves for

    P(n,Q)=(n+1)(3/2)^Q, rho=27/28

that every leaf weakly decreases P and every branch multiplies it by at most rho. The active lower bound is P>=3. Define, using exact rational comparisons,

    B(P)=min{k>=0 : rho^k P < 3},
    R(n,Q)=3 B(P(n,Q))+Q  (Q>=1); R(n,0)=0.              (3)

B exists by geometric decay. If P'<=P then B(P')<=B(P). If P'<=rho P, B(P')<=B(P)-1: the exponent B(P)-1 already gives a product below 3. Thus a leaf lowers R by at least 1, while a branch lowers 3B by at least 3 and adds only 2 to Q, again lowering R by at least 1. A transition directly into Q=0 also lowers the positive rank. Every auxiliary return therefore stops in at most its initial R actual accelerated steps.

For example the native active run `(13,1)->(5,3)->(1,2)->(1,1)->(1,0)` has valuations `(3,4,2,2)` and exact ranks `163->162->38->4->0`. The rank is conservative, but every decrease is literal.

No floating-point logarithm is needed. This is the October11 energy/credit-rank mechanism specialized to the Session100 auxiliary controller; it is not a new abstract ranking theorem. It also produces an explicit native-step bound, unlike a bound only on completed action receipts.

**The reset debt remains.** Restarting Q from 0 to 1 multiplies P by 3/2. It can increase (3); termination of each return does not imply termination of infinitely many returns or descent below an immutable original source. No reset is silently authorized as a decreasing edge.

For b=5, the actual 13->33->83->13 cycle has valuations (1,1,5). With the (5,5) rule, starting at 83 with any active Q follows (Q,Q+2,Q+1,Q), so the native state repeats. Strict rank inequalities around the three-edge loop would give R<R. An unbounded history sidecar cannot repair strict well-founded decrease along the resulting infinite run either.

## 4. A carry debt that says precisely when an original source is paid

For immutable odd source n0, multiplier b, and chronological exact valuations a_1,...,a_J, keep D=sum a_i and C initialized at zero with

    C' = b C + 2^D; D' = D+a; J'=J+1,
    2^D n_J = b^J n0 + C.

The exact source cylinder is n0 == (2^D-C)(b^J)^(-1) mod 2^(D+1). Keep it or replay every exact native valuation; mere divisibility is weaker.

Define signed integer payment

    L_J=(2^D-b^J)n0-C = 2^D(n0-n_J).                    (4)
    L_(J+1)=b L_J + 2^D[(2^a-b)n0-1].                 (5)

L_J>0 iff the displayed endpoint is smaller than the ORIGINAL source. It is not a positive rank, and is not monotone. Its role is an exact debt ledger preventing a slope-only argument from declaring payment. First-return ranks must be coupled to an independently proved eventual L_J>0 or to another genuine terminal.

If every odd n>1 has a finite native word ending below n, strong induction proves all reach1. Equivalently the accelerated record-low sequence is a strictly decreasing positive-integer rank at successful block boundaries. This standard first-descent criterion is useful only when every block has a separately proved finite completion and true source payment. Searching until a descent appears assumes exactly the missing all-source obligation.

For b=5 and w=(1,1,5), J=3,D=7,C=39. Its affine slope125/128<1 but its fixed point is39/(128-125)=13. At n0=13, L_3=0. For b=3 the same word has C=19 and fixed point19/101; its rational orbit is 19/101 ->79/101 ->169/101 ->19/101. Positivity in the real numbers does not make these integral sources. Session100 proves this denominator persists through inverse prehistories; Session101 shows appending and reclosing is a different operation, with denominators35/175/7 depending on carry and attachment. A denominator is not a descent rank on these rational cycles, and reclosure is not an actual native step.

## 4a. Compile paid resets into one natural-number rank

Here is a concrete way to join locally decreasing episodes without appealing only to a lexicographic order. Label each macro source by a positive integer n, held immutable throughout its episode. Suppose independently supplied data give:

1. A finite nonnegative integer bound B(n), and a local integer rank 0<=r<=B(n).
2. Every internal native edge keeps n fixed and strictly lowers r.
3. At local completion r=0, any nonterminal reset is authenticated to a new source n'<n and a refreshed local rank 0<=r'<=B(n').
4. Every nonterminal state has a permitted internal edge or such a paid reset; terminal labels are genuine target certificates, not unexplained dead ends.

Set

    F(n)=sum_(1<=m<n) (B(m)+1),    G(n,r)=F(n)+r.           (6)

Every internal edge lowers G. At a reset,

    G(n',r') <= F(n')+B(n') < F(n) = G(n,0),

because the sum F(n)-F(n') includes B(n')+1 and all other terms are nonnegative. Thus G is one ordinary nonnegative-integer rank across BOTH native steps and resets, not just an ordinal-valued lexicographic rank. Adding one makes every value positive if desired. A terminal state need not have G=0; it simply has no outgoing proof obligation. The displayed assumptions exclude infinite runs. If B is computable (respectively primitive recursive), then so are F and G; no closed-form evaluation of the finite sum is needed. Finite availability and genuine-terminal condition4 are essential: a decreasing proof procedure with an unresolved dead end is not a target proof.

For positive odd Collatz sources use the ordered labels j=(n+1)/2, with B(j)=R(2j-1,1) from (3). An episode's exact local controller rank is bounded by its initial B(j), since internal edges decrease. **If** every completed return from a nonterminal source n>1 had endpoint n'<n, these local ranks would compile by (6). But this missing premise is false for the raw controller: starting (n,Q)=(3,1), the exact valuation is1, so the first return is (3,1)->(5,0). Its endpoint is larger than its immutable source. The original-source payment in (4) is L_1=2(3-5)=-4, not positive. Resetting at5 cannot be declared paid. The repaired compiler therefore applies only to an independently authenticated family of source-decreasing macro episodes (possibly grouping several returns with their own termination account), not automatically to every (3,3) return.

This supplies a precise engineering target for the analogy: establish a bounded local completion rank, retain immutable source and carry through the whole macro, and prove positive source payment before every reset. The arithmetic source-payment requirement is not removed by the finite sum F.

## 5. Rounding and weighting obstructions

A real potential that strictly decreases need not yield an integer rank: 1+1/(t+1) decreases forever. If instead V>=0 and V(Tx)<=V(x)-epsilon for a fixed epsilon>0 on every active edge, floor(V/epsilon) decreases by at least1 (use terminal stopping); ceil also works. Variable steps tending to zero do not justify that rounding. Rationality alone does not help, because rational denominators can grow.

The multiplicative construction (3) works because it combines a uniform ratio, a positive active lower bound, and a credit that pays neutral edges. An averaged contraction does not give that pointwise statement. Likewise polynomially decaying full-support source priors do not themselves become contracting flow weights: the inherited actual word of d ones sends 2^(d+1)-1 to 2*3^d-1, obstructing bounded positive corrections of those priors for universal discounted incoming flow. Source measurement and transition payment are distinct accounts.

## 6. Source table and evidence scope

| Source actually inspected | Input used | Status here |
|---|---|---|
| `math/05-knowledge/results/collatz_effective_prefix_mass_20261005.md` | positive atomic coverage; missing mass estimate; discounted flow caveat | inherited |
| `math/05-knowledge/results/collatz_atomic_rank_interface_20261011.md`, A1–A4 | conditional drift; atomic integrability; exact energy/credit integer rank | inherited template, not rediscovery |
| `math/05-knowledge/results/collatz_atomic_run_bridge_20261011b.md` | native common-future receipts, action versus physical clocks, terminal-child obligation | inherited boundary |
| `math/05-knowledge/results/collatz_composed_weight_floor_20261011.md`, C3–C4 | target floor versus alias tail; fixed source law versus flow; expanding-word hostile | inherited obstruction |
| `session100_reframes/controllers/RETURN_CONTROLLER_TRAP.md` | rho=27/28, active bound3, 5-cycle trap, reset cost | inherited proved controller; (3) specializes it |
| `session100_reframes/drift/DRIFT_TARGET.md` | source-marked clock barrier; Haar-null target problem | inherited boundary |
| `session100_reframes/bijections/VALUATION_CONJUGACY_AND_INTEGER_SIDECAR.md` | 19/101 hostile; inverse denominator persistence | inherited proof inspected; not rederived here |
| `session99_collatz_poles/transfer/CARRY_CLOCK_POLES.md` | native guards; exact carry; raw versus normalized moments | inherited |
| `session101_denominator_connections/continuations/CONTINUATION_STATE_AND_SYMMETRY.md` | affine concatenation and reclosure attachment | inherited |
| Cook, Kroening, Rümmer, Wintersteiger, *Ranking Function Synthesis for Bit-Vector Relations*, Definition1/Lemma2, [Microsoft Research PDF](https://www.microsoft.com/en-us/research/wp-content/uploads/2010/01/ckrw10.pdf) | classical well-founded ranking implication | external primary-source confirmation, not a novelty claim |

All local sources were read as filesystem files, not fetched through lazy git. The scripts and proofs below are not a claim that historical receipts have been rerun. `check_descent.py` gives independently coded finite exact controls for the new specialization and interfaces, including malformed-input rejections. The infinite statements are proved above, not inferred from finite samples.
