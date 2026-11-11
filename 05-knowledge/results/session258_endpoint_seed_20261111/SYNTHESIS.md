# From gcd one to endpoint-realizable seed cores

**Session258, 11 November 2026. New proofs and finite exact constructions; universal odd-slope seed existence and Chan–Pak remain OPEN.**

The main advance is a structural classification, not another common-divisor test: **a unit-offset tail seed has exactly one selected strong core, with only unselected transitive padding outside it.** For every fixed exact multiplier q, seed existence is now a finite, complete decision problem over bounded strong cores, with literal constructors when it succeeds.

## 1. The exact gap after C_k(T)=1

The inherited definition is C_k(T)=gcd{H(U):U extends T by k named vertices}, and C_(n−1)(T)=1 is inherited as proved. It says a common divisor eventually disappears. It does not supply a prescribed orientation, an exact count, or a unit-offset ray seed.

Three parameters must remain separate:

- **N:** actual Hamilton-path count. N=7 and21 are excluded.
- **c:** odd ray slope in N=1+c2^s. Every odd c is meaningful, including7 and21.
- **q=c2^j:** exact multiplier U of a selected-set seed. The seed equations are H(B)−W_B(O)=1 and U_B(O)=q.

The previous bank covers all odd slopes c≤4097 using a mixture of selected-set tail roots and doubling vertices. Those mechanisms must not be conflated. A seed gives the ray from s=j+1 onward; its earlier points are a separate obligation.

## 2. Endpoint constraints become a positive cut, not an integer ideal

For distinct old vertices a,b, let A_ab sum H(T−V(P)) over all directed simple a-to-b paths P having an even number of vertices. If the newborn beats exactly X, then

    H(T+v_X)=H(T)+2 sum_(a∈X,b∉X) A_ab.

The [endpoint proof](endpoints/PROOF.md) derives this by splitting the odd-cycle expansion at the unique new-vertex cycle. Every original arc has positive capacity, but nonarcs can have capacity through longer paths and must be retained.

This gives an exact constructor and its obstructions: a target half-gap r must be an actual directed-cut weight; every capacity>r becomes an endpoint implication; heavy-arc SCCs cannot be split. Zero-gap births are exactly singleton-SCC insertions.

A sharp consequence closes an entire naive recursion:

> A single vertex can increase H by2 only as1→3 or3→5.

The proof forces a unique capacity-one crossing, transitive shores, and at most one old directed triangle. It does not restrict edge flips, two births, root changes, or graph replacements. The H=5 hostile already has gcd-one extension image{5,9,11,13,15}, but no+2 extension.

For ray seeds, the same capacities give exact coupled budgets:

    sum_(a,b∈O) A_ab=(H(B)−1)/2,
    sum_(a∈O) [e_a+2 sum_b A_ab]=q.

These are a selected induced-edge budget and an endpoint-adjusted vertex budget, imposed on the **same** base and selection. Neither congruences nor separately attainable totals discharge them.

## 3. Strong-core rigidity removes the unbounded padding problem

The [core classification](seed_decision/CORE_CLASSIFICATION.md) proves:

    B=TT_a → S → TT_b,       O⊆V(S),
    H(S)−W_S(O)=1,
    U_B(O)=2^b U_S(O),

where S is strong, including the singleton case. Every exterior vertex is unselected. Conversely every such core and padding is a seed.

The mechanism is positive. If two selected vertices lie in different SCCs, their two ordered endpoint-cover couplings already exceed the entire H(B), contradicting W=H−1. With all selections in one SCC, every exterior strong factor divides the offset1 and therefore must equal1. Splitting each Hamilton path before each selected vertex also proves

    U_S(O)≥|O|H(S).

Every nontransitive seed has |O|≥2. Hence a core for exact multiplier q has H(S)≤q/2. The inherited Busch strong-order floor bounds its size logarithmically in q. Arbitrarily long padding has been eliminated exactly, not truncated experimentally.

The resulting finite strong-core classifier decides every exact q. Throughq29 it finds precisely

    q∈{1,2,4,6,8,11,12,16,20,22,24,29}.

All other q in1..29 have no selected-set tail seed at any order. Literal core witnesses and padded constructors are in the receipt. For example q6 supplies c3 from H=3,W=2,U=6; q11 supplies c11 from H=5,W=4,U=11. The small positive seeds are constructions, not new claims of numerical coverage beyond the inherited bank.

The q7 absence is a useful joint hostile:15 and29 are individually realized, yet they cannot be the first two tail counts of one unit-offset seed with U=7. Likewise no q9 tail seed coexists consistently with the inherited doubling-vertex ray of slope9.

## 4. Positive continuation and native Zeckendorf inverses

An independent [nine-state construction](seed_decision/PROOF.md) retains the first path's empty/end-selected state and the second path's empty/start-selected state. Ordinal block concatenation is a positive linear transfer. Capping coordinates at q+2 preserves the exact seed acceptance predicate and yields a finite all-length BFS, shortest block-word witnesses, and finite absence certificates. It independently confirms the all-order q1..13 classification.

The [historical bridge](grammar/HISTORY_BRIDGE.md) uses the repository's actual inverse65 and native Zeckendorf advances:

- Singular positive block maps are inverted relationally by27-coordinate suffix profiles, rather than pretending to have scalar or rational inverses.
- The native two-register Fibonacci reader retains the shifted value and previous-digit guard.
- The three physical U channels require an explicit carry cocycle with four possible carries; disjoint path classes do not mean disjoint Fibonacci digit supports.
- Forward states and backward target guards meet in a literal reconstructed graph. Finite graph distance supplies the “path of least resistance” rank.

The earlier append-only spectral obstruction remains valid. The present compiler changes the model: it reads the entire target, uses a target-dependent alphabet, retains endpoint channels, and may rebuild the source. The exact TM-to-CA interpretation preserves its halting predicate and grammar quantifier; it is not an argument that computational universality implies mathematical coverage.

The [two-positive-probe lemma](grammar/PROOF.md) gives a complementary interface: on the same rooted base, actual tail totals1+2q and1+4q are equivalent to the unit-offset equations. Two unrelated tournaments with those totals are insufficient.

## 5. Carmichael, Catalan, and99·2^k+23

The [arithmetic note](arithmetic/RESEARCH.md) recovers Session200's physical Carmichael twins and joint descent, the Catalan ordered-carry denominator, and the latest actual ready-return family. It proves that every99·2^k+23 lies outside the two-run Carmichael chamber2+(2^r−1)(2^s−1), despite having the inherited ready representatives. Multiplicity at a represented target and filling an empty target fibre remain distinct.

The25 primes below100 give an exact phase sieve on k, with valid period428468040. They do not prove primality. The first surviving composite is

    k16:6488087=1301·4987.

All terms are2 modulo3 and5 modulo9, so none is a square or cube. More decisively for the seed task, for k≥2 their canonical odd slopes are99·2^(k−1)+11 and their2-adic heights are1: varying k is not one fixed-c unit-offset ray. The inherited full response contains3^k even though H is dyadic; no polynomial in H recovers that full continuation coordinate on the whole family.

## 6. What remains to prove

The structural result reduces universal selected-set seed existence to a precise **strong-core odd-part surjectivity** problem:

    Every positive odd c is oddpart(U_S(O)) for some strong S
    with W_S(O)=H(S)−1.

The fixed-q compiler decides each candidate exponent, and trailing unselected sinks preserve c while increasing that exponent. What is missing is a uniform nonemptiness invariant or bound guaranteeing a successful exponent for every c. Then one still needs uniform construction of each line's below-start points to prove Chan–Pak. No such theorem is claimed here.

## 7. Incoming strong zero digit: a universal repair obstruction

The final incoming checkpoint5f24689c49a9 supplies repeatable strong H→3H and the repaired2673·2^k+645 family. Our+2 theorem upgrades its finite5991→5993 newborn-cut test: **no realizer of5991, at any order, has a one-vertex5993 child**, and the analogous statement holds for every k in that family. The incoming two-birth zero digit remains valid; the nonzero repair must change its interface. [Full integration and content distinction](INCOMING_SYNTHESIS.md).

## 8. Evidence and publication

The new endpoint identities have exact controls over1,099 labelled bases and33,866 births. The global nine-state compiler has an independent permutation/split audit of all146 transfer vectors and complete q1..13 closures. The stronger core controls classify q1..29. The Zeckendorf bridge has33,462 pullback,11,154 charge-interface, and8,000 cocycle checks. All tests retain their explicit checks under optimized Python.

See [reproduction](REPRODUCE.md), [source pins](SOURCE_PINS.json), and [open questions](OPEN_QUESTIONS.md). Session257 and the subsequent concurrent checkpoints, including5f24689c49a9, were inherited and preserved. Classical inputs are clearly separated from new proofs and exact experiments; no Lean formalization, external priority search, or universal conjecture closure is asserted.
