# Exact bridges, retained interfaces and the remaining target theorem

## 1. What this session establishes

Let H(T) be the number of directed Hamiltonian paths of a tournament. The intended conclusion remains that every positive odd integer except7 and21 is H(T) for some tournament. **This session does not establish that conclusion or an eventual threshold.** Its advances are exact structural descriptions and constructive selectors inside rigorously specified families.

The organizing idea is **observable-preserving reduction with an explicit defect account**: export a tournament family to a simpler mathematical object only after recording which interface coordinates, baseline and admissibility constraints survive. This is a proposed synthesis principle, not a claim that an unnamed universal theorem has been proved. Older repo results already supply the fusion laws, root packets and several statistical identities used here; each lane credits its sources.

## 2. A complete fiber, not merely a convenient sublanguage

For the inherited five-vertex deletion core with H=13 and Q=104, all32 root cuts form a finite physical alphabet. Equality of the full central interface with a constant baseline is recognized by a bounded rational carry graph. Exhausting its mathematically bounded states is a proof of the full fiber, not inference from a finite word census.

Exactly four constant baselines have nontrivial fibers; the other28 are rigid. One nontrivial fiber is the regular language `(A | B E* C)*`. Pinned binary spins describe its words, with normalized energy

    j = sum(s_i) - 9 sum(s_i s_(i+1)).

For length k>=2, actual counts differ from the all-A count by `2*13^(k-2)*j`. The attainable normalized values are exactly

    {0} union {m-8r : m>=1, r>=0, 2m+r<=k}.

Translation reduces the range to the numerical semigroup <8,17>. From k=17 its56 boundary holes persist, while an explicitly described interior is gapless. The proof gives a target inverse, multiplicities, moments and actual rooted-class distinctions. The k=2 minimum boundary is handled separately.

**Power gained:** a complete infinite family has an exact grammar, support and inverse. **Information still missing:** actual mesh2*13^(k-2) and the absolute baseline are not eliminated. A long normalized interval does not become a long interval of odd path counts.

For a general finite carry language, connected integer edge flows encode exact length and energy; Euler trails reconstruct physical words. Connectivity is essential: disconnected balanced cycles give false selectors. This is exact feasibility, not a polynomial-time complexity claim.

## 3. A two-sided semigroup theorem isolates the physical supply problem

Suppose actual full-central banks have counts `b_i+2E_i`, with pairwise coprime deletion counts h_i>1. Write D=product(h_i), w_i=D/h_i. Fusion gives a weighted bounded sum with its actual baseline retained.

If E_i=[0,m_i] and m_i>=h_i-1, then

    sum_i w_i[0,m_i] = {z in [0,S] : z in Gamma and S-z in Gamma},
    Gamma=<w_1,...,w_r>, S=sum_i w_i m_i.

CRT residues give both an exact membership test and a constructive digit inverse. Thus the desired interval is controlled by two reflected finite lists of holes, rather than an uncontrolled search over compositions. The theorem is conditional only in its application to tournament banks: its arithmetic statement is proved, but this session does not construct all those banks.

The obstruction is also exact. Projection modulo h_i detects a missing digit independently of other banks. The named mixed13/15 alphabet never saturates the needed residues at any length. An all-order content argument, with a finite short-word base, explains this failure. A generalized finite alphabet with K_i=0 and integral Q_i/h_i also has an explicit deletion-product threshold beyond which content prevents saturation. These restrictions do not cover all tournaments or arbitrary reset operations. Unit-deletion padding supplies a positive finite example: six actual14-vertex full-central values6689,6691,6695,6697,6701,6703 share (h,Q,U,V,d,e)=(15,30720,24834,8259,0,15). Both-sided unit contexts retain their offsets. Their normalized digit set is {0,1,3,4,6,7}, not a complete interval. An attempted repeated slide pump fails because its endpoint defects cancel the later scalar gain; no unbounded complete-bank supply is proved.

**Power gained:** the missing constructive theorem can be stated as physical complete-digit supply plus baseline control. **Limitation:** arithmetic completion cannot be substituted for physical realizability.

## 4. Resetting the mesh creates a new interaction geometry

Reroot the preceding spin family at a specified vertex of its final A block, then prepend a rooted triangle. For k>=3 this gives an actual5k+3 vertex tournament family with2^(k-2) spin choices and count

    F_k(s)=c_k + sum_i a_i s_i + sum_(i<j) b_(j-i) s_i s_j.

All pair coefficients are strictly negative and are given exactly in the reset proof. Quadraticity is structural: only arcs incident with one old root vary, their indicators depend affinely on spins, and any Hamiltonian path uses at most two such arcs. No higher-order interaction is being discarded.

An18-state exact transfer produces coefficients involving13,104 and471; exact matrix annihilation and finite response identities propagate them to every length. The unique minimizing word has ones in its first floor(k/2) free positions and zeros afterwards; the maximizing word is its complement. A minimum-cut representation supplies an independent optimization lens.

A further all-order theorem recovers every free spin from the single reset count by outside-in13-adic decoding. Thus the reset is injective on all2^(k-2) states: an O(k)-value central collapse becomes exponentially many distinguishable scalar counts. The decoder includes an exact final equality test and gives polynomial-time membership in this specified family, not membership in the full tournament spectrum.

This is a genuine transfer of statistical-mechanical tools to literal tournament counts. However a submodular optimizer finds extrema, not all intermediate energies. Resetting changes the old central interface. Moreover this particular family has only O(X^(log(2)/log(13))) possible counts below X, so it cannot supply all odd targets. The result explains an infinite family fully while proving the limit of that route alone.

## 5. Homometry turns nonlinear changes into linear differences

Old repo work on cyclic marker homometry and phase reconstruction, rather than an obvious Chan–Pak source, supplies the wildcard bridge. In a cyclic Cayley tournament, equal-cardinality cuts with equal cyclic autocorrelation yield identical full interface packets.

For two such cuts, toggle any common subset of sites on which they agree. Each resulting count retains quadratic interactions, but those interactions cancel in the difference between the two replicas. A positive-coefficient Laurent polynomial records the exact possible differences together with the cardinality change. Thus subset-sum methods act on genuine physical central switches, with admissibility retained.

The concrete13-vertex cyclic core example gives14-vertex counts24,389,195 and24,389,151 with the same non-H packet. Their difference44 is the sum of local changes−2390 and2434. This is a useful cancellation witness, not a count44 or a unit switch. Paley symmetry supplies a sharp failure case: a constant interaction matrix makes the agreeing-site phase changes powerless. Prime-order local-field injectivity does not remove that restriction.

The labelled-layer interaction inverse is also exact. Some mean and variance identities are inherited, and are credited rather than republished as discoveries.

## 6. How the bridges can compose, and what is still needed

The four views separate three distinct tasks:

- **Internal control:** classify all states and exact support within a fixed interface fiber. The automata/semigroup lane does this for a nontrivial infinite physical family.
- **Resolution change:** produce switches not trapped in that fiber's growing mesh. Rerooting and homometric-replica cancellation provide lawful examples, but do not yet supply a universal digit bank.
- **Absolute anchoring and renewal:** place a prescribed target in the actual translated range while retaining a usable interface for the next stage. This is still open.

A precise candidate next theorem is: construct a sequence of physical full-central banks with explicit baselines and complete bounded digits, or an equivalent correlated selector, whose two-sided semigroup intervals overlap at every necessary scale and whose interfaces admit the next construction. A reset is useful only if its output can be fed into that theorem with all coordinates accounted for.

A possible roundabout route is to use homometric replicas to engineer cancelling local fields inside a carry-recognizable family, then use the exact flow selector to enforce the common interface and the two-sided semigroup criterion to prove interval filling. This is a research proposal, not a proved compatibility of the current constructions. The current examples have different interfaces; identifying their scalar energies alone would lose the very information needed for composition.

The session therefore moves from “many operations and encouraging gaps” to several complete mathematical descriptions with explicit interfaces, while leaving the decisive interoperability theorem explicitly visible.
