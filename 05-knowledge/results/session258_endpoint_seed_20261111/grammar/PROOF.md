# Target-saturated seed compilation: a finite reachability bridge

Status: **PROVED two-positive-probe seed equivalence and finite-reachability reduction**, **FINITE-EXACT implemented one-apex controls**. No universal odd-count coverage claim. The 83-state cutwidth-two compiler is inherited from THM-4878; this packet proves its target-dependent finite quotient and implements the complete procedure for the four-register one-apex subgrammar. It does not execute the potentially enormous 83-register search.

## Inheritance and source map

Closest positive mechanism: `01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md`, items 1 and 6: positive integer cut-state matrices, and cutwidth-two attainment through 8191. The maintained source is `05-knowledge/results/spine_grammars_cutwidth_apexes_hankel_seams_20261104.md`. We use its positive **83-state** representation for G(0,2), not its potentially signed **rank-17** compression.

Hostile control: 13 is an actual tournament count but is absent from the one-apex grammar. The universal tournament exclusions 7 and 21 remain exclusions. These distinguish grammar completeness from tournament completeness.

Closest return source: `05-knowledge/results/berggren_target_return_20261111/reset/TARGET_RETURN.md`, §§1–5. Its six-coordinate funnel compiler, ready barriers at 1997, and growing continuation defect are accepted at their stated scopes. Our compiler avoids a requirement that every scalar intermediate be 1- or 2-ready: it reads an entire physical spine word, retaining all cut fragments until final closure.

Natural-recursion source: `01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`; `05-knowledge/hypotheses/HYP-1822-natural-operation-cospans.md`; historical `07-reflections/summand-graph-fermat-zeckendorf.md` (its warning remains authoritative). The map below transfers labelled-operation fibres and shortest ancestry, rather than identifying the additive shadow with tournament realizability.

Zeckendorf source: the preceding historical reflection, §2, and `05-knowledge/results/zeckendorf_tournaments.out`. We use ordinary Zeckendorf uniqueness to give canonical target addresses, not that reflection's speculative forbidden-value numerology.

Machine source: `01-canon/theorems/THM-4858-exact-translations-of-the-goldbach-machines-a-53-symbol-radius-one-cellular-automaton-and-a-568-fraction-fractran-program-from-lengs-25-state-turing-machine.md`, §(a). Its general TM-to-CA construction applies verbatim; its special 53-symbol size belongs to Leng's machine and is not claimed here.

Concept board: target saturation; labelled operation fibres; cut-fragment endpoints; canonical Zeckendorf address; shortest physical witness; halting predicate.

## 0. Main seed upgrade: two positive probes replace the signed equation

**Notation firewall.** In HYP-9658(A), c is an odd **ray slope**, not a Hamilton-path target. Every odd slope is requested, including slopes 7 and 21. In §§1–5 below, the letter c in the general saturation lemma denotes a single numerical output target; equivalently rename it N. None of those output-target exclusions is an exclusion of slopes. The existing bank already covers all odd slopes through 4097; we do not claim to extend that bank here.

For a rooted L with B=L-v and O=N+(v), inherit THM-5126(A):

    A_t(L)=H(L o 0^t)=D+2^t U,  t>=1,
    D=H(B)-W(O), U=U_B(O).

Fix an odd slope c and an exponent j>=0, and put M=c2^j. Then the following is an exact equivalence:

    D=1 and U=M
      iff A_1(L)=1+2M and A_2(L)=1+4M.                 (Seed-2)

Proof: subtract the two equations to get 2U=2M, then substitute in the first. The reverse direction is immediate. Crucially both A_1 and A_2 are actual nonnegative Hamilton-path counts, even if D is obtained by large cancellation. This gives a positive-observable replacement for the signed endpoint constraint.

**Fixed-exponent seed compiler.** Let a regular physical rooted grammar have positive finite-state path transfers for its one-tail and two-tail decorations. Run the two transfers in direct sum (or use one shared transfer and two final contractions when available). Cap all coordinates at B=2+4M. Accept precisely when the two outputs equal 1+2M and 1+4M. The saturation theorem below proves that exhaustive BFS decides and constructs a unit-offset seed of this exact slope and exponent, over every word length in the grammar. A shortest witness is bounded by the finite product-state count. No bound on H(B) or W(O) is needed.

For bounded-cutwidth spine bases, this positive decorated transfer exists by the same finite frontier construction: retain the marked root and the at most two additional tail vertices as finitely many named boundary vertices, and use their fixed orientations in local edge choices. The cut slots and named boundary vertices are bounded, so the set of fragments is finite. This is the general frontier construction, **not** a claim that the inherited 83-state G(0,2) dump already contains these additional boundary labels. Its joint dimension d_seed must be instantiated before a numerical state bound or computation is reported. The seed machine is mathematically specified here, not implemented by our one-apex control.

This resolves a real obstacle: separately capping H and W could map H=1000,W=999 to the same overflow symbol and erase their difference 1. Seed-2 instead compares two observable counts, and positive saturation is exact. The physical output graphs contain the same rooted base and the stipulated tails; two unrelated graphs attaining those two numbers do **not** certify a seed.

The finite problem is indexed by (grammar,c,j). Enumerating j gives a semidecision procedure for a seed of slope c in the grammar, not a terminating absence test unless an exponent bound is supplied. If J(c) is a proved exponent bound, the finite disjunction j<=J(c) is a complete terminating slope compiler. Thus the new uniform proof obligation is explicit: a physical grammar plus an exponent bound and a nonemptiness invariant for these two positive target constraints. Zeckendorf can canonically encode c and j or M; it does not supply the bound.

Finally, a seed produces N=1+c2^s only for s>=j+1. Its finitely many smaller head points require actual path witnesses separately. Even universal seed existence alone does not settle Chan–Pak. The line exceptions N=7 and N=21 belong respectively to slopes 3 and 5, not slopes 7 and 21.

## 1. The target saturation theorem

Let a physical construction language L be regular, with legality controller Q, finite alphabet A, and a positive integer linear representation

    H(w) = u M[a1] ... M[at] f,

where u, f and every M[a] have nonnegative integer entries and dimension d. Each accepted word must decode to an actual tournament, and the displayed H must count its Hamilton paths. These are exactly the inherited cut-state hypotheses, not an assumption that arbitrary nonnegative matrices come from tournaments.

For a fixed target c >= 1 put B=c+1 and pi(x)=min(B,x). Define capped addition and multiplication by

    x (+) y = min(B,x+y),     x (*) y = min(B,xy).

Then pi is a semiring homomorphism N -> {0,...,B}. In particular, zero times B is zero: a large partial-path channel may be killed by a later letter, and it must not poison channels it does not reach.

Coordinatewise saturation commutes with every matrix transition and final contraction. Therefore the deterministic state

    (q, pi(u M[w]))

has at most |Q| (c+2)^d possibilities and recognizes EXACTLY the physical words with H=c, using final test pi(v f)=c. Proof is induction on word length using the homomorphism; because c<B, equality to c is never confused with overflow.

**Constructive corollary.** BFS starting at (q0,pi(u)) terminates, either returning an actual shortest accepted word of value c or proving that this fixed grammar never realizes c at any order. Store the incoming letter and predecessor for each state. Every minimum-length witness has length strictly less than the number of reachable states, since a repeated deterministic state allows removal of the intervening loop while preserving legality and target acceptance. This is a derived witness bound, not a chosen enumeration cutoff.

**Least-resistance corollary.** Reverse graph distance from the accepting states is an integer rank. Whenever a target state is reachable, choose a successor whose distance decreases by one. This gives a terminating target-directed recursion. With positive integer letter costs, Dijkstra gives the least-cost word and the corresponding cost-to-go ranking. The rank is an actual certificate; a heuristic promise of decreasing H is unnecessary and generally wrong.

For the inherited G(0,2) cut-state machine this gives a finite decision procedure with 83 positive coordinates plus a finite legality controller. It permits arbitrarily many reversed arcs. Thus it escapes THM-4878's bounded-reversed-arc counting obstruction. It is not a claim that its worst-case state bound is practical.

## 2. What the quotient preserves and what it loses

Source: physical labelled construction histories with all endpoint fragments. Target: capped endpoint multiplicity vectors and legality state. Preserved predicate: existence and shortest length of an accepted completion with absolute target c. Lost data: exact multiplicities above c, graph isomorphism type, and identities of previously inserted vertices. Restoration sidecar: predecessor word and deterministic physical decoder. All emitted witnesses retain their full actual adjacency data.

Two equal capped states have the same future target language. This is a **target-dependent continuation equivalence**, not equality of tournaments or permission to splice the retained vertices of one graph into another. If a name-retaining surgery is demanded, its marked identities must be added to the controller and its finite-state status proved separately.

The cap must be applied to positive path-fragment coordinates. Applying it to the signed six-state formula in TARGET_RETURN §3, or to the rank-17 linear representation, is not justified: large positive and negative terms can cancel back below c. Also, pruning a whole state because one coordinate exceeds c is unsound; that coordinate might subsequently be multiplied by zero. Our generic theorem keeps the overflow marker per coordinate.

This is the exact extra information missing from C_k(T)=1: a finite target machine remembers reachable endpoint channels and their nonnegative multiplicities, not their gcd or only their scalar sum.

## 3. An implemented four-register physical control

Take an increasing chain v1,...,vm and one apex a. Encode vertex vi by bit bi=1 when a->vi, and 0 when vi->a. Every directed cycle uses a, so two directed cycles cannot be vertex-disjoint. By the inherited odd-cycle collection identity,

    H=1+2 C,

where C counts odd directed cycles. An odd cycle is uniquely an even-length increasing chain subsequence whose first bit is 1 and last bit is 0, closed through a.

Let o,e count all odd/even positive-length increasing subsequences starting at a 1 (without a last-bit restriction). Scan bits left to right. The physical response (h,o,e,1), initially (1,0,0,1), has updates

    b=0: (h,o,e,1) -> (h+2o, o+e,   e+o, 1),
    b=1: (h,o,e,1) -> (h,    o+e+1, e+o, 1).

The newly counted cycles are precisely the old odd subsequences followed by a zero. These equations provide explicit nonnegative 4x4 matrices. The constant register is fixed, so capped BFS has at most (c+2)^3 nonconstant states. In this particular grammar h is monotone, so a state with h>c may safely be discarded; this extra pruning is not part of the general theorem.

`compile.py` exhausts the capped reachable graph for every odd target through 63, and independently counts each emitted tournament by subset DP. It also compares the transfer with subset DP on every binary word of length 0..8 (511 tests). Normal and optimized Python runs agree. Across all 32 target searches it visits 2135 states. The all-length misses in this one-apex language are

    7,13,21,29,35,41,43,49,55,57,59.

These are not merely misses through length eight. Length eight is only the independent transfer control. The finite saturated closures prove absence for every finite chain length. In particular the 13 control rejects a false inference from successful digit parsing to physical seed realization.

An elementary positive family is the word 1 followed by r zeros, r>=1. Here o after j zeros equals 2^(j-1) for j>=1, and summing the cycle increments gives H=1+2^r. Thus this compiler emits the known infinite dyadic family and general mixed words without any ready-return assumption. This family is inherited, not claimed as a new spectrum result.

## 4. Zeckendorf input as a literal finite-state interface

Write c=sum_i d_i F_(i+2), with F2=1,F3=2, d_i in {0,1}, and d_i d_(i+1)=0. The least-significant-first address is unique after removing leading zeros at the high end. A finite controller (p,r,z) stores previous digit p, position r mod3, and parity z. Reading d requires pd=0 and performs

    (p,r,z) -> (d, r+1 mod3, z XOR (d * epsilon_r)),
    epsilon_0=1, epsilon_1=0, epsilon_2=1.

This is an explicit state map because Fibonacci weights have parity pattern 1,0,1. A fixed-target compiler first decodes this word by exact integer arithmetic and then constructs its capped physical graph. `compile.py` implements both maps and checks the parity interface on every target.

The arithmetic address has O(log c) digits, but the endpoint machine still depends on c. This distinction is important: uniqueness of a short number representation does not produce a uniform finite-state tournament selector. The open synthesis target is now precise: exhibit a uniform construction policy on digit addresses, with a proved invariant relating its endpoint registers to the remaining Fibonacci weight and with a physical decoder. A register invariant or a well-founded ranking may be enough; a finite unweighted state set is not assumed.

## 5. An exact machine question, rather than a universality analogy

Define SeedGap_G as follows. Enumerate canonical Zeckendorf words by numeric value (or enumerate odd integers and encode them), skip 7 and 21, build the finite capped graph for target c, and halt with that address iff the graph has no accepting state. Every individual iteration terminates by §1. Thus

    SeedGap_G halts iff G misses some allowed odd target.

For G=G(0,2), nonhalting would prove this fixed-grade coverage conjecture, hence Chan–Pak. Halting would refute only the grade's completeness, not Chan–Pak; 13 in the one-apex control demonstrates why the distinction is necessary. This construction is an effective algorithm and a mathematical halting equivalence, not a supplied small explicit TM transition table.

Compile the algorithm to any standard single-head TM M with alphabet Gamma and states Q_M. Apply the inherited radius-one rule on Gamma union (Q_M x Gamma) union {HALT}: the head cell writes M's symbol; the neighbour in M's move direction acquires the new head/state; HALT persists. A valid single-head configuration maps step-for-step to M. Consequently this CA emits HALT exactly for a missing target in the chosen physical grammar. The output includes the Zeckendorf target address and the exhausted finite graph certificate, so the obstruction is recoverable rather than merely a flashing cell.

This connects the repository's CA/TM translations to **the exact seed predicate** while preserving the grammar quantifier. It does not use a universal CA as evidence that any prescribed H can be built.

## 6. Concrete next obligation

Instantiate the inherited 83-state dump, its legality controller, and physical word decoder, then target the first difficult seed/ready barriers (1997, 5993, and successors). For modest c, reduce by target-language equivalence or backward coaccessibility before BFS. The positive 83-state coordinates must survive any minimization, or the capping proof must be replaced.

The theorem already guarantees that a stalled fixed-grade target search has a finite certificate to seek. What remains open is a uniform argument that every allowed c has an accepting state. An independently verified digit-policy invariant would close that obligation; further finite censuses alone would not.

## Reproduce

    python3 compile.py
    python3 -O compile.py
    python3 seed_probe_check.py
    python3 -O seed_probe_check.py

The seed-probe control checks 1,681,440 admissible integer affine-law cases, including a large-cancellation hostile. It does not enumerate physical seeds.

The two compile runs write the deterministic `receipt.json`. No third-party packages are required. Source pins in this note use repository paths at the session's read snapshot 13268e7459bb; historical scope warnings are retained, and this packet does not modify those sources.
