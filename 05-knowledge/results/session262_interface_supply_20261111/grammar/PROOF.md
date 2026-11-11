# Nonlinear reusable sources: a Noetherian certificate and finite initial-case compiler

2026-11-11. **PROVED:** the physical polynomial register grammar, terminating equality-guard compiler, finite inverse-cover reduction, and a one-initial-case nonlinear source service. **FINITE-EXACT:** the supplied new controls. **CONDITIONAL:** completeness from a finite physically inhabited controller certificate. **OPEN:** such a certificate for every allowable Chan–Pak count. No inherited audit or claim of historical novelty is made.

## 1. What changes, and what does not

Session261 `wildcard/PROOF.md` §7 and `FOURTH_WAVE.md` §4 certify linear finite controllers by closing defect rowspaces. `FIFTH_WAVE.md` proves severe recurrent length restrictions when each edge appends a fixed word. The new question is different: **can a physical operation copy its current source block, rather than append a fixed block?** The operation is then nonlinear in the complete context matrix, and a fixed-word spectral argument no longer describes it.

The result here is a finite all-continuation certificate for that larger grammar. Its certificate terminates by the ascending-chain condition for polynomial ideals, not by a fixed-dimensional rowspace bound. Separately, a residue compiler turns affine target coverage into finitely many explicit initial cases. Neither step substitutes “every missing integer has a source” for the actual obligation.

Closest mechanism: the complete positive two-sided context carrier. Hostile: an algebraic point satisfying the response equations need not have a physical source. Corrected near miss: immediate correct digits do not imply reusable output interfaces. Least-used sidecar: the literal word-register program and its nonempty-block address. Live board: physical copying / future equality guards / inverse target cover / initial sources / domination budget.

The marked reselection in Session261 `constructors/MARKED_RESELECTION.md` remains important: equal H/W/U tensors can conceal different future port menus. No arbitrary reselection is installed as a polynomial map here. It must come with a sufficient rooted response carrier and a proved coordinate map, or be kept as a distinct literal operation with its own certificate.

## 2. A nonlinear operation with literal tournament semantics

Use the two-apex spine grammar with core a_1→a_2 from Session261. A word w specifies the vertices outside the core, their forward spine arcs, and their two core incidences. Let M_w be its 36-dimensional nonnegative transfer matrix, in column convention. Thus

    M_(uv)=M_v M_u,        H(w)=L M_w v_0.

A physical **copy operation** replaces w by ww: take a disjoint copy of its spine vertices, retain their core incidences, and orient all first-copy to second-copy arcs forward. The core is shared, not duplicated. Its context update is exactly

    X → X².

This is not a formal squaring of H, nor an ordinal sum of the whole tournament. Usually H(ww) is not H(w)². It is a quadratic operation on the complete context.

More generally retain finitely many word registers w_1,...,w_r. Each instruction replaces a register by a specified concatenation of constants, old registers, and their converses. All right-hand sides use the old register values. Every occurrence of a register means a disjoint copy of its spine, with the prescribed forward cross arcs. Words are actual finite objects even when their descriptions are shared as a directed acyclic graph.

Let w*=phi(reverse(w)), with phi(0)=3, phi(3)=0, phi(1)=1, phi(2)=2. Retain X_i=M_(w_i), Y_i=M_(w_i*). Concatenation obeys

    X_(uv)=X_v X_u,       Y_(uv)=Y_u Y_v;

converse swaps X and Y; constants have known matrices. Thus every instruction is a **nonnegative integer polynomial map** on a fixed finite carrier, of dimension 2r·36². The physical word program is the realizability sidecar. It distinguishes matrices with identical responses from distinct words, without requiring the words to be reconstructed from the matrices.

For a nonnegative output and a fixed target N, saturation at N+1 commutes with all these polynomial maps. Therefore the language of instruction sequences producing exactly N is regular and has an effective finite state graph. This extends the earlier positive-linear exact-target decision theorem to physical nonlinear copying. It does not make simultaneous polynomial equalities decidable by the same capping method: signed differences of two large coordinates are not recovered by saturation.

## 3. Noetherian future-guard theorem

Let Q be a finite control graph. State q has carrier x in Q^d, polynomial output h_q(x), and a finite list G_q of mandatory polynomial equality guards. Each edge e:q→p has a specified physical operation whose carrier map F_e is polynomial over Q, and a desired polynomial scalar law P_e. Set

    δ_e(x)=h_p(F_e(x))−P_e(h_q(x)).

Start with I_q^(0)=ideal(G_q, {δ_e:e leaves q}). Iterate simultaneously

    I_q^(j+1) = I_q^(j) + sum_(e:q→p) F_e^*(I_p^(j)).

Here F_e^* substitutes F_e(x) into a polynomial. A finite generating set for I_p suffices: substitution is a ring homomorphism, so images of its generators generate the pulled-back ideal as an ideal in the source ring.

**Theorem (PROVED).** This process terminates. A physical seed x in state q annihilates the terminal ideal I_q if and only if every finite legal control path from it satisfies every equality guard and every prescribed scalar law at every step.

**Proof.** Each I_q^(j) is an increasing chain of ideals in the Noetherian ring Q[x_1,...,x_d]. Each chain stabilizes; finitely many states give a common stabilization stage. At that stage the ideals are closed under every edge pullback. By induction, all future guard and defect polynomials belong to the relevant initial ideal, proving sufficiency. Conversely the iterates are generated by finite-path pullbacks of the original guards and defects. If every future guard and defect vanishes at a seed, every such generator, and hence the whole ideal, vanishes there. This proves necessity. Gröbner bases give an exact terminating implementation of each equality/stabilization test. No useful complexity bound is claimed.

A shorter proof certificate need not run this full closure. Supply ideals J_q by generators and exact membership identities showing:

1. every local guard and defect lies in J_q;
2. g∘F_e lies in J_q for every generator g of J_p and edge e:q→p;
3. each literal initial source annihilates its J_q.

These three finite checks certify all continuations. The theorem is stronger than a rowspace certificate because copying produces higher-degree defects and pullbacks. For linear homogeneous maps/linear guards it specializes to the inherited rowspace mechanism.

**Guard boundary.** Positivity, nonempty modules, designated-root incidence, and port accessibility are structural physical guards, not silently encoded by the polynomial ideal. The instruction must preserve them by its literal construction, or carry a separately proved invariant. Arbitrary inequalities are not covered by the Noetherian equality theorem. Nor does a rational zero of I_q provide a literal initial source.

## 4. A finite arithmetic compiler that exposes every initial case

This second theorem is independent of how the physical guard certificate was obtained. Fix B≥1, a modulus M≥1, and residue sets R_q⊆Z/MZ. The intended labelled targets in state q are

    D_q={n≥B: n mod M belongs to R_q}.

Each physical edge e:q→p has certified scalar law

    h → a_e h+b_e,       a_e≥2 integer, b_e integer.

Assume the preceding ideal certificate preserves its source class and all future equalities. Define

    C=max(B+1, max_e(a_e B+b_e), max_e(1+floor(−b_e/(a_e−1)))),
    L=lcm(M, {a_e M:e an edge}).

Perform the following **finite** test. For every state p and residue t modulo L with t mod M in R_p, require an incoming edge e:q→p such that

    t≡b_e (mod a_e),
    (t−b_e)/a_e mod M belongs to R_q.

Division here is after the stated divisibility test. Its result modulo M is independent of the representative t because L is divisible by a_e M.

**Finite initial-case theorem (PROVED).** If that residue test passes and a literal source in the required certified interface is supplied for each labelled n∈D_q∩[B,C−1], then every labelled target in every D_q has a literal certified source.

**Proof.** For n≥C, the finite residue certificate chooses an incoming edge. Its parent m=(n−b_e)/a_e is an integer in the correct residue set. The definition of C gives m≥B and m<n. Strong induction first supplies the parent in state q; execute its physical edge to produce n in state p. Targets below C are exactly the listed initial sources. Replaying the inductive derivation backwards is an explicit construction algorithm. Since slopes are at least two and the offsets are fixed, sufficiently large targets halve up to a bounded additive error, giving O(log n) edges before the finite initial region.

The seed obligation is finite and explicit, not a hypothesis over every target. Labels matter: a seed of scalar count n in the wrong interface is not an initial case. The residue test says nothing about those physical initial sources or the actual edge identities; they must be supplied separately.

For the desired binary rules 2h±1, B=23, M=2, R={1}, the compiler gives C=47 and L=4. The complete initial list is

    23,25,27,29,31,33,35,37,39,41,43,45.

This reproduces the earlier binary interval reduction as a special case, but allows mixed integer slopes, state-dependent residue domains, and multiple incoming source interfaces. With only the plus edge, residue 1 modulo 4 fails. An infinite plus ray cannot be mislabeled uniform odd-count supply by this certificate.

**CONDITIONAL Chan–Pak consequence.** A finite physical controller satisfying this test for all odd n≥23, plus literal witnesses for each allowable smaller count, would settle the spectrum. The arithmetic certificate is supplied here; the requisite physical digit controller and its twelve initial sources are NOT supplied. More general mixed-radix controllers can now be tested without guessing a seed interval.

## 5. A complete nonlinear source service with its actual initial case

A useful nonvacuous example is the inherited literal family

    w_m=(1,2,0,2^m),        m≥1,
    H(w_m)=15·2^m−1.

The all-size count identity is inherited from Session261 `wildcard/PROOF.md` §6; it is not re-proved from sample values here. Regard the final block of letter2 vertices as a designated nonempty source register. Define physical operations C_0,C_1 by replacing that block by two copies, then appending e additional letter2 vertices:

    C_e: m→2m+e,           e∈{0,1}.

They preserve the actual nonempty designated register and its core incidence. Starting with the **single initial source m=1**, word (1,2,0,2), every m≥1 is constructed by reading the binary digits of m after its leading1. This supplies every member of the family, including its initial case, by a reusable physical program with O(log m) instructions. Materializing the tournament still requires its m+5 vertices; compressed certificates do not erase output size.

Write x=2^m and h=H(w_m). The exact scalar-carrier maps are

    F_e(x,h)=(2^e x², 2^e(h+1)²/15−1).

Their denominators are harmless on the physical class. The one-generator ideal g=h−15x+1 has the exact membership certificate

    g(F_e(x,h)) = [2^e(h+1+15x)/15]·g(x,h).

The literal initial source has (x,h)=(2,29), hence g=0. This is a genuinely nonlinear all-continuation certificate with both operations, not a fitted recurrence and not merely an ambient algebraic line. It also has its complete target-index initial-case proof, with no m=0 supplied by wishful extrapolation.

The forbidden zero-module point is a sharp hostile: (x,h)=(1,14) satisfies g=0, but the actual word (1,2,0) has H=13. Thus the nonempty register address and literal initial source do real mathematical work. This is the same distinction between a formal orbit and its addressed initial component that appears in the randomly recovered polynomial and Collatz sources below.

The nonlinear service gives **no new scalar coverage beyond that inherited ray**. It establishes that a physical source can be reused by operations outside the fixed-word linear controller, and that the stronger finite certificate is inhabited. Copying a variable block avoids the hypotheses of the fixed-appended-word spectral theorem; it does not refute that theorem or magically provide ± binary count digits.

## 6. Where the universal domination-defect inequality enters

Take as inherited the universal inequality

    d(G_z)−1 ≤ 2P_even(z) ≤ H(T)−H(T−z).

For a source of count h and a proposed affine count operation h→ah+b, realize the copying/insertion instruction as successive additions z_1,...,z_r. The original source remains induced throughout this sequence. Telescoping gives the necessary exact budget

    sum_i [d(G_(z_i))−1] ≤ sum_i 2P_even(z_i) ≤ (a−1)h+b.

Count-preserving converse steps can be separated from insertion segments and contribute zero to the scalar telescope. A destructive replacement that does not preserve the previous tournament is not licensed to use this particular telescope without additional work.

This changes the research objective from “does growth produce enough counts?” to “which physical source classes can afford, realize, and return the exact target budget?” The inequality is a one-sided source filter; the polynomial identity is the equality obligation; the invariant ideal is the reusability obligation; the residue test is target coverage; the finite literal list is initial supply. None substitutes for another. In particular, domination lower bounds do not prescribe an exact increment, and a large algebraic invariant set does not supply inhabited interfaces.

The practical next search should use budget violations to discard impossible physical digit candidates before ideal closure, then compute guard ideals only for survivors, then search for the finitely many labelled initial sources required by the inverse-cover compiler. For a more general source operation that is not polynomial on this spine carrier, a new sufficient carrier remains necessary.

## 7. Random historical recovery and transfer ledger

The exact sorted pool, seed262310, pool size6175, and three uniform-without-replacement draws are recorded in `random_sources.json`. All three sources were read; none of their older theorem claims is independently audited or used as a logical dependency.

* **THM-3897**, `01-canon/theorems/THM-3897-f-zero-residual-all-degree-global-emptiness.md`: polynomial solutions need both the equation and the source address. Transfer: an invariant variety is not an inhabited physical source interface. The m=0 point above is the decisive new hostile.
* **THM-1570**, `01-canon/theorems/THM-1570-one-sided-elimination-and-laplace-layer.md`: finite polynomial ideal computations can certify all coefficient values rather than a coefficient grid. Transfer: replace linear future-defect spans by polynomial pullback ideals. The present proof uses Hilbert's basis theorem directly; no GMC claim is imported.
* `05-knowledge/results/inverse_ray_ternary_addresses_20261004.md`: retain an exact parent pointer, a source guard, and the quotient that a residue projection loses. Transfer: a finite residue test chooses an actual smaller target parent, while a literal source program certifies the terminal interface. Exact equality and a finite initial bank replace density or full residue support alone.

Connection map: literal word registers → paired context matrices → polynomial future-guard ideals. Preserved: every specified future output identity. Lost: vertex identities and physical reachability. Restoration sidecar: literal word programs and marked nonempty register addresses. Cheapest test: copy a tiny word and compare X² with its literal doubled-word count.

Connection map: labelled affine edges → residues modulo lcm(a_e M) plus threshold C. Preserved: exact inverse-parent integrality, interface residue, and strict descent. Lost: actual graph sources. Restoration sidecar: every listed literal initial source and the all-continuation interface certificate. Cheapest hostile: remove the minus digit and observe the missing 1 mod4 class.

## 8. Reproduction and result boundary

From repository root:

    python3 05-knowledge/results/session262_interface_supply_20261111/grammar/certificate.py

`receipt.json` records 85 new square-context checks, six independent literal subset-DP path counts for copied blocks, the exact one-generator polynomial identities, 1000 binary-index decodings, the complete twelve-entry arithmetic initial list, the plus-only failure, and the zero-module phantom. These finite controls check the new operations; the displayed arguments establish the all-size statements.

The main new deliverable is a **terminating nonlinear interface certificate combined with an explicit finite initial-case compiler**. Its nonvacuous physical example supplies one whole indexed family from one genuine initial source. A uniform odd-count physical controller, its inhabited seed bank, and the full Chan–Pak theorem remain OPEN.
