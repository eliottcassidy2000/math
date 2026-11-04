# Complete termination of the singleton–triangle return grammar

2026-11-04. **PROVED all-length inverse and normalization; FINITE-EXACT complete return certificate.** The natural ordered-component enlargement of Session183's terminal-C3 construction does **not** close under both ternary digits. Its entire clock-uniform return tree from the common TT3 entry has seven normalized states and maximum depth four, even allowing arbitrary fixed integer clock shifts. This is a stronger scoped obstruction than the previous transitive-only theorem, not a global spectrum obstruction.

## 1. Sources, scope, and the actual grammar

Use Session182 `induction/PROOF.md` for the seven-coordinate ordered-sum law and root attachment, Session183 `mining/PROOF.md` for its terminal-triangle +2 construction, and `multiport_defect_grammar_20261103/transport/TWO_DIGIT_ENTRY.md` for the two original siblings. These sources were checked against main at `0e8bdcdb77f49045e65ece0334b06a5a17d6d8b1`.

The older `session163_ternary_triples_20261103/nonuniform/NOTE.md` already proves the general exceptional-interface bound 3^max(k−2,0) for k triangle modules and one nonuniform exceptional vertex. We explicitly **do not claim that coarse divisibility bound as new**. Section 5 gives its exact leading coefficient for this ordered grammar; Sections 2–4 give the new complete inverse and return-tree termination.

A word consists of blocks S_a (one vertex, a=0 or1 ports) and C_a (a directed triangle, a=0,1,2,3 ports), with every earlier block beating every later block. The marked vertices are actual subsets. All subsets of the same cardinality in a directed triangle are equivalent by its cyclic automorphisms; thus these six symbols cover every word of singleton/C3 strongly connected components, up to marked isomorphism. Retain any concrete labelled representatives when constructing the graph.

Let B be such a word and W its port set. Attach a nonempty TT_t prefix A and a root r with A→B, W→r, and r→A and r→(B\W). The full tournament is strong precisely when the last block has at least one port. Write its inherited signature as (h,p,F,G,J,e,s), and c=h−J. Then

    H_t(B,W)=c+2^t F.

A target with a fixed integer clock shift q has H_(t+q)(B')=3H_t(B)+d at two distinct legal clocks, hence all legal clocks, exactly when

    2^q F'=3F,          c'=3c+d,          d in {2,4}.       (1)

Only this clock-uniform notion is classified. A success at a single chosen clock, a suffix with a larger strong component, a different attachment architecture, or multiple nonuniform root interfaces lies outside the result. No native Hamilton-complement program or old-name-preserving operation is asserted.

## 2. A complete inverse at fixed F, with no word-length cutoff

The literal component signatures are

    S_a: (1,2,a,a,0,a,a),
    C_a: (3,12,3a,3a,a(a−1),a,a).

These follow by enumeration inside one vertex or a directed triangle. The inherited sum law, for signatures (h,p,F,G,J,e,s) and (H,P,f,g,j,E,S), is

    (hH, pP, FH+pf, GP+hg, JH+hj+eS+Gf, hE, sH).       (2)

Let a word contain k triangles. Its h is 3^k. Define each block's binary width w_i to be1 for S and2 for C, and put q_i=sum_(j<i) w_j. Dividing the F-coordinate recurrence by h gives

    F=3^k u,            u=sum_i a_i 2^(q_i).             (3)

Indeed p_i/h_i equals2 for S and4 for C, while F_i/h_i=a_i in both cases. Each S digit fits one binary position and each C digit fits two; all allowed digits, including zero, are present. Thus (3) is a mixed-radix expansion with no carries between blocks.

**Inverse theorem.** For a prescribed positive F, every strong-attachment word with that F is obtained exactly once by the following procedure:

1. Choose 0<=k<=v_3(F), and put u=F/3^k.
2. Read the binary expansion of u from low to high, partitioning it into width-one and width-two chunks, using exactly k width-two chunks.
3. Replace a width-one chunk of value a by S_a and a width-two chunk by C_a.
4. The last chunk must be nonzero. It may be width two with a single extra high zero bit; no other high zero padding is allowed.

Proof: (3) recovers all digits of any actual word. Strongness excludes a zero final block, so its effective binary length is bitlength(u) or bitlength(u)+1; the latter occurs only for final C_1. Conversely every such partition constructs the actual ordered sum, has F by (3), and has a terminal port. Chunk positions and values recover the word uniquely. ∎

In particular, this is a finite, complete all-order decision procedure for any prescribed pair (F,c): enumerate these words, then use (2) to test c. A returned word has at most bitlength(u)+1+k vertices. The search is not stopped at an experimentally chosen graph order.

## 3. Removing arbitrary clock shifts

For a digit output c'=3c+d, c' is nonzero modulo3. Consider any candidate target word. Its initial zero-port blocks form an ordered prefix P. By (2),

    c_(P→D)=h_P c_D.

Consequently P cannot contain a triangle: that would make c' divisible by3. Its entire initial zero-port prefix therefore consists of S_0 blocks.

Delete those initial S_0 blocks to obtain a normalized word whose first block is marked. Deleting each such block divides F by2 and leaves c unchanged. In a normalized word, (3) shows that v_2(F) is0 or1: a first S_1, C_1, or C_3 makes u odd; a first C_2 makes u exactly2 modulo4.

Let F_0 be the odd part of3F_source. Equation (1) forces the same odd part for F_target. Thus **all normalized targets for either digit occur in the two finite inverse lists F_target=F_0 and F_target=2F_0**, restricted to first block marked. Restore any number a of initial S_0 blocks afterward, compensating by the integer clock shift required by (1). All sufficiently large input clocks remain legal.

This proves completeness over every word length and every fixed positive or negative clock shift. It also justifies normalizing after each step of an entire lineage: extra initial S_0 blocks alter only the dyadic scale of F, not c or the two normalized successor lists. The normalized state retains its literal marked graph; this is not a claim that arbitrary equal-count graphs can be interchanged in a native program.

## 4. The complete common-entry tree terminates

Start at the actual marked TT3 word S_1 S_1 S_1, for which (F,c)=(7,−5). Let a digit word record successive requested d values. The exhaustive inverse certificates in `candidates.json` and `receipt.json`, computed by (2)–(3) and the normalization theorem, give exactly these states:

| Digit word | Normalized component word | (h,F,J,c) | H at t=1 | Permitted next digits |
|---|---|---|---|---|
| empty | S1 S1 S1 | (1,7,6,−5) | 9 | 2,4 |
| 2 | S1 S1 C1 | (3,21,16,−13) | 29 | 2 |
| 4 | S1 S0 S1 S0 S1 | (1,21,12,−11) | 31 | 2 |
| 22 | S1 S0 S1 C2 | (3,63,40,−37) | 89 | none |
| 42 | S1 S0 C1 S1 | (3,63,34,−31) | 95 | 2 |
| 422 | S1 S0 C1 C1 | (9,189,100,−91) | 287 | 2 |
| 4222 | C1 S1 S1 S1 S1 S0 S1 | (3,567,274,−271) | 863 | none |

For candidate odd F values21,63,189,567,1701, the normalized inverse lists have respectively5,10,18,29,50 words, including possible doubled-F words. Each list is finite by proof, not by a search cutoff. At each state both target constants3c+2 and3c+4 are compared against **every** entry. Repeated source F values lead to280 comparisons across the fourteen digit tests, with112 distinct actual candidate words.

All accepted successors in this tree have F'=3F with no clock shift in their displayed representatives. Each accepted normalized successor is unique within this grammar. Therefore the complete address language is

    {empty, 2, 4, 22, 42, 422, 4222}.                    (4)

In particular, neither original sibling has a +4 successor anywhere in the singleton–triangle grammar, even after arbitrarily many added components, arbitrary port changes, or arbitrary fixed clock shifts. The two longer literal all-clock +2 chains at t=1 are

    9 → 29 → 89,
    9 → 31 → 95 → 287 → 863.

They exhibit actual strong representatives, but neither chain is a recursively reusable two-digit construction. The fact that triangle count decreases from two at422 to one at4222 is allowed by this broad census; imposing preservation of old induced triangle components could only remove options. The theorem already allows arbitrary replacement within the entire grammar.

The finite arithmetic certificate, together with the all-order inverse theorem, proves termination at every depth. It does not require induction from an unproved finite-menu heuristic. Larger strong suffix components or non-clock-uniform changes remain viable routes and are not excluded.

## 5. Exact leading 3-adic obstruction and one-triangle eligibility

For a word with k>=2 triangles define

    A=sum a_i a_(i+1),

where the sum is over **adjacent component positions both occupied by triangles**. Then the stronger exact leading law is

    3^(k−2) divides H_t,
    H_t / 3^(k−2) = −A (mod3).                         (5)

In particular v_3(H_t)=k−2 exactly when A is nonzero modulo3; if A=0 modulo3 there is at least one additional factor3. This coefficient is independent of the clock.

Proof. Resolve the physical joint count J by the components containing its first endpoint and second start. If both are in the same component i, its contribution is J_i times the product of all other h_j. If the first endpoint is in an earlier component i and the second start in a later component, these components must be adjacent, with contribution e_i s_(i+1) times all other h_j; otherwise an intervening component could belong to neither path. In the opposite order, for endpoints in j<i, the contribution is

    G_j F_i (product_(j<l<i) p_l)(product_(l<j or l>i) h_l).

These cases partition the actual two-slot covers. Every opposite-order term is divisible by3^k because F_i,G_j,p_l,h_l each retain the triangle factor3. A same-component term is divisible by3^(k−1). An adjacent-component term loses two factors3 only when both components are triangles, giving exactly3^(k−2)a_i a_(i+1). All other adjacent terms are divisible by3^(k−1). Finally h and F are divisible by3^k, so H_t=h−J+2^tF has leading residue (5). ∎

For k=1, let a be the triangle's number of ports, and let b,c be the port bits of its immediately preceding/following singleton components, using zero at a missing boundary. The same decomposition gives

    H_t = −a(a−1+b+c) (mod3).                          (6)

Thus the triangle must be partially marked for a count coprime to3. With a=1, at least one adjacent singleton must be marked. With a=2, the two adjacent singleton positions must not both be marked. For k=2, (5) says coprimality to3 holds **if and only if** the two triangles are adjacent and each has one or two ports. With k>=3 every count is divisible by3. These are complete mod3 eligibility statements within the grammar, not sufficient conditions for any prescribed prime or7p value.

The generic k−2 lower bound is inherited from S163. The exact coefficient and its cancellation criterion specify when the old capacity bound is sharp for this ordered grammar. For example, k consecutive C_1 blocks have A=k−1, so the bound gains at least one extra factor precisely when k=1 modulo3 (for k>=2). No higher exact valuation is asserted in the cancellation case.

## 6. Reproduction and verification boundary

Run `python verify.py` and `python -O verify.py`. All checks use explicit exceptions, and output agrees in both modes.

- All882 component words through suffix order7 have all seven signature coordinates checked by independent directed subset DP, not by assembly.
- All492 of these words with terminal ports round-trip through the complete inverse.
- All72 terminal-port words through suffix order5 have actual TT1 root attachments independently counted and checked strongly connected.
- All112 distinct normalized candidates used in the complete tree have their full signatures independently recomputed on their literal graphs.
- Each of the seven accepted states is checked as a literal strong full attachment at clocks1 and2.
- The leading coefficient is checked on112 small words and87 additional words with2 through30 triangles.

`candidates.json` retains every normalized candidate, its component word and actual signature, so a reviewer can inspect rejected constants rather than trusting a hidden search. `receipt.json` retains the whole return tree and candidate counts. None of these computational controls replaces the all-order inverse and normalization proofs.

**Global status:** all odd Hamilton-path counts other than7 and21 are still not proved realizable. This result closes one natural nontransitive return language and supplies explicit longer +2 chains. A complete proof must escape this language or abandon clock-uniform returns, and separately prove entry/coverage and any requested native-operation legality.
