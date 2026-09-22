# The extended Collatz graph E: one backward braid step, a 3-adic Terras theorem, and the cycle census

**Status: PROVED (scoped) — leaf identity as an index-shift of the inherited inverse fibre, transient forest, inverse-move equivalence for Q2, residue drift table, 3-adic hostile family, greedy-word determination mod 3^(J+1), exact residue Markov chain, E[2^{K_J}] = 3(7/3)^{J-1}, density-1 greedy stopping time, the pn+1 hierarchy with its general criticality dichotomy and the p = 7, 17 dead-exit probabilities; FINITE-EXACT — Q1/Q2 to 10^6, SCC data to 10^5, the 74-cycle bounded census and its GLOBAL completion for every cycle length <= 26, the signed side to 10^6, the p = 5 Chernoff certificate; OPEN — Q2 itself, the one-SCC conjecture C_E, h = ceil(a log2 3) beyond length 26, Q1_5; REFUTED — "greedy suffices" on the minus side (witness m = 4); SCOPE — no map found to THM-4139/4146, THM-3341/3333/1745.** No Collatz convergence, cycle-uniqueness, or Goldbach statement is claimed; the drift table, the transient forest and the Terras-type mechanism carry no novelty claim (classical or elementary). Session collatz-mod6-20260917 (machine mac-mini), lane `extended_collatz_scc`; script recovered from the agent transcript on 2026-09-21, audited, extended (section S8) and regenerated on that date.

## Inheritance and concept board

This lane takes the user's question "what arrows would appear if evens also went to 3n+1?" as the definition of an object rather than a remark. The closest proved mechanism is the inherited inverse-fibre theorem of [arithmetic_braids_20260917_collatz.md](arithmetic_braids_20260917_collatz.md) (sections 1–2: three-row typing, fibre formula (B1) n_j = (2^(h0+1) 4^j u - 1)/3 with R(n) = 4n+1, the exact word condition (B5), the 2-adic hostile block n = 2^(L+1)-1, and the fact that the source (4u-1)/3 of the lowest target is 1 mod 8) together with Terras' stopping-time theorem for the forward map (CITED below). The canonical hostiles are the 3-adic family m = 3^j+1 (section 2), the start m = 4 on the minus side (section 6), and the 5n+1 control E_5 whose greedy chain has the same density mechanism but extra cycles (section 7). The corrected near misses are, first, the recovered draft's carry formula (exponent of 3 reversed; section 4), its k-word for the m = 4 rescue, and its "peak" that silently mixed E-nodes with compound values (section 3); second, the bounded census figure "74", both of whose caps bind (section 5). The least-used sidecar is the tilt eigenvalue rho_p of the greedy chain, which is what separates 3n+1 from 5n+1 when the density theorem cannot (section 7).

Also inherited and cited, not re-derived: [arithmetic_braids_20260917_summand.md](arithmetic_braids_20260917_summand.md) section 6 (the "unhalved extension": every target v >= 7 with v = 1 mod 6 has exactly one even 3n+1-predecessor (v-1)/3, and path-existence versus all-path convergence are different statements in the nondeterministic graph) and section 5 (the 5n+1 cycle 13->33->83); [arithmetic_braids2_20260917_signed_cycles.md](arithmetic_braids2_20260917_signed_cycles.md) section 5 (the b = -1 positive odd cycles (1), (5,7), (17,25,37,55,41,61,91), which are the odd skeletons of the three E_- cycles of section 6, and the minimal parameter q = |Delta|/gcd(B,|Delta|)); [arithmetic_braids2_20260917_inverse_completion.md](arithmetic_braids2_20260917_inverse_completion.md) (every finite halving word is realized in every basin — the deterministic-inverse counterpart of the greedy-word determination of section 4, and the reason a finite residue word cannot decide Q2); and the blueprint audit [collatz_blueprint_20260921_synthesis.md](collatz_blueprint_20260921_synthesis.md) section 4 / [collatz_blueprint_20260921_energy.md](collatz_blueprint_20260921_energy.md) (a natural density of sources is not a time average on each orbit; a density-1 statement leaves a quantifier gap before an every-integer statement — exactly the gap between Theorem 4.4 and Q2 here). Canon [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md), [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md), THM-3341, THM-3333 and THM-1745 were read; section 8 records that no map was found (SCOPE).

The object is the nondeterministic digraph on the positive integers

    E:   n -> n/2   (n even),        n -> 3n+1   (all n).

Its arrows are the reversed doubling forest D of the summand note together with the strict summand arrows n -> n+(2n+1). E contains the deterministic Collatz map C as a sub-relation; the added arrows are exactly n -> 3n+1 for even n. E_- is the same with 3n-1; E_p with pn+1 is the hostile control.

## 1. The new arrows are one backward step of the inverse-fibre braid (PROVED; mostly inherited)

**Theorem 1.1 (leaf identity).** Let v = 6j+1 >= 7. The E-predecessors of v are exactly 2v (halving) and 2j = (v-1)/3 (new). Let n0 = (4v-1)/3 be the least odd T-predecessor of v (inherited (B1) with h0 = 1). Then

    n0 = 8j+1 = 1 (mod 8),     3 n0 + 1 = 4v,     R^{-1}(n0) = (n0-1)/4 = 2j.

Hence E extends the fibre formula n_j = (2^2 4^j v - 1)/3 of (B1) to the index j = -1, and to no further index: (2^2 4^{-2} v - 1)/3 = (v/4-1)/3 is not an integer for odd v. Targets v = 5 mod 6 and v = 3 mod 6 acquire no new predecessor.

*Attribution.* The unique even predecessor (v-1)/3 of each v = 1 mod 6, v >= 7, is the inherited summand note, section 6; that the source (4u-1)/3 of the lowest target is 1 mod 8 is the inherited collatz note, section 1. What is new here is only the reading R^{-1}(n0) = 2j, i.e. that the new arrow is the fibre formula (B1) evaluated at index -1 and that index -2 is never integral.

*Proof.* An n with 3n+1 = v exists iff v = 1 mod 3, and it is even iff v = 1 mod 6 (for even v the predecessor (v-1)/3 is odd and the arrow is already in C). Compute (4(6j+1)-1)/3 = 8j+1, 3(8j+1)+1 = 24j+4 = 4v, (8j+1-1)/4 = 2j. For v = 3 or 5 mod 6, v is not 1 mod 3. The fibre formula at j = -1 is (v-1)/3 = 2j, and at j = -2 it is (v/4-1)/3, not integral since v is odd. QED

So the user's "1 mod 6 chain" is precisely where E adds arrows; the "3 mod 6 chain is never visited" statement is Theorem 1.2.

**Theorem 1.2 (transient forest; elementary, no novelty claim).** No arrow of E (or E_-) enters 3Z from outside; inside 3Z only halving arrows exist. Hence the multiples of 3 are singleton SCCs, and every path starting at 3a leaves 3Z after at most v_2(3a) halvings (a 3n+1 arrow leaves immediately) and never returns.

*Proof.* 3n+1 = 1 mod 3, 3n-1 = 2 mod 3, and n/2 = 0 mod 3 iff n = 0 mod 3. The halving forest is acyclic. QED

Script S1 cross-checks the predecessor sets for all v <= 10^5 (16666 new arrows, one per v = 1 mod 6 in [7, 10^5]) and the transience for n <= 10^5.

## 2. Two reachability questions and the residue drift (PROVED)

Define **Q1**: every n reaches 1 in E (implied by Collatz, since the C-path is an E-path). **Q2 (new)**: 1 reaches every m with 3 not | m.

**Theorem 2.1.** Q2 holds iff every m with 3 not | m can be reduced to 1 by the inverse moves m -> 2m and m -> (m-1)/3 (when integral, any parity) staying inside the non-multiples of 3; equivalently by compound moves

    m -> (2^k m - 1)/3,     k >= 0,     2^k m = 4 or 7 (mod 9).

*Proof.* Reverse the arrows. Forward paths from 1 never enter 3Z (Theorem 1.2), so a reverse path from m to 1 must avoid 3Z, which kills the branch (m-1)/3 = 0 mod 3. A reverse path to 1 cannot end with a doubling (2x = 1 is impossible), so it is a word of compound moves. (2^k m-1)/3 is an integer not divisible by 3 iff 2^k m = 1 mod 3 and 2^k m != 1 mod 9, i.e. 2^k m in {4,7} mod 9. QED

**Lemma 2.2 (drift table, minimal admissible k; elementary).** Since ord_9(2) = 6:

| m mod 9 | k_min | admissible k mod 6 | factor 2^k/3 | result mod 3 |
|---|---|---|---|---|
| 1 | 2 | 2, 4 | 4/3 | 1 |
| 2 | 1 | 1, 3 | 2/3 | 1 |
| 4 | 0 | 0, 2 | 1/3 | 1 |
| 5 | 3 | 3, 5 | 8/3 | 1 |
| 7 | 0 | 0, 4 | 1/3 | 2 |
| 8 | 1 | 1, 5 | 2/3 | 2 |

The parity of k is [m = 2 mod 3]; of the three classes mod 6 with that parity exactly one is excluded (it gives 2^k m = 1 mod 9). Residues {1,2,4,5} produce a result 1 mod 3, residues {7,8} a result 2 mod 3. So 4,7 shrink by 1/3; 2,8 by 2/3; 1 grows by 4/3; 5 by 8/3.

**Theorem 2.3 (3-adic hostile family).** Let v_3(m-1) = j >= 1 exactly (m = 1 mod 3^j, m != 1 mod 3^(j+1)). The greedy word starts 2^{j-1} 0 and the j-th greedy image is m_j = 4^{j-1}(m-1)/3^j. Thus m_j < m for j <= 4, and for j >= 5 one has m_j > m as soon as m > 4^{j-1}/(4^{j-1}-3^j), which every m >= 244 satisfies. The net factor 4^{j-1}/3^j exceeds 1 iff j >= 5.

*Proof.* If j = 1 the residue is 4 or 7 mod 9, k = 0 and m_1 = (m-1)/3. If j >= 2, m = 1 mod 9 forces k = 2 and m_1 - 1 = 4(m-1)/3 has 3-valuation j-1; iterate. At valuation 1 the residue is 4 or 7 mod 9, so k = 0 and m_j = (m_{j-1}-1)/3 = 4^{j-1}(m-1)/3^j. Finally 4^{j-1} < 3^j exactly for j <= 4 (1 < 3, 4 < 9, 16 < 27, 64 < 81, 256 > 243). QED

| j | m = 3^j+1 | greedy values to first shrink | net | 4^{j-1}/3^j |
|---|---|---|---|---|
| 1 | 4 | 4, 1 | 1/4 | 1/3 |
| 2 | 10 | 10, 13, 4 | 2/5 | 4/9 |
| 3 | 28 | 28, 37, 49, 16 | 4/7 | 16/27 |
| 4 | 82 | 82, 109, 145, 193, 64 | 32/41 | 64/81 |
| 5 | 244 | 244, 325, 433, 577, 769, 256 | 64/61 | 256/243 |
| 6 | 730 | 730, 973, 1297, 1729, 2305, 3073, 1024 | 512/365 | 1024/729 |
| 7 | 2188 | ..., 12289, 4096 | 1024/547 | 4096/2187 |
| 8 | 6562 | ..., 49153, 16384 | 8192/3281 | 16384/6561 |

Note the exact images 2^{2j-2} of m = 3^j+1: the growth block lands on a power of two. Greedy from 244 continues and descends at step 6: compound values 244, 325, 433, 577, 769, 256, 85; the E-node peak of this segment is 2308 = 4*577 (a doubling intermediate), the compound-value peak is 769 (script S2, S8.2).

**Typed analogy (2-adic block <-> 3-adic block).** Source: the inherited Collatz block n0 = 2^{L+1}-1 with L consecutive k = 1 steps, T^L(n0) = 2*3^L-1 (inherited (B5)). Target: m with v_3(m-1) = j and j-1 consecutive k = 2 greedy steps, m_{j-1} = 4^{j-1}(m-1)/3^{j-1}+1. Map: reverse arrows and exchange the roles of 2 and 3 (parity word mod 2^L <-> greedy word mod 3^{j+1}). Preserved: forced word prefix determined by a single prime-power congruence; unbounded expansion factor. Lost: determinism (E's inverse is nondeterministic; greedy is one strategy) and the criticality class (section 4). Sidecar: the free choice of k within the admissible classes. Cheapest decisive test (script S8.8): 244 -> 256 versus the accelerated T-orbit 31 -> 47 -> 71 -> 107 -> 161 = 2*3^4-1 (n0 = 2^5-1, L = 4).

## 3. Finite-exact verification and SCC data

**Q2 to 10^6 (FINITE-EXACT).** For every m <= 10^6 with 3 not | m the greedy strategy alone reaches a value below m (statuses: 666664 'below', 2 'one' — m = 2 and m = 4 hit 1 directly — 0 'cycle', 0 'cap'). Chaining descents, the maximal number of compound moves to 1 is 74 and of E-arrows 172, both at m = 984104. **Peak semantics (corrected).** Two peaks must be distinguished. The E-node peak ranges over all nodes of the E-path 1 -> m, including the doubling intermediates 2^k x; its maximum is 150994948 = 4*37748737 = 9*2^24+4 at m = 797162 (peak/m = 189.416). The compound-value peak (the values m_i only) at the same m is 50331649; the printed inverse path of m = 797162 (39 compound moves) climbs to 50331649 -> 16777216 = 2^24 and then falls. The recovered draft quoted the first number while displaying the second kind of path.

**Q1 to 10^6 (FINITE-EXACT).** Every n in [2, 10^6] reaches a smaller value under C; induction.

**SCC data (FINITE-EXACT, iterative Tarjan on E|[1,N]).**

| N | #SCC | nontrivial SCCs | giant size | non-multiples of 3 outside giant | of those <= N/2 | smallest outsiders |
|---|---|---|---|---|---|---|
| 10^3 | 854 | 1 | 147 | 520 of 667 | 202 | 31, 47, 55, 71, 73 |
| 10^4 | 8249 | 1 | 1752 | 4915 of 6667 | 1737 | 383, 511, 575, 608, 667 |
| 10^5 | 82924 | 1 | 17077 | 49590 of 66667 | 17764 | 1535, 2047, 2207, 2287, 2303 |

All multiples of 3 are singleton SCCs at every N (Theorem 1.2), and #SCC = #(multiples of 3) + #outsiders + 1 in each row. The finite restriction underestimates because n belongs to the giant SCC of E|[1,N] iff both a forward path n -> 1 and an inverse path 1 -> n stay inside [1,N]. For the five smallest outsiders at N = 10^5 the binding constraint is the forward Collatz peak (118096, 1276936, 190996, 250504, 118096, all > N) while their greedy inverse segments peak at 12280, 2047, 4414, 9148, 4606, all below N. By the Q1/Q2 verification each fixed n <= 10^6 lies in the giant SCC for all sufficiently large N.

**Conjecture C_E (OPEN).** In E the non-multiples of 3 form one strongly connected component and every multiple of 3 is a transient singleton feeding it. The Q1 half is the Collatz conjecture; the Q2 half is new.

## 4. A 3-adic Terras theorem for the greedy inverse strategy (PROVED)

Write m_0 = m and m_{i+1} = (2^{k_{i+1}} m_i - 1)/3 with k_{i+1} = k_min(m_i mod 9); K_i = k_1+...+k_i. Then

    3^i m_i = 2^{K_i} m - B_i,      B_i = sum_{l=1}^{i} 3^{l-1} 2^{K_i - K_l} > 0.        (4.1)

(The recovered draft displayed B_i = sum_{l=0}^{i-1} 3^{i-1-l} 2^{K_i-K_{l+1}}, with the exponent of 3 reversed; script S8.1 checks (4.1) for all m <= 3000, 3 not | m, i <= 6 — 12000 instances — and finds the reversed form wrong in 8679 of them, e.g. the word (0,2) from m = 4 has B_2 = 4 + 3 = 7, 9*1 = 16 - 7, whereas the reversed form gives 13. Only B_i > 0 is used below, so Theorems 4.3–4.4 were never at risk.)

**Theorem 4.1 (word determination).** The letters k_1,...,k_J and the residue m_J mod 3 are functions of m mod 3^{J+1}. The modulus is sharp: m = 1 and m = 1+3^J have different J-th letters. The finer residue m_J mod 9 is a function of m mod 3^{J+2} and not of m mod 3^{J+1} (script S8.4: witnesses m = 1 versus 10, 28, 82 for J = 1, 2, 3, whose m_J mod 9 are 1 versus 4).

*Proof.* Induction on i. Suppose m mod 3^{i+1} fixes k_1..k_i and m_i mod 3. By (4.1) m_i is affine in m with slope 2^{K_i}/3^i, so replacing m by m + 3^{i+1} t changes m_i by 3*2^{K_i} t; hence the three lifts t = 0,1,2 send m_i mod 9 bijectively onto the three residues congruent to m_i mod 3, and m mod 3^{i+2} determines m_i mod 9, hence k_{i+1} and m_{i+1} mod 3. The base i = 0 is trivial. Sharpness: by Theorem 2.3 the word of 1+3^J is 2^{J-1}0 while the word of 1 is 2^J. QED

**Theorem 4.2 (exact chain law).** Under the uniform measure on the units mod 3^{J+2}, the residues r_i = m_i mod 9 (0 <= i <= J) form a Markov chain: r_0 is uniform on {1,2,4,5,7,8}, and given r_i, r_{i+1} is uniform on {c, c+3, c+6} with c = 1 if r_i in A = {1,2,4,5} and c = 2 if r_i in B = {7,8}. Consequently the class sequence (A/B) is i.i.d. with P(A) = 2/3, the residue law after every step i >= 1 is 2/9 on {1,4,7} and 1/9 on {2,5,8}, and E_pi[k] = 2(2/9)+1(1/9)+3(1/9)+1(1/9) = 1 < log_2 3.

*Proof.* The uniform measure on units mod 3^{J+2} is a uniform choice of m mod 9 followed by independent uniform choices of the base-3 digits of order 2,...,J+1. By the proof of Theorem 4.1 the digit of order i+1 bijects with the position of m_i mod 9 inside its coset. The class law follows from the table: {1,4,7} has classes A,A,B and {2,5,8} has A,A,B. QED (Script S4 checks the law at steps 1..4 on units mod 3^6; S8.4 rechecks it for J = 1,2,3 on units mod 3^{J+2}.)

**Theorem 4.3 (mean ratio).** For every J >= 1, E[2^{K_J}] = 3 (7/3)^{J-1}.

*Proof.* Given class A the next residue is uniform on {1,4,7} with (k, class) = (2,A),(0,A),(0,B); given B it is uniform on {2,5,8} with (1,A),(3,A),(1,B). The tilted matrix M_{c,c'} = E[2^k 1{next class = c'} | c] is [[5/3, 1/3],[10/3, 2/3]] with det M = 0 and tr M = 7/3, so M^{J-1} = (7/3)^{J-2} M for J >= 2. The initial vector u = (E[2^{k(r_0)} 1{r_0 in A}], E[... in B]) = ((4+2+1+8)/6, (1+2)/6) = (5/2, 1/2) satisfies u M 1 = 7. Hence E[2^{K_J}] = u M^{J-1} 1 = 7 (7/3)^{J-2} = 3 (7/3)^{J-1}; for J = 1 directly (4+2+1+8+1+2)/6 = 3. QED

**Theorem 4.4 (density-1 greedy stopping time).** Let sigma(m) = min{i >= 1 : m_i < m} (infinity if none). For all J, X >= 1,

    #{m <= X : 3 not | m, sigma(m) > J} <= (7/9)^{J-1} (2X/3 + 2*3^J).

Hence the set of non-multiples of 3 with sigma(m) = infinity has natural density 0: the greedy strategy descends on a set of density 1.

*Proof.* If 2^{K_J(m)} < 3^J then m_J < m by (4.1) since B_J > 0, so sigma(m) <= J. Thus {sigma > J} is contained in {2^{K_J} >= 3^J}, which by Theorem 4.1 is a union of residue classes mod 3^{J+1}. By Markov's inequality and Theorem 4.3 the number of such unit classes is at most sum_r 2^{K_J(r)}/3^J = 2*3^J * 3(7/3)^{J-1}/3^J = 2*3^J (7/9)^{J-1}, and each class has at most X/3^{J+1}+1 members in [1,X]. Letting J -> infinity gives upper density 0 for {sigma = infinity}. QED

Exact enumeration (script S4) gives, with f_J = fraction of units mod 3^{J+1} having some prefix i <= J with 2^{K_i} < 3^i and g_J = fraction with 2^{K_J} < 3^J:

| J | units | f_J | g_J | E[2^{K_J}] |
|---|---|---|---|---|
| 1 | 6 | 2/3 | 2/3 | 3 |
| 2 | 18 | 8/9 | 5/6 | 7 |
| 3 | 54 | 25/27 | 43/54 | 49/3 |
| 4 | 162 | 26/27 | 73/81 | 343/9 |
| 5 | 486 | 236/243 | 214/243 | 2401/27 |
| 6 | 1458 | 239/243 | 686/729 | 16807/81 |
| 7 | 4374 | 241/243 | 2126/2187 | 117649/243 |
| 8 | 13122 | 2173/2187 | 12647/13122 | 823543/729 |
| 9 | 39366 | 19609/19683 | 19339/19683 | 5764801/2187 |
| 10 | 118098 | 58868/59049 | 6413/6561 | 40353607/6561 |

f_J is monotone and 1-f_J <= (7/9)^{J-1} holds with much room. All word counts for J <= 4 equal 2*3^J P_chain(word). For J = 6 and m <= 10^5 the residue prediction "descends within 6 steps" agrees with the actual behaviour with zero mismatches. The bound is loose by a factor of about forty (script S8.3): #{1 <= m <= 10^6, 3 not | m, sigma(m) > 10} = 2044 (this includes the greedy fixed point m = 1, for which sigma(1) = infinity; 2043 starts m >= 2), against the Markov bound (7/9)^9 (2X/3 + 2*3^10) whose integer part is 81740, 40.0 times the truth.

**Typed analogy (Terras 1976; CITED: R. Terras, "A stopping time problem on the positive integers", Acta Arith. 30 (1976), pages 241–252 from recollection — verify the page range before external use).** Source: Collatz forward map, parity vector determined by n mod 2^J, stopping time finite on a density-1 set. Target: greedy inverse map on non-multiples of 3, k-word determined by m mod 3^{J+1}, Theorem 4.4. Map: reverse arrows, 2 <-> 3. Preserved: residue determination of the word, tail bound, density-1 stopping time, and the fact that neither implies the full conjecture (the blueprint audit's quantifier gap, cited above). Lost: the criticality class. The Collatz forward step has E[ratio] = (3/2+1/2)/2 = 1 (mean-critical), so Markov's inequality gives nothing and Terras needs a binomial tail; the 3-adic greedy inverse is mean-subcritical with ratio 7/9 per step after the first, so a one-line Markov inequality suffices. Sidecar restoring the difference: the tilt eigenvalue rho (7/3 here versus 3 for the mean-critical case). Decisive test: section 7 shows rho_p < p only for p = 3.

## 5. Cycle census of E: bounded (FINITE-EXACT) and global for every length <= 26 (FINITE-EXACT)

**Theorem 5.1.** Every E-cycle with a arrows of type 3n+1 and h halvings satisfies 2^h > 3^a (E_-: 2^h < 3^a); every E-cycle other than (1,4,2) uses at least one even->3n+1 arrow, provided its minimum is <= 10^6.

*Proof.* Going once around, 2^h n_0 = 3^a n_0 + B with B a positive combination of powers of 2 and 3 (sign reversed for 3n-1). A cycle avoiding even->3n+1 arrows is a cycle of C; by section 3 every n in [2, 10^6] descends under C, so no C-cycle has minimum in [2, 10^6]. QED (The script's S5 remark about the published verification bound 2^68 is UNCITED-RECOLLECTION — Barina, J. Supercomputing 2021, from memory — and is not used anywhere.)

**Bounded census (script S5).** The canonical DFS from the minimum node (iterative, length cap 40) finds exactly **74** simple cycles on nodes <= 2000, none through 3Z. Length histogram: 3:1, 8:1, 13:6, 16:1, 21:2, 26:11, 34:9, 39:43. The (a,h) pairs are exactly (1,2),(3,5),(5,8),(6,10),(8,13),(10,16),(13,21),(15,24); in every case h = ceil(a log_2 3), i.e. 3^a < 2^h < 2*3^a. The full list of the 74 cycles (length, a, h, e = number of even->3n+1 arrows, cycle from its minimum) is printed in the .out; the shortest ones are

```
 3 | 1| 2|0| [1, 4, 2]
 8 | 3| 5|1| [4, 13, 40, 20, 10, 5, 16, 8]
13 | 5| 8|2| [13, 40, 20, 61, 184, 92, 277, 832, 416, 208, 104, 52, 26]
13 | 5| 8|1| [14, 43, 130, 65, 196, 98, 49, 148, 74, 37, 112, 56, 28]
13 | 5| 8|2| [16, 49, 148, 74, 37, 112, 56, 28, 85, 256, 128, 64, 32]
13 | 5| 8|1| [19, 58, 29, 88, 44, 22, 67, 202, 101, 304, 152, 76, 38]
13 | 5| 8|1| [19, 58, 29, 88, 44, 133, 400, 200, 100, 50, 25, 76, 38]
13 | 5| 8|1| [20, 61, 184, 92, 46, 23, 70, 35, 106, 53, 160, 80, 40]
16 | 6|10|1| [2, 7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4]
21 | 8|13|1| [5, 16, 8, 25, 76, 38, 19, 58, 29, 88, 44, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10]
21 | 8|13|1| [5, 16, 49, 148, 74, 37, 112, 56, 28, 14, 7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10]
```

**Both caps of "74" bind (script S8.5).** Raising the length cap to 45 on nodes <= 2000 gives 100 cycles, the 26 new ones all of length 44 with (a,h) = (17,27); raising the node cap to 2500 at length <= 40 gives 104 cycles (30 more than 74); nodes <= 4000 and length <= 45 give 410 cycles. In all 410, h = ceil(a log2 3).

**Global completion (FINITE-EXACT, script S8.5).** A simple E-cycle rotated to start with a 3n+1 arrow is a composition (h_1,...,h_a) of h into a parts >= 0, and going around once gives n_0 (2^h - 3^a) = B with B = sum_{i=1}^{a} 3^{a-i} 2^{h_1+...+h_{i-1}} determined by the word; so n_0 = B/(2^h-3^a) must be a positive integer and the walk must halve only even nodes. Enumerating every such word of every length L <= 26 (6720473 compositions, no node bound) gives:

| L | #E-cycles in all of E | (a,h) | largest node | with a node > 2000 |
|---|---|---|---|---|
| 3 | 1 | (1,2) | 4 | 0 |
| 8 | 1 | (3,5) | 40 | 0 |
| 13 | 6 | (5,8) | 832 | 0 |
| 16 | 1 | (6,10) | 52 | 0 |
| 21 | 2 | (8,13) | 148 | 0 |
| 26 | 22 | (10,16) | 7168 | 11 |

and no E-cycle of any other length <= 26. Hence: every E-cycle of length <= 25 has all its nodes <= 2000, so the bounded census is complete there; there are exactly 1, 1, 6, 1, 2 E-cycles of lengths 3, 8, 13, 16, 21 in all of E; the user's cycle 2->7->22->11->34->17->52->26->13->40->20->10->5->16->8->4 is the unique E-cycle of length 16 globally (a finite exhaustive check, not a bounded observation); at length 26 the bounded census sees 11 of 22 (the 11 others have a node > 2000, largest 7168, e.g. [5, 16, 49, 148, 445, 1336, 668, 334, 167, 502, 1507, 4522, 2261, 6784, 3392, 1696, 848, 424, 212, 106, 53, 160, 80, 40, 20, 10]); and h = ceil(a log2 3) holds for every E-cycle of length <= 26. Whether h = ceil(a log2 3) is forced for all E-cycles remains OPEN beyond length 26.

## 6. The signed side E_- (3n-1)

**Theorem 6.1 (mirror leaf identity, PROVED).** In E_- the arrows absent from the deterministic 3n-1 map are 2j -> 6j-1 (targets exactly 5 mod 6). For v = 6j-1 the least odd T_- predecessor is n0 = (4v+1)/3 = 8j-1 = 7 mod 8, with 3n0-1 = 4v and (v+1)/3 = 2j = R_-^{-1}(n0) for R_-(n) = 4n-1. Targets 1 and 3 mod 6 get nothing new; 3Z is again a transient forest. (Script S6 checks v <= 10^5.)

**Theorem 6.2 (conjugate drift, PROVED).** The inverse move in E_- is m -> (2^k m + 1)/3 with 2^k m in {2,5} mod 9, and the minimal k satisfies k_-(r) = k_+(9-r): (r, k_-) = (1,1),(2,0),(4,3),(5,0),(7,1),(8,2). Negation m -> -m conjugates the two tables, so Theorems 4.1–4.4 hold verbatim for E_- (hostile family m = 3^j-1, e.g. 242 -> 323 -> 431 -> 575 -> 767 -> 256).

**FINITE-EXACT results (script S6, S8.6).** Q2_- holds for all m <= 10^6 with 3 not | m (E-node peak 150994940 = 9*2^24-4 at m = 797161; max 82 compound moves at m = 919795), but the greedy map G_-(m) = (2^{k_-} m+1)/3 has the 2-cycle {4, 11} (the reversal of the E_- cycle 4 -> 11 -> 32 -> 16 -> 8 -> 4, which uses the new arrow 4 -> 11); m = 4 is the only start <= 10^6 that needs a non-greedy k, via 4 -> 11 -> 59 -> 20 -> 7 -> 5 -> 2 with k-word (3, 4, 0, 0, 1, 0) (the recovered draft wrote (3, 4, 0, 1, 1, 0); 20 -> 7 is (20+1)/3 with k = 0). **REFUTED:** "greedy suffices" for E_- (witness m = 4); the plus-side greedy map G has no cycle with minimum in [2, 10^6]. Q1_- holds for all n <= 10^6: the deterministic 3n-1 path descends below n except for the cycle minima n = 5 and n = 17, each rescued by an even->3n-1 arrow (DFS expanding 13 and 47 nodes, peaks 383 and 3407). The three known 3n-1 cycles (1,2), (5,14,7,20,10), (17,50,25,74,37,110,55,164,82,41,122,61,182,91,272,136,68,34) — the odd skeletons (1), (5,7), (17,25,37,55,41,61,91) of the inherited signed-cycles note — lie in one SCC of E_- together with 1, by the BFS-shortest (hence simple) paths

    1 -> 2 -> 5   (2 arrows);   5 -> 14 -> 7 -> 20 -> 10 -> 29 -> 86 -> 43 -> 128 -> 64 -> 32 -> 16 -> 8 -> 4 -> 2 -> 1   (15 arrows);
    1 -> 2 -> 5 -> 14 -> 41 -> 122 -> 61 -> 182 -> 91 -> 272 -> 136 -> 68 -> 34 -> 17   (13 arrows);
    17 -> 50 -> 25 -> 74 -> 37 -> 110 -> 55 -> 164 -> 82 -> 41 -> 122 -> 61 -> 182 -> 91 -> 272 -> 136 -> 68 -> 203 -> 608 -> 304 -> 911 -> 2732 -> 1366 -> 683 -> 2048 -> 1024 -> 512 -> 256 -> 128 -> 64 -> 32 -> 16 -> 8 -> 4 -> 2 -> 1   (35 arrows).

(The recovered draft displayed a 23-arrow walk 1 -> ... -> 17 that repeats 11, 32, 16, 8 and a 53-arrow 17 -> 1 path; both are superseded by the BFS paths above.) Tarjan on E_-|[1,10^5]: 82979 SCCs, one nontrivial of size 17022, all multiples of 3 singletons, 49645 of 66667 non-multiples of 3 outside the giant (17825 of them <= N/2; smallest 1889, 2305, 2401, 2551, 2561). The E_- census on nodes <= 2000, length <= 40, has 70 cycles (lengths 2:1, 5:2, 15:2, 18:12, 23:4, 28:1, 31:2, 36:46), each with 2^h < 3^a, and its deterministic members are exactly the three known 3n-1 cycles. Q1_- and Q2_- are OPEN in general: Q1_- reduces to "no divergent 3n-1 orbit and every 3n-1 cycle minimum has an E_- path to 1".

## 7. Hostile control: the same chain for pn+1 (PROVED, with exact certificates)

Let p be a prime and o = ord_p(2); write 2^o = 1 + pe with e != 0 mod p (this excludes only Wieferich-type primes; it holds for every prime tested). For s in the subgroup <2> of (Z/p)^* let k0(s) in [0, o-1] be the least k with 2^k s = 1 mod p. For a unit r mod p^2 with s = r mod p in <2> put a(r) = (2^{k0(s)} r - 1)/p mod p.

**Theorem 7.1.** The minimal admissible k for the inverse move m -> (2^k m-1)/p (integral, result not 0 mod p) is k(r) = k0(s) + o[a(r) = 0], and the next class is c(r) = a(r) if a(r) != 0, else e. Over each coset s+pZ the value a runs through Z/p exactly once, so the class sequence of the greedy chain is i.i.d. with law P(c) = (1+[c=e])/p on (Z/p)^*. If 2 is a primitive root (o = p-1) every class is in <2>, the tilted matrix M_{s,s'} = (2^{k0(s)}/p)(1+2^{p-1}[s'=e]) has rank one with eigenvalue rho_p = (1/p)(sum_s 2^{k0(s)} + 2^{p-1+k0(e)}), and E_pi[k] = (p-1)/2 + k0(e)/p. If 2 is not a primitive root, a class outside <2> has no inverse move at all, and the per-step probability of this dead exit is exactly (p-1-|<2>|+[e not in <2>])/p.

*Proof.* 2^{k0+ol} r = 2^{k0} r (1+pe)^l = 1 + p(a + e l) mod p^2, so the result (2^k r-1)/p is a+el mod p; it is nonzero for l = 0 iff a != 0, and for a = 0 the first admissible l is 1 with result e. As t varies, a = a(s) + 2^{k0} t covers Z/p once, which gives the class law on every coset, including the cosets of non-primitive p. The stationary law of the i.i.d. class chain is its one-step law; for primitive 2, k0 is a bijection of (Z/p)^* onto [0,p-2], so E_pi[k] = (1/p) sum_s k0(s) + k0(e)/p + (p-1)/p = (p-2)(p-1)/(2p) + (p-1)/p + k0(e)/p = (p-1)/2 + k0(e)/p, and sum_s 2^{k0(s)} = 2^{p-1}-1; the rank-one matrix has eigenvalue v.(2^{k0(s)})_s = rho_p with v_{s'} = (1+2^{p-1}[s'=e])/p. The dead-exit probability is the mass of the classes outside <2>: p-1-|<2>| classes of mass 1/p, plus one more unit of mass if e itself lies outside <2>. QED

**Corollary 7.2 (general criticality dichotomy, PROVED).** For every prime p >= 5 with 2 a primitive root and e != 0: rho_p = (2^{p-1}-1+2^{p-1+k0(e)})/p >= (2^p-1)/p > p (since 2^p-1 > p^2 for p >= 5), so the chain is mean-supercritical; and for p >= 7, E_pi[k] >= (p-1)/2 >= log_2 p (since 2^{p-1} >= p^2), so it is log-supercritical. Hence among all such primes, mean-subcriticality rho_p < p holds only for p = 3 and log-subcriticality E_pi[k] < log_2 p only for p in {3,5}. (Script S8.7 verifies the formulas residue by residue mod p^2 for p in {3,5,7,11,13,17,19,23,29,37}.)

| p | ord_p 2 | e | E_pi[k] | < log_2 p | rho_p | rho_p < p | dead exit |
|---|---|---|---|---|---|---|---|
| 3 | 2 | 1 | 1 | yes | 7/3 | yes | - |
| 5 | 4 | 3 | 11/5 | yes (2^11 < 5^5) | 47/5 | no | - |
| 7 | 3 | 1 | - | - | - | - | 3/7 (<2> = {1,2,4}) |
| 11 | 10 | 5 | 61/11 | no | 66559/11 | no | - |
| 13 | 12 | 3 | 86/13 | no | 1052671/13 | no | - |
| 17 | 8 | 15 | - | - | - | - | 8/17 (|<2>| = 8) |
| 19 | 18 | 3 | 176/19 | no | 8650751/19 | no | - |
| 23 | 11 | 20 | - | - | - | - | 12/23 (|<2>| = 11, e not in <2>) |
| 29 | 28 | 1 | 14 | no | 536870911/29 | no | - |
| 37 | 36 | 1 | 18 | no | 137438953471/37 | no | - |

For p = 7 the nodes reachable from 1 in E_7 within [1, 10^4] have residues {1,2,4} mod 7 only, so the Q2 analogue is REFUTED outright for p = 7 (witness m = 3). Thus the Markov-inequality proof of Theorem 4.4 is special to p = 3. For p = 5 a rational Chernoff certificate still proves density-1 greedy descent: with lambda = 21/20, rho(lambda) = (1+lambda+lambda^2+lambda^3+lambda^5)/5 = 17876501/16000000 satisfies rho^10 = 3.031275 < lambda^23 = 3.071524 and 2^23 < 5^10, so P(2^{K_J} >= 5^J) <= P(K_J >= 2.3J) <= C (rho^10/lambda^23)^{J/10} (rate 0.99868 per step versus 7/9 for p = 3); the rank-one identity E[lambda^{K_J}] = (1/4) 1^T M^J 1 is checked exactly for J <= 3 and lambda in {2, 21/20}. Script S7 gives the E_5 greedy table directly: k_min by residue mod 25 has mean 23/10 = 2.300 < log2 5 = 2.322, E[2^{K_J}]/5^J = 3, 141/25, 6627/625, 311469/15625, 14639043/390625 and prefix-descent fractions f_J = 3/5, 17/25, 91/125, 489/625, 2514/3125 for J = 1..5. Meanwhile the E_5 greedy map has the cycles {13,83,33} and {17,27,43}, the reversals of the deterministic 5n+1 cycles through 13 and 17 (for m <= 10^5, 5 not | m: 79997 descend, 2 cycle, 0 cap), and a forward search from 13 in E_5 does not reach below 13 within 200001 nodes while 33 and 83 do (7 and 4 nodes); for n <= 2000 with a 20000-node budget the unresolved forward starts begin 7, 9, 13, 17, 21, 23, 29, 31, 37, 39, 41 (full list in the .out). Q1_5 is unresolved, not refuted. The density theorem therefore does not distinguish 3n+1 from 5n+1; the mean-ratio eigenvalue rho_p does.

## 8. Typed analogies and non-connections

- **Doubling forest <-> halving arrows.** Source: inherited strict summand shadow with complement D = {x -> 2x}. Target: E. Map: reverse D and add the summand arrows n -> n+(2n+1). Preserved: every E-halving is a reversed doubling; every 3n+1 arrow, including 1 -> 4 = 1 + 3, has the distinct parents (n, 2n+1). Lost: which of the two summand parents is "chosen" (E keeps only the affine one). Corrected decisive test (script S8.8): the inherited diagonal parent pair (1,1) -> 2 is the shortcut arrow 1 -> 2, which is not an E-arrow; the recovered draft's claim that (1,4,2) "uses the diagonal" was false.
- **Inverse fibre (B1) <-> new arrows.** Theorem 1.1: the map is index shift j -> j-1; preserved: 3 n0 + 1 = 4v and the target; lost beyond j = -1: integrality.
- **Terras (2-adic) <-> Theorem 4.4 (3-adic).** See section 4; the sidecar is the tilt eigenvalue.
- **Finite-word completion (braids2) <-> word determination (Theorem 4.1).** Source: every finite halving word of the deterministic inverse occurs in every basin (inverse_completion note). Target: the greedy k-word is a function of m mod 3^{J+1}. Map: the same affine-slope computation on the other prime. Preserved: a finite residue word cannot decide the global question. Lost: the deterministic-inverse completion is unconditional, the greedy word is one strategy among the admissible k-classes.
- **SCOPE (no map found).** THM-4139/4146 (x^2-29/16, mod-63 lift) and the user's {63, -7/4, -29/16} systems share with E only the word "cycle": E-cycles are integer solutions of 2^h n_0 = 3^a n_0 + B, the quadratic three-cycle is a rational preperiodic point of a degree-two map. THM-3341/THM-3333 (Gaussian squaring of triples, Pell hypotenuses) and THM-1745 (h-spectrum holes {7,21}) were read for the session-wide inheritance list; no predicate shared with E was found. These are stated non-connections, not theorems (the recovered draft labelled the first one PROVED; that label is withdrawn).

## 9. Reproduction

    cd /tmp/math-wt-collatz-mod6-b && python3 04-computation/experiments/collatz_mod6_20260917_extended_collatz_scc.py > 05-knowledge/results/collatz_mod6_20260917_extended_collatz_scc.out

Runs in about 11 s with a peak of about 373 MB; the `python3 -O` replay is identical modulo timing lines (all checks use explicit `raise`). Source sha256 85379a7e840da021fffd5a80edb0f65ac695a9232a489327dc6cc454da18ce66 (printed in the .out); output sha256 0d1e56342410312f92f93a2a3c9d42f66fad5e70ef08f05dcfc2af55f643f77b. Universes: v, n <= 10^5 (S1, S6 leaf identities), m, n <= 10^6 (S3, S6, S8.3), N in {10^3, 10^4, 10^5} (Tarjan), residues mod 3^{J+1} for J <= 10 (S4) and mod 3^{J+2} for J <= 3 (S8.4), nodes <= 2000 and length <= 40 (S5) plus the cap variants 2000/45, 2500/40, 4000/45 and the global word census for all lengths <= 26 (S8.5), p in {3,5,7,11,13,17,19,23,29,37} (S8.7), E_5 to m <= 10^5 and n <= 2000 (S7). Hostile controls: 3^j+1, 3^j-1, m = 4 in E_-, E_5, E_7. History: the script was recovered from the 2026-09-17 agent transcript on 2026-09-21; the recovered draft passed all its own checks and its output matched the frozen .out modulo timing, but it predated the agent's pn+1 section, so section S8 (audit-driven additions) was written and the .out regenerated; the recovered draft note's errors listed above (carry formula, rescue k-word, peak semantics, the 1 -> 17 walk, 31 -> 121, the (1,4,2) diagonal claim, the PROVED label on a non-connection) were corrected against the two audit reports.

## 10. Stopping boundary / next question

The greedy inverse strategy is now understood exactly at the level Terras understood the forward map: residue-determined words, an explicit Markov chain, and a density-1 stopping theorem with an elementary tail (loose by a factor of about forty at J = 10). What is not available is any control of the composition beyond the residue word, i.e. the 3-adic version of the missing global coordinate in the inherited notes; the blueprint audit's quantifier gap applies verbatim. The next question is whether the free choice of k (two classes mod 6 at every step) can be used to build a non-greedy strategy with a provable Lyapunov function; the E_- witness m = 4 shows that at least one non-greedy step is sometimes necessary, and section 7 shows that any such argument must use more than the drift table, since the drift alone is favourable for p = 5 where extra cycles exist. Two sharper finite questions are left open by the census: whether h = ceil(a log2 3) is forced for every E-cycle (true for all lengths <= 26 and all 410 cycles on nodes <= 4000 of length <= 45), and the true exponential rate of 1-f_J (the Markov rate 7/9 is far from the observed decay).
