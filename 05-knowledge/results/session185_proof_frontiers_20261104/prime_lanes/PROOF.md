# The unique 7p channel in the completed seven-state return tree

2026-11-04. **PROVED all-clock strong family and exact valuation distribution; FINITE-EXACT prime-cofactor witnesses.** This extracts an arithmetic consequence of Session184's completed return tree. It is not a new return, not a claim of a new globally attained count, and not a proof of infinitely many prime cofactors or either complete prime lane.

## 1. Inherited interfaces and a sharply scoped increment

The construction and seven-state classification are inherited from Session184 `spectrum/PROOF.md`, Sections 1–4, with the same tree for all suffix SCCs of order at most four by `spectrum_four/PROOF.md`. We use the actual named marked component word, not an arbitrary graph with the same Hamilton count. The relevant word at address 422 is

    S_1 S_0 C_1 C_1.

A symbol S_a is a singleton with a ports. C_1 is a directed triangle with one specific marked vertex. Earlier components beat later components. Attach a transitive prefix TT_t, t≥1, and a root r: prefix beats suffix, every marked suffix vertex beats r, and r beats prefix and unmarked suffix vertices. The terminal marked triangle makes the full tournament strongly connected.

The inherited seven-coordinate signature is (h,p,F,G,J,e,s), where p counts ordered two-slot path covers including empty slots, F/G prescribe the first endpoint/second start to lie in the marked set, J prescribes both, and e/s count full Hamilton paths ending/starting in the marked set. The ordered-sum law is

    (hH, pP, FH+pf, GP+hg, JH+hj+eS+Gf, hE, sH).

The root attachment has Hamilton count h−J+2^t F. Applying this exact interface law gives

    (h,p,F,G,J,e,s) = (9,576,189,333,100,3,9),
    H_t = 189·2^t−91 = 7(27·2^t−13),       |V|=t+9.       (1)

Thus (1) is an unbounded family of actual strong tournaments. The prime-cofactor subfamily remains a primality question, not a consequence of strongness.

Older repository families already give infinitely many strong counts divisible by 7, for example Session175's 27·2^r−17. This packet does not claim that existence as new. Its increment is the unique 7-divisibility channel in the completed tree, its all-depth 7-adic clock classification and exact valuation densities, and literal/reproducible witnesses on this precise return lineage.

## 2. Exactly one of the seven clock families can support 7p

The seven normalized addresses have (F,c=h−J):

- empty: (7,−5)
- 2: (21,−13)
- 4: (21,−11)
- 22: (63,−37)
- 42: (63,−31)
- 422: (189,−91)
- 4222: (567,−271)

All slopes F are divisible by 7. Exactly one constant c, that at 422, is divisible by 7. Therefore every count on the 422 ray is divisible by 7, and no count on any other listed ray is divisible by 7. This remains true after any integer clock shift and after adding unmarked singleton prefixes. In the larger SCC≤4 tree, the inherited normalization produced no further representatives, including no unmarked-Q prefix branches.

Consequently, within this completed tree, 422 is the **only** possible 7p lane, and it supplies no prime count. The other six rays can supply prime counts but never 7p. This classifies the arithmetic roles of the inherited tree; it makes no statement about other entry states, one-clock substitutions, or larger SCCs.

## 3. Complete 7-adic clock filtration

Put q_t=27·2^t−13, so H_t=7q_t. Modulo 7, q_t=6(2^t−1). As the order of 2 modulo 7 is 3,

    v_7(H_t)=1  iff  t is not 0 modulo 3.             (2)

For each j≥1 there is a unique residue a_j modulo L_j=3·7^(j−1) satisfying

    7^j divides q_t  iff  t=a_j modulo L_j.          (3)

The residues are effectively determined by a_1=0 and the unique k in {0,...,6} for which

    a_(j+1)=a_j+k L_j,
    27·2^(a_j+k L_j)=13 modulo 7^(j+1).             (4)

Proof: induction gives 2^L_j=1+7^j modulo 7^(j+1). Thus replacing t by t+kL_j changes q_t modulo 7^(j+1) by k·27·2^t·7^j. Its coefficient is a unit modulo 7, so exactly one k lifts each root. Conversely, reducing any higher root gives the previous root. The same induction shows the order of 2 modulo 7^j is exactly L_j: the order starts at 3, and each new modulus multiplies it by 7. This proves uniqueness over a complete period, not merely along a sampled list. ∎

The first lifts are a_1=0 (mod 3), a_2=6 (mod 21), and further exact residues are stored in `receipt.json`. For e≥2,

    v_7(H_t)=e
    iff t=a_(e−1) mod L_(e−1), but t≠a_e mod L_e.   (5)

Every e≥1 therefore occurs on an infinite arithmetic progression of legal clocks. The exact natural density among positive integer clocks is

    density(v_7(H_t)=1)=2/3,
    density(v_7(H_t)=e)=2/7^(e−1),  e≥2.             (6)

These densities sum to 1; they are valuations of counts in this specified family, not statistical claims about all strong tournaments.

## 4. What this proves and does not prove about prime cofactors

If q_t is prime, necessarily t is not 0 modulo 3. The exceptional possibility q_t=7 cannot occur because 27·2^t=20 has no legal solution. Also t cannot be 2 modulo 4: then q_t is divisible by 5, and q_t≥41 excludes q_t=5. Thus a necessary elementary clock sieve is

    t mod 12 in {1,4,5,7,8,11}.                       (7)

It is not sufficient; e.g. t=5 gives 851=23·37. The exact trial-division prime clocks through 35 are

    1, 4, 8, 13, 19, 28, 32, 35.

They give prime cofactors

    41, 419, 6899, 221171, 14155763,
    7247757299, 115964116979, 927712935923.

For example, the t=35 graph has 44 vertices and exactly 6493990551461 Hamilton paths, equal to 7·927712935923. This count follows from the proved exact interface formula, not a direct enumeration of its 44-vertex paths. The checker certifies the displayed primalities by trial division, without probabilistic primality assumptions.

The existence of infinitely many primes of the form q_t is not proved here. Finite witness lists, the density of eligible clocks, and arbitrarily high 7-adic valuations do not imply prime infinitude. They do not establish all 7p values, the ordinary prime lane, or global spectrum completeness. Any use of this ray as a complete constructor must also overcome the already documented finite-ray barriers.

## 5. Reproduction and provenance

Run `python verify.py > normal.out` and `python -O verify.py > optimized.out`; outputs agree. The script is self-contained and uses only the standard library. It recomputes all seven inherited signatures, verifies the modular lane partition, independently counts literal full graphs by directed subset DP at clocks 1,2,3,4,5,8 and checks strongness, computes twelve exact lifted residue classes, checks the valuation equivalences at clocks 1 through 1000, and trial-divides all candidate cofactors through clock35. Explicit exceptions remain active under optimization.

The all-clock claims rest on the physical ordered-sum identity and Sections 2–3, not extrapolation from these finite controls. The S184 sources were read from the existing shared verified packets. A targeted current-main search at ad09a24c8fd25eca9ea320d1b240c04990ada96b located the older 7-divisible family and finite-ray barriers, notably `05-knowledge/results/session175_tournament_spectrum_20261104/synthesis/NOTE.md` and `05-knowledge/results/session151_spectrum_metagraph_20261103/provenance/AUDIT.md`. No comprehensive novelty claim is made.
