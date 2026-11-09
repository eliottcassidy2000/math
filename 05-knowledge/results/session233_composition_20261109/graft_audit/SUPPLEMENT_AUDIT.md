# Independent audit: weighted nonunit loops and direct CRT composition

Verdict: PASS for the constructions and formulas below. This supplements, without replacing, the original stable theorem audit. No arbitrary-target placement or Chan–Pak completeness follows.

## Independently verified physical inputs

Permutation-and-split enumeration recovers the nonunit13.json state (13,832,73,584,46,0,11) and nonunit29.json state (29,1160,133,668,60,0,21). Literal graph closure counts agree with 2A+h-T. These are saturated deletion cores, despite their cores not needing to be strong themselves.

For h13 the 12-bit words have last(X)=1, first(Y)=0. The exact flanked identity state(s1=>X=>C=>s0) minus state(s1=>C=>Y=>s0) is (0,0,0,0,-10,0,0). For h29 the 28-bit words have last(X)=0, first(Y)=1; the analogous s0/s1-flanked difference is (0,0,0,0,-6,0,0). The inherited h45 12-bit identity has flanks0/1 and difference (0,0,0,0,+2,0,0).

Define F13(K,i)=s1=>X13^i=>C13=>Y13^(K-i)=>s0, F29(J,j)=s1=>s0=>X29^j=>C29=>Y29^(J-j)=>s1, and the inherited protected F45(L,k). Their deletion orders are 12K+10, 28J+11, and 12L+11. Every first deletion vertex is selected and dominates later blocks. An adjacent movement retains the flanked finite difference, because the word boundary bits match those flanks. Induction gives pure-T changes -10i,-6j,+2k for every nonnegative length, not only tested lengths.

## Weighted two-core theorem

For a pure-T change -2g in a deletion h-family and the h45 change +2, ordinal composition and reindexing k=L-u produce count increments 2(45g i+h u). When gcd(45g,h)=1, the audited bounded coin-box proof applies with a=45g,b=h,K>=h-1,L>=45g-1. This proof requires no unit loop in the first family.

At fixed K12 the h13 pair therefore supports 13L-2675 consecutive odd values on 12L+167 vertices for every L>=224. At fixed K28 the h29 pair supports 29L-3723 consecutive odd values on 12L+808 vertices for every L>=134. The offsets from the zero/reversed-corner baseline are respectively the coin intervals [2688,225K+13L-2688] and [3752,135K+29L-3752]. Strongness follows from the leading selected singleton and the common v,r closure.

## Direct three-core CRT theorem

Take F13(12,i)=>F29(28,j)=>F45(L,L-k), with 0<=i<=12,0<=j<=28,0<=k<=L, then add the common closure v,r. Let B_L be the closure count at i=j=k=0. Exact composition gives

    H = B_L + 2(6525i + 1755j + 377k).

Indeed half of the first drift is 5*29*45=6525, half of the second is 3*13*45=1755, and the last multiplier is 13*29=377. All non-T interface coordinates are fixed. The h value is 16965, q is divisible by it, and no unsaturated-fiber claim is implied.

Modulo 13, 6525 is invertible and 1755 vanishes. Modulo 29, 1755 is invertible and 6525 vanishes. Thus the 13*29 choices of (i,j) give each residue modulo377 exactly once. Set R(i,j)=6525i+1755j. The smallest representative is0 and largest is127440. For each t in

    [127064, 377L+376],

select the unique (i,j) with R(i,j)=t mod377 and set k=(t-R(i,j))/377. If k were negative, t<=R-377<=127063, a contradiction. If k>L, t>=R+377(L+1)>=377L+377, a contradiction. Therefore this selector is physical for every target in the interval. The interval is nonempty for every integer L>=337 and has 377L-126687 integers, giving that many consecutive odd counts.

The selected predecessor127063 cannot occur: it has the same residue as the unique representative127440 and would need k=-1. The successor377L+377 has the same residue as representative0 and would need k=L+1. Thus the displayed central interval has sharp adjacent boundaries in this fixed rectangle; this does not claim a globally maximal interval over other constructions.

The total graph order is (12*12+10)+(28*28+11)+(12L+11)+2=12L+962. All witnesses are strong by the first-selected-singleton argument. Examples: L337 gives362 consecutive odds on5006 vertices; L400 gives24113 on5762 vertices. The asymptotic width/order ratio is377/12, but the baseline remains exponential and no target-placement argument is supplied.

## Independent executable evidence and scope

The separate audit_supplement.py imports only independently authored elementary helpers from audit_common.py, never the producer verifier. It recovers all three input core states by explicit permutation enumeration, checks their literal star closures, and computes the finite flanked identities exactly. It tests every family position through k100 at selected lengths. For six pair parameter choices it checks every physical rectangle position and every target in the claimed window. It checks all377 CRT residues, and every target and its full interface in the direct triple intervals at L337,338,400. Both immediately adjacent omitted offsets are independently checked.

Run `python3 graft_audit/audit_supplement.py` in an assembled packet with sibling graft and graft_audit directories; otherwise give `--source PATH`. Normal and -O outputs must agree. The separate `--sentinel` intentionally raises RuntimeError before loading data, under -O as well. Receipts are SUPPLEMENT_CHECK.json, SUPPLEMENT_CHECK.optimized.json, and SUPPLEMENT_SENTINEL.txt. No literal graph DP at the large witness orders is claimed.

## Final producer text and nested comparison

The stable WEIGHTED_PROOF.md, verify_weighted.py and WEIGHTED_CHECK.json were subsequently reviewed. Their stated concrete constructions and receipt values agree with the independent calculations. The general theorem applies the earlier lemma verbatim when a,b>=2. If a=1, the hypotheses L>=0 and K>=b-1 make [bj,bj+K] adjacent for successive j, giving all integers [0,K+bL]; the b=1 case is symmetric. Thus the general statement is also true in these cases, but this one-sentence extension should be supplied or a,b>=2 explicitly required.

The optional nested construction is sound. For h13/h29 the half-count coefficients are145 and39, with K38,J149. The coin conductor is5472, the maximum offset is11321, and the closed interval [5472,5849] has378 integers. The deletion order is12*38+10+28*149+11=4649. Independent recomposition checks all378 selected inner states, their common interface, and T differences -2r. Retain r0 through376 and compose with h45. The coefficients are45 and377, the conductor is16544, and the outer maximum is16920+377L. Hence the width is377L-16167 and order12L+4662 for L>=44. Protected first vertices retain strongness. Independent physical endpoint and midpoint controls pass at L44,100,1000. The L44 case gives421 consecutive odd counts on5190 vertices.

The direct CRT L1000 lower endpoint was independently recomputed by repeated word composition, not producer exponentiation: bit length12966, SHA-256 fe78f2c31a532384fc58a3329368bcdf8bcb3cc22d382b7f8ba3bb9f1fc22ac1, as claimed. Five physical target controls at that size also pass. The supplement normal/-O receipts remain identical after these additions.

Recommended editorial qualifications: call the weighted result a sufficient escape from the repeated-core divisor obstruction rather than a converse (necessity was not proved). Replace “does not prove cofinal occupancy” by “does not prove eventual coverage of all admissible targets,” because unbounded attained counts already constitute a cofinal subset in the usual order-theoretic sense. Neither issue affects the explicit positive theorem. Discovery census counts are producer provenance, not independently audited here and not needed for validity. No literature-priority claim is certified.
