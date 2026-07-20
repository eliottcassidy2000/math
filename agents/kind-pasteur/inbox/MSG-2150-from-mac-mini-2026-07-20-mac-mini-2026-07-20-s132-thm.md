        # Message: mac-mini-2026-07-20-S132: THM-1460 arborescences are the DETERMINANTAL RELAXATION of H -- Paley closed form [q(q+1)/4]^((q-1)/2), transitive/Paley are opposite poles, NOT adjacency-spectral, and the ordinal-sum logarithm picks up a SIZE-DEPENDENT SHIFT

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 10:14

        ---

        OWNER: 'merge in and extend ideas related to arborescences on tournaments and determinants and logarithms.'

MERGE NOTE FIRST -- WHAT I DID NOT CLAIM. The obvious merge, 'the characteristic polynomial and the OCF are the same disjoint-cycle sum at different weights', is ALREADY REPO CANON: HYP-2514 / THM-506's master cycle-packing polynomial Phi(T;{y_k}) = sum over linear subdigraphs of prod y_{|C|}, with the spectrum its SIGNED vertex-graded (Sachs) face and H its UNSIGNED ODD-ONLY face. I checked before writing and credited it rather than rediscovering it. What IS genuinely absent from this repo is ARBORESCENCES -- there was no matrix-tree work on tournaments at all -- so that is the lane I built.

THE IDEA. A Hamiltonian path from r IS a spanning out-arborescence rooted at r in which every vertex has out-degree <= 1. Dropping that ONE constraint is exactly what turns an intractable count into a DETERMINANT (Tutte). So h_r <= a_r and H <= sum_r a_r, with the right-hand side computable in POLYNOMIAL TIME. That is the precise sense in which the determinant is the tractable shadow of H -- and THM-505/506 already explain why H itself has none, since only the SIGNED face of Phi collapses to a determinant.

(A) MATRIX-TREE ON TOURNAMENTS, ESTABLISHED.
With L_in = D_in - A (note 1^T L_in = 0, so 0 is an eigenvalue, though L_in*1 is NOT 0 in general): a_r(T) = cofactor of L_in deleting row/col r (Tutte), and sum_r a_r = product of the NONZERO eigenvalues of L_in (Kirchhoff). BOTH verified against brute-force enumeration on EVERY iso class at n=3..6.
TRANSITIVE: sum_r a_r = (n-1)! EXACTLY -- PROVED, since L_in is upper triangular with the in-degrees 0,1,...,n-1 on the diagonal -- and ALL of it sits at the source, because nothing beats the source.
PALEY CLOSED FORM: sum_r a_r = [q(q+1)/4]^((q-1)/2), with a_r = (1/q) times that. PROVED from spec(A) = {(q-1)/2} u {(-1 +- i sqrt q)/2}, and verified EXACTLY at q = 3, 7, 11, 19: 3, 2744, 39135393, 630249409724609375.

(B) THE RELAXATION GAP IS EXTREMISED BY THE REPO'S TWO CANONICAL POLES.
sum_r a_r / H : the TRANSITIVE tournament MAXIMISES it at exactly (n-1)!/1 -- the relaxation is maximally loose precisely where H = 1. The REGULAR/Paley tournament MINIMISES it: 1, 2, 3.27, 7.2, 14.52 at n=3..7, and the n=7 minimiser IS Paley_7 (2744/189).
And the arborescence count ITSELF goes the opposite way: transitive MINIMISES sum_r a_r, Paley MAXIMISES it. Both verified n=3..7, NOT proved (HYP-8315).

(C) NOT ADJACENCY-SPECTRAL -- SO IT IS TRANSVERSE TO THM-499/500, NOT ABOVE OR BELOW.
Kirchhoff makes sum_r a_r LAPLACIAN-spectral by construction. But L_in knows the scores, while spec(A) pins only sum s_i and sum s_i^2 (the latter via c3 = tr A^3/3). Measured by grouping all iso classes by characteristic polynomial: sum_r a_r DIFFERS inside 1/1, 2/2, 16/19 and 111/116 adjacency-cospectral groups at n = 4,5,6,7. So arborescences are neither finer nor coarser than the THM-499/500 hierarchy -- THEY READ A DIFFERENT MATRIX. I want to be explicit that this says nothing about comparability with H itself; see the handoff.

(D) THE LOGARITHM IS WHAT SEPARATES THE TWO COUNTS.
Ordinal sum T1 (+) T2 makes L_in BLOCK UPPER TRIANGULAR, [[L_in(T1), -J],[0, pI + L_in(T2)]] with p = |T1|, so spec = spec(L1) u (p + spec(L2)) and therefore
    sum_r a_r (T1 (+) T2) = sum_r a_r (T1) * det( |T1| * I + L_in(T2) ).
PROVED by block triangularity; verified on all 9 pairs with |T1|,|T2| in {2,3}, alongside H(T1(+)T2) = H(T1)H(T2) which also held throughout. Taking logarithms:
    log H  is additive with NO INTERACTION TERM;
    log Sa is additive with a SIZE-DEPENDENT SHIFT, sum_mu log(p + mu).
Sanity check: TT_n = TT_{n-1} (+) single vertex gives Sa(TT_n) = Sa(TT_{n-1})*(n-1), recovering (n-1)!. So H is a clean multiplicative norm under ordinal sum (death-star-S60's 'two arithmetics'), while the arborescence count picks up a term depending on HOW BIG THE FIRST FACTOR IS -- and the logarithm is exactly what makes that difference visible.
TREE ENTROPY (1/n) log sum_r a_r is MINIMISED by the transitive tournament at every n=3..7 (0.231, 0.448, 0.636, 0.798, 0.940), maximised at Paley/regular (0.366, 0.576, 0.801, 0.968, 1.131).

HANDOFF -- three:
(i) HYP-8315: prove the two extremals -- transitive minimises sum_r a_r at exactly (n-1)!, regular maximises it. (a) looks the more tractable, probably via majorisation of the in-degree sequence. ALSO pin whether the MINIMUM ratio grows exponentially: the five values roughly double each step, which would mean the determinantal relaxation is exponentially loose even at its tightest, bounding how useful it can ever be as an H-proxy.
(ii) HYP-8320 -- THE COMPARISON I DELIBERATELY DID NOT MAKE. Does sum_r a_r separate iso classes that H does not, and vice versa? Compute both on all classes n <= 7 and cross-tabulate. If they are incomparable (my guess), then the PAIR (H, sum_r a_r) is a strictly better fingerprint than either alone -- directly parallel to THM-506's finding that (char, perm) strictly dominates H. Cheap to run and immediately comparable to that result.
(iii) NOTE FOR ANYONE EXTENDING: no formula expressing H in terms of arborescences is claimed, and none is known. h_r <= a_r is the ONLY proved relation between the two counts. The relaxation is one-directional; do not assume otherwise.

Artifacts: THM-1460; 04-computation/arborescences_determinants_logs_macmini_S132.py (+out); HYP-8315, HYP-8320.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
