# Repunit recovery and a measure-aware rarity theorem

## Deliverables and scope

This packet searches present tracked filenames and locally present UTF-8 content, without Git lazy fetching. `catalogue_initial.json` freezes the first pass at dca0424aed4cf8a7a1b5397add40b3266789d51c: 113,749 tracked paths, 72,868 text files searched, 40,455 absent tracked files, 2 oversized files skipped, 2,640 matching entries and 6,720 matching lines. `catalogue.json` is the reproducible expanded pass after the shared checkout advanced to 54bf9664293, with exact HEAD, hashes, line numbers, and the regex recorded. `discovery_summary.txt` records its totals. This is a discovery index, NOT a claim that every hit concerns Hamiltonian-path counts or that every historical blob was searched. Broad Mersenne hits include Collatz and inverse-polynomial work. Formula regexes cannot capture every typography or unnamed use of 31/63/127. The 700-character match previews may truncate very long lines; source paths are retained.

Reproduce: `python3 discover.py /path/to/math catalogue.json`.

Four targeted GitHub connector reads recovered sources missed or only partially surfaced locally; copies under `fetched/` preserve the returned text. They are the exact named paths below. Connector default-branch reads were on 2026-10-22; they are not silently equated with the initial checkout snapshot. No repository files were changed in this lane.

## Inheritance and concept board

Closest proved mechanism: multiplicativity over the ordered strong components (Session24; THM-4749). Hostile example: 255 is a repunit yet ultimately more common than 257. Corrected near miss: minimum feedback-arc complexity does not explain the observed fixed-order residue depletion (MISTAKE-657). Underused sidecar: automorphism weights distinguish uniform isomorphism classes from uniform labelled tournaments.

Board: (1) ordered SCC factors; (2) one-apex cancellation; (3) binary carry/odd-cycle coefficients; (4) automorphism group action; (5) trace versus Hamiltonian-count typing; (6) Collatz truncation versus tournament construction. The rational generating function links (1) and (4), while the fan links (2) and (3). Shared Mersenne integers alone do not identify (5) or (6) with either.

## Curated provenance routes

Each path is relative to the repository. Status below distinguishes a theorem's valid core from historical glosses.

1. **2026-03-07, Mersenne ORDERs, PROVED.** `01-canon/theorems/THM-067-c1-formula-mersenne.md`: c1=2^(f+1)-d-2 vanishes at tournament order n=2^(f+1)-1 for even f. This produces n=7,31,127,511, not path counts H=7,31,127,511. It must not be used as a rarity theorem.
2. **2026-03-14, gap repunits, HISTORICAL / partly false.** `04-computation/repunit_gap_synthesis.py`, corresponding `05-knowledge/results/repunit_gap_synthesis.out`: 7=111 base2 and21=111 base4 are exact; generalized realizability claims are unsupported and even the proposed coefficient algebra fails (see corrections).
3. **2026-03-16, k-nacci trace, PROVED core.** `01-canon/theorems/THM-227-knacci-mersenne.md`, `04-computation/prove_knacci_mersenne_s90ap.py`: Tr(M_k^r)=2^r-1 for1<=r<=k by Newton identities; first departure Tr(M_k^(k+1))=2^(k+1)-k-2. The current Connections section still says 'forbidden value'; this does not follow and is false as an H assertion.
4. **2026-03-16, trace/Cayley exploration, HISTORICAL.** `04-computation/mersenne_connections_s90aq.py`, matching `.out`, and `04-computation/synthesis_tanh_mersenne_s90au.py`: valid trace and rational identities coexist with invalid spectral AM-GM and forbidden-value commentary. Do not promote the whole exploratory output to theorem status.
5. **2026-05-29, 63 unlocked, FINITE-EXACT.** `01-canon/theorems/THM-344-n8-h63-complete-omega.md`: exactly two classes at n8, both asymmetric, both complete odd-cycle conflict graph K31, cycle lengths {3:8,5:17,7:6}. Reproduction: `04-computation/h63_n8_isoclass_census_s10.py` and `h63_complete_omega_s10.py`. Supersedes HYP-1754 and the 7*3^j permanent-gap tower.
6. **2026-06-10, binary coefficient tower, PROVED.** `01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md`: H=sum alpha_j2^j, alpha_j counts j mutually disjoint directed odd cycles. These are redundant coefficients, not binary digits individually. Carries matter.
7. **2026-07-20, permanent gaps, VERIFIED proof / prior art.** `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`: SCC multiplicativity, insertion monotonicity, Moon strong deletion, finite small spectra exclude7 and21 at all orders. All other odds being attainable remains OPEN globally; original MathOverflow2016 attribution is explicitly restored.
8. **Family and factorization routes.** `05-knowledge/variables/single-core-cycle-count.md` already contains the two H63 apex words; `THM-1355-proth-unification-observer-hypotenuse-two-axes.md`, `THM-1360-figurate-stagger-atlas-pisot-ladder-proth-mersenne-relational.md`, and `THM-1425-the-ocf-is-the-weight-two-shifted-pascal-identity.md` supply earlier Proth, Fermat, and Jacobsthal relations. HYP-2198/2199/1762 and Session16 feedback-word sources appear in the catalogue/provenance synthesis, but a full content recovery of absent files is NOT asserted.
9. **DRT/lattice order tower, a different Mersenne carrier.** `01-canon/theorems/THM-448-drt-mersenne-doubling-tower.md`, THM-447/450/451/455/480 nearby: doubling tournament ORDER and Hadamard/code constructions, not a theorem about low multiplicity of a fixed H. The modular/lattice lane supplies the detailed map.
10. **2026-09-24, digits and Collatz, PROVED core / corrected boundaries.** `THM-4473-collatz-digit-chains-rotation-repunits.md`, `05-knowledge/results/procgen_repunit_20260924_digits_rotation_repunits.md`: explicitly retain the map's convention, k>=2 in multidigit rotation-fixed claims, and the correction that -1/3 is fixed for x->4x+1, NOT the Collatz odd-step map.
11. **2026-10-06 to19, Mersenne Collatz line.** THM-4555/4556, THM-4605, THM-4712/4713, HYP-9213/9217/9243/9245/9332/9333; `mersenne_line_barriers_20261007.md`, `mersenne_horizon_20261019_clusters_orphans_reanchoring.md`. These concern coalescence, debt, clocks, exits and finite horizons, not H-count rarity.
12. **2026-10-21, actual multiplicity law.** `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md` and `h_spectrum_multiplicities_20261021.md`; THM-4744/4745; MISTAKE-652 repairs explanations based merely on primality.
13. **2026-10-21, repunit-specific synthesis.** THM-4748 and THM-4749 (full filenames in catalogue); `05-knowledge/results/repunit_rarity_hyp9357_20261021.md`; experiments `repunits_hyp9357_20261021/` including `fan_uniqueness.py`, `repunit_degrees.py`, `repunit_rarity.py`, `class_residues.py`, `refuted_leads.py`; HYP-9359 remains a scoped empirical hypothesis, not a universal rarity proof.
14. **2026-10-22, one-apex scope and measured boundary.** `05-knowledge/results/session26_random_bridges_20261022/SYNTHESIS.md` and its audited subpackets classify one-apex induced wiring, not all tournament structure.

### Correction lineage that overrides older statements

MISTAKE-024/050:63 is attained. MISTAKE-028:31 and127 are attained. MISTAKE-653:2^29-1 is composite233*1103*2089, withdrawn forbidden-Mersenne gloss, and stale HYP-1154. MISTAKE-652: primality and factor7 alone do not explain finite multiplicity. MISTAKE-657: subtraction is NOT restricted to shared-endpoint back arcs outside the one-apex setting; low-feedback-arc residue depletion is refuted across n6..9; full strong spectrum remains conjectural; 'blind predictions' need the precise partial data exposure history. The finite inventory counts were independently reproduced, but this does not restore those causal explanations.

## A proved measure-aware explanation

Let S be the set of isomorphism classes of nontrivial strongly connected tournaments. Write h(A)=H(A), m(A)=|A| and a(A)=|Aut(A)|. For fixed h, only finitely many A with h(A)<=h occur, by the known exponential lower bound on strong Hamiltonian counts used in Session24. Set

N_h(n)=number of ordinary isomorphism classes of order n with H=h,
W_h(n)=sum_{[T], |T|=n,H=h}1/|Aut(T)|.

**Theorem (weighted extension of Session24).** For n large enough,

W_h(n) = sum_{r>=0} sum_{(A_1,...,A_r): product h(A_i)=h} binom(n-M+r,r)/product a(A_i),

where M=sum m(A_i). For all n use zero for terms with n<M. The same formula without the automorphism denominator gives N_h(n). For attainable h, consequently both have degree w(h), the maximum number of nontrivial strong factors in a factorization of h by actually attained strong values. Their leading coefficients are positive sums over maximizing sequences divided by w(h)!. For unattainable h both counts vanish, and no maximum over an empty factorization set is asserted. For h=1 the empty sequence gives w(1)=0 and N_1(n)=W_1(n)=1.

**Proof.** Every tournament has a unique totally ordered condensation of strong components. Nontrivial components form the displayed ordered sequence. Insert the n-M singleton components in its r+1 gaps in binom(n-M+r,r) ways. Hamiltonian paths concatenate and hence their counts multiply. An automorphism fixes the ordered condensation pointwise, so it acts independently inside each nontrivial component; its order is product a(A_i). All summands are positive. The degree and leading coefficient follow. This uses the actual strong spectrum, not its conjectured completion. QED.

The number of labelled tournaments is L_h(n)=n!W_h(n); uniform-labelled probability is n!W_h(n)/2^(n(n-1)/2). Thus fixed-h comparisons have the SAME polynomial-degree distinction under uniform classes and uniform labelled tournaments, although the coefficients change. Weight H/|Aut| likewise multiplies W_h by h and preserves its degree for fixed h. This is not a claim that the finite-n histograms coincide.

The generating function is a finite sum of z^M/(1-z)^(r+1), with the corresponding positive component weights. This establishes a concrete bridge to the modular lane's rational-GF obstruction without numerical analogy.

### Application to31,63,127

31 and127 are prime and attained, so w=1. The number63 has only proper nontrivial factorizations through7 or21, both impossible, so w(63)=1 despite its compositeness. In contrast33=3*11,65=5*13,129=3*43 each admit two actual strong factors, so w>=2. Therefore both N_31/N_33 and W_31/W_33 tend to0, and similarly63/65 and127/129. This is a proof of relative rarity as ORDER n grows at FIXED requested values, not an explanation of every small-order outlier.

For63, the strong floor at n9 exceeds63 and the complete n<=8 census gives exactly the two asymmetric order8 atoms. Hence exactly
N_63(n)=W_63(n)=2(n-7) for n>=8,
L_63(n)=2(n-7)n!, and P_labelled(H=63)=2(n-7)n!/2^(n(n-1)/2).
At n8 this is80640/2^28. At n<=7 there are no realizations. The structural gap at strong order7 is still a census fact, not explained by the polynomial asymptotic.

### Hostile controls and limits

15=3*5 and255=3*5*17 have degrees2 and3, respectively, whereas17 and257 are prime with degree1. Repunits can be MORE common than their 2^k+1 counterparts. Minimum feedback arc number3 versus1 is a construction-complexity distinction, not a counting theorem. The Fibonacci fan proves actual strong realizations for every2^k-1,k>=4, but it does not prove uniqueness or distributional rarity.

For the different measure S(h)=number of strong classes across all orders at exact H=h, singleton padding is absent. The SCC degree argument says nothing about S(31)/S(33) or any large-h density. Strong classes bounded by H<=B, classes at fixed n, and automorphism-weighted ensembles must be reported separately. The S51 arithmetic-support fit does not substitute for a theorem in that regime.

## Reproducible verification

`python3 verify_weighted.py > verify_weighted.json` independently enumerates every labelled tournament at n1..6 (33,867 objects), computes H by subset DP, tests strong connectivity directly, and verifies the ordered-SCC weighted generating identity using exact Fractions. No prior H census is imported. It also independently verifies the Fibonacci fans by Hamiltonian DP for k4..12. All checks passed. These finite checks support the code/interpretation; the all-order result above follows from the proof, not extrapolation.

## Proposed repairs, not applied

`corrections.patch` repairs THM-227's currently surviving 'forbidden' gloss without changing the true trace theorem. The old exploratory `repunit_gap_synthesis.py` and `.out` should be given a prominent historical-warning route, not treated as a current proof. Its purported alpha=(T-13,3,1) gives1+2alpha1+4alpha2+8alpha3=2T-5, not2T+1; and alpha=(0,3,1) is unrealizable since one independent triple already requires at least three vertices and three pairs in the conflict graph. Thus its claimed global completeness proof fails before any deep tournament question. No new theorem ID or mistake ID is reserved here.

### Final expanded discovery totals

At 54bf96642934c6cdf01f02e364f3bba00892ae22: 113,808 tracked paths; 72,894 present text files searched; 40,488 absent tracked files; 424 nontext/read-error files; 2 oversized skips. The expanded broad-term index has 7,290 matching entries and 12,379 matching lines. Literal term file counts overlap: repunit 249, Mersenne 1237, Wagstaff 30, Fibonacci fan 20, trailing ones 26. The increased total versus the first pass principally reflects broad 'all ones' coverage, not newly discovered relevant theorems. Use the curated provenance routes before interpreting raw hits. `fetch_manifest.json` stores connector blob SHAs and verified URLs.
