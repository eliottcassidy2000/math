# A response cone and constructive low-cycle arithmetic inverse

Session 197. Status: all-order response inequalities and exact target inversion; finite-exact coverage. No cofinite or all-order scalar-coverage theorem is claimed.

## 1. Inherited physical grammar

Use the definitions and proved one-seam law in `ROOTED.md`, its complete integer gadget compiler `unique_start.py`, and the complete C=1 classification in `C1_COMPLETE.md`. For a rooted input write x=(h,e,R), where h counts all Hamilton paths, e counts paths ending at the root, and R counts ordered covers (P,Q) with P possibly empty and Q starting at the root. Always R>=h>=1 and 0<=e<=h.

The allowed gadget has transitive root-deletion and exactly one path ending at its root. It is indexed by L>=1, with constants g,Q,B and t=L mod2:

    h'=g h-t e+L R,  e'=h-[L=1]e,  R'=Qh+BR.       (1)

Here g(1)=Q(1)=B(1)=1. For ell>=1, even L=2ell has (g,Q,B)=(g_ell,Q_ell,2B_ell). For odd L=2ell-1>=3 these are (g_ell+ell-(ell mod2),Q_ell+B_ell,2B_ell). Every such integer L is an actual gadget, not a formal transition.

The closest inherited mechanism is the complete bounded rooted image in `ROOTED.md`. The new work below sharpens its response domain and gives an exact congruence-based final-gadget inverse, with actual positive witnesses. It does not change the one-seam law.

## 2. An all-order output cone

**Theorem 1.** Applying any allowed gadget to any nonempty rooted tournament gives

    h' <= R' <= 2(h'-e').                            (2)

In particular this holds for every nontrivial output of the rooted C=1 grammar. It is not asserted for its singleton seed (1,1,1).

First, Q>=g and B>=L. These follow inductively from the displayed recurrences: in the odd step the increment of Q-g is B_ell-ell+(ell mod2)>=0. Also B(L) is the least power of two at least L, so for L>=2, B<=2L-2. Hence R'-h'=(Q-g)h+t e+(B-L)R>=0 for L>=2; for L1 it equals e.

For the upper bound set d=2g-Q. Even steps preserve d. Odd L=2ell-1 gives

    d(L)=d(ell)+2ell-2(ell mod2)-B(ell).

Induction gives d>=1 for every L and d>=3 for odd L>=3. If ell is even the increment is at least2; if ell is odd>=3 it is nonnegative and d(ell)>=3. The base ell2 gives d(3)=3.

For L>=2,

    R'-2(h'-e')=(Q-2g+2)h+2t e+(B-2L)R.

If L is even this is at most h-2R<=0. If L is odd it is at most -h+2e-2R<=0. For L1 it is h-R<=0. This proves (2).

**Corollary 2.** From any nontrivial grammar state, every next output satisfies h'>=3R/2>=3h/2. Indeed L1 gives h'=R+(h-e)>=3R/2; L>=2 gives h'>=LR>=2R. Thus a construction of target H<=N has every nontrivial preterminal response R<=floor(2N/3), and its number of subsequent gadget steps is logarithmic in N. This is a necessary physical domain, not a sufficient realization criterion.

In coordinates (a,b,c)=(e,h-e,R-h), the singleton gadget is the unimodular map

    (a,b,c) -> (b,a+b+c,a).

Its formal inverse is (a,b,c)->(c,a,b-a-c); nonnegativity of that inverse is exactly the upper cone inequality. **It is not a physical inverse-closure theorem.** The actual L3 output from the singleton is (h,e,R)=(5,1,7); its formal singleton inverse is (3,2,4), which violates (2) and is not a nontrivial rooted-grammar state. This is a concrete obstruction to replacing graph inversion by a cone calculation.

## 3. Exact final-gadget arithmetic inversion

Fix a physically reachable prefix (h,e,R) and a target N. Check L1 separately. For L>=2, write L=2k+t, t in{0,1}. Equation(1) forces

    2Rk = N+t e-Rt  (mod h).                         (3)

Let d=gcd(2R,h). If d does not divide the right side, that parity is impossible. Otherwise k lies in one explicit residue class modulo h/d, found by a modular inverse (modulus1 is handled separately).

The exact finite candidate interval is contained in

    ceil(N/(h+R)) <= L <= floor((N-h+t e)/R).         (4)

The lower bound uses g(L)<=L, proved from the recurrence and oddness of g; the upper bound uses g>=1. For each candidate in (3)-(4), evaluate g recursively in O(logL) and check equation(1). This is an exact inverse for the fixed prefix, not just a necessary filter. The parity and interval tests omit no solution; the final equality discards every false candidate.

A prefix bank with R<=B is complete within that bound: start at(1,1,1), use every L<=B/R and retain outputs with R'<=B. Since R'>=h+R, every earlier prefix of a retained state also satisfies the cap. Deduplicate the full triple, not just h. Applying the exact final inverse to every state in this finite bank decides precisely its one-final-gadget image.

## 4. A failed fixed bank and adaptive repairs

The R<=500 bank has2,717 states. Its one-final-gadget image covers all odd N<=1,000,001 except the29 known C1 omissions listed below. However through10,000,001 it additionally misses

    5,099,741; 5,477,461; 8,231,473; 9,551,177.

These are failures of this fixed-prefix strategy, not omissions from the C1 spectrum. Increasing the bank cap to1000, with7,497 states, yields the following actual recipes (listed from the singleton outward):

- 5,099,741: [4,17,1,1,4500], order26.
- 5,477,461: [96,3,4061], order25.
- 8,231,473: [21,2,1,2,4869], order26.
- 9,551,177: [3,1,18,1,9628], order27.

`arithmetic-hard-witnesses.json` records the full graphs, rooted profiles and two-port words. Besides the scalar recurrence, the independent physical two-port compiler checks every graph arc and its uncapped27-coordinate cover transfer gives the stated H and C=1. No large adjacency subset-DP claim is made for these orders.

## 5. Finite positive coverage and all26 C2 repairs

The complete all-order rooted closure through N=5001 from `ROOTED.md` has exactly these29 odd omissions:

    7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,
    277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647.

The R<=1000 prefix bank followed by one arbitrary final gadget constructively covers every other odd value from3 through10,000,001. This is a finite positive-coverage assertion. Its failure to produce the29 values alone would not exclude other C1 constructions at larger targets; their exclusion here already follows from the complete closure through5001.

`arithmetic-c2-witnesses.json` supplies actual C2 graphs for all26 entries above other than7,21,67. Each recipe starts from a displayed C1 graph and named pair of start-one vertices. Add a new root beating that pair and no other old vertex, followed by the displayed transitive gadget tail; then apply the listed ordinary one-seam gadgets. The verifier builds every arc, counts all Hamilton paths and cycles directly, checks both source roots and the displayed (h,E,R,K) response, and never assumes a scalar pair-transfer identity. Twenty-four witnesses have order<=12. The last two have order13:

- 1861: nonadjacent pair response(9,5,25,4), pair gadgetL1 produces(25,4,34), then ordinaryL35.
- 2647: adjacent pair response(5,2,12,4), pair gadgetL2 produces(25,5,29), then ordinaryL63.

Here E is the sum of the two ending counts, R the sum of the two rooted two-cover counts, and K the sum of the two ordered nonempty-cover counts ending at one named vertex and starting at the other. These response data are provided for subsequent arithmetic work; positive coverage depends only on the directly verified graphs.

Consequently every odd integer 3<=N<=10,000,001 other than7,21,67 occurs in the strong C<=2 class. The absence of those three from this class is separately certified by the frozen complete two-port compiler/control packet through501. Without using that external negative certificate, the present packet asserts only positive coverage and C1 exclusion for the29 listed values.

## 6. Reproduction, scope and next obstruction

Run `python3 verify_arithmetic_inverse.py` and `python3 -O verify_arithmetic_inverse.py`. The verifier compiles `arithmetic_coverage.cpp`, checks both fixed prefix caps through10,000,001, hashes the complete packed value images, checks small cone/inverse controls, verifies all26 C2 recipes and the four large C1 recipes, and writes `arithmetic-inverse-certificate.json`. Every check uses explicit exceptions, not Python asserts. The exact source dependencies and hashes are recorded in the receipt.

The all-order open question remains whether C<=2, or C<=3, covers every sufficiently large odd integer, or exactly all odd integers other than7,21,67. A finite bank surviving ten million does not prove it survives larger targets. The formal unimodular inverse also fails to preserve physical reachability. A successful uniform descent must retain enough root/interface information to prove that its smaller predecessor is physically realizable.
