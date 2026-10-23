# Session41: three marked units, zero deletion, and coherent Fibonacci rows

Status: PROVED elementary statements under explicit conventions; FINITE-EXACT executable checks. No tournament or Collatz theorem is inferred.

## Sources and corrected inheritance

Read `05-knowledge/results/zeckendorf_guard_automaton_20261003.md` first for current truth: it already proves the compressed-three-unit representation fibres (E1), distinguishes the owner's stripe from the auxiliary F2^2 charge, and records the old row35/36 continuation obstruction. That source routes to `reset_20260926_colours.md` and `crossroads_crossing_20260926_colour.md`; the principal earlier candidate is `collatz_crossings_20260926_potential_and_seeds.md`, §1.1. `fibonacci_two_copies_collatz_20260926.md` provides the zero-removed two-sided context, with explicit audit repairs.

The older authoritative stripe translated R→C,K→A,B→B puts B at35. The canonical word below puts C there. The current 16-symbol prefix agrees with the canonical word but does NOT retract or explain the older35-symbol data. Consequently the following is a rigorously reconstructed canonical model of the current prefix, not proof that it uniquely determines the user's intended infinite coloring. Block recurrences and compressed fibres are inherited; the retained-zero comparison and all-size proper-expansion count sharpen their boundary analysis.

## 1. Two grammars, not two sets of lengths

Let F0=0,F1=F2=1 and Fn=F(n−1)+F(n−2). Mark A=F(-1)=1, B=F1=1, C=F2=1. All representations use distinct slots with no adjacent selected slots. The zero slot is always unselected; selecting it would create irrelevant zero multiplicities.

R (retained-zero): slot indices are −1,0,1,2,3,... and adjacency is the original integer-index adjacency. Because0 is absent, A is compatible with B.

D (deleted-zero): reconnect the surviving slots into the path A,B,C,F3,F4,... . Now A and B are adjacent and incompatible. No other positive-slot adjacency changes.

Merely ignoring zero length does not justify replacing R by D. The first disagreement is2=A+B, legal in R but not D;2=A+C is legal in both.

Write Z(n) for standard Zeckendorf support on indices≥2, Z(0)=empty, and r(n)=1 when2∈Z(n), otherwise0. Negative arguments mean no support.

**Theorem1 (complete fibre comparison).** For n≥0, D has precisely the applicable supports

- Z(n);
- A+Z(n−1), n≥1;
- B+Z(n−1), n≥1 and r(n−1)=0.

R has exactly these plus A+B+Z(n−2), n≥2 and r(n−2)=0.

Proof. Partition a support by which of A,B it uses. Its remaining support is a unique ordinary Zeckendorf support. B excludes C; A excludes B only in D. These are all restrictions. This also proves the converse for every listed support. The D part is inherited E1; the fourth R sheet is the zero-sensitive addition.

For n≥1, |D(n)|=3−r(n−1). For n≥2,
|R(n)|=4−r(n−1)−r(n−2).
Since consecutive integers cannot both have a lowest Zeckendorf term1, the latter is always3 or4. One proof of nonconsecutivity: adding1 to a support ending inF2 triggers the ordinary carry chain, whose normal form has noF2. Thus deletion is a surjection onto numerical values but removes a specific marked sheet; it is not a bijection on supports.

For finitely many slots, let Z_tail(x) count independent subsets of F3,F4,... with exponent their numerical sum. Then the exact difference of weighted representation polynomials is
P_R(x)−P_D(x)=x² Z_tail(x).
This follows by stripping the forced pair A+B; C is forbidden and F3 onward is unrestricted except consecutive adjacency. It retains incidences, not just lengths.

## 2. Proper self-decompositions require a selection rule

A proper decomposition of Fn uses only slots below n. Nonconsecutivity alone does not select one: F4=F3+A=F3+B in both grammars.

**Theorem2 (all-size ambiguity count).** For every n≥3, Fn has exactly2 proper R decompositions. It has exactly1 proper D decomposition when n is odd and exactly2 when n is even.

Proof. The fibre includes Fn itself and every other support uses smaller slots. Ordinary Zeckendorf expansion of Fn−1 has lowest index2 for n odd,3 for n even. Expansion of Fn−2 has the opposite presence of index2. The explicit expansion of Fn−1 is F(n−1)+F(n−3)+... (ending inF2 for odd n, F3 for even n). For odd n, Fn−2=F(n−1)+F(n−3)+...+F4 (empty for n=3). For even n, Fn−2=F(n−1)+F(n−3)+...+F5+F2 (the displayed odd-index sum is empty for n=4). Each identity follows by telescoping and has nonadjacent support. Thus r(Fn−1)+r(Fn−2)=1. Apply Theorem1 and remove singletonFn.

Choose the lexicographically largest support, comparing descending slot indices. This is an additional canonical rule, not a consequence of legality. Its proper support is
E_n={n−1,n−3,n−5,...≥2} ∪ {−1 if n odd, 1 if n even}.
Its weights sum toFn by telescoping the Fibonacci recurrence; it is legal in both models. Greedy descending choice forces these indices and the terminal marked unit, so it is the unique lexicographically largest proper support. In particular canonical proper expansion never uses bothA andB. The forgotten zero changes legal alternatives, but not this selected expansion.

## 3. Exact coherent row theorem and substitution

Set W_A=A,W_B=B,W2=C. Expand each Fn recursively using E_n, concatenating summand blocks in ascending slot order. Then
W3=AC, W4=BAC, and W_n=W_(n−2)W_(n−1) for n≥5.
Proof: E_n below its largest slot n−1 is E_(n−2); its expanded block is W_(n−2). Base n=5 follows from A+C+F4. Length isFn.

Reverse blocks: V_n=reverse(W_n). Thus V3=CA,V4=CAB and V_n=V_(n−1)V_(n−2). They are nested prefixes and define an infinite wordV beginning
CABCACABCABCACAB...

**Theorem3 (the added-left-unit property).** For every N≥1 define rowR_N by taking ordinary Z(N), expanding every summand throughW, and concatenating in ascending summand order. Then reverse(R_N)=V[1..N]. Consequently R_N=V[N] R_(N−1): the new unit appears at the left and all earlier colored units remain unchanged.

Proof. Induct onN. For Fn≤N<F(n+1), put m=N−Fn<F(n−1). Greedy Zeckendorf yields reverse(R_N)=V_n reverse(R_m). For n≥4 the limiting word begins V_n V_(n−1), and by induction reverse(R_m) is its prefix of lengthm. This gives exactly the firstN symbols. HandleN=1,2 directly.

An exact three-letter substitution is
σ(C)=CA, σ(A)=B, σ(B)=CA.
It sends V3 toV4 andV4 toV5, and therefore sends everyV_n toV_(n+1) by the concatenation recurrence. The limit is the fixed point beginningC. Its incidence matrix in order(C,A,B) is
[[1,0,1],[1,0,1],[0,1,0]],
with characteristic polynomial λ(λ²−λ−1). Hence this explicitly selected word's limiting frequencies are C:A:B=φ:φ:1, normalized to (φ^−2,φ^−2,φ^−3). For a proof without invoking a substitution theorem, use block count vectors and the same Fibonacci recurrence; normalization converges by Binet's formula and the finite-prefix remainder decomposition.

Equivalently, V[N]=C if minZ(N)=2, A if that minimum is odd≥3, and B if it is even≥4. Induct using N=Fn+m: for m>0 the lowest index is that ofm, and for m=0 the final letter ofV_n isA orB according to n's parity. This recovers the earlier Wythoff candidate by an explicit map from rows rather than by fitting a prefix.

### All-size fibre counts read directly from the word

Every C inV is followed byA and every A is preceded byC. This follows from the substitution images CA,B,CA: a C occurs only at the beginning of CA, and an A only at its end. Consequently for N≥1,

|D(N)|=3−1_{V[N]=A}.

For N≥2, the adjacent pair r(N−1),r(N−2) sums to0 exactly when V[N−1]=B, and otherwise to1. Thus

|R(N)|=3+1_{V[N−1]=B}, with |R(1)|=3.

The limiting mean fibre sizes are therefore 3−φ^−2=φ² forD and 3+φ^−3=2φ forR. The mean number of sheets removed by deleting zero is 2φ−φ²=φ−1=φ^−1. These are exact asymptotic statements for these two grammars, not assertions about arbitrary chosen colored representations.

### Bounded discrepancy, not just limiting density

For n≥3, V_n has exactly F_(n−2) copies each ofC,A and F_(n−3) copies ofB. This follows from the base blocks and their Fibonacci recurrence. Binet gives the exact error
F_(n−2)−φ^−2 F_n=−(−φ)^−n.
Every prefix is a concatenation of distinct Zeckendorf blocks, possibly including the exceptional blockV2=C. Summing absolute errors over all possible indices n≥3 is bounded by sum φ^−n=φ^−1. Thus every prefix has uniformly bounded discrepancies from the frequencies, for example absolute error below2 for each color. Consequently

sum_(1≤n≤N)|D(n)|=φ²N+O(1),
sum_(1≤n≤N)|R(n)|=2φN+O(1).

These bounds follow from the proved block model, not finite density fitting. They sharpen the limiting means while remaining conditional on the explicitly chosen canonical expansion convention.

## 4. Shared-port interpretation, with limits

Source object: weighted independent supports on a slot path, with distinguished boundary portsA,B,C. Operation: delete the zero vertex and reconnect its neighbors. Preserved: numerical weights, ordinary tail, canonical greedy word. Lost: compatibility of the pairA,B and therefore a whole representation sheet. The cheapest decisive test is supportA+B at value2.

This is an exact example of why an apparently empty cell can retain incidence information. It is analogous to keeping endpoint ownership before gluing tournament cells. No map from these Fibonacci supports to tournament Hamiltonian paths is supplied; in particular the two hidden units do not prove a {7,21} obstruction. Likewise the three letters do not establish an invariant for3n−1 or Collatz. Those are open connection problems, not consequences of Theorems1–3.

## 5. Reproducibility

`check.py` independently enumerates both legal-support grammars through all valuesN<610, checks the complete fibre formulas, checks recursive rows and lowest-index colors forN<6765, verifies the requested16-symbol prefix and the old35 discrepancy, and tests the corresponding Beatty positions with integer square roots. `receipt.json` and `receipt-optimized.json` are identical normal/python−O outputs. No test depends on Python assertions. Beatty position tests are FINITE-EXACT here; the all-size statements proved above do not rely on them.
