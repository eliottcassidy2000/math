# Raw source energy under re-pairing: exact identities, sharp obstruction, and a compatible-tree survivor

**Status: PROVED elementary algebra / FINITE-EXACT original-source controls / OPEN native ROOT debt target.**

This note uses the original R3 raw energy ledger, not conditional-expectation norms and not positive pair contributions. The source is pinned to commit `96e35642110104cc547516a63f37eebcf930f707`, `05-knowledge/results/paid_root_cycle_budget_20261026/README.md`, R1–R3. The source-charge formula is E1 in `05-knowledge/results/root_eisenstein_first_arrival_20261018j.md`. No theorem ID or external priority claim is assigned.

## 0. Inheritance, scope, and the recovered numerical configuration

- Closest proved mechanism: R3 allocates every same-terminal pair of original labelled sources to its first differing low-bit level exactly once.
- Canonical hostile: entrywise positive/absolute interactions erase signed cancellation; reduced expected energy does not bound the positive part of each realized energy difference.
- Corrected near miss: the three four-child pairings have mean middle energy `(P+2C)/3`; this removes debt only in the **averaged energy ledger**, not necessarily in any actual pairing or in expected actual debt.
- Least-used sidecar: the complete labelled membership of every temporary bin, with original source weights, native integer clock, arrival cutoff, terminal and entry port retained separately.
- Live concepts: raw quadratic energy; whole-level positive variation; equal-size block design; compatible balanced trees; deterministic matching selection; native low-bit chronology.

The trusted user history mentions positive pair contributions `433` versus prescribed whole-level positive debt `79` at a finite cutoff. **A fully specified native configuration producing exactly these numbers is now FINITE-EXACT verified:** `N=512`, all 768 original sources at complete finite first arrival, maximum native depth 42. It has `S_N=1244`, terminal energy `1313`, whole-level debt `79`, positive parent-bin interactions `293`, and positive labelled-pair interactions `433`. The complete labelled source ledger is supplied. Historical provenance of the user-quoted computation is still unidentified; numerical agreement does not prove this was its identical historical cutoff. Independent additional controls include `N=1024`, and `N=4096` at depth cutoffs `1,3,10` and complete finite arrival.

The generic hostile source configurations below are honest integer Eisenstein source configurations. They are not asserted to be the complete prescribed family `A_N(s)` at any `N`. A native-scope failure is claimed only where the actual native source construction is run.

## 1. The exact native ledger

Write `omega^2+omega+1=0`, `Q(a+b omega)=a^2−ab+b^2`, and

`beta(x,y)=Q(x+y)−Q(x)−Q(y)`.

For each original primitive source `s>0`, odd and not `5 mod 8`, retain

`(s,u_s,p_s,t_s,d_s,A_N(s))`,

where `u_s` is the first frontier terminal, `p_s` its exterior entry port (zero for a source already in the frontier), `t_s` the complete integer braid clock, and `d_s` the native edge depth. The literal source is

`A_N(s)=sum_{j:B^j(s)<4N} omega^(B^j(s)) * (omega−omega^2 if B^j(s)<2N, otherwise omega)`, `B(x)=4x+1`.

In particular zero-charge labels are retained. Set `v_s=(−1)^t_s A_N(s)` and `m_s=(s−1)/2`. For the sources arrived by a specified finite cutoff,

`W_(u,r)^(q)=sum_{u_s=u, m_s=r mod q} v_s`,

`E_q=sum_{u,r mod q} Q(W_(u,r)^(q))`, `q=1,2,4,...`.

Every residue bin exists, including empty bins and bins whose nonempty source sum vanishes. Multiplicities of equal charge vectors are not identified. With `q` separating the source labels, `E_q=S_arrived`; the full denominator `S_N=sum_s Q(A_N(s))` continues to include moving/unarrived sources.

The native identity is

`Delta_q=E_q−E_(2q)=sum_(u,r mod q) beta(W_(u,r)^(2q),W_(u,r+q)^(2q))`.

For distinct same-terminal source labels `s,z`, their interaction appears precisely at `q=2^v2(m_s−m_z)`. No cross-terminal interaction occurs. The native debt is

`D_native^+=sum_q max(Delta_q,0)`.

The positive part is taken **after summing the entire level over all terminals and all parents**. In general

`D_native^+ <= sum_(q,u,r) max(beta(child sums),0) <= sum_(same-terminal labelled pairs) max(beta(v_s,v_z),0)`.

Both inequalities can be strict. These three quantities have different meanings and the checker prints all three. It verifies the first-differing-bit allocation by direct labelled-pair enumeration independently of the bin-energy computation.

No normalized martingale energy is substituted for `E_q`. First arrival, the full original source and the low-bit source tree are three distinct coordinates.

## 2. Perfect-matching and equal-block averaging

Let `x_1,...,x_n` be labelled bin sums in a real inner-product space, with `Q=||.||^2`; the Eisenstein norm is a case. **When there are several terminals, use the orthogonal direct sum `H=direct_sum_u H_u`.** For example, a residue-slot vector is `x_r=direct_sum_u W_(u,r)`. Then `P=||sum_r x_r||^2=sum_u Q(sum_r W_(u,r))`, never `Q(sum_(u,r) W_(u,r))`. Terminal tags remain present through all additions, and no cross-terminal interactions are introduced. Labels can denote complete source sub-bins, not just individual charges. Put

`P=Q(sum_i x_i)`, `C=sum_i Q(x_i)`.

### Proposition 2.1 (uniform perfect matching)

For even `n>=2`, choose uniformly among the `(n−1)!!` perfect matchings of the **labelled** children. If

`M=sum_{pairs {i,j}} Q(x_i+x_j)`,

then

`E[M] = C+(P−C)/(n−1) = [P+(n−2)C]/(n−1)`.

**Proof.** A fixed unordered pair is matched with probability `1/(n−1)`. Expand `M=C+sum_matched beta(x_i,x_j)` and use `sum_(i<j) beta(x_i,x_j)=P−C`. This proof retains empty bins and equal-valued but distinct labels. QED.

For `n=4` the three perfect matchings give exactly `sum M_i=P+2C`. For native four-grandchild bins at modulus `4q`, use the ordered addresses `r,r+q,r+2q,r+3q`; the actual modulus-`2q` grouping is `(0,2)|(1,3)`, not `(0,1)|(2,3)`.

### Proposition 2.2 (uniform equal-size blocks)

If `b|n` and a uniformly random partition has blocks of size `b`, then its raw bin energy `E_b` satisfies

`E[E_b] = [(b−1)P+(n−b)C]/(n−1)`.

A fixed pair is in the same block with probability `(b−1)/(n−1)`, which proves the formula by the same expansion. This includes `b=1` and `b=n`.

### Compatible balanced hierarchy

For `n=2^h`, take one uniformly random permutation of the labelled children and use consecutive blocks of sizes `n,n/2,...,1`. These partitions are nested; hence they form an actual compatible binary tree. All their marginal expectations satisfy Proposition 2.2. The expected energy sequence lies on the segment between `P` and `C` and is monotone. Consequently

`D^+(the expected energy sequence)=(P−C)_+`.

This is not `E[D^+(the realized tree)]`. Jensen gives the reverse-direction lower bound

`E[D^+(tree)] >= (P−C)_+`.

Sampling different partitions independently at each size preserves the marginal formulas but need not produce a tree. For example, on eight labels, coarse blocks `{0,1,2,3}|{4,5,6,7}` and fine pairs `{0,4}|{1,5}|{2,6}|{3,7}` are not nested. One needs the shared permutation, a bottom-up matching process, or another explicit compatibility sidecar.

## 3. Insertion, refinement and the exact obstruction

For three nonnegative energies `P,M,C`, the two-step debt satisfies

`(P−M)_+ +(M−C)_+ = (P−C)_+ + dist(M,[min(P,C),max(P,C)])`.       (3.1)

This follows by the three cases `M` below, in, or above the endpoint interval. Inserting an intermediate refinement never decreases debt. It leaves debt unchanged exactly when the inserted energy lies between its endpoints. Deleting a level never increases debt.

For fixed endpoints, replacing one intermediate partition changes debt by exactly the change in the distance term. In particular its absolute effect is at most `|M−M'|`. For any two finite energy sequences of the same length and endpoints,

`|D^+(E)−D^+(F)| <= sum_(interior k) |E_k−F_k|`.

Indeed `D^+(E)=[TV(E)+E_0−E_h]/2`, and triangle inequalities bound the variation difference by twice the displayed sum.

Averaging the three middle energies places their mean in the endpoint interval, but does not force one of the actual middle energies into that interval. The first failed implication is moving a positive-part operator through expectation.

### Sharp four-child counterexample

Take

`(x_1,x_2,x_3,x_4)=(2−omega,2−omega,omega,omega)`.

Then

`P=Q(4)=16`, `C=7+7+1+1=16`,

while the three middle energies are `32,8,8`. Their mean is `16`, but their actual two-step debts are `16,8,8`. Thus:

- the averaged energy path has debt zero;
- no actual pairing has debt zero;
- the minimum actual debt is `8`;
- expected actual debt under uniform pairing is `32/3`.

This is a source-preserving counterexample. Use eight labelled unit charges in the four bins:

`x_1: 1,1,−omega; x_2: 1,1,−omega; x_3: omega; x_4: omega`.

Every named source remains present and the full source norm is `8`. The child-bin energy `16` is deliberately distinguished from full source norm `8`; the counterexample is about the debt inserted between the stated parent and child partitions. It does not pretend those child bins are already the separated-source leaves. The scalar integer alternative `(5,3,−1,−1)` has endpoints `36,36` and middle energies `68,20,20`.

### An original-native minimal-support counterexample, including an empty bin

The same failure occurs inside the prescribed source family. At `N=512`, terminal `435`, parent residue `0 mod 1`, the four `mod 4` child bins have sums

`(omega, omega, 0, −1−omega)`.

They contain respectively source `1833`, source `435`, no source, and source `1375`. Their original charges have the displayed values and all their clocks are even. The complete records are:

| source | terminal | entry port | integer clock | native depth | odd-index m | original A |
|---:|---:|---:|---:|---:|---:|---|
|1833|435|4643|10|5|916|omega|
|435|435|0|0|0|217|omega|
|1375|435|4643|8|4|687|−1−omega|

No additional source belongs to this terminal at the complete `N=512` cutoff. The endpoints have `P=C=3`, and the three middle energies are `5,2,2`: every actual pairing has positive local debt, `2,1,1`. The averaged path has zero debt. The actual low-bit pairing is the second pattern, `(0,2)|(1,3)`, with middle energy `2`.

There is also an example at a scale covered by the native-open target. At `N=4096`, terminal `443`, parent residue `0 mod 4`, the four `mod 16` child bins have sums `(omega,−1−omega,−1−omega,0)`, supported respectively by sources `12609`, `13993`, `9457`, and the empty set. Their clock/depth pairs are `8/2`, `24/12`, and `6/1`, all with entry port `9457`. Again `P=C=3`, and middle energies are `2,2,5`. This local statement holds once these sources have arrived, in particular at depth 12, and remains so at complete finite arrival. Other parent bins may cancel this local debt at a global level; no global lower bound is inferred from it.

Three nonzero child bins are minimal for this obstruction: with at most two nonzero child charges, every middle energy equals either `P` (the two charges paired together) or `C` (paired separately with zeros), so no pairing inserts excess. The examples use only three actual nonzero original source charges and preserve the genuinely empty fourth child bin.

At the whole-level scale `q=16` for the same complete `N=512` data, endpoints are `P=1255`, `C=1252`. Using each one of the three pairing patterns uniformly at every parent produces global middle energies `1251,1242,1266`; none lies in `[1252,1255]`. Their debts are `4,13,14`, whereas the averaged middle energy `1253` gives endpoint debt `3`. The second pattern is native. This is a compatible *three-member correlated ensemble at one inserted level*. It does not exclude an alternative using a different pairing at each parent.

### Sharp local bounds and deterministic selection criterion

For four children, every middle energy obeys

`P/2 <= M_i <= 2C`, and `sum_i M_i=P+2C`.                       (3.2)

The lower bound is `Q(y)+Q(p−y)>=Q(p)/2`; the upper bound is the two-source Cauchy inequality inside each pair. Define `L=min(P,C)`, `U=max(P,C)`.

**Proposition 3.2.** Some pairing has excess over endpoint debt at most

`B_det(P,C)=min(P/2, max(C−P/2,0))`.

The uniform-pairing expected excess is at most

`B_av(P,C)=[dist(2C,[L,U])+2dist(P/2,[L,U])]/3`,

namely

- `(P+C)/3` if `P<=C`;
- `2(2C−P)/3` if `C<=P<=2C`;
- `0` if `P>=2C`.

**Proof of deterministic bound.** If `P<=C` and no `M_i` lies in `[P,C]`, at least one must be below `P` because their mean lies inside the interval. Such an energy is at least `P/2`, so its excess is at most `P/2`. If `C<=P<=2C`, the same mean argument supplies a value below `C` unless one is already in the interval (a value above `P` forces a compensating below-`C` value); its distance is at most `C−P/2`. If `P>=2C`, (3.2) puts every `M_i` inside `[C,P]`. QED.

**Proof of expected bound.** The feasible triple lies in the simplex `M_i>=P/2`, `sum M_i=P+2C`; its vertices are permutations of `(2C,P/2,P/2)`. The sum of distances to `[L,U]` is convex, so its maximum occurs at a vertex. QED.

These bounds are sharp as real-inner-product inequalities. Given `p,a`, set `x_1=x_2=(p+a)/4` and `x_3=x_4=(p−a)/4`. Then `P=Q(p)`, `C=(Q(p)+Q(a))/4`, and the middle energies are `(2C,P/2,P/2)`. Every admissible real pair `C>=P/4` is realizable by choosing the norms of `p,a`. The integer Eisenstein witness above attains both sharp bounds at `P=C`.

In particular at equal endpoint energies `P=C`, some pairing always has excess at most `P/2`, but that constant cannot be improved; expected excess is at most `2P/3`, also sharp. If `P>=2C`, **every** pairing incurs precisely the endpoint debt. If `P=0`, at least one pairing incurs zero debt. For arbitrary prescribed data, a zero-excess deterministic pairing exists **if and only if** some actual candidate middle energy belongs to `[L,U]`; the mean identity alone is not such a certificate.

### Whole-level aggregation is essential

For many parent bins, the candidate global middle energy is the sum of one pairing energy from each parent. Zero excess is equivalent to this attainable sum meeting the **global** endpoint interval. Demanding zero excess independently in each parent is stronger and is not the R3 target.

For `B` disjoint terminal copies of the sharp witness, global endpoints are `16B` and attainable middle energies are `8B+24k`, `0<=k<=B`. Hence zero global excess exists exactly when `3|B`, despite every individual parent having unavoidable positive excess. For other `B`, the minimum global excess is `8`. This also shows why taking positive parts before the whole level is not harmless.

## 4. A genuine global deterministic survivor (but for a changed tree)

### Theorem 4.1 (minimum-matching compatible tree)

For `n=2^h` labelled bin charges in the terminal-tagged direct-sum Hilbert space and terminal-sum energy `P=sum_u Q(total_u)`, there is a compatible balanced binary tree, preserving every leaf charge, terminal tag and identity, such that

`D_tree^+ <= (1−1/n) P`.                                    (4.1)

The construction repeatedly chooses, bottom-up, a perfect matching minimizing the next raw energy. The bound is sharp in the class of all source configurations: equal charges attain it.

**Proof.** At a stage with `m` labelled current bins of total energy `E`, the total charge is unchanged, so its squared norm is still `P`. By Proposition 2.1 the minimizing matching has next energy

`E' <= [P+(m−2)E]/(m−1)`.

Its positive gain in the bottom-up direction, which is exactly the debt of this edge in root-to-leaf orientation, satisfies

`(E'−E)_+ <= (P−E)_+/(m−1) <= P/m`,

using `E>=P/m`. Sum over `m=n,n/2,...,2` to obtain `(1−1/n)P`. Each new bin is literally the union of its two predecessor bins; thus the selected matchings are compatible by construction. For `n=1`, there are no edges and the bound is zero. For `n` equal charges `v`, `P=n^2Q(v)` and every tree has debt `P−nQ(v)=(1−1/n)P`. QED.

Consequences and boundaries:

1. If the total charge is zero, some compatible balanced tree has debt zero, irrespective of the original source norm or the debt of a supplied tree.
2. For several terminals, either use the terminal-tagged direct-sum vectors in the theorem directly (one selected grouping shared by the residue slots), or apply the construction separately within each terminal, padding to a common power-of-two leaf count with empty bins if needed. In the latter construction the global whole-level debt is at most the sum of terminal debts. Both procedures preserve terminals; neither computes the square of a sum over different terminals.
3. The procedure requires the actual charge sums to select matchings; it is not a universal fixed pairing independent of the input.
4. Minimum matching is a finite selection certificate. No efficient global algorithm or native dynamical realization is asserted here.
5. Most importantly, (4.1) controls an **alternative** source-partition tree. It neither bounds the original low-bit `D_native^+` nor makes source refinement into native trajectory chronology. It is also circular as a prospective bound on `P` from source norm: its right side is already the target terminal energy.

### Unbounded fixed-tree obstruction with unchanged unit-source norm

Let `n=2^h`, use one terminal, label leaves by the odd-source coordinate `m=0,...,n−1`, and assign charge `v_m=(−1)^m`. Under the fixed low-bit partitions:

`E_1=0`, `E_(2^k)=n^2/2^k` for `1<=k<=h`, `S=n`.

Therefore

`D_fixed^+=n^2/2−n`, so `D_fixed^+/S=n/2−1`.

The terminal sum and all original unit charges are unchanged. An alternative tree pairing opposite signs first has debt zero. Hence terminal sums plus source norm do not bound the debt of a fixed source tree by an absolute multiple of source norm. This is a rigorous generic partition obstruction, **not** a counterexample to the prescribed native `A_N(s)` family, which has additional arithmetic and clock constraints.

### A native marked-observer sidecar

At `N=512`, source `799` is already at terminal `799`, with `A=v=−2−omega`, clock/depth `0/0` and port `0`. Source `1065` has `P(1065)=799`, clock/depth `2/1`, exterior port `1065`, and `A=v=omega`. Their source addresses are `399` and `532`, so their interaction belongs to native level `q=1`. The two charges satisfy `beta_Q=0`, but for the actual marked readout `c=2a−b`, their marked cross term is `2*(-3)*(-1)=6>0`. Thus a vanishing norm-polarization pair need not be invisible to the marked observer. Both original source records occur in the provided ledger; this is not merely a freely chosen vector example.

## 5. What was established, and what remains native-open

**PROVED:** higher-arity matching and equal-block means; compatible random hierarchy means; exact insertion-excess formula; sharp deterministic and average four-child excess bounds; the integer Eisenstein no-zero-debt witness; a sharp compatible-tree bound obtained by minimum matching; unbounded fixed-tree debt for generic labelled unit sources.

**FINITE-EXACT:** the checker retains original native labels, zero charges, full source denominator, first-arrival clocks, entry ports and actual low-bit bins for its declared scales/cutoffs. Direct labelled-pair allocation agrees with energy differences. The check uses explicit exceptions rather than Python `assert` and is replayed under normal Python and `python -O`.

**FINITE-EXACT recovered instance:** `433` positive labelled-pair contributions versus native whole-level debt `79`, at `N=512`, complete first arrival, maximum depth 42. The historical computation itself remains unlocated; no provenance identity is claimed.

**OPEN:** `D_native^+(N,d)<=S_N/2` for every `N>=4096` and every finite native first-arrival depth cutoff. The present averaging and selection theorems do not prove it. They identify the required missing sidecar: a signed, whole-level, original-tree argument or a quantitative cost for transferring from an alternative compatible tree back to the original low-bit tree.

## 6. Reproduction and artifacts

Run in any directory containing `verify.py` and the pinned `frontier.txt`:

```bash
python3 verify.py > receipt.json
python3 -O verify.py > receipt.optimized.json
cmp receipt.json receipt.optimized.json
```

The script generates `native_sources_N512.json`, `native_sources_N1024.json`, and `native_sources_N4096.json`, containing all original source labels and their full routes' terminal/port/clock/depth/charge data. The receipts give SHA-256 digests of those exact files. The frontier SHA-256 is `df88f03b576b13c98c0c5378fc1cf61685eeb6817073681266d2a781dc329b26`; it has 627 states and is checked forward-closed.

The finite algebra universe is all `5^4=625` ordered four-tuples from the stated five-vector palette, plus all perfect matchings for declared samples at arities `2,4,6,8`. Compatible random-tree means and debt are enumerated over all `8!=40320` permutations of eight named leaves. Minimum-matching construction is checked on all `2^8=256` scalar sign words and coherent sharpness controls. Empty and zero-sum bins are not silently deleted: the native receipt states conceptual bin counts and unoccupied counts at every level.

The elementary proofs, not finite tests, support the universal algebraic claims. No external literature claim is made.
