# Complete fixed-module returns for34·2^k+5

**PROVED reduction plus independently recertified FINITE-EXACT head.** The
classification concerns a fixed tournament R with one vertex replaced by
TT_k, retaining every outside comparison uniformly. It permits arbitrary
fixed seed order and neutral SCC padding. Variable-incidence prefixes and
multiple growing modules are outside its scope.

## 1. Why an infinite-family classification has a finite head

For a pointed tournament(R,v), let A_j count spanning words on the other
vertices with j separated occurrences of the v-placeholder. Replacing v
by TT_k gives the inherited exact nonnegative expansion

    H(R[v←TT_k]) = sum_(j>=1) A_j j! S(k,j),       k>=1.              (1)

Each ordered nonempty assignment of the k module vertices to j runs
contributes one path, ordered internally by transitivity. Conversely every
path has a unique run decomposition. All A_j are nonnegative integers.

Suppose this count equals34·2^k+5 for infinitely many unbounded k. Any
positive A_j for j>=3 would grow as j^k and contradict the equality.
Thus only j1,2 occur. Using2!S(k,2)=2^k−2 forces

    A1=73, A2=34, A_j=0 for j>=3.                                  (2)

In particular the equality extends to every k>=1 and H(R)=73. By the
published sharp strong-order bound f(9)=75, a strong such R has order<=8.
We use this bound exactly as in the inherited
[ready inverse](../../ready_return_pentagon_20261103/OBSTRUCTION.md#o2-all-order-ready-value-inversion-has-a-finite-constructive-head),
which cites [Busch2006 Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

For a weak R, its ordered SCC counts multiply to the prime73. Exactly one
component has H73 and all others have H1, hence are singleton components.
If v lies outside the H73 component its inflation leaves H constant73,
contrary to(2). Otherwise the family is an active strong pointed template
with arbitrary fixed transitive padding before and after it. This proves
the claimed reduction for every fixed seed order.

## 2. The complete answer is three families and their converses

The complete order<=8 universe has19 strong H73 classes: eight at order7
and eleven at order8. Exactly six pointed classes satisfy(2). They are
the following three and their converses, with free prefix point0:

* A: T_a(1), the C4 tail01,11,10;
* B: T_b(1), the C4 tail01,01,10;
* C: (126,236,120,16,226,72,8,109), the preceding uniform replacement
  family from [RESET section5](../../chanpak_five_core_bridges_20261110/alternatives/RESET.md).

The receipt retains every row and prefix-fixing isomorphism, so this is
a classification of pointed structures rather than equality of scalar
responses. None of the six is a new source-return family.

The independent recertification uses all7412 inherited representatives
of orders1..8. For EVERY representative it recomputes H and strongness,
canonicalizes by all score-compatible permutations, and counts the exact
automorphism stabilizer. At each order the distinct class masses sum to
2^(n choose2). Distinctness plus that full labelled mass proves exhaustivity;
the inherited census label is not the sole justification. Every selected
module response is also checked on literal inflations k1..4.

## 3. The complete family list has an all-order source barrier

At k2 the three active strong sources have H141. Fresh complete strong
frontiers under425 have successive sizes

    A:25,4,0;      B:27,4,0;      C:27,11,0.

None contains425. The [low-ratio strong-ear theorem](../../chanpak_five_core_bridges_20261110/algebra/LOW_RATIO_FRONTIER.md)
applies since425<(n+2)141. Every proper divisor of425 is<141, so the
strong-source exclusion extends to arbitrary target tournaments. Converse
preserves induced inclusion and H; neutral padding cannot remove a retained
blocked component. Thus ALL six templates, with arbitrary neutral padding,
exclude425 at every possible target order when k2.

**Consequence for fixed-template programs.** No finite surgery outside an
intact uniformly seen TT_k can restore34·2^k+5 and then support a fixed
induced-extension template of count102·2^k+17. This remains impossible if
both identities are required only for an unbounded set of k: positivity
in(1) extends each identity to k2, where the second template would contain
an excluded141 source. There is no bound on the number of fixed auxiliary
vertices in this statement.

The conclusion does not forbid k-dependent incidence, several growing
modules, a different count family, or the successful smaller-count return
in [UNIFORM_RETURN](UNIFORM_RETURN.md). That construction deliberately
returns to BA_k at14·2^k+3 and proceeds directly to the desired target.

## 4. Additional complete and bounded probes

The older deleted-port parent W(k) has H=7·2^k+5. At k1 its cyclic two-ear
cost H(W)+Q_W(2)=275 exceeds221. Both H and Q(2) are monotone under induced
inclusion: insert each new vertex in the first nonempty bin, or into an
empty bin under a fixed convention. Thus its entire strong frontier under
221 uses single-vertex strong births, by the same relative-ear proof.

Prune a prefix template only if A_j>0 for j>=3 or A2>102, since such a
template cannot lie inside a102·2^k+17 family for every k. Canonicalize only
the newly added names, keeping all six W(1) labels fixed. The complete
frontier has sizes34,572,906,219,0 and no target template. This excludes
every fixed-template direct branch retaining precisely that deleted-port
parent. It does not exclude deleting a different old vertex; the positive
middle-tail deletion is the explicit surviving option.

Splitting the two prefix vertices of each of the22 old k2 port-deletion
reset hits into TT_m and TT_m (or TT_m and TT_(m+1)) also gives no strong
uniform return. `two_clock.py` computes the full two-module run polynomial,
expands it in the independent integer exponential bases, and checks literal
inflations. This is a complete22-hit probe, not a classification of every
two-clock construction.

The complete inherited order9 bank supplies two converse-merged pointed
templates for102·2^k+17, and none for102·2^k+19. The new positive family
is one of them. This order9 assertion does not exhaust the latter ray:
H223 permits strong base order10 (f(10)=125, f(11)=225). The five-core
normalized base223 inverse likewise yields no fixed single-module target,
but this grammar screen is not a complete order10 census.

Reproduction: `python -B classify_modules.py`, `python -B direct_branch.py`,
`python -B two_clock.py`, and `python -B scan9.py`. The last uses C++17 and
temporary files. The receipts state exact universes and retain witnesses.
