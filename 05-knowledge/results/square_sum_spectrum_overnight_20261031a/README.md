# A prescribed square-sum path at its exact threshold, and terminal-controller rigidity

**Status: PROVED with exact finite certificates for the two scoped theorems below. The unrestricted raw tournament spectrum remains OPEN.** Same-agent mathematical review uses independently implemented consumers, including a full-vector check of every finite square-sum cycle. There is no external priority claim. Ordinary interval square-sum Hamiltonicity is Gerbicz's known result, already recovered in the [previous packet](../square_sum_spectrum_overnight_20261030/README.md).

## Inheritance and live concept board

| Lane | Proved mechanism inherited | Hostile / corrected near miss | Preserved information and new test |
|---|---|---|---|
| Anchor: raw tournament spectrum | The exact 17-dimensional continuation representation and integer lattice | Integral positive states need not be actual words; 7 and 21 remain forbidden | Actual prefix orbit, terminal suffix, and all continuation responses |
| Anchor: changing terminal controllers | Fixed scalar contexts already fail | A terminal suffix might have restored an affine recurrence; the identity below excludes this for the unrestricted prefix language | Restrict the source language explicitly before trying another controller |
| Niche: prescribed square-sum forests | All 294 pair transitions and 392 final glue states; the 2303-path theorem | Optional edges need not be jointly compatible | All named singleton ports, cut directions, parity, and one connected glue cycle |
| Wildcard: reuse a hostile | The forced proper cycle in Q_24 | Opening that cycle avoids the immediate obstruction, but degree tests alone suggest the false onset 68 | Three forced incidences at vertex 36 explain why the exact onset is 72 |
| Shared mechanism | Exact quotient plus finite certificates | Linear span loses actual word entry; a degree cover loses one-cycle connectivity | Certify membership in the original object class, not just the projected statistic |

The least-used sidecar here is the explicit Q_24 obstruction in [session32_square_monster_20261023/square-sum](../session32_square_monster_20261023/square-sum/README.md). Its ordered cycle, rather than just its non-Hamiltonicity conclusion, supplies the new prescribed path. The tournament dependencies include THM-002 / [odd-cycle formula](../../../01-canon/theorems/THM-002-ocf.md), THM-1370 / [permanent raw exclusions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md), and THM-4094 / [two prime lanes](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md). Raw H and presentation count H/Aut remain separate throughout.

## 1. An exact prescribed-path threshold

Let Q_N have vertices 1,...,N, with two **distinct** vertices adjacent precisely when their sum is a square. Set

    P = (1,8,17,19,6,10,15,21,4,12,24).

**Theorem.** Q_N has a Hamiltonian cycle containing P as a consecutive path, in either direction, **if and only if N >= 72**.

Additionally, for every N >= 128453 the cycle can simultaneously contain the inherited path 482–3362–6242. At N >= 6294221 it can also retain all 2303 paths in the inherited [connected-path theorem](../square_sum_all_size_connected_paths_20261019/ALL_SIZE_CONNECTED_PATHS.md), including its four-vertex seed extension.

The ten sums along P are 9,25,36,25,16,25,36,25,16,36. Its unused closing edge 24–1 also has square sum. Requiring that closing edge as well would force a proper closed component for every N>24. The theorem deliberately prescribes the **open path**, retaining its two available endpoint incidences.

### Necessity: why 68 does not suffice

The internal vertices of P already use both Hamiltonian-cycle incidences. Thus an outside vertex cannot use any of them as a cycle neighbor.

- N<24: the required vertex 24 is absent.
- 24<=N<=29: vertex 18 has only neighbor 7. Its next distinct square-sum neighbor is 31.
- 30<=N<=31: the only neighbors of 30 are 6 and 19, both internal vertices of P.
- 32<=N<=67: the neighbors of 32 that could occur are 4,17,49. The first two are internal to P, leaving at most one available neighbor.
- 68<=N<=71: the available neighbors of 28,45,64 are respectively {36,53}, {36,55}, {36,57}. All three vertices therefore force an edge to 36, contradicting degree two at 36.

This last obstruction survives the individual degree tests. The new vertex 72 makes the edge 28–72 available because 28+72=100; it releases precisely one member of that forced three-edge star. This explains the observed threshold rather than just reporting solver infeasibility.

### Sufficiency through 128452: every integer has a certificate

The directory `native_n24_threshold/` contains:

1. All 57 explicit cycles for N=72,...,128 in `finite_head_cycles.json`.
2. A complete insertion/reversal trace from 128 through 128452. It introduces exactly the next vertex at each step. Every reversal preserves one cyclic ordering; no move removes a prescribed edge of P. The trace uses 531148 reversals, 20544 direct insertions, and at most 45 reversals in any one step. Search success is **not** assumed at any unsaved index.
3. `verify_native_trace.py`, a standard-library consumer using an implicit sequence and exact bad-edge induction. It also directly verifies all 57 head cycles and every elementary negative case above. Its 5857676 checks pass.
4. `verify_native_trace.cpp`, an independent full-vector consumer. It checks the complete vertex set, every cyclic sum, and all ten P edges at **every** N=128,...,128452. It passes 24757389843 checks, including 8250014250 square-edge checks. It also independently matches the stored final cycle.

The union covers all 128381 integers from 72 through 128452. The two consumers use different sequence representations; the C++ consumer does not use the incremental bad-edge proof. Neither consumer uses a Hamiltonicity solver. The discovery program is supplied for provenance, not as a proof oracle.

### All larger integers: a complete replacement final table

The inherited [eventual-cycle theorem](../seeded_eventual_cycle_theorem_20261019/SEEDED_INTERVAL_THEOREM.md) supplies B_n pairs for every n>=53 and both protected A families for every n>=2621. Its base bundle and all 294 transitions were replayed in the previous session; their exact dependency hashes and receipts are retained in [INHERITANCE.json](../square_sum_spectrum_overnight_20261030/INHERITANCE.json).

Given N>=128453, write

    n=floor((N-24)/49),   d=N-49n,   24<=d<=72.

Choose the same A family and the same lifted input members as in that theorem. For every offset -24<=c<=24, the affine block is

    49 S_i[j] + (-1)^j c,       positions j are one-based.

The common input positional parity makes these blocks partition [25,N]. Singletons 1,...,24 complete the interval. Cut only the certified 2–7 edge in each block. The resulting 98 segments have numerical endpoint pairs determined by n modulo two, d, and the direction of 2–7 in each of the two input paths.

`n24_forest_table.json` replaces **only** the final glue table. It supplies all

    2 parities * 49 residues * 2 directions * 2 directions = 392 states.

For each state, `verify_native_forest.py` checks all 98 segments and 24 singletons exactly once, the chosen input index and orientation, distinct named ports, every cyclic connector sum, all ten P edges, and the connector 482–3362. A single cyclic atom word, rather than an arbitrary degree-two cover, establishes connectedness. The inherited initial pair vertices 69,127 give the internal edge 3362–6242 at offset 19; the cut at 2–7 cannot disturb it. This proves simultaneous seed retention.

Internal square sums are multiplied by 49=7^2. All new external sums are explicitly checked. Consequently the old all-index induction with this replacement final table proves every N>=128453. This infinite implication is independent of the three large expanded controls in the receipt. Those controls at N=128453,132063,132128 additionally verify actual full cycles and both prescribed paths.

Finally, the old 2303-path argument concerns internal segments untouched by the 2–7 cuts. The new table uses the identical segments and cuts, so that entire proof carries over at its existing threshold 6294221. The new P lies in the singleton range [1,24], disjoint from those internal paths. No new density-to-universality inference is involved.

## 2. No nontrivial universal affine terminal insertion

Use the four-letter two-feedback convention of the previous packet: marked vertices x->y, a transitive chain, and bit r of a letter set when marked vertex r beats that chain vertex. H(w) is the raw Hamiltonian-path count. Empty words are allowed.

**Theorem.** Fix words A,B and constants a,b. The identity

    H(w B A) = a H(w A) + b       for every word w

holds exactly in the following cases:

    a=1, b=0, and either B is empty,
    or both B and A consist entirely of the letter 3.

The all-3 case simply adds sinks. In particular a fixed terminal suffix cannot support a universal nontrivial doubling/addition controller on the unrestricted prefix language.

### An all-length signed prefix identity

For every word z,

    2H(z) -2H(1z) -3H(2z) +3H(3z)
    -2H(10z) +2H(11z) +H(12z) -H(13z)
    +H(22z) -H(23z) +2H(102z) -2H(103z)
    -H(122z) +H(123z) = 2.                       (1)

The coefficients sum to zero. `terminal_rigidity.py` proves (1) in the original 36-state gap model, before the 17-state quotient. Let u,M_0,...,M_3,v be that exact representation, and D the signed sum of prefix rows in (1). The script constructs a basis C of the complete observable column space and verifies

    Dv=2,   (DM_i-D)C=0,   M_i C contained in span(C),  i=0,1,2,3.

The closed space has dimension 20 and contains v, so induction on the suffix length proves (1) for every suffix. All arithmetic is exact. Independently, a vertex-subset Hamiltonian-path algorithm checks (1) on all 21 suffixes through length two; these controls supplement, rather than replace, the closed-space proof.

Apply the proposed affine identity at each prefix appearing in (1) and take that signed sum. The left side is 2, the nonconstant right side is 2a, and the b terms cancel. Hence a=1.

### Every genuine insertion has unbounded positive debt

The odd-cycle-collection formula, THM-002 / [ocf](../../../01-canon/theorems/THM-002-ocf.md), says H(T) is the sum of 2^{number of cycles} over collections of vertex-disjoint directed odd cycles. Old collections embed unchanged when new vertices are added.

Assume B is nonempty. Put w=3^m, m>=1.

If some letter of B is not 3, its vertex b beats a marked vertex x. Choose any odd-cardinality subset of the m old chain vertices. In increasing order it makes a new directed odd cycle

    x -> chosen chain vertices -> b -> x.

There are 2^{m-1} such subsets. Each new singleton collection contributes two and uses the inserted vertex b. Therefore

    H(3^m B A) - H(3^m A) >= 2^m.                (2)

If B consists of 3s but A contains a non-3 letter, choose b in B and a vertex q in A that beats a marked vertex x. Every **even**-cardinality subset of the old chain, including the empty subset, gives a new odd cycle

    x -> chosen chain vertices -> b -> q -> x.

The same bound (2) follows. These are genuinely new collections because they contain b; no subtraction or cancellation is hidden in the argument. Thus the difference cannot be the fixed constant b in the proposed identity. The remaining all-3 padding case and the empty insertion are immediate. This completes the classification.

The checker also tests (2) on 1680 bounded positive/hostile cases. The all-order lower bound itself is the explicit cycle injection, not an extrapolation from those tests.

### A stronger rerouting test: all six closed binary source alphabets

**Theorem.** Let Sigma be any two-element subset of {0,1,2,3}, let A be any fixed terminal word, and let B belong to Sigma*. If

    H(w B A)=a H(w A)+b  for every w in Sigma*,

then a=1 and b=0. For nonempty B, the complete list of cases where the identity holds is:

| Sigma | Exactly the neutral insertions |
|---|---|
| {0,1} | B contains only 1s, and A contains only 1s and 3s |
| {0,2} | None |
| {0,3} | Both B and A contain only 3s |
| {1,2} | None |
| {1,3} | Both B and A contain only 3s |
| {2,3} | Both B and A contain only 3s |

Empty B is always neutral. This theorem concerns a fixed source alphabet **closed under the inserted block**; it is not a prohibition on switching source charts.

Here is the exact proof certificate for the new slope restriction. In each alphabet the checker supplies a signed sum D of actual prefix rows, whose coefficients sum to zero and which satisfies D M_s=D for both source letters s. For alphabets {0,3} and {1,2}, D M_A v is respectively the nonzero constant 24 or 4 for every suffix A. The original 36-state consumer independently verifies those constant identities before quotienting.

For the other four alphabets, the integral 17-state representation gives

    D M_A v = 2 (mod 4)  for every suffix A.

`binary_controller_rigidity.py` verifies this by enumerating the **entire closed row orbit modulo four**, eight states for each alphabet. This is an integral finite-state proof, not an invalid reduction of a rational basis modulo four. Every state has output 2 modulo four. The signed prefix expressions and all orbit states are recorded in `binary_controller_rigidity.json`; 126 independent vertex-subset evaluations check their literal tournament meaning.

Summing the proposed affine identity at those prefixes gives

    D M_B M_A v = a D M_A v.

Since B is a source-alphabet word, D M_B=D. The displayed nonzero integer forces a=1.

For completeness, all nonneutral rows of the table have an exponentially growing insertion debt, proved as follows. For alphabets containing 3, use w=3^m and the earlier cycle injection. For {1,2}, use w=2^m if B contains 1, and w=1^m otherwise; the chosen inserted letter beats the relevant marked vertex. For {0,1}, use w=1^m: a 0 in B gives the odd-subset cycle, while B all 1s followed by a 0 or 2 in A gives the even-subset cycle with one inserted vertex. In the remaining neutral case, T(wA) is the order join of the one-apex tournament on x and w with the one-apex tournament on y and A; inserting initial 1s in the latter adds sources only. Finally, for {0,2}, put w=2^m. An inserted 0 gives y -> odd subset -> b -> y. If B is all 2s, an inserted b gives y -> even subset -> b -> x -> y, using the marked arc x->y. Each construction supplies 2^{m-1} new odd cycles, so the debt is at least 2^m. It cannot equal a fixed b. The checker verifies 3024 bounded cases in addition to these all-order injections.

The closed-source condition has a concrete hostile. For every w in {0,1}*,

    H(w 3 1)=3 H(w 1).

Here A=1 and B=3 lies **outside** the source alphabet. The same order-join decomposition proves the formula: the terminal marked-y triangle has three paths instead of one. This does not iterate within the declared binary source chart, but it demonstrates precisely what changing charts can recover. A next controller should retain this domain transition rather than erase it.

## 3. What survives and the next decisive questions

The first tournament obstruction applies to **every prefix** with a fixed terminal A and insertion B. The binary refinement also excludes nontrivial affine controllers on each fixed two-letter chart when the inserted block remains inside that chart. Neither excludes controllers selected from the actual state, changing source languages, multiple reported continuation coordinates, or non-affine inverse procedures. The next useful candidate must specify its domains and prove that they cover the intended prime and seven-prime targets. Another unrestricted terminal search would repeat a now-closed route.

For square sums, a concrete failed object has become a positive universal construction: remove one edge of the old forced proper cycle, retain its endpoint ports, and use a complete final table. **OPEN:** Does every fixed finite square-sum linear forest embed unchanged in a Hamiltonian cycle of Q_N for every sufficiently large N? A proper prescribed cycle remains an immediate hostile. The present single-path theorem and the inherited many-path theorem do not prove the general statement.

The common exact map is restriction to boundary data. In the square-sum construction, named ports, parity and the complete connected atom word suffice for the claimed embedding. In the tournament construction, passing from an actual word to a positive lattice point loses reachability, and signed identities detect a limitation of affine controllers. These are analogous information boundaries, not an asserted isomorphism between the two graph classes.

## Reproduction

From this packet directory:

```bash
python3 terminal_rigidity.py
python3 -O terminal_rigidity.py
python3 binary_controller_rigidity.py
python3 -O binary_controller_rigidity.py
python3 verify_native_forest.py
python3 -O verify_native_forest.py
cd native_n24_threshold
python3 verify_native_trace.py
python3 -O verify_native_trace.py
c++ -std=c++17 -O3 verify_native_trace.cpp -o /tmp/verify_native_n24
gzip -dc full_incremental_trace.bin.gz > /tmp/native_n24_trace.bin
gzip -dc final_cycle.txt.gz > /tmp/native_n24_final.txt
/tmp/verify_native_n24 /tmp/native_n24_trace.bin /tmp/native_n24_final.txt
```

The universal final-table consumer and both trace consumers use no optimizer. SymPy is used for the exact tournament column closure; OR-Tools is used only by optional discovery. Checks are explicit failures, not removable Python assertions. `MANIFEST.json` hashes this packet and its declared inherited proof inputs. The weaker preliminary short-path table and solver onset probe are not proof dependencies of the stronger theorem.
