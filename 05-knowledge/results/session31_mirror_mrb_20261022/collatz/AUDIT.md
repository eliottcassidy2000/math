# Audit of the supplied Collatz manuscript

Session31, 2026-10-22. This is an audit of Motohiro Hamaguchi's supplied 31-page manuscript, *A Complete Proof of the Collatz Conjecture: Deterministic Liquidation Mechanism via Binary State Transitions and Hierarchical Collapse*, dated May 23, 2026, including its two-page appendix. Page references below are the printed PDF pages, which agree with the PDF page numbers. The manuscript has section labels but no numbered theorem/lemma system.

**Verdict: the proposed proof is invalid.** The first explicit false algebraic transition is in section 2.3, page 9. Independent failures remain after correcting it: the induction loses its numerical baseline, the asserted run/cascade contraction is false, and the appendix's proposed post-cascade fuel exhaustion is false. These are counterexamples to the manuscript's intermediate claims, **not counterexamples to the Collatz conjecture**.

**Status ledger.** REFUTED: the indicated manuscript claims. PROVED, narrowly scoped: the exact local-run formulas and adversarial finite-word constructions below. FINITE-EXACT: the explicitly bounded checks and trajectory census in `receipt.json`. OPEN: universal Collatz convergence, the repository's uniform native 153N gate, and arbitrary-source entry. No literature-priority claim is made for these elementary calculations.

## 1. Sources, inheritance, and audit boundaries

The primary source is the supplied PDF with SHA-256 `6b5c2408391e7880efef17de205ed3b5bda9a001b533f6558213cdd076829ae4`. Its supplied extracted text was checked against rendered pages 9, 17, and 31 for the decisive passages. The audit uses the ordinary map C(n)=n/2 for even n and C(n)=3n+1 for odd n; it never silently substitutes a different Collatz convention.

The inheritance pass used the repository operating contract, START-HERE, relevant CURRENT-FRONTIER entry, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL, and the following targeted sources:

- `05-knowledge/results/session18_native_annulus/README.md`, especially the exact finite-horizon bill and source/clock/terminal requirements.
- `05-knowledge/results/scale_credit_supply_20261021a/README.md`, especially sections 4-5: source-disjoint linear relation supply is not timely arrival or a uniform credit budget.
- `05-knowledge/results/tail_geometry_20261021c/README.md`: linked-source future bounds, including the exterior and arrival-boundary obligations.
- The Session30 local `session30-collatz/REPORT.md` and `CONTRACTION_EXTENSION.md`: completed-bank contraction is not contraction of the entire native energy.

Closest proved mechanism: an exact finite branch/arrival identity. Canonical hostile: a large local decrease can end above the original induction threshold. Corrected near miss: an algebraically exact peak formula does not bound the peak's 2-adic valuation. Least-used sidecar: the post-cascade odd quotient, which can already contain an arbitrarily long trailing-one run. Concept board: exact map / induction baseline / run length / collapse valuation / post-collapse quotient / source-owned debt.

The bounded read-only startup reported a clean cached checkout at `9860c2977d36`. The targeted message peek could not run because `.machine-id` was absent; it requested registration. No registration, repository mutation, Git mutation, or broad sparse Git search was performed. The audit did not need that queue to decide the explicit arithmetic claims.

## 2. First false transformation: the shifted valuation is not a halving coordinate

### Manuscript claim: section 2.3, page 9

After writing N=a·2^n−2 with a odd, the manuscript states that repeated division of N by 2 decreases n by one at every step and leads to 2a−2, then a−1.

In fact,

    (a·2^n−2)/2 = a·2^(n−1)−1,

not a·2^(n−1)−2. For n≥2, the actual result is odd, so it leaves the even-number grid immediately. A small nondegenerate witness is

    10 = 3·2²−2;   C(10)=5, whereas 3·2¹−2=4.

The grid in section 2.1 is a valid parametrization: for every positive even N,

    n = v₂(N+2),    a = (N+2)/2^n

is unique and a is odd. The pair a=1,n=1 produces 0 and must be excluded if the domain is positive even integers. The parametrization does not make n a Lyapunov function for C. For example,

    12 = 7·2¹−2  ->  6 = 1·2³−2,

so a single ordinary halving can increase n from 1 to 3.

The qualitative inevitable-collapse assertions in section 2.2, pages 8-9, are already unsupported. Section 2.3 is the earliest displayed dynamical rule that is directly false, rather than merely announced without proof.

## 3. The induction uses the right principle with the wrong landing condition

### Sections 3.1-3.4, pages 10-13

The identity P/2=Q for P=a·2^n−2 and Q=a·2^(n−1)−1 is correct. Q is odd only for n≥2, not for n=1. This is the ordinary fact that an even integer follows the trajectory of its half. It proves no new descent for an odd starting integer.

The valid strong-induction statement is:

> If all positive integers below K are known to converge and an orbit reaches an actual integer less than K, then that orbit converges, regardless of later excursions.

What is missing is a proof of that actual numerical inequality. A lower layer index does not suffice when the odd core changes: 100 has layer 1 and is larger than 6, whose layer is 3. Likewise, landing in one of the residue classes 0, 2, or 4 modulo 8 does not by itself place the value below an earlier reference threshold.

A concrete baseline test starts at 27:

    27 -> 82 -> 41 -> 124 -> 62 -> 31.

82 is 2 modulo 8 and 124 is 4 modulo 8. Even after their advertised local reductions, the displayed values are still above the original 27. This is not a counterexample to eventual descent; it pinpoints the missing implication from a local class label to descent below the fixed induction baseline.

The positive local survivors are elementary: an 8m term permits at least three successive halvings, an 8m+4 term permits exactly two before the next odd term, and an odd x≡1 mod4 with x>1 satisfies C³(x)=(3x+1)/4<x. Their validity is relative to the **current** value, not an arbitrary earlier K. Section 3.4's promise of later descent refers forward to precisely the unproved and false contraction mechanism audited below.

## 4. The peak formula survives, but its use as a global envelope does not

### Exact repaired local formula

For odd x=a·2^r−1 with a odd and r=v₂(x+1)≥1, define the half-step T(x)=(3x+1)/2 while x is odd. Then

    T^j(x)=a·3^j·2^(r−j)−1,    0≤j≤r.

The values for j<r are odd, and the value for j=r is even. Therefore the raw peak at the end of this first increasing odd run is

    P = 2(a·3^r−1).

Writing s=v₂(a·3^r−1)≥1, the complete following halving cascade ends at

    y = (a·3^r−1)/2^s,

an odd integer. The raw peak has valuation s+1. The run length r and cascade length s+1 are different coordinates.

For the manuscript's even N=a·2^n−2 with n≥2, take x=N/2 and r=n−1. This proves the displayed formula in section 4.1, page 14:

    P_n = 2(a·3^(n−1)−1).

It is the first peak divisible by 4 after the initial rising run, not necessarily the global maximum of the entire trajectory. At n=1, the formula simply returns the already-divisible-by-four initial N.

### The strict envelope has omitted boundary cases

For n≥2,

    P_(n−1) − P_n/4 = (a·3^(n−2)−3)/2.

Thus P_n/4<P_(n−1) holds exactly when a·3^(n−2)>3. It is equality for N=6 (a=1,n=3) and N=10 (a=3,n=2). For N=2 it fails, with the purported lower peak equal to 0, outside the positive domain. Even where the inequality is strict, comparison to a hypothetical other layer's peak does not exhibit a trajectory landing below the original K.

### A second peak label is used inconsistently

Section 4.2, pages 15-16, correctly computes two odd updates from 8R+6 to 18R+16. It then calls this the peak P_n. For R=1, N=14 reaches 34 after those two updates, but the first peak divisible by 4 is 52. The manuscript itself subsequently finds that 18R+16 may be 2 or 6 modulo 8, inconsistent with the earlier P_n, which really is always divisible by 4. The branch arithmetic survives; the identification of the intermediate number with the completed-run peak does not.

## 5. Long runs do not force compensating cascades

### Sections 4.3-4.4, pages 17-18

At the elementary binary level, x≡15 mod16 implies

    3x+1 ≡ 46 ≡14 mod16,

so this immediate result has exactly one trailing zero, not a long string of zeros. Also, 3x+1 is not multiplication by 1.5 followed by addition of 1, as stated on page 17; even the combined half-step is 1.5x+0.5.

More importantly, the claimed compensation fails even after following the entire increasing run and complete following cascade. Let

    N=2^n−2,    n≡3 mod4.

The run has n−1 odd updates. Since n−1≡2 mod4,

    3^(n−1)≡9 mod16,

and therefore v₂(3^(n−1)−1)=3 exactly. Every corresponding first peak has **exactly four** factors of 2, although the rising run can be arbitrarily long. The cascade ends at

    y=(3^(n−1)−1)/8.

Its ratio to N is unbounded as n increases through this progression. Thus no run-length-funded compensating cascade is available.

Explicit hostile with a genuine 8m cascade:

    N=2046,   P=118096,   P/16=7381>2046.

The four-step cascade is

    118096 -> 59048 -> 29524 -> 14762 -> 7381.

This refutes even net contraction for the natural cycle consisting of one whole increasing run and its following complete cascade. The weaker example N=14 ends its first such cycle at 13, which is greater than (3/4)·14=10.5. It directly contradicts the displayed section 4.4 bound at i=1 under that cycle definition.

The manuscript never supplies a different precise cycle definition. If a “cycle” is instead chosen to end only once the required contraction has occurred, the inequality is a stopping condition whose universal attainment still requires proof, not a descent theorem.

On a fixed run/valuation branch the endpoint is

    y = [3^r x + (3^r−2^r)] / 2^(r+s).

This is an injective affine function of its admissible starting integers. The upper remainder has not been erased: it remains in the nonzero coefficient of x. Branch changes can merge different integers, but do not justify a uniform reset or a decrease across every branch.

## 6. The appendix's proposed fuel exhaustion is false

### Appendix steps 1-2 and conclusion, pages 30-31

The appendix asserts that after a large cascade the remaining integer is depleted of trailing ones, cannot immediately rebuild a strong ascent, and must fall below its original value during a cooling period.

For any odd r≥3 and even k≥2, let

    y=2^r−1,     x=(2^k y−1)/3.

Because y≡1 mod3 and 2^k≡1 mod3, x is a positive odd integer. Exactly

    3x+1=2^k y.

The k-fold cascade lands immediately at the all-ones integer y. Both k and r are arbitrarily large. For example,

    2645 -> 7936 -> 3968 -> 1984 -> 992 -> 496
         -> 248 -> 124 -> 62 -> 31,

and 31=11111₂ immediately starts its long increasing odd run. There is no forced post-cascade absence of trailing ones. A smaller witness is 37 ->112 ->56 ->28 ->14 ->7, followed immediately by 7 ->22 ->11 ->34 ->17 ->52.

### Stronger survivor: all three local parameters can be prescribed

For arbitrary positive integers r,s,t, choose a positive odd b solving

    2^(s+t)b ≡ 2^s−1  (mod 3^r).

This is possible because 2 is invertible modulo 3^r; adding 3^r switches parity, and adding 2·3^r preserves it. Put

    a=[2^(s+t)b−2^s+1]/3^r,
    x=a·2^r−1,
    y=2^t b−1.

The integer a is positive and odd. Then

    a·3^r−1 = 2^s y,

and y is odd with v₂(y+1)=t exactly. The first increasing odd run has length r, its final raw peak has valuation s+1, and its post-cascade odd endpoint immediately has a new run of length t. Thus a long prior ascent, a large collapse, and an immediate long new ascent are arithmetically compatible, for any prescribed finite sizes. Example r=s=t=20 is recorded and literally checked in `receipt.json`.

This theorem makes a finite-word obstruction precise. It does not assert a divergent orbit, an infinite prescribed valuation sequence represented by a positive integer, or uniform density of such words as the parameters grow.

### The “approximately 5 to 6 rebounds” boundary cannot be a uniform finite-word bound

The manuscript gives no exact definition or specified census for that empirical bound. Under the natural run/cascade definition above, arbitrary finite strings of net-growing cycles exist. For any m≥1 let

    x₀=2^(3m+1)−5,
    x_j=9^j·2^(3(m−j)+1)−5,  0≤j≤m.

For j<m, x_j≡11 mod16. Its odd-step valuation word is (1,2), and its complete cycle map is

    x_(j+1)=(9x_j+5)/8>x_j.

All intermediate raw values in that cycle are at least x_j. Hence the first m complete cycles stay above x₀ and strictly increase their endpoints. Taking m=7 already rules out a universal six-cycle limit in this precise meaning, and m is arbitrary. This finite construction does not rule out eventual descent after the constructed prefix.

### Guarded hostile: unit conditions modulo 15 and 35 do not bound this prefix

Read the two guards directly on the ordinary starting integer x₀=2^(3m+1)−5. Since 2¹²≡1 mod105, the residue pair depends only on m modulo 4. The exact table is:

| m modulo 4 | x₀ modulo 15 | x₀ modulo 35 | x₀ modulo 105 | Unit modulo both 15 and 35? |
|---|---:|---:|---:|---|
| 0 | 12 | 32 | 102 | No |
| 1 | 11 | 11 | 11 | Yes |
| 2 | 3 | 18 | 18 | No |
| 3 | 14 | 4 | 74 | Yes |

Here m≥1; the row m≡0 applies beginning at m=4, not the negative initial value produced by m=0. The joint modulus is lcm(15,35)=105, and the two residues agree modulo their shared factor 5. Thus

    gcd(x₀,15·35)=1  if and only if  m is odd.

Odd m is unbounded. Consequently, even the conjunction of these two starting-value unit guards permits arbitrarily long prefixes of complete run/cascade cycles with increasing endpoints. In particular, no constant depending only on these two pass/fail guard outcomes bounds the number of such cycles before the orbit falls below x₀. The same obstruction holds within either passing residue pair separately, by taking m≡1 or m≡3 modulo 4.

This is a statement about the guards on **x₀ itself**. No identification with the repository's source-address registers or guards has been supplied, and none is claimed. Those coordinates need their actual map, full native clocks, and ownership data; the table provides neither native 153N coverage nor a counterexample to it. It also does not prescribe the residues of every later endpoint. The checker verifies the table directly for 1≤m≤512, including the product/lcm unit equivalence, while the preceding prefix replay checks m≤40; the displayed periodicity and prefix formulas prove the unbounded statement.

## 7. Base-state and final-gate claims

The classifications on page 19 are correct: n=1 means N≡0 mod4, and n=2 means N≡2 mod8. But their union is not forward closed. Already 12 ->6 changes n=1 to n=3. Another segment is 82 ->41 ->124 ->62, with the even layers 2,1,6. No irreversible odd-core or bit-length descent is established by calling the first two layers “base states.”

The claims on pages 20-21 that every non-power-of-two start must pass through 3 or 5 are false:

    21 ->64 ->32 ->16 ->8 ->4 ->2 ->1.

More generally, (2^(2j)−1)/3 for j≥3 goes directly to a power of two and avoids both proposed gates. Conditional on convergence, every non-power-of-two starting integer does pass through 16: its last odd integer above 1 is (2^k−1)/3 with even k≥4. This conditional observation is not a proof of convergence and does not imply the false 3-or-5 gate.

## 8. Finite empirical claims: preserve the data, correct the accounting

### The showcase trajectory, section 5.1, pages 22-24

The named start 14693280766 really reaches 1 in 445 ordinary C-steps. Its actual maximum is 564859072960, reached after 18 steps. The text's 564859000000 is a rounded value, not the exact integer.

The maximum has v₂=6. Its complete descending block is

    564859072960 ->282429536480 ->141214768240 ->70607384120
      ->35303692060 ->17651846030 ->8825923015.

The first four terms are each divisible by 8. Those labels describe **overlapping** three-halving windows. They do not authorize adding four separate three-bit contractions. The passage has six halvings in total, not the at-least-twelve claimed on page 23. This is a literal double-counting of operations, not a disagreement about metaphors.

### The below-10000 census, sections 5.2-5.3, page 25

The checker exhaustively follows every integer 1≤n≤10000 to 1. In this exact universe:

- 3997 orbits pass through 1822: 39.97%, not approximately 16%.
- 4011 pass through 9232.
- 1579 have their global trajectory maximum equal to 9232: 15.79%.

The counts are unchanged for 1≤n<10000, so the endpoint convention does not explain the discrepancy. The near-16% number is numerically consistent with the maximum-9232 statistic, but this audit does not claim to know why the manuscript reports it for the 1822 hub.

These are finite observations, not evidence that the paper's unproved uniform inequalities hold. The suggested 9232 “energy concentration” mechanism is not supplied as a mathematical derivation. Also, 911 and 1367 mentioned as backward routes from 1822 are successive forward odd stages: 1822 ->911 ->2734 ->1367. The underlying inverse integrality restriction modulo 3 is real, but does not imply a universal catchment mechanism.

## 9. Precisely what transfers to the source-owned debt task

There is a useful methodological connection, not a new native-energy theorem.

1. **Source:** actual ordinary Collatz trajectory segments. **Target:** a run/cascade word. **Map:** retain its exact odd-update count, every halving, and endpoint. **Preserved:** literal orbit edges and total arithmetic multiplier. **Lost by residue labels alone:** original baseline, upper quotient, total distinct operations, and continuation. **Restoring sidecar:** complete exponent word and actual endpoints. **Cheapest test:** 2046 ->118096 ->7381 and the overlapping showcase cascade.

2. The showcase's four overlapping 8m labels cannot each spend their own three-halving reduction. This is the same accounting defect, at the level of literal operations, that the native source ledger forbids for reused source norms or duplicated pair credits. The correspondence is a bookkeeping analogy only: numerical decrease of an ordinary orbit is not the Eisenstein energy Q of a source bank.

3. The Session30 completed-bank theorem uses original primitive sources, disjoint ownership, actual common terminals, full native P clocks, and an exterior-plus-current-boundary term. None of those quantitative ingredients is bounded by the supplied manuscript. Its ordinary odd-run formulas are not a conjugacy to the native P map, which also strips B-braids. They do not authenticate new source-owned relations or pay the missing uniform 153N remainder.

4. The new three-parameter construction supplies hostile finite words for proposed local potential or reset arguments. Any future native transfer must explicitly compile each word to native P edges, preserve full clocks and source identities, and prove the relevant source measure and exterior/arrival bounds. A family existing at arbitrary height is not automatically a linear-density source bank in each annulus.

The strongest immediate conclusion for the uniform relation task is therefore negative and sharply scoped: the manuscript supplies no missing arrival or debt lemma. The existing exact bank identities remain intact, and their previously stated uniform exterior-plus-boundary obligation remains open.

## 10. Reproduction and stopping boundary

Run from this directory:

    python3 verify_audit.py --out receipt.json
    python3 -O verify_audit.py --out receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The two receipts agree byte-for-byte. The checker records 20,870 explicit checks: the unique grid; 1,581 literal first-run peak cases; 1,581 envelope identities; 20 fixed-four-halving hostile family cases; 400 inverse-recharge cases; 1,728 small three-parameter constructions and one (20,20,20) case; 40 repeated-growth-prefix constructions; 512 direct starting-value guard cases and their modular periodicity; the decisive small witnesses; the full named large trajectory; and all 10,000 starts in the finite census. Validation uses explicit exceptions, not Python assertions, and integer arithmetic only. Sequence hashes pin the finite data.

The all-parameter results rest on the displayed proofs, not on these bounded checks. The proof audit is complete at the level needed to reject the claimed complete proof. No search for a Collatz counterexample, no claim of universal convergence, and no unsupported replacement Lyapunov theorem is made.
