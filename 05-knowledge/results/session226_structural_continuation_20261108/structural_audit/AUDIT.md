# Independent audit: Session226 structural extensions

**Verdict: PASS for both submitted packets, 2026-11-08.**

This audit reviewed the proofs and verifiers independently and did not modify
producer files. Exact input hashes are in `audited_inputs.sha256`. The parent
owns publication. No completeness, priority, or unrestricted finite-line claim
is certified here.

## 1. Exact strong grafts: PASS

Reviewed `session226_structural_continuation_20261108/strong_grafts/PROOF.md`
and its verifier, with Session225 graft definitions and formulas.

- **Sufficiency is universal in the gadget.** When the old root reaches every
  old vertex and there is an old nonroot, it reaches that nonroot, then any
  gadget vertex. Every gadget vertex can return to the old root. This also
  gives gadget-to-gadget paths even when the gadget itself is disconnected.
- **Necessity is genuine.** A shortest path from the old root to another old
  vertex is simple. If it enters the gadget, its first return to the old
  graph must repeat the root. Therefore old-root reachability cannot be
  repaired by the gadget. For a singleton source the old root is a sink
  and the graft has at least two vertices, so strongness fails.
- **Iteration quantifiers are sound.** Every nonempty iterate is strong and
  therefore its newly marked root reaches every vertex. Gadget roots, sizes,
  internal strongness, and matrix invertibility are irrelevant.
- **The H9 claim is a finite separator only.** Both five-vertex witnesses
  have nine Hamilton paths, and the five output counts for each are disjoint.
  All ten output graphs are strong. No statement follows about all possible
  source graphs or scalar coverage by other constructions.

The inherited first-birth, nonisomorphism classification and self-converse
claims remain attributed to Session224; this audit does not repeat that full
classification. The separator itself is independently verified without it.

Normal and optimized producer runs agree byte for byte: all 5,625 pairs
pass, with 3,225 strong grafts, including 2,064 using nonstrong gadgets.
The independent script replaces bitset reachability by Boolean Warshall
closure and tests 2,175 pairs with every source root through order three and
all gadgets through order four. It replaces the producer subset path DP by
literal permutation enumeration on each of the ten nine-vertex grafts:

- A: 423, 605, 357, 185, 261.
- B: 315, 541, 551, 217, 207.

## 2. Physical seed and finite-line extension: PASS

Reviewed `session226-root-descent/PROOF.md` and its verifier, with Session225
root-descent Section 1 and graft operation conventions.

### Physical cone and readout

The cut immediately before the root is an injection into ordered covers.
Its image has root-singleton second bin exactly when the original path ended
at the root. Subtracting these matched subsets justifies the stronger
`D-E <= U-d`, not merely `D <= U`. Deleting the terminal root also gives
`E <= d`, including the empty residual graph of a singleton. The implicitly
used `E <= D` follows immediately because E counts a subset of Hamilton
paths. Root sinks and singleton seeds satisfy the nonstrict inequalities.

The endpoint after OJ is `(h,D,h,D)` with `h=D+2U`, `g=2U`, and `h-g=D>=1`.
The scalar `N=D+4U` is the Hamilton count after one **extra** J, not the
unclosed endpoint's D. The packet preserves that distinction. No unjustified
strongness assertion is made for arbitrary seeds.

### Bounded intercept and distinct-value count

D never decreases. Splitting at the **last** non-O letter handles every finite
word without any bound on its length. If the final D is at most C, then just
before this letter `U_0=D_1-D_0+E_0 <= C`, using `E_0 <= D_0`. Both Z and J
leave a positive integer U at most 2C. Trailing O letters then only double U.
The O-only case separately retains the finitely many seed U values, which may
exceed 2C. These cases exhaust all words, including empty words.

Consequently the proof bounds **distinct values**, not words, paths, or graph
orders. Arbitrarily many root-sink Z steps can preserve a state and do not
contradict the count. Finite S and fixed C are essential. Since each positive
base dyadic sequence contributes O(log(X+2)) values, the claimed finite union
has the stated size. Boundary intercepts at most zero are impossible.

### Residues and discarded slopes

Solving the line equation gives `g=(aN-c)/(a+b)`; integer even g is equivalent
to `aN=c mod 2(a+b)`. No illegal division by a modulo a noncoprime modulus is
used. A common even multiple M and an odd residue r give an odd progression
with X/M+O(1) positive members. The residue-gap hypothesis excludes every
retained strict-slope line on that progression. Bounded boundary intercepts
remove at most O(log X) of its members.

For discarded slopes below 2/3 the inequality bounds h; since `0<g<h`, it
also bounds N. Thus they contribute finitely many scalar values even if
infinitely many words represent those values. Negative, zero, and positive
intercepts are handled by the same inequality. For `|c|<a`, the difference
`a-c` lies strictly between zero and `2a<2(a+b)`, making r=1 a valid choice.

The hostile five-line bank correctly covers every odd residue modulo ten.
For each matching residue, `(h,g)=((3N+c)/5,(2N-c)/5)` eventually has
`h>g>0` and `3g>=2h`. It does not establish physical reachability. This
explicit failure of the residue-gap hypothesis prevents an overclaim about
unrestricted finite nonunit line banks.

## 3. Reproducibility and independent controls

Producer runs, normal and optimized, are saved as `strong_normal.json`,
`strong_optimized.json`, `root_normal.json`, and `root_optimized.json`.
Each pair is byte-identical. Root controls cover 5,405 labelled rooted
profiles through order five, 46 distinct profiles, 13,120 words, and 167,895
congruence tests. These are finite controls, not proofs of universal claims.

`independent.py` uses no imports from either producer verifier. In addition
to the strong-graft controls above, it enumerates two-path covers by choosing
the two vertex sets and then their path orders separately. On every rooted
labelled tournament through order four, it independently checks the cone and
855 literal graph O/Z/J transitions. It also checks 812 state/cap cases
against the **exact** finite dyadic envelope, using caps 1,2,3,5,9 and six
transition generations from every distinct seed profile through order four.
It tests all 500 odd N below 1000 against the hostile residue partition and
checks positive cone chambers from N=13 onward.

Run from this audit directory:

```
python3 independent.py > independent.replay.json
python3 -O independent.py > independent.optimized.replay.json
cmp independent.replay.json independent.json
cmp independent.optimized.replay.json independent_optimized.json
```

The saved independent outputs are byte-identical and report PASS. No
assert-only guards are used. The mathematical proofs, not any finite
experiment, establish the universal quantifiers.

## Publication note

The root-descent proof header was promoted from pending audit to PASS, and an
explicit prior endpoint-defect attribution was added after this review. Its
mathematical statements and proofs are unchanged. The input hash ledger
records the exact pre-publication source files reviewed; the package checksum
manifest records the final published bytes.
