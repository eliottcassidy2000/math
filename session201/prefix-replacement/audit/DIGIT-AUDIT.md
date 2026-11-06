# Independent audit of the separate-index doubling stopping rule

**Verdict: PASS in its explicitly formal, coefficient-wise scope.** No defect was found in the sidecar proposition, its index formulas, or its finite verifier. This is a distinct audit addendum. It does not alter or enlarge the main maximal-prefix repair audit.

## Frozen object

The original sidecar inputs are byte-preserved under `frozen-inputs/`:

- DIGIT-ANSATZ.md: `f9a7b0537003b23b6db15423da25f8fd576461ebb58b39c210e90b037fe6fe59`.
- check_digit_ansatz.py: `c4bf47af37fd6839309d9d4def09b7ec7dee01e147c64bdd6e0979a864db56b8`.
- DIGIT-normal.json and DIGIT-optimized.json: both `caeda63b8c707e6d8b4fcbe15fac125f8cca1df5e9acb2bf200938348a46ee1c`.
- The imported source-operator implementation remains the frozen main verifier at `2868922f90f5ba2cc8d6122408d87df692bfd64d1afa5747b814ca4f070af8dc`.

No original input or main-audit source/receipt was changed. No GitHub writes or Git-index changes were made.

## Proof review

### Separate index maps

The O and Z rules follow directly from the definition of g. For even u=2h, the sibling rule is `g(u-1)=g(h)+h-(h mod 2)=g(u)+u/2-(u/2 mod 2)`. For odd u=2h-1, the inverse sibling identity is `g(u+1)=g(u)-h+(h mod 2)`. Restricting all intermediate intended indices to at least two correctly avoids the excluded singleton boundary.

On each original residue class modulo four, O/Z/S has an affine index rule independent of A=g(u), and its g update has coefficient one on A. Composition preserves both facts. Modulo four suffices even for arbitrarily many compositions: the index remains `a*u+b` with integral a,b, and all branch choices require only the current index modulo four. The g correction is affine in u with fixed rational coefficients. The same reasoning applies independently to v.

### The nonempty-word obstruction

With chronological column-vector convention, `M_W` is the left product of the source matrices. The coefficient of AB in the pairing is `(Q M_W)[0,0]=M_W[3,0]`. The first J sends the formal vector e_D to `(1,1,1,1)`.

The parity invariant in the proof is exact. In raw order `(D,E,U,d)`, the full residue set satisfying D=d=1 and U=E modulo two consists of just

- `(1,0,0,1)`, and
- `(1,1,1,1)`.

Each of O, Z, and J maps each of these two states back into this set. Thus d remains odd after every later letter, proving that the AB coefficient is odd for every considered nonempty word. The proof does not require a finite word-length assumption. In fact, this invariant needs only that the first letter is J; the stated narrower class beginning and ending in J is valid.

Separate substitutions `A'=A+c*u+d`, `u'=a*u+b` and their right-index analogues preserve the AB coefficient. There is no other term that can introduce AB: neither new index contains A or B, and neither side's variables depend on the other side's variables. Exchanging the two indices transposes the pairing and leaves its AB coefficient unchanged. Therefore an odd coefficient cannot equal 2^k times an odd coefficient for k>=1.

### Empty cases

The empty expansion is exactly

`A+B+u*v-u-v+p*q-p-q+1`.

Its AB coefficient vanishes. This rules out both empty/nonempty cross cases, since the nonempty AB coefficient is odd and hence nonzero even after multiplication by 2^k. If both words are empty, the individual A and B coefficients remain one under the allowed separate substitutions and exchanges, while the scaled right-hand side has both coefficients 2^k. Adding a constant cannot change any of these coefficients.

These arguments cover all stated cases.

## Independent exact controls

`check_digit_independent.py` does not import the sidecar verifier or the author's arithmetic functions. It uses explicit 4-by-4 integer matrices and an iterative evaluation of g. Its affine controls infer each residue-class map from two independently evaluated index samples, rather than reusing the author's coefficient-update formulas.

The independent checks passed:

- A complete six-transition parity-closure certificate on the two residue states. This verifies the finite induction mechanism used in the all-word proof, rather than treating bounded word enumeration as the proof.
- 9,842 considered nonempty central words through length ten. Matrix multiplication verifies AB coefficient extraction, oddness, and the full parity invariant.
- 15,625 exact empty-pairing expansion controls.
- 181,515 legal O/Z/sibling composition controls through composition length five and starting index 500. Every affine index and g correction agrees with direct integer evaluation.
- The author verifier's exact smaller universe count, 11,947 legal controls, was independently reproduced.
- 18,928 substituted-matrix controls using independently reconstructed affine maps, including exchange and the individual A/B coefficients for the empty case.

The independent normal and optimized receipts are byte-identical. Both normal and optimized replays of the frozen author verifier are byte-identical to the supplied sidecar receipts. All checks use explicit exceptions rather than removable assertions.

## Reproduction

From the common workspace parent:

```bash
python3 session201-prefix-repair-audit/check_digit_independent.py > session201-prefix-repair-audit/digit-independent-normal.json
python3 -O session201-prefix-repair-audit/check_digit_independent.py > session201-prefix-repair-audit/digit-independent-optimized.json
cmp session201-prefix-repair-audit/digit-independent-normal.json session201-prefix-repair-audit/digit-independent-optimized.json
python3 -B session201-prefix-repair-audit/frozen-inputs/check_digit_ansatz.py > session201-prefix-repair-audit/digit-author-replay-normal.json
python3 -B -O session201-prefix-repair-audit/frozen-inputs/check_digit_ansatz.py > session201-prefix-repair-audit/digit-author-replay-optimized.json
cmp session201-prefix-repair-audit/frozen-inputs/DIGIT-normal.json session201-prefix-repair-audit/digit-author-replay-normal.json
cmp session201-prefix-repair-audit/digit-author-replay-normal.json session201-prefix-repair-audit/digit-author-replay-optimized.json
```

All commands and comparisons passed.

## Scope boundary

This result excludes a coefficient-wise polynomial identity after fixing residue classes and treating A,u,B,v as independent formal variables. Actual target coverage still requires A=g(u) and B=g(v). Special relations among the two indices can make formal coefficients cease to be independent; the proof correctly makes no statement excluding functional identities on such subsets.

Accordingly this is a stopping rule for one separate-index coefficient-matching approach. It does not refute coupled-index constructions, target-dependent replacement, finite central-menu occupancy, or universal C1/C2/C3 coverage.
