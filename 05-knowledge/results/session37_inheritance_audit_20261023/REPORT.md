# Session 37: independent inherited-claim audit

Status: COMPLETE. Two false inherited statements are corrected, with independent exact certificates. No Git operation or repository edit was performed. This work audits the two named internal claims only; it does not assess external attached papers or reopen excluded publication topics.

## Results first

1. **THM-506's advertised first failure at eight vertices is false.** The pair of characteristic and permanental polynomials determines the Hamiltonian-path count for every tournament on at most **eight** vertices, not merely seven. This is a finite-exact, independently cross-checked result. At nine vertices an explicit same-polynomials/different-path-count pair exists. Adding sources extends that pair to every order at least nine.
2. **The old rank computation was not itself the problem.** The n=8 pooled rational rank really is 3 before H and 4 after H. Nevertheless every actual equal-coordinate fiber has one H. There are even 59 individual characteristic-polynomial classes in which H is not affine in permanent coordinates, while arbitrary functional recovery still holds.
3. **Bosonic THM-1810's unsigned Vandermonde identity is false already at n=3.** The all-tournament score polynomial is `prod_(i<j)(x_i+x_j)`. The Vandermonde permanent counts only transitive tournaments. At n=3 their difference is exactly `2*x_0*x_1*x_2`, and their values at all ones are eight and six.
4. **“No Wick minus signs, so nothing cancels” is false for polynomial moments.** For `L(s^k)=k!`, `L(s-1)=0`. The strongest relevant detection-depth statement survives through the already-repaired THM-1790, using homogeneous projective existence and eventual nonvanishing, not a signlessness argument.

## Audit scope and correction lineage

Snapshot: `989b11f2ad12def609f2680967fdf90d1b1bc299` in `session35-publication/clean-repo`, observed clean before and after the audit. Exact file targets:

- `01-canon/theorems/THM-506-permanental-companion.md`.
- `01-canon/theorems/THM-1810-the-bosonic-permanent-side-is-why-gmc2-is-hard.md`.

There is also a different `THM-1810-binary-forms-tournaments-subquestions.md`. The correction of its Q1 by `THM-1820-intransitivity-is-schur-convexity-not-covariant-vanishing.md` is unrelated to the bosonic filename; do not apply it as an existing correction of the present error.

Targeted exact-ID, exact-statement, synonym and mistakes-ledger searches found no prior correction of the THM-506 threshold or the unsigned Vandermonde equality in this snapshot. Navigation repeats the wrong threshold in `00-navigation/OPEN-QUESTIONS.md`, around line 6284. The old computational scripts cited by these two files are absent from the entire accessible shared checkout tree. Consequently this audit does not claim to have rerun those historical scripts. It instead supplies new independent scripts and preserves the old rank observation as a separately typed statement.

Existing correction lineage that must be retained: `THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one.md`, repaired 2026-07-24, and `MISTAKE-246`. The original affine-equation-count proof was invalid; the projective homogeneous repair is the correct surviving mechanism.

## 1. The sharp functional-recovery boundary is nine

Let `p(x)=per(xI+A)` and `c(x)=det(xI-A)`. At coverage degree m set `O_m=(p_m-c_m)/2`, `E_m=(p_m+c_m)/2`, where p_m and c_m denote coefficients of x^(n-m). These count cycle packings with odd/even numbers of cycles, respectively. Define D_ab as the number of disjoint cycle pairs of lengths a,b, and T_333 as the number of disjoint directed-triangle triples. Cycles are counted modulo rotation but not reversal.

The established odd-cycle packing identity gives:

- n<=7: `H=1+2(O_3+O_5+O_7)+4E_6` (missing indices are zero).
- n=8: `H=1+2(O_3+O_5+O_7)+4(E_6+D_35)`, with `E_8=D_35+D_44`.
- n=9: `H=1+2(O_3+O_5+O_7+O_9)+4(E_6+D_35)+6T_333`, with `O_9=c_9+T_333`.

The formal n=8 coordinate map merges D_35 and D_44. It does **not** follow that two actual tournaments realize the same retained coordinates and different D_35. The n=8 census establishes the opposite: the realized fibers fix both D_35 and D_44 individually.

### Complete n=8 universe and two independent methods

Every tournament has a Hamiltonian path: inductively insert a new vertex immediately before the first vertex it beats, or at the end if it beats none. Relabel a Hamiltonian path as `0->1->...->7`. Its seven adjacent arcs are then fixed; the other 21 arcs are arbitrary. Thus the `2^21=2,097,152` presentations with that fixed path cover every isomorphism type. No isomorphism algorithm or random-coverage assertion is needed.

The producer's computation (`session37-exclusion/exhaustive_n8.cpp`) uses row-assignment DP for determinant/permanent polynomials and endpoint DP for H. Independent audit `exhaustive_n8_cycle_audit.cpp` uses minimal-vertex-rooted cycle counts, all 651 possible disjoint cycle-pair placements, and the odd-cycle packing formula. Both give exactly **3,217 joint-polynomial fibers and zero unequal-H fibers**. Audit checksum: sum of H over presentations is 815,136,768; minimum H=1, maximum H=661. All counters fit signed 32 bits except the presentation sum, which is signed 64 bit.

The audit keys `(c3,c4,c5,c6,c7,c8,D33,D34,D35+D44)` are exactly equivalent to `(p,c)` at n=8, because no three cycles can fit and the signed/unsigned coefficients separate single-cycle and two-cycle counts. The TSV contains one representative for each of the 3,217 keys. The universe is presentations, not 3,217 isomorphism classes.

### A three-tournament witness to affine failure without functional failure

Upper-pair-lex edge masks use bit 1 for `i->j` when i<j. The masks

`172254133`, `172256153`, `172238829`

all have characteristic polynomial coefficients, in ascending powers of x,

`[-2,-9,-20,-28,-22,-12,0,0,1]`.

Their permanent polynomials have low three coefficients `(12,31,40)`, `(16,39,44)`, `(18,43,46)` and the same remaining coefficients `[28,22,12,0,0,1]`. Their H values are 177, 197, 209. The two permanent difference vectors are `(4,8,4)` and `(6,12,6)`, in ratio 2:3, but H differences 20 and 32 are not in that ratio. Hence no affine function recovers H even on these three points. Their permanent coordinates are distinct, so there is no obstruction to arbitrary functional recovery.

All three values were independently checked by direct enumeration of all 8! permutations, both as path orders and as cycle covers. The full n=8 rank audit yields 1,523 characteristic-polynomial classes, pooled rank 3->4, and 59 classes with their own affine-rank increase. The original rank table's n=8 row therefore survives only as an affine obstruction.

More generally, on the abstract finite carrier set `(x,y)=(0,0),(1,1),(2,4)`, both coordinates are linear projections of the ambient plane and `rank[Dx|Dy]>rank[Dx]`, yet y is the function x² of x. The earliest failed implication is replacing a function on a finite realized image by an affine map on its linear span. Pooling classes can impose the further unnecessary restriction of a common slope across distinct characteristic classes.

### Explicit n=9 collision

Masks `29759988540` and `18875161993` share:

- `per(xI+A)` ascending: `[384,692,629,380,154,68,28,0,0,1]`.
- `det(xI-A)` ascending: `[-20,-82,-167,-192,-154,-68,-28,0,0,1]`.

Their H counts are respectively **2843** and **2841**. Independently enumerating all 9! permutations for each matrix confirms the two polynomial arrays and both path counts. Full matrices and packing counts are in `verify_tournament_witness.json`.

The decisive sidecars are, first/second:

- D_35 = 209/210, D_44 = 96/95, so E_8=305 on both.
- T_333 = 11/10, c_9 = 191/192, so O_9=202 on both.
- D_33=94 on both, and c_3=28,c_5=154,c_7=398 on both.

Thus the H difference is `4*(-1)+6*(1)=2`. This is a genuine realized collision, not a kernel vector of an ambient coordinate map.

### Extension and minimality

Adjoining a new source multiplies both p and c by x and leaves H unchanged: the source must be first in any Hamiltonian path, and all cycle covers keep the source fixed. Repeating gives collisions at every n>=9. Conversely any collision at a smaller order would extend to n=8, contradicting the exhaustive result. Therefore the universal recovery property is true exactly for n<=8. The n<=7 portion also has the direct formula above; the n=8 step is computer-assisted/finite-exact.

## 2. The two unsigned tournament polynomials are different

Use `V_ij=x_i^(n-1-j)` for 0<=i,j<n so `det V=prod_(i<j)(x_i-x_j)` without a global sign ambiguity. If ascending exponent columns are used instead, multiply the signed identity by `(-1)^(n choose 2)`; the permanent is unchanged.

Choosing one endpoint variable from each factor proves

`U_n(x)=prod_(i<j)(x_i+x_j)=sum_(all tournaments T) x^score(T)`.

Expanding the permanent proves

`per V=sum_(permutations pi) prod_i x_i^(n-1-pi(i))`.

The exponent vectors here are exactly the permutations of `(0,...,n-1)`. A tournament with such distinct scores is transitive: remove the score-(n-1) winner and induct. Conversely each total order gives exactly one such score vector. Therefore

`per V=sum_(transitive tournaments T) x^score(T)`.

Consequently `U_n-per V` is the nonnegative-coefficient generating polynomial of intransitive tournaments. It vanishes identically exactly for n<=2. At n=3 it is `2*x_0*x_1*x_2`. At all ones, U_n is `2^(n choose 2)` and per V is n!. For n>=3 the former is strictly larger: the base case is 8>6 and the successive factors satisfy 2^n>n+1.

The signed tournament identity survives. The strongest useful comparison is that the signed product cancels nontransitive score fibers whereas the unsigned edge product retains them; replacing det by per in the *Vandermonde matrix* is a different unsigned lift. All three objects should stay distinct.

Independent exact coefficient dictionaries verify these identities for n=1..5, with all-tournament totals 1,2,8,64,1024 and permanent totals 1,2,6,24,120.

## 3. What survives of the Gaussian interpretation

For a normalized circular complex Gaussian Z, `E[Z^a conjugate(Z)^b]=0` unless a=b, and `E[|Z|^(2a)]=a!=per(J_a)`. Also det(J_a)=0 for a>=2. These algebraic moment/repeated-state facts survive. They do not identify the all-tournament score polynomial with per V.

All-positive *combinatorial* Wick signs do not prevent coefficient or phase cancellation: `L(s-1)=1!-0!=0`, whereas `L((s-1)^2)=1`. Already-corrected THM-1790 gives the stronger two-charge witness `H=s(s-2)`, `P=W+Z(ZW-2)`: its first three Gaussian moments vanish and its fourth is 48. The degree-uniform lower-bound mechanism is projective homogeneous equation solving plus eventual nonvanishing.

No reduction from permanent computation to GMC(2), or converse, is provided in the inherited note. Hence permanent hardness cannot establish GMC(2) hardness, a necessary proof method, or impossibility of all fermionic methods. Constant-term projection is not automatically a fermionic determinant functional. Keep these analogies explicitly heuristic, and route detection-depth claims to their precise proved statements rather than to a sign slogan. This audit makes no new claim about the current status of GMC(2).

## Files and reproduction

Run from this directory:

```
python3 audit_polynomial_identities.py
python3 verify_tournament_witness.py
python3 verify_tournament_witness.py --n 8 --masks 172254133 172256153 172238829 --output n8_affine_obstruction_direct.json
g++ -O3 -std=c++17 -Wall -Wextra -pedantic exhaustive_n8_cycle_audit.cpp -o /tmp/session37-n8-audit
/tmp/session37-n8-audit 2097152 n8_fibers.tsv
python3 audit_n8_ranks.py
```

`proposed-canonical/` contains proposed complete replacements and a unified diff, outside the repository. `MISTAKES-proposed-entries.md` preserves the first failed implications, minimal witnesses, strongest survivors and repair lineage; the integrating editor should allocate final ledger IDs. No historical data or unrelated work is deleted by this audit.
