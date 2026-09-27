# Independent audit of the all-k return cylinders

**VERIFIED proof audit / FINITE-EXACT independent arithmetic controls.**
Date: 2026-09-27. No mathematical error found in the audited claims of
[the recursive-entry note](entry_20260927_recursive.md). This audit does
not assert universal coverage or completeness of its exploratory policy.

The independent [script](../../04-computation/experiments/entry_20260927_recursive_review.py)
imports no functions or certificate implementation from the note's script.
It uses repeated integer division to compute valuations and actual U steps,
and reads the older bank's explicit frozen rows only for the bank comparison.
The [output](entry_20260927_recursive_review.out) records **1,218,256 explicit
checks** with identical normal and optimized-Python outputs.

## 1. Universal theorem and actual clock

Put A=8^(k+1), B=9^(k+1), for k>=0. Let t>=1 satisfy

```text
2^t(A-5)>B-5,
```

and let b be positive with `bB=5 mod2^t`. Then b is odd. For
`n=bA-5`, the actual first k pairs have valuations `(1,2)` and reach

```text
8b*9^k-5.
```

The following step has valuation1 and reaches `12b*9^k-7`.
Its numerator is exactly `4(bB-5)`, so its next valuation is
`2+v2(bB-5)` and its next endpoint is `oddpart(bB-5)`.
The final two steps thus consume at least t+3 divisions by2; their U
step count is exactly2, not t+3. This resolves the potential clock ambiguity.

The repeated pairs increase the shifted value by9/8 at each pair, and
every intermediate grows from its preceding boundary. The penultimate
step also grows. Thus all first2k+1 actual U iterates exceed n, including
the k=0 case, where the first value is12b-7>8b-5.

For the final comparison write D=2^t A-B. The hypothesis gives
`D>5(2^t-1)>0`, and hence

```text
2^t n-(bB-5)=bD-5(2^t-1)>0.
```

The endpoint is no larger than `(bB-5)/2^t`, proving its strict descent.
Therefore the asserted first-descent time2k+2 is exact. The assumption
works for every admissible b, with no large-source qualifier.

The least residue beta of `5B^(-1) mod2^t` is positive odd. Its source
residue `beta*A-5` is positive and smaller than A*2^t. Consequently every
positive integer in the claimed cylinder is represented with positive b;
there is no unverified fringe. Choosing the least t for the displayed
inequality is an exact sufficient prescription. The note correctly does
not claim optimality after exploiting the minimum allowed beta.

## 2. Hostiles and density comparison

The two weakened-budget witnesses replay exactly:

| k | b | source | endpoint after2k+2 | v2(bB-5) |
|---:|---:|---:|---:|---:|
|5|3|786427|797159|1|
|11|1|68719476731|70607384119|2|

Every earlier iterate also exceeds its source. These are failures of
the weakened selected-return bound, not claims of divergence.

Different k cylinders have distinct exact values `v2(n+5)=3k+3`.
Their disjointness is therefore symbolic, not inferred from a search.
The tail k>=L lies in `n=-5 mod2^(3L+3)`, which justifies the natural
density as the sum of the cylinder measures. Moreover

```text
(B-5)/(A-5)>B/A,
2^t A>B,
1/(A*2^t)<9^(-k-1).
```

Summing the geometric majorant gives the stated strict tail bound
`1/(8*9^L)`. It is not an orbit-frequency assertion.

For the old bank, all171 frozen residue rows disagree with `-5 mod4096`
at their common binary precision. This checks the entire containing
cylinder, so the all-k conclusion for k>=3 is justified even though only
finitely many new rows are generated. The first three new cylinders are
exactly the existing rows at q=7,51,83. Their measures sum to73/1024.
No other overlap is possible, by the preceding containing-cylinder test.

The independent calculation reproduces the first20 added density
`2553380107527241/2^64` and augmented density
`6990313556829423891/2^65`. The note correctly distinguishes first descent
on this set from convergence of every member to1.

## 3. Nested rise in the original fixed-coefficient27 family

For `n_k=4*8^k-5`, the first2k U steps reach `m=4*9^k-5`.
Factoring powers of9 gives

```text
v2(m+1)=2+v2(9^k-1)=5+v2(k).
```

For an odd x with `v2(x+1)=H`, its next H-1 valuations are all1.
They all increase x. Therefore `l=4+v2(k)` more U steps stay above
the original n_k, proving the strict lower bound `tau(n_k)>2k+l`.
The closed endpoint follows from the exact affine word:

```text
3^l(m+1)/2^l-1=3^l(9^k-1)/2^(l-2)-1.
```

For odd k this becomes `(81*9^k-85)/4`, with next exponent
`a_k=v2(243*9^k-251)-2`. Both the k=1 failure and k=3 success in the
note replay: the latter has first descent at step11, ending at691.

For the sparse-gate bound, success requires

```text
2^a_k > (243*9^k-251)/(16*8^k-20)
       > (243/16)*(9/8)^k.
```

The second strict inequality is valid: cross multiplication reduces to
`1215*9^k>1004*8^k`. No unsafe omission of the affine constants is needed.
For two odd indices k<j, subtracting their exact divisibility relations
and using `v2(9^(j-k)-1)=3+v2(j-k)` yields

```text
v2(j-k)>=min(a_k,a_j)-1.
```

If the two a values differ, equality holds here. For two successful gates,
the preceding necessary budgets imply exactly the claimed separation
`j-k>(243/32)(9/8)^k`. This does not prove that only finitely many gates
succeed, and it does not control other later return words. Both limitations
are correctly retained in the note. Reduction modulo32 also verifies
`a_k=2` for k=1 mod4 and the failure of that entire quarter-class gate.

## 4. Certificate types and finite controls

The affine composition `(M,D,C)->(pM,dD,pC+Dc)` and terminal comparison
`(D-M)N>C` are correct with N retained as the immutable source. The
Repeat1 and Repeat12 guards ensure actual valuation words. ToFixed uses
an odd target and a positive power-of-two quotient, so it certifies the
actual next odd integer. A child call preserves its current integer as
its own source and cannot discharge a different parent's target without
the separate terminal comparison.

The all-k family uses at most three macro nodes; k=0 correctly omits
the empty repeat. Literal Steps can represent a known finite descent,
but the note makes no inference that the exploratory generator terminates
for every source. Likewise27 reaching47 in the completed-family control
is correctly rejected until its fixed suffix goes below27.

Independent control universes:

* k=0..256, with nine admissible coefficient lifts each, including a
  lift indexed by `10^25+37`:2,313 actual source trajectories.
* Both weakened-budget hostiles, checking every preceding iterate.
* All171 old-bank rows, the three exact overlap rows, and both densities.
* Fixed27-family k=1..512, checking every growing iterate and nested clock.
* All pairwise comparisons among the256 odd gate indices in that range,
  including the sharper equality when their valuations differ.

Run from the repository root:

```text
python 04-computation/experiments/entry_20260927_recursive_review.py
python -O 04-computation/experiments/entry_20260927_recursive_review.py
```

The universal scope comes from the algebraic proof audit above; the finite
checks independently test the formulas, clocks, guards and equality boundaries.
