# Independent audit of the row-conditioned exponential theorem

**Verdict: PASS — analytic proof accepted; no correction requested.**
The reviewed report is `overnight11_20260906_no3line_rowfreeze.md` in this
directory. I derived the conditional normalization, collision bounds,
finite sum, and exponential constant before reading the complete prose.
The producer source was neither read nor imported.

The inherited tenth theorem supplies the positive-part expansion, the
sharp transposition bound four, and its elementary exponential-moment
lemma. The new source is one uniform column permutation after an arbitrary
deterministic row ordering has been fixed. The retained data are the
prefix/suffix row-induced matching polynomials. The explicit n5 change in
conditional mean is the hostile preventing reuse of the old unconditional
mean. The audit accepts a conditional probability theorem, not an extremal
existence conclusion or an independence statement about line events.

## 1. Conditional normalization and the capacity boundaries

For a diagonal with fixed source row set S, let L=|S| and let m_k(S) count
unordered k-edge matchings in G[S]. Each such matching determines k source
columns and prescribes a distinct target column label for each. Its
containment probability is `1/(n)_k`; ordering the selected cells gives

```
E_sigma(Y)_k = k!m_k(S)/(n)_k,        0<=k<=n.
```

There is no second injection denominator and no additional automorphism
factor. For k>L the matching number is zero; for k>n the moment is directly
zero without evaluating the quotient. The exact conditional mean therefore
uses the induced matching polynomial on each prefix or suffix, with the
full row set appearing once. A full skeleton cycle profile is insufficient
to determine this conditional quantity.

To verify the uniform estimate, choose k ordered distinct rows in S and
one of the two incident edges independently at each row. If q is the
probability that the selected columns are distinct, then
`k!m_k(S)=(L)_k 2^k q`. Let c(S) count columns with both neighbors in S.
The exact collision probability for a selected pair of rows is

```
c(S)/(2L(L-1)) <= 1/(2(L-1)),        L>=2.
```

Counting columns rather than distinct row pairs retains the C4 case,
where two different columns can share the same neighbor pair. The union
bound on k positions and the target-row sampling bound give

```
0 <= (2theta)^k-E_sigma(Y)_k
 <= 2^k k(k-1) theta^(k-1)/n,       theta=L/n.       (1)
```

The simplification uses `theta^k/(L-1)<=2theta^(k-1)/n` only when L>=2.
For k>L the moment is zero and `k(k-1)/2>=L` bounds its deficit already
by the row-sampling term. This covers k>n as well. L=0 is direct and L=1
has only the preceding case for k>=2. Thus no capacity boundary or zero
denominator is hidden in (1). The first moment is exactly 2L/n.

## 2. Uniform mean error and the finite fourth-order bound

The factorial moment bound supplies uniform absolute domination for the
complete positive-part series. Separating its odd and even terms gives

```
-(2/n)lambda^2 cosh(lambda)
 <= E_sigma(Y-2)_+-g(lambda)
 <= (2/n)lambda^2 sinh(lambda),
g(lambda)=lambda-2+(lambda+2)e^-lambda.
```

This is valid for every row set, not only for a random row set. The
diagonal-length trapezoid identity is unchanged by conditioning. Its
integral is `alpha=1-5e^-2`. Integrating the error terms, with n>=4,
reduces the two worst constants to

```
3e^2-9e^-2 <21,
3e^2+13e^-2-4 <20.
```

The latter includes the trapezoid excess. Independent rational Taylor
intervals certify both inequalities. Hence
`alpha*n-21<=mu(rho)<=alpha*n+20` uniformly in G and every fixed rho.

The finite lower bound also checks. With three selected edges, two
different column collisions cannot coexist: either they would require
four edges or all three edges at one column of degree at most two.
Therefore the exact coefficient is

```
m3(S)=(4/3)(L)_3-2c(S)(L-2).
```

Ignoring all column collisions gives `m4(S)<=(2/3)(L)_4`. Apply the
pointwise lower truncation `f(y)>=binom(y,3)-2binom(y,4)`. The bound
c(S)<=L may be substituted with the indicated direction only for L>=2;
the L=1 terms must be kept zero. The producer does this explicitly.
Independent summation gives

```
sum_d E binom(Y_d,3) >= 2(n-3)/3+2/[n(n-1)],
sum_d E binom(Y_d,4) <= 2(2n-3)/15.
```

Thus the finite conditional lower bound is exactly
`2(n-9)/15+2/[n(n-1)]`. Its positivity starts at n=9. The theorem
correctly combines it with zero and the asymptotic lower estimate; it
does not silently square a negative lower bound.

## 3. Concentration with one random permutation

Fixing rho does not alter the column-transposition bound: removing at
most four old points and inserting their replacements changes the excess
by at most four in absolute value. Expose n-1 column labels. Swapping two
available labels pairs the completion sets while preserving every previous
reveal. Conditional expectations consequently lie in an interval of length
at most four at each step.

The inherited elementary exponential-moment argument now has n-1 steps,
so its right side is `exp[2(n-1)s^2]`. Optimization gives

```
P_sigma(X=0 | rho) <= P_sigma(F=0 | rho)
                   <= exp[-mu(rho)^2/(8(n-1))].
```

Together with the uniform mean estimate, this proves the claimed rate
`alpha^2/8` for every fixed row ordering. It doubles the coefficient
supplied by the tenth two-permutation exposure argument. It does not
reduce the sharp universal transposition constant or assume the conditional
mean is invariant under row reorderings.

The scope restrictions are correct: the column permutation must remain
uniform conditional on the selected skeleton and row ordering. One cannot
choose rows after examining the realized columns. Uniform consequences
may be averaged over skeletons or row orders, but concentration around an
arbitrary mixture's own mean is not inferred. The original model may take
the minimum with the inherited tenth finite bound; that older finite bound
is not asserted here for every fixed row ordering. Finally, F=0 equals
selected-direction success, while only X=0=>F=0 holds for the full target.

## 4. Independent exact controls

The referee uses a different rook engine. It enumerates every matching
once by its increasing edge IDs, accumulating its exact occupied-row
mask. A subset zeta sum then gives all row-induced matching polynomials.
This avoids both the producer's row-by-row column-mask dynamic program
and its component matching-polynomial traversal.

The complete small universe is all **6,764 row subsets** of all **29 cycle
profiles** with shore size 2,...,9. Every subset checks the exact m2/m3
identities, the m4 bound, column multiplicities, and every factorial deficit
through order n+4, including zero/full subsets and k>L or k>n. For sizes
through five, literal column permutations also verify every factorial
moment of a canonical target matching on each row subset.

Separately, **296 row orderings** and **29,992 column permutations** cover
every cycle type and every row ordering at n=2,...,5. Their literal defect
means equal the induced-rook formula. All row orderings are mean-invariant
through n4. For C10 at n5, the means 11/10,17/15,7/6 each occur on 40 of
the 120 row orderings. The explicit two orderings in the producer report
give 11/10 and 7/6 as claimed. The C4+C6 histogram is also retained in
the transcript, providing an additional carrier control.

Exact conditional-prefix averages on **383 nodes** check the one-permutation
interval bound independently of the mean calculations. These small trees
support the separate general completion-coupling proof; no observed small
maximum replaces its proved universal bound.

The four larger targets are reconstructed without either small rook engine:
for n96, selected paired rows contribute factors 1+2t or 1+4t+2t^2; for
the n97 cycle, selected cyclic runs of q rows contribute path factors
`sum_k binom(2q+1-k,k)t^k`. The full selected cycle uses its cycle matching
polynomial. This computes the identity and multiplier-five row orderings
independently. All four exact rational means are printed, have floor 30,
and satisfy both positive lower bounds and the uniform envelope.
The producer subsequently compared these four complete rational numbers
against its component-polynomial engine and confirmed exact equality.

## 5. Reproduction and frozen identities

```
python -B 04-computation/overnight11_20260906_no3line_rowfreeze_audit.py
python -B -O 04-computation/overnight11_20260906_no3line_rowfreeze_audit.py
```

Both runs pass **134,492 always-active gates**, with byte-identical LF
output. All numerical gates use integers or exact rational arithmetic.

```
referee source SHA256:
8996356d3ef8010895685ac55de1001f85267e3fc0d216bf90842473442f0178
referee output SHA256:
430d336aa115951c4f2b6e1494b0009d56175e7e3b1a92eefd840766ed834331
semantic SHA256:
af91373bc0528b3de729083023cec06f53c7247c8779e5fc8774e78a79f59f3d
producer source identity supplied for review, not imported:
12e962f404d849ecba8b2acb1db0dd252d276bb7df2d7dec9c5f1d70dfa360f8
producer output identity:
c48ea4c4a0755e08b69eee212390fc7ca75fb6742c60c90bc120ca60a92ca94d
```

The accepted all-size conclusions follow from the analytic conditional
collision bound, summation, and coupling. The finite controls check the
normalization and failure boundaries; no large-size claim is extrapolated
from them. No producer file or repository truth surface was edited.

**Filing:** root integrated these audited artifacts in the eleventh checkpoint;
reproduction commands above are relative to the repository root. Outside-worktree
locations preserve author provenance, not the present reproduction location.
