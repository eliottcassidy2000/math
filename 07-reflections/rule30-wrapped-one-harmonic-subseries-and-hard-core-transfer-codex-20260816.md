# Wrapped Rule-30 ones form a convergent harmonic subseries

**Status: COROLLARY OF PROVED
[THM-3493](../01-canon/theorems/THM-3493-rule30-dyadic-wrap-atlas.md).**
The argument below is exact from its independently audited dyadic wrap-prefix
atlas.  This reflection is not a separate canon promotion or a Rule 30 prize
claim.

## 1. The marked subset of the harmonic series

Let `c_k` be the distinguished Rule 30 center bit.  Split the positive depths
into THM-3493's wrapped and hard sets `W,H`, and mark their one-supports:

```text
C_1={k>=1:c_k=1},
W_1=C_1 intersect W,
H_1=C_1 intersect H.                                  (1)
```

For the dyadic block

```text
B_m={2^m,...,2^(m+1)-1},                              (2)
```

THM-3493 says either the entire block is hard or its wrapped part is an
initial prefix whose center word is

```text
0...01.                                                (3)
```

Therefore

```text
|W_1 intersect B_m|<=1                                (4)
```

for every `m>=0`.  Any marked depth in that block is at least `2^m`, so

```text
sum_(k in W_1) 1/k
 <=sum_(m>=0) 2^(-m)
 =2.                                                   (5)
```

Thus `W_1` is a convergent subset of the harmonic series.  More quantitatively,
for `M>=0`,

```text
sum_(k in W_1, k>=2^M) 1/k <=2^(1-M).                 (6)
```

The estimate is uniform in the unknown innovation valuations.  It uses only
the one-mark-per-dyadic-block theorem.

## 2. Harmonic divergence transfers to the hard core

Write

```text
S_A(N)=sum_(k in A, k<=N)1/k.                         (7)
```

Since `C_1=W_1 disjoint_union H_1`, equation (5) gives

```text
0<=S_(C_1)(N)-S_(H_1)(N)<=2                           (8)
```

for every `N`.  Consequently,

```text
sum_(k in C_1)1/k diverges
 iff sum_(k in H_1)1/k diverges.                      (9)
```

This is stronger than saying wrapped ones have zero natural density: the
entire reciprocal mass they can ever contribute is bounded.

If the center-one counting function satisfies

```text
|C_1 intersect [1,N]|=delta N+o(N),                  (10)
```

partial summation gives

```text
S_(C_1)(N)=delta log N+o(log N).                      (11)
```

Equation (8) then transfers the same logarithmic main term to the hard ones:

```text
S_(H_1)(N)=delta log N+o(log N).                      (12)
```

In particular, the still-open center-balance conjecture would force

```text
S_(H_1)(N)=(1/2)log N+o(log N),                       (13)
```

while the wrapped correction remains bounded by two.

## 3. Recurrence meaning and the exact boundary

In a nonempty wrap block, its unique marked wrapped depth is

```text
v_m=nu_2(R_(2^m)-1)=kappa_(m+1),                     (14)
```

the next period-doubling innovation.  Hence `W_1` is not an arbitrary sparse
set: it is the visible subspine of the innovation recurrence, with at most one
selected node at each dyadic scale.  The harmonic weight turns that scale
separation into summability.

The unmarked wrapped set `W` is different.  A wrap prefix may occupy a large
fraction of its block, and

```text
sum_(k in W)1/k                                      (15)
```

need not converge from the current theorem.  The full-prefix hostile would
contribute nearly `log 2` in every dyadic block.  It is the terminal `1` in
the word `(3)`, not wrappedness alone, that creates `(5)`.

Likewise, `(5)` does not prove center nonperiodicity or balance.  It says that
any unbounded harmonic signal in the center support must live on the hard
pointed-terminal depths.  Those are precisely the depths for which the
staircase/phase machinery still carries its large coordinate burden.

## 4. What “every subset is a subset of the harmonic series” preserves

The map

```text
A subset N  |->  sum_(n in A)1/n                    (16)
```

is useful only after the native recurrence and scale partition are retained.
For `W_1`, it preserves the dyadic sparsity strongly enough to prove
convergence.  It forgets the actual innovation labels, center word, phase
profile, and ancestry of the selected depths.  Two very different recurrence
subsets can therefore have the same convergence status or even the same sum.

No intrinsic pairwise orientation appears in `(16)`, so forcing these subsets
into a tournament would add a gauge rather than information.  The faithful
object is the graded subset together with its dyadic-block address and marked
recurrence role.
