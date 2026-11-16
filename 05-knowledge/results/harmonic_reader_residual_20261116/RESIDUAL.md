# A sharp signed product residual for binary readers

Let D=2^n and write R=rev_n(b), with J as in [the previous reader formula](../variable_reader_coverage_20261116/PROOF.md). The
one-final-near-sink reader has the exact product form
```
H = I + (Z+g b/D)(2D−R) + g E(b),
E(b) = b R/D − J(b).
```
The remainder is nonnegative, and its sharp upper bound is
```
0 ≤ E(b) ≤ 2k/3 − (2/9)(1−4^(−k)),   k=ceil(n/2).
```
The lower equality holds only for the empty mask. The upper equality
holds only for fully alternating masks with the maximal number of
ones (one phase for odd n; two for even n). Thus the exact product
estimate always undercounts, and the worst signed error term is
asymptotic to gn/3, not an unspecified two-sided O(n).

## Positive run decomposition

Expand bR/D over all ordered pairs of occupied bit positions.
The 1/2 diagonal weight and opposite-direction geometrically
decaying pairs are omitted from the geometric part of J. Subtract
the adjacency surcharge. If b has w ones in λ one-runs, that gives
```
E = λ−w/2 + Σ_(i>j both marked) 2^(−(i−j+1)).
```
For a contiguous run of length ℓ its internal pair sum is
ℓ/2−1+2^(−ℓ). Hence
```
E = Σ_(one-runs) 2^(−ℓ)
      + Σ_(set pairs in different runs) 2^(−(i−j+1)).
```
This already proves exact positivity for every nonempty word.

## A telescoping potential and its exact optimum

Read the bits low to high. Just before a run starts, let x be the
backward density Σ_(j<i) ε_j 2^(−(i−j)) of all previous bits.
It lies in [0,1/2] because the immediately previous bit is zero,
or because the word has just started. Group each one-run with all
following zeros, appending a single dummy final zero if necessary.
The dummy zero changes E by zero, as seen from the positive run
formula. Ignore unused leading zeros. For a group with ℓ ones and
d≥1 zeros, put p=2^(−ℓ), q=2^(−d). Its contribution to E and the
density at the next run are exactly
```
reward = p+(1−p)x,
x' = q(1−p+px).
```
If ℓ=d=1 then
```
reward = (ℓ+d)/3 + (2/3)(x−x').
```
For any other group the right side exceeds the reward by at least 5/12.
To see this, the difference is affine in x. When ℓ≥2 its slope is
−1/3+p(1−2q/3)<0, so its minimum is at x=1/2 and equals
```
(ℓ+d)/3 − 1/6 − p/2 − (2q/3)(1−p/2).
```
The expression subtracted from (ℓ+d)/3 is increasing in p and q, so
the minimum is attained at ℓ=2,d=1 and is 5/12. For ℓ=1,d≥2 the
slope is positive and the minimum at x=0 equals
(1+d)/3−1/2−q/3≥5/12.

With all groups 10 the density iteration is x↦1/4+x/4 starting
at zero, so after k groups x=(1−4^(−k))/3. Telescoping gives the
claimed exact value. To prove optimality, for odd n a nonunit group
would lose at least 5/12 from the generic (n+1)/3 bound, worse than
alternating, which loses less than 2/9. For even n with no dummy zero,
use the same argument with n/3. If n is even and the dummy zero
was needed, the final x is at least 1/4 (the last zero immediately
follows a one), hence the loss is at least 5/12+1/6 from (n+1)/3.
This is at most n/3−1/4, again worse than alternating. If there
are only unit groups, there are at most ceil(n/2) and equality
forces precisely the alternating phases stated above.

## Signed target cut and an exact rounding sublanguage

Put P(b) = I+(Z+gb/D)(2D−R). A target H=N can only be in the
downward product slice
```
N−g E_max(n) ≤ P(b) ≤ N.
```
In particular P>N is a genuine absolute overshoot, independent
of any unknown correction. The rule does not prove that the
slice contains any b. It gives a product-space cut of width
O(log N) at fixed reader scale, rather than a signed uncertainty.

There is also an exact rounding certificate. If gE(b)<2 then
H is the smallest odd integer at least P(b), since H is odd.
For example, if b is one run of ℓ ones then E(b)=2^(−ℓ)
regardless of its position or any zeros around it, so the
certificate applies once 2^ℓ>g/2. More generally, for λ runs
whose one lengths are ≥L and whose internal zero gaps are ≥G,
the initial density of any run after the first is <2^(−G), and
```
E(b) ≤ λ·2^(−L) + (λ−1)·2^(−G).
```
This supplies physically decoded sparse arithmetic families but
does not turn the product inverse into full coverage. In
particular large signed residuals are real: alternating masks
saturate the sharp bound at each length.

## Average correction and concentration

If b is uniform on all n-bit words, the run-start indicator gives
E[λ]=(n+1)/4, while E[w/2]=n/4. Every pair of bits is jointly set
with probability 1/4, so the geometric sum in the decomposition has
mean (n/2−1+2^(−n))/4. Hence
```
mean E = n/8 + 2^(−n−2).
```
Changing a single bit from zero to one changes E by
1/2−c+Σ_(other present bits at distance d)2^(−d−1), where c is
the number of present adjacent positions (0,1,2). With c=0 the
geometric sum is <1/2, with c=1 it lies in [1/4,3/4), and with c=2
it lies in [1/2,1). Thus every one-bit change has magnitude at most
one. The standard bounded-differences estimate gives
P(|E−mean E|≥u)≤2 exp(−2u²/n). This concerns uniform physical
programs, not programs selected by proximity to a prescribed target;
it does not justify substituting the average correction into a decoder.
