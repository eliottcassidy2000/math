# A one-arc unit shield with an arbitrary binary head

This is an exact two-apex seam different from the previous dyadic
fixed-slope switches. It both deletes to a factor-1 transitive spine
and reads a signed factor equal to 1 from a growing strong guard.
The head is completely arbitrary (subject to a leading and trailing 1);
the exact length of a repeated tail is chosen from its signed residual.

Use the [reverse bridge and word coordinates](SHIELD.md). Take the
two short left words

    u0 = 111010             (h0,P0) = (39,59)
    u1 = 111110             (h1,P1) = (33,63).

Exactly one apex–spine arc differs, at zero-based site 3. Their
profile displacement is 2(−3,2). Thus the physical change in H
is 2 if the right strong word v has j=h(v) and Q=A(rc(v))+1
with 2Q−3j=1. Both tournaments are strong and deleting the two
apices leaves the same transitive spine. In fact

    H0 = 85Q−14,                 H1 = 85Q−12.             (1)

## The exact tail language

Temporarily write w=rc(v) and L(w)=2(A(w)+1)−3h(w).
Reading 10 sends the ordinary word state to

    h ← h+2A+2,     A ← 4A+2,      last ← 0,

independently of the previous last bit. Consequently **L←L−2**.
This is a true fixed-step drift, not an approximate change.
If the head has a positive odd L, just repeat the tail the
number of times needed to spend L−1.

Precisely, all strong solutions of L(w)=1, without any width
cutoff, are

    w = x (10)^t,
    x any binary word of length ≥2 beginning and ending in 1,
    t = (2(A(x)+1)−3h(x)−1)/2.                           (2)

For the forward construction, if x ends in 1, its last vertex is a
sink so its h is unchanged from the preceding prefix. The elementary
bound h≤A+1 on that prefix gives h(x)≤(A(x)+1)/2. Thus L(x) is odd
and at least (A(x)+1)/2, so t is an integer ≥1. The word in (2)
has L=1 by the drift identity and is strong.

For the converse, strip the maximal terminal (10)^t from w.
It cannot leave the empty word: L(empty)=−1 and the repeated tail
only decreases it. If the remainder ends in 00, the two zero appends
added at least 3A/4+1 to the count with final address A (including
the positive count of the earlier prefix). Thus
L≤2(A+1)−3(3A/4+1)<0, an impossibility. No remainder ending 10 is
allowed by maximality. So the remainder ends in 1 and begins in 1.
It cannot be the singleton word 1, because L(1)=1 then t=0 and w is
not strong. Its length and t are now exactly (2). Stripping is
a literal, unique inverse, so the free heads never collapse.

For example the first count pairs (ordered by their baseline) are:

| head x | t | guard w before reverse-complement | Q | H0,H1 |
|--|--:|--|--:|--|
| 101 | 1 | 10110 | 23 | 1941,1943 |
| 11 | 2 | 111010 | 59 | 5001,5003 |
| 1011 | 7 | 101110101010101010 | 191147 | 16247481,16247483 |

This is not the analogous false claim that a fixed mixed periodic
word has zero drift. Its nonzero drift is essential: the head is
chosen separately and its signed error fixes the repetition count.

## Sharp growth of the set of anchors

The distinct lower counts in (1) are in bijection with the heads:
distinct heads produce distinct w by the maximal-tail inverse,
the leading 1 makes the numerical word address injective, and
H0 is strictly increasing in Q. This seam has exactly

    (2/3) ln X + O(sqrt(ln X))

such lower counts at most X, with natural logarithm. Counting both
members of each adjacent pair doubles the main term; the pairs
are disjoint since distinct Q differ by at least one.

Let A be the ordinary address of the head, so A runs over every odd
integer >=3. Let n0 be its bit length, D0=2^n0, and let
ν(A)=rev_(n0)(A)/D0 be its binary radical inverse; since A is odd,
1/2≤ν(A)<1. The geometric pair kernel J of the
[signed reversal defect](../../harmonic_reader_residual_20261116/RESIDUAL.md)
satisfies E=Aν−J=O(n0). Summing the elementary recurrence shows
h(x)=1+A−(last bit)−J=A(1−ν)+E, since the head ends 1:
each 1 adds its geometric interactions against the later 0 bits,
with the adjacent surcharge, and before discarding later 1 bits
its total future weight is its address weight, except for the last bit.
Thus the length of w is

    n = n0 + 2t = A(3ν(A)−1) + n0 + 1 − 3E.

Also log₂ H0=n+O(1), since Q is the address of a leading-one word,
plus one. Put T=log₂ X. The condition H0≤X differs from
A(3ν(A)−1)≤T only by an O(log T) boundary band: trivially the condition
or its band forces A=O(T), because t≥(A−1)/4 and 3ν−1≥1/2.

Choose dyadic M≈sqrt(T) and partition the odd A by A mod M.
For a residue a, use v=rev_(log₂ M)(a)/M. Then ν(A) lies in
[v,v+1/M), so replacing A(3ν−1) by A(3v−1) creates uncertainty
O(T/M) on the relevant range. In this one residue class the number
of A in the approximate threshold is T/[M(3v−1)]+O(1), and
the entire uncertainty contributes O(1+(T/M+log T)/M) terms.
The v values run once through the evenly spaced grid in [1/2,1).
The corresponding Riemann sum has error O(T/M), because
v↦1/(3v−1) is Lipschitz on this interval. Summing all errors
gives O(sqrt(T)), and the integral is

    T ∫_(1/2)^1 dv/(3v−1) = T ln(4)/3 = (2/3) ln X.

This shows exactly what the negative drift buys: an arbitrary
binary head can survive in the guard, and one physical arc
closes the residual 2, but the linear drift tariff makes
the baselines logarithmically sparse. It does not close arbitrary
absolute targets.

## Why a different fixed suffix cannot reduce the drift tariff

Take any fixed suffix y of length ℓ, applied to arbitrary zero-ended
heads. As above, A'=2^ℓ A+a and h'=h+C A+c, where C is the
reverse address of the zero sites in y. For the change in L=2(A+1)−3h
to be independent of the head, the coefficient of A must vanish:

    2(2^ℓ−1)−3C=0.

This forces ℓ even and C=2+2³+...+2^(ℓ−1). The actual suffix
is then necessarily (10)^(ℓ/2), and it spends exactly ℓ units of
the residual. Thus another fixed translating suffix on this
same one-apex projector cannot accelerate the correction. More
ports or a changed projector would be needed to compress the
linear tariff.
