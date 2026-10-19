# A matched-block height test beyond expansion powers

19 October 2026. **PROVED extension conditional only on the proved primary
height and overlap lemmas.** Independent audit: PASS.

Use its expansion-first first-crossing language and notation: N is Terras
length, L is expansion-bit length, lambda=3^19/2^30, and B=10^10. All matched
blocks below are actual chronological parity strings, not just equal
multipliers or carry-blind grammar labels.

## 1. Separated matches and the exact return alternative

Suppose identical r-bit parity blocks occur at two distinct macroblock
boundaries i<j of a lawful word W; both occurrences lie wholly in W and may
overlap. Let y_i,y_j be the actual orbit values at their starts. If

    100r>=L,                                     (1)

then either the endpoint strictly descends or y_i=y_j.

To prove this, suppose the endpoint does not descend and y_i!=y_j. Two
sources of the same r-bit word belong to its unique residue modulo 2^r, so

    |y_j-y_i|>=2^r.                              (2)

For every macroblock prefix P before the final first crossing, the primary
carry identity and multiplier envelope give

    y_P=S_P x+e_P
       <=S_P(x+N/6)
       <=lambda^(L/30)(x+N/6).                   (3)

All values are positive. Hence their nonzero difference is strictly less
than their maximum, and (2)--(3) imply

    x+N/6 > 2^r/lambda^(L/30).                   (4)

The primary non-descent height ceiling is

    log(x+N/6)<B+(B+1)log N.

Its proof is even easier than the augmented-height ceiling used there.
Combining (1), (4), log2>1/2, and log lambda<1/12 gives

    log(x+N/6)>L/450.

The same short/long overlap (N<15L; L0=2^48; the exact N<=2^52 lemma) is
incompatible with these two bounds. Thus y_i=y_j is necessary.

Equivalently, a non-descending endpoint with distinct checkpoint values
obeys the quantitative necessary matching budget

    r log2 < (L/30)log lambda+B+(B+1)log N.        (5)

Only the expansion bits before the later start are needed in the first
term, if the prefix-specific bound is retained.

## 2. Three conditions that remove the return alternative

Under (1), any of the following suffices for strict endpoint descent:

1. The bank word separating the starts has affine multiplier at least one.
   Its carry is positive, so y_j>y_i.
2. The earlier start is the original source. Before the final crossing,
   the cumulative multiplier is at least one, so condition 1 applies.
3. Some later aligned parity bits differ while both lie inside W. More
   precisely, there is a t>=r such that i+t and j+t are both valid bit
   positions of W and their parity bits differ. Equal starting values would
   have identical entire forward trajectories, contradicting that bit.

Thus sufficiently long repeated prefixes, or long matched factors followed
by a witnessed mismatch, are direct endpoint exclusions even when the
matched strings contain resets. An arbitrary matched factor without any
of these conditions has a genuine return alternative; it cannot simply be
declared descending.

## 3. Relation to the primary power theorem

If W contains U^k for a nonempty expansion word U, compare the occurrences
of U^(k-1) starting at the first and second copies of U. Their length is
r=(k-1)|U|, and their separating word U expands. The primary theorem follows
as a special case. This supplies an independent cylinder-difference reading
of the repeated-word fuel argument.

For a general pair of equal bank strings, literal bit lengths and aligned
starts must be checked. Equal slope or equal odd-bit count is insufficient.
A grammar cancellation that changes the actual guarded word supplies no
such pair.

## 4. Sharper remaining hypothesis

A tractable sufficient hypothesis for eliminating a family of actual S_d
witnesses is:

> Every sufficiently long word in that family has two matching parity
> blocks of length at least L/100 at macroblock starts, and either their
> separator has multiplier at least one or a later aligned mismatch occurs
> before the endpoint.

This is strictly broader than finding consecutive pure-expansion powers.
It is still a hypothesis requiring proof for the realized words. The
low-description greedy family in the companion hostile avoids the latter
powers but is not claimed to avoid these more general matches.

## 5. Remove the macroblock-start restriction

The same one-percent conclusion holds for two equal parity factors starting
at arbitrary Terras bit positions. Here is the explicit cost of crossing a
bank boundary, so this generalization does not assume a carry lemma outside
its valid macroblock scope.

Every bank letter has at most 122 bits. At an arbitrary position t before
the endpoint, let P be the last completed macroblock prefix and let y_t be
the current state. For any positive z, T(z)<=(3z+1)/2, so a partial block of
at most 122 bits gives

    y_t < 2^122 (y_P+1)
        <= 2^122 lambda^(L/30)(x+N/6+1).          (6)

The first inequality follows by iterating
T(z)+1<=3(z+1)/2<2(z+1). The second uses the macroblock carry bound; the
factor lambda^(L/30) is at least one. This argument does not require the
partial block itself to be expanding.

Put Y=x+N/6+1. For a non-descending endpoint, x>=r4>1 and the primary
height proof gives

    log Y < B+(B+1)log N.

Two unequal starting values sharing r parity bits differ by a nonzero
multiple of 2^r. Equations (1) and (6) therefore give

    log Y > r log2-(L/30)log lambda-122 log2
          > L/450-122.                           (7)

For L<L0=2^48, N<15L0<2^52 still invokes the exact short-range descent
lemma. For L>=L0, the same derivative argument as in the primary proof
works, with the slightly stronger checked base inequality

    L0/450>53B+174.

Thus (7) is incompatible with non-descent. Equal long factors at arbitrary
bit starts force either endpoint descent or an exact return y_i=y_j. If
their aligned future parity bits differ before either occurrence exits W,
that return is impossible and strict descent follows.

## 6. A precise factor-complexity obstruction

Let p_W(r) be the number of distinct length-r contiguous parity factors of
W. These are literal bit factors at all starts, not bank-letter factors.
Let

    r=ceil(L/100),  P=p_W(r).

If

    P < N-r+1,                                   (8)

then every non-descending lawful source of W enters an exact periodic orbit
by time at most P-1, with a positive period at most P.

Indeed (8) says P+r<=N, so the first P+1 starts 0,1,...,P all support r-bit
factors. Two have the same factor. Section 5 forces their starting orbit
values to agree if the endpoint does not descend. If those starts are i<j,
then i<=P-1 and 1<=j-i<=P. Determinism gives an actual periodic orbit from
time i onward. No assertion is made that the original source x is itself
periodic; it may have a preperiod.

A convenient sufficient numerical hypothesis is

    p_W(r)<=C r+D,  (C+1)r+D<N+1.                 (9)

For fixed C<99 and fixed D, (9) holds for all sufficiently large L, since
N>=L and r<=L/100+1. Thus a uniform linear literal factor-complexity bound
with coefficient below 99 forces an actual internal cycle at every
sufficiently long non-descending first crossing. This is an obstruction,
not a general convergence theorem: a contracting separator can genuinely
have a positive fixed point, and the cycle alternative must be retained.

A sharper right-continuation test is available. Under non-descent, every
pair of equal length-r factors has equal starting values. Consequently
whenever both occurrences have a next bit in W, those next bits must agree;
and all their later aligned bits inside W must agree too. Therefore one
length-r factor with two observed different continuations already proves
strict endpoint descent. A family for which long low-complexity words are
necessarily continuation-branching is excluded outright; low complexity
alone proves only the periodic-orbit alternative.

The fixed regular grammar or O(log m) description in the greedy hostile is
not itself a factor-complexity bound. Section 6 of that companion note
separately proves p_W(r)<=(62/61)r+498 by an explicit rotation partition and
variable-length encoding argument, and obtains this exact cycle obstruction
for every member with m>=5. The proof counts literal factors rather than
substituting program length or automaton size.
