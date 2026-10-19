# Every greedy first-crossing word has strict endpoint descent

19 October 2026. **PROVED; fresh independent audit requested.**
This theorem closes the entire greedy V/F family, not only its
upper-convergent subsequence. It uses the audited matched-block height
theorem, including its cited logarithmic-form input. It does not establish
arbitrary bank entry, general Collatz convergence, or the full S_d boundary.

## 1. Precise statement and language

Use T(x)=x/2 for even x and (3x+1)/2 for odd x, with chronological bank

    K=110110110110110101011011101010,
    V=K^4 F, F=10.

V has 122 bits, slope v=3^77/2^122>1, and every nonempty literal prefix of
V has multiplier greater than one. Put mu=3/4 and

    theta=log(v)/log(4/3),  0<theta<1.

Starting with cumulative multiplier S=1, repeat m times: append V; then
append F if and only if its output cumulative multiplier would still be
at least one. Finally append one additional F. Call the resulting word
W_m and its Terras length N. Its expansion-bit length is L=122m.

The exact invariant 1<=S<4/3 after each completed iteration proves that
W_m is expansion-first, all iteration checkpoints precede the first crossing,
and the final extra F is its first explicit-reset coefficient crossing.
The previous notes prove the invariant by rational comparisons. In particular
all these claims concern the fixed specified parse, not a reordered word.

**Theorem.** For every m>=1 and every positive integer source following the
entire lawful word W_m,

    T^N(x)<x.

Thus no W_m witnesses x in S_d for any d>=0. Every W_m nevertheless has an
infinite positive lawful cylinder; the claim is a strict endpoint exclusion
on that cylinder, not the absence of admissible sources.

## 2. The two arithmetic inputs

### Equal long futures force descent or an exact return

The audited arbitrary-start matched-block theorem says: in an
expansion-first first-crossing bank word, two equal r-bit future factors
with 100r>=L force either strict endpoint descent or equality of the actual
states at their starts. They can overlap. They must both be wholly present
in the actual finite word. This is a theorem about lawful integer orbit
values; it is not a formal pumping rule.

### The two encoded symbols have disjoint output images

Encode a reset indicator by

    0 -> V,  1 -> VF.

Every lawful V output is 82 modulo 243, and every lawful VF output is 1
modulo 243. Indeed V has b=103 and M=22 modulo 243, with 22^(-1)=232.
Thus its output is 103*232=82 modulo 243. Applying F changes that residue
to (3*82+1)/4=1 modulo 243. Equivalently the exact values of v_3(y-1) are
4 after V and 5 after VF.

This uses an output-image guard, not an assumption about the source's
3-adic distribution. It holds for every actual integer input on the stated
lawful branch. The fixed-point exception F(1)=1 causes no problem: a V
output has finite type 4 and therefore cannot be 1.

## 3. A finite minimal-left-occurrence lemma

Let s_1...s_m be a finite word over an alphabet, and suppose some length-k
factor occurs twice. Among all equal-factor pairs choose one whose earlier
zero-based start i is as small as possible; write its later start as j>i.

Then either i=0, or the immediately preceding symbols differ.

For if i>0 and s_i=s_j in one-based notation, shifting both starts one
place left gives two equal length-k factors: their first symbols are equal
by assumption, and their remaining k-1 symbols agree by the original
match. The earlier start becomes i-1, contradicting minimality. The shifted
factors remain inside the same finite word. The argument also applies to
overlapping factors and to k=1.

The point of minimizing the earlier start is essential. An arbitrary
repeated pair may have equal preceding symbols and retain a genuine cycle
alternative.

### Stronger equivalent structural proof: reverse-injective coding

In fact all completed-iteration checkpoint values y_0,...,y_m are pairwise
distinct, without any factor-complexity or endpoint assumption beyond
lawfulness and the pre-crossing invariant.

Suppose y_i=y_j with 0<=i<j<=m. If i>0, their equal output values determine
the same preceding code symbol, since the V and VF output-image classes
are disjoint. That same branch is an injective affine map, so equality of
its outputs implies y_(i-1)=y_(j-1). Repeat this argument i times. The result
is y_0=y_(j-i), where j-i>=1. But the nonempty completed-iteration prefix
has multiplier at least one and positive carry, so y_(j-i)>y_0. This is a
contradiction.

This cancels identical actual past branches after their types have been
identified. It neither deletes an unmatched reset nor assumes that the
full Collatz map is globally invertible. The initial source anchor and the
separated output images are essential.

Consequently any repeated reset-symbol factor can be used in Section 5,
without minimizing the earlier occurrence: its starting values are already
known to be different. The minimal-left-occurrence proof remains a direct
finite-word version of the same obstruction.

## 4. Exact factor complexity of the reset-symbol prefix

Let f_n be the indicator of the conditional F after the nth V, before the
final extra F. The multiplier invariant gives

    f_n=floor(n theta)-floor((n-1)theta).

Theta is irrational by unique prime factorization: theta=p/q would imply
3^(77q+p)=2^(122q+2p). A length-k factor is a function of the phase
z={(n-1)theta}. Its coordinates change only at the k points {-j theta},
1<=j<=k. The floor convention makes this a partition into k+1 half-open
intervals, including their assigned endpoints. Hence the number p_f(k) of
distinct length-k factors in the finite prefix f_1...f_m satisfies

    p_f(k)<=k+1.                                 (1)

This proof only bounds factors that actually occur in the finite prefix.
It does not claim that every infinite rotation factor appears there, or
use an unproved finite recurrence estimate.

## 5. Proof for every m

For m>=3 set k=ceil(m/100). The exact inequality

    m>2k

holds: k<=(m+99)/100, and 2(m+99)<100m for every m>=3.
There are m-k+1 length-k occurrence positions. By (1),

    m-k+1>k+1>=p_f(k),

so a repeated k-factor exists. Choose a pair i<j by Section 3.

Let t_i be the literal Terras length of the first i complete code images;
t_0=0. Equal k-symbol factors encode identical literal parity strings of
some length r. Variable image lengths 122 and 124 cause no ambiguity:
the same finite symbol sequence is encoded on both sides. Moreover

    r>=122k, 100r>=122m=L.

Since j+k<=m, both encoded futures lie inside the encoded prefix
f_1...f_m, before the final extra F. Both starts are genuine macroblock
checkpoints strictly before the final first crossing.

There are two cases.

- If i=0, the earlier state is the original source x. The later completed
  iteration prefix has cumulative multiplier at least one and strictly
  positive carry, so its actual value is greater than x. The two states
  cannot be an exact return.
- If i>0, Section 3 makes the preceding symbols different. One checkpoint
  is therefore a V output and the other a VF output. Section 2 gives
  their different residues 82 and 1 modulo 243, so they cannot be equal.

The matched-block theorem now forces strict endpoint descent in both cases.
The initial-source case is why the coded prefix is anchored at the original
source; silently inserting an earlier arbitrary context would invalidate
that part of the proof.

For m=1,2, the words have at most 250 Terras bits. The audited exact
short-range theorem for expansion-first first crossings, N<=2^52, proves
strict descent. This completes the all-m proof.

## 6. Reusable terminal-image criterion

The same finite argument works for a coded word over any finite alphabet
under the following explicit hypotheses:

1. The coded word begins at the original positive source, and every factor
   start used is before the first coefficient crossing. A later nonempty
   coded prefix has multiplier at least one and positive carry.
2. Distinct code symbols have disjoint lawful output-image classes,
   independent of the integer input entering that code symbol.
3. Every code image has literal length at least ell, some k satisfies
   100 ell k>=L, and the number of distinct length-k symbol factors is
   smaller than the number m-k+1 of occurrence positions.
4. The ambient lawful word is within the proved matched-block theorem's
   expansion-first bank scope; any final tail and every matched future
   are explicitly retained.

Each fixed code branch is injective, as every Collatz affine branch has
nonzero multiplier. Reverse decoding therefore also proves that all coded
checkpoints are pairwise distinct. Its endpoint strictly descends. A minimal
earlier occurrence alternatively either
begins at the source or has distinct preceding output classes, exactly as
above. Injectivity of terminal images is the required arithmetic property;
unique decipherability or prefix-freeness of the literal code is not needed.
The entire original bank alphabet does not automatically satisfy condition
2: F's image type depends on its input and cannot be given a false constant
type just to apply this lemma.

## 7. Significance and remaining boundary

The greedy words have floor(m theta)+1 explicit resets, tending to infinity.
Their expansion runs have length at most ten V's. Thus all sufficiently
large W_m avoid the previous one-percent pure-expansion power criterion,
even after K^4F<->V rewrites. The present theorem closes their actual endpoint
realizability, rather than raising a finite episode cap.

This closes the endpoint-descent goal behind the typed-left-special question
recorded in the frozen image-guard packet. The guaranteed useful pair is
either source-anchored or has different incoming types; finite pigeonhole
plus minimality proves exactly this disjunction. The original occurrence
hypothesis requiring different incoming types in every sufficiently large
prefix is not established here. No irrational-rotation recurrence bound is
required.

Other first-crossing bank families may have high factor complexity,
nonseparated terminal images, or no anchored coded representation of this
kind. Their realized S_d sources, the complete 2-adic limit boundary, and
arbitrary bank entry remain open.

## 8. Fixed-template generalization, with the one-copy boundary explicit

Fix any nonempty expansion word U over {K,V}; its length ell=|U| may be
arbitrarily large. After each copy of U, append as many F resets as possible
while keeping the cumulative multiplier at least one. After m such
iterations append one more F, giving the first crossing. Call this word
W_(U,m).

**Corollary.** Every positive lawful source of W_(U,m) has strict endpoint
descent for every fixed U and every m>=2.

Write s_U=3^h/2^ell>1 and theta_U=log(s_U)/log(4/3). The completed iteration
multiplier lies in [1,4/3). Its number j_n of conditional resets is exactly

    j_n=floor(n theta_U)-floor((n-1)theta_U).

As before, theta_U is irrational by unique prime factorization. If
a=floor(theta_U) and beta=theta_U-a, then j_n is either a or a+1 and
j_n-a=floor(n beta)-floor((n-1)beta). The binary reset-choice sequence has
at most k+1 length-k factors by the same rotation-partition proof.

The two code images are U F^a and U F^(a+1). The last expansion letter of U
is K or V, so U's output type is a fixed t in {3,4}. The two coded outputs
have distinct finite types t+a and t+a+1, independent of their lawful
integer inputs. Each branch is injective. All complete-iteration
checkpoints therefore are pairwise distinct by reverse decoding.

For m>=3 use k=ceil(m/100), exactly as above. The expansion length is
L=m ell and each matching k-symbol future has at least k ell bits, so it
meets the same height threshold. This proves descent.

For m=2 no factor-complexity estimate is needed. Compare the two literal
U futures beginning at the original source y_0 and at y_1, the start of the
second iteration immediately after the first complete iteration. They have
the same length ell=L/2, and y_1 is strictly greater than the source because
its multiplier is at least one and its carry is positive. The matched-block
theorem proves descent.

**No universal m=1 assertion for arbitrary U is made.** Such a word can be
arbitrarily long, so the short-range theorem cannot automatically supply
that missing case. The earlier all-m theorem for U=V handles m=1 separately
because its actual length is small. Nor does this construction grant an
entry theorem for an arbitrary integer or an arbitrary interleaved word.
