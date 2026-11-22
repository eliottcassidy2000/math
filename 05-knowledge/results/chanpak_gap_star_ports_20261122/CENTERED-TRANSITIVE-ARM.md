# Transitive arms at a prescribed perfect-number center

**PROVED.** Let `M=2^p−1` be prime, with `p≥3`, and put `E=2^(p−1)M`. Consider the two-arm star of width `m≥1`. Prescribe its **growing term** to equal `E` exactly. No such star with transitive right arm has count `E−1` or `E+1`, for any left arm or either tag set. If the left arm is transitive, any surviving count `E±1` has an empty-left star representation on the same right arm at a possibly greater width. This includes exponent 61. General stars can reach the same target with a different growing term and a compensating intercept.

## Reading a free Boolean cut

For a tagged arm let `h` count Hamiltonian paths and `k` count them after insertion of an apex beating precisely the tags. Let `P` count nonempty paths on subsets starting tagged, each weighted by the Hamiltonian count on its complement. Let `e` restrict this sum to even orders, and `J` further restrict it to paths ending tagged. The signed endpoint identity gives `k=h+2(e−J)`. In particular `k≤h+2P`. A tournament has count one exactly when transitive, by the odd-cycle packing identity since every nontransitive tournament contains a directed triangle. Writing `C=h_A+P_A` and `g=h_B`, the star output is

\[
                  2^m CP_B + I,\qquad I=gk_A-2J_B h_A.
\]

An empty left arm has `h_A=k_A=C=1`. [Star construction](../chanpak_moonshine_prime_blocks_20261118/STAR_SEAM_PRODUCT.md) · [Signed endpoint identity](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md)

On a transitive arm, list the tags as bits in increasing tournament order, with an arc directed from each earlier vertex to each later one. The integer with those bits from most to least significant is precisely `P`. In fact, a path whose first vertex has `d` successors can independently take any of them, contributing `2^d` when that vertex is tagged. If there are `t≥1` tags, their nonempty even subsets already give `J≥2^(t−1)−1`. A single tag gives `J=0`.

## A transitive right arm is excluded

Since `2^m CP_B=E`, exactly one of `C,P_B` contains the prime `M` and the other is a power of two.

* If `P_B=M 2^j`, the tag word has `p` ones, `C≤2^(p−2)`, and `g=1`. As `k_A≤2C−h_A`,
  \[
       I\le 2C-(2J_B+1)h_A
         \le 2^{p-1}-(2^p-1)h_A
         \le 1-2^{p-1}\le-3.
  \]
* If `P_B` is dyadic, it has a single tag. Thus `I=k_A≥1`. At equality the left arm with its inserted apex is transitive. Deleting the apex gives a transitive arm whose tags form a suffix, so `P_A=2^s−1` and `C=2^s`, including the empty suffix `s=0`. This contradicts `M|C`.

The hypothesis `p≥3` matters. At `p=2`, take an empty left arm, a transitive right arm with two tags on two vertices, and width one. Then `P_B=3, J_B=1`, and the resulting count is `6−1=5`.

## A transitive left arm reduces to the empty arm

Suppose `A` is transitive. If `C=M 2^a`, its tag word is the binary expansion of `C−1`, namely `p−1` ones, one zero, and `a` ones, up to leading zeros. Those leading zeros and trailing ones are ordinal sources and sinks in the apex tournament and can be stripped. The remaining tournament is a transitive order with a single reversed arc from the last vertex to the first, enclosing `p−1` interior vertices. Its unique increasing Hamiltonian path persists. Each path using the reversed arc is determined by the subset of the interior vertices placed before that arc; the others follow it in increasing order. Thus `k_A=1+2^(p−1)`. Since `P_B≤2^(p−2)` and `J_B≤P_B`,

\[
                  I\ge (1+2^{p-1})g-2^{p-1}\ge1 .
\]

Equality forces `g=1` and hence a transitive right arm, already excluded.

The remaining possibility is `C=2^a`. The tag word for `P_A=C−1` is a suffix of ones, whence `h_A=k_A=1`. The entire ray becomes

\[
                    H=2^{m+a}P_B+g-2J_B .
\]

This is the empty-left formula at width `m+a`. Any new centered two-arm representation beyond this one-arm reduction therefore requires two nontransitive arms. In the reduction the necessary and sufficient right data are `g−2J_B=±1` and `P_B=M2^j` for some `0≤j≤p−2`, with empty-left width `p−1−j`. Existence of such an arm at exponent 61 remains open here.
