# Logarithmic cutoff for induced Hamilton-count gaps

Follow-up to the [induced-gap and three-slot descent packet](../chanpak_induced_gap_three_slot_descent_20261118/README.md). A background follow-up found the full linear cutoff and the first three gap rows; this note adds a logarithmic cutoff and independently replays the rows. Nothing here assumes Chan–Pak completeness.

Write H(T) for the directed Hamilton-path count, H(empty)=1. Remove every singleton strongly connected component (SCC) of T and delete the same vertices from the marked induced subtournament S, allowing an empty remainder. Ordinal multiplication shows that neither count changes.

**Theorem.** If H(T)-H(S)=2g>0, the resulting normalized host has at most
`min(2g+2, floor(10+10 log_2 g))` vertices. In particular, at any fixed positive gap the *complete* induced endpoint language can be computed by a bounded search. Logarithmic order in g is optimal in general.

## The new step: park the parity choice in a Hamilton cycle

Use the [Grinberg–Stanley packing formula](https://arxiv.org/abs/2307.05569): packings of vertex-disjoint directed odd cycles of length at least three have weight 2 per cycle, and unused vertices are singletons. An induced deletion retains precisely the packings avoiding the deleted set. Thus every omitted directed odd cycle alone uses one unit of half-gap.

Let C be a changing SCC of order n (its original count b exceeds its retained count a), put `d=(b-a)/2`, and take a deleted vertex v. By vertex-pancyclicity, there is a directed triangle D through v. Packings that contain D alone as a prescribed cycle contribute `2 H(C-D)` to the loss, so `H(C-D) <= d`. Suppose C-D has t singleton SCCs. For each other SCC of size s≥3, [Busch's strong lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) gives `log_5 H >= (s-1)/3 >= 2s/9`. Consequently the total order of these other SCCs is at most `(9/2) log_5 d`.

Now take a directed Hamilton cycle of C. Cut its cyclic word at the three vertices of D. Each of the at most three surviving segments lies in C-D and runs forward in the SCC order, with possible repetitions within a nonsingleton SCC. Keep the first and last vertex of every nonempty segment, and keep all nonsingleton-SCC vertices. Every other singleton-SCC vertex can be dropped independently: any new adjacency crossing one or more dropped singletons jumps forward between distinct SCCs, where the needed arc is forced. At most six singletons were made ineligible. All resulting cycles retain D's three vertices, but do not have to use the arcs of D. If there are k≥1 eligible singletons, exactly half of their 2^k deletion choices make the resulting cycle odd. These are distinct cycles through v, so `2^(k-1) <= d`. The inequality `t <= 7+log_2 d` now follows, also when k=0. Therefore
`n <= 10+log_2 d+(9/2)log_5 d <= 10+10 log_2 d`.

The same changing SCC separately satisfies `n <= 2d+2` by vertex-pancyclicity: every odd length from 3 through n uses one unit of d. The proof of vertex-pancyclicity in the prior packet also supplies the Hamilton cycle used above.

## Multiplicative bookkeeping

Let P be the product of the unchanged SCC counts, r≥1 the number of changing SCCs, and `g'=g/P`. For the changing factors write `b_i=a_i+2d_i`, all `a_i,d_i >= 1`. The half-difference of the products expands to
`g' >= 2^(r-1) product(d_i)`.
Thus the changing order is at most `10r+10 sum(log_2 d_i) <= 10+10 log_2 g'`. The Busch bound gives the unchanged nonsingleton order at most `(9/2)log_5 P <= 10 log_2 P`. This proves the logarithmic cutoff.

For the linear cutoff, the same expansion gives `g' >= sum(d_i)+2 binom(r,2) >= sum(d_i)+r-1`. Changing order is at most `2g'+2`. An unchanged nonsingleton SCC of count b has order at most `2(b-1)` by pancyclicity and b≥3, so their total order is at most `2(P-1)`. Since `(P-1)(g'-1)>=0`, the sum is at most `2g+2`.

For sharpness of the logarithmic order, reverse only the extreme arc of a transitive n-tournament, n≥3, and delete the first vertex. The host is strong; its directed odd cycles each use that reversed arc and an odd-size subset of the n-2 middle vertices, and they all intersect. Hence the endpoint pair is `1 -> 1+2^(n-2)`, with `g=2^(n-3)`. Conversely the theorem implies `H(T) <= (floor(10+10 log_2 g))!`, so a sequence of unbounded upper endpoints requires `g >= exp(Omega(log H / log log H))`. This rules out bounded or polylogarithmic exact repairs by induced addition, not by rewiring.

## First three exact rows

| Gap | All endpoint pairs |
|---|---|
| +2 | 1→3, 3→5 |
| +4 | 1→5, 5→9 |
| +6 | 3→9, 5→11, 9→15 |

Two changing factors alone force half-gap at least 4. For g≤3 there is therefore one changing factor, and the only nontrivial unchanged multiplier possible is P=3 at g=3, multiplying the +2 row. The strong host has order ≤8. At orders ≤6, the replay checks every labeled host and every induced subset. At orders 7 and 8, a deleted vertex v must already account for at least three distinct odd cycles by pancyclicity. Hence any gap ≤6 is exact on deleting just v, and its count pair is the same as for the original deletion. Fix a Hamilton path of C-v as the label order, enumerate all remaining arcs, and retain only strong hosts with exactly three odd cycles through v. Compute the actual deletion difference by directed subset DP; the odd-cycle test alone is only a filter. The strong-host row witnesses are printed as n, retained subset bitmask, and outgoing adjacency bitmasks; when the mask is zero, one can retain any singleton instead.

Run:
```sh
g++ -std=c++17 -Wall -Wextra -Werror -O3 replay.cpp -o /tmp/chanpak-gap-replay
/tmp/chanpak-gap-replay | diff -u replay.out -
```
