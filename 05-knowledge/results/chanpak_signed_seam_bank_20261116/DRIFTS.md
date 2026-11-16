# A modular clock with finite signed moves

An exact mirrored +2 repair requires the boundary pair z=(L,T) to obey `N(z')=N(z)−2 mod p` at every prime p dividing the core count, where `N(x,y)=x²+y²`. This note studies the resulting finite selection problem, not whether that selection lifts to actual masks or fixes the absolute count. It is a necessary test to apply *before* searching for a bank of physical moves of a specified drift pattern.

## Periodicity and minimum bank size at an inert prime

Let p≡3 mod 4 and d=z'−z. The condition for the signed direction d is the affine line

    z·d = −1 − N(d)/2.                         (1)

Every direction of a valid step is nonzero. If a word of k signed directions is run twice, with total displacement D over the first pass, subtract (1) for each letter between the two passes. This gives d·D=0 for each letter, hence N(D)=0 on summing. Since p≡3 mod 4, D=0. Telescoping the norm drops in the first pass gives p|k. Thus **no periodically repeated word shorter than p** can underlie a physical mirrored +2 script at such a prime. One successful pass is not enough. The bound is sharp in the unrestricted modular graph: every field value is a sum of two squares (the two sets of size (p+1)/2 must intersect after translation), so choose one z per successive norm value and close after p. This sharpness does not imply a physical lift.

In particular two letters cannot provide four successive steps. Equal consecutive letters are impossible by (1) because N(d)≠0. Four steps on two letters would be ABAB, a double pass of period 2, which is forbidden. For split primes p≡1 mod 4 there is instead a one-letter modular clock: choose nonzero d with N(d)=0, then a starting point z with z·d=−1. This is the arithmetic contrast to the restrictive transitive carry law and to the particular nonisotropic drift used in the 13-path example.

### Stronger bound for a full cycle with few directions

In a closed norm-step cycle of length p, all starting points of steps are distinct, because their norms are distinct. Partition the direction alphabet of size m into parallel classes of sizes m_i. Two successive nonparallel letters (d,e) fix a unique starting point: intersect z·d at its value from (1) with (z+d)·e at the value for e. Each ordered pair of nonparallel letters can therefore occur at most once in the cycle.

A maximal run using a single parallel class i lies on a single spatial affine line, and each particular letter d in that class is usable at most once in this run: along the spatial line the quantity z·d is injective, since N(d)≠0. Hence the run length is at most m_i. A cycle cannot use only one parallel class, because N restricted to an affine spatial line is a quadratic with at most (p+1)/2 distinct values. Charge every run to its outgoing ordered pair (last letter, first letter of the next run). There are at most m_i(m−m_i) such pairs for class i. Consequently every feasible p-cycle obeys

    p ≤ Σ_i m_i² (m−m_i).

In particular **at least four directions** are needed for p≥7: for three letters the possible nontrivial partitions 1+1+1 and 2+1 both have bound 6. If all directions are pairwise nonparallel the sharper condition is p≤m(m−1). This restriction concerns modular cycles or twice-repeated words, not a finite noncyclic path with changing drift.

## Why this is the right archive transfer

- [THM-2682](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-2682-central-arrival-clock-trap-and-three-event-dilation-nilpotence.md) has a physical two-event product that stops at the third event because the fixed arrival clock forces equal shallow labels. Here the exact analogue is the repeated-letter subtraction in (1); checking individual steps or a single pass is not enough.
- [THM-3145](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3145-bass-serre-two-three-tree-and-tetrahedral-congruence-quotient.md) separates local quotient edges, cycles in a finite quotient, and closed lifts in a Bass–Serre tree. The direction clock above is deliberately only a quotient: a modular cycle cannot by itself be a closed +2 trajectory of physical tournament states, since a genuine closed return has net change zero.
