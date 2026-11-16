# A golden five-adic clock inside a literal rank-four return
The earlier [Gray selector](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/gray_selector_recursion_20261018e.md) distinguishes the full state of the Fibonacci resistance ladder from its projective readout. The [golden rigidity notes](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md) give the ladder and Pisano comparison. The [decimal-window work](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/powers_two_decimal_windows_20261007.md) retains exponent phases when lifting a suffix. The [power-gap notes](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/prime_shells_20260921_elliptic.md) interpret 1024 = 1000 + 24 as an exact identity, not a reusable bounded correction. Here all three mechanisms appear in one actual tournament return: a projective clock, a five-adic lift and an absolute-count obstruction.

## One fixed ribbon word, all powers of five

Use the exact rank-four return lattice of the earlier eight-mask packet. Starting state (s=(-43,-5,4,53)ᵀ), readout (f=(32,14,10,27)), (H_0=25). Repeat the literal eight-letter word (R=(0,6,5,0,3,5,2,4)). It is one of the 168 authenticated returns: it acts invertibly on the same exact rank-four image when iterated. The old integer lattice compiler gives
```
M = [ 5632  2235  1222  4221
       672   267   149   506
      -680  -272  -147  -511
     -6896 -2736 -1498 -5169 ].
H_n = f M^n s.
```

* The characteristic polynomial is (t^4-583t^3+1025t^2-1047t-424). Modulo 5 it factors as ((t-3)^2(t^2-2t-1)). The vector (v=e_4) has (u=(M-3I)v=(1,1,4,3)ᵀ ≠ 0) mod 5 and ((M-3I)u=0). The plane spanned by (u,v) is exactly a nontrivial size-two Jordan block at 3. This is the same conjugacy type as the Fibonacci matrix ([[1,1],[1,0]]) mod 5 because its characteristic polynomial (t^2-t-1) reduces to ((t-3)^2). Thus on this plane the full-state order is 20 and the projective order is 5; for non-eigenvectors both are exact, by (M^n v=3^n v+n3^{n-1}u) mod 5. No √5 eigenvalue over the reals or global similarity of the four-dimensional operator to a Fibonacci ladder is claimed.
* The complementary quadratic is irreducible mod 5 and its root (θ) satisfies (θ^3=2). Its order is 12: its cube has order 4, while (θ^4=2θ
e1). The whole mod-5 operator therefore has order 60. The seed does too, by checking the three maximal proper divisors 12, 20 and 30. The observability matrix with first four rows (fM^i), (0 ≤ i < 4), has determinant 1 mod 5, so scalar and vector periods coincide at every (5^k).
* Put (M^{60}=I+5K). The initial displacement is (Ks=(1,1,4,3)ᵀ) mod 5, verified by multiplication modulo 25. For every a ≥ 0,
  ```
  (I+5K)^{5^a} ≡  I+5^{a+1}K (mod 5^{a+2}).
  ```
  Binomial expansion inductively proves this, so the exact period of (H_n mod 5^k) is (60·5^{k-1}). For the lower bound, a period must already be a multiple of 60 at mod 5; writing the multiplier (5^a d) for (5 ∤ d) shows the corresponding state difference has valuation exactly (a+1).

### Four fixed binary phases, every five-adic suffix

Literally checking the mod-32 orbit gives (H_n mod 32=(25,29,17,1)) in phases (n mod 4), from the start. For each phase in the 60-term mod-5 cycle, every scalar mod-5 residue occurs **three** times. Moreover all 60 initial exponents (j) give nonzero (fM^jKs mod 5) (indeed the values (2,1,3,4) by (j mod 4)). For fixed such (j), the samples (H_{j+60t}) as (t) ranges mod (5^{k-1}) give a **bijection** onto the (5^{k-1}) lifts of that sample's mod-5 value. Induction: the five lifts (t+d5^{a}) change the output at next precision by (d5^{a+1} fM^jKs) with a unit coefficient. Consequently in a full period (60·5^{k-1}), every odd CRT class with one of the four stated mod-32 phases and any specified residue mod (5^k) occurs exactly three times. One may realize any infinite sequence of such requests by waiting at most one period per request, always running the same literal block and keeping rank four. The other twelve dyadic lanes are impossible for this fixed block. This is not steering all odd classes of (32·5^k); that needs a genuine mixed-block controller.

Replay: `python check_golden_decimal.py path/to/return_automaton.json`. Generate that JSON by the earlier eight-mask lattice compiler if needed. The checker verifies the literal matrix and seed against the authenticated bank, and uses integer arithmetic for the finite starting checks.

## Why this excellent clock is not an absolute-digit recursion

The historical identity (2^{10}=10^3+3·2^3) says ((2^{10})/(10^3)=128/125). Iterating it alone does not supply arbitrarily small relative corrections: the relative excess after (m) steps is ((128/125)^m-1), and (1024^m-1000^m) has exactly (2)-valuation (3m) and (5)-valuation zero. A state with hidden scale cannot absorb this discrepancy merely by changing its phase. The old suffix work's order (4·5^{k-1}) of 2 gives order (2· 5^{max(k-2,0)}) for 1024 mod (5^k), while 1000 itself is a nonunit. Thus even approximately aligned scales behave differently in the actual five-adic phase.

There is an unconditional **absolute sparsity obstruction** for the previous rank-four bank, not just for the single golden orbit. The verified [isolated triple minimum](https://github.com/eliottcassidy2000/math/blob/3a2423f0bb11b8a63dc81b7b3616a881c101b3f6/05-knowledge/results/chanpak_eight_braid_20261116/RETURN_GROWTH.md) is (B=13,460,883>168^3). For a fixed seed, an (m)-block word gives (H ≥ 25 B^{⌊ m/3⌋}) by disjoint-support odd-cycle packing. Yet there are only (168^m) block words. Summing over the permitted lengths gives
```
#{H ≤ X}=O(X^{3 log(168)/log(B)}),
3 log(168)/log(B) ≈ 0.93644.
```
Distinct values are counted at most by the number of words, not assumed to coincide. A uniformly steerable mixed decimal controller using only this bank would still omit density one of absolute targets. Every single bank block has (H′ ≥ 55H), so such a block cannot implement the usual one-decimal-digit update (H'=10H+d), (0 ≤ d ≤ 9). Two blocks have (H′ ≥ 27143H), so cannot implement a three- or four-decimal-digit update either.

The bound survives a small fraction of escape moves. On the fixed eight-letter grid allow at most (e ≤ m/100) blocks outside the bank, each unconstrained among (D=8^8) words. Consecutive good runs have at most two leftover blocks each after forming triples. There are at most (e+1) such runs, so (H ≥ 25 B^{(m-3e-2)/3}), regardless of whether the escape moves retain rank. With (A=168, z=A/(100D)), the number of permitted words is at most
```
z^{-m/100}(A+D z)^m .
```
Indeed for each term of ((A+Dz)^m) with (e ≤ m/100), its factor (z^{e-m/100}) is at least 1. The exact integer comparison ((A+Dz)^{300}z^{-3}<B^{97}) shows that this enlarged language still has a power-saving sparse absolute spectrum (exponent below 0.998). This is only a fixed-grid, fixed-seed obstruction; frequent excursions, other literal words, other seeds and Chan–Pak itself remain open. It gives a quantitative meaning to “least resistance” here: making the next residue cheap in a finite quotient is insufficient if the absolute height grows faster than the number of available words.
