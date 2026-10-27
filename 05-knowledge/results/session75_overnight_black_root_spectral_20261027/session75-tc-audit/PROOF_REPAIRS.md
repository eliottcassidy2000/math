# Minimal publication corrections for the frozen THM-4812

1. **Greedy sum:** explicitly restrict to `k` for which every factor `c3-i(D+1)` is positive. Stop before the first nonpositive factor; do not sum subsequent re-positive products.

2. **Replace the false pointwise-increasing tail margin:** for `k=floor(m/20)>=20`,
   `c3-k(D+1)>=(11m^3+37m)/480`, hence
   `2^k i_k>=((11/12)m^2)^k>=((1100/3)k^2)^k`.
   Also `23*3^((m-1)/2)<=23*3^(10k+9)`.
   The exact integer inequality `440000^20>23*3^9*177147^20` establishes the k20 base, and `(1100k^2/3^11)^k` increases for integer k>=20. This proves the tail.

3. **Explain the infinite order bound rather than only citing the n<200 script:** once n>=10, every allowed exceptional-orbit selection fits and `B_v(n+3)=3B_v(n)`, while `f(n+3)=5f(n)`. Verify the three residue bases13/14/15 for v7 and23/24/25 for v21, and exclude21 separately.

4. **Explicitly eliminate index5:** among the24 configurations the sole index-cap>=5 configuration is n13 with twoC3 blocks and oneP7 block; its product group has order189, not divisible by5. All configurations with5 dividing the group order have cap<5. Thus the odd index can only be1 or3.

5. **Independent order12 certificate:** add an alternate proof using the same orbit-product/VT reduction and six configurations listed in AUDIT.md. Their caps are1,4,1,3,2,1, so index1/3 suffices. Their45 jobs enumerate46,968 strong generated tournaments, with10,104 exact automorphism calls and45,384 module/Held–Karp checks, no targets. This replaces the unrerun original large order12 computation in the independent acceptance path; do not label the original scan independently reproduced.

6. **Orbit wording:** say "with prescribed G-orbits" rather than implying the scanner enforces equality with full-Aut orbits. Larger automorphism groups cause safe over-enumeration and are checked exactly.

Optional dependency simplification: replace the15-class regular9 input by the short regular-subgroup argument in AUDIT.md and the complete CayleyC9/C3² census. All50 subgroups ofC3 wr C3 were independently checked;13 transitive, all with one of those regular subgroups.

Part(g) added at019882 is not included in this audit's acceptance. Prime-untouchable relative-density-zero is existing S73 work, not a new uniqueness result.
