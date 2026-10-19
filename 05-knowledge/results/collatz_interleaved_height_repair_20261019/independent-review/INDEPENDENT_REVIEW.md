# Independent review: 24 interleaved bank episodes

19 October 2026. **Verdict: PASS for the stated restricted first-crossing theorem.** No substantive proof gap or numerical discrepancy was found. This is a computer-assisted descent-or-smaller-cycle theorem under explicit bank-language and episode-count hypotheses, not a proof of general Collatz convergence.

Audited document: `INTERLEAVED_HEIGHT_REPAIR.md`, SHA-256
`744401a10e3ce93bf3a62c5f5efdd6212f349670cc41675ee2ccf74a64a74d5c`.
The predecessor `ORIGINAL_HEIGHT_REPAIR.md` and its independent review were also consulted. The new checker imports neither author's verifier nor the predecessor's checker.

## Proof review

1. **First crossing is defined at the right granularity.** The definition examines each explicit F/H endpoint, including intermediate endpoints of a consecutive reset run. Thus every earlier explicit reset endpoint has cumulative multiplier at least one. Every K/V block and each of its nonempty prefixes has multiplier greater than one, so all macroblock endpoints before N also have cumulative multiplier at least one. V's internal F is consistently part of V, and all 122 V prefixes were checked. Using only the final endpoint of each reset run would not justify the same argument.

2. **The carry estimate is valid for arbitrary episode lengths.** Writing `S_i=s_1...s_i`, the total normalized carry is `e=sum(e_i S/S_i)`. The bank envelopes give `e_i<C s_i+1`, and the first-crossing hypothesis gives `S_i>=1` for `i<m`. Consequently

       e < CmS+(m-1)S+1 = 1+(m(C+1)-1)S < m(C+1).

   The last inequality holds for every `m>=1`, since `C>0` and `0<S<1`. No individual episode must contract. For `m<=24`, a nonnegative integer endpoint displacement is therefore in `0,...,12798`. The exact bound is independent of the distance of S from one.

3. **The forced source cylinder is sound.** If the initial expansion run is K, K², or K³, the first F already contracts and descends for every lawful source `x>=815464443`; H is pointwise smaller. Otherwise the literal word begins with K⁴. Independent dyadic lifting recovered

       r4=588943162011782273059357097228697595 modulo 2^120.

   Since `0<r4<2^120`, every positive source in this cylinder has `x>=r4>466400`.

4. **All 12,798 positive offsets are certified on all lifts.** The independent checker computes affine carries from the closed sum over odd positions, rather than the author's forward recurrence, and reconstructs short source guards by dyadic lifting and literal replay. Every first-contracting prefix has length at most 92; its modulus divides `2^120`; and its exact threshold `(a d+b)/(2^j-a)` is below 466400. The residue congruence propagates the complete parity word to every `x+d` with `x=r4+u2^120`. The uniform affine inequality then proves actual descent below x. Sample lifted replays are controls, not the proof of the universal quantifier.

   All offset rows reproduce the published digest. Maximum depth is 92, uniquely at `d=4`. The largest threshold is `810692053755261919/1738366812781`, uniquely at `d=11808`. An independent prefix-free-cylinder argument and exact arithmetic-progression counts verify the 718 disjoint classes cover the entire offset interval.

5. **Zero displacement proves the stated distinct-cycle alternative.** Before N, complete macroblock endpoints have multiplier at least one and strictly positive carry. Interior K/V prefixes increase their input. Every proper F/H prefix is strictly above that block's endpoint at inputs greater than one; this includes prefixes of the final reset whose endpoint is x. Thus `T^t(x)>x` for every `1<=t<N`, establishing minimum x and least period N. The independently lawful companion yields `T^N(x+1)=x+1`, and its initial 00 gives the periodic integer `(x+1)/4`, strictly between 2 and x. Its cycle is nontrivial and distinct from the lower cycle. This does not turn companion descent into lower-orbit descent. These minimum and periodicity facts also establish exclusivity of the two alternatives and the least-counterexample corollary.

6. **Exact normalization preserves the relevant crossing.** A replaced K⁴F segment before N starts with cumulative multiplier at least one. Its full V multiplier, and every V prefix multiplier, exceed one. The removed explicit F therefore cannot be the first crossing. Parity words, guards, affine coefficients, and N remain unchanged; the episode count cannot increase. Finite controls over all parses of at most eight macroblocks that have a crossing include 106 actual rewrites and 36 episode mergers. The preceding symbolic argument, not that enumeration, covers arbitrary lengths.

7. **The two limitations are accurately separated.** At `d=12852`, all 120 inherited prefix multipliers are at least one. With positive offset and nonnegative carry, none can descend below x. This disproves that particular universal inherited-prefix lemma for the 25-episode carry budget; it does not prove realizability of that offset, failure of the 25-episode theorem, or sharpness of 24. The irrational-rotation construction correctly gives first-crossing words with unbounded normalized carry `b/2^N`. For sufficiently large L, `j_L>=2`, so the allowed K⁴F rewrites still leave resets between episodes; exact normalization does not remove this obstruction. Large carry alone makes no claim of actual non-descent on lawful lifts.

## Exact controls and reproduction

`check_interleaved.py` uses only Python's standard library, exact integers/rationals, and explicit exceptions instead of assertions. It reports **267,283 checks**. In addition to the offset certificates, it verifies every one of the 72 earlier reset endpoints of `(K^11 F^3)^24 F`, reconstructs its carry by an independent closed sum, and replays a positive lawful source through all 8,066 bits. Further exact controls use `L=69, j=19` for 24, 100, and 1,000 episodes; they satisfy the sufficient inequalities for the unbounded-carry family.

From this directory:

```sh
python3 check_interleaved.py > /tmp/interleaved-review-normal.json
python3 -O check_interleaved.py > /tmp/interleaved-review-optimized.json
cmp /tmp/interleaved-review-normal.json /tmp/interleaved-review-optimized.json
cmp /tmp/interleaved-review-normal.json receipt.json
sha256sum -c SHA256SUMS
```

Both independent executions reproduce `receipt.json`. The author's normal and optimized executions were also run and matched its published receipt byte for byte: PASS, 102,714 explicit checks. No general-source census, external literature/priority review, or proof that arbitrary orbits enter a bounded-episode bank prefix is supplied by this audit. That pointwise coverage question remains open.
