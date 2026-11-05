# Independent audit: unit determinants and physical sign switching

**Status: PASS.** This audit checks the all-order arguments and supplies independent finite controls. It makes no global Hamilton-path spectrum or literature-priority claim.

## Frozen target

The reviewed final source is copied under `frozen/`. Its proof SHA-256 is `1e3a7915f014b129918aff5a47b24e1d02abf5b72b4025356143adeaa0939fcf`; the verifier SHA-256 is `eadbd6b5ac018845c3800a4b7ea1aa6e521ca55b2547805ae4631d0de6a3e799`. The companion synthesis proof matches the previously audited SHA-256 `31982897f0092fafa040c9d60452fb2b69c128d4b7173b5d642d3ff56a622499`.

Independent normal and optimized replays of the author's frozen verifier pass and are byte-identical to both supplied receipts. Those four JSON files share SHA-256 `b597a95ed38f7667c2a80a8f4d3c4739472afcf9a2d30e06a7f7270652bd71b0`. The final independent C++ source was compiled warning-free and replayed, with byte-identical receipts.

## Reviewed statements

1. **Deletion-one recursion.** A deletion-one exterior is transitive. In its source-first order, the bits `[c -> a_i]` form a binary integer q. Counting the first vertex after c gives `R=q+1`, while the last exterior vertex gives `e=1-(q mod 2)`, including the empty-word convention. Appending zero gives `U'=2U`, `e'=1`, `sigma'=sigma`; appending one gives `U'=2U+1`, `e'=0`, `sigma'=sigma+U+1-e`. These yield precisely the stated recursions. Leading zero padding preserves `(h,e,R,d,U,sigma)`. It need not preserve the unused left response L of X.
2. **Fixed slack.** Removing terminal zero bits reduces positive U to its odd core v. The odd recursion implies `v<=2 sigma+1`; conversely each admissible odd core and every terminal-zero extension are actual apex words. The zero-slack and slack-two classes are exactly the stated powers-of-two rays. The additional U=0 is not a positive amplifier. The finite-partner corollary follows by solving `sigma=(g-beta)/t` for each partner and sign, then taking `k=j+m>=1`; it does not silently use deletion m=0.
3. **First missing slack.** The eleven odd cores up to 21 have slacks `0,2,2,6,4,8,8,14,8,14,12`, so none has slack ten. Combined with the all-order odd-core bound this is a genuine exclusion at every order, rather than an unbounded claim inferred from a raw search. The actual partner with `(g,t,L)=(11,1,13)` is verified directly. The alpha-eleven positive contrast has deletion count five and does not contradict that exclusion.
4. **U=1.** Nonnegative even slack forces `sigma=0,h-e=1`. Zero-slack geometry gives a unique out-neighbor a. The singleton-birth identity gives `R_a(X-c)=h-e=1`, forcing deletion count one. The physical odd count then forces exactly `(h,e,R,d)=(1,0,2,1)`. Arbitrarily many neutral leading vertices do not evade this argument.
5. **L=1 and odd prime amplitudes.** `L>=g>=1` forces transitive Y, and its exact left response `L=2^j` forces c to be the source, including singleton Y. The alternative factorization U=1 would force beta=+1 and transitive Y with power-of-two L, so cannot produce odd prime amplitude. Thus all stated necessary conditions hold. At deletion count one, the only odd-prime possibilities are amplitudes 3 and 5 with beta=-1. The theorem restricts amplitudes, not all tournament Hamilton counts.
6. **Physical sign switch.** Doubling the marked vertex makes a literal tournament Y with `g=2H-beta` and `t=H`, for either choice of its two clones as c. The fixed three-vertex X has `U=3,d=1,sigma=2`, hence the synthesized new clone ray has offset `-beta`. The graphs add four vertices and retain an induced original input. The cycle through c and all three new vertices proves strongness when Y is strong. Repetition is legitimate because the newly marked vertex is again an actual positive binary amplifier.
7. **Same-count distinction.** The switched operation lives at `N=6L-beta>=12H-7beta`, not at the old count H. The source explicitly preserves this distinction. It neither supplies two fixed commuting clone buttons nor proves all-target reachability. Its amplitude is composite and divisible by three.

One scope correction was required in the newly added finite-partner corollary: it must restrict to U>0, or separately include the constant count one. The literal singleton pair X=Y gives U=0, beta=+1 and H(F_m)=1, so it is not a positive dyadic ray. The corrected source restricts the corollary to U>0. No change to the main theorems was needed. The proofs cover arbitrary order and positive clone multiplicity; the finite tests below do not supply the all-order quantifiers.

## Independent computation

`audit.cpp` was written independently and does not import or call the author's verifier. It constructs every tested graph at adjacency level, counts paths by subset endpoint DP, checks all profiles through four vertices by literal permutations and cuts, and uses a recursive bipartite matching search to identify binary-amplifier inputs. Its clone builder uses an ordered owner list rather than the author's exterior-first indexing. The fixed sign compiler retains Y's vertices in place and adds its three new vertices directly.

The exact finite universes are:

- All 33,867 labeled tournaments of orders one through six and every root: 202,013 rooted graphs. Physical inequalities, slack parity, all 872 U=1 cases, and all 873 L=1 cases pass.
- All 273 distinct relevant X profiles from that universe paired with all 46 distinct relevant Y profiles through order five. All seven odd-prime/unit-offset profile pairs satisfy the theorem. Both amplitude-two/unit-offset pairs have positive offset and deletion count one.
- All 2,047 incidence words of lengths zero through ten checked as actual apex graphs, including arbitrary leading zeros. Another 255 checks add two leading zeros explicitly.
- All integers `0<=q<2^20` checked against the slack recursions, odd-core bound, and exact slack-zero/slack-two classifications.
- Every labeled unit-offset amplifier through order four, plus one representative of every new `(H,beta,L_first,L_second)` signature occurring at order five. Both choices of clone c are tested, producing 332 physical compiler cases and 996 actual counts at new clone multiplicities one, two, and three. The synthesized graphs agree exactly with direct cloning of the newly marked vertex. The 678 applicable strongness checks pass.
- All eleven odd cores in the exact slack-ten exclusion, the displayed blocked partner profile, and the alpha-eleven positive-contrast graphs at multiplicities one through five. The neutral U=0 singleton pair is retained as a boundary control.
- Two successive sign switches from the directed triangle for each of the four ordered choices of first/second marked clone. All actual second-step counts, signs, and strongness checks pass.

The larger search finds 2,790 unit-offset rooted inputs through order five, represented by 22 compiler signatures. Only the explicitly described subset is used for the more expensive physical compiler check; no complete-order-five compiler census is claimed.

## Reproduction

From this directory:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic audit.cpp -o audit
./audit > receipt.json
python3 frozen/verify.py > replay.normal.json
python3 -O frozen/verify.py > replay.optimized.json
cmp replay.normal.json replay.optimized.json
cmp frozen/receipt.json replay.normal.json
```

All audit writes are confined to this dedicated scratch packet. The parent owns any publication.
