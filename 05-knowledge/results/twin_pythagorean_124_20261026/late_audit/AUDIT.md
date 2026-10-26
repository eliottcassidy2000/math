# Bounded independent audit of the late S60 corrections

**PASS for the four specified repairs and the added power-of-two-center
synthesis. This is not a full audit of THM-4805.** No source or earlier
audited artifact was edited by this reviewer. Unrelated assertions,
computational censuses, C'(1)–(4), and broader open hypotheses retain their
existing status, including the pending full-theorem audit.

## Reviewed claims

1. **Small-source domain and divisor count.** The classifier now concerns
   `2<=v<=14`; `s(1)=0` and the prime fibre at target 1 are separate.
   It counts proper divisors, `t=tau-1<=3`. The repaired source shapes
   `p²,p³,pq` are exactly the composite shapes with that bound. A source
   with four proper divisors is `p⁴` and has sum at least 15; with five it
   is `p⁵` or `p²q`, of sum at least 31 or 16. Six or more distinct proper
   divisors have sum at least 21, already outside this target interval.
   Thus the exclusion needed for the repair is valid. The witnesses
   `s(16)=15`, `s(12)=16`, and `tau(6)=4` verify the stated boundaries.

2. **Mersenne singleton scope.** S57's source statement lists exactly the
   odd singleton targets `3,7,9,11,39` below `10^8`. Among `2^k-1` for
   `1<=k<=26`, only 3 and 7 occur in that list; target 1 has infinitely
   many prime sources. This is correctly labelled CITED FINITE-EXACT,
   with no claim that the present repair recomputed that census.
   For every `k>=4`, the distinct-prime Goldbach antecedent applied to
   `2^k-2` gives an odd semiprime source `pq` different from `2^k`.
   This proves precisely the stated conditional global iff, not an
   unconditional global classification. Exponent zero is rightly excluded.

3. **The missing unit in the ray equality.** With `j>=0` and `k>=1`,
   `D(1)={1} union P(2)` and `D(2)=P(2)`. If unrestricted positive bases
   satisfy `D(n)=P(m)`, then `m>=2` because the row is infinite. Their
   least and second elements give `n=m` and `2n=m²`, hence `n=m=2`.
   The uniqueness concerns equality of the two sets, not uniqueness of
   operation addresses. The segment `1,2,4` and its sum seven survive.

4. **Positive abc domains.** At `n=1`, the almost-perfect identity holds
   but `(1,s(1),1)=(1,0,1)` is not positive. For `n>=2`, the additive
   condition `1+s(n)=n` is equivalent to `sigma(n)=2n-1` and automatically
   gives pairwise coprimality. Positive Mersenne exponent `k>=1` likewise
   keeps `(1,2^k-1,2^k)` in the positive abc domain. No abc census or claim
   about unknown almost-perfect integers is certified by this repair.

The theorem title/status, source body and typing summary, MISTAKES entry,
and the two S60 navigation mirrors all express these repaired scopes. The
literal mirror replacements affect the S60 summaries, not preceding S59
material. The maintained MISTAKES purpose banner remains at its original
position; the correction follows it.

## Added exact synthesis

For `m=2^r`, `r>=1`, both `m-1,m+1` are prime exactly for `r=2`.
The case `r=1` has lower neighbour 1. If r is even, divisibility of
`2^r-1` by 3 forces that prime to equal 3; if r is odd and at least 3,
`2^r+1` is a larger multiple of 3. Thus the survivor is `m=4` and the
triple `(15,8,17)` with twin radii 3 and 5.

For the **even positive exponent** k with `2^k=m²`, this selects exactly
`k=4`. Its radical is `rad(1*15*16)=30>16`, so it is not an abc hit.
The first **Mersenne-family** hit, `1+63=64`, instead has `m=8`, triple
`(63,16,65)`, radii 7 and 9, and `Omega(63)=3`. Its radical is 42, less
than 64. The first five Mersenne values are `1,3,7,15,31`, all squarefree,
which verifies the stated first occurrence. The qualifiers “even k” and
“Mersenne-family” are essential and are present in the reviewed text.
This synthesis does not infer twin-prime infinitude from abc.

## Exact byte lineage

Independently replayed every recorded replacement from
`git show 54617edb71:<path>` with lazy fetching disabled. All old anchors
occurred exactly once at their application point. The 30 replacements
across the five files reconstructed current bytes exactly; each before
and after SHA-256 matched `LATE_CORRECTIONS.json`.

| Repaired file | Final SHA-256 |
|---|---|
| `01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md` | `94f853fea0ad204d010839afd1f64a948e43e64f3c1af20ca9ad750bc880de94` |
| `05-knowledge/results/aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md` | `872e7cf8fa6377ee73f99e3738342402d621ee166d9e1c675b985922d8160010` |
| `01-canon/MISTAKES.md` | `f79d694f3a3b2d746e5cab35de6c18e4e97b220b9e3e8a78852034dfb88529c4` |
| `00-navigation/PROBLEM-LEDGER.md` | `4adf53da0297549afb67cca93ad2c2049517965ca7a87206e4014fb128bc26e5` |
| `05-knowledge/hypotheses/INDEX.md` | `59e8d09f2446771f13a8266a027fb4f4655010125ac923cc74276229574507d6` |

The reviewed correction plan has SHA-256
`126cee7c0e2dc9c74cee98137d74c2380bc758ee4b23ceada6e5e890dc13737f`;
the reviewed `INCOMING_S60.md` has SHA-256
`fad3849c736a5a5090891b8a5275cb64a7d92b1d4cfd1c2a2ee31790323e5325`.
The exact before hashes and replacement text remain in that plan. No
full-theorem or independent-census approval should be inferred from this
bounded PASS.
