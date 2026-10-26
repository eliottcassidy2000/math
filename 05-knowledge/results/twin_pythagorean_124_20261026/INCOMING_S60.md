# Incoming S60: bounded domain and quantifier corrections

**BOUNDED REVIEW / REPAIRED. Full THM-4805 audit remains pending.**
The incoming revision `54617edb71` supplies the actual
[THM-4805 / aliquot-sources-of-2-5-7-21](../../../01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md)
and [source note](../aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md).
The root packet's frozen reference to the S60 reservation describes its
historical research base; the incoming theorem is now populated. This
continuation does not adopt its unrelated unaudited claims as dependencies.

The bounded review found and repaired exactly these four mathematical issues
in the theorem, its title/status, the source synthesis/body/typing summary,
the mistakes ledger, and the two S60 summaries in PROBLEM-LEDGER and the
hypothesis INDEX (the preceding S59 text is unchanged):

1. **A'(i), unit and prime fibres.** Targets 0 and 1 cannot obey the three
   composite-source-shape classification: `s(1)=0`, and every prime has
   aliquot sum 1. Restrict the classifier to 2<=v<=14 and state those two
   fibres separately. The regime counts proper divisors, `t=tau-1<=3`,
   not tau<=3; source 6 already has tau=4. The source's first-handful prose
   is correspondingly restricted to targets 2..5.
2. **A'(ii), finite versus conditional.** For positive exponents 1<=k<=26,
   the single-source iff k=2,3 is CITED FINITE-EXACT from S57's list below
   10^8. That range was not independently recomputed by this review.
   Extending to all k>=1 uses the explicit distinct-prime Goldbach
   antecedent. At k=1 the target 1 has infinitely many prime sources;
   exponent 0 is excluded because target 0 has the unique source 1.
   The former unconditional global conclusion is unsupported by its proof,
   not disproved by a counterexample.
3. **C'(5), the missing unit.** The inherited and implemented domains are
   j>=0 for doubling and k>=1 for powers. Thus D1={1} union P2 and D2=P2.
   If D(n)=P(m), infinitude gives m>=2, and least/second elements give
   n=m and 2n=m²; hence n=m=2. This is the unique exact equality among
   unrestricted positive row bases. The 1,2,4 segment of D1 remains valid.
   The source now says 2 is the only prime on D1, not on an arbitrary row.
4. **D', positive abc.** Require n>=2 for the positive additive triple
   (1,s(n),n), and k>=1 for the Mersenne triple. The n=1 case is (1,0,1),
   despite satisfying sigma(1)=2*1-1. No assertion about the separate abc
   census or unknown almost-perfect numbers was re-audited or altered.

The compact proofs and minimal witnesses were checked directly. S57's cited
finite range and all unrelated producer computations remain at their inherited
scope. The remaining B', C'(1)–(4), numerical censuses, literature claims and
open hypotheses have not received a full new audit here. Existing broader
pending-audit labels remain in place.

A concrete bridge to this packet also survives the repair. For a power-of-two
centre `m=2^r`, `r>=1`, the two neighbours `m-1,m+1` are both prime if and only
if `r=2`, so `m=4`. Indeed, `r=1` has lower neighbour 1. For even r,
`3 | (2^r-1)`, so primality forces `2^r-1=3`; for odd r>=3,
`3 | (2^r+1)` and `2^r+1>3`, excluding primality. The survivor gives twin radii 3 and 5
and the primitive triangle `(15,8,17)`. Thus among even positive k with
`2^k=m²`, the Mersenne identity `1+(2^k-1)=2^k` is the centre identity of
a twin-radius triangle exactly when `k=4`. Its radical is
`rad(1*15*16)=30>16`, so it is not an abc hit.

By contrast, the first Mersenne abc hit `1+63=64` has centre 8 and the primitive
triangle `(63,16,65)`, with inradius 7 and corresponding exradius 9. Its radical
is `rad(1*63*64)=42<64`, but 9 is composite and `Omega(63)=3`. This is exactly
the inherited multiplicity obstruction to the twin-radius image. The literal
factorizations `15=3*5`, `16=2^4`, `63=3^2*7`, `64=2^6` verify the two radical
comparisons. For k=1,...,5 the Mersenne values 1,3,7,15,31 are squarefree,
so the repaired abc criterion also verifies directly that k=6 is the first
hit. This bridge is an elementary consequence, not an assertion that abc
implies infinitely many twin primes.

[LATE_CORRECTIONS.json](LATE_CORRECTIONS.json) stores source revision,
repository-relative paths, before/after SHA-256, and every exact old/new
replacement. Applying each replacement once in order to `git show
54617edb71:path` reconstructs the current file exactly. MISTAKES has one
insertion immediately below its unchanged purpose banner. No root README,
connection board, diagram, verifier, manifest or other audit file was changed.

| Repaired path | Before SHA-256 | After SHA-256 |
|---|---|---|
| `01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md` | `0c2a4f22368dea54475e1f31addbb1bc959fa8f6a6afa42b97f7289b3d689e0f` | `94f853fea0ad204d010839afd1f64a948e43e64f3c1af20ca9ad750bc880de94` |
| `05-knowledge/results/aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md` | `2bb8e01d22c3905200b324a5f0e78105a3d87cac7665b21882acdb765a22912c` | `872e7cf8fa6377ee73f99e3738342402d621ee166d9e1c675b985922d8160010` |
| `01-canon/MISTAKES.md` | `d6b9d9569e6e1c7e5e61f957f6ea355ea9c9c9f28936da84d066967cd9568af1` | `f79d694f3a3b2d746e5cab35de6c18e4e97b220b9e3e8a78852034dfb88529c4` |
| `00-navigation/PROBLEM-LEDGER.md` | `09c33a0d568f5426bc21ef0e07a37ca6bb46df08069b659071cac1cd0cf72298` | `4adf53da0297549afb67cca93ad2c2049517965ca7a87206e4014fb128bc26e5` |
| `05-knowledge/hypotheses/INDEX.md` | `9d3e64ef58cb45c9b8dbee6c7d7e85e208c7f6a003349b6b77e7d26aff24c51a` | `59e8d09f2446771f13a8266a027fb4f4655010125ac923cc74276229574507d6` |
