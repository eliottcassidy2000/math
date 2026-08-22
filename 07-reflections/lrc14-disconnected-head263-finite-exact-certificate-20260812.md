# Disconnected-low finite raw head below 264

**Status: PROVED finite reduction + FINITE-EXACT + independently
hostile-audited.** This closes the finite raw-level head in the generalized
Dirichlet reduction. It does not, by itself, close the affine tails, the
disconnected-low branch, or LRC(14).

## Statement

Put

```text
Dmax/5 = 186636088362/58865718786875.
```

For every `g in {1,2,3}`, every coprime `P,Q` with

```text
P<Q<8P,  P+Q>=8,  gP<264,
```

and every one of the 2,530 upper-median physical contexts `(L,j,e,f)` with
`L<4592`, the exact reflected overlap at raw levels `(gP,gQ)` is strictly
larger than `Dmax/5`.

The task universe contains 201,377 channels and therefore
509,483,810 exact channel-context comparisons. Its ordered semantic digest
is

```text
2771f2f901f2f052952343fd77412114ae1d1d99543f42539bc28d0f0f1948af.
```

The unique global minimum is

```text
g=1, (P,Q)=(3,5), (L,j,e,f)=(336,174,12,3),
I=158/46397,
I-Dmax/5=49340690507272/210091750350356875>0.
```

At `g=2,3`, the weakest primitive channel remains `3:5`, with masses
`1258/186873` and `3254/421429`.

## Exact context and mass objects

The context generator imports the canonical THM-3352 connected-head
compiler, pinned at SHA-256
`aebfe98ab72f7eb0dc1718dfb17529e5b3f9288c6ed97d57f048bf3b29281f19`,
restricts to `L<4592`, and freezes 29 rulers, 62 `(L,j)` cells, and 2,530
ordered label contexts. Their bank digest is

```text
efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4.
```

For every frozen row it verifies

```text
L/14 <= e*j mod L,  (e*j mod L)+e <= 13L/14,
L/14 <= f*j mod L,  (f*j mod L)+f <= 13L/14.
```

Thus every nominal tooth is full and the endpoint correction in the general
THM-3352 mass engine is identically zero on this universe.

The C scanner is an integer-only port of that bulk evaluator. It orients by
the actual tooth counts `z=Lp-e` and `w=Lq-f`, computes three Euclidean floor
moments, compiles exact residue-prefix counts and sums, and evaluates every
complete tent turn plus its two tails. All moment, tent, and target
cross-products use signed `__int128`; the final reduced numerator and
denominator are narrowed only after explicit checked bounds. A full
undefined-behavior/signed-overflow sanitizer replay had no diagnostic and
generated the same output and ledger.

## Independent audits

The generated 201,377-row ledger has SHA-256

```text
1caa410d5e57d6085e730270091da3d8433f14d0bd74d98f6283f2ec8a4ca7a0.
```

It is reproducible and intentionally not stored as an 8 MB repository blob.
The independent audit checks every reported minimizer with the slower
literal-audited reference engine (SHA-256
`da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e`).
It then recomputes all 2,530 contexts for the 20 weakest and 200 seeded-spread
channels, another 556,600 exact masses. The hostile audit checks 20,241
C/reference equalities, including the equal-level orientation regression

```text
(L,j,e,p,f,q)=(168,33,3,1,8,1),  I=8/55,
```

the three weak `3:5` dilations, `1:7`, cone endpoints, maximal raw levels,
and reversed orientation. All checks pass; ordinary and optimized Python
audits are byte-identical.

An independently developed screened compiler gives a second proof on the
non-`3:5` subuniverse. It uses a downward-rounded `2^56` midpoint screen, an
unrelated C++ floor-moment port, and 99 three-route literal controls. It
proves all `509,476,220` non-`3:5` rows exceed `1/294` and finds their exact
minimum `92/7645` at `(g,P,Q;L,j,e,f)=(1,4,5;168,90,12,6)`. See
`lrc14_disconnected_low_finite_head_20260812.py`; this is an independent
audit, not a replacement universe and not an affine-tail argument.

## Reproduction

From the repository root:

```bash
python3 04-computation/lrc14_disconnected_head263_contexts_20260812.py
clang -std=c11 -Wall -Wextra -Wconversion -Wshadow -Werror -O3 -pthread \
  04-computation/lrc14_disconnected_head263_exact_scan_20260812.c \
  -o /tmp/lrc14-head263-scan
/tmp/lrc14-head263-scan \
  04-computation/lrc14_disconnected_head263_contexts_20260812.txt \
  /tmp/lrc14-head263.ledger 8
python3 04-computation/lrc14_disconnected_head263_independent_audit_20260812.py \
  /tmp/lrc14-head263.ledger
python3 04-computation/lrc14_disconnected_head263_hostile_audit_20260812.py \
  /tmp/lrc14-head263-scan
```

Frozen SHA-256 values:

| Artifact | SHA-256 |
|---|---|
| C scanner | `498ef114fbf7b3d54e62de4556cd0f669b0300f882d61d16c39f1566f5efb23f` |
| context generator | `468a88781de94cb6d0d49d500371bb234b1ddb2615e8e8a8c6203c827d1ca298` |
| context bank | `efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4` |
| scan output | `44befca1fff5c832d52658b221f070930d683f621dbfe8d1ce612911c72d9cb5` |
| independent audit | `f82d7bde225253f4a3262bc94a2dac59ff6639c80a18da62a7ffd5781531333d` |
| independent output | `ffb7670550505cc81ebad0e88fbbb2104649fb91fb9fd38f77eb66a14409c08f` |
| hostile audit | `4b6ae438ef91ae88358cf945d441a51a4465e133033b65caf6952b51a052920f` |
| hostile output | `099363056fef772749f4eb9e077d8e56b455bb7a5fe66da65c5dfc4cf009410b` |

## Consequence and remaining boundary

This discharges the `p<264` gate in the generalized Dirichlet reduction.
The live object is its affine tail. The refined turn-band argument and exact
bridge computations are being audited separately; in particular, the
equality face `T_infinity=1` must be retained rather than silently absorbed
into the strict `T_infinity<1` chamber.
