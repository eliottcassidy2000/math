# THM-2907 intrinsic-child-H2 strengthening: theorem-ready handoff

Status: exact independent strengthening/referee of the canonical hybrid
THM-2907 proof.  The landed hybrid proof is sound and should remain the
canonical primary unless its theorem text and locked counts are deliberately
updated together.

## Structural mechanism

Fix one of THM-2901's 52 pair-cap-exception marked carriers `C`, an actual
pair `L` in its finite `H4(C)` core, and the literal child
`R=C\D_L`, of mass `hR`.  Any remaining triple in a parent five-cover spans
a triangle in the inherited binary link

```text
U_R({y,z}) >= hR-q1(C).
```

The strengthening does not use the inherited threshold to decide whether a
child high core is finite.

1. Seal the exact three largest allowed singleton coverages of `R`.  If their
   sum is strictly below `hR`, subadditivity excludes every triple.
2. Otherwise compute the attained exact child cap `q1(R)`.  Apply THM-2893
   afresh with `(k,s)=(3,2)` and `B=q1(R)`.  Every triple cover contains at
   least two vertices of

   ```text
   H2(R)={w:c_R(w)>=(hR-q1(R))/2}.
   ```

3. On all 4,247 children reaching this step,

   ```text
   5hR/7-q1(R)>0,
   ```

   so THM-735 gives a finite exact `H2(R)`.  Enumerate its heavy pairs, then
   enumerate the final singleton through THM-2893's longest-component cutoff.
   Retain all three inherited binary-link inequalities.

This removes the canonical hybrid engine's optimization-dependent
finite-H2/generic split.  Failure boundary: ten H4 flags have literal
closed-danger completions; they are not empty residuals, so they require the
endpoint exit below.

## Exact finite census

```text
actual H4 flags                                      18,290
exact top-three closures                             14,043
fresh child-H2 routes                                 4,247
nonfinite fresh H2 cores                                  0
H2 routes with no completion                          4,237
boundary flags                                           10
rank-three singleton scans                       12,821,172
largest rank-three scanned label                      3,028
H2 labels scanned                                 1,922,445
H2 cutoff range                                     167..1,395
total actual H2 memberships                           13,235
maximum actual H2 size                                    33
H2 core-pair tests                                    18,649
actual heavy H2 pairs                                  9,984
terminal singleton checks                            170,275
maximum terminal cutoff                                  120
raw heavy-pair/leaf cover hits                            32
deduplicated flag/tail incidences                         16
deduplicated full tails                                    2
```

The ten flags are exactly the edges of `K5` on
`{16,18,20,26,48}`.  The 32 raw leaf hits deduplicate to 16 flag/tail
incidences and exactly two full tails:

```text
(16,18,20,24,26),
(16,18,20,26,48).
```

All 48 inherited binary-link inequalities on these incidences are strict;
their minimum margin is `155/6552`.

## Analytic boundary discharge

The two reconstructed thirteen-speed rows are

```text
2*{1,...,13},
2*{1,...,11,13,24}.
```

They are discharged, not left as obligations.  THM-369's Lean-checked master
lemma says:

> for natural `q,n` with `0<q<=n`, integer `a` coprime to `q`, and every
> speed not divisible by `q`, `t=a/q` is `n`-lonely in the non-strict
> sense.

Take `q=n=14` and `a=1` after dividing either row by its gcd 2.  Neither
normalized row contains a multiple of 14.  Therefore `t=1/14` is weakly
lonely for the normalized row, equivalently `t=1/28` for the raw row, and
the exact minimum distance is `1/14`.  Equality is allowed.  THM-407's
`M(cS)=M(S)` is a compatible scaling statement, but the explicit time
transport already proves the needed claim.

## Locked independent replay

The final intrinsic-v4 ordinary and optimized runs are byte-identical for
the result, stdout, and complete 18,290-row ledger.

```text
source SHA-256
928332c6c68054219e8008747f26fafa44c58d47b419a03f77b40afef0c09ebd

result/stdout SHA-256
93c0cd1b962102e4798ba1706fd51fce37da1c1a6b46fa29b7189ba6d377c586

ledger SHA-256
5e5671b867eef5e43359a5cd58c41b3233c2d84bf554379d4380648f9a0c9e56

semantic digest
7ab6db7e83beb02338f7a490af7cc3ffa9364036c346fe0b196f96ebac45f831

52-key full-prefix interchange digest
ae895885371a12b3986e9cc4d47f3bc1540bd519052cae4ac1b456d1c21a2d12
```

The source is preserved in the named path-limited git stash
`h4-intrinsic-v4-after-origin-promotion`; outputs are
`canonical_locked_v4.{ordinary,optimized}.*` in this directory.

## Exact branchwise compositions

The hardened THM-2904 pivot source hash is
`99f1938f264d90c2b34ec3c64566605cc8fd12520424ad2f5cd0957342202ba0`;
its regenerated full ledger hash is unchanged at
`bec35518329b5d9e6ba2c9a8c87bfb20234a0c07dc1a5c5f2babec21888d452a`.

Full `(body,gate-size,rank,apex,prefix)` branch joins give:

```text
G5                                      2,964 branches, 16 hard roots
G5 + THM-2904 pivots                    3,243 branches, 26 hard roots
G5 + pivots + H4                        3,295 branches, 26 hard roots

finite-H1 + G5                          6,054 branches, 129 hard roots
finite-H1 + G5 + pivots                 6,118 branches, 133 hard roots
finite-H1 + G5 + pivots + H4            6,170 branches, 133 hard roots
```

Thus H4 adds no whole root against the proved THM-2904 union:

```text
88 -> 88, residual 3,344.
```

The older finite-H1/G5 candidate `179` becomes `181` after the pivot
certificates, but H4 again adds no root:

```text
181 -> 181, residual 3,251.
```

This zero addition is now checked through the three-way/four-way branch
unions; it is not inferred merely from the 52 exception body names.  A future
route can consume the `CLOSED_BRANCH` interchange rows without losing the
excluded prefix.
