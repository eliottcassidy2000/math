---
id: THM-3102
title: "Projected k3 z233 exact screen and complete-cell cardinality descent"
status: >
  PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT HOSTILE AUDIT REQUIRED.
  The candidate closes all 62 rows in the pinned projected k=3, z1=233
  necessary layer.  If promoted, it updates the projected ledger to 374325
  and lowers the cap to z1<=232.  It makes no LRC(14) claim.
source: codex-thm3094-hostile-audit-2026-08-02
depends_on:
  - THM-3098-z234-final-two-body-height-and-complete-cell-closure
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z233_exact_screen_complete_cell_cardinality_descent_thm3102.py
output: 05-knowledge/results/lrc14_j7_k3_z233_exact_screen_complete_cell_cardinality_descent_thm3102.out
script_sha256: 3db29155f4c1332c91605e5589dfdc19319eb81dc309a33ff8af161abb39a036
output_sha256: f358eae14c52783b561b8b799b02fb07f2988ef6b1b9f6cd33f3030ce177727d
semantic_sha256: 6cf01affce52a5dcc67a8634da815780e3d357e72b295a7c4951211c6f12b0da
hash_basis: LF-normalized bytes
---

# THM-3102 -- projected k3 z233 exact screen and complete-cell cardinality descent

**PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT HOSTILE AUDIT REQUIRED.**

## 1. Candidate statement

In the pinned THM-2941 projected `k=3` necessary atlas, the complete
`z_1=233` layer consists of

```text
62 rows = 45 wall + 17 order.                              (1)
```

A fresh evaluation of every row through the promoted THM-3078 screen gives

```text
1,642 states = 570 crude + 1,039 exact-status + 33 residual. (2)
```

The `33` residual states occupy four wall bodies.  Every one has a positive
duplicate-permitting two-high gap.  The inherited wall condition supplies at
least one high suffix label, so every possible packet in these states has
exactly one high label.  The exact terminal then closes all `33` one-high
cases by complete-cell translated-band cardinality.  Thus the candidate
conclusion is

```text
all 62 projected k=3, z_1=233 necessary rows are empty.     (3)
```

No ledger or cap change follows before independent audit and promotion.  If
the candidate is promoted, composition with THM-3098 gives the exact disjoint
layer subtraction

```text
374387-62=374325,
projected k=3 cap: z_1<=232.                               (4)
```

The next layer is retained, not discarded.  It has

```text
z_1=232: 3 rows = 3 wall + 0 order,

(1,2,6,9,12,14),    L=3528,   H=348;
(1,5,6,9,12,14),    L=17640,  H=1738;
(1,9,10,11,12,14),  L=194040, H=19111.                    (5)
```

Here and below `H=floor(13L/132)+1` is the inherited high floor.

## 2. Pinned universe and screen direction

The companion pins the promoted source, output, and semantic hashes of
THM-3078 and THM-3098.  It reparses all `6,060` structured rows in the
THM-2941 atlas, whose LF-normalized hash is

```text
cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda. (6)
```

The independently recovered neighboring census is

```text
z_1=234: 381 = 330 wall + 51 order;
z_1=233:  62 =  45 wall + 17 order;
z_1=232:   3 =   3 wall +  0 order.                        (7)
```

The exact ordered tuple of the `62` selected rows has digest

```text
c503a3402cb2243ab586284d9d0773d07d0d0bddc95aae941c00ff7b2bf348ac. (8)
```

The logical direction of the inherited screen is load-bearing.  Its ray
quotient relaxes global suffix order, while preserving the fixed first label
`233`, distinct maximizing labels, their denominators, and the wall
at-least-one-high gate.  Therefore every actual projected packet maps into
the enumerated state set.  A crude upper exclusion or an exact Farkas
exclusion of that enlarged set is sufficient; no converse from the quotient
is used.

All `1,039` status exclusions pass the promoted legacy full-table exact
Farkas verifier.  The narrowly repaired THM-3078 zero-good direct branch is
never invoked:

```text
direct Farkas certificates:    0;
legacy Farkas certificates: 1039.                           (9)
```

The complete screen-record digest, committing to every row, state count,
status audit, and residual tuple, is

```text
2d24fe53a76095d32e7dfe0667ac1247fc01a6301a76f91dd8e5645be5c92313. (10)
```

The order subbank is independently visible inside `(2)`:

```text
17 rows, 85 states = 53 crude + 32 exact-status + 0 residual. (11)
```

Thus no order row enters the wall terminal.

## 3. Exact residual bank

The four residual rows and the complete residual-mask bank are:

| body `E` | `L` | `H` | masks | residual-tuple SHA-256 |
|---|---:|---:|---:|---|
| `(1,2,5,9,12,14)` | `17640` | `1738` | `2` | `c9bfec7a7749796239d67be94074926cd2ac89e37515ea4753ac0e63357c4258` |
| `(1,4,5,9,12,14)` | `17640` | `1738` | `2` | `1532111bfb5cafba03c14ab3e116d1ed09a1fd194e453009fe03f4198311f2d9` |
| `(1,5,6,9,12,14)` | `17640` | `1738` | `3` | `6203d6d71cedb53284763e70b2cf82fbdb17288eac6c3c5c6eff7ccbb9b40a6a` |
| `(1,5,9,11,12,14)` | `194040` | `19111` | `26` | `b4ed40d958f3692bf6bbbd8b0122c8770ee7da4dff0016f90a2743329bed4cfb` |

For completeness, the first three tuple banks themselves are

```text
(1,2,5,9,12,14):
  (2,980,1960,17640), (8,980,1960,17640);

(1,4,5,9,12,14):
  (2,980,1960,17640), (3,980,1960,17640);

(1,5,6,9,12,14):
  (2,490,980,17640), (2,980,1960,17640),
  (2,980,2940,17640).                                    (12)
```

The `26`-tuple fourth bank is serialized without abbreviation in the
companion output and committed by both `(10)` and its row digest above.

## 4. Two-high exclusion and exhaustive one-high terminal

Let `u_d` be the rigorous high-ray supremum for denominator `d`, as in
THM-3078.  For every residual mask the terminal grants repeated suffix labels
and both relevant suffix slots their independent suprema.  The resulting
duplicate-permitting two-high gap is therefore a necessary upper relaxation,
not a physical packet assertion.  The exact minima are all positive:

| body `E` | masks | minimum two-high gap | zero-high hostiles | one-high cases | coarse / exact closes | case SHA-256 |
|---|---:|---:|---:|---:|---:|---|
| `(1,2,5,9,12,14)` | `2` | `21113341/6171345180` | `0` | `2` | `0 / 2` | `deb3af778b1d9d0103a013bf9a1a02bb64f10d88c4640c1d4272eb331663c7fe` |
| `(1,4,5,9,12,14)` | `2` | `2583751/748041840` | `0` | `2` | `0 / 2` | `af4732c86fc10e42088c6f929c0024afd4438bbcd33c5735af39ee1ae8cd5cc3` |
| `(1,5,6,9,12,14)` | `3` | `233977/93059967` | `0` | `3` | `0 / 3` | `b07ec9fc8f3b7b00cd5f00f9cb4c2593661c9a78317c0b640c66f8e5a87526a4` |
| `(1,5,9,11,12,14)` | `26` | `39679901180/9741776933889` | `23` | `26` | `23 / 3` | `f6288fd9f7eceac090acee225c6a6d9f7e654f1bfa3d0ca1d80b5af62f3ef5d3` |

Hence even the enlarged state set has no packet with at least two high
labels.  The wall grammar already forces at least one high label, so exactly
one remains.  The one-high enumerator retains every finite low label,
including negative scalar contributions, enforces distinct literal labels,
and grants the unique high slot its ray supremum.  Its `33` cases are a
superset of actual one-high assignments.  Closing all `33` therefore closes
the entire residual bank.

The terminal totals are

```text
4 residual rows; 4 positive two-high gaps; 4 closed rows;
33 one-high cases = 23 coarse-cardinality + 10 exact-cardinality;
0 max-gap cases; 0 failures; 0 unit-phase rescues.          (13)
```

The complete terminal-record digest is

```text
c1447cdbc2e1aabdf94f082a5bab5ac256480d25b6f236dd870f73240a989e0e. (14)
```

The `23` zero-high scalar passes in `(1,5,9,11,12,14)` are hostile controls,
not closures.  They are excluded only because the inherited wall condition
requires a high suffix label.  They are neither silently discarded nor
counted among the `33` one-high certificates.

## 5. Complete-cell cardinality and exact implication

Fix one one-high case, its body `E`, its two finite low suffix labels, and
the high denominator `d|L`.  Let `J` be the complete `1/L` cells safe from
the body, the first drift `233`, and both low labels, and project their cell
indices modulo `d`:

```text
S_d={j mod d:j in J},                 kappa(d)=ceil(d/7). (15)
```

THM-2984 shows that an arbitrarily translated strict high-danger band has
capacity `kappa(d)`.  This is the translated count `ceil(d/7)`, not the
smaller centered quantity `beta(d)`.  Consequently

```text
|S_d|>kappa(d)                                             (16)
```

gives, at every local coordinate and for every high label on that exact
denominator ray, a complete cell missed by the high comb.  Thus

```text
P_(E,Z)=T.                                                 (17)
```

For `C=|J|`, a residue class contains at most `L/d` cells, so

```text
|S_d|>=ceil(C/(L/d)).                                      (18)
```

Inequality `(18)` proves `23` cases without constructing the support; the
other `10` supports are enumerated exactly.  Every one of the `33` cases has
strict support slack at least one in `(16)`.  An independent local audit
rebuilds each vectorized cell set with the scalar weak-endpoint predicate,
recomputes every exact residue support, and recovers the split `23+10`, the
minimum slack `1`, and record digest

```text
6db491bfaf2a70f0c12f0a690198cac8ab7a9adedfa70cb132606596ce6afa79. (19)
```

Only six distinct fixed-safe carriers occur:

| body `E` | low labels | `|J|` | least endpoint control |
|---|---|---:|---|
| `(1,2,5,9,12,14)` | `(234,243)` | `3516` | `(0,1350,14,L)` |
| `(1,4,5,9,12,14)` | `(234,243)` | `3388` | `(0,1350,14,L)` |
| `(1,5,6,9,12,14)` | `(234,243)` | `3556` | `(0,1350,14,L)` |
| `(1,5,6,9,12,14)` | `(234,324)` | `3648` | `(0,1350,14,L)` |
| `(1,5,6,9,12,14)` | `(234,402)` | `3654` | `(0,1350,14,L)` |
| `(1,5,9,11,12,14)` | `(234,243)` | `36550` | `(0,14850,14,L)` |

The zero slack in each last column occurs at the left endpoint for body
label `14`.  These are strict-open controls: the weak complete-cell endpoint
inequalities are necessary and were used.  Replacing them by strict
inequalities would delete genuine safe cells.

Finally, the direction from THM-2941 `(25g)` is

```text
literal completion  ==>  P_(E,Z) subset U_A.              (20)
```

It is not the converse and is not inferred from a scalar intermediate.
THM-1166 bounds the union for the three distinct aligned labels by

```text
mu(U_A)<=36/91<1.                                         (21)
```

Equations `(17)`, `(20)`, and `(21)` contradict one another.  This proves
the candidate conclusion `(3)` once the exact computation is independently
audited.

## 6. Complete row-digest ledger

The following ledger records all `62` rows in the exact order committed by
`(8)`.  `B` is `W` for wall and `O` for order.  The last column is the
residual key.  Its full digest is `0` for an empty residual and `A--D` for
the four banks in Section 3:

```text
0  2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d
A  c9bfec7a7749796239d67be94074926cd2ac89e37515ea4753ac0e63357c4258
B  1532111bfb5cafba03c14ab3e116d1ed09a1fd194e453009fe03f4198311f2d9
C  6203d6d71cedb53284763e70b2cf82fbdb17288eac6c3c5c6eff7ccbb9b40a6a
D  b4ed40d958f3692bf6bbbd8b0122c8770ee7da4dff0016f90a2743329bed4cfb
```

| # | B | body `E` | stage SHA-256 | residual |
|---:|:---:|---|---|:---:|
| 1 | O | `(1,2,3,8,10,12)` | `774fb5909228070027b14c101abeacd114afa02975861cb8a8493b99cec9e05b` | 0 |
| 2 | W | `(1,2,3,9,12,14)` | `37970975dfa8fa175f813d0883fc3df74f2b6aedab7c95ffa1982f5bf40130ab` | 0 |
| 3 | O | `(1,2,4,8,10,12)` | `d2d22e6f082060c4125d25eaf6cd0133d6ba837da1fe773d34b3bbe8a73446c8` | 0 |
| 4 | W | `(1,2,4,8,10,14)` | `8ee89817f5a6c70bb9c42ceb0036ad03a8182df1c59eae3cd1c1d892a5ffa767` | 0 |
| 5 | W | `(1,2,4,9,12,14)` | `3ec71c334e3c393e6c7822a34eba1268eb2cfce24a6f361bc8e53a974ca3d3f8` | 0 |
| 6 | W | `(1,2,5,6,9,12)` | `1dbd9ab2bfef2f529bc736a7cd5b5504e51fd6d49afc8daddabc8292df4ebc26` | 0 |
| 7 | W | `(1,2,5,9,12,14)` | `1829cd674cd3542204498db93abb7232f4c67d5733e49349b6d808176fa1026e` | A |
| 8 | W | `(1,2,6,7,9,12)` | `f07f25a8791dcc45e6c156545f60d098bce1582a2a5359bc9541d30a15b7aedb` | 0 |
| 9 | O | `(1,2,6,8,10,12)` | `c50aee1cb4d1cfe5ca0c7d2f47521dec0386656e8a8df0b28aa1204406db5385` | 0 |
| 10 | O | `(1,2,6,8,12,14)` | `fd16685b6fddcd780c83fa7273f0c11630cf56270c65627c6025ad2ed4f2e75b` | 0 |
| 11 | W | `(1,2,6,9,10,12)` | `c8b2e501cf4b8aef4e5b334472763d1698ed1d17907d644f0c2133357fdb1e0b` | 0 |
| 12 | W | `(1,2,6,9,12,14)` | `292d110bf9469557baf60a8646519a53291549ac1b0f656b78a126559740f8c4` | 0 |
| 13 | W | `(1,3,4,5,9,12)` | `7d17dfdb06c371b3667f6be6aad3edb7993073fb1d86bddd45e7a204b971ff00` | 0 |
| 14 | O | `(1,3,4,6,8,10)` | `e3f0cdb84fc7d17cc5c02e0badd3fbdeb202c90a47874769728eeb73c5ecdfbb` | 0 |
| 15 | O | `(1,3,4,8,10,12)` | `182a62189612134fa80fecf5d9e421907d797b8282595c2d25dada0cf386f1d7` | 0 |
| 16 | W | `(1,3,4,9,10,12)` | `40c7fa3d3cab31293eeeae509c9e99a19bee7aeac101ad1772ecc6916a569a93` | 0 |
| 17 | W | `(1,3,5,6,9,12)` | `5774b6e018f6edb0a9595d69181e0fb10d01d62e05823465de71f0d512aa08f4` | 0 |
| 18 | W | `(1,3,5,9,12,14)` | `c6800896ef572c271c5d5cb379e95af6277e890a1d7f529f69bb1ee72ad11e3e` | 0 |
| 19 | W | `(1,3,6,7,9,12)` | `9af2e9f8314b4b920d2e1a36b12ad6bdd7b3d4cd51acf8e2d3cf334386f03b0f` | 0 |
| 20 | O | `(1,3,6,8,12,14)` | `2204db1507619a7ec0c734404a73227a4bd619847a7942d6688606fa4c1d2d74` | 0 |
| 21 | W | `(1,3,6,9,10,12)` | `6d3f49e36b28eabb1afe8c36ae67ca32e1a1849360ec15729b0acdaa1a9097b2` | 0 |
| 22 | W | `(1,3,6,9,12,14)` | `b2892d47be3eb5acf0fbe5f1d3f23a9b57d0236aa0b31abbdbf98f7751e406d5` | 0 |
| 23 | W | `(1,4,5,6,9,12)` | `1de587f587194d8f9e3ae4e49e0a713aa9e012c95ab69acb1d4184957034ecbe` | 0 |
| 24 | W | `(1,4,5,6,9,14)` | `2d49ce8ba6c6845f591b436f3486aed32cf1a421a15ca2f1991bff82e586e6e4` | 0 |
| 25 | O | `(1,4,5,8,10,12)` | `3d81b869ba1a979b7e19fdc4de2ee58554ce5d063b7ce37f2b549d6ab6773b4b` | 0 |
| 26 | W | `(1,4,5,9,10,12)` | `ce02db171df41ca51d2e1d5299266cc5d5197b4cd4cb06a752c456d21122d9bb` | 0 |
| 27 | W | `(1,4,5,9,12,14)` | `756d8e52f5755f5a57f94781ef7b55fb6d07ec78b6add7980bacf09e8bb6f792` | B |
| 28 | W | `(1,4,6,7,9,12)` | `cf086dcc13dace1dd08eae1dbfe00695f6702eef95b5d09ec68aa164fc723616` | 0 |
| 29 | O | `(1,4,6,8,10,12)` | `fb112b63abbfa6197ebf0c435c7f9aedb40511a34b2c3989cae7360239c3b2eb` | 0 |
| 30 | O | `(1,4,6,8,12,14)` | `c2bc10ce77d5511385a8463556078842f6297526e8617cf9a15ee9b70f261cb6` | 0 |
| 31 | W | `(1,4,6,9,10,12)` | `545ca0759d91cf0d3cf738dad177fbd5e45c05c47649776d2a8289178622d3b7` | 0 |
| 32 | W | `(1,4,6,9,12,14)` | `6c871a919b8c0c2467f0ddd8534dc9d5a117d869ed1551246d23e6d221e74451` | 0 |
| 33 | W | `(1,4,8,10,12,14)` | `892dacc1753bd07d9756b306912d42292324e550f5fb43a8cbd39de0f3cc8630` | 0 |
| 34 | W | `(1,4,9,10,12,14)` | `392aeed00d9d22f137010076a3d5a778239853d8754b6ef16154521be9c88e6e` | 0 |
| 35 | W | `(1,5,6,7,9,12)` | `58e44526f0cd6516696954f98bb24c564c7abbdc1b778b5ff24204405d4689f1` | 0 |
| 36 | W | `(1,5,6,9,10,12)` | `a0cf8e411d3a38d6953ee159f462bdcf107c8c5c6784a46ff66c7f6da01d7701` | 0 |
| 37 | W | `(1,5,6,9,12,14)` | `c421860e1f2a57a8ae5e15a63f975455641f427f767a211bb7929a6646921479` | C |
| 38 | W | `(1,5,9,11,12,14)` | `878edf844a96494972537bb1ed3a5dd8b602a3e5c95d3477cadbfde2c1bb1fac` | D |
| 39 | O | `(1,6,7,8,12,14)` | `0398ce4c25422a5cbae722bc426fbc2152683854bb1f71a5e9eb19f86ab85c5a` | 0 |
| 40 | W | `(1,6,7,9,12,14)` | `7a27c2bb44b4ff9d3664051afb798681086f6e826233e0ff2297b65e09e4539c` | 0 |
| 41 | W | `(1,6,8,9,12,14)` | `9df7490ea32d627023d7331abd9c465dda7bb00a601ae3c4e52acb2f60f327a2` | 0 |
| 42 | W | `(1,6,9,10,12,14)` | `19a83b2fdc4477e1227671fb1d5ddba5cb9336be05d50b58553f74554487d625` | 0 |
| 43 | W | `(1,8,10,11,12,14)` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | 0 |
| 44 | O | `(2,3,4,8,10,12)` | `2db2fab32cd4527adcfe546fcecbd12f9babefc2e22e448d74a1dc05551c31f5` | 0 |
| 45 | O | `(2,3,6,8,12,14)` | `b052a40f99f50636e900aff3b888845fe3dbf15c6daf50a757781b538ed1ce1e` | 0 |
| 46 | W | `(2,3,6,9,12,14)` | `25d44b3b902a0348c741e26afa9589d7c1f16a1318118cba12a7b3c7246a0124` | 0 |
| 47 | W | `(2,4,5,9,12,14)` | `377d56aaf6a0db18fb4b3abc8e84d7414d0172a53f8fd83d69a7b692fedb8abc` | 0 |
| 48 | O | `(2,4,6,8,10,12)` | `5bf7855b567f76860fa76979140d5c49ee3c91386d6503ebb739e714bdd53a09` | 0 |
| 49 | O | `(2,4,6,8,12,14)` | `8074d508ad91d38d626db803ce1260eb63cd46357b328bf0f3ff86136c91591a` | 0 |
| 50 | W | `(2,4,6,9,12,14)` | `4448b577d5c51c5729ab6296c5ebf0b767598f13ef2bf137b70fedea3baffce7` | 0 |
| 51 | W | `(2,5,6,9,12,14)` | `01594e666eafdb945d9788f9f94828ddd632e0818e723735f7a93b0647a0ab94` | 0 |
| 52 | W | `(2,5,9,10,12,14)` | `fb5c7cb9af9942856af400dffb08cb67ba57978fad6237e77646aef7a77295b9` | 0 |
| 53 | O | `(2,6,7,8,12,14)` | `079b2c28995efee7556656f1a5ec49232eaabe8a9650c9e926a0ea33d3b057d0` | 0 |
| 54 | W | `(2,6,7,9,12,14)` | `b20b3148c0da9db5f7fd080cd80563c9d301a4cfbdb796ed02bba911a9c3dbba` | 0 |
| 55 | W | `(2,6,8,9,12,14)` | `1e6bd40e73574cbaf10674a323095e64b1147c00700f7cbdd6b1a25addbe0891` | 0 |
| 56 | W | `(2,6,9,10,12,14)` | `1fae028d5de507bc0dd5c470d607dd45edc3ff2aaa91ca379141dd48a64e3717` | 0 |
| 57 | W | `(3,4,5,9,12,14)` | `fa517d459321c294f8f07fbf208fc9902fc74f7a5f66e74c96a35ecc87bd64ee` | 0 |
| 58 | O | `(3,4,6,8,10,12)` | `91435fb68f9750ab25cd8db604f4b2f037f6005b92c0fbbbe88ef38121c429a4` | 0 |
| 59 | W | `(3,4,6,9,12,14)` | `7a63358b75aa78cae198d04169a7148085bfd3a9d9b717ac31ec4742df92e907` | 0 |
| 60 | W | `(3,4,8,10,12,14)` | `da32775c253cec5a10fd4a675c8423ac869d9b4a2eedcd75d00df3ecfe6d87c0` | 0 |
| 61 | W | `(3,5,6,9,12,14)` | `ded8df4c22545c45f1b8e3dd7520d4b8e9044885fb6fa327d706ef45305a82fa` | 0 |
| 62 | W | `(3,6,7,9,12,14)` | `968a8a940f957b310e2d5937681fb56fe0040fcf1a724a396026cc1703beddcc` | 0 |

## 7. Exact evidence and scope

Run

```text
python 04-computation/lrc14_j7_k3_z233_exact_screen_complete_cell_cardinality_descent_thm3102.py --processes 8
python -O 04-computation/lrc14_j7_k3_z233_exact_screen_complete_cell_cardinality_descent_thm3102.py --processes 8 --output <optimized-output>
```

The companion contains no truth-bearing Python `assert`.  A fresh normal run
and an optimized run both pass the frozen semantic gate; their transcripts
are LF-byte-identical to the stored output and end in
`all_exact_controls=PASS`.  The evidence hashes are

```text
script:   3db29155f4c1332c91605e5589dfdc19319eb81dc309a33ff8af161abb39a036
output:   f358eae14c52783b561b8b799b02fb07f2988ef6b1b9f6cd33f3030ce177727d
semantic: 6cf01affce52a5dcc67a8634da815780e3d357e72b295a7c4951211c6f12b0da
```

**Scope.**  This candidate acts only in the pinned projected `k=3`
necessary atlas.  It does not classify physical covers outside that
projection, say anything about arbitrary `k<=1` packets or the final rung,
or prove LRC(14).  Until independent hostile audit and explicit promotion,
THM-3098's proved ledger `374387` and cap `z_1<=233` remain canonical.

QED (candidate; audit required).
