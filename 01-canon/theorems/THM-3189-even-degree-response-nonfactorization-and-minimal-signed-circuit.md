---
id: THM-3189
title: "Even-degree response nonfactorization and minimal signed circuit"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On THM-3184's complete 1,820-state depth-seven bank, none of its ten
  selected degree-N+2 upset rows is a linear or affine function of the
  complete degree-N partition response.  Exact ranks are certified modulo
  1,000,003 and lifted to Q by the zero-mass upper bound.  A circuit-minimal
  twelve-state affine hostile splits into two rational probability laws with
  identical complete degree-six response and different selected degree-eight
  upset response.  No nonlinear, Markov, or enriched-carrier obstruction is
  asserted.
source: root/multiscale-newton-flag/product-gamma-width3/2026-08-02
audit: >
  The integer/modular companion pins THM-3184 and its two transitive source
  helpers, rebuilds all 1,820 states and 7,280 source balance identities,
  verifies the four maximal zero-mass ranks, all ten linear and affine rank
  jumps, response injectivity at each source degree, and the exact primitive
  twelve-state circuit including every proper deletion.  Normal, optimized,
  and stored replay agree exactly.  An independent immutable audit rederived
  the modular-to-rational rank bounds, affine factorization criterion,
  probability split, target separation, and circuit minimality; replayed the
  full 1,820-state normal computation; matched both hashes; and accepted the
  nonlinear/Markov scope boundary.
depends_on:
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3184-depth-seven-degree-fourteen-farkas-death
related:
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
script: 04-computation/gmc_even_degree_response_nonfactorization_thm3189.py
output: 05-knowledge/results/gmc_even_degree_response_nonfactorization_thm3189.out
script_sha256: 6d86e3f56c09075c5bde6fd4b56c0ca8d3f122dd60d6b7a583ed30e793628dfa
output_sha256: 7e54730c27eb337e23f6f115c13ae2e05a5d079afc874ab0efa0cc0c853841bc
hash_basis: LF-normalized bytes
---

# THM-3189 -- even-degree response nonfactorization and minimal signed circuit

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3184's separator uses degrees `8,10,12,14`.  A tempting explanation is
that a single two-step response operator generates each higher-degree facet
from the complete response two degrees below.  This is false in the strongest
linear sense visible on the state bank: every one of the ten selected higher
rows adds a new rational covector, even after adjoining probability
normalization to the lower response.

## 1. Complete response matrices

Retain THM-3184's support-`(1,3)`, bank-`I2` physical state bank

```text
S=S_<=7,                         |S|=1,820.              (1)
```

For `N in {6,8,10,12}`, let

```text
A_N(mu,sigma)=G_N^sigma(mu),
mu in Par(N),                    sigma in S.             (2)
```

Thus `A_N` has `p(N)` rows and `1,820` columns.  Every column has total mass
zero, so

```text
rank_Q A_N <= p(N)-1.                                      (3)
```

Adjoin the normalization row

```text
A~_N=[A_N; 1_S].                                           (4)
```

For an upset `U` of `Par(N+2)`, write its response row as

```text
t_(N+2,U)(sigma)=sum_(mu in U)G_(N+2)^sigma(mu).           (5)
```

## 2. Ten new covectors

The target rows are exactly THM-3184's ten selected upsets, grouped by their
putative source degree.  Generators are the full minimal antichains.

| source -> target | `|U|` | minimal generators of `U` |
|---|---:|---|
| `6 -> 8` | 21 | `(2,1^6)` |
| `8 -> 10` | 40 | `(3,1^7),(2,2,1^6)` |
| `10 -> 12` | 76 | `(2,1^10)` |
| `10 -> 12` | 74 | `(2,2,1^8)` |
| `12 -> 14` | 130 | `(3,2,1^9),(2,2,2,1^8)` |
| `12 -> 14` | 128 | `(4,1^10),(3,2,1^9),(2^6,1^2)` |
| `12 -> 14` | 121 | `(7,1^7),(3,3,2,1^6),(2^4,1^6)` |
| `12 -> 14` | 132 | `(2,2,1^10)` |
| `12 -> 14` | 129 | `(5,1^9),(3,3,1^8),(2^3,1^8)` |
| `12 -> 14` | 127 | `(5,1^9),(4,2,1^8),(3,3,1^8),(2^4,1^6)` |

Let `q=1,000,003`.  Exact row reduction over `F_q` gives

| `N` | `p(N)` | `rank A_N` | `rank A~_N` |
|---:|---:|---:|---:|
| 6 | 11 | 10 | 11 |
| 8 | 22 | 21 | 22 |
| 10 | 42 | 41 | 42 |
| 12 | 77 | 76 | 77 |

For **each** target row `t` in the preceding table,

```text
rank_Fq [A_N;t]       =p(N),
rank_Fq [A~_N;t]      =p(N)+1.                           (6)
```

Consequently

```text
t notin Row_Q(A_N),
t notin Row_Q(A~_N).                                     (7)
```

In particular no selected degree-`N+2` facet is a linear or affine function
of the complete degree-`N` response on `S`.

### Proof

Reduction modulo a prime cannot increase rational rank.  The first modular
rank in the table attains the rational upper bound `(3)`, so

```text
rank_Q A_N=p(N)-1.                                       (8)
```

The normalization row raises the modular rank once, while its row count gives
the matching rational upper bound; hence `rank_Q A~_N=p(N)`.  Finally the two
augmented modular ranks in `(6)` attain the respective next row-count upper
bounds.  They therefore also hold over `Q`, proving `(7)`.  QED.

This is stronger than failure relative to the other nine certificate rows:
the target remains new after receiving **every** partition coordinate two
degrees below.

## 3. A circuit-minimal probability hostile

The first affine failure admits a compact exact witness.  In the state order

```text
sigma=((1),(2),(3),(4),(5),(6),(7),(8),
       (1,1),(1,2),(1,3),(1,4)),                         (9)
```

put

```text
c=(-1047277306414574923871547479704943667621575439892091,
    7504006088796364741516689429373777425450644438777969,
   -23076608092958877009276603275374222862197557645092271,
    39446840155066170582433740489847512927510928685260765,
   -40481016835146433236610937331785735797691389653391585,
    24940417689554928500780094388795198827192237041614851,
   -8542034917490033954128648539285677801015113585173749,
    1254725791744637777888128377808698770041845615340271,
    4534007100750039738438417212526753888724041742080,
   -8125859585348642127289323848561206793955708255360,
    6660162235070956799665346800833144498192277589760,
   -2120882902654833141730499839406513262980068520640).
                                                               (10)
```

This integer vector is primitive, has no zero coordinate, and satisfies

```text
sum_j c_j=0,                 sum_j c_j A_6(.,sigma_j)=0. (11)
```

Its positive and negative masses agree exactly:

```text
m=73157183894497922599156756449838547848582572100325696.
                                                               (12)
```

Let `lambda_+` put mass `c_j/m` on positive coordinates and let `lambda_-`
put mass `-c_j/m` on negative coordinates.  They are rational probability
laws with

```text
A_6 lambda_+=A_6 lambda_-.                                (13)
```

For the degree-eight upset generated by `(2,1^6)`, however,

```text
t lambda_+-t lambda_-
 =2277546515253459787947043677351991747527514706519412855776760
  /1143080998351530040611824319528727310134102689067589
 >0.                                                       (14)
```

The witness is circuit-minimal in the affine response matroid: after deleting
any one of the twelve columns in `(9)`, the remaining eleven columns of
`A~_6` are independent modulo `q`, hence over `Q`.  Thus no proper sub-support
of `(10)` carries an affine dependence.

This probability pair is more informative than an arbitrary null vector.
It shows that complete lower-degree response data can agree exactly on two
honest laws while the selected higher-degree facet separates them.

## 4. What the rank jump means

The degree axis is not a two-step linear functor on these response
coordinates.  Each selected `N+2` covector supplies information transverse to
the entire degree-`N` response space.  In the holotopy language, there is no
edge label determined solely by `A_N` which transports these ten facets; an
additional coordinate or a state-level lift is necessary.  This is a linear
rank statement, not a claim of a topological loop.

There is a formal resemblance to THM-3183's transverse exterior coordinate:
both mechanisms expose data lost by a lower-dimensional projection.  No map
between the factorial Gauss--Manin lattice and the partition-response bank is
known, so that resemblance is not used in the proof.

## 5. Exact evidence and scope

Run

```text
python 04-computation/gmc_even_degree_response_nonfactorization_thm3189.py
python -O 04-computation/gmc_even_degree_response_nonfactorization_thm3189.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
integer response rows and modular row reduction only.  It pins both transitive
source helpers and THM-3184's exact script/output, independently regenerates
the target upsets, and verifies all twelve proper circuit deletions.

The individual complete response columns at each of degrees `6,8,10,12` are
pairwise distinct on `S`; hence a hostile consisting of two individual states
cannot exist.  The affine circuit is the correct mixture-level replacement.

This theorem is restricted to the fixed support, invariant bank, state bank,
and ten named targets above.  It does not rule out a nonlinear map, a Markov
lift using enriched or state-level data, a map using the physical state itself,
or a relation involving different facets.  It proves no positivity-cone
feasibility or infeasibility beyond THM-3184, no Gaussian-moment counterexample,
and no LRC consequence.

QED.
