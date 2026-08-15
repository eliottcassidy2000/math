---
id: THM-3451
title: "Target-free prime-quarter literal half-twist cap-seven classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  After excluding
  every previously proved lower supported divisor, strict literal half-twist
  covers on Q=4p, p an odd prime, by at most seven distinct transverse owners
  have exact support p=17,37.  No composite-even, arbitrary-time, decrement,
  or LRC(14) consequence is claimed.
source: root-prime-quarter-half-twist-rank-seven typed four-sheet closure, 2026-08-15
audit: >
  independent immutable full-package audit reconstructed the type-size,
  central-fibre, profile, joint-period, line/coset, large-a plateau,
  positive-weight, raw-zero, rational-descent, 96-class graph, finite-DFS,
  witness, multiplicity, legacy-scope, and boundary-stitch gates; clean-room
  direct-mask and graph controls, dependency/hash, AST/security,
  normal/optimized/stored, ID/routing, scope, diff, and documentation gates
  clean
depends_on:
  - THM-3434-seventeen-fibre-two-sided-mass-closure
  - THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists
  - THM-3445-prime-even-half-twist-cap-seven-classification
related:
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3421-prime-half-twist-rank-seven-classification
script: 04-computation/lrc_prime_quarter_half_twist_cap7_classification_thm3451.py
output: 05-knowledge/results/lrc_prime_quarter_half_twist_cap7_classification_thm3451.out
script_sha256: edda04033f24aa4275e19b05b56857db9570d9e78e02912c8a167445591b0422
output_sha256: 00af2309eb53f3e8937d2604a3db6ea1a7a57c2a7f39f7873915655a2d511aa6
semantic_sha256: b455ca48ce8454d754b1a7d4f816ecd76bb431a5b9202757152badd5e0993482
hash_basis: LF-normalized bytes
---

# THM-3451 -- target-free prime-quarter literal half-twist cap-seven classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  The analytic proof,
exact companion, predecessor interfaces, direct-mask controls, finite search,
and boundary stitch passed a separate immutable full-package audit.

## 1. Exact scoped statement

Let `p` be an odd prime, put `Q=4p`, and use the strict literal masks

```text
B_r={ell in Z/(4p)Z: ||r(2ell+1)/(8p)||<1/14}.          (1)
```

Choose the sign representatives

```text
1<=r<4p.                                                (2)
```

The transverse owner condition excludes the zero class.  The three
transverse representatives `p,2p,3p` have empty masks and may be omitted
without loss.  All selected owners are distinct.

Call the row **post-THM-3445 target-free** when

```text
p not in {3,5,7,11,13,19,23,29}.                        (3)
```

These are exactly the prime-quarter rows containing one of the previously
proved lower supported divisors

```text
8,9,10,11,12,13,14,15,23,25,29,38,51.                  (4)
```

The candidate classification is

```text
some at-most-seven family covers Z/(4p)Z
                    iff p is 17 or 37.                  (5)
```

The complete positive-profile controls are as follows.  The profile records
the counts `(a,b,u,v)` of the four owner types `(A,B,E_o,E_e)`, defined in
Section 3, and `Omega=sum_r |B_r|-4p`.

| `p` | profile | `Omega` | direct witness |
|---:|---:|---:|---|
| `17` | `(4,0,0,3)` | `8` | `(8,11,23,24,45,56,57)` |
| `17` | `(4,0,1,2)` | `4` | `(1,24,32,33,35,36,67)` |
| `17` | `(4,0,2,1)` | `0` | `(8,11,12,23,28,45,57)` |
| `17` | `(4,2,0,1)` | `0` | `(8,13,21,30,38,47,55)` |
| `37` | `(4,0,2,1)` | `8` | `(8,33,41,100,107,115,140)` |

There are two target-free conventions in the predecessor literature.  The
older THM-3435 `A_0` convention omits only

```text
p in {3,5,11,13,23,29}.                                (6)
```

In that legacy universe, THM-3445 supplies the additional proper-period
pullbacks at `p=7,19`.  Thus the corresponding legacy support is exactly
`{7,17,19,37}`.  Statement `(5)` always means the strengthened convention
`(3)`.

## 2. Inheritance and the joint-period gate

The closest proved mechanism is
`01-canon/theorems/THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md`.
Its four-sheet chart supplies the owner gates and the affine fibre sidecar.
The odd proper-period input is
`01-canon/theorems/THM-3434-seventeen-fibre-two-sided-mass-closure.md`, and
the even proper-period input is
`01-canon/theorems/THM-3445-prime-even-half-twist-cap-seven-classification.md`.

For an owner define its quotient order

```text
m_Q(r)=Q/gcd(Q,r),
L=lcm_r m_Q(r).                                         (7)
```

The standard period descent writes every mask as the pullback of a literal
mask on `L`; transversality, distinctness after deduplication, and coverage
are preserved.  Since `L|4p`, a proper joint period lies in

```text
1,2,4,p,2p.                                             (8)
```

The periods `1,2,4` are empty by direct strict-mask evaluation.  If `L=p`,
THM-3434 would force `p` to be one of `11,13,23,29`, all excluded by `(3)`.
If `L=2p`, THM-3445 would force `p=7` or `19` after its own lower-supported
rows are removed, again excluded by `(3)`.  Therefore every hypothetical
cover in `(5)` has joint period `4p`, and all joint-period gates of THM-3435
apply.

This proper-period step explains why `p=7,19` are exclusions rather than new
negative computations.  It also prevents a full-period finite search from
silently rediscovering a lower theorem.

## 3. Four owner types, exact sizes, and overlap budgets

Write each nonempty representative in one of the forms

```text
A:   r=q,   q odd,                 f=0, alpha=8;
B:   r=2q,  q odd,                 f=1, alpha=4;
E_o: r=4q,  q odd,                 f=2, alpha=2;
E_e: r=4q,  q even,                f=2, alpha=2.         (9)
```

Here `alpha=2^(3-f)`.  On the THM-3435 fibre over `Z/pZ`, an `A` owner
selects one of four sheets, a `B` owner selects one two-sheet coset, and both
`E` types select all four sheets of every active base fibre.  The parity of
`q` is still essential: at the central base fibre

```text
j_c=(p-1)/2,                                             (10)
```

every `A,B,E_o` mask is empty, whereas every `E_e` mask contains all four
central sheets.  This is the coordinate lost by the naive projection to the
odd base.

Write `p=14k+s`, where `s` is one of `1,3,5,9,11,13`.  Exact strict odd-word
counting gives

| `s` | `|A|` | `|B|` | `|E_o|` | `|E_e|` |
|---:|---:|---:|---:|---:|
| `1` | `8k` | `8k` | `8k` | `8k+4` |
| `3` | `8k+2` | `8k` | `8k` | `8k+4` |
| `5` | `8k+2` | `8k+4` | `8k` | `8k+4` |
| `9` | `8k+6` | `8k+4` | `8k+8` | `8k+4` |
| `11` | `8k+6` | `8k+8` | `8k+8` | `8k+4` |
| `13` | `8k+8` | `8k+8` | `8k+8` | `8k+4` |

For a target-free joint-period cover, THM-3435 gives exactly seven owners and

```text
a>=2,             v>=1,             a+2b>=4.            (11)
```

The last inequality is the two-sheet capacity gate.  The central fibre in
`(10)` costs

```text
Omega>=4(v-1).                                          (12)
```

Combining the size table with `(11)--(12)` yields the complete profile ledger.
The `raw` column imposes `(11)` and `Omega>=0`; `central` also imposes `(12)`.

| `s` | exact `Omega` | raw profiles | central profiles | `Omega_max` | maximizers |
|---:|---|---:|---:|---:|---|
| `1` | `4(v-1)` | `26` | `26` | `12` | `(2,1,0,4)` |
| `3` | `2a+4v-12` | `19` | `10` | `8` | `(4,0,0,3)` |
| `5` | `2a+4b+4v-20` | `13` | `5` | `4` | `(2,3,0,2),(2,4,0,1)` |
| `9` | `20-2a-4b-4v` | `19` | `13` | `8` | `(2,1,3,1),(4,0,2,1)` |
| `11` | `12-2a-4v` | `13` | `10` | `4` | `(2,b,4-b,1)`, `1<=b<=4` |
| `13` | `4-4v` | `13` | `13` | `0` | all survivors have `v=1` |

In particular,

```text
Omega<=12.                                               (13)
```

If a covered sheet has multiplicity `mu`, its contribution to the total
excess is `mu-1`.  Hence every pair in a cover satisfies

```text
w(r,t):=|B_r intersect B_t|<=Omega<=12.                 (14)
```

The central invoice is not merely a search optimization.  Over the `77`
bounded target-free primes it reduces `1,312` raw mass-feasible profiles to
`959` genuinely central-admissible profiles.

## 4. The fully typed line-and-coset compiler

This section treats every pair except `E_e/E_e`.  Order the pair so the first
owner `X` is non-`E_e` and no deeper than the second owner `Y`.  Thus its
effective coefficient is odd.  The nine ordered sectors are

```text
AA, AB, AE_o, AE_e, BB, BE_o, BE_e, E_oE_o, E_oE_e.    (15)
```

This ordering is load-bearing: an `E_e` effective coefficient is even and
must never be silently used as the odd normalized root.

Put `N=2ell+1 modulo 8p`.  Multiplication by an odd unit permutes these sheet
words, so normalize the first effective coefficient to one.  Let `x` be its
centered odd coordinate modulo `alpha_X p`, and let `y` be the centered target
coordinate modulo `alpha_Y p`.  The mask conditions are

```text
|x|<alpha_X p/14,             |y|<alpha_Y p/14.         (16)
```

The fourteen-gap argument produces coprime positive integers

```text
1<=b_0<=13,       1<=a_0<=floor((p-1)/14),
b_0 q_Y-epsilon a_0=m p.                                (17)
```

Put

```text
c_A=4, c_B=2, c_Eo=c_Ee=1,
eta_A=eta_B=eta_Eo=1, eta_Ee=0.                         (18)
```

The coefficient range and parity conditions in `(17)` are exactly

```text
0<=m<=c_Y b_0,
gcd(m,b_0)=1,
m = b_0 eta_Y-a_0 (mod 2).                              (19)
```

Every intersection point lies on a line

```text
b_0 y-epsilon a_0 x=Jp,                                 (20)
```

but the affine lift in `(17)` also forces the exact sheet-coset congruence

```text
J=m x (mod b_0 alpha_Y).                                (21)
```

For a fixed `J`, count odd residue classes `x` modulo `2b_0 alpha_Y` that
satisfy `(21)`; call this count `n_J`.  The compiler table is

| sector | `alpha_X` | `alpha_Y` | congruence modulus | odd-word period | literal factor |
|---|---:|---:|---:|---:|---:|
| `AA` | `8` | `8` | `8b_0` | `16b_0` | `1` |
| `AB` | `8` | `4` | `4b_0` | `8b_0` | `1` |
| `AE_o,AE_e` | `8` | `2` | `2b_0` | `4b_0` | `1` |
| `BB` | `4` | `4` | `4b_0` | `8b_0` | `2` |
| `BE_o,BE_e` | `4` | `2` | `2b_0` | `4b_0` | `2` |
| `E_oE_o,E_oE_e` | `2` | `2` | `2b_0` | `4b_0` | `4` |

The centered line section has scaled length

```text
G_J=max(0,
  min(alpha_X a_0,14J+alpha_Y b_0)
 -max(-alpha_X a_0,14J-alpha_Y b_0)),
g_J=p G_J/(14a_0).                                      (22)
```

An open interval of length `g` contains at least `ceil(g/d)-1` points of any
fixed arithmetic progression of step `d`.  Consequently the literal pair
weight has the exact phase-free lower bound

```text
L_XY=2^f_X sum_J n_J max(0,
       floor((pG_J-1)/(28a_0 b_0 alpha_Y))).             (23)
```

It is monotone in `p` for fixed `(X,Y,a_0,b_0,m)`.  A row is disjoint exactly
when no line has both `G_J>0` and `n_J>0`.  Formula `(21)`, not the base ratio
alone, decides this.

## 5. Uniform positive-weight floor for `p>=449`

The companion evaluates `(22)--(23)` by exact integer arithmetic.  For
`1<=a_0<=20`, it uses `p=449`:

| sector | rows | active | no-line | active lines | residue classes | least active lower bound |
|---|---:|---:|---:|---:|---:|---:|
| `AA` | `2030` | `2023` | `7` | `17668` | `45040` | `16` |
| `AB` | `1020` | `1015` | `5` | `7148` | `17276` | `16` |
| `AE_o` | `515` | `514` | `1` | `3694` | `7388` | `16` |
| `AE_e` | `515` | `512` | `3` | `3683` | `7366` | `16` |
| `BB` | `1020` | `1011` | `9` | `4646` | `11228` | `16` |
| `BE_o` | `515` | `510` | `5` | `2170` | `4340` | `16` |
| `BE_e` | `515` | `507` | `8` | `2165` | `4330` | `16` |
| `E_oE_o` | `515` | `499` | `16` | `1368` | `2736` | `16` |
| `E_oE_e` | `515` | `500` | `15` | `1415` | `2830` | `16` |

For `21<=a_0<120`, evaluate at the smaller admissible boundary

```text
p_0=max(383,14a_0+1).                                   (24)
```

There are `34,460` rows, none disjoint.  Their sector counts and least bounds
are

```text
AA      9768 / 16,       AB      4909 / 20,
AE_o    2479 / 22,       AE_e    2479 / 22,
BB      4909 / 16,       BE_o    2479 / 20,
BE_e    2479 / 20,       E_oE_o  2479 / 16,
E_oE_e  2479 / 16.                                      (25)
```

For `a_0>=120`, the plateau

```text
|14J|<=alpha_X a_0-alpha_Y b_0                          (26)
```

has `G_J=2alpha_Y b_0`.  The admissible `J` residues in `(21)` have cyclic
gap at most `alpha_Y`.  Checking the two starting parities `a_0=120,121`
gives at least

```text
15,32,66,67,15,32,32,15,15                             (27)
```

plateau lines in the nine sectors of `(15)`.  Because `(17)` implies
`p>=14a_0+1`, every such line contributes at least one progression point.
After the literal factors, the least sector bounds are

```text
15,32,66,67,30,64,64,60,60.                            (28)
```

It remains to treat `E_e/E_e`.  Put `h=floor((p-1)/14)` and normalize the
effective ratio to `lambda`.  The base intersection lattice is

```text
L_lambda={(x,y) in Z^2:y=lambda x (mod p)},
det L_lambda=p.                                         (29)
```

For its successive `l_infinity` minima `sigma_1<=sigma_2`, Minkowski gives
`sigma_1 sigma_2<=p`.  If `2sigma_1<=h`, the five points
`0,+/-v,+/-2v` lie in the strict box.  Otherwise

```text
sigma_2<=p/sigma_1<2p/h<=h,                             (30)
```

because `h^2>=2p` from `p=409` onward (`h=29` and `h^2=841` at the boundary).
Then `0,+/-v_1,+/-v_2` are five distinct base points.  Each base point has all
four literal sheets, so

```text
w(E_e,E_e)>=20.                                         (31)
```

Combining `(23)--(31)` proves the tail dichotomy

```text
p>=449: every pair is disjoint or has weight at least 15>12.   (32)
```

## 6. Raw zero atlas and the modulo-420 graph

The no-line rows in Section 5 form a finite exact atlas.  Their counts in the
nine sectors are

```text
7,5,1,3,9,5,8,16,15,          total 69.                (33)
```

Every row has raw height `a_0+b_0<=7`.  Normalize an `A` owner to coefficient
one.  Its complete root bank, with entries `(a_0,b_0,m)`, is

| target type | no-line tuples |
|---|---|
| `A` | `(1,1,2),(1,1,4),(1,3,4),(1,5,4),(1,5,12),(3,1,4),(5,1,4)` |
| `B` | `(1,1,2),(1,3,2),(1,5,2),(1,5,6),(3,1,2)` |
| `E_o` | `(1,2,1)` |
| `E_e` | `(1,1,1),(1,2,1),(1,3,1)` |

For each tuple, the signs in `(17)` are instantiated only when the coefficient
is integral, in range, and has parity `eta_Y`.  Every prime `p>=449` therefore
has exactly `15` root-zero neighbors, with type counts

```text
7A+4B+E_o+3E_e.                                         (34)
```

This graph is finite for an arithmetic reason, not by extrapolation.  A root
candidate and a candidate-pair zero relation use rational numerators and
denominators at most six.  If a modular equality did not already hold over
the rationals, clearing the three denominators would give a nonzero integer
of absolute value at most

```text
2*6^3=432<p.                                            (35)
```

Thus modular equality implies rational equality for `p>=449`.  Candidate
existence, type, and affine lift depend only on `p modulo 60`.  Combining this
with the profile/excess class modulo `14` gives the safe full period

```text
lcm(60,14)=420.                                         (36)
```

The companion enumerates all `96` unit classes modulo `420`.  It evaluates
direct literal masks at two primes in every class: first primes range from
`449` to `3319`, and second primes from `877` to `3739`.  It finds `16`
labelled graph signatures, each occurring in six classes, and in every class

```text
root degree=15,             rooted clique number=6.     (37)
```

The rooted maximum is unique and has the coefficient presentation

```text
{1,4p-1,2p-2,2p-1,2p+1,2p+2}.                          (38)
```

It consists of three `Q`-complementary pairs and has type `A^4E_oE_e`.
Direct full-root controls at `p=449,1009` recover exactly the bank `(34)` and
minimum positive weights `16,36`.

There is also an independent root-`E_e` hostile sidecar.  Direct controls at
`p=409,449,1009` have exactly `35` zero neighbors, of types
`12A+12B+11E_o`, with `109` edges among the neighbors and rooted clique
number six.  Their two maximum presentations are

```text
{8,p-4,p+4,3p-4,3p+4,4p-8},
{8,2p-16,2p-8,2p+8,2p+16,4p-8}.                        (39)
```

This sidecar is a hostile control rather than a dependency of the tail proof.

## 7. Tail contradiction and exact finite stitch

Suppose `p>=449` supported a cover.  By `(14)`, every pair weight is at most
`12`; by `(32)`, every pair must therefore be disjoint.  Gate `(11)` supplies
an `A` owner.  Normalize it to one.  The seven owners would form a rooted
seven-clique in the exact graph of Section 6, contradicting `(37)`.

The bounded range is decided without a node cap.  The companion searches all
`77` post-THM-3445 target-free primes through `443`.  It first enumerates the
`1,312` raw mass profiles, applies the proved central invoice `(12)`, and runs
the exact literal-mask DFS on the remaining `959` profiles.  A mandatory
`E_e` owner is odd-unit normalized to `8`.  The DFS tally is

```text
states=14,839, branches=13,966,
budget prunes=0, capacity prunes=12,054, reach prunes=1,220.   (40)
```

It finds precisely the five supported profile rows listed in Section 1:
four at `p=17` and one at `p=37`.  Every other bounded prime is negative.
The last finite prime is `443`; `445` and `447` are composite, so the next
prime is the analytic boundary `449`.

The positive controls also freeze the equality geometry:

- for the `p=17`, `A^4E_e^3` witness, `64` sheets have multiplicity one and
  the four central sheets `{8,25,42,59}` have multiplicity three; the only
  positive pairs are the three pairs among `E_e` coefficients `8,24,56`, each
  of weight four;
- for `p=17`, the `A^4E_oE_e^2` witness has `64` singleton sheets and four
  double sheets, while the other two profiles are partitions; and
- for the `p=37` witness, `140` sheets have multiplicity one and eight have
  multiplicity two.  The only positive pairs are `E_o=100` with
  `A=33,41,107,115`, each of weight two; `E_e=8` and `E_o=140` are isolated.

These profiles are positive controls, not assumptions in the exclusion.
Together with Sections 2--7 they complete the provisional proof of `(5)`.

## 8. Source, target, loss, and scope ledger

| field | exact content |
|---|---|
| source | labelled strict literal masks on all `4p` sheets |
| target | a typed weighted raw-ratio/coset graph |
| map | `ell -> N=2ell+1`, odd-unit normalization of the non-`E_e` endpoint, the reduced gap relation `(17)`, then the line/coset compiler `(20)--(23)` |
| preserved with sidecars | strict endpoints, owner type, coefficient label, exact pair weight, literal multiplicity factor, and central-fibre excess |
| destroyed by bare base projection | the `A` mod-four sheet, `B` mod-two coset, `E_o/E_e` parity, affine lift `m`, line class `J`, higher sheet multiplicities, and the union predicate |
| required sidecars | `(A,B,E_o,E_e)` type, `m`, congruence `(21)`, factor `2^f`, the `E_e/E_e` lattice, and the profile budget `Omega` |
| canonical positive hostile | `Q=68`, where three `E_e` owners coincide on the same four central sheets |
| oriented positive hostile | `Q=148`, where one `E_o` owner meets four `A` leaves with weight two |
| corrected near miss | projecting only to `Z/pZ` loses the dyadic coset and effective full-fibre parity |
| equality boundary | `Omega_max=12`; tail positive floor `15`; rooted zero clique number `6` |
| finite boundary | exact DFS through prime `443`; uniform proof from prime `449` |
| scope loss | no composite-even row, arbitrary centre/time, decrement, runner-current certificate, or physical LRC cover |

The pair observable is symmetric and has ties, so no intrinsic tournament is
present.  The faithful finite object is the labelled weighted fibre-grid
clutter together with its affine coset sidecars.

## 9. Reproduction and audited status boundary

Run from the repository root:

```bash
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_prime_quarter_half_twist_cap7_classification_thm3451.py
PYTHONHASHSEED=1 python3 -B -O 04-computation/lrc_prime_quarter_half_twist_cap7_classification_thm3451.py
```

Both transcripts must equal
`05-knowledge/results/lrc_prime_quarter_half_twist_cap7_classification_thm3451.out`
byte for byte.  The companion pins the three proved dependencies by
LF-normalized SHA-256, checks its own AST against `assert`, dynamic execution,
and unexpected imports, and freezes the semantic payload independently of the
rendered transcript.

The package passed a separate immutable audit of every analytic and finite
gate, including clean-room direct-mask and graph reconstructions.  It is a
proved dependency in exactly the scope stated above.  LRC(14) remains open.
