---
id: THM-3445
title: "Target-free prime-even literal half-twist cap-seven classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Target-free literal
  half-twist covers on Q=2p, p an odd prime, by at most seven distinct
  transverse owners have exact support p=7,19.  No composite-even,
  arbitrary-time, decrement, or LRC(14) consequence is claimed.
source: even-prime-rank7 typed interval and weighted-clique closure, 2026-08-15
audit: >
  independent theorem-universe, target-free period-descent, typed-size,
  overlap-budget, line-section, parity/AP, A-orientation, tangent-threshold,
  successive-minima, rational-descent, mod-240 graph, clique, finite-DFS,
  boundary-stitch, and witness audits; clean-room graph and direct-mask
  reconstruction; normal/optimized/stored and predecessor-atlas replays;
  dependency, hash, AST/security, ID, routing, scope, and documentation gates
  clean
depends_on:
  - THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists
related:
  - THM-3423-odd-interval-ratio-complement-and-dyadic-clique-law
  - THM-3434-seventeen-fibre-two-sided-mass-closure
script: 04-computation/lrc_prime_even_half_twist_cap7_classification_thm3445.py
output: 05-knowledge/results/lrc_prime_even_half_twist_cap7_classification_thm3445.out
script_sha256: ca061487acdf56f3d0507ad55d9b33c5e9724ddb9d7fdba815d46a1b7deb4382
output_sha256: bcdef4ca8e5ba5a786be31aa6757f704f4aae0e7ae0953e504701acc7a4290f2
semantic_sha256: 45e7e843b1767ef7254be729add71ea6ef501160503b46a03cbb43e3e53517c3
hash_basis: LF-normalized bytes
---

# THM-3445 -- target-free prime-even literal half-twist cap-seven classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  The analytic argument,
exact finite companion, and predecessor atlas passed independent derivation,
clean-room graph, direct-mask, hostile-boundary, and replay audits.

## 1. Exact scoped statement

Let `p` be an odd prime, put `Q=2p`, and use the literal strict masks

```text
B_r={ell in Z/(2p)Z: ||r(2ell+1)/(4p)||<1/14}.          (1)
```

Choose sign representatives

```text
1<=r<2p,             r!=p,                              (2)
```

so the owners are transverse and distinct.  Call the prime-even row
**target-free** when

```text
p not in {5,11,13,23,29};                              (3)
```

the excluded values are precisely the prime-even rows already containing a
lower supported literal base in the finite-atlas universe.

The candidate classification is

```text
some at-most-seven family covers Z/(2p)Z
                    iff p is 7 or 19.                   (4)
```

The two direct witnesses are

```text
p=7,  Q=14: (1,3,4,5,9,11,13),
p=19, Q=38: (1,9,17,20,21,29,37).                      (5)
```

Both have type profile `A^6E`.  The `Q=14` row is a partition.  At `Q=38`,
`34` sheets have multiplicity one and `4` have multiplicity two.

## 2. Inheritance, owner types, and the overlap budget

THM-3435 is the proved dyadic source.  On `Q=2p` it gives three exact owner
types.  Write

```text
d_A=1, d_B=2, d_E=4,       alpha_A=4, alpha_B=2, alpha_E=1,
alpha_X d_X=4.                                          (6)
```

For an owner of type `X`, put `q=r/d_X modulo p`.  In the centered base
coordinate its support is

```text
A: x odd, |x|<4p/14, one oriented literal sheet;
B: x odd, |x|<2p/14, both deck sheets;
E: x arbitrary, |x|<p/14, both deck sheets.             (7)
```

Every `E` owner contains the reflection-fixed base point; `A` and `B` avoid
it.  An `A/even` intersection equals its base intersection, an even/even
intersection is twice its base intersection, and only `A/A` needs an
orientation sidecar.

If `p>7`, write `p=14k+s`, and let a seven-owner row have counts `(o,b,e)`
of `(A,B,E)`.  The exact sizes and total excess
`Omega=sum_r |B_r|-2p` are

| `s` | `|A|` | `|B|` | `|E|` | `Omega/2` |
|---:|---:|---:|---:|---:|
| `1` | `4k` | `4k` | `4k+2` | `e-1` |
| `3` | `4k` | `4k` | `4k+2` | `e-3` |
| `5` | `4k+2` | `4k` | `4k+2` | `2-b` |
| `9` | `4k+2` | `4k+4` | `4k+2` | `b-2` |
| `11` | `4k+4` | `4k+4` | `4k+2` | `3-e` |
| `13` | `4k+4` | `4k+4` | `4k+2` | `1-e` |

For a target-free cover, THM-3435 proves that there are exactly seven owners,
at least two are `A`, and at least one is `E`.  The common fixed fibre of the
`e` E-owners already spends `2(e-1)` units of excess.  Comparing this with
the table makes `s=3` impossible and gives the sharp surviving maxima

```text
s:                    1    5    9   11   13
maximum possible Omega: 8    4    4    4    0.          (8)
```

Thus every unresolved cover has `Omega<=8`.  At every covered sheet with
multiplicity `m`, the excess contribution is `m-1`.  Therefore, for every
pair of owners in a cover,

```text
w(r,t):=|B_r intersect B_t| <= Omega <=8.               (9)
```

This is why the proof needs exact pair multiplicity, not merely a
disjointness heuristic.

## 3. The typed cross-interval complement

For Sections 3--7 assume `p>=607`.  Take an ordered pair of types `X,Y`,
excluding `A/A` and `E/E`, and put

```text
lambda=q_Y/q_X modulo p.                                (10)
```

Among the `14` points `0,lambda,...,13lambda` on the `p`-circle, two
consecutive points have circular gap at most `p/14`.  Their difference,
after removing a common divisor, gives

```text
b lambda=epsilon a (mod p),
gcd(a,b)=1,  1<=b<=13,  1<=a<=floor((p-1)/14).          (11)
```

An intersection point lies on a line

```text
b y-epsilon a x=j p.                                   (12)
```

Put

```text
A=alpha_X a,       C=alpha_Y b,       D=A+C.            (13)
```

The parity and deck data are exact:

| ordered sector | direct parity mismatch | admissible `j` | x-AP step | literal factor |
|---|---|---|---:|---:|
| `AB,BA` | `a+b` odd | odd; otherwise even | `2b` | `1` |
| `BB` | `a+b` odd | odd; otherwise even | `2b` | `2` |
| `AE` | `b` even | odd; otherwise all | `b` or `2b` | `1` |
| `BE` | `b` even | odd; otherwise all | `b` or `2b` | `2` |
| `EA` | `a` even | odd; otherwise all | `b` or `2b` | `1` |
| `EB` | `a` even | odd; otherwise all | `b` or `2b` | `2` |

In a mixed row the first step is for the mismatch branch and the second for
the compatible branch.  Consequently

```text
the two cross-type masks are disjoint
 iff the direct parities mismatch and D<=14.            (14)
```

The strict equality case belongs to the disjoint side.  Indeed, an
admissible nonzero line has `|j|>=1`, while a line section is nonempty only
when `14|j|<D`.  This proves the forward implication of `(14)`.  The reverse
implication, including a quantitative lower bound, follows from the exact
section certificate below.

The complete reduced **effective-relation** obstruction pairs `(a,b)` have
ordered counts

```text
AB 4, AE 6, BA 4, BB 12, BE 11, EA 6, EB 11.           (15)
```

These are counts before affine coefficient-lift/sign multiplicity.  For an
`A` target the unique odd lift survives both signs; for a `B` target, the two
signs have opposite effective parity and exactly one survives; for an `E`
target, canonical sign normalization retains exactly one half-range lift.
Thus the literal cross contributions at normalized roots are

```text
A root: AB 4, AE 6;
B root: BA 8, BB 12, BE 11;
E root: EA 12, EB 11.                                  (15a)
```

With the `A/A` calculation below, the normalized literal zero degrees are
therefore

```text
degree(A),degree(B),degree(E)=(19,31,23).                (16)
```

Finally, the raw coefficient ratio in `(11)` is

```text
r_Y/r_X=epsilon alpha_X a/(alpha_Y b) (mod p).          (17)
```

In every mismatch row the numerator and denominator on the right share a
factor at least two.  Hence every reduced cross zero ratio has height at most
seven.

## 4. Exact line sections and the uniform positive-weight floor

For any line index `j`, the open x-section in `(12)` has length

```text
g_j=p G_j/(14a),
G_j=max(0,min(A,14j+C)-max(-A,14j-C)).                  (18)
```

An open interval of length `g` contains at least
`ceil(g/(qb))-1` points of any fixed arithmetic progression of step `qb`.
Thus each admissible line supplies the integer lower bound

```text
floor((p G_j-1)/(14a q b)),                             (19)
```

multiplied by the literal factor in the table.  Because `(11)` implies
`p>=14a+1`, this bound is monotone in `p` at fixed `(a,b)`.

For `1<=a<100`, evaluation at

```text
p_0=max(607,14a+1),       1<=b<=13,                    (20)
```

gives the complete phase-free certificate:

| sector | cells | least nonzero lower bound | rows at most `8` |
|---|---:|---:|---|
| `AB` | `834` | `14` | none |
| `AE` | `834` | `13` | none |
| `BA` | `834` | `14` | none |
| `BB` | `834` | `12` | none |
| `BE` | `834` | `8` | `(3,10),(5,6)` |
| `EA` | `834` | `14` | none |
| `EB` | `834` | `8` | `(6,5),(10,3)` |

These four exceptions all have `D=16` and raw ratios `3/5` or `5/3`; they
are the only boundary tangencies needing their exact AP phase.  The other
`D=16` rows already have a lower bound above eight.

For `a>=100`, one does not enumerate.  Since `C<=52`, the centered plateau
`|14j|<=A-C` contains enough complete admissible lines.  Evaluating only the
two starting parities `a=100,101` gives uniform plateau lower bounds

```text
AB 52, AE 55, BA 40, BB 48, BE 54, EA 28, EB 44.       (21)
```

The plateau radius only increases when `a` advances in the same parity
class, so `(21)` proves the whole `a>=100` tail.

The four mixed tangencies have exact one-sign counts

```text
N_35(p)=#{x:     2p/21<x<p/7,  3x=p (mod 10)},
N_53(p)=#{x odd: 4p/35<x<p/7,  5x=p (mod 6)},           (22)

w_BE(3/5)=4N_35(p),       w_BE(5/3)=4N_53(p),           (23)
```

and reversal gives the `EB` rows.  Strict integer endpoints are respectively

```text
floor(2p/21)+1, floor((p-1)/7),
floor(4p/35)+1, floor((p-1)/7).                         (24)
```

Both counts increase by one under `p -> p+210`.  The last odd arguments with
count below three are `601` and `609`; the latter residue is not coprime to
`210`.  A complete residue period gives `N_35,N_53>=3` for every prime
`p>=607`.  Therefore every cross-type positive intersection has

```text
w(r,t)>=12>8.                                           (25)
```

## 5. The A/A orientation sidecar

Normalize the first `A` coefficient to `1`, and let `s` be the resulting odd
second coefficient modulo `4p`.  Apply `(11)` to `s modulo p`, and retain

```text
kappa=(b s-epsilon a)/p modulo 4.                       (26)
```

For centered odd phases `x,y`, literal intersection requires

```text
b y-epsilon a x=j p,       j=kappa x (mod 4).           (27)
```

This is the bit destroyed by a bare ratio modulo `p`.

- If `a,b` have opposite parity, then `kappa` and `j` are odd.  The phase
  congruence selects a compatible progression of step `4b`.  The masks are
  disjoint exactly for `(a,b)=(1,2),(2,1)`, giving
  `+/-2,+/-1/2`.
- If `a,b` are odd and `kappa=0`, then `j=0 (mod 4)`, the step is `2b`, and
  the direct point `(x,y)=(b,epsilon a)` prevents disjointness.
- If `a,b` are odd and `kappa=2`, then `j=2 (mod 4)`, the step is `2b`, and
  disjointness holds exactly for

  ```text
  (1,1),(1,3),(3,1),(1,5),(5,1).                       (28)
  ```

For each odd pair in `(28)`, the two coefficient signs have complementary
`kappa` values, so exactly one sign has `kappa=2`.  Together with the four
opposite-parity neighbors, this gives the `A/A` zero degree nine used in
`(16)`.

Formula `(18)` still applies with `A=4a,C=4b`.  The finite certificate has

| branch | cells | least regular bound | rows below `12` |
|---|---:|---:|---|
| opposite parity | `544` | `16` | none |
| odd/odd, `kappa=0` | `290` | `13` | none |
| odd/odd, `kappa=2` | `290` | `10` | `(3,5),(5,3)` |

For `a>=100`, plateau-only bounds are `48,52,48`.  Thus every nonzero `A/A`
intersection already exceeds eight.  The two middle rows admit the sharper
exact formula

```text
N_AA(p)=#{x odd: 4p/21<x<2p/7, 3x=2p (mod 5)},
w_AA=2N_AA(p).                                          (29)
```

Here `N_AA(p+210)=N_AA(p)+2`.  Crucially, the last odd `p` with
`N_AA(p)<6` is `605`, not the earlier weight-at-most-eight boundary `521`.
A full period `607<=p<=815` has minimum six.  Hence `(29)` is at least `12`
for every odd `p>=607`.

## 6. The E/E successive-minima floor

Put `h=floor((p-1)/14)` and, for the effective ratio `lambda`, consider

```text
L_lambda={(x,y) in Z^2:y=lambda x (mod p)},    det L_lambda=p.   (30)
```

Let `sigma_1<=sigma_2` be its successive minima for the unit `l_infinity`
square.  Minkowski's second theorem gives

```text
sigma_1 sigma_2<=p.                                     (31)
```

If `2sigma_1<=h`, the five lattice points
`0,+/-v,+/-2v` lie in the strict E-box.  Otherwise

```text
sigma_2<=p/sigma_1<2p/h<=h,                             (32)
```

where `h^2>=2p` holds from `p=607` onward (`h=43` at the boundary).
Then `0,+/-v_1,+/-v_2` are five distinct base intersection points.  Every
base point has two literal deck sheets, so

```text
w_EE>=10>8.                                             (33)
```

## 7. The finite typed zero graph

Let

```text
D_7={epsilon a/b:gcd(a,b)=1, a,b>0, a+b<=7}.           (34)
```

Sections 3 and 5 prove that every zero neighbor of a normalized owner has a
raw ratio in `D_7`, with endpoint type and `A` lift retained.  Conversely,
the exact parity and orientation tests decide every labelled candidate in
`D_7`.

This graph is genuinely finite and periodic.  If two candidate fractions and
their pair ratio agree only modulo `p`, clearing the three denominators gives
two integer terms of size at most `6^3`; their difference has absolute value
at most `432<p`.  Hence the modular relation is an equality over the
rationals.  A candidate coefficient solving

```text
d c=epsilon a r+t p                                    (35)
```

has its type and `A` orientation determined by `p modulo 4d`.  Since
`d<=6`, the entire labelled graph depends only on

```text
p modulo lcm(4,8,12,16,20,24)=240.                     (36)
```

The companion enumerates all `64` unit classes modulo `240`, using the
formula graph and checking a second prime in every class.  In every class it
finds

```text
root degrees (A,B,E)=(19,31,23),
maximum clique size=6,
rooted maximum presentations=(2,1,1).                  (37)
```

There are four normalized presentations of the equality orbit.  Every
maximum clique is exactly three complementary pairs `r,2p-r` and has type
profile `A^4BE`.  Direct literal-mask comparisons at `p=607` and `1009`
agree on every candidate edge and on every normalized zero neighbor; their
root-wise minimum positive weights are `(12,12,12)` and `(18,20,20)`.

## 8. Tail contradiction and finite stitching

Suppose `p>=607` supported a target-free cover.  By `(8)--(9)`, every pair
would have weight at most eight.  Sections 4--6 show that every positive pair
has weight strictly above eight.  Thus the seven distinct owners would be a
seven-clique in the exact graph of Section 7, contradicting `(37)`.

The remaining bounded range is stitched without extrapolation:

- the pinned FINITE-EXACT atlas
  `lrc_prime_even_half_twist_cap7_finite_atlas_20260815.py` searches every
  target-free odd prime `p<=599`, with no node cap, and finds only `p=7,19`;
- the same exact profile DFS decides `p=601` negatively in `15` profiles,
  `65` states, and `50` branches; and
- there is no prime strictly between `601` and `607`.

The direct controls `(5)` prove the positive direction.  This completes the
frozen candidate proof of `(4)`, subject to the required independent audit.

## 9. Source, target, loss, and boundary ledger

| field | exact content |
|---|---|
| source | labelled strict literal masks on all `2p` sheets |
| target | the typed weighted coefficient-ratio graph |
| map | odd-unit normalization, effective ratio modulo `p`, then the reduced gap relation `(11)` |
| preserved | endpoint types, exact pair weight, strict endpoints, and coefficient labels |
| lost by the bare base ratio | the `A` orientation/lift, literal sheet locations, the fixed E fibre, higher multiplicities, and the union predicate |
| required sidecars | ordered `A/B/E` types, `kappa mod 4`, deck factor, complementary partner, and fixed-fibre excess |
| equality boundary | strict `D=14` mismatch is disjoint; `D=16` contains the four mixed height-eight tangencies |
| finite boundary | exact predecessor atlas through `p=599`, then direct `p=601`; the uniform proof begins at `607` |
| scope loss | no composite-even row, four-sheet `Q=4R` row, arbitrary centre/time, decrement, or physical LRC cover |

The graph has symmetric zero/positive ties and no intrinsic orientation, so a
tournament would be cosmetic.

## 10. Reproduction and status boundary

Run from the repository root:

```bash
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_prime_even_half_twist_cap7_classification_thm3445.py
PYTHONHASHSEED=1 python3 -B -O 04-computation/lrc_prime_even_half_twist_cap7_classification_thm3445.py
```

The companion pins both predecessor scripts and outputs by LF-normalized
SHA-256, replays the `p=601` bridge and positive controls, and checks that it
contains no Python `assert`.  Normal and optimized transcripts must equal the
stored output byte for byte.  The exact dependency and semantic hashes are in
the frontmatter and transcript.

An independent audit checked the analytic quantifiers, AP existence,
orientation lifts, periodic graph proof, finite universe, and both execution
modes.  Thus the target-free prime-even classification is proved in its stated
scope.  In particular, LRC(14) remains open. **QED.**
