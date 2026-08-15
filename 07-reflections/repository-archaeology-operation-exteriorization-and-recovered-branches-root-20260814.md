# Repository archaeology: recovered proofs, operation exteriorization, and next probes

**root/repository-archaeology-2026-08-14.  Reflection, not canon.**  Current
truth lives in the cited theorem and correction files.  `PROVED CANDIDATE`
below remains outside the proof graph until independently audited and promoted.

## Outcome

This was not merely a citation-count scan.  It inspected live remote refs,
reflogs, unreachable objects, stale reserved namespaces, orphan computations,
first untested predictions, and operation columns around current tournament,
GMC, JC, and LRC work.

The session produced:

1. [THM-3079](../01-canon/theorems/THM-3079-laguerre-pf-row-transform-and-strict-integer-mesh-terminal-minus-one.md),
   a genuinely unmerged Newton--PF theorem, recovered from one remote branch,
   independently audited, and promoted to **PROVED**;
2. [THM-3369](../01-canon/theorems/THM-3369-skew-deletion-response-and-ordered-join-orientation-current.md),
   a **PROVED** skew exterior response which restores one ordered-join bit lost
   by the symmetric deletion Gram;
3. [THM-3372](../01-canon/theorems/THM-3372-multiaffine-deletion-transform-variance-and-skew-join-current.md),
   a **PROVED** multiaffine deletion transform, sharp transitivity variance,
   and scalar skew join current extending dormant THM-026;
4. a **REFUTATION** of THM-300's quadratic-signature conjecture at its first
   untested size `n=9`, with the extra negative direction localized to the new
   layer Schur complement;
5. a repair and all-cycle-length strengthening of **PROVED** THM-309, replacing
   false ordered 2-transitivity by sharp unordered-pair transitivity;
6. correction of THM-101's missing `rk(d4)` term and closure of stale reserved
   THM-1279 as a superseded scalar target; and
7. three coherent branch-only JC proof candidates, a forgotten exact LRC
   phase toolkit, and several small typed transfers ranked below.

None of this proves LRC(14), JC(2), DC(2), FC(3), or a new GMC statement beyond
the declared theorem scopes.

## Search coverage and negative archaeology

At the branch auditor's frozen initial snapshot `origin/main=d9db6b53146a`, the
repository had `74` actual local/remote refs: `38` were ahead of main and `36`
were contained.  There were `396` live-ref-reachable commits absent main and
`1,078` commits in the union of live refs and reflogs absent main.

The genuinely loose-object layer was mostly noise rather than treasure:

- `334` commits were reachable by neither refs nor reflogs;
- `329` were stash scaffolding (`WIP`, `index`, `On`, or untracked commits);
- the five substantive-looking clones were already integrated by patch or
  exact file content; and
- no new theorem came from the true loose-object set.

That negative result changes future search strategy: live divergent branches,
reserved stubs, stale exact outputs, and the first untested prediction are much
higher-yield than another indiscriminate `fsck` sweep.

Major nonrecoveries were also classified.  The reflected THM-2990 branch is
superseded by THM-3349/3350/3352/3355; old THM-2595 is current THM-2604; the K4
Hafnian scratch became THM-3058; modular Farey became THM-3035; the fixed
finite-denominator LRC bank is refuted by THM-2072; and several product-ring or
solvable-radical routes lose field or globalization hypotheses.  These should
not be resurrected as theorems.

## Inheritance pass and live concept board

The anchor remained **LRC(14), OPEN**.  The niche was branch-proof recovery;
the wildcard was relation exteriorization of commutative responses.  The
inheritance pass used:

- closest proved mechanisms: THM-3056/3065 for Gamma/Hankel signs,
  THM-3324 for deletion responses, THM-026/002 for deletion OCF, THM-3356/3359
  for folded arithmetic support, and THM-3068/3072/3076/3077 for the C3/JC
  valuation tower;
- canonical hostile: `K1 triangleright C3` versus its reverse;
- corrected near misses: THM-3065's zero-prefix H3/H4 failures, THM-300's
  extrapolated inertia, THM-309's false symmetry premise, and THM-1279's
  gcd-forgotten decay; and
- least-used sidecars: skew relation signs, full deletion Taylor levels,
  folded gcd pairs, endpoint owner/phase, and the initial-coefficient module
  of a Laurent key tower.

The board stabilized at six concepts:

1. **operation response**, rather than one more scalar invariant;
2. **controlled forgetting**, with an explicit kernel and smallest sidecar;
3. **branch-only proof chains**, separated from canon until hostile audit;
4. **new-layer Schur complements**, rather than fitted inertia sequences;
5. **symmetry matched to the actual incidence type**; and
6. **polynomiality debt after Laurent/valuation success**.

Each successful pull changed another lane: the skew-response law suggested the
deletion transform; the deletion transform exposed the transitivity variance;
the Paley repair demonstrated that weaker symmetry can prove a stronger
correctly typed conclusion; and the THM-300 failure promoted Schur-layer tests
for every nested matrix conjecture.

## The reusable abstraction: response--relation exteriorization

Suppose a response `F_T` is commutative under ordered join,
`F_(X triangleright Y)=F_X F_Y`, while marked responses `r_(T,v)` transform
factorwise.  Contract the marked row against the intrinsic lost relation
`K=A-A^T`:

```text
Xi_T(z,w)=r_T(z) K_T r_T(w)^T.                         (A)
```

The diagonal blocks of `K_(X triangleright Y)` reproduce the factor currents;
the cross blocks `+J,-J` give one exterior product.  Composition closes once

```text
q_T=sum_v r_(T,v)                                      (B)
```

is derived or carried.  Three genuinely different responses now test the
schema:

- THM-3369: the split spectral response `s=(P,zN)`, with `q=ns-zs'`;
- THM-3372: the Hamiltonian deletion transform `D`, with `q=D'`; and
- THM-3166's complete path-colour polynomial `Q`, where exact order-six masks
  `16` and `83` have the same `Q` but different `q`, proving that `(B)` is a
  genuine extra sidecar there.

This is now recorded in
[META-PATTERNS.md](../00-navigation/META-PATTERNS.md).  Counterindications are
essential: a nonconstant/cyclic cross block, unlawful deletion response,
target-bearing ties, or a target inside the contraction kernel stops the
transfer.  `Xi` is not automatically injective or chronological.  THM-3372's
order-five masks `8,10` have the same `(D,xi)` and distinct marked decks; every
self-converse object forces the skew current to vanish.

Two noncanonical probes sharpen the next boundary.  First, keep each Boolean
deletion layer as a multiset rather than only its sum:

```text
M_T(t)=sum_X t^|X| [H(T-X)],       [a][b]=[ab].        (C)
```

This is an exact commutative semiring compiler,
`M_(X triangleright Y)=M_X M_Y`; `D` is its additive-moment quotient.  A
**FINITE-EXACT** census through all `456` order-seven classes finds exactly the
same collision fibres for `M` and `(Gamma,D)`, while `(M,xi)` is injective.
Order eight is the first honest test of whether that equivalence is structural
or a small-order coincidence.  Second, replacing the OCF weight `2^|S|` by a
formal fugacity `x^|S|` preserves every deletion, product and skew-current law
coefficientwise.  Its first possible separation of packing count from covered
support is order nine, where full support can have cycle types `{9}` and
`{3,3,3}`.  These are research probes, not claims of all-order completeness.

## Ranked recoveries and extension probes

| Rank | Object and status | Source -> target map | Preserved / destroyed information | Needed sidecar and cheapest decisive test |
|---|---|---|---|---|
| 1 | THM-3074/3080/3081 C3 toric key chain — **PROVED CANDIDATES + VERIFIED-EXACT**, branch-only | coordinate-line C3 escape -> primitive first key -> finite depth partition/gcd descent -> terminal Möbius decoder | preserves completed-field valuation, depth budget, residue field, and root torsors; loses global polynomial realizability | retain the initial-coefficient ring/module inside `C[x,y]`; test low budgets `h=E±1`, where exact packets polynomialize exactly one target and force the other reciprocal |
| 2 | Forgotten LRC exact contact/phase toolkit — **STALE FINITE-EXACT**, not a theorem | old contact graph/F2 rank/phase-labelled Fourier vector -> current `110` walls / `12` families at `z1=216` | phase distinguishes translated collision pairs that pair sums and Fourier magnitudes merge; old runner conventions and owner legality are lost | port conventions and retain phase plus owner; hostile bank: AP13, V*, stored q=9/q=25 translate collisions, then THM-2058's globally blind pair |
| 3 | THM-3166 path-colour exterior current — **PROVED-DERIVED + FINITE-EXACT probe**, not canon | multiplicative `Q` and marked deletions -> `psi=dKd^T` ordered-join current | keeps one orientation channel and richer path-cover grading; destroys marked ownership and full order | `q^Q=sum_v Q_(T-v)` is independent: masks `16,83` share `Q=t^6+8t^4+8t^2` but have `q` ending in `8t` versus `4t`; freeze this hostile before theorem filing |
| 4 | THM-3356/3359 folded matching colour — **RECOVERED PROVED ROUTING**, not novelty | antipodal root classes -> unordered gcd pair `kappa={gcd(M,s-r),gcd(M,A(r+s)+2h)}` | restores the three Boolean XOR perfect matchings lost by common harmonic density; destroys lift order, owner, phase, orientation, and ancestry cost | at `M=1105`, verify colours `{13,85}`, `{5,221}`, `{17,65}` on the six K4 edges; add owner/endpoint words only if physical transport is attempted |
| 5 | Terminal-prefix `-2` Newton--PF extension — **OPEN** | THM-3079's one reciprocal-Gamma baseline -> two baselines / compound PF transform | would preserve all-order checkerboard signs beyond terminal `-1`; generic Hadamard closure and shallow-strip signs are unsafe | search exact rational minors of orders `3..6` using THM-3079's 128 strict-prefix controls and THM-3065's H3/H4 hostiles before seeking a convolution proof |
| 6 | THM-300 layer-defect classification — **OPEN after REFUTATION** | nested quadratic coefficient matrices -> exact old/new Schur inertia word | preserves the true new-layer contribution; destroys the misleading global fitted sequence | extend Schur complements past `n=12`; first ask whether the `(5+,2-)` seven-tile `n=9` layer is sporadic, periodic, or a partition-threshold event |
| 7 | THM-1279 decorated multi-tooth replacement — **OPEN; scalar target SUPERSEDED** | located lcm bank -> joint `(gcd, endpoint-owner pair, chronology)` packet | retains correlations erased by `1/lcm=O(1/hj)` scalarization; no scale-free projective floor survives without them | use a primitive fixed-ratio near-coprime scaling family and measure the full distinct-owner joint invoice, not the sum of marginal lcm credits |
| 8 | Modular repeated-current word — **PROVED corollary if the chosen current is fixed; low priority** | repeated joins -> C-finite modular support and THM-3359 harmonic scar | exact word retains current sign; density/scar forget it, diagonal `z=w` kills it, bad moduli erase coefficients | for the THM-3372 hostile, test `xi_(T^r)(1,2)=42r60^(r-1)` mod 11; compare the signed word, not only density `10/11` and scar `60/121` |

## Branch-only JC signal in more detail

The branch `origin/codex/jc-resolvent-bridge-20260801` contains:

- THM-3074 at `cb3167b192`: a two-pole leading wedge forces a primitive
  binomial and a first depth key;
- THM-3080 at `18ca246002`: strict steps spend a finite budget with
  `B_next=B-e`, `g_next=gcd(g,e)`, so terminal depths partition the budget and
  the last lattice is primitive; and
- THM-3081 at `a11b9bb805`: the terminal residue ratio generates `C(u)`, hence
  is Möbius, and earlier stages form exact `mu_d` torsors.

All three companions replay byte-identically in normal and optimized modes,
and no current correction or supersession was found.  They remain candidates
because the completed Laurent-field argument has not yet been independently
audited into a polynomial Keller chart.

The cheapest new low-budget packet makes the obstruction concrete.  With
terminal budget `D=E`, one-stage key

```text
R=u s^h,
M=1+[E/(h-E)]u s^E,
x=R^(-1), y=MR,
```

one has `dx wedge dy=du wedge d(s^E)`, `T=xy-1`, and
`theta=T^h x^E=c^h u^(h-E)`.  Residue degree one forces `h=E±1` in this family.
For `h=E-1`, the `Q` target is polynomial while `P=c^h/theta` is reciprocal;
for `h=E+1`, `P=theta/c^h` is polynomial while `Q=T/(cP)` is reciprocal.  At
`E=3`, these are the explicit `h=2` and `h=4` dual cells.  Depth and gcd are
already correct; simultaneous polynomiality is the missing coordinate.

There is also an exact arithmetic sieve on the depth ledger.  If the positive
depths `(e_0,...,e_N)` sum to `D`, terminal primitivity is
`gcd(h,e_0,...,e_N)=1`.  Since the coordinate-line branch has `D=E (mod h)`,
only primes dividing `gcd(h,E)` can survive to the terminal gcd.  The number
of primitive compositions is therefore

```text
sum_(d | gcd(h,D)) mu(d) 2^(D/d-1).                  (D)
```

For `E=3`, all compositions pass when `3` does not divide `h`; otherwise the
counts at `D=3,6,9` are `3/4`, `30/32`, and `252/256`.  Thus generic partition
enumeration is the wrong next expense: retain the first depth not divisible
by three as the unique `mu_3`-discharge event and test whether the three
`A4/C2` charts synchronize that event inside their polynomial initial modules.

## Corrections found by treating failure as data

- [THM-300](../01-canon/theorems/THM-300-quadratic-signature.md) is
  **REFUTED**: exact inertia at `n=9` is `(20,8,0)`, and the new layer has two
  negative directions.  The observed `negative=n-1` for `9..12` stays OPEN.
- [THM-309](../01-canon/theorems/THM-309-five-cycle-design-paley.md) remains
  **PROVED and strengthens**: square-affine maps act sharply on unordered
  pairs, so every simple directed cycle length forms a pair-balanced
  multidesign.  The old long trace table counted closed walks, not cycles.
- [THM-101](../01-canon/theorems/THM-101-dt-cancel-filling.md) retains its
  finite `beta_2=0` conclusion, but surplus is
  `rk(d4)+beta_3-beta_2`, not `beta_3-beta_2`.
- [THM-1279](../01-canon/theorems/THM-1279-complementary-near-top-cluster-closure.md)
  is no longer a misleading live reservation: the ratio-only target is
  superseded, with decorated owner/gcd/chronology retained as the open repair.

The mechanisms and minimal witnesses are frozen in MISTAKE-378--380.

## Recommended next session

Use the portfolio, not a single tunnel:

1. **Anchor — LRC:** port the exact phase/contact toolkit only far enough to
   test AP13, THM-2058, and the twelve current families; stop if phase is merely
   a coordinate gauge or owner legality fails.
2. **Niche — JC:** independently audit THM-3074/3080/3081 and turn the
   `h=E±1` polynomial/reciprocal dichotomy into an initial-coefficient-module
   obstruction or a realizability criterion.
3. **Wildcard — response exteriorization:** freeze the THM-3166 `psi` law with
   the masks `16/83` sidecar hostile; in parallel, test the deletion-layer
   compiler `(C)` at order eight and the packing fugacity at order nine.  Then
   try a genuinely nontransitive substitution whose cross block is not rank
   one.  That is where the abstraction either grows or honestly stops.

The main methodological conclusion is simple: forgotten work becomes valuable
when it is pulled through a current operation, boundary, or lost-coordinate
question.  Branch existence alone is not progress; an exact map, an honest
kernel, and a cheap hostile are.
