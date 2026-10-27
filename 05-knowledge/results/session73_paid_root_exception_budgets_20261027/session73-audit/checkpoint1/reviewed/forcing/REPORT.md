# Native ROOT forcing: the phase pays, the frontier still owes

2026-10-27, Session73. **PROVED scoped elementary deductions from inherited
ROOT weights and occupation identities; FINITE-EXACT controls.** This packet
does not prove all-source ROOT positivity, a uniform spatial ROOT-current
bound, or a uniform finite-policy zero-extension bound. It uses no random
orbit model and replaces no source measure.

## Findings before the research notes

1. For the inherited signed ROOT profile, forgetting its actual shortcut-clock
   phase makes the passive forcing bill infinite under the original native
   source probability `lambda`, already at the first odd transition.
2. Retaining the native phase repairs that bill exactly. The all-depth,
   original-source-weighted paid forcing is
   `E = sum_(m>1 odd) W(m)/(6m+1) <= 13/57`.
   The same upper bound applies to every finite collection of authenticated
   first-hit paths with their original `lambda` weights and native guards.
3. This is not a bound for zero-extending an arbitrary finite certified profile:
   every sharp finite first-hit horizon has infinite frontier forcing in this
   metric. Infinite inverse sibling families impose a genuine trace condition.
4. The native probability has finite expected **odd-step** lifetime (inherited)
   but infinite expected **literal** first-hit clock (new deduction). Modulus
   refinement or renaming the clock cannot erase that difference.

The scalar payment does not yet control the inherited spatial current
`R_f(N)`; its source aggregation, covariance and inverse-mass observation costs
remain separate. A separately audited additive supplement may sharpen this scalar constant;
it does not change these foundational statements.

## Inheritance and source pins

All eight source snapshots and their exact Git blobs are in `SOURCE_PINS.json`,
read at commit `92e7f12f3b95840cff92c5e6cf40b2ae8952cae4`.

- `05-knowledge/results/weighted_renewal_moving_modulus_20261026/README.md`,
  **R3**: the positive-branch passive identity with `b=3n/(3n+1)` and
  `d=1-b^2`; **R1/R4**: supplied-weight Hardy criterion and adaptive aliasing.
  These are inherited, not new claims here.
- `05-knowledge/results/collatz_adaptive_mixture_flow_20261005.md`, **P2/P3/P5**:
  the exact beta-mixture `W`, its incoming fibres, native probability `lambda`,
  and **already proved** occupation identity and mean odd-step bound. In
  particular this packet does not claim to discover its Green representation.
- `05-knowledge/results/marked_contacts_root_current_20261017e.md`, section 2:
  the genuine stopped ROOT shortcut-clock sign and the spatial folding kernel.
  `root_current_ground_truth_contact_rotations_20261017.md` identifies the
  uniform spatial energy estimate as a substantial still-open obligation.
- `collatz_source_refinement_floor_20261005.md`, definitions and SF4: positive
  support is an additional obligation; a weighted Hilbert quotient may erase it.
- `collatz_floor_transport_deadlines_20261005.md`, F3: finite approximate values
  of `W` do not establish a positive floor.
- `native_return_necklaces_20261026/README.md`, N6: branch guard and clocks.
- `collatz_poisson_source_dual_20261005.md`, P1/P2: independent signed dual
  inequalities, not mere norm existence, yield positive source floors.

Closest mechanism: R3's passive identity plus P5's **inherited** occupation law.
Canonical hostile: a deleted source atom, and the unsigned gauge's odd-valuation
inverse siblings below. Corrected near miss: an eventual-profile estimate is
not automatically a finite controller's total zero-extension estimate.
Least-used sidecar: the parity of the actual shortcut clock, together with the
complete named inverse fibre sharing one target port.

Concept board: **source mass / native clock phase / incoming fibre / passive
energy / certified frontier / spatial signed covariance**. The hostile changes
all six: it selects the original `lambda`, exposes a phase cost, uses infinitely
many siblings sharing a target, separates two energy domains, distinguishes
finite guards from total profiles, and leaves the spatial cross terms intact.

## 1. Types and conventions

Let `U(n)=oddpart(3n+1)` on positive odd integers. Write
`a(n)=v2(3n+1)`. ROOT is the **first** visit to 1; no outgoing root loop is
executed. A strict root word `(a_1,...,a_tau)` has

```
L=tau-1,  K=sum floor((a_i-1)/2),
W(n)=2 K! (L+1)!/(L+K+2)!    (n>1 rooted),
W(1)=1,  W(n)=0             (unrooted).
```

No finite word or counters are assigned to an unrooted source. On odd `n>1`,
`lambda(n)=W(n)` if `3|n`, and is zero otherwise. Inherited P2/P5 gives

```
sum W <= 16/3,   sum_(n>1) W <= 13/3,   sum lambda=1.
```

`lambda` is the native ROOT-supported source probability, not the unrelated
positive-at-every-source prior `nu_r`. For a rooted odd source define

```
f(n)=(-1)^(a_1+...+a_tau),   f(1)=1;
f(n)=0 on unrooted sources.
```

This is the eventual ROOT-supported extension of the profile in the inherited
current. A finite authenticated controller snapshot is a different object.
One odd transition uses `a(n)` **shortcut** steps but `1+a(n)` literal standard
steps (one `3n+1`, followed by individual halvings). Therefore

```
sigma(n)=(-1)^a(n),             f(U(n))=sigma(n) f(n).       (1)
```

At `n=21`, the word is `(6)`: shortcut sign is `+1`, while the seven-tick
literal-clock sign is `-1`. Both clocks must be retained; neither is replaced.

## 2. The exact native phase-aware invoice

Put

```
b=3n/(3n+1),  d=1-b^2=(6n+1)/(3n+1)^2,  c=sigma(n)b.
```

For arbitrary real `e,F`, set `e'=ce+F`, `r=de-cF`. Since `c^2+d=1`,

```
(e')^2+r^2/d = e^2+F^2/d.                                (2)
```

This is R3's elementary identity with the actual signed gain; it remains
positive because `n>0`. All statements also hold for complex signals with
absolute squares. On the genuine native profile, `e=f(n), e'=f(U(n))`, (1)
gives `F=(1-b)f(U(n))`. Consequently

```
F^2/d = r^2/d = 1/(6n+1)       on rooted n>1,              (3)
F=r=0                         on unrooted n.
```

No source was replaced and no phase chosen retrospectively: `sigma(n)` is
determined by the actual branch valuation. The sign convention is essential.
Using `b` instead of `c` gives instead

```
[f(U(n))-b f(n)]^2/d = 6n+1          if a(n) is odd,
                         1/(6n+1)  if a(n) is even.       (4)
```

### A genuine native first-step divergence

For every integer `k>=0`, let

```
n_k=(10*64^k-1)/3 = S^(3k)(3),   S(x)=4x+1.
```

Then `3|n_k`, the strict word is `(6k+1,4)`, `U(n_k)=5`,
`f(n_k)=-1`, `f(5)=1`, and

```
lambda(n_k)=4/[(3k+2)(3k+3)(3k+4)].                      (5)
```

The unsigned-gauge invoice contributed by this source is exactly
`lambda(n_k)(20*64^k-1)`, whose terms tend to infinity. Its sum is infinite.
The phase-aware invoice is instead `lambda(n_k)/(20*64^k-1)`.
This refutes the **unsigned-gauge** budget for this specific native current,
not all paid carriers. It is not merely a beta-weight-score surrogate.

## 3. The original native weights pay every authenticated path

Let `V={3,5,...}` and let `K` push source masses along `U`, killing an edge
when its target is ROOT. The inherited P5 identity is

```
W|V = sum_(t>=0) K^t lambda.                             (6)
```

Its proof iterates `W-KW=lambda`; the remainder goes to zero in `ell1`
by finite rooted support paths and summable `W`. Thus for every nonnegative
vertex cost `q`, Tonelli gives the **same-source** formula

```
sum_s lambda(s) sum_(0<=t<tau(s)) q(U^t(s))
     = sum_(m in V) W(m)q(m).                            (7)
```

The source label `s` and mass `lambda(s)` remain unchanged throughout a path.
Regrouping visits by target is not a resampling or a replacement by `W(U^t s)`.
Equation (7), not new here, is the exact payment interface.

Apply it to (3):

```
E := sum_s lambda(s) sum_(t<tau(s)) F_(s,t)^2/d_(s,t)
   = sum_(m in V) W(m)/(6m+1)
   <= (1/19) sum_(m in V) W(m) <= 13/57.                  (8)
```

This proves a uniform all-depth forcing budget for the **eventual native
profile with the correct branch phase**. It also bounds the sum over any
finite collection of distinct authenticated ROOT sources and any retained
prefixes of their first-hit paths: all summands are nonnegative. Distinct
named sources count once with their given mass; copying a source into several
records does not authorize duplicating its probability mass. A source split
must split its mass, and every alias/port stays named.

A finite guarded path certificate therefore has a paid forcing account.
But it does not establish the ROOT status of a missing source, and it must
not silently impose a total forcing equation on uncertified frontier edges.
Also, diagonal temporal payment is not the spatial covariance estimate for
`R_f(N)` in the inherited folding problem.

### The exact incoming contraction is on W, not on lambda alone

On all positive odd vertices, define `Cf(n)=f(U(n))` for `n>1`, and `Cf(1)=0`.
The root incoming row has mass `sum_(j>=1)W(S^j1)=1`; all other rows are P3.
For every `h in ell2(W)`,

```
||h||_W^2-||Ch||_W^2 = sum lambda(n)|h(n)|^2.             (9)
```

Thus this native incoming carrier retains `W` as vertex mass and `lambda`
as the boundary defect. It cannot be defined on `ell2(lambda)` alone:
`delta_ROOT` has zero lambda-norm but `C delta_ROOT` has positive lambda-norm
(at source 21 already). The `ell2(W)` quotient likewise ignores all `W=0`
sources. Contractivity cannot establish their missing positive mass.

## 4. Finite profiles have an unpaid frontier

Let `A` be a subset of the authenticated ROOT component, contain ROOT,
and be forward closed until ROOT. Put `f_A=f 1_A`. Its phase-aware total-profile forcing is

```
F_A(n)=f(U(n))[1_A(U(n))-b(n)1_A(n)].                    (10)
```

It pays `(6n+1)^(-1)` on `A`, zero when neither endpoint is in `A`, and
`1/d_n` when `n notin A` but `U(n) in A`. Those are three distinct domains.

**Every finite A has infinite first-step lambda-weighted frontier bill.**
Indeed it omits all but finitely many sources

```
z_k=(64^(k+1)-1)/3,
lambda(z_k)=2/[(3k+3)(3k+4)],   U(z_k)=1.                (11)
```

They have `d_z <= 2/(3z+1)`, so the individual frontier invoices
`lambda(z_k)/d_z` grow exponentially rather than tending to zero.

**Every finite sharp first-hit horizon also fails.** Put
`A_T={n:tau(n)<=T}`, with ROOT included. Choose a rooted unit target `m`
with `tau(m)=T`. Such targets exist at every T: starting at 1, choose a
larger unit inverse child among legal exponents at most seven; a full proof
is below. Let `b0(m)` be its minimal inverse base, and take the infinitely
many siblings `S^j(b0(m))` divisible by three. Their first-hit length is T+1,
so they are frontier sources. Their native weights decay only polynomially
in j (the factorial ratio with fixed target counters), whereas `1/d` grows
like `4^j`. The first-stage forcing sum is infinite.

For the existence assertion, a unit target has one parity of legal exponents.
For `m=1` use child 5. For `m>1`, among exponents `{2,4,6}` if `m=1 mod3`
and `{3,5,7}` if `m=2 mod3`, exactly one of the three children is divisible
by three; the remaining two are units and exceed m. Each prepends one actual
first-hit edge, producing the next depth. This is the inherited SF5 inverse
palette, used only for this scoped existence statement.

Therefore a pointwise or weighted-L2 limit of certified horizon profiles
may have the finite energy (8), although every sharp horizon has infinite
energy in the total-profile forcing metric. This is not a contradiction and
is not a claim that no other approximation scheme exists.

### Incoming trace condition

For the root family (11), write `alpha_k=lambda(z_k)/d_(z_k)`. It tends to
infinity. If a signal g has finite first-stage phase-aware energy
`E_g=sum lambda(n)|g(U n)-sigma(n)b(n)g(n)|^2/d_n`, then

```
|g(z_k)-g(1)| <= [sqrt(E_g/alpha_k)+(1-b(z_k))|g(1)|]/b(z_k) -> 0. (12)
```

This follows by bounding one summand and rearranging; `sigma(z_k)=1`.
The error is `O(k*8^(-k))` for fixed E_g and root value. Thus a nonzero
ROOT value requires an actual trace on arithmetically distant inverse
siblings. Finite-support nonzero-root profiles are excluded from this
energy domain. For any rooted unit target, its lambda-positive sibling
subsequence similarly has trace `sigma*g(target)`. No condition is asserted
at an unresolved zero-weight target. A finite symbolic description with
infinite support is not excluded; the complete named target fibre is the
natural port to retain.

## 5. Odd-step payment does not pay the literal clock

On the root family (11), the odd-step length is one, shortcut clock is
`6k+6`, and literal standard clock is `6k+7`. Direct arithmetic gives

```
lambda(z_k)(6k+7) >= 1/(k+1).
```

Hence the native lambda probability has **infinite expected literal
first-hit clock**, despite inherited `E_lambda[tau_odd]<=13/3` and the finite
forcing budget (8). Its expected shortcut clock is infinite as well, since
subtracting the integrable odd-step count from the literal clock cannot
make both finite. The branch-precision consumption `sum a_i` is the shortcut
clock on these odd trajectories, so it too has infinite mean.

Every individual certified source still has a finite integer clock. An
infinite mean is not nontermination, and increasing an observation modulus
does not change these operation counts. Full source-weighted energy and full
source-weighted execution cost are different accounts.

## 6. Connection contracts and boundaries

| Source -> target | Map and preserved data | Loss / required sidecar |
|---|---|---|
| R3 scalar contraction -> native signed profile | `c=(-1)^a b`, same source and `d` | Branch parity must be actual shortcut parity, not literal-clock parity |
| Native lambda paths -> vertex occupation | Equation (7), retaining every original source and visit | Aggregation loses provenance unless source/clock are retained as sidecars |
| Finite guarded certificates -> paid forcing | Restrict (8) to their literal first-hit paths | Missing frontier equations remain obligations, not zero values |
| Eventual profile -> sharp finite horizon | Pointwise restriction | Does not preserve finite passive energy; complete inverse target ports are needed |
| Native source energy -> spatial current | Only a weighted observable map is available | Inverse source masses and signed spatial covariance are unpaid |
| Clock -> address precision | Retain `a` and exact word/cylinder | Precision and execution are distinct; means can diverge under native lambda |

**OPEN:** independent positivity for every native source; a useful finite
controller approximation in the strong energy domain with its frontier paid;
uniform spatial current bounds; transfer to a different supplied all-positive
source measure; the uniform efficiency of any full tournament/arithmetic
compiler. No LRC(14), Collatz, or tournament-classification closure follows.

## 7. Reproduction and audit scope

```
python3 check_core.py > core.normal.json
python3 -O check_core.py > core.optimized.json
cmp core.normal.json core.optimized.json
```

Checks use exact integers/Fractions and explicit exceptions, not removable
assertions. The universe is all 501 positive odd sources through 1001 with
literal first-hit replay; all their 11,088 odd transitions; exact analytic
incoming-row tails for 51 targets at 12 cutoffs; all 3,125 value assignments
in `{-2,-1,0,1,2}^5` on a stated target set; 20 terms in each infinite-family
hostile; and 13 first-hit horizons with eight explicit frontier siblings each.
The all-order assertions rely on the written proofs, not extrapolation.
Finite lambda-source occupation is independently regrouped by named visited
vertex and compared against each exact native W. All four clocks/roles
(source index, odd-step index, shortcut ticks, literal ticks) remain typed.
The separate agent audit is delivered by the parent integration packet; this
file itself claims only producer-authored exact controls.
