# AMM 12592 -- hostile audit of the golden floor chain (boxeph, 2026-08-03)

Chain audited: THM-3009 -> THM-3017 -> THM-3027 -> THM-3024.
Referee: `04-computation/amm12592_golden_floor_hostile_audit_boxeph.py`,
output `05-knowledge/results/amm12592_golden_floor_hostile_audit_boxeph.out`.
Exact arithmetic only (int / Fraction / sympy 1.9); every float below is
display-only.

## 0. Verdict summary

| Theorem | Verdict | Scope the result deserves |
|---|---|---|
| THM-3009 | **SOUND** (floor direction) | `C* > 1.5970` for balanced-block schemes is PROVED-exact; the asymptotic floor `C*_block >= log_5(5 phi^2)` is rigorous **modulo one routine, unwritten Stirling-transfer lemma** (sec 10.3 of the file, failure direction only). Threshold *exactness from above* still rests on an uncertified numeric scan. |
| THM-3017 | **SOUND within scope** | Threshold of a necessary criterion for the `H=1` checkpoint-closure program. Its 20-digit numerics are superseded by the symbolic collapse (THM-3027 + this audit). |
| THM-3027 | **SOUND** | Tangency collapse re-derived fully symbolically here; two upgrades delivered (inner concavity, gamma-monotonicity) make the floor direction scan-free; `sigma*` interiority upgraded to an exact certificate. Laplace/Stirling rate passage remains ASSUMED (self-declared). |
| THM-3024 | **GAP FOUND -- HIGH severity** | The promotion `C*_general = log_5(5 phi^2)` is **not established**. Within the theorem's own transportation model the degree-resolved Hall cuts do not bind below golden once the shell window extends past the binding shell; the reported binding cuts are truncation edge artifacts. Demote; general-class floor is OPEN again. |

The honest bottom line for the problem: `C* in [1.5970..., 1.6]`-style claims
survive only in the balanced-block class; `C*_general = golden` reverts to
hypothesis (HYP-9061 territory). The **upper** bounds and the THM-3029
attainment work are untouched by this audit.

## 1. Classification of every load-bearing step

### THM-3009 (archimedean floor, balanced blocks)

| Step | Status |
|---|---|
| Deviation form (D), exact centring (hockey-stick) | PROVED-symbolic; identity re-verified exactly (part 3.0) |
| Reduction A (no 2-adic obstruction) | PROVED-symbolic (Lucas + F2); re-verified for m = 8..64 |
| Reduction B ((ARCH) necessity), `sum_i C(a,i)C(i,r) = C(a,r)2^(a-r)` | PROVED-symbolic; identity re-verified exactly |
| (ARCH) monotone in `a_k` (validity of the binary search) | PROVED (Pascal) **and** hostile full-ladder sweeps at m = 8,16,32 are monotone (part 3.1) |
| Finite-m certified bounds m <= 4096 | PROVED-exact. Independently recomputed from scratch for m = 8..1024 (all match committed values); committed refuted candidates at m = 2048, 4096 re-verified to fail and their successors to pass (part 3.4) |
| Continuum limit (T) of (ARCH) | Heuristic scaling as stated; the **failure direction** (all that the floor needs) is routine Stirling but the transfer lemma of sec 10.3 (floor perturbations `O(log m)`) is sketched, **not written**. Severity: LOW. |
| Closed form of the tangency system (sec 3.1, (CLO) lemmas) | PROVED-symbolic; re-verified in sympy (part 2) |
| Interiority (I1)-(I3) | PROVED-symbolic structure; the two strict constants (`ell*` margins 0.34 and 0.0617) verified at 60 dps against closed log-forms |
| Global minimality of `delta*` at `gamma*` (sec 10.2 scan) | VERIFIED-NUMERIC, no Lipschitz certificate (the file admits this). **After this audit it is no longer load-bearing for the floor** (see 2b below); it remains load-bearing only for "the criterion's threshold is exactly gamma*" (upper direction). Severity: MODERATE for exactness, NONE for the floor. |
| Sec 11 fold constant / convergence-rate analysis | Heuristic diagnostic, not load-bearing |

### THM-3017 (golden threshold for checkpoint closure)

Variational conditions (S)/(T)/(V): PROVED (calculus; derivative formulas
re-verified symbolically in part 2). The closed-form solution
`rho = sqrt5, sigma = phi` was only VERIFIED-NUMERIC (20+ digits) in this
file; it is now PROVED via the THM-3027 collapse. Its part-B "interior
minimiser" check samples 8 points -- weak, but superseded. Scope statement in
the file (necessary criterion for one sufficient program, not `C*`) is
correct and should be preserved.

### THM-3027 (capacity threshold = log(phi)/log(sqrt5))

| Step | Status |
|---|---|
| Rate reduction of (*) by Laplace/Stirling | **ASSUMED** (self-declared). Needs an error-term lemma. For the *barrier/floor* direction this is routine (term-wise Stirling upper bound on the sum, lower bound on the target, floors absorbed in `O(log R)/R`); for the *ampleness* direction (`gamma > gamma*` closes all large epochs) it additionally needs uniformity in `tau` near the boundary and is not written anywhere. |
| (S),(T),(V) are the threshold conditions | PROVED; derivative formulas re-verified in sympy (audit part 2) |
| Key identity (K), multiplier cancellation, `(V) <=> (1-tau)^2 = tau` | PROVED-symbolic; re-derived independently in positivity-safe coordinates, residuals identically 0 |
| `gamma* = log(phi)/log(sqrt5)` given the rate reduction | PROVED-symbolic (surd back-substitution re-verified: `u* = 1/sqrt5`, `(1-tau*)/tau* = phi`, `2+phi = phi sqrt5`) |
| `sigma*` interior | **UPGRADED to exact**: `sigma*` is monotone decreasing in `gamma`; with the certified rational bracket of gamma* (below) and exact Q(sqrt5) arithmetic, `0 < sigma* in (0.038635396253, 0.038635396266) < tau*`. |
| Inner max over sigma is the stationary point | **UPGRADED to PROVED**: `d2psi/dsigma2 = -[(1+g)/(1-rho)+1/rho](1+tau)/(g(1+sigma)^2) < 0` verified symbolically -- psi is strictly concave in sigma (perspective of the entropy), so no scan is needed for the inner problem. |
| Global worst tau = tau* | VERIFIED-NUMERIC (3000-point scan) only. **Not load-bearing for the floor** (see below); load-bearing only for threshold exactness from above. |
| b-alphabet universality of `tau* = phi^-2` | PROVED-symbolic (identity (K_b) re-verified with symbolic b) |

**New lemma delivered by this audit (floor direction made scan-free).**
`dpsi/dgamma = -(1+sigma) log(1-rho) > 0` (verified symbolically), so
`Psi(gamma, tau)` is strictly increasing in `gamma` at fixed `tau`. Combined
with concavity-in-sigma, interiority of `sigma*`, and the symbolic value
`Psi(gamma*, tau*) = H(tau*)` from the collapse, one gets **symbolically**:
for every `gamma < gamma*`, `Psi(gamma, tau*) < H(tau*)` -- the continuum
criterion fails at `tau* = phi^-2`. Hence the numeric global-tau scan is NOT
needed for the lower-bound direction; the only analytic debt of the floor is
the finite-`R`/finite-`m` Stirling transfer lemma.

### THM-3024 (degree-resolved Hall, general-class promotion)

| Step | Status |
|---|---|
| Reduction of the general class to the forward transportation model (per-shell demand `binom(m-1,d)`, supply = (ARCH) boxes of the extremal profile, forward-only routing) | **MODELLING PREMISE, not proved.** Nowhere derived from THM-2966/THM-3008 word-level structure. This is the load-bearing step and it is missing. |
| (G1) "the degree-blind cut is a valid Hall cut for any routing rule in d, so opus's sign flip proves infeasibility unconditionally" | **UNSOUND as an unconditional claim.** (i) opus's S4 script computes the *per-shell* continuum (ARCH) margin (mpmath 25 dps, delta-grid 1/100, x-grid 1/4000) -- it never computes any cross-shell cut; the one cross-shell sentence is asserted, not computed. (ii) When the aggregate/tail cuts ARE computed (audit part 4, exact integers), they do **not** flip at golden: see (G3) row. |
| (G2) decoupling into per-(M,d) tail cuts | The Hall argument (disjoint forward chains, closed sets = unions of tails) is a genuine proof **given** the degree-preservation premise. But the premise ("forward routing preserves absolute degree") is a plausibility statement about the carry mechanism, and it **contradicts** the scale-invariance ("the band rebinds at the deeper shell's own delta*") that opus's aggregate argument implicitly needs. Under preservation the floor argument fails (see next row); under rescaling the decoupling proof fails. Neither transport law is derived. Answer to the audit question: (G2) has a full proof of a conditional statement whose hypothesis is unproved -- i.e. the decoupling *argument* is fine, the *decoupling premise* is a plausibility argument. |
| (G3) "the cuts add nothing" (float computation, shells 8..512) | The equality it reports is real, but the reading is wrong. In exact integers (audit part 4): at `gamma = 71/125` (certified below gamma*), per-shell (ARCH) fails at m = 64, 128, 256, 512, yet of the **819** violated degree-resolved tail cuts over the window, **all 819** live at degrees `d > 255` where only the deepest shell (m = 512) has any demand or supply (`demand_m(d) = supply_m(d) = 0` for `d >= m`), and **zero** have genuine multi-shell content. The "binding degree-resolved cut" is the deepest shell's own per-shell (ARCH) constraint -- a **truncation edge artifact**. Any per-shell deficit interior to the window is absorbed by the next shell at the same absolute degree with exponential room (audit 4.2: at (m, 2m) = (64, 128), d = 36, deficit ~2^55 vs slack ~2^143; the two-shell tail cut is satisfied with ~2^88 to spare). At fixed absolute `d`, demand `binom(m-1,d)` grows polynomially in `m` while supply grows exponentially, so in the model as stated (unbounded forward routing) **every** tail cut with a deeper shell available is satisfied for any `gamma > 0`: the transportation relaxation yields no general-class floor at all. (The exact integer-profile per-shell verdicts fail earlier in m than THM-3024's float continuum control -- floors lower supply -- which changes nothing above.) |

**Consequence.** `C*_general = C*_block = log_5(5 phi^2)` is unsupported.
The gap is not cosmetic: the model itself, taken at its word, is *feasible
below golden*, so no amount of extra computation inside this model can
restore the theorem. What is needed is a genuinely new ingredient -- e.g. a
deadline-bounded routing window derived from the extractor axioms (a deficit
at shell m must be settled within finitely many shells, with the window width
tied to C), or the joint (all-degrees-simultaneously) constraint that the
flow relaxation discards. Until then the correct statement is:

* balanced-block class: `C* > 1.5970` (exact, m = 4096), asymptotically
  `>= log_5(5 phi^2)` modulo the transfer lemma;
* `H = 1` checkpoint-closure programs: barrier `log_5(5 phi^2)` (same caveat);
* general exactly fair extractors: floor OPEN (only the trivial bounds).

## 2. Answers to the three audit questions

**(a) Is opus's degree-blind Hall-cut sign flip at gamma = golden proved
symbolically or only computed in floats?** Floats only -- and it is not even
a Hall-cut computation. `amm12592_cross_shell_capacity_opus_S4.py` evaluates
the *per-shell* continuum margin `min_delta [supply - H(delta)]` in mpmath
(dps 25) on a delta-grid of step 1/100 and an x-grid of 1/4000; its printed
margin at golden is `+1.17e-5` (grid artifact; at the exact binding point it
prints `-1.09e-7`). The location of the sign flip is pinned symbolically only
through THM-3009 sec 3.1 / THM-3027's tangency collapse, which this audit
re-derived and confirmed. The cross-shell content of opus's file is one
uncomputed sentence.

**(b) Is the Laplace/Stirling rate passage in THM-3027 rigorous?** No -- it
is explicitly assumed ("standard... ASSUMED, not reproved"), and it does need
an error-term lemma. After this audit the state is: the inner max over sigma
is concave (proved, so no uniformity-over-sigma issue survives at the rate
level), `sigma*` is interior (now exact), and the floor direction needs only
the *failure* half of the passage, which is a routine finite-Stirling
sandwich (`2^{nH} / (n+1) <= binom <= 2^{nH}`, at most R terms, floors cost
`O(log R)`) -- sketched in THM-3009 sec 10.3 but never written as a lemma
with constants. The *ampleness* half (`gamma > gamma*` => criterion holds for
all tau and all large R) needs in addition a certified global-tau statement
(currently a 3000-point scan with no Lipschitz bound) and boundary uniformity
as tau -> 0, 1. Severity: LOW for the floor, MODERATE for any claim that the
criterion's threshold is *exactly* gamma*.

**(c) Does THM-3024's (G2) forward-routing decoupling have a full proof?**
The combinatorial half (given the model: disjoint chains => the Hall family
is exactly the per-(M,d) tail cuts) is a correct standard proof. The premise
that routing is forward-only at *preserved absolute degree*, with those
supplies and demands, is a plausibility argument about the carry mechanism --
the file itself calls it a "modelling premise" while the status line says
PROVED (structural). Worse, the fallback sentence in its Scope ("if degree
could shift... (G1) already bounds that from below unconditionally") fails
because (G1) is itself unsound as an unconditional claim (see table). So:
no, (G2) does not close the general class.

## 3. The exact-arithmetic re-referee (decisive computations)

**gamma* bracket (pure integer certificates).** `p/q > gamma*` iff
`5^p > phi^{2q}` iff `X := 2*5^p - L_{2q} > 0` and `X^2 > 5 F_{2q}^2` -- a
Fibonacci/Lucas integer comparison. Certified:

```text
115939/193882  <  gamma*  <  105183/175895        (width 1/34102874390 ~ 2.9e-11)
1/2, 2457/6592, 2457/4135, 149/250, 597987/10^6   all BELOW gamma*  (exact)
3/5, 299/500, 597988/10^6                         all ABOVE gamma*  (exact)
```

**(ARCH) sign-flip ladder, exact integers, independently recomputed.**
Candidate ladders verified monotone by exhaustive sweep at m = 8, 16, 32
(no pass-then-fail anywhere), validating the binary search. Then:

```text
m      largest refuted C     gamma_m      vs gamma* (exact)   gap
8      3/2                   0.500000     BELOW               0.09799
16     14/9                  0.555556     BELOW               0.04243
32     64/41                 0.560976     BELOW               0.03701
64     115/73                0.575342     BELOW               0.02264
128    239/151               0.582781     BELOW               0.01521
256    224/141               0.588652     BELOW               0.00933
512    508/319               0.592476     BELOW               0.00551
1024   1992/1249             0.594876     BELOW               0.00311
2048   3890/2437 (verified)  0.596225     BELOW               0.00176
4096   6709/4201 (verified)  0.597001     BELOW               0.00099
```

All ten refuted rates certified `< gamma*` by the integer test; the ladder is
strictly increasing; every committed value reproduced (m <= 1024 by full
independent binary search; 2048/4096 by exact re-verification of the refuted
candidate, its successor candidate passing, and `C = 8/5` passing). Exact
finite-m flip brackets: `(3890/2437, 423/265]` at m = 2048 (first failure at
d = 1262, d/m = 0.6162) and `(6709/4201, 1704/1067]` at m = 4096 (first
failure at d = 2527, d/m = 0.6169 vs 1/phi = 0.618034); the m = 4096 flip
sits `0.00099` below gamma*, consistent with the fold prediction
`~ c log2(m)/m`, and the binding-degree fraction converges to `1/phi`.
D0-robustness: at m = 2048, `gamma = 59/100` fails (ARCH) with the floor
profile `a_k = floor(gamma(m+k)) + D0` at D0 = 0 (d = 1165) **and** D0 = +2
(d = 1173), so the refutation is not an artifact of the profile offset
convention.

**Tangency algebra (sympy 1.9).** Eliminating the multipliers from
(S)+(T)+(V) collapses, with residuals identically zero, to
`(1-tau)^2 = tau`; the root in (0,1) is `tau* = (3-sqrt5)/2 = phi^-2`;
back-substitution gives `u* = 1/sqrt5` (via `2+phi = phi sqrt5`, exact) and
`gamma* = log(phi)/log(sqrt5)` exactly. Confirmed b-universal.

## 4. Demotions / edits requested

1. **THM-3024**: demote status from "PROVED (structural) + VERIFIED-NUMERIC"
   to "MODEL-CONDITIONAL; floor promotion REFUTED within the stated model
   (see audit)". The headline `C*_general = log_5(5 phi^2)` should move back
   to hypothesis status. Any downstream text (including the AMM submission
   writeup, if it cites THM-3024 for the general class) must re-scope the
   lower bound to balanced-block schemes.
2. **THM-3009**: no demotion; annotate that (i) the sec 10.3 transfer lemma
   is the one unwritten analytic step of the asymptotic floor, and (ii) after
   this audit the delta-scan is only needed for exactness-from-above.
3. **THM-3027**: no demotion; annotate the concavity + gamma-monotonicity
   upgrades and the exact sigma* certificate from this audit.
4. **THM-3017**: no demotion; mark "AWAITING INDEPENDENT HOSTILE AUDIT" as
   discharged by this audit (within its stated scope).
