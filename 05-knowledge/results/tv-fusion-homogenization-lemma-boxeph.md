# TV-fusion / homogenization lemma at n=2: closed form, fusion inequality, rigid face, transfer probe

- **Agent:** boxeph, 2026-08-03.
- **Script:** `04-computation/tv_fusion_bin2_verification_boxeph.py`
- **Output:** `05-knowledge/results/tv_fusion_bin2_verification_boxeph.out` (ALL CHECKS PASS)
- **Source of the block statement:** Kontorovich, *TV Homogenization Inequalities*, arXiv 2601.04079v3.
- **Method discipline:** exact only — sympy 1.9 polynomial identities with sign-resolved closed
  regions, `fractions.Fraction` scans; no floats in any decision, no numpy.

## Statuses (literal, per claim)

| claim | status |
|---|---|
| Closed form `TV(Bin(2,x),Bin(2,y)) = \|x-y\|(1+\|x+y-1\|)` on `[0,1]^2` | **VERIFIED** (symbolic region-resolved identities + exact 33×33 Fraction grid vs raw definition) |
| Fusion inequality `TV(Bin(2,pbar),Bin(2,qbar)) <= delta_1+delta_2-delta_1*delta_2` on `[0,1]^4` | **VERIFIED** (two-case polynomial certificates, manifestly nonnegative term-by-term, + exhaustive 17^4 exact grid + 20000-pt ragged-denominator hostile scan: zero violations) |
| Rigid equality face = `F0 ∪ F1 ∪ F2 ∪ F3 ∪ F4` (stated below) | **VERIFIED** (forced by vanishing of every certificate term; grid equality set = predicted set pointwise, 351 = 351, both directions; cross-checked by inclusion–exclusion 289+68−6) |
| General block statement `delta_N <= delta_I + delta_J - delta_I*delta_J` for fusing blocks `I ⊔ J` (arbitrary sizes) | **OPEN — this is Kontorovich's stated CONJECTURE, not his theorem** (scope correction 2026-08-03 after full-text read, see `gvc-tv-provenance-hunt-boxeph.md`: the paper PROVES only Lemma 1.4 `delta_N <= 2(delta_I+delta_J)` and remarks the product-form sharpening is conjecturally optimal with "no pathway" via its methods; "rigid face" is not the paper's term). **The `\|I\|=\|J\|=1` certificate above is therefore the FIRST PROVED CASE of that open conjecture**, strictly stronger than Lemma 1.4 at n=2. |
| Transfer to AMM 12592 exact-extraction deadline lane | **REFUTED as a transfer** — verdict NO TRANSFER (six-field spec below) |

## 1. Closed form (proof sketch)

`Bin(2,x)` has masses `((1-x)^2, 2x(1-x), x^2)`. With `s = x+y`, `d = x-y` the three
absolute-value arguments factor as polynomial identities (verified by `expand`):

```text
x^2 - y^2           =  d*s
x(1-x) - y(1-y)     =  d*(1-s)
(1-x)^2 - (1-y)^2   = -d*(2-s)
```

On the box, `s ∈ [0,2]`, so the factors `s` and `2-s` are nonnegative throughout and only
`sgn(d)`, `sgn(1-s)` vary. On each of the four closed regions `{d≥0 | d≤0} × {s≤1 | s≥1}`
the sign-resolved sum is a polynomial identity:

```text
TV = (1/2)|d| ( s + 2|1-s| + (2-s) ) = |d| (1 + |s-1|) = |x-y| (1 + |x+y-1|).
```

All four region identities verified in sympy; independently the closed form equals the raw
TV definition on the exact grid `{0,1/32,...,1}^2`.

## 2. Fusion inequality (proof sketch, n=2)

Write `d_i = p_i - q_i`, `delta_i = |d_i|`, `S = pbar + qbar ∈ [0,2]`. By the closed form,

```text
LHS = |d_1 + d_2|/2 * (1 + |S-1|),      RHS = delta_1 + delta_2 - delta_1*delta_2.
```

**Symmetry reduction (verified symbolically).** The box automorphisms
`sigma` (block swap), `pi` (`p ↔ q`), `kappa` (complement `x → 1-x`) each preserve LHS and
RHS; their action on `(d_1, d_2, S)` is `(d_2,d_1,S)`, `(-d_1,-d_2,S)`, `(-d_1,-d_2,2-S)`.
A finite orbit enumeration (verified) maps every open sign/branch region into one of two
canonical closed cases; boundaries lie in region closures and the inequalities are weak.

**Case B** (`d_1,d_2 ≥ 0`, `S ≤ 1`; parametrize `q_i = c_i ≥ 0`, `d_i = t_i ≥ 0`):
here `1+|S-1| = 2-S ≤ 2 - (t_1+t_2)/2` (since `S ≥ (t_1+t_2)/2`, the "`2 - something`"
step), and with `T = t_1+t_2`, AM-GM gives `LHS ≤ T - T^2/4 ≤ T - t_1 t_2 = RHS`.
Exact certificate (sympy `expand == 0`):

```text
RHS - LHS = (t1+t2)(c1+c2)/2 + (t1-t2)^2/4        >= 0 term by term.
```

**Case A** (`d_1 = a ≥ 0 ≥ d_2 = -b`, WLOG `a ≥ b` via `sigma∘pi`, `S ≤ 1`;
parametrize `q_1 = u ≥ 0`, `p_2 = v ≥ 0`, so `a ≤ 1`):

```text
RHS - LHS = b(1-a) + b + (a-b)(u+v)/2 + (a-b)(a+b)/4   >= 0 term by term.
```

Both certificate identities verified exactly in sympy; nonnegativity of each term is
structural on the closed case region (products of nonnegative box quantities — no floats,
no numerical optimization). Exhaustive exact scan over `{0,1/16,...,1}^4` (83521
quadruples, LHS computed from the *raw* TV definition, Fractions): **zero violations**;
plus 20000 random quadruples with denominators up to 64: zero violations.

Equivalent complement-product form: `1 - TV(Bin(2,pbar),Bin(2,qbar)) >= (1-delta_1)(1-delta_2)`.

## 3. The rigid equality face (exact characterization)

Equality `LHS = RHS` holds **precisely** on

```text
F0 = { p1=q1, p2=q2 }                                   (both sides 0)
F1 = { p1=p2=t, q1=q2=0 },   F2 = { p1=p2=1, q1=q2=t },
F3 = { p1=p2=0, q1=q2=t },   F4 = { q1=q2=1, p1=p2=t },   t ∈ [0,1].
```

Equivalently: equality is nontrivial (`delta_1 delta_2 > 0`) **iff the two blocks carry
identical pairs** (`p1=p2`, `q1=q2`) **and one of the two fused distributions is
deterministic** (`Bin(2,0) = δ_0` or `Bin(2,1) = δ_2`), where `1 - TV = (1-delta)^2`
exactly. Derivation: the certificates are sums of nonnegative terms, so equality forces
every term to zero — case A collapses to `a=b=0`; case B forces `t_1=t_2` and
(`t=0` or `c_1=c_2=0`), i.e. the face `{p1=p2=t, q1=q2=0}`, whose `G`-orbit is exactly
`F1..F4`. Grid confirmation: the 17^4 scan's equality set equals the predicted face set
pointwise in both directions (351 = 289 + 68 − 6 points, inclusion–exclusion checked).

This is the "rigid face" of the homogenization step: fusing two blocks is TV-lossless
against the product bound only at these boundary configurations.

## 4. Transfer probe (HYP-9061 sec 2e "couple two biases"; THM-3024 Hall floor) — six-field spec

**Candidate transfer:** use `1 - TV >= (1-delta_1)(1-delta_2)` (fusion of
indistinguishability across independent flip blocks) to lower-bound stopping data /
deadlines for exact fair-coin extraction (AMM 12592 lane).

1. **Map.** Flip blocks ↦ `Bin(2,·)` fused blocks; the two certificate biases
   `q_A = 896/2181`, `q_B = 2974400/11821757` (HYP-9061 sec 2e) ↦ a pair
   `(p_i, q_i)` with `delta = |q_A - q_B| = 4105127872/25783252017 ≈ 0.1592`;
   "couple two biases" ↦ the complement-product fusion bound.
2. **Preserved predicate.** Only this: *two independent biased blocks, homogenized to
   their average bias, gain no distinguishability beyond the product law* — an
   archimedean/L1 statement about distributions up to outcome relabeling.
3. **Loss.** Three independent kills. (i) *Exactness:* the lane's fairness constraint is
   an exact polynomial identity over Q (ledger of THM-3002/THM-3024); TV is blind to it —
   leaf-probability denominators can change arbitrarily at TV distance 0. (ii) *Valuation
   data:* sec 2e's genuine dual must cancel archimedean mass between `N^(A), N^(B)` while
   keeping 2-adic leading terms **misaligned** (`s_A = 7` vs `s_B = 6`); TV is invariant
   under relabeling and carries no p-adic structure, and the fusion RHS is symmetric in
   `(delta_1, delta_2)` while the needed coupling is intrinsically asymmetric.
   (iii) *Direction/typing:* the lemma **upper-bounds** a TV between two four-parameter
   distribution pairs; sec 2e needs a **lower bound** on a two-evaluation-point dual
   functional of one generating object — an arity/typing mismatch: the lane has two
   biases, the lemma consumes two *pairs*.
4. **Sidecar that would be needed.** A valuation-aware fusion law — a "2-adic TV" on
   leaf-probability ledgers with a complement-product inequality. No such object exists
   in canon, and TV's relabeling invariance suggests none factors through it.
5. **Cheapest hostile test (run, exact).** Instantiate at the actual certificate biases:
   `TV(Bin(2,q_A), Bin(2,q_B)) = 141573722928218785664/664776084572134568289 ≈ 0.2130`,
   fusion bound `1-(1-delta)^2 ≈ 0.2931`. Both are dimensionless O(1) rationals with **no
   n-dependence**: the lane's frontier is a per-degree *rate* (inequality (27) shape,
   `C*` envelope), and no ledger inequality in THM-3002/3009/3017/3024 mentions any
   TV-like quantity. The numbers make contact with nothing. Also, THM-3024's floor is
   already closed by exact Hall/transportation cuts (combinatorial, unconditional); a
   statistical indistinguishability bound could at best re-derive strictly weaker
   approximate-fairness statements.
6. **Verdict: NO TRANSFER** to the exact-extraction deadline lane. Honest residue worth
   recording: the lemma is the `|I|=|J|=1` case of a genuine product-type family (cf. the
   standard `1 - TV(P⊗P', Q⊗Q') >= (1-TV(P,Q))(1-TV(P',Q'))`), homogenized through the
   bias average — file as an external CITED landmark. Its only conceivable in-repo
   consumer is a hypothetical *approximate*-fairness side question (none currently open);
   it cannot strengthen, simplify, or shortcut the Hall-cut floor or the sec-2e 2-adic
   coupling.

## Reproduce

```bash
python3 04-computation/tv_fusion_bin2_verification_boxeph.py
# expected tail: "OVERALL: ALL CHECKS PASS"; runtime ~13 s.
```
