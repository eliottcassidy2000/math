---
id: HYP-2909
status: PROOF-TEMPLATE AUDIT / apex-shell correction integrated; residual sharpened
source: codex-2026-06-22-S119
tags: [lrc14, thm079, tight-locus, apex7, goddyn-wong, covering-sets, tournaments, open-q-108]
depends_on:
  - HYP-+2878
  - HYP-2877
  - HYP-2879
  - HYP-2905
  - HYP-2906
  - HYP-2907
  - HYP-2908
  - HYP-2896
  - THM-079
  - THM-201
  - THM-343
  - THM-523
  - THM-524
  - THM-527
  - THM-560
  - THM-565
  - THM-566
  - THM-568
related:
  - KPS-S31y
  - KPS-S31z
  - KPS-S31aa
  - KPS-S31ab
  - HYP-+2909
  - HYP-2911
  - HYP-2929
  - mac-mini-S47
  - mac-mini-S48
  - OPEN-Q-108
results:
  - 04-computation/lrc14_thm079_star_audit_codex_s119.py
  - 05-knowledge/results/lrc14_thm079_star_audit_codex_s119.out
  - 04-computation/lrc14_gap_M_exact_mac-mini-2026-06-16-S3.py
  - 05-knowledge/results/lrc14_gap_M_exact_mac-mini-2026-06-16-S3.out
  - 04-computation/lrc14_tight_locus_census_kps.py
  - 05-knowledge/results/lrc14_tight_locus_census_kps.out
  - 04-computation/lrc14_M_bounded_sweep_floor_mac-mini-2026-06-16-S3.py
  - 05-knowledge/results/lrc14_M_bounded_sweep_floor_mac-mini-2026-06-16-S3.out
---

# HYP-2909: the THM-079 template for LRC14 reduces to STAR, then THM-568 sharpens the residual

The owner's THM-079 analogy is the right proof shape, but it should be stated
with one extra guardrail.  The analogy does not by itself prove the bounded
atom; it identifies exactly what theorem about the atom is still missing.

Post-pull integration plus S120 correction: incoming THM-568 proves the local
apex-shell half of the STAR target for tight rows (`14|D` and
`D|(active pair sum)`), but not by itself primitive denominator `D=14`.  So this
note should now be read as an audit of the THM-079 template plus a handoff to
the shell-collapse / covering-strictness residual:

```text
14-covering primitive reduced atom  =>  M(S) > 1/14.
```

The structural statement "primitive tight implies apex denominator `D=14`" is
therefore still a theorem target unless supplied by an additional shell-collapse
or covering-strictness argument.  The remaining analytic/combinatorial burden is
covering strictness, especially when at least seven runners are multiples of
14.

Final-rebase integration: incoming mac-mini S49 / HYP-+2909 proves the weaker
forward binding-pair statement

```text
M(S)=1/14  =>  some active pair satisfies 14 | (s_i+s_j).
```

This is exactly the safe local conclusion of THM-568 as corrected in HYP-2929
and formalized in `LRCApexShell.lean`.  The stronger denominator-`14` conclusion
requires the separate shell-collapse step.

## The template

THM-079 forbids `H=21` by reducing the OCF value to one hard component and
then forcing that component to generate extra odd-cycle mass.  The LRC14
analogue is:

```text
Move A: peel/reduce every non-atomic row
        by omit-prime witnesses, dilation, and large-speed boundary-state
        induction.

Move B: prove the remaining atom cannot be a covering counterexample
        because the only sharp rows are the AP and Goddyn-Wong boundary
        tilers, both apex-denom-14 and non-covering.
```

HYP-2905/HYP-2906 make Move A theorem-shaped, with THM-524/THM-527/THM-565
providing the binding-pair, density, and apex-ruler machinery around the
existing reductions.  HYP-+2878 is the older strong-atom covering route, and
HYP-2877/HYP-2879 are the tournament strong-component/strong-ear templates.
THM-523 says a counterexample must be covering, hence must contain a multiple
of `14`.  Therefore every denom-14 apex point is blocked in a covering row:

```text
s = 14m in S  =>  ||s * (a/14)|| = 0.
```

So if sharpness is known to occur only at denom-14 apex witnesses, covering
rows cannot be sharp.

## STAR before and after THM-568

The exact theorem needed is not merely a slogan that AP/GW are tight.  Before
THM-568, the template needed one of the following two forms.

### STAR+ (direct sharp extremality)

For every primitive reduced 13-speed atom,

```text
M(S) >= 1/14,
```

with equality only for the AP/Goddyn-Wong tight locus, up to the normalizations
allowed by the reduction, and with an apex denom-14 maximizing point for every
equality row.

This is the cleanest statement: it is exactly the bounded-atom theorem.  Then
THM-523 removes non-covering rows by q-witnesses, and covering rows are strict
because the only equality rows are non-covering.

### STAR0 plus boundary forcing

If one proves only

```text
M(S) = 1/14  =>  S has an apex denom-14 optimum
                 and S lies in the AP/GW tight locus,
```

then one still needs a separate reduction/compactness statement:

```text
any hypothetical sub-tight covering family with M(S)<1/14
  produces a reduced boundary atom with M=1/14.
```

Without that boundary-forcing step, equality classification alone does not
logically rule out a non-apex sub-tight covering row.

This is the main leak to avoid when using the THM-079 analogy.  THM-079's
component proof has the forcing theorem inside the tournament category.  LRC14
still has to prove either STAR+ directly or prove the boundary-forcing bridge
that makes STAR0 sufficient.

Incoming THM-568 proves the local apex-shell half of STAR0:

```text
M(S)=1/14 at t=a/D with opposite binders  =>  14|D and D|(u+v).
```

Thus the pasted `M(S)=1/14 => apex denominator` clause should be split into the
proved shell statement and the still-open shell-collapse statement:

```text
primitive tight shell D=14h  =>  h=1.
```

The remaining covering statement is:

```text
14-covering primitive reduced atom  =>  M(S) > 1/14.
```

For 14-free rows, THM-523 already supplies a denom-14 q-witness, so tightness
ensures at least one apex-denom-14 optimum.  For 14-covering rows with at most
six multiples of 14, the S31v comb-teeth union bound over the 14-free core's
margin handles the row.  The live residual is therefore the case with at least
seven multiples of 14 over a small 14-free core, or equivalently an apex-shell
second-moment/equidistribution problem.

## Exact audit

`lrc14_thm079_star_audit_codex_s119.py` checks the local arithmetic:

```text
AP {1..13}               covering=False has14=False M=1/14 tau=5/14
GW {1..11,13,24}         covering=False has14=False M=1/14 tau=5/14
cover C_1={1..11,13,84}  covering=True  has14=True  M=7/89 tau=37/89
cover C_2                covering=True  has14=True  M=14/173 tau=72/173
cover C_6                covering=True  has14=True  M=42/509 tau=212/509
```

The covering rows block every denom-14 unit point, but they are still safely
lonely off the apex at binding-pair denominators.  Therefore raw denom-14
blocking is not the bounded-atom theorem.  After THM-568, the structural
"tight implies apex denominator" part is proved; the remaining theorem is
covering-strictness.

The same audit also verifies the graph atom used by HYP-2908:

```text
connected I(G,2)=7  <=>  G=K_3
```

through `n<=5`, matching the elementary proof in HYP-2908.

## Routes from here

### Route 1: THM-568 residual / multiples of 14

Use THM-568 to avoid full tight-locus classification.  The target becomes:

```text
S = R union Q,  Q subset 14*Z,  |Q|>=7,  R 14-free and |R|<=6
  =>  M(S) > 1/14.
```

The core `R` has a proven `1/13` margin by smaller LRC.  The multiples of 14
must fail to cover the corresponding margin interval.  This is the localized
Node-3 / second-moment residual singled out by S31aa.

### Route 2: tight-locus / three-gap

Prove STAR+ directly by classifying the reduced tight locus.  This is stronger
than what THM-568 now requires, but remains a clean independent route.

Known pieces:

- THM-560 proves the difference-closed exact tilers are exactly AP dilates.
- HYP-2893 identifies the Goddyn-Wong acceleration mechanism and the LRC14 row
  `{1,...,11,13,24}`.
- Exact single-swap and bounded scans keep finding only AP/GW equality rows:
  `lrc14_tight_locus_census_kps.py`,
  `lrc14_tight_locus_and_true_inf_kps.py`, and
  `lrc14_M_bounded_sweep_floor_mac-mini-2026-06-16-S3.py`.
- `lrc14_gap_M_exact_mac-mini-2026-06-16-S3.py` is the exact rational `M(S)`
  engine behind the local equality and second-spectrum checks.

Missing theorem:

```text
No reduced non-AP/GW atom has M <= 1/14.
```

Equivalently, every non-AP/GW reduced atom has a positive off-apex interval or
binding-pair slack.  This is the finite Node-2 theorem behind the pasted
`M(S)=1/14 => apex` statement.

### Route 3: forbidden-H7 state lift

Use HYP-2908 instead of a direct tight-locus classification.  Prove that a
reduced covering atom trying to beat or exactly saturate the apex bound
constructs a tournament-conflict-realizable packet graph with connected
`I(.,2)=7`.  Then the graph is `K_3`, and THM-201/THM-343 forbid it as a
tournament OCF component.

Missing theorem:

```text
LRC14 apex over-cover / sub-tight boundary packet
  -> tournament-conflict-realizable connected activity-2 packet graph
     with I=7.
```

HYP-2907 is the guardrail: generic binary digraphs and incomplete orientations
realize `7`, so the lift has to land in the tournament OCF category or a
packet category satisfying the same closure rule.

## THM-079 analogy: exact and inexact parts

Exact transfer:

- multiplicative/component thinking becomes boundary-state induction;
- a scalar obstruction (`7`) becomes useful only in the right category;
- the final step is a single-atom forcing theorem, not another raw induction.

Inexact transfer:

- THM-079 has actual OCF multiplicativity and proved forcing lemmas;
- LRC14's Move A still carries finite-comb and scale-balance hypotheses;
- after THM-568, LRC14's Move B still needs covering-strictness for the
  multiples-of-14 residual, or a stronger STAR+/HYP-2908 closure.

So the proof should now be advertised as one of:

```text
LRC14 = HYP-2905/HYP-2906 reductions
        + THM-568 apex-denominator lemma
        + ">=7 multiples of 14 over a 14-free core are not tight"
```

or the stronger direct route:

```text
LRC14 = HYP-2905/HYP-2906 reductions
        + STAR+ bounded-atom theorem / tight-locus classification
```

or

```text
LRC14 = HYP-2905/HYP-2906 reductions
        + HYP-2908 tournament packet state lift.
```

The pasted argument is therefore a strong proof template, not yet a closed
proof.  After THM-568, its main contribution is even sharper: it locates the
residual at multiples of the apex denominator rather than across the whole
tight locus.

## Tournament Analysis / assumption challenge

Candidate vertex sets considered:

```text
runners
denom-14 apex points
q-witness denominators
binding-pair denominators
covering obligations q=2..14
AP/GW tight-locus atoms
proof obligations in the THM-079 template
tournament OCF conflict packets
connected graph atoms with I(.,2)=7
```

Chosen carrier for this audit: proof obligations plus graph atoms.  This
preserves the predicate "a counterexample is impossible after the reduction"
and exposes dependency direction:

```text
peel_exhaustiveness -> THM-568 apex denominator
                    -> covering_strictness
                    -> state_lift
```

It destroys raw speed identity, which is acceptable only because the exact-M
rows above are kept as counterchecks.  The challenged assumption is that
denom-14 obstruction alone is the finite atom bound.  The covering family
`{1,...,11,13,84m}` shows otherwise: it blocks the apex, but its witness simply
moves to an off-apex binding pair.
