# Independent referee: the denominator-23 phase bank

**Status: independent analytic/source audit PASS, with independent exact
controls.** The final producer report and input identity are recorded below.
LRC(14) remains open; the result is
a sufficient residue-class condition. It does not require or establish
decoder entry for its generic controls.

## Exact statement and supplier scope

Let `V` contain six positive integer speeds whose residues modulo 23 belong
to `A=+/-{1,3,5,6,10,11}`. Let `U` contain seven positive integer speeds and
put `L=max U`. Let `t,g` be positive coprime integers. Define

```
B={3,5,6,10,11,12,13,17,18,20}/23,
rho(gB)=m(g)/23.
```

If `28 m(g) L <= 69 t`, the full physical row `tV union gU` has a phase
with every clearance **strictly greater** than `1/14`. The arithmetic proof
also works with repeated speeds, although the recorded thirteen-speed
controls are distinct. No primitive-gcd or parity assumption on `V` is used
beyond the displayed residue hypothesis; coprimality of the two scales is
used in the lift permutation.

The only non-elementary mathematical input is a seven-speed phase of
clearance at least `1/8`. This is **CITED**, not reverified globally here:
Sungkawichai and Trakulthongchai,
[Eleven, twelve, and thirteen lonely runners, Theorem 1.3](https://arxiv.org/html/2604.23906v2)
states LRC for at most twelve nonzero integer speeds. I checked the primary
September 1, 2026 version and its convention counting nonzero speeds.
The paper's computational proof is not repeated by the finite bank below.
The repository's maintained larger-unit and nonunit component theorems use
this same proper-component supplier.

## Independent proof

For every `a in A` and `b in B`, literal reduction modulo 23 gives
`||ab||>=2/23`. In fact `B` is the complete set of residues with this property:
zero is forbidden, and the excluded nonzero residues are exactly the
inverses of `A`, because `A=-A`. The common strict margin over the target
is `2/23-1/14=5/322`.

Every lifted phase `x=(k+b)/t`, `0<=k<t`, protects `tV` with that margin.
Write `gb=l_b+gamma_b`, with integral `l_b` and `0<=gamma_b<1`. Its image
in the `U` phase coordinate is `(gk+l_b+gamma_b)/t`. Since multiplication
by `g` permutes residues modulo `t`, the union of all such images is exactly

```
{(j+gamma)/t mod1 : 0<=j<t, gamma in gB mod1}.
```

The largest circular gap of this set is `rho(gB)/t`, including the gap
crossing zero. If `23|g`, all ten residues collapse to one and the gap is
`1/t`, so that residue must not be handled as a ten-point permutation.
More generally, if `d=gcd(t,g)`, write `t'=t/d`, `g'=g/d`. The original
index list repeats the reduced list exactly `d` times, giving the exact
gap `(d/t)rho(g'B)`. The unreduced expression `rho(gB)/t` is false in
general; `t=g=2` is an explicit control. The main theorem uses the coprime
normalization, so this additional identity changes no stated sufficient bound.

Choose a proper-component phase `beta` with every `U` clearance at least
`1/8`. The minimum clearance is `L`-Lipschitz on the circle. Consequently
the closed arc centered at `beta` with radius `3/(56L)` is weakly safe at
`1/14`, and its interior is strictly safe. Its length is `3/(28L)`.
Any closed circular arc at least as long as the largest gap contains a bank
point. Thus the displayed sufficient condition is exactly

```
m(g)/(23t) <= 3/(28L)  iff  28m(g)L <= 69t.
```

If the hit is in the interior, both parts are already strict. If it lies
at an endpoint and one `U` constraint is exactly weak, move the physical
phase slightly in the direction that moves its `g`-image toward `beta`.
The `U` Lipschitz bound becomes strict. A sufficiently small displacement,
less than `5/(322t max V)`, preserves the strict `V` margin. This proves
strict safety **also at equality**. There is no claim that a closed-arc
endpoint itself must already be strict.

The complete independent residue profile, for `g=0,...,22`, is

```
m(g)=23,6,4,10,5,5,6,4,5,7,8,6,6,8,7,5,4,6,5,5,10,4,6.
```

In particular `g=1 mod23` gives `m=6`, and the sufficient condition simplifies
to `56L<=23t`. This simplification includes equality.

## Independent exact controls and scope of the gain

The referee source imports no producer or repository program. It enumerates
the residue bank directly and also obtains it by modular inversion. It
compares literal lifted phase unions against the effective-order gap formula
for all **420** pairs `1<=t<=14`, `1<=g<=30`, including both noncoprime
scales and the collapsed residue.

For every one of the 120 seven-subsets of `{1,...,10}`, it finds a proper
phase from the literal rational walls `(8k+/-1)/(8u)`. This is a finite
existence check on those controls only. With each scale
`g in {1,2,3,7,10,23}`, it chooses the least admissible coprime `t` above
the bank bound, computes a bank lift by modular inversion, and verifies
the resulting full phase directly. The fixed six-shape
`V_i=a_i+23(100+i)`, `a=(1,3,5,6,10,11)`, `i=0,...,5`, is primitive
and has both parities. These **720** actual distinct thirteen-speed safety
controls do not carry a decoder-entry assertion.

A dedicated equality control has

```
V=(1,3,5,6,10,11), U=(1,2,3,4,5,6,23), t=56, g=1.
```

Here `56L=23t`, the proper center is `beta=1/8`, and the actual bank endpoint
`x=41/322=beta+3/(56*23)` has full clearance exactly `1/14`, attained by
speed 23. The phase `50507/396704`, obtained by a small perturbation toward
the center, is strictly safe for all thirteen speeds. This checks the
endpoint argument on actual inequalities, not only on an abstract arc.

Both inherited large balanced rows also satisfy the bank condition. However,
**both primitive shapes in those rows are all odd**, including the row with
an even physical scale `t=36883259192`. Hence they already admit a stronger
actual-clearance single-phase proof. Both shapes have clearance `1/2` at
their half phase. The `U` interval then has length `6/(7L)`, and `7L<=6t`
suffices. For the mixed physical row a direct single-phase witness is

```
x=36883259191/73766518384,
clearance=24879646471/73766518384 > 1/14.
```

The odd physical row is safe at `1/2`. These computations prevent either
inherited row being described as inaccessible to actual-clearance
single-phase methods. The denominator-23 lemma's surviving scope is its
arbitrary-parity residue class and scale bound, not optimality of the chosen
bank for the inherited odd star. Denominator-12/14 odd-only improvements
require their own stated assumptions and are not needed in this audit.

The producer's separate nonunit seven-shape is also checked independently:
`U=(2,3,4,5,6,7,12003612719)`, with the inherited `V`, `g=1`, and
`t=36883259192`. It is primitive and has both parities. At phase
`26509842545/212078740354` the full row has clearance exactly
`17507133005/212078740354>1/14`. This is a physical safety control, with
no inferred decoder equality or claim that every other phase method fails.

The final producer also supplies a distinct **actual** mixed, unitless
`6+7` control. It changes the fifth inherited six-shape coordinate to
`201676672560` and uses `U=(2,3,6,331,109561,36264691,12003612721)`.
For that row the referee reconstructs every atlas edge from the strict
primitive-sum factor criterion in **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
There are eleven edges, exactly the declared components, and the literal
bounded relation matrix has rank eleven by integer row elimination.

All 105 two-core/one-outside supports require a distinguished coefficient
above `Q=91^6`; the exact minimum is `36446911794129044`. For all 126 reverse
supports, even the unamplified target exceeds the maximum two-term budget:
the minimum raw ratio is `778001386374677670/778001386051180109>1`.
The producer's stronger ratio after retaining the least permitted
distinguished coefficient is also reproduced exactly:
`3112005545498710680/3102632035080254131`.
These checks cover support-two crossings through zero pair coefficients.
Thus every bounded support-three relation is internal; the atlas rows
already span the full componentwise kernel. This proves the claimed
`W=V_dec` without confusing graph rank with full entry.

The row is primitive and lies inside the physical `Q^2` box, with exact sum
`43300016482122890820283`. At phase `26509842545/212078740354`, its clearance
is exactly `26509924715/212078740354>1/14`. Both primitive shapes have both
parities and lack a unit. The independently enumerated scalar gcd maxima
for subset sizes 12 down through 7 are `(1,1,1,1,2,3)`. Neither producer nor
referee claims the full joint-shadow profile bank was checked.

The final report's odd-only denominator-12 comparison is correct and is
also checked literally. Its recovered seven-clock corollary is elementary:
if `7` does not divide `t`, select `k` with `tk=3 mod7`; the nearest half
lift differs from `k/7` by exactly `1/(14t)`. Odd `V` remains at clearance
`1/2`, while `U` units modulo 7 with `max U<t` remain strictly above `1/14`.
This preserves the favorable translate rather than replacing it by a
worst-gap bound. The referee checks these identities for `1<=t<=99` and
all admissible positive `u<t`. The named incoming `5+8` numerical row is
an inherited application, not an additional actual-entry audit here.

## Reproduction and frozen referee evidence

```
python -B 04-computation/continuing1_20260906_phase_bank_audit.py
python -B -O 04-computation/continuing1_20260906_phase_bank_audit.py
```

Both runs pass **20,205 always-active gates** with byte-identical LF output.
The script is standalone and requires only Python's standard library.

```
source 4986aa7f5a23b1ac69e550dcb22f39a3fb84ed97906a95902835a251b1dd25a5
output b55d2d450bd2004bcf41f2f8f6f1eb5206c447073fd66b3af20ae1aed51a1594
```

No producer files, repository files, or Git state were changed.

## Producer review record

I read the final main report in full, and the producer source after completing
the independent core verifier. The proof, strict endpoint, noncoprime
extension, actual mixed entry, and novelty boundaries are accepted.
The producer is frozen in `C:/w/continuing1_20260906_lrc` with stem
`continuing1_20260906_lrc_phase_bank`:

```
source 5f4cf70ec86a5b65b392b5b2994afd5b1b0821672113de5eace9d2643d712872
output 7223e6ab65192653bb492dbe77b38ce894541b261243951d293caaa6a8ef2006
```

The report's `+/-` residue descriptions have a harmless encoding artifact
in a few prose lines; the literal residue list, proof and source are
unambiguous. Correcting that typography changes no audited mathematics.
