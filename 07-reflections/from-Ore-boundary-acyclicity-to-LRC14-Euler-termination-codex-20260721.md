# From Ore boundary acyclicity to LRC14 Euler termination

**Session:** codex-2026-07-21-DC2-LRC14-termination  
**Status:** three proved no-go/structure theorems, one exact carrier audit, and
an explicit finite-atlas target; neither DC(2), planar JC(2), nor LRC(14) is
claimed solved.

## Outcome first

The attempted DC(2) obstruction at boundary grade six disappears when both
members of the Weyl pair are allowed to move. THM-2049 proves that the
associated-graded correction map is onto in every grade, and the exact Ore
calculation raises the residual through grades `6,7,...,14`. A formal
beta-adic solution of `[S,T]=1` exists; the real gate is finite polynomial
termination and compatibility with `D`.

The same distinction is load-bearing in LRC14. THM-2050 proves that AP13 and
the strict row `12->26` have identical complete local phase-height germs near
all six period-14 unit points. Local exactness cannot decide global loneliness.
The companion audit shows what repairs that loss and what does not:

```text
local unit germ       preserves period-14 tight contact; loses magnitude/exit
G_{1/14} topology     exact weak-LRC carrier; sees AP's isolated Euler points
peel quantization     sufficient strict-volume pruning; not a complete sidecar
resolved first exit   complete finite strict certificate once M>1/14
```

After a fresh mainline pull, the audit also independently reproduced
THM-2048's genuine Cover14 gain at peel `93`: core
`mu=35517/280280`, `theta=19801/40040`, `r=50`, and positive tax excess
`2413467317/235670635200`. This is the right success case for the valuation
sidecar, while the missed smaller rows remain the right controls.

## The transfer, without identifying the problems

The common architecture is operation-level:

```text
DC(2) Ore chart                       LRC14 phase-height chart
-------------------------------       --------------------------------
localize at x!=0                      localize at period-14 top vertices
standard Weyl pair [ell,t]=1          identical binding-tent germs
beta=v_x-2deg_ell sidecar             owner magnitude/resolved-q sidecar
associated-graded map is onto         local germ does not determine M(S)
finite support still open             global Euler/first-exit closure open
```

No algebra map between the Weyl algebra and the runner circle is asserted.
The transferable lesson is narrower: a localization proof needs a valuation
sidecar and a termination theorem. In both settings, producing corrections at
every grade/local chart is weaker than producing one finite global object.

The incoming GMC(2)/NC2 work adds a second precise division. The constant-term
functional controls the angular/volume side and can certify strict exits, but
tight LRC witnesses can have zero Haar volume. THM-2047's Euler characteristic
sees those isolated points. Thus the honest synthesis is

```text
GMC/CT moments + peel discrepancy  -> strict-volume branch,
signed phase-height Euler word     -> tight boundary branch,
Noetherian labelled deletion       -> missing glue.
```

## Assumptions challenged by exact controls

1. **"The first invariant boundary grade is the obstruction."** False.
   THM-2049 makes grade six exact and every subsequent graded equation
   solvable.
2. **"A complete top germ should determine the global row."** False.
   THM-2050 gives identical germs with maxima `1/14` and `1/12`.
3. **"The new peel tax can be the scalar descent height."** False on the
   named controls. It fires for `12->96` and `12->84`, but not for `12->26`,
   `12->36`, or P10+K33.
4. **"More local jet depth repairs magnitude loss."** False by THM-2043's
   infinite Hasse-indistinguishable `q=41` family.
5. **"GMC-style positive norm/volume sees weak LRC."** False at equality:
   AP has six isolated safe points and zero safe volume.

## Exact finite strict branch

THM-2047 gives a useful termination fact. If `M(S)>1/14`, a maximizing reduced
phase `a/q` is supported on an opposite-sign pair and

```text
q | v_i+v_j <= 2 max(S).
```

Therefore the first strict primitive exit satisfies `q_exit<=2 max(S)`. The
strict branch needs no conjectural denominator cutoff. What remains is the
zero-margin branch: prove an Euler point survives when no positive pair-sum
margin is present.

## Next decisive target

For a primitive pair-sum phase set

```text
C_(q,a)(S)=min_v (14 min(av mod q,q-av mod q)-q).
```

The exact Wall-A target is a labelled deletion trichotomy for
`S=C union {w}`:

1. a THM-2048 fiber-tax inequality fails, giving positive safe volume;
2. some `C_(q,a)>0`, giving a resolved strict exit; or
3. an owner-labelled endpoint of `G_{1/14}(C)` is not deleted by the danger
   arcs of `w`.

The first two clauses are exact finite certificates. The third is the only
load-bearing open clause. It should be attacked as positivity of a signed
cyclic endpoint word/Euler valuation, with magnitude retained on every owner.
An unlabelled toric layer, scalar energy, or fixed residue packet cannot state
the required survival event.

The extensive repo pull identifies HYP-2108 as the exact prior-art interface.
For a peel `S=C union {w}` and a lifted core-safe component of midpoint `m_i`
and length `l_i`, containment in one open `w`-danger arc is equivalent to

```text
||w m_i||+(w/2)l_i<1/14.
```

Thus weak endpoint survival is exactly
`P_w(C)=max_i(||w m_i||+(w/2)l_i-1/14)>=0`. The active-owner circuit theorem is
no longer a vague incidence request: it must force `P_w(C)>=0` on THM-2051's
bounded circuit hyperplanes after the tax and positive-margin branches are
removed. HYP-3117/HYP-3120 independently name endpoint owner as the missing
proof-circuit input.

THM-2051 now proves the structural cut in stronger form: after paying pair
covariances exactly, absence of a bounded support-`3..5` relation of height
`2^21` forces positive safe volume. THM-2052 then uses the finite-height theorem
and a pigeonhole code to show every hypothetical counterexample already has
eleven independent bounded three-support relations, hence lies in a finite
atlas of rational subspaces of dimension at most two. Its two-anchor refinement
makes these finitely many one-projective-parameter star families.

THM-2053 closes the infinitary part of that branch. If `u,z` is a saturated
basis of a two-plane, `L=max_i||(u_i,z_i)||`, and `v(a,b)=a u+b z`, then

```text
M_T-max_i|a z_i-b u_i|/(2(a^2+b^2)) <= M(v(a,b)) <= M_T.
```

The centrally symmetric set of signed columns has a nonradial secant by
Ungar's planar direction theorem. Perpendicular projection therefore has full
support and a repeated absolute speed. Settled LRC through thirteen total
runners gives `M_T>=1/13`, so the exact safe gate is
`max_i|a z_i-b u_i|<=(a^2+b^2)/91`. Its round corollary removes every
primitive `||(a,b)||>=91L`. Each rank-eleven star consequently leaves only a
finite anisotropic residual; a twelfth relation is a shortcut, not the unique
finiteness mechanism.

Completing the square makes that residual concrete. Gate failure is the union
over `i` and `sigma=+-1` of the open disks

```text
||(a,b)-(91 sigma/2)(z_i,-u_i)||^2
    < (91^2/4)(u_i^2+z_i^2).
```

All `26` boundary circles pass through the origin. The faithful next vertex is
therefore a primitive lattice point in a tangent disk intersected with the
star's positivity cone, not a raw relation matrix or even the whole round
parameter ball.

The one-tail plane calibrates the scale exactly. With
`u=(1,...,13)`, `z=e_12`, and replacement `w=12+b`, the determinant maximum is
`13|b|` for `b!=0`; the gate is `b^2-1183|b|+1>=0`, hence every integer
`|b|>=1183` is certified. HYP-2896's divisibility fan then discharges the
bounded core (and in fact the whole family) by exact `q=12`, `q=14`, or affine
binding phases. This validates the proposed geometry/arithmetic handoff.

The enlarged audit immediately challenges an overoptimistic reading of that
reduction. Every control already has a height-one support-three circuit; AP and
the six AP-derived strict controls all share `1+2=3`. The genuine Cover14 tax
row has `1+11=12`. Raw circuit existence therefore merges tight and strict
routes even more badly than the local germ. A useful circuit must meet the
peeled speed or an active endpoint owner, or change the signed cyclic wall
word. The next target is an active-owner circuit lemma, not a census of more
relations.

The resulting decisive target is now **finite fan-or-Euler discharge**. Compute
a reduced saturated basis for each two-anchor template, enumerate only the
primitive parameters failing the quadratic gate inside `||(a,b)||<91L`, and
quotient only by data preserving pair-sum
margins and endpoint ownership. HYP-2896 supplies the worked model: its
one-tail family splits at the divisibility walls `12|w` and `14|w` into cells
with a `q=12` witness, a `q=14` witness, or one affine binding phase. HYP-8845
adds a clean equivariant dividend on covering rows: one half-circle survivor
automatically has a mirror partner, so `chi>=2`.

This also repairs a tempting misuse of HYP-4346. An unpointed rigidity theorem
that finds *some* safe direction in a plane says nothing about the specified
counterexample direction. The uniform geodesic estimate has the right
quantifier: it handles every specified direction outside the finite disk.

A final incoming pull supplies the natural analytic partner, THM-2054
(currently claimed/in progress): relative Fejer whole-product decorrelation
along a character line. If completed with explicit constants, it should erase
whole off-resonance cells of the tangent-disk atlas, leaving exact pair-sum and
Euler work only on resonance fans. HYP-8860's Paley-prime table can choose the
useful modulus (`7` for the period apex, `3` for the smallest resonance atom,
`11` for the rank scale), but its tournament orientation is only a lens and
cannot replace signed coefficient heights or endpoint owners.

## Tournament methodology

The audit deliberately uses proof carriers as vertices. Pairwise comparison
uses `(cross-route merges,-fiber count,cost)`, and the switch quotients the
named row bank by each carrier signature. The tournament is transitive with
score histogram `0,...,7`, no directed triangles, singleton SCCs, and one
Hamiltonian path. The ranking is diagnostic rather than a proof: threshold
topology and the fused termination carrier lead, while raw unit germs finish
near last and raw unit-circuit existence finishes last. The quotient preserves
certificate discrimination on the named bank and
destroys the circle geometry, which remains in the signed phase-height word.

## Artifacts

- THM-2049: Ore boundary correction acyclicity.
- THM-2050: period-14 local-germ blindness.
- THM-2053: every rank-two parameter plane has an explicit finite geodesic
  terminal.
- HYP-8841: Noetherian first-exit/Euler termination target.
- HYP-8846: pointed finite-plane completion program.
- `dc2_ore_descent_codex_20260721.py` and matching output.
- `lrc14_termination_sidecar_codex_20260721.py` and matching output.
