# LRC14 Density-Denominator Bridge

Codex 2026-06-18.

The prompt asked to push the equivalence

`rho*(P,E) >= c0 > 0`

and

`covering forbids binding denominator D = 14j - r with small r`.

The useful correction is that this bridge has two stages, not one.

## Three Reservoirs

For a mixed S3 set `S=P union {V-e:e in E}` there are three related but distinct objects.

1. `rho*(P,E)`: the limiting shape reservoir

   `x in G_P` and `maxgap({e*x}) > 2/7`.

   This is the THM-527 sufficient route. It says a suitable fast phase exists.

2. Grid opportunities at `q=14V`: points `a/q` sampling the same shape reservoir, with the finite fast phase visible. These track the Weyl/discrepancy tail, but they are not automatically actual CRT witnesses because the finite phase `V*a/q` may miss the available cluster gap.

3. Actual CRT witnesses: residues `a mod 14V` for which every runner of `S` is level-`1/14` safe. These are the literal finite covering-system survivors.

The new script keeps these separate. That separation matters: in the bounded examples the `rho*` grid-opportunity density and the actual CRT witness density are both positive, but neither is a subset of the other at the same residue. The missing proof layer is therefore the integer-vs-real fast-phase alignment, exactly the remaining THM-527-G/OPEN-Q-108 issue.

## Denominator Remainder

At the exact optimum, `M(S)=j/D` and `M>1/14` is the same inequality as

`D = 14j - r`, with `r=14j-D > 0`.

The new scout checked this for the hard bounded families:

- k=5 exact rho* minimizer: `rho*=95/2548`. Actual covering lifts at `V=28,70,140` have global CRT densities about `0.097,0.102,0.102`; M-binding remainder debts `r=36,124,260`.
- k=6 exact rho* minimizer: `rho*=3488/63063`. Covering lifts at `V=112,280` have global CRT densities about `0.078,0.081`; M-binding debt `r=173` at `V=112`.
- Angle-B broad-scan `rho*=1/90` shape: covering lifts at `V=84,140,280` have global CRT densities about `0.073,0.072,0.070`; M-binding debts `r=125,236` on the exact-M rows.
- In an 80-row bounded random S3 census: no LRC breaks; min global CRT density `9/203`; min rho-grid opportunity density `11/728`; smallest M-binding remainder seen was `r=11`, with no `r<=10` rows in the sample.

This does not prove that covering forbids small `r`. It says that small `r` is not a one-pair property: a near-`1/14` binding pair becomes dangerous only if all other runners clear at the same crossing. So the proof should target the whole covering system: conditional blockers cannot erase every rho* opportunity.

## Tournament Lens

I used runner vertices with a conditional marginal-blocking gauge:

for each pair `(v,w)`, hold all other runners fixed modulo `14V`, and orient `v -> w` if `v` blocks more of the remaining candidate residues than `w`; ties follow the increasing-speed Hamiltonian path.

This quotient preserves the CRT covering obstruction but destroys the interval-component geometry. In the named hard families it was transitive with one Hamiltonian path. That is itself information: runner-level blocking dominance is too simple. The next useful tournament vertices should likely be safe components, denominator events, four-window packets, or proof obligations rather than runners.

## New Proof Target

Replace:

`private pair has j >= D/14`

with:

`for every bounded-spread shape and small part P, the family of conditional blocker sets cannot cover the rho* grid reservoir at all sufficiently large V`.

This is a finite/discrepancy covering lemma. It is the cleanest merger I see between the uniform density floor and the small-remainder binding-denominator language.
