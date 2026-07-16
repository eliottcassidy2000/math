# HYP-7084 — the owner packet as a parallel-class current

**Status:** EXACT REDUCTION + SYNCHRONIZED-WALL THEOREM PROVED; RAW
CROSSING-ENERGY CLOSURE REFUTED (codex-2026-07-16-S19).  The noncommon
gcd-lattice wall sum remains open.

This note reserves the exact class-circle reformulation of `HYP-7083`.  On the
one-miss set of `F=E union {t}`, write

```text
q=floor(14{tx})=2c+eta,   eta in {0,1},   u=2c mod 7.
```

Then `sec(tx)=c` while `sec(2tx)=u+eta`.  The two possible target classes are
therefore adjacent vertices `u,u+1` of the `THM-913` parallel-class circle:
the even half-sector is the degenerate pair `(c,c)` of sum class `u`, and the
odd half-sector is the boundary pair `(c,c+1)` of sum class `u+1`.

Let `mu_s^eta(u)` be the exact mass on which `F` has unique miss `s` and the
owner half-sector has class coordinate `(u,eta)`.  Put

```text
M_s(u)=mu_s^0(u)+mu_s^1(u),
J_s(u)=mu_s^0(u)-mu_s^1(u).
```

The exact decomposition is

```text
Delta_F(2t)
 = 1/2 sum_s (M_s(s)+M_s(s-1)) - p1(F)/7
   + 1/2 sum_s (J_s(s)-J_s(s-1)).
```

Thus the parity remainder is the same-label trace of the discrete divergence
of an oriented current on `C_7`.  For the `THM-913` crossing profile

```text
xi=(0,0,2,3,3,2,0),
```

let `L_xi` be its circulant Laplacian and set

```text
E_pc=sum_s [M_s^T L_xi M_s + J_s^T L_xi J_s].
```

The adjacent-current and adjacent-average dual norms are exactly `6/29` and
`22/203`, giving the proved comparison

```text
Delta_F(2t)^2 <= (16/29) E_pc.                                (3)
```

The exact `HYP-7083` bank verifies the decomposition, the polarized energy
identity, and (3) on all `6,900` doubling rows.  It also refutes raw energy as
the closure: only `5,346/6,900` rows meet the energy threshold needed for
`0.097`; the worst energy-only upper bound is `0.24485`.  More revealingly,
the current term itself stays between `-0.02132` and `0.02222`, while the
symmetric term reaches `0.07171`.  At the dangerous row
`E=(0,3,4,5,6,7),t=8`, the correction splits as

```text
2173/27440 = 5903/82320 + 11/1470.
```

So the finite crux is mostly symmetric adjacent-endpoint aliasing, not raw
half-sector parity.  This quotient preserves the exact correction and its
miss-label incidence, but destroys slow-wall order, the additive relation
lattice of `E`, and individual chord crossings.  The next sharpening restores
those sidecars by subtracting the closed `THM-891` residue-two limit and writing
the residual as an exact slow-wall endpoint potential on the fourteen-cycle.

## Exact endpoint palette

Let `m_E(x)` be the slow core's missed-sector set.  For a fixed pattern `m`,
define the fourteen-cell waveform `H_m(q)` by adding a source in sector
`c=floor(q/2)` to `m`: if one miss `s` remains, set

```text
H_m(q)=1_[q mod 7=s]-1/7;
```

otherwise set `H_m(q)=0`.  Let `P_m` be the periodic primitive of
`H_m-mean(H_m)`.  Exact cell integration gives

```text
Delta_F(2t)
 = C_2(E)
   + (1/t) sum_(slow walls p)
       [P_(m before p)({tp})-P_(m after p)({tp})].              (4)
```

Here `C_2(E)=F_2(E)/2` is the already closed residue-two limit of `THM-891`.
Equation (4) is the `THM-727` endpoint sum with a finite parallel-class
palette; it is proved by telescoping `P_m(tx)/t` across the slow-wall cells.
All `6,900` bank rows verify (4) fraction-exactly.  The finite wall remainder
ranges from `-0.03610` to `0.06704`; it is the genuine remaining alias.

## The synchronized deck closes sharply

At `x=k/7`, every slow runner crosses a sector wall simultaneously.  The sum
of these six nontrivial wall contributions depends only on the multiset of the
five slow residues modulo seven and on `t mod 7`.  Exhausting the `462*7=3,234`
residue rows gives coefficient range

```text
-23/98 <= t G_sync <= 229/686.
```

The all-zero residue multiset is incompatible with primitivity.  Every other
five-residue multiset has a primitive distinct-speed realization by diameter
`35`; minimizing the feasible diameter and then the congruent `t>D` proves the
sharp universal form

```text
-23/784 <= G_sync <= 229/5488.                                (5)
```

Both extrema occur already at `D=7,t=8`.  The upper equality is the dangerous
row `E=(0,3,4,5,6,7),t=8`; its synchronized deck supplies `229/5488`, while
the closed limit supplies `2209/144060` and the noncommon wall remainder is
`797/36015`.

After removing (5), the bounded-bank noncommon remainder lies in
`[-0.03275,0.03922]`, with finite maximum `269/6860`; this is evidence, not a
universal bound.  Its owner sets are precisely gcd-sheet intersections.  This
is the real-space finite-alphabet face of the general-cluster `THM-887`
per-owner comb law and its lcm/CRT compound resonances, and it realizes the
concurrent Opus `HYP-7100` proposal that resolution into parallel/torsion
classes is the shared LRC grammar.  A scalar class energy loses this owner-set
sidecar; the remaining proof must retain it.

Reserved verifier/output:
`04-computation/lrc14_owner_packet_parallel_class_current_codex_S19.py` and
`05-knowledge/results/lrc14_owner_packet_parallel_class_current_codex_S19.out`.
The remaining task is a universal signed bound for the noncommon gcd-lattice
part of (4); a scalar crossing energy cannot supply it by itself.
