# HYP-7083 — the finite-`t` owner-doubling packet

**Status:** OPEN after an exact reduction and finite rational census
(codex-2026-07-16-S18).

Let `E={0,e_1,...,e_5}`, let `t>max(E)`, put `F=E union {t}`, and write `p_j(F)`
for the mass on which `F` misses exactly `j` sectors.  The `THM-727` finite correction
has the exact increment form

```text
Delta_F(w)=p_0(F union {w})-p_0(F)-p_1(F)/7.
```

This immediately removes the first owner tooth:

```text
Delta_F(t)=-p_1(F)/7 <= 0,
```

because the duplicated `t` phase cannot fill the unique sector already missed by `F`.
Thus `a=1` is never an upper-tail propagation obstruction.

For `w=2t`, set

```text
q_t(x)=floor(14{tx}),  c_t(x)=floor(q_t(x)/2).
```

Then `sec(tx)=c_t(x)` and `sec(2tx)=q_t(x) mod 7`.  If `s_F(x)` is the unique
miss label on the one-miss set, the remaining correction is exactly

```text
Delta_F(2t)
 = integral_[p1(F)] (1_[s_F(x)=q_t(x) mod 7]-1/7) dx.          (1)
```

Equivalently, with `epsilon=(-1)^q`,

```text
1_[q mod 7=s]
 = (1_[s=2c]+1_[s=2c+1])/2
   + epsilon(1_[s=2c]-1_[s=2c+1])/2.                          (2)
```

Equation (2) isolates the finite wall remainder as one signed half-sector parity
packet.  The correct vertices are therefore `(miss label, t-sector, half-sector
parity)`, not runners, arcs, or the six limiting residues.  This quotient preserves the
finite LRC correction and destroys slow-wall chronology and the additive relation
lattice; those sidecars must be restored in a universal proof.

The exact verifier scans all `96,600` rows with primitive five-speed slow core through
diameter `10`, `D<t<=4D`, and `1<=a<=14`.  Every row is below `0.097` in absolute value.
The largest positive correction is

```text
2173/27440 = 0.07919...
```

at `E=(0,3,4,5,6,7)`, `t=8`, `a=2`.  All `6,900` `a=1` rows verify the exact
negative identity above.  This makes the owner-doubling parity packet, rather than the
full unspecialized `O_E(1/t)` term, the next finite-`t` target.

Tournament Analysis on multipliers is only telemetry: scalar positive-risk and
absolute-risk gauges are transitive with singleton SCCs and unique tie Hamiltonian
paths.  Switching gauges flips edges, but neither ordering retains the proof-bearing
miss/parity incidence.

Verification:
`04-computation/lrc14_finite_t_owner_packet_codex_S18.py` ->
`05-knowledge/results/lrc14_finite_t_owner_packet_codex_S18.out`.
