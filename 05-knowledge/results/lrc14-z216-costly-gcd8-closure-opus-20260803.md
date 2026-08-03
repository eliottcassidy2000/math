# LRC(14) projected `z1=216` costly gcd-eight closure

**Status: FINITE-EXACT + INDEPENDENTLY AUDITED.**  This closes two rows in a
necessary projected `k=3` atlas.  It is not a physical-cover classification,
does not close the remaining gcd `24/36/72` wall families, and does not prove
LRC(14).

## Exact statement

In the corrected THM-2941 projected atlas at `z1=216`, the two rows omitted
from THM-3264 solely by its intrinsic work cutoff are empty:

```text
row 238: E=(1,8,10,11,13,14), L=560560, r=34,
         high_floor=55207, wall=True;
row 370: E=(2,8,10,11,13,14), L=560560, r=34,
         high_floor=55207, wall=True.
```

Each has invoice `L*r=19,059,040`, above THM-3264's declared `2,000,000`
cutoff.  The exact THM-3139 screen and complete THM-3078 terminal mechanism
give

| row | states | crude | exact status | residual | terminal split | min slack |
|---:|---:|---:|---:|---:|---:|---:|
| 238 | 1035 | 237 | 683 | 115 | `131` coarse + `6` exact | 1 |
| 370 | 354 | 213 | 122 | 19 | `16` coarse + `3` exact | 1 |

Every residual mask has `D/first_d=8`.  The duplicate-permitting two-high
gaps are respectively

```text
66973458766897 / 29716657283313000 > 0,
3668224417997 / 1297249202749500 > 0.
```

The wall gate requires at least one high slot, the positive gaps forbid two,
and the complete enlarged one-high banks violate the necessary support count
in every case.  The direction used is only

```text
enlarged necessary system empty  =>  original projected row empty.
```

No converse from quotient feasibility is used.

## Independent audit

The primary probe imports the promoted screen/terminal engine.  A separate
wrapper-free audit therefore:

- reparsed all `6060` atlas rows and the `480` rows at `z1=216`;
- recovered target addresses `238,370` and the layer digest
  `53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649`;
- checked all `683+122=805` Farkas contradictions with a separately spelled
  exact rational verifier and rejected a negative-multiplier mutation;
- independently rebuilt the signed low/high ray tables, making `524599` and
  `70839` primitive ray checks;
- rebuilt every complete terminal cell on the full `Z/LZ` grid before
  comparing with the inherited range representation; and
- reproduced both screen and terminal semantic hashes, with no max-gap
  fallback and no failure.

The audit also caught and repaired one non-load-bearing metadata error in the
first scratch transcript: `55207` is `high_floor`, not the component count;
the correct component count is `r=34`.  The invoice and all downstream tuple
indices were already correct.

The logical audit checked the wall consequence against THM-2941, translated
`kappa(d)=ceil(d/7)` against THM-2984/MISTAKE-334, nonpositive ray retention
against MISTAKE-338, and the noncanonical Farkas-witness digest boundary
against MISTAKE-331/333.

## Composed consequence

The targets are disjoint from THM-3264's 17 low-cost gcd-eight rows,
THM-3270's 33 order rows, and THM-3281's 47 natural-family rows.  Hence the
current exact ledger changes by precisely two:

```text
projected ledger: 373186 -> 373184,
z1=216 wall rows:     382 -> 380,
cap:                         216 (unchanged).
```

This completes the gcd-eight stratum at `z1=216`.  The next economical probe
is a complete intrinsic ruler family among the remaining gcd `24/36/72`
rows, ranked by total `sum L*r`, rather than another externally chosen cost
slice.

## Reproduction and artifacts

Primary inherited-engine replay:

```text
python .scratch/lrc14_z216_costly_gcd8_probe_lrc_probe_agent.py --processes 2
python -O .scratch/lrc14_z216_costly_gcd8_probe_lrc_probe_agent.py --processes 2
```

Independent reconstruction:

```text
python .scratch/lrc14_z216_gcd8_independent_audit_gcd8_independent_audit.py --processes 2
```

Artifacts and current SHA-256 hashes:

```text
primary script  1eacc80c17baa70f9b11d99bade90b3cd470e9c11aee0b15e8426e7e4ceb0234
primary output  936b5871ac8cc616b1d52a1bd1771c2cb796e8baea593b8021533c6544fa6956
audit script    e6d6c8428c6216e27d3e0581231884bbfa74dc9bf5c1311314d2c37b7ee6e35e
audit output    1432ffaca0025137da91ed392313c0a9da48d6d657c7b6faefcc157660c5e92a
```

The primary wrapper remains useful execution redundancy; the independent
audit is the load-bearing promotion check.
