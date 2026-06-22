# LRC additive-energy majorization

S103 attacks HYP-2885 by separating the true Fejer carrier from false monotone
orders.

The exact computation gives the clean correction:

- scalar additive energy is not monotone for `p0` or `L_y`;
- pairwise difference-profile majorization is not monotone for `p0`;
- one-step compression toward an interval is not monotone for `p0`;
- nevertheless the AP difference profile majorizes every tested row and the AP
  is the exact `p0` and `L_y` maximizer in the same banks.

The theorem shape should therefore be AP-facing:

```text
AP difference-profile majorization
  + signed sector/Fourier remainder control
  => L_y(E) <= L_y(AP_k).
```

It should not be stated as `L_y <= G(A(E))` for a scalar monotone `G`.

After rebasing over concurrent S103/S39 work, this note is attached to
HYP-2889/T1003.  The octahedral-current result owns HYP-2887/T1002, and
mac-mini's HYP-+2888 adds the boundary endpoint: the strict-threshold AP
extremizer covers exactly to measure `1`, with a rational boundary witness.
So the Fejer branch should prove no non-AP row over-covers past that boundary
value; it should not try to produce a positive-measure safe floor at the strict
threshold.

## Exact evidence

Script/output:

- `04-computation/lrc_additive_energy_majorization_codex_s103.py`
- `05-knowledge/results/lrc_additive_energy_majorization_codex_s103.out`

Banks:

```text
k=8:  1716 rows, E=(0 plus 7 from 1..13)
k=9:  1287 rows, E=(0 plus 8 from 1..13)
k=10: 220 rows,  E=(0 plus 9 from 1..12)
```

All three have:

```text
AP diff-profile majorization failures = 0
rows with p0 > p0(AP)                 = 0
rows with L_y > L_y(AP)               = 0
```

The false routes fail sharply:

```text
scalar A -> p0 monotone inversions: 3137 total
profile-majorization -> p0 monotone violations: 556752 total
profile-up one-step compression with p0 down: 1177 total
```

The worst scalar inversion is not small noise.  At k=9:

```text
lower-A  (0,1,2,3,4,5,6,7,12) A=401 p0=1093/2940
higher-A (0,2,3,4,5,7,8,9,10) A=405 p0=689/5880
gap = 499/1960
```

This rules out a clean scalar theorem of the form "higher additive energy gives
higher cover probability."

## Proof implications

The likely rigorous sequence is:

1. Prove the AP difference-profile majorization theorem.  This is a
   Pollard/Karamata-style finite additive-combinatorics lemma:

```text
d_AP^down = (k-1,k-2,...,1) majorizes d_E^down
```

where `d_E(h)=#{ unordered pairs {a,b}: |a-b|=h }`.

2. Decompose the THM-534 `L_y` functional into an AP-facing Fejer part plus a
   labelled signed sector remainder.

3. Prove the remainder cannot reverse the AP advantage.  Low relation-depth
   packets should go to finite AP/Freiman atlases; high relation-depth packets
   should go to the HYP-2636/HYP-2884 cancellation machinery.

This connects HYP-2885 to HYP-2886 cleanly: the cap branch gets the AP-facing
Fejer proof, while the exact-period packet branch handles finite witness atoms
when denominator/residue labels are the sharper carrier.
