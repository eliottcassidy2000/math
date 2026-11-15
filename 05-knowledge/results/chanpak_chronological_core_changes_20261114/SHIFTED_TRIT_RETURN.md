# A physical shifted-trit return on one entire dyadic arm

The [induced-extension theorem](../chanpak_collatz_exact_addresses_20261114/TRIANGLE_TARIFF.md) prohibits getting H→3H+2 by keeping an arbitrary H>1 source and a disjoint directed triangle as induced subtournaments in ordinal order. The [one-outlet transfer](ONE_OUTLET.md) bypasses that restriction for one exact infinite family: it replaces the exterior host rather than retaining that triangle bank.

The lower arm of the consecutive triple is h=9·2^m−5 for every m≥2. Its witness uses the exterior rows (20,1,2,7,14) with outlet 3 and a transitive replacement of size m−1. Change the exterior tournament to rows (34,0,3,7,15,30) with outlet 3, and use a transitive replacement of size m. The exterior path-count and parity profile is (A,E,O)=(9,11,16), in the notation of ONE_OUTLET. The exact new count is

    A−2E+(E+O)2^m = 27·2^m−13 = 3h+2.

Both templates are physical tournaments. No new triangle bank is retained, so this does not contradict the induced-extension barrier. For example the construction sends 67 to 203, and the same switch works at every width in the arm, with the width known exactly from h. The output is 5 mod 9 while the input arm is 4 mod 9, so these two templates do not give closure under iteration or an inverse for arbitrary targets.

Replay with Python from this folder:

```python
from replay_one_outlet import gate, attach, transitive, table
R = (34,0,3,7,15,30)
assert gate(R,3) == (9,11,16)
for m in range(2,7):
    assert sum(table(attach(transitive(m),R,3))[-1]) == 27*2**m-13
```
