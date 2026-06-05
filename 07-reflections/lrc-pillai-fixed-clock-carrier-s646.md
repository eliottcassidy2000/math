# LRC Pillai Fixed Clock Carrier S646

The perfect-number prompt paid off in a very specific way.  Perfect numbers are
not just solutions to a scalar equation; they are fixed points of a divisor
carrier.  For LRC the analogous divisor carrier is not `sigma`, but the
gcd-shell mass of the pair-sum clock

```text
C = 2n - 1.
```

For odd `C`, define

```text
A(C) = sum_{1 <= a <= (C-1)/2} gcd(a,C).
```

This is half of Pillai's gcd-sum after removing the zero residue:

```text
A(C) = (P(C)-C)/2,
P(C) = sum_{r=1}^C gcd(r,C).
```

The fixed equation is therefore

```text
A(C) = C
```

or `P(C)=3C`.  That is the LRC shell-clock analogue of
`sigma(N)=2N`.

The delightful part is that the odd fixed clocks are exactly

```text
C = 15 and C = 27.
```

So the corresponding LRC rows are

```text
n = 8 and n = 14.
```

This is not a numerology echo.  It is a short divisor-carrier classification.
The local factor at `p^a` is `1+a*(1-1/p)`.  Solving the product equation
`P(C)/C=3` leaves only `3*5` and `3^3` in the odd world.

That makes `n=14` feel different.  The modulus `27` is not merely composite,
not merely `3^3`, and not merely the THM-401 pair-sum clock.  It is a fixed
clock of the gcd-shell divisor carrier:

```text
9 unit shells     contribute 9
3 gcd-3 shells    contribute 9
1 gcd-9 shell     contributes 9
total             contributes 27
```

The AP row is shell-fixed.  `Vstar` is the interesting one:

```text
AP:    1,2,3,4,5,6,7,8,9,10,11,12,13
Vstar: 1,2,3,4,5,6,7,8,9,10,11,13,24
```

Modulo `27`, the replacement moves shell `12` to shell `3`.  So shell support
detects a defect.  But both shells are gcd-`3`.  The divisor carrier sees no
defect at all:

```text
AP gcd counts     {1:9, 3:3, 9:1}
Vstar gcd counts  {1:9, 3:3, 9:1}
```

That is progress because it tells us what the known non-AP tight row is.  It is
not a random leaked shell; it is a gcd-carrier fixed row.  In the AP
single-swap scan through `n=14`, every tight non-AP row preserves this weighted
gcd-shell mass.  The tight rows are:

```text
n=5:   2 -> 7
n=6:   2 -> 9
n=8:   6 -> 12
n=14: 12 -> 24
```

The first two are shell-reflection fixed.  The latter two are gcd-stratum fixed
but not shell-fixed.  That splits the sporadic floor rows into two mechanisms:

```text
same shell        -> reflection sporadic
same gcd stratum  -> divisor-carrier sporadic
```

The next proof target is now sharper than "classify all weird rows."  It is:

```text
mass-changing row -> loose
mass-fixed row    -> AP/Vstar or killed by pair-pinch/carry-owner labels
```

This is exactly the kind of side-channel jackknife the repo keeps rediscovering.
The scalar `M=1/n` is the public quotient.  The gcd carrier, pair-pinch owner,
and carry cocycle are the retained channels.  Perfect numbers taught us to
look for a fixed divisor carrier, and `C=27` obliged.

The caution is also clear.  Fixed gcd mass is necessary in the small
single-swap scan, not sufficient.  There are many same-mass loose rows.  So the
carrier is a filter, not the proof.  The proof would be the no-leak lemma that
same carrier mass plus the THM-401/THM-369 pair-pinch ledger leaves only AP and
Vstar on the floor.
