# The exact factor bank and all-modulus clock fit in one sparse grammar

In either selective quotient allow both modules to be arbitrary ordinal
sums of singletons, arbitrary strong one-apex cores, and the seven
[dyadic factors](DYADIC_INVERSE.md). This includes the fixed-source
prime-cut reservoirs and all modules used in the
[modular aperture](RESIDUE_APERTURE.md): the extra cyclic wrap there
is itself one apex over a chain, with cut 10...0. Hence this single
enlarged grammar still realizes every odd-compatible residue modulo
every modulus.

It nonetheless has only O(X^(1999/2000)) absolute outputs through X.
Let L(s) be the sum of Q(2)^(-s) over this alphabet (counting duplicate
letter types if needed). The [one-apex bound](ONE_APEX_CAPACITY.md)
and the seven dyadic costs give
```
L(1) <= 1189691/1242864 + Σ_(m in [5,8,12,12,13,14,14]) 2^(-m)
     = 315979709/318173184 < 199/200.
```
Likewise L(3/4)<8, by the inherited bound <7 and adding a sum <1.
At s=1999/2000, Hölder gives
L(s) < (199/200)^(499/500) · 8^(1/500) < 1, since
8·199^499 < 200^499. These are exact integer comparisons.

Ordinal multiplication of Q(2) makes the total weight of all finite
sequences at most 1/(1−L(s)); Q(3)>=Q(2) makes the same estimate
valid for the source side. The output bound
N>=2 Q_U(3) Q_V(2) then bounds the number of admissible input pairs
and hence the distinct absolute values by O(X^s).

Thus combining the physical factor decoder with unrestricted serial
reuse of the entire one-apex clock, even with all the verified new
dyadic cores, cannot cover all primes or all 7p targets. One needs
additional physical capacity (e.g. genuinely different strong cores
or changing the selective quotient), not a larger CRT modulus.
