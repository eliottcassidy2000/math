# Vieta Perfect Aliquot Carriers S644

The new useful object is not "perfect numbers explain LRC."  The useful object
is a carrier stack:

```text
A = C(n,2)
8A+1 = (2n-1)^2
s(A) = sigma(A)-A.
```

The first line is the pair-count scalar.  The second line is the Vieta square
root that recovers the LRC shell clock.  The third line is the divisor-sum side
channel.

This makes the `n=14` row unexpectedly crisp:

```text
C(14,2) = 91 = 7*13
s(91) = 1+7+13 = 21
sqrt(8*91+1) = 27.
```

So the same triangular carrier gives both the hard LRC clock `27` and a
visible `21` shadow.  That is a better statement than "14 and 21 are related."
It says exactly which observer sees 21: the aliquot observer on the pair-count
row.

Perfect numbers fit as fixed controls.  The even perfect numbers are triangular:

```text
6=C(4,2)
28=C(8,2)
496=C(32,2)
8128=C(128,2)
```

and each is a fixed point of `s`.  Thus perfect numbers are where the aliquot
side channel returns the pair count itself.  The `n=14` row is different: it is
deficient, but its shadow is the durable problem-scalar `21`.

The user's `6` seam is also exact.  If `2p=3q` for primes `p,q`, then `p=3`,
`q=2`, and the common value is `6`.  That same `6` is the first perfect number,
the triangular row `C(4,2)`, the shell-root row `sqrt(8*6+1)=7`, and the
unique distinct positive product-sum resonance from THM-361:

```text
6 = 2*3 = 1+2+3 = 1*2*3.
```

This makes `6` a base crossing of additive, multiplicative, triangular, and
divisor-sum observers.

The most promising extension is the semiprime family:

```text
n=2p, q=2p-1 prime
C(n,2)=p*q
s(C(n,2))=1+p+q=3p.
```

The `p=7` instance is exactly:

```text
n=14
C(n,2)=91
s=21
```

and it matches S642's diagonal packet `(2p,3p)=(14,21)`.  Incoming
HYP-2219/S643 sharpens that diagonal into a companion-graph duplicate branch:
the unordered pair `{p,p}` has only one odd companion instead of two.  So S642,
HYP-2219/S643, and S644 now meet in a shared prime-pair/aliquot carrier:

```text
Goldbach/Lemoine diagonal: p -> (2p,3p)
Bridge-graph duplicate:    (E,O)=(2p,3p) is a branch locus
Triangular aliquot family: C(2p,2) -> 3p, when 2p-1 is prime.
```

The warning is important: none of this proves LRC `n=14`.  The proof still
lives in the `C=27` lift/CRT/owner channel.  But this gives a new cheap
arithmetic observer that may help classify which rows are carrier-crossings and
which are merely scalar echoes.
