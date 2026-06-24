# LRC14 Toeplitz-PSD Fourier Dual

Hypothesis: HYP-2974

This pushes the Fourier idea behind HYP-2971 and complements HYP-2973's
count-dual route.  If the danger arcs cover, then

```text
f_S(t)=C_S(t)-1 >= 0.
```

Every nonnegative function has PSD Toeplitz moment matrices:

```text
T_d(S) = (fhat_S(i-j))_{0<=i,j<=d}.
```

For LRC14 the coefficients are explicit:

```text
fhat_S(0)=6/7,
fhat_S(k)=sum_{s|k} sin(pi*(k/s)/7)/(pi*(k/s)).
```

So a negative eigenvalue of `T_d(S)` is a Fourier dual certificate.  If `c` is
the eigenvector and `p(t)=sum c_j e^{2*pi*i*j*t}`, then

```text
int (C_S(t)-1)|p(t)|^2 dt < 0.
```

That forces `C_S=0` on positive measure, hence a strict lonely interval.

The scan over named hard rows through degree `160` found:

```text
AP, GW 12->24: no PSD failure through 160
all 13 positive hard rows tested: PSD failure
```

First failure degrees:

```text
drop(2,6)->add(17,42): 32
13->26:                 37
12->84 / 12->168:       51
12->48:                 53
6->69 / 6->75:          56
drop(4,6)->add(19,42):  57
12->26:                 59
10->20:                 70
12->36:                 101
P10+GW:                 160
```

This is not simply better than HYP-2971.  HYP-2971 often certifies by a lower
scalar multiplicity barrier.  The Toeplitz route buys localization: the
negative eigenvector is a trigonometric square that should point toward the
safe interval and can be reattached to endpoint-owner, twist-ladder, or K33
state-lift labels.

The theorem target:

```text
If every finite Toeplitz section of C_S-1 is PSD,
then S is AP/GW boundary or state-lift-owned.
```

Tournament Analysis uses Fourier/Toeplitz proof carriers as vertices, not
runners.  The transitive carrier order is:

```text
Toeplitz certificate
> eigenvector localization
> Fourier coefficient ledger
> multiplicity histogram
> endpoint packet
> raw row.
```
