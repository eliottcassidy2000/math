# THM-2928 septimal-peel audit and the located `d=6` closure

**Scope:** independent audit, scratch only.  LRC(14) remains open.

## Audit verdict

The following steps are sound.

1. If `7` does not divide `D`, every `D`-fibre contains exactly `M/7`
   cells from each one-body danger word.  Indeed, for `P=L/f`, the orbit in
   the fibre has length `h=P/gcd(P,D)`.  Here `7|h`, `h|M=L/D`, and a cyclic
   half-open block of length `P/7=gcd(P,D) h/7` contains exactly `h/7`
   points of the orbit.  Repetition gives `M/7`.  Six bodies delete at most
   `6M/7`, proving `S_D=Z/DZ`.
2. From `d|7q`, coprimality after dividing by `gcd(d,q)` gives
   `d/gcd(d,q) in {1,7}`.  In the second case the seven lifts in a
   `q`-fibre are distinct and a strict interval of length `d/7=gcd(d,q)`
   contains at most one.  In the first case all seven lifts coincide.
3. The optional-bit marginal is exactly `r/7` by Fubini.  No independence is
   used.  The fractional-cover optimum is an upper relaxation of the actual
   joint law.  If a positive-demand cover exists, the relevant integer
   superlevel event is open and contains the compact aligned-safe set.  A
   killed case has fractional optimum at most `u_2<1`, so the event is
   automatically proper and the comparison is strict.  `N_c=0` must and does
   pass separately.
4. The exact-lcm coefficient follows by Mobius inversion on `E|D` after
   partitioning the divisor alphabet, relative to the fixed `q=D/7`, into
   transverse, positive-remainder vertical, and zero-remainder vertical
   symbols.  These alphabets are disjoint; fixing the positive multiset with
   lcm `ell` leaves exactly

   ```text
   sum_(ell|E|D) mu(D/E) Mult(U_D(E),c) Mult(Z_D(E),z).
   ```

   The C++ state engine agrees with the Python Mobius count for every queried
   `D,c`, and the Python referee literally enumerates all 52,925 exact-lcm
   shapes in the 18 support-hard cases through `D<=300`.

The one correction is at strict endpoints.  The displayed identity

```text
Y_d(u)=w(a+X_d(u))
```

is an almost-everywhere identity.  If `r=0`, both ends of the open interval
can be lattice points, and the count falls to `w(a-1)` at equality phases.
For example `d=7`, numerator `1`, `u=1/2` selects no residue although
`a=1`.  The implemented `wa` is nevertheless a pointwise upper envelope,
so granting it is a valid necessary relaxation.  Canonical prose should say
this explicitly and call the zero-remainder term an a.e. floor / pointwise
upper grant.

## Uniform located closure of the lone `d=6` spike

Let the reduced numerator of the vertical denominator six be
`s in {1,5}`.  The high optional bit is

```text
V_s={u: ||s u||<3/7}.
```

Indeed multiplication by `s` permutes the six residue classes.  Its compact
low set is

```text
C_s=T\V_s={u: ||s u-1/2||<=1/14},       mu(C_s)=1/7.
```

For one aligned multiplier `a`, put `g=gcd(a,s)`, `A=a/g`, `S=s/g`.
Scaling by `g` and Parseval give

```text
mu(C_s cap D_a)
 =1/49 + 2 sum_(t>=1)
   (-1)^(At) sin(pi S t/7) sin(pi A t/7)/(pi^2 A S t^2)
 =1/49 + [B_2({(S+6A)/14})-B_2({(S+8A)/14})]/(AS),
```

where `B_2(x)=x^2-x+1/6`.  Write the Bernoulli difference as `h/49`.
For `S=1`, the exact table by `A mod 14` is

```text
h=(0,-1,5,-3,3,-5,1,0,-1,5,-3,3,-5,1).
```

For `S=5`, it is

```text
h=(0,-5,4,-8,8,-4,5,0,-5,4,-8,8,-4,5).
```

Consequently

```text
mu(C_s cap D_a)<=1/14,
```

with equality only at `a=2s`; for every `a!=2s` it is at most `1/28`.
For two distinct aligned multipliers, at most one can equal `2s`, hence

```text
mu(C_s cap (D_a union D_b))<=1/14+1/28=3/28<1/7=mu(C_s).
```

Thus `C_s` is not covered by the aligned danger union, equivalently the
aligned safe set is not contained in `V_s`.  Whenever `c=4`, the sole
vertical denominator is `6`, and `N_4>0`, a cover would force precisely
this impossible containment.  This closes that subfamily uniformly.

The exact residual pass removes `8,998,004` global denominator shapes and
`248,154,771` body/denominator occurrences.  The current floor/exception
ledger becomes

```text
26,899,164,786 shapes; 200,141,092,521 occurrences;
4,354 rows; 2,966 bodies; 147 resolving moduli.
```

Its `c=4` column becomes

```text
3,199,399,598 shapes; 39,514,128,418 occurrences;
3,185 rows; 2,812 bodies; 130 resolving moduli.
```

The former minimum `D=168`, `F=(1,2,3,4,6,12)` survivor is removed.  The
cheapest next computation is to add the exact Bernoulli overlap filter to the
canonical GF and print the new minimum survivor overall and by `c`; then scan
the other one-spike denominators by their translated low-comb overlap tables.
