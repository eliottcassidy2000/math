# A dyadic factor gives a narrower single flip centered at the 61st-exponent perfect number

**PROVED.** Let \(E=2^{60}(2^{61}-1)\). A pair of tournaments on **125 vertices** differing in precisely one arc has Hamilton path counts
\[
 2658455985209121452320167624203544159,\qquad
 2658455997930542036989217607704140193 .
\]
Their midpoint is exactly \(E\) and their halfgap is \(6360710292334524991750298017\). The halfgap improves the [uniform strong pair on 122 vertices](CENTER.md) by a factor greater than 78 million. The 125-vertex tournaments are ordered joins of two strong components. The halfgap-one straddle stays open.

## The physical flip and factor

For a transitive spine with one apex, a word records 1 when the apex beats a spine vertex. Set \(u=102580081,\;w=5619616816\). Take the two contact words
\[
 \mathrm{bin}(u)\;b\;\mathrm{reverse\mbox{-}complement}(\mathrm{bin}(w)),
 \qquad b\in\{0,1\},
\]
using ordinary binary representations. Their two 62-vertex tournaments differ only at the central contact. The two flanks, in physical order, are

```
110000111010011111101110001
111100111100111011010000100011010
```

Their respective counts are \(1152921501848328609\) and \(1152921507365365343\), exactly \(2^{60}\pm2758518367\). An independent direct count from either word \(s\) uses
\[
 H(s)=1+\sum_{i<j:\,s_i=1,\,s_j=0} (2^{j-i-1}+[j=i+1]).
\]
Adjoin after each tournament the same 63-vertex component with contact word \(1\,0^{29}\,1\,0\,1\,0^{29}\) and orient all joining arcs forward. This component has count \(2^{61}-1\) by the same direct formula, or by the symmetric three-block rail with \(a=29,b=1\). Every spanning path in the join is the concatenation of one path in each component. Multiplying both endpoint counts thus preserves one physical flip and gives the stated endpoints. All three words start with 1 and end with 0, so their apex and entire spine form a directed cycle. [Symmetric Mersenne-offset rail and direct apex formula](../chanpak_heegner_absolute_targets_20261123/MERSENNE-RAIL.md)

## Optimality for the dyadic contact step

**FINITE-EXACT.** Among all single-apex contact flips with a transitive spine and midpoint \(2^{60}\), the smallest possible halfgap is \(2758518367\). The first-stage optimizer is unique in the two integer addresses up to exchanging them. This quantifier permits arbitrary spine length. It does not cover flipping a spine edge, other tournament architectures, or centered flips whose midpoint takes a different factorization.

Use the zero-based reader \(f(0)=0,\;f(2q+1)=f(q),\;f(2q)=f(q)+q+[q\text{ odd}]\). With \(c(x)=\lceil x/2\rceil\), the inherited exact contact formula gives
\[
 M(u,w)=2uw+f(u)+f(w)+1+c(u)+c(w),\qquad
 d=|c(w)-c(u)|.
\]
Here \(u,w\) are arbitrary nonnegative addresses. Symmetry permits \(u\le w\). The recurrence gives \(0\le f(x)\le x\): at an odd argument the reader drops to its half, while at twice an odd \(q\) the bound \(f(q)\le(q-1)/2\) absorbs the extra unit. Every \(f(x)\) is even. [One-site contact identity](../chanpak_apex_flip_local_global_20261122/README.md)

Let \(D=2758518367\). If \(u<50,000,000\) and \(d\le D\), then \(w\le u+2D+1\le5,567,036,735\). Even substitution of these upper bounds into \(2uw+u+w+1+c(u)+c(w)\) gives less than \(2^{60}\). At the other extreme \(2u^2\le 2uw<M\) forces \(u\le759,250,125\). In the remaining interval set \(R=2^{60}-f(u)-1-c(u)\). Every candidate has
\[
 \left\lceil\frac{2R-1}{4u+3}\right\rceil \le w\le
 \left\lfloor\frac{2R}{4u+1}\right\rfloor,
 \qquad c(w)\not\equiv c(u)\pmod2.
\]
The bounds follow by applying \(0\le f(w)\le w\) to \(R=2uw+c(w)+f(w)\). Testing every such integer pair leaves exactly \((102580081,5619616816)\), with \(f(u)=45284624\) and \(f(w)=5351979710\).

The [replay](dyadic_lift.cpp) defaults to the three direct, quadratic-time word counts. Its optional exhaustive mode computes the recurrence in 12-bit blocks. Appending a low block \(s\) to a prefix \(h\) adds \(h\) times the reversed bitwise complement of \(s\), then \(f(s)\), then one if \(h\) is odd and the high bit of \(s\) is zero. This follows by summing the recurrence at all zero bits of the appended block. All scanned addresses fit in three blocks. Compile and run in this directory:

```sh
c++ -O3 -std=c++17 dyadic_lift.cpp -o /tmp/chanpak-dyadic-lift
/tmp/chanpak-dyadic-lift
/tmp/chanpak-dyadic-lift 50000000 759250125
```
