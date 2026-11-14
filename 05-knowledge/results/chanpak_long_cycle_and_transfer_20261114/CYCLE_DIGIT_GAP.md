# A cycle digit is already a signed absolute path-count correction

Use the [binary cycle encoder](BINARY_CYCLE_GATE.md) with
the ordinary orientation \(b\to v\) and a declared digit
length \(n\). It realizes \(K=\sum_{j=1}^n d_j2^{n-j}\).
Let \(T\) be that tournament, \(Q=T-a\), and put
\[
 B(d)=\sum_{j=1}^n d_j\,2^{j-1}3^{n-j},
 \qquad K^\mathrm{rev}=\sum_{j=1}^n d_j\,2^{j-1}.
\]
There is an exact, whole-count identity
\[
             H(T)-H(Q)=2K+2B(d)+K^\mathrm{rev}-d_1 .
\]
This is a literal count gap, not a congruence or the number
of paths at one endpoint.

## Explain both signs by the physical root

In \(T\), the removed vertex \(a\) beats only \(v\).
It is a zero-width bridge added to \(Q\) at \(v\), as in the
[ordinary root transfer](https://github.com/eliottcassidy2000/math/blob/bd8f3030aee3e75af8b67695a76a857e0094a18a/05-knowledge/results/chanpak_literal_prefix_compiler_20261114/ARBITRARY_ROOTED_BRIDGE.md).
Hence \(H(T)-H(Q)=R_+(Q,v)+R_-(Q,v)\).
By the [signed endpoint identity](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/READY_CLASSIFICATION.md),
the second term is minus the ordinary number of full paths
ending at \(v\), which is easily counted here.

Delete \(v\) too: on \(u,x_1,\ldots,x_n,b\) the only back
arc is \(b\to u\). There is one path without it, ending at
\(b\). Other paths put an arbitrary subset of the \(x_i\)
before \(b\to u\) and the rest after \(u\). The empty final
subset ends at \(u\); for a nonempty final subset with maximum
\(j\) one has \(2^{j-1}\) choices, and it can join \(v\) iff
\(d_j=0\). Thus the full paths ending at \(v\) number
\(1+2^n-K^\mathrm{rev}\).

For \(R_+(Q,v)\), the singleton selected set \(\{v\}\)
leaves \(1+2^n\) full paths outside. Every other selected
root-start path must have a first nonempty ascending run
\(S\) of \(x_i\)'s with allowed earliest digit. Cases:

* If neither \(b\) nor \(u\) is selected, the complement
  has \(1+2^{n-|S|}\) paths, except when \(|S|=n\) and it
  has only the two-vertex path, contributing 1.
* If \(b\) is selected but \(u\) is not, the selected path
  ends at \(b\) and the complement is transitive.
* If both are selected, the path uses \(b\to u\); each
  remaining \(x_i\) can be in its later ascending run or
  the complement, giving \(2^{n-|S|}\) choices for this case.

Selecting \(u\) alone without \(b\) is impossible after \(v\).
The number of nonempty allowed \(S\) is \(K\).
The sum of \(2^{n-|S|}\) over those \(S\) is \(B(d)\):
for earliest \(j\), each earlier vertex contributes 2,
and each later vertex contributes \(1+2=3\).
The all-selected exceptional correction is \(d_1\).
Thus
\(R_+=1+2^n+2K+2B(d)-d_1\), which proves the displayed gap.

## Four ordinary counts recover any prescribed even digit pulse

Pad the same digit word with \(j\ge1\) leading zeros.
The cycle count stays \(K\); both \(B(d)\) and
\(K^\mathrm{rev}\) are multiplied by \(2^j\), and the first
digit is now zero. Thus the literal path-deletion gap \(G_j\)
at that width is
\[
 G_j=2K+2^j(2B(d)+K^\mathrm{rev}),\qquad
                          2G_j-G_{j+1}=2K .
\]
This compresses the retrieval of the binary cycle target on
this special family to **four ordinary path counts**, all on
graphs of logarithmic order and with explicitly marked vertices.
It does not contradict the no-fixed-color-readout theorem,
which concerns a fixed formula on *all* strong odd-order inputs.
It still requires physically implementing a signed difference
of gaps to become one desired ordinary output, rather than
reading off an integer from known outputs. The actual excess
channel is isolated here, not replaced by a claim of absorption.

Run cycle_digit_gap.py for a literal subset-DP replay of the
formula and leading-zero cancellation.

## The whole ordinary total, and three concrete terminal phases

There is also a direct formula for the absolute output.
Let \(D\) be the one-apex graph on \(v,x_1,\ldots,x_n\)
with the same digit contacts, write \(h_D=H(D)\), \(f_D=F_2(D)\),
and retain the original bits \(d_1,d_n\).
Consider \(Q=T-a\) as \(u\to D\) with a last vertex \(b\)
beating \(u,v\), i.e. two contacts at width zero.
For a general contact set at width zero the added count is
\(R_++R_--2E\), where \(E\) counts old even-order paths with
both endpoints contacted, weighted by the complement paths.
This follows from the same odd-parity argument as the
[positive-width two-contact law](https://github.com/eliottcassidy2000/math/blob/3c4ef16d160baa9aca5ccbb7bbb86e8d52b37905/05-knowledge/results/chanpak_literal_prefix_compiler_20261114/TWO_CONTACT_REPAIR.md):
without any follower, the chosen old path must be even and
its last endpoint cannot be another contact.
Here the \(u\)-contact has \(R_+=f_D,R_-=0\), and the
\(v\)-contact has \(R_+=1+K,R_-=-[d_n=0]\).
The disallowed end-to-end path must be from \(u\) to \(v\).
The intervening \(x_i\)'s form an even subset, either empty or
ending at a zero digit. Counting by the largest selected index
gives \(2E=2^n-K^\mathrm{rev}+d_1\). Thus
\[
\begin{aligned}
 H(Q)&=h_D+f_D+K+d_n-2^n+K^\mathrm{rev}-d_1,\\
 H(T)&=h_D+f_D+3K+2B(d)+2K^\mathrm{rev}+d_n-2^n-2d_1 .
\end{aligned}
\]
The first formula also handles \(n=1\).
For direct computation, let \(P(z)=\sum_{i<j}
d_i(1-d_j)z^{j-i-1}\) be the inversion polynomial of the
contact word. Then the elementary one-apex path formula and
its weighted odd-cycle expansion give
\[
 h_D=1+P(2)+P(0),\qquad
 f_D=2^{n+1}+2^{n-1}\bigl(P(3/2)+P(1/2)\bigr).
\]
In the second identity all cycles meet \(v\);
an inversion with \(m=j-i-1\) interior sites contributes
\(2^n\sum_{t\ \mathrm{even}}\binom mt2^{-t}\), which
explains both arguments of the polynomial.

After \(j\ge1\) leading zeros, the exact **ordinary** total
is therefore a single affine ray
\[
 H(T_j)=h_D+3K+d_n+
       2^j\bigl(f_D+2B(d)+2K^\mathrm{rev}-2^n\bigr).
\]
The intercept comes from the original digit word: its one-apex
count stays fixed under prepending leading zeros, while its
\(F_2\) doubles at each such step. This gives an exact inverse
in the padding width whenever the underlying cycle target
is fixed. Also any two successive ordinary counts recover
the intercept by \(2H(T_j)-H(T_{j+1})\); together with \(h_D\),
three ordinary counts recover \(3K\) after subtracting \(d_n\).

For the previous focal numbers, these particular binary sources are:
\[
\begin{array}{c|c|c}
 K&H(T_j)&7\mid H(T_j)\ \text{iff}\\\hline
 11&37+152\cdot2^j&j\equiv0\pmod3\\
 13&45+176\cdot2^j&j\equiv2\pmod3\\
 26&97+564\cdot2^j&j\equiv1\pmod3 .
\end{array}
\]
So all three terminal phases actually appear in the arithmetic
of this cycle construction. For example at \(K=13\) the same
cycle count accompanies the ordinary prime 397 (\(j=1\))
and \(749=7\cdot107\) (\(j=2\)).
This is an exact physical relation, not a claim that these
three sparse rays fill every cofactor.
