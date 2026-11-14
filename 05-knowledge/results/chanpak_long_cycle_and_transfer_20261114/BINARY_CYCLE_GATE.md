# Every cycle count has a logarithmic physical encoder; ordinary counts still need a switch

For distinguishability, this is about directed Hamiltonian **cycles**
modulo cyclic rotation, not the ordinary count of Hamiltonian paths.
It gives an exact way to encode every positive integer in strong
tournaments of logarithmic order, and to read that integer back from
a short list of ordinary counts. It does not turn the signed readout
into one ordinary path count.

## 1. Choose the binary digits of the cycle count directly

Given \(K>0\), write it in \(n\) binary digits \(d_1\cdots d_n\);
leading zeroes are allowed and \(n\) may be taken odd. The vertices
are \(x_1,\ldots,x_n,a,b,u,v\). Use all of these arcs:

- \(x_i\to x_j\) for \(i<j\); all \(x_i\) beat both \(a,b\);
  the last pair is ordered \(b\to a\).
- \(u\) beats all \(x_i\), \(a\), and \(v\); \(b\to u\).
- \(a,b\to v\), and \(v\to x_i\) exactly when \(d_i=1\).

A Hamiltonian cycle rooted at \(u\) has the form
\[
                            u\,P\,v\,Q\,u .
\]
Here \(P,Q\) split the spine vertices \(x_i,a,b\). The only old
vertex that can close the cycle into \(u\) is \(b\), so \(Q\) ends
there. In the spine alone \(a\) is a sink, so \(a\) cannot be in
\(Q\) and must end \(P\). Conversely every assignment of the \(x_i\)
to the two runs respects their increasing order; the only remaining
rule is that \(Q\) must begin at an \(x_i\) with \(d_i=1\)
(it cannot begin at \(b\) since \(b\to v\)).
For a fixed earliest \(x_i\) in \(Q\), every later \(x_j\) is
independently in \(P\) or \(Q\), and all earlier ones are in \(P\).
Therefore the literal cycle count is
\[
                  \operatorname{HC}(T)=\sum_{i=1}^{n}d_i2^{n-i}=K.
\]
The explicit cycles also prove strongness when \(K>0\).
Deleting the two hubs \(u,v\) leaves a transitive tournament;
the single reversed final pair in the display is just a
different order of the last two vertices. Thus two hubs
suffice for the complete cycle spectrum. A single hub over a
transitive remainder gives at most one directed Hamiltonian
cycle, since starting at the hub fixes the entire remaining
vertex order. This is a sharp boundary for cycle counts by that
description. At each **fixed** order \(n+4\) this fills
\(1,\ldots,2^n-1\), and the endpoint \(2^n\) is obtained
by making all the digits 1 and changing \(b\to v\) to
\(v\to b\): the only extra choice is \(Q=(b)\).
Conversely any tournament with a transitive remainder after
deleting two hubs has at most \(2^{N-2}\) cycles at order
\(N\), since each is determined by the assignment of those
\(N-2\) vertices to its two ascending runs.
Taking an extra leading zero if needed realizes
\(K\) in **odd** order at most \(\lfloor\log_2 K\rfloor+6\).

## 2. The ordinary-count readout is exact but signed

For a tournament \(T\) on odd \(N\ge3\), define the polynomial
\(F_m(T)\) as in the [color-channel formula](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/WRAP_CHANNELS.md).
A set of disjoint directed odd cycles \(\mathcal C\) contributes
\(2^{|\mathcal C|}m^{N-|V(\mathcal C)|+|\mathcal C|}\). The term
of degree one must use a single spanning odd cycle, so its
coefficient is \(2\operatorname{HC}(T)\). The remaining powers
have the same parity; write \(F_m(T)/m=P(m^2)\) with
\(\deg P=(N-1)/2\).

Let \(L=(N+1)/2\). Lagrange interpolation at \(1^2,\ldots,L^2\)
gives the explicit signed readout
\[
 2\operatorname{HC}(T)=
 \sum_{i=1}^{L}
 \frac{2(-1)^{i-1}(L!)^2}{(L-i)!(L+i)!}\,\frac{F_i(T)}{i}.
\]
For instance the coefficient is found by taking the product
\(\prod_{j\ne i}(-j^2)/(i^2-j^2)\); there is no numerical
approximation.

These color counts can in turn all be recovered by **ordinary**
Hamiltonian path counts of explicit tournaments of order \(O(N)\).
Let \(W\) be the single outer wrap. Use the proved transfer
\[
   F_m(W(T))=m^2F_m(T)+m(F_{m+1}(T)-F_{m-1}(T)).
\]
In \(H(W^j(T))\) the highest-index original channel is \(F_{j+1}(T)\)
with coefficient \(j!\), since every one of the \(j\) steps must
increase the channel to reach it. Every other term uses only
indices at most \(j\). Hence the \(L\) ordinary counts
\(H(T),H(W(T)),\ldots,H(W^{L-1}(T))\) recover all \(F_i\) by exact
triangular inversion. The largest graph has \(N+2(L-1)=2N-1\)
vertices, logarithmic in \(K\) for the binary encoder.
What has been transferred is the **readout**, not a physically
signed gluing of the ordinary counts.

## 3. Why a fixed finite linear readout cannot become the missing gate

Even on strong tournaments of odd order, no fixed affine linear
combination of finitely many channels \(F_m\) can equal
\(2\operatorname{HC}\). Test on the nearly-transitive strong
tournament \(U_N\), obtained by reversing only the first-last arc
of the transitive order, for any odd \(N\ge3\). It has
\(\operatorname{HC}=1\). All directed cycles use the same
backward arc; summing its subsets of intermediate vertices
in the odd-cycle expansion gives, for each integer \(m\ge1\),
\[
 F_m(U_N)=m^N+m\bigl((m+1)^{N-2}-(m-1)^{N-2}\bigr).
\]
Its dominant exponential base in \(N\) is \(m+1\), strictly
different for each \(m\). In a finite fixed linear combination,
the largest positive-index channel present cannot be cancelled.
So all coefficients at \(m\ge1\) would vanish, leaving a constant.
The binary construction with \(K=2\) supplies a strong odd-order
counterexample to that constant. A fixed finite affine combination
of ordinary counts after bounded repeated wraps also fails, by
their expansions in finitely many fixed channels.

This isolates a specific missing proof step: to send this universal
cycle spectrum to ordinary counts, a transformation must either
perform a genuinely nonlinear/size-dependent removal of the
nuisance channels, or find an entirely different way to make the
desired binary choice contribute at absolute height. Fixed signed
color sampling alone cannot be that transformation.

Run `python binary_cycle.py` to replay the literal binary encoder,
odd orders, color coefficient, and inverse wrap readout.
