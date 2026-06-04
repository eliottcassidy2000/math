# Pi/e Lonely Shadow Fallout S638

The new useful phrase is:

```text
an algebraic e+pi or e*pi would be lonely.
```

S635/S636 already made the trace/norm point: `S=e+pi` and `P=e*pi` cannot both
descend, and any two of `S,P,D=e-pi` reconstruct the hidden pair.  S638 adds the
transverse fallout.  If one elementary coordinate descended, many nearby
symmetric coordinates would be forced transcendental.

The cleanest surprise is the `log(pi)` gate.  If `L=log(pi)` were algebraic,
then `pi=e^L`.  Since `pi` is neither `1` nor `e`, the exponents `0,1,L` are
distinct, and Lindemann-Weierstrass would make `1,e,e^L` linearly independent
over algebraic coefficients.  So `e+pi` and `e-pi` could not be algebraic.
Hermite-Lindemann would also make `e*pi=e^(1+L)` transcendental.  Therefore:

```text
any algebraicity among e+pi, e*pi, e-pi would prove log(pi) transcendental.
```

That does not solve the original problem, but it moves the dependency graph.  A
proof that `e+pi` is algebraic would not be a small oddity; it would immediately
settle a separate log-commensurability question in the hard direction.

The power-sum fallout is the other handle.  Newton gives

```text
p_k = e^k + pi^k = f_k(S,P).
```

If `S` is algebraic and some `p_k` with `k>=2` is algebraic, then `P` becomes
algebraic over `Qbar`, and Vieta collapses `e,pi` into an algebraic quadratic.
Impossible.  If `P` is algebraic, the same argument forces every `p_k` with
`k>=1` to be transcendental.

This is the repo pattern again: do not ask only whether a scalar descends.  Ask
which transverse packet is forced to stay wild if that scalar descends.

The incoming H=21 closure is the perfect neighbor.  `H=21` is not excluded by
staring at the number `21`; it closes when the strong-component and `c3<=10`
side condition is retained.  Likewise, `e+pi` and `e*pi` are not understood by
raw decimal or rationality checks.  They become a carrier with:

- elementary coordinates `S,P`;
- branch coordinate `D`;
- log-commensurability coordinate `log(pi)`;
- transverse symmetric shadows `p_k`;
- theorem-machine labels: Vieta, Lindemann-Weierstrass, Schanuel.

For Tournament Analysis, I did not use numbers as vertices.  I considered
shadows, assumptions, implication arrows, theorem machines, branch sheets, and
repo carrier analogies; the selected vertices were proof routes.  This preserves
field-descent obstruction data and discards decimal near-miss noise.

The next finite analogue should be a scalar-fixed tournament fiber: hold `H`
constant and vary `c3`, SCC, beta, OCF packet, and rooted-perspective data.  The
pi/e result says that a descended scalar often makes the transverse packet more
important, not less.
