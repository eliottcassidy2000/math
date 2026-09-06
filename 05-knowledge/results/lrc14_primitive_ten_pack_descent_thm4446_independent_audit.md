# Independent audit of THM-4446

**Verdict: PASS.** The primitive ten-pack proposition and integer
dilation-ray corollary are theorem-safe with their stated dependencies.
No hypothesis, inequality, endpoint convention, or conclusion required
repair. LRC(14) remains open.

## Proof audit

For \(C=pC'\), the \(3p\) points \((y+j)/(3p)\) are the complete fibre of
the circle map \(x\mapsto3px\) above \(y\). They all preserve the body. A
tail \(t\) produces \(\gcd(t,3p)\) copies of an order
\(3p/\gcd(t,3p)\) grid. The referee rederived the open-arc bound
\(\lceil m/7\rceil\), including the exact \(7\mid m\) endpoint case.

For a ternary-unit tail, the only cases are:

    p|t, p!=3:  order 3,  at most p labels killed;
    p does not divide t: order 3p, at most ceil(3p/7) killed.

Primitivity excludes three \(p\)-divisible tails. The independent word
census verifies \(kp+(3-k)\lceil3p/7\rceil<3p\) for every admissible
\(k\le2\) and every prime through 1000. At \(p=3\), \(k=0\), so the actual
margin is three.

The ten-speed citation is used only to supply a phase for \(C/p\) at
clearance \(1/11\ge1/14\). For the dilation corollary, the referee separately
checked \(\gcd(C_0)=1\), \(g\mid h\), preservation of ternary-unit tails, and
primitivity after division. Scale one invokes THM-4442; larger normalized
scale invokes the proposition.

## Sharp hostile

The one-label \(p=2\) margin is attained. Put

    C'=(1,2,3,4,6,7,8,9,10,11),
    T=(1,4,10),
    y=23/56.

Then \(y\in G_{C'}\), the primitive thirteen-speed row is \(6C'\cup T\),
and exactly label \(j=3\) among the six common lifts survives all tails.
Thus the theorem's worst arithmetic boundary is real rather than an artifact
of summing three separate caps.

## Engine separation and coverage

The clean-room script imports neither the primary verifier nor a repository
geometry engine. It performs:

- event-cell sweeps for every translated grid order \(1\le m\le210\);
- every prime and admissible tail-divisibility word through \(p=1000\);
- every distinct ternary-unit tail triple through 20 at \(p=2,3,5,7\);
- all 286 bounded-body gcd checks; and
- 292,500 normalization rows.

This totals 2,189,781 explicit checks. Normal and optimized transcripts
match.

    python -B 04-computation/lrc14_primitive_ten_pack_descent_thm4446_independent.py
    python -B -O 04-computation/lrc14_primitive_ten_pack_descent_thm4446_independent.py
