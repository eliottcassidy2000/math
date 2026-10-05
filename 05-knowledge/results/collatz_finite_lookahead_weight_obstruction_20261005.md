# Finite lookahead weights on the sibling quotient: an exact obstruction and a local repair

**Status: PROVED scoped observer obstruction, variable-radius bound, and
guarded weight comparison; FINITE-EXACT implementation controls.** A global
positive summable paid base weight, and universal Collatz, remain **OPEN**.
The obstruction concerns the specified observations, not arbitrary finite
algorithms or a controller retaining the full source/cofactor.

## 1. Inheritance and the target inequality

The closest proved mechanism is the exact sibling-fibre elimination in
[three-bit sibling flow, C3](collatz_three_bit_sibling_flow_20261005.md).
Let \(S(n)=4n+1\), and let \(\mathcal B\) be the positive odd integers
with \(v_2(3b+1)\in\{1,2\}\). Every positive odd source decodes uniquely
as \(S^j(b)\), \(b\in\mathcal B\). At a nonroot base define

\[
 U(b)=S^{d(b)}(G(b)),\qquad G(b)\in\mathcal B.
\]

The decoder repeatedly removes \(S\) exactly when its current source is
five modulo eight. It terminates by ordinary size, independently of
Collatz. With fixed \(0<r<1\), \(0<\rho\le1\), a positive summable base weight
\(g\), extended by \(f(S^jb)=r^jg(b)\), satisfies the full killed
incoming-flow inequality if and only if

\[
 g(b)\le\rho(1-r)r^{d(b)}g(G(b))
 \qquad(b\in\mathcal B\setminus\{1\}).                \tag{1}
\]

The canonical earlier hostile was the actual path \(9\to13\) in
[the neighbour macros](collatz_2_10_42_neighbour_macros_20261005.md).
It cannot simply be reused for the quotient: its intermediate endpoint
13 is the sibling \(S(3)\), so quotienting changes that path. The present
proof checks every sibling depth instead of assuming it is zero.
The corrected near miss is that longer fixed lookahead might automatically
repair a register-only weight. The least-used sidecar is the remaining
precision around the rational cycle anchor, beyond the observed window.

Incoming commit `69f6b3903`,
[critical flow, C7–C8](collatz_three_bits_critical_flow_20261005.md),
already strengthens the class-obstruction mechanism: the plateau
\(107\to161\to121\to91\), with an additional incoming predecessor,
rules out a positive critical flow depending only on its specified current
height, climb register, and golden-reader state. It also extends the
sibling criterion to \(\rho=1\), the endpoint allowed in (1) here.
The present extension varies the observation window and its residue
precision; it does not claim that the golden reader is included in this
different observer. C8 independently supplies a computable nonnegative
summable canonical flow whose support is exactly the rooted component.
Its unresolved strict positivity is consistent with the obstruction below.

The inherited negative anchor \(-5\) of the word \((1,2)\) is already
used in [adaptive switch payment](adaptive_switch_payment_20261004.md)
and in [guarded pumping memory, section 5](collatz_guarded_pumping_memory_20261004.md).
Here its role is to construct positive guarded aliases, not to posit a
positive cycle. The board is **sibling depth / current height / future
climb register / cofactor residue / negative-cycle precision / refuel bill**.

## 2. Precisely what the observer retains

For fixed integers \(L\ge0\), \(M\ge1\), consider the following
observation of a base \(b\):

* its complete current bit length;
* at every \(G^i(b)\), \(0\le i\le L\), the exact climb register
  \(h_i=v_2(G^i(b)+1)\), the odd cofactor
  \(t_i=(G^i(b)+1)/2^{h_i}\) modulo \(M\), and the state modulo \(M\);
* the actual valuation and sibling-removal depth of every intervening edge;
* a flag if ROOT is reached, at which point observation stops.

Call the combined record \(\mathcal O_{L,M}(b)\). The state residue
is redundant once the exact register and cofactor residue are known,
but is included to make the comparison explicit. No outgoing ROOT edge
is observed. The full integer cofactor is **not** retained; together with
\(h_i\), it would reconstruct the source itself.

**Theorem.** There is no positive function of this observation alone
that is a base weight satisfying (1). This remains true if an arbitrary
positive correction of this observation multiplies either of the two
inherited source weights

\[
 \mu(n)=2^{-(2\lfloor\log_2((n+1)/2)\rfloor+1)},\qquad
 \nu(n)=\frac8{3\,4^{\operatorname{bitlength}(n)}}.
\]

The proof constructs infinitely many exact positive witnesses for each
fixed \((L,M)\).

## 3. An expanding square with indistinguishable observations

Choose

\[
 K\ge3\lfloor L/2\rfloor+4,\qquad
 s\ge\operatorname{bitlength}(M)+2,
 \qquad t=M\left\lceil\frac{3\,2^s}{M}\right\rceil.
\]

Then \(M\mid t\) and
\(3\,2^s\le t<(13/4)2^s\). Define

\[
 \boxed{n=2^{K+3}t-5,\qquad x=9\,2^Kt-5.}            \tag{2}
\]

They are primitive positive odd bases, and the two actual quotient edges
are

\[
 n\xrightarrow{1}3\,2^{K+2}t-7\xrightarrow{2}x.
\]

Both sibling depths are zero, so \(G^2(n)=x>n\) exactly. The complete
current bit lengths agree:

\[
 \operatorname{bitlength}(n)=\operatorname{bitlength}(x)=K+s+5. \tag{3}
\]

Indeed the two sources respectively lie between \(24\,2^{K+s}-5\)
and \(26\,2^{K+s}-5\), and between \(27\,2^{K+s}-5\) and
\((117/4)2^{K+s}-5\). All four bounds are strictly inside the same
binary height shell.

For a starting value \(2^bu-5\), the shadow formulas are

\[
 X_{2j}=9^j2^{b-3j}u-5,\qquad
 X_{2j+1}=3\,9^j2^{b-3j-1}u-7.                       \tag{4}
\]

Apply these separately with \((b,u)=(K+3,t)\) and \((K,9t)\).
The bound on \(K\) ensures that every displayed perturbation through
index \(L\) is divisible by eight. Therefore even-index states are
three modulo eight, odd-index states are one modulo eight, and all are
primitive bases. Their exact registers alternate \(2,1\); their
intervening valuations alternate \(1,2\); no sibling is removed.
Every state is positive and exceeds one. Thus (4) is an actual positive
quotient trajectory on the required finite window.

The difference between the two trajectories at index \(i\) is

\[
 D_{2j}=9^j2^{K-3j}t,\qquad
 D_{2j+1}=3\,9^j2^{K-3j-1}t.                         \tag{5}
\]

These differences are divisible by \(M\). They remain divisible by
\(M\) after division by the common register power, four or two,
because the displayed powers of two have at least the required spare
precision. Hence both state residues and odd-cofactor residues agree.
Together with (3), this proves

\[
 \mathcal O_{L,M}(n)=\mathcal O_{L,M}(x).              \tag{6}
\]

The two baseline atoms also agree. They have the same full height and
both are three modulo eight, so neither gamma index is at the exceptional
power-of-two boundary. A direct index bit-length calculation gives the
same conclusion without the measure-comparison theorem.

For any proposed positive weight factoring through this observer,
\(g(n)=g(x)>0\). Multiplying (1) along the two checked depth-zero
edges would instead give

\[
 g(n)\le[\rho(1-r)]^2g(x)<g(x),
\]

a contradiction. The identical argument applies to baseline times
observer-correction weights because the baseline atoms agree.
This also applies at the critical boundary \(\rho=1\), since the full
sibling-fibre factor \(1-r\) still forces strict progress along a base edge.
It does not depend on computability or summability of the proposed
function. Increasing \(s\) yields infinitely many distinct witnesses.

There is also a direct critical-flow corollary, without assuming a geometric
sibling extension. Let \(y=U(n)\), so \(U(y)=x\), and let \(f\) be
strictly positive on every nonroot odd source with \(\mathcal Kf\le f\).
The extra source \(S(n)\) is distinct from \(n\) and has the same
successor \(y\). Thus

\[
 f(x)\ge f(y)\ge f(n)+f(S(n))>f(n).
\]

Any observation-based equality \(f(n)=f(x)\) is impossible. This applies
the incoming C7 external-predecessor argument to the quantified aliases
above; it is not a new general criterion for critical flows.

Small control: \(187\to281\to211\) has the same current bit length
eight and register two at its endpoints. This witnesses the \(L=0,M=1\)
case with smaller integers than the deliberately uniform construction.

## 4. Some unbounded lookahead rules still fail

The construction also quantifies the necessary precision for this specific
observer family. Fix \(M\), put \(s=\operatorname{bitlength}(M)+2\),
and keep the corresponding \(t\) fixed. For each height \(N\ge s+9\),
choose \(K=N-s-5\). The pair (2) has exactly that current height and
aliases every radius satisfying

\[
 \boxed{L\le2\left\lfloor\frac{N-s-9}{3}\right\rfloor+1.} \tag{7}
\]

Consequently, if the observer uses a height-dependent radius \(L(N)\)
meeting (7) at even one such height, the corresponding pair blocks a
universal strict weight comparison. If it meets (7) infinitely often,
there are infinitely many witnesses. In particular every sublinear
radius \(L(N)=o(N)\) fails eventually; so do radii at most \(cN\)
eventually for a fixed \(c<2/3\).

Here \(N\) is bit length, approximately \(\log_2 n\). This is a
sufficient alias range, not a proved optimal threshold. The argument
does not rule out larger radii, an increasing residue modulus, the full
cofactor, additional rational-anchor precision, or arbitrary finite
computations on the exact source.

## 5. A constructive two-anchor correction, with its unpaid boundary

The missing precision in (4) can pay the displayed negative-cycle phases
locally. Define on positive odd integers

\[
 \zeta(b)=
 \begin{cases}
  v_2(b+5),&b\equiv3\pmod8,\\
  v_2(b+7),&b\equiv1\pmod8,\\
  0,&\text{otherwise},
 \end{cases}
 \qquad g_*(b)=\nu(b)8^{-\zeta(b)}.                   \tag{8}
\]

This is an explicit positive summable weight, since \(0<g_*\le\nu\).
The following **guarded** comparisons hold:

1. If \(b\equiv3\pmod8\) and \(\zeta(b)\ge4\), the next
   actual valuation is one, no sibling is removed, and its phase is one
   modulo eight with precision \(\zeta(b)-1\). The endpoint's bit
   length increases by at most one. Hence
   \(g_*(b)/g_*(G(b))\le4/8=1/2\).
2. If \(b\equiv1\pmod8\) and \(\zeta(b)\ge5\), the next
   valuation is two, no sibling is removed, and the next phase is three
   modulo eight with precision \(\zeta(b)-2\). The actual positive
   endpoint is smaller, so its \(\nu\)-atom is at least the source's.
   Hence \(g_*(b)/g_*(G(b))\le1/64\).

For example, \(11\to17\) and \(25\to19\) attain the respective
bounds. Both guarded families satisfy (1) with \(r=1/16\),
\(\rho=2/3\), since \(\rho(1-r)=5/8\) and their depths are zero.
The guards are precisely \(b\equiv11\pmod{16}\) and
\(b\equiv25\pmod{32}\). Their union has source mass \(5/256\)
under \(\nu\), and \(15/256\) under \(\mu\), by the inherited
exact cylinder formulas. This is guaranteed local weight payment, not
global positivity of a paid solution.

The same globally defined candidate fails at the explicit refuel edge

\[
 7\longrightarrow11,\qquad
 \frac{g_*(7)}{g_*(11)}=4\,8^4=16384.                \tag{9}
\]

Its arrival creates four units of the new anchor register from a source
assigned zero. Thus the correction repairs the specified phase edges
while locating a new incoming bill exactly. A longer negative-cycle
episode consumes this precision and eventually leaves the guard; neither
that exit nor ROOT is silently certified by (8).

## 6. Reproduction and finite scope

The [standalone exact script](../../04-computation/experiments/collatz_finite_lookahead_weight_obstruction_20261005.py)
implements the quotient independently by literal sibling removal, the
closed shadow formulas, the observer, the source construction, and the
guarded pole correction. It performs no computation on import and writes
no files by default. ROOT observations stop without adding a self-edge.

```text
python 04-computation/experiments/collatz_finite_lookahead_weight_obstruction_20261005.py
python -O 04-computation/experiments/collatz_finite_lookahead_weight_obstruction_20261005.py
```

The [saved output](collatz_finite_lookahead_weight_obstruction_20261005.out)
records 186,630 checks: 975 source pairs from all \(L=0,\ldots,64\)
and 15 stated moduli; 1,215 variable-radius pairs; all 6,144 guarded
pole edges among positive odd sources below 65,536; and 13 malformed,
nonbase, root-edge, and insufficient-height hostiles. Each displayed
shadow state is independently compared with actual quotient iteration.
The finite census verifies the implementation; the infinite statements
are proved above.

The constructive next obligation is the missing incoming refuel inequality,
with a summability bound that survives such insertions. The observer
theorem identifies a class of insufficient inputs to that search. It does
not establish that global paid weights cannot exist; the inherited
equivalence says that their existence is exactly the unresolved universal
coverage problem.
