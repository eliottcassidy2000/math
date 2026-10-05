# The 2, 10, 42 sequence: rooted centres, paid neighbours, and refuel debt

**Status: PROVED exact macro, payment, source-mass, and obstruction statements;
FINITE-EXACT implementation controls.** Universal positive Collatz remains
**OPEN**. The positive macro below lies inside inherited first-descent
coverage. Its contribution is an exact compressed family and certificate
transport, not a new convergence basin.

## 1. Inheritance and the coordinate that must survive

The closest proved mechanism is the fixed-target sibling ray in
[ternary Berggren recursion, sections 3–4](ternary_berggren_20260925.md).
There the sequence

\[
 a_k=\frac{2(4^k-1)}3,\qquad a_{k+1}=4a_k+2,
 \qquad k\ge0,
\]

already appears as the Berggren **height** sequence
\(0,2,10,42,170,\ldots\). Its map from the parameter \(k\) to height
is a 3-adic isometry and conjugates addition by one to \(a\mapsto4a+2\).
Every finite ternary quotient is covered, while the ordinary height image
is sparse. Those are different coverage predicates.

The canonical hostile here is the refuel family in
[effective prefix mass, P8](collatz_effective_prefix_mass_20261005.md).
The repaired near miss is that paying a valuation-one climb does not pay
the preceding event that replenishes its climb register. The least-used
sidecar is the signed **offset from the even height**. The elementary
binary identities retain it exactly:

\[
 a_k=(10)^k_2,\quad a_k-1=((10)^{k-1}01)_2,
 \quad a_k+1=((10)^{k-1}11)_2 \qquad(k\ge1).
\]

Thus almost identical compressed descriptions can have opposite
source-relative payment behaviour. No tournament orientation is needed:
the intrinsic directed relation in this package is the actual odd Collatz
edge, together with checked common-future receipts.

The working concept board is: sibling root rays; ordinary versus ternary
height; offset-dependent word macros; source atoms versus climb-register
statistics; immutable-source payment; and supplied certificate transport.

Write \(U(n)=(3n+1)/2^{v_2(3n+1)}\) on positive odd integers. ROOT has
an empty first-hit word; the formal edge \(1\to1\) is never stored.
A word lists successive exact halving exponents.

## 2. The rooted even centre and its exact atomic mass

For \(k\ge1\), put

\[
 r_k=a_k/2=(4^k-1)/3.
\]

The case \(r_1=1\) is ROOT. For \(k\ge2\), the word \((2k)\)
is an exact one-odd-edge first-hit receipt for \(r_k\), since
\(3r_k+1=4^k\). The even source \(a_k=2r_k\) first halves to
that source. Its ordinary first-hit time is one for \(k=1\), and
\(2k+2\) for \(k\ge2\). These are the inherited inverse ROOT ray.

Use the two previously proved full-support atomic measures on odd sources:

\[
 \mu(n)=2^{-(2\lfloor\log_2((n+1)/2)\rfloor+1)},\qquad
 \nu(n)=\frac8{3\,4^{\operatorname{bitlength}(n)}}.
\]

The first is the gamma-source measure from
[atomic prefix certificates](atomic_prefix_certificate_20261005.md);
the second is the mixture measure from
[effective prefix mass](collatz_effective_prefix_mass_20261005.md).
For \(k\ge2\), \(r_k\) has bit length \(2k-1\), while
\((r_k+1)/2\) has bit length \(2k-2\). Therefore

\[
 \mu(r_k)=\frac2{16^{k-1}},\qquad
 \nu(r_k)=\frac2{3\,16^{k-1}}.
\]

With the exceptional ROOT atoms \(\mu(1)=1/2\), \(\nu(1)=2/3\),
the entire already certified ray has exact mass

\[
 \boxed{\mu(\{r_k:k\ge1\})=\frac{19}{30},\qquad
        \nu(\{r_k:k\ge1\})=\frac{32}{45}.}
\]

After retaining \(r_1,\ldots,r_K\), the respective ray tails are
\(2/(15\,16^{K-1})\) and \(2/(45\,16^{K-1})\). This supplies an
infinite symbolic ROOT receipt with a computable mass remainder.
The full-source residual masses \(11/30\) and \(13/45\) remain positive;
large mass of a sparse rooted family is not full source coverage.

## 3. The plus neighbour pays, and transports a smaller child receipt

Set \(p_k=a_k+1\), and let
\(\epsilon_k=1\) for even \(k\), \(\epsilon_k=2\) for odd \(k\).
The following all-height macro is exact:

\[
 \boxed{p_k\xrightarrow{\ (1,\,2^{[k-1]},\,2+\epsilon_k)\ }
        y_k=\frac{3^k+1}{2^{\epsilon_k}}.}                 \tag{1}
\]

Here \(2^{[k-1]}\) means a run of \(k-1\) letters equal to two, not
one exponentiation-valued letter. First \(U(p_k)=4^k+1\).
For \(j=0,\ldots,k-1\), the intermediate integer after the next
\(j\) exponent-two steps is

\[
 3^j4^{k-j}+1.
\]

For \(j<k-1\), the next numerator is exactly four times an odd
integer. At \(j=k-1\), it is \(4(3^k+1)\), and reduction modulo
eight gives \(v_2(3^k+1)=\epsilon_k\). This proves every valuation
in (1). The final endpoint is positive and strictly below \(p_k\):
indeed \(y_k\le(3^k+1)/2<(2\,4^k+1)/3\) for \(k\ge1\).

The earliest descent is even shorter. The cases \(k=1,2\) are
\(3\xrightarrow{(1,4)}1\) and \(11\xrightarrow{(1,2,3)}5\).
For every \(k\ge3\), the first three letters are \((1,2,2)\),
the first two checkpoints exceed the immutable source, and

\[
 U^3(p_k)=\frac{27p_k+23}{32}<p_k.
\]

This is the existing exact cylinder \(n\equiv43\pmod{64}\), whose
least positive source already exceeds the threshold \(23/5\).
So (1) does not enlarge the generic first-descent atlas.

Its useful longer-range identity is the common future

\[
 y_k=U(c_k),\qquad c_k=3^{k-1}<p_k\qquad(k\ge2).        \tag{2}
\]

If a supplied, checked first-hit ROOT receipt for \(c_k\) is
\((\epsilon_k,\mathrm{tail})\), substitute the source word (1) to obtain
the ROOT receipt \((1,2^{[k-1]},2+\epsilon_k,\mathrm{tail})\) for
\(p_k\). None of the displayed source intermediates is ROOT; for
\(k\ge2\), even the joining endpoint exceeds one. The supplied tail
therefore preserves first-hit safety. The source receipt has \(k\)
more odd edges and \(3k+1\) more ordinary edges than the child receipt.
For \(k=1\), use the explicit word \((1,4)\); do not pad the empty
ROOT receipt with an outgoing root edge.

This is a conditional proof compiler. It checks the supplied child receipt
and source identity; it does not search the child's orbit or assume that
every power of three is rooted.

## 4. The minus neighbour refuels, grows, and restores the register

For \(k\ge2\), put \(q_k=a_k-1\). It is exactly the inherited
refuel source with \(h=2k-1\). The entire return to climb register one is

\[
 \boxed{q_k\xrightarrow{\ (2,\,1^{[2k-2]},\,2)\ }
       z_k=\frac{3^{2k-1}-1}{2}>q_k.}                    \tag{3}
\]

The first edge reaches \(2^{2k-1}-1\). Its next \(2k-2\) exact
valuation-one steps end at \(2\,3^{2k-2}-1\). The next numerator is
\(2(3^{2k-1}-1)\), which has valuation two because \(2k-1\) is odd.
Both endpoints have \(h(n)=v_2(n+1)=1\). Strict growth follows from

\[
 z_k-q_k=\frac{9^k-4^{k+1}+7}{6}>0 \quad(k\ge2).
\]

The word has \(2k\) odd edges, halving cost \(2k+2\), and ordinary
cost \(4k+2\). One further exact step has exponent \(2+v_2(k)\)
and endpoint

\[
 U(z_k)=\frac{9^k-1}{2^{3+v_2(k)}}.                    \tag{4}
\]

Equation (4) follows by factoring \(9^k-1\): its valuation is
\(3+v_2(k)\), since \(9=1+8\). It is a compression identity, with
no uniform payment or ROOT conclusion appended to it.

The omitted boundary \(k=1\) has \(q_1=1\). Its strict first-hit
word is empty; the formal word \((2,2)\) is not an admissible ROOT
receipt. The program treats this exception explicitly.

For every \(k\ge2\), the two neighbours \(q_k,p_k\) have identical
\(\mu\)-atoms and identical \(\nu\)-atoms. Their shared description
length and mass do not determine whether the selected block pays.

## 5. Same climb distribution, different height concentration

The relation between the two source measures is inherited from
[algorithmic Collatz measure](algorithmic_collatz_measure_20261005.md):
\(\nu(2m-1)/\mu(2m-1)\) is \(4/3\) when \(m\) is a power of two,
and \(1/3\) otherwise.

There is nevertheless an exact common climb-register law:

\[
 \boxed{\mu\{h(n)=k\}=\nu\{h(n)=k\}=\frac3{4^k},
        \qquad k\ge1.}                              \tag{5}
\]

For gamma coordinates, this is the cylinder
\(m\equiv2^{k-1}\pmod{2^k}\), of mass
\(4^{-k}+2/4^k\). For the mixture, it is
\(n\equiv2^k-1\pmod{2^{k+1}}\), of mass
\(8/(3\,4^k)+1/(3\,4^k)\). Both reduce to (5), so the common
mean register is \(4/3\).

What the register forgets is the ordinary representative. Conditional on
the cell in (5), its least source \(2^k-1\) has probability \(2/3\)
under gamma and \(8/9\) under the mixture. That source is exactly the
maximal initial valuation-one climb in its dyadic cell.

For completeness, the Shannon source entropies are three bits for gamma
and \(\log_2(3)+1/3\) bits for the mixture. The latter follows from
the ROOT atom \(2/3\) and the bit-length shell masses
\(2^{1-\ell}/3\), \(\ell\ge2\); their normalized conditional mean
bit length is three. The same climb law is compatible with different
source entropies. None of these entropy identities establishes payment.

## 6. A sharp obstruction to every climb-only correction

Let \(b\) be either \(\mu\) or \(\nu\), and let
\(f:\mathbb N_{>0}\to\mathbb R_{>0}\) be arbitrary. Put
\(w(n)=b(n)f(h(n))\). There is **no** \(0<\rho<1\) for which every
surviving odd edge satisfies \(w(n)\le\rho w(U(n))\).
Consequently no such correction satisfies the stronger incoming-flow
inequality \(\mathcal Kw\le\rho w\) from the inherited flow theorem.

The smallest displayed refuel return proves this directly:

\[
 9\longrightarrow7\longrightarrow11\longrightarrow17\longrightarrow13.
\]

Both endpoints have register one. Moreover
\(\mu(9)=\mu(13)=1/32\) and \(\nu(9)=\nu(13)=1/96\).
Multiplying the four proposed edge inequalities would give
\(w(9)\le\rho^4w(13)=\rho^4w(9)\), impossible.
This does not require summability or an asymptotic estimate. It also
rules out a weight depending only on the pair (bit length, climb register),
since both endpoints have bit length four.

The repair must retain an address-dependent refuel cost beyond that pair,
or use guarded blocks with a separate actual-source payment proof. Section
3 provides the latter on its explicit guard; section 4 prevents extending
it by ignoring the offset. The all-source positive summable flow remains
open, exactly as in the inherited discounted-flow equivalence.

## 7. Implementation, controls, and next obligation

The [standalone script](../../04-computation/experiments/collatz_2_10_42_neighbour_macros_20261005.py)
exports a two-integer `Neighbour(k, offset)` record, an exact integer-source
membership decoder, guarded literal replay, and `lift_plus_receipt`.
No code runs on import. The source decoder tests the exact power-of-four
identity; matching a residue alone never substitutes a source.

Reproduction:

```text
python 04-computation/experiments/collatz_2_10_42_neighbour_macros_20261005.py
python -O 04-computation/experiments/collatz_2_10_42_neighbour_macros_20261005.py
```

The [saved output](collatz_2_10_42_neighbour_macros_20261005.out) records
11,254 explicit checks: all formulas for \(k=1,\ldots,128\), complete
two-offset membership on 4,096 odd sources below 8,192, 64 register cells,
64 independent finite ray sums plus their exact tails, and 14 malformed,
wrong-valuation, and ROOT-padding hostiles. Sixteen finite supplied-child
receipts are generated independently for test purposes, then consumed by
the transport API; their resulting source words contain 590 odd edges.
Those finite tests do not assume that all powers of three converge.

The concrete next obligation is to price a refuel event using its exact
incoming address and retained terminal certificate. Repeating a climb-only
weight search cannot succeed by section 6. An inverse-fibre mass formula
may simplify the incoming operator, but a source must still be discharged
by the actual inequality or receipt; a compressed recurrence alone is
insufficient.
