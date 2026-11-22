# The 2–3–11 bracket, the three negative cycles, and the missing binary clock

**PROVED** (with small exact replays). This session distinguishes three
objects: an archimedean bracket, a modular multiplication orbit, and
a legal integer orbit. Identifying these without their extra coordinates
was the failure mode in the older bracket-only dynamics. None of the
results below classifies all negative-integer cycles or decides Chan–Pak.

## Two exact appearances of the same prime triple

In the [inherited odd-square theorem](../procgen_numerology_20260923_odd_square_brackets.md)
the first question is whether \(p\) and \(2p\) lie in a common bracket
\(((2m-1)^2,(2m+1)^2]\). The primes are exactly \(2,3,11\), with
\(2,4;\ 3,6;\ 11,22\). For \(m\ge3\) twice the least integer in the
bracket already exceeds its upper endpoint. For these primes the
inequalities are strict at the top; the open/closed convention therefore
does not affect the prime list. See the [Session 74 inheritance note](../session74_eleven_mechanisms_20261027/session74-inheritance/INHERITANCE.md)
for the full integer-multiplier and real-multiplier classifications.
This is an order/size statement, not a statement about residues.

A second exact prime set comes from the Collatz multiplier \(3/2\):
take its two nonunit primes, \(2\) in the denominator and \(3\) in
the numerator, and add the **least odd prime outside these two**
where the order of the numerator 3 does not grow from \(p\) to \(p^2\).
That least prime is 11: at 5 and 7 the orders are \(4,6\) and grow,
whereas \(3^5-1=242=2\cdot11^2\), with \(\operatorname{ord}_{11}(3)=5\).
But \(2^{10}-1=1023=93\cdot11\), so the base-2 clock *does* grow at 11.
In particular the order of the **ratio** \(3/2\) does grow from
10 to 110. The first-bracket proof and the first-order-anomaly proof
are independent. The useful connection is that the exceptional order
is visible in the **relative** clock of the actual long negative
cycle, not that it predicts the bracket or the cycle.

## Local return periods of the three known negative cycles

Negating inputs conjugates ordinary \(3x+1\) on negative integers to
\(3x-1\) on positive integers. Use the shortened step \(T(x)=x/2\) for
even \(x\), \((3x-1)/2\) for odd \(x\). The three known negative cycles
correspond to positive anchors \(1,5,17\) in this convention.
[Known cycles, Kaneda](https://www.fq.math.ca/Papers1/53-2/KanedaM2172015.pdf)
uses the unshortened version, whose periods are \(2,5,18\).
The odd-to-odd valuation words for the shortened version are below;
their sums count **all** halvings, including the halving on the
odd step. This is why the shortened periods are \(1,3,11\), not
\(2,3,11\).

| Anchor \(x_*\) | Valuation word | Number of odd steps \(k\) | Shortened period \(A\) | \(3^k-2^A\) |
|---:|---|---:|---:|---:|
|1| \(1\) |1|1|1|
|5| \(1,2\) |2|3|1|
|17| \(1,1,1,2,1,1,4\) |7|11|139|

Prescribe one cycle's actual parity word. Its affine return fixes
\(x_*\); on displacement \(d=x-x_*\) it is multiplication by
\(\lambda=3^k/2^A\). This is a permutation modulo every odd modulus
coprime to 3. At 11, \(2\) generates the units and \(3=2^8\bmod11\),
so \(e=\operatorname{ord}_{11}(\lambda)=10/\gcd(10,8k-A)\). Put
\(u=v_{11}(\lambda^e-1)\). The equality
\[
 v_{11}\bigl(\lambda^{10}-1\bigr)
   =v_{11}\bigl(3^{10k}-2^{10A}\bigr)
\]
and the binomial lifting law give \(u=1\) for the first two rows
and \(u=2\) for the third. Indeed the valuations of the two terms
after subtracting 1 separately are \(2+v_{11}(k)\) and
\(1+v_{11}(A)\). They differ for the first two rows. For the long
row the valuations agree at 2, and the residues modulo \(11^3\)
of the original powers are \(727\) and \(606\), whose difference is
\(121=1\cdot11^2\). Since \(10/e\) is an 11-unit the same depth holds
at exponent \(e\).

For an 11-unit displacement, the exact period of the affine
return map modulo \(11^j\) is therefore
\[
\begin{array}{c|c|c|c}
x_* & e & u & r_j\\ \hline
1&10&1&10\cdot11^{j-1}\\
5&10&1&10\cdot11^{j-1}\\
17&2&2&2\cdot11^{\max(0,j-2)}
\end{array}
\]
The formula follows because any return exponent must first be a
multiple of \(e\) and raising to the 11th power adds exactly one
valuation to a principal unit. If the displacement has valuation
\(v\), replace \(j\) by \(j-v\) when positive; otherwise the
displacement is already zero modulo \(11^j\). The long sheet has
a genuine **extra square of phase stability**: every unit displacement
returns in two copies of that word modulo \(121\), not just an
isolated numerical congruence for the anchor.

## Exact odd/even return and its absolute cost

Fix one row, and take \(x=x_*+d>0\) with \(d>0\), so the anchor
itself is not counted. Let \(M=2^s m\) with \(m\) odd and \(s\ge1\).
The following is an exact equivalence:
\[
\begin{gathered}
\text{the next \(rA\) shortened steps follow the prescribed cyclic parity}\\
\text{word and return to \(x\bmod M\)}
\\[2pt]\Longleftrightarrow\\[-2pt]
v_2(d)\ge rA+s,\qquad
m\mid (3^{kr}-2^{Ar})d .
\end{gathered}
\]
Indeed, through step \(i\) the displacement is \(3^{k_i}d/2^i\),
where \(k_i\) is the number of odd inputs seen so far. To match
all \(rA\) input parities requires and suffices that
\(v_2(d)\ge rA\). Equality of start and finish modulo \(2^s\)
then requires \(s\) more bits, since \(3^{kr}-2^{Ar}\) is odd.
For the odd modulus one can invert the denominator and obtain the
stated divisibility. If \(s=0\) and we require the final odd input
as well (an actual completed *accelerated* word), replace \(s\) by
1 in the binary bound; the odd condition is unchanged. On a bare odd
modulus a formal affine return does not check this final bit.

For \(m=11^j\) and unit \(d\bmod11\), choose the minimum possible
return count \(r=r_j\). The **least positive displacement** for the
actual residue return modulo \(2^s11^j\) is exactly \(2^{Ar_j+s}\):
the bound follows from the binary gate and is attained by this power
of two. At \(j=1\) the clock costs \(Ar_j\), in anchor order, are
\(10,30,22\). At \(j=2\) they are \(110,330,22\). For every \(j\ge2\)
the \(17\)-sheet beats both short sheets. For example, for
\(M=2\cdot121\) it has a real positive starting integer
\(17+2^{23}\) that follows the long word twice and returns modulo
\(M\), without returning to its original integer value.
This crossover compares **legal word shadows**, not total stopping
times or newly constructed Collatz cycles. At any fixed nonzero \(d\)
the binary condition eventually fails as \(r\) grows. This explains
why adding an odd-prime period computation is not a route around the
[finite valuation obstruction](../../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md).

More generally if \(m\) is coprime to 3, the odd condition is the
multiplicative period modulo \(m/\gcd(m,d)\). If \(3\mid m\) then
\(3^{kr}-2^{Ar}\) is a 3-unit, so **all** of the 3-part of \(m\) must
already divide \(d\). Thus 2 measures precision spent by division,
3 is the contracting odd factor, and 11 is the first odd unit prime
where the numerator clock sticks on passage to the square. These
are three mathematically distinct roles.

This is the useful transfer to the [Farey absolute inverse](../chanpak_farey_two_three_clocks_20261122/README.md): preserve
legality and absolute size separately from any local phase return.
The modular Collatz period is not a new exact \(E_{61}\) witness.

## A second triple with the same size but not the same dynamics

The [inherited \(\sigma^{-1}(144)\) classification](../borwein_ladder_owner_numbers_vsc_twelve_20261028.md) gives exactly
\(\{2,3,11\}\) and \(\{2,5,7\}\) among squarefree triple supports.
A reason to compare them beyond the product identity is that for
squarefree \(n\), the projective line over \(\mathbb Z/n\) has
\(\prod_{p\mid n}(p+1)=\sigma(n)\) states. Indeed CRT turns its
primitive homogeneous pairs, modulo unit equivalence, into
independent projective points over the prime fields.

On the odd parts \(33\) and \(35\), the map
\([u:v]\mapsto[2u:v]\) is a genuine permutation. Over each prime
field it fixes 0 and infinity, and every other point has period
\(\operatorname{ord}_p(2)\). The order lists are \((2,10)\) and
\((4,3)\), respectively. A pair of cycles of lengths \(a,b\)
splits into \(\gcd(a,b)\) cycles of length \(\operatorname{lcm}(a,b)\).
Thus the two 48-state odd quotients have cycle spectra
\[
 \begin{array}{c|l}
 33 & 4\times1,\ 2\times2,\ 4\times10\\
 35 & 4\times1,\ 4\times3,\ 2\times4,\ 2\times12
 \end{array}
\]
where \(c\times d\) denotes \(c\) cycles of length \(d\). For instance,
the second power fixes eight points at 33 and only four at 35.
To lift the comparison to even moduli, use the odd unit \(35\) modulo
66 or \(37\) modulo 70 (both are 2 on the odd part and 1 modulo 2).
The 3 projective points modulo 2 give **three copies** of the odd
spectrum, with both ambient sets having size 144. Literal multiplication
by 2 on the even projective space is *not* a global projective map:
\([1:0]\) at 2 would go to the invalid pair \([0:0]\). On ordinary
affine residues multiplication by 2 instead collapses the 2-coordinate
and permutes the odd coordinate. Thus the equal sigma fibre fixes the
number of projective states, not the dynamics or the needed lift.

Run `python3 replay.py` in this directory for the prime-power
orders, explicit binary return witness, and the two projective spectra.

## Further consequences in this packet

- [Clock-controlled rational cycles](RATIONAL-CYCLES.md) independently
  prescribe a binary residue, a ternary unit residue, and an arbitrary
  11-adic ball, in primitive rational cycles of either sign. The
  construction keeps the denominator sidecar. This proves local
  flexibility, not the existence of additional integer cycles.
- [The canonical Farey clock at 11](FAREY-ELEVEN.md) has a one-digit
  delay: fixing its binary coefficient modulo \(11^r\) forces its
  canceled ternary coefficient modulo \(11^{r+1}\). This is a
  non-product restriction even though each coordinate has a clean
  individual range.
- [All-depth sparse signed Farey search](FAREY-SPARSE.md) excludes
  any Collatz relation between the two positive E61 flanks that
  uses at most two odd steps, even with independent signs. This
  extends the bounded-depth search without treating Collatz
  dynamics as a necessary condition on a general pair.
