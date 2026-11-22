# The exact unit-gap Collatz mirror

**PROVED.** Define the accelerated signed map on odd integers by \(C_\epsilon(n)=(3n+\epsilon)/2^{v_2(3n+\epsilon)}\), \(\epsilon=\pm1\). This section concerns ordinary integer legality, not just the affine maps on an odd-prime residue ring. Write \(x,x+2\) for the ordered inputs.

For a prescribed two-step script \((a,b)\), the integer branch expression is
\[
     C_\epsilon^2(x)=(9x+\epsilon(3+2^a))/2^{a+b}.
\]
The first exponents at \(x\) and \(x+2\) differ, and exactly one is 1, since the two numerators differ by 6. Suppose both scripts have the same sum \(S\) of division exponents, and let \(c\ge2\) be the larger first exponent. The output difference has numerator \(18\pm(2^c-2)\) over \(2^S\), with \(c\le S-1\). For \(S\ge4\) its absolute value is strictly less than 2, since \(16+2^{S-1}<2^{S+1}\). Hence the two outputs, which are both odd, cannot differ by two. For \(S=3\), necessarily \(c=2\), and the signed numerator has to be 16 (the other is 20). This leaves exactly the following two mirrored legal cases:
 
| sign | initial lower residue | scripts at \(x\), \(x+2\) | image pair |
|---|---|---|---|
| \(+\) | \(x\equiv9\pmod{16}\) | \((2,1),(1,2)\) | \(G_+(x),G_+(x)+2\), \(G_+=(9x+7)/8\) |
| \(-\) | \(x\equiv5\pmod{16}\) | \((1,2),(2,1)\) | \(G_-(x),G_-(x)+2\), \(G_-=(9x-5)/8\) |

For example, on the positive sheet the lower first exponent is 2 for \(x=1\bmod8\), and its second is 1 precisely when \(9x+7=8\bmod16\); this also forces the upper second exponent 2, since its predivision affine numerator is \(9x+23=8\bmod16\). On the minus sheet \(x=5\bmod16\) yields valuations \((1,2)\) and \((2,1)\) by substituting in the respective formulae. If the script sums differ, the difference as a function of \(x\) has nonzero slope; such a fixed pair of scripts cannot be a unit-preserving identity on a progression. This statement does not classify isolated coincidences for unequal totals.

The pair maps have fixed lower values \(x_*=-7,5\) respectively and satisfy \(G_\epsilon(x)-x_*=(9/8)(x-x_*)\). At least \(k\) consecutive applications are possible if and only if \(x\equiv x_*\pmod{2^{3k+1}}\): successive starts must have \(x_i-x_*\) divisible by 16 and each step loses exactly three binary digits unless \(x=x_*\). Thus the density of starts following \(k\) paired moves is \(2^{-(3k+1)}\) among all integers. The only permanent positive start for these two scripts is \(5\) on the minus sheet; the orbit \(5\leftrightarrow7\) is the constituent two-cycle, not the different greedy inverse cycle from \(\{4,11\}\).

## The base-3 return and the binary clock must be compared separately

For an odd prime \(\ell\ne3\), let \(q_\ell(a)=(a^{\ell-1}-1)/\ell\bmod\ell\) for units. The branch comparison on a displaced coordinate \(d=x-x_*\) is multiplication by \(9/8\), so after \(\ell-1\) paired steps the *unscaled* multiplier equals \(1+\ell(2q_\ell(3)-3q_\ell(2))\pmod{\ell^2}\). If instead the observation puts back the \(8^{\ell-1}\) removed by the binary clock, the multiplier is \(9^{\ell-1}\). Its return to 1 modulo \(\ell^2\) is exactly Mirimanoff's base-3 condition \(3^{\ell-1}\equiv1\pmod{\ell^2}\), since squaring is injective on the kernel of reduction mod \(\ell\). To see this on legal trajectories rather than formal maps alone, use CRT to start in the legal \(2^{3(\ell-1)+1}\)-cylinder with \(d\) an \(\ell\)-unit. Such starts are unbounded.

This separates **carry cancellation, clock removal, and legality**. The affine and reversal viewpoint inherits [the reciprocal Collatz carry](../collatz_five_mirrors_20260929.md) and the classical [Bernstein–Lagarias parity coding](https://websites.umich.edu/~lagarias/doc/bernstein.pdf). The tournament metagraph analogy is precise at the level of a marked unit seam, not an identification of the odd Collatz graph with its Hamiltonian-count projection. In particular \(E_p-1\equiv15\bmod16\) for \(p\ge5\) is in neither legal branch above.

## A direct, clock-independent Mirimanoff test at the Collatz preimage of zero

There is a second formulation that does not require removing a binary clock. Let \([u]\) be the Teichmüller representative of a residue \(u\) at an odd prime \(\ell\ne3\). In the plus branch, \(a=-1/3\bmod\ell\) is the unique residue mapping to zero. Then \([a]=-[3]^{-1}\), so for every \(j\ge1\)
\[
  v_\ell\bigl((3[a]+1)/2^j\bigr)
  =v_\ell(3-[3])
  =v_\ell(3^{\ell-1}-1).
\]
The last equality holds since \(3/[3]\equiv1\bmod\ell\), and raising to \(\ell-1\) preserves the depth. On the minus branch take \(a=1/3\). Thus **a one-step Collatz branch carries the canonical lift of its zero-preimage to zero modulo \(\ell^2\) exactly at a base-3 Mirimanoff prime**. It is independent of the branch's power of two. At \(\ell=11\), the plus branch has \(a=7\bmod11\) and \([a]\equiv7^{11}=40\bmod121\). The positive input \(x=161\) lies in this lifted class, and its actual first numerator \(3x+1=484\) is divisible by \(121\).

This is an exact extension to arbitrary prescribed scripts. Write their affine map \(F(x)=(3^b x+C)/2^S\), with \(b\) odd steps, \(S\) total halvings and integer carry \(C\). If \(\ell\nmid C\), put \(a=-C/3^b\bmod\ell\) and \(U=3^b/C\in\mathbb Z_\ell^\times\). Then \([a]=-1/[U]\) and hence
\[
       v_\ell F([a])=v_\ell(U-[U])
                    =v_\ell(U^{\ell-1}-1).
\]
The first-depth obstruction is the Fermat-quotient difference \(bq_\ell(3)-q_\ell(C)\), not merely the number of appearances of 3. It can be seen on positive actual integers: the backwards formulas \(x_i=(2^{a_i}x_{i+1}-\epsilon)/3\) for odd \(x_{i+1}\) show that the script is realized on one class modulo \(2^{S+1}\). At each backwards step the denominator is a 2-adic unit and the computed \(x_i\) is odd. Use CRT to choose positive inputs in this class and in \(x=[a]\bmod\ell^k\), or in any finite simultaneous collection of such prime-power classes. The resulting trajectory's last state has exactly the predicted truncated valuation. This is a lifting statement about a prescribed finite script, not about unbounded integer trajectories or tournament attainability.
