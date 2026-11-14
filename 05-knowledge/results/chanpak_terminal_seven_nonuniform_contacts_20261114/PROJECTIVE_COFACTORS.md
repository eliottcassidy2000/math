# Every fixed terminal-seven axis is locally complete on its cofactor

The [recent simultaneous congruence theorem](https://github.com/eliottcassidy2000/math/blob/11d40d1d10bef0d270e7aa3fb8edad094110a116/05-knowledge/results/chanpak_mrb_leech_clock_20261114/FULL_CONGRUENCE_CLOCK.md) for the five-clock
\((4,2,5,3,7)\) can tolerate a new left readout. Replace the
ordinary \((1,1)\) by any primitive integer row \((u,v)\) with
\(u\) odd and \(v\not\equiv0\bmod3\). The same physical clock
states still make \(ua+vb\) hit every locally compatible odd
residue modulo every modulus. Positivity of the readout is only
needed if interpreting it as a realized positive count.

The good-group and dyadic proofs require essentially no change:
\((u,v)\) is primitive on every good part; after setting all
terminal widths large at 2, the binary correction multiplies
the terminal baseline by the odd unit \(u\). At the bad part, the
compressed clock symmetry and all the on-line factors are
unchanged. The sole new requirement is the finite mod-63 terminal
guard test with the readout \((u,v)\).

It suffices to check the primitive projective rows modulo 9 and
7. Their counts are 12 and 8: mod 9 use \((1,j)\) for every
\(j\), plus \((j,1)\) for \(j=0,3,6\); mod 7 use \((1,j)\) plus
\((0,1)\). Keep the 9 mod-9 rows whose second coordinate is a
unit mod 3, making 72 combined cases. For **every** such row
and each of the 63 desired output residues, the same six width
variables (seed plus five terminal widths) have a choice among
the six powers of 2 mod 63 with unit derivative for each of
3 and 7. The attached replay gives an exact finite certificate.
The excluded projective rows cannot pass: the first state
coordinate is always zero mod 3 on completing a full five-clock,
as is seen by the product of its five matrices, for any odd width
choices. The derivative lifting argument and all synchronizations
therefore apply exactly as for ordinary readout to the 72 good rows.

For every \(y=4\cdot8^j\), the primitive positive row
\((A_y,B_y)=((57+68y)/7,(68+81y)/7)\) from the nonreset
terminal-seven rule meets these conditions: \(A_y\) is odd
and \(B_y\equiv2\bmod3\). Thus one can prescribe every finite
compatible congruence of the **cofactor** while enforcing literal
integer divisibility by 7 in the ordinary total. For example,
setting an arbitrary nonzero cofactor residue mod 7 ensures
that the total has exact 7-valuation one. This does not mean the
cofactor is prime. The same \(o(X^{0.406})\) ordered-clock
bound applies to each fixed positive row: its output is at
least a constant times \(a+b\) for every underlying word.
So not even this literal seven-divider overcomes the absolute
height problem. What it supplies is a robust interface for
the distinct \(7p\) frontier, without using an impossible
seven-path component or an exact rank-one reset.
