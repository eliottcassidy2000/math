# Two independently removable local factors at a fixed count

The four literal cores from the main note can be combined without interfering with one another, although gcds of integer-valued polynomials do not generally multiply.

Type 15 has two same-order versions A and B with f(1)=15 and f(2)=120,108, respectively. Type 4151 has two same-order versions U and V with f(1)=4151=7·593 and f(2)=88256,93104, respectively. A's uniform-context gcd contains 5, U's contains 7, and all nonempty evaluations of A are divisible by 5, of U by 7. The last assertions follow from the independently verified raw gcds (or directly from their automorphisms acting freely on every ordered path cover). Both 15 types have every evaluation divisible by 3.

Form an ordinal word with m blocks of type 15 and n blocks of type 4151, not both zero. Choose r of the m to be A, and s of the n to be U. The order and closed count are fixed at

    order = 5m+10n,       H=15^m · 4151^n.

The exact raw and uniform-context gcds of the word both equal

    3^m · 5^r · 7^s.

Proof: f multiplies under ordinal sum. The universal evaluation divisibilities just stated give the lower bound on the gcd. At j=1, v₃ is exactly m. At j=2, v₅ equals r, v₇ equals s, and v₅₉₃ is zero, by the four values 120,108,88256,93104. Both j=1 and j=2 have factorial weight 1, so they prove the sharp upper bound for both gcds. No unlisted prime can occur since H itself is an evaluation.

Thus for fixed block positions there is an exact rectangular family of local levels: replace A with B to remove one 5, or U with V to remove one 7, without changing H, order, or the other local level. The unavoidable power 3^m remains. This is a literal version of selective normalization in this construction class, not an equivalence relation of Galois representations, nor a guarantee of unchanged continuation ports. It works because the crucial exact valuations occur at common low readings; replacing the 15 blocks by arbitrary models of the same count would lose that proof.
