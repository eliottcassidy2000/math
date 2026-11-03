# Prime-lane obstruction for square-graph transfer

## Scope and main result

Use the transfer in NOTE.md, with an actual Hamilton path p in the square-sum graph on [n] and at least one selected interval for each path gap. For all n≥4, all three exact source counts are divisible by 3, irrespective of how the square-sum interval embeddings and transitive padding are chosen. They exceed 21. Therefore none is prime, and none has the form 7p with p prime. This is a family obstruction, not an obstruction for tournament path counts generally.

Proof: If every gap of p were odd, the adjacent endpoints would have opposite parity and their sum would be an odd square, congruent to 1 modulo 4. Hence each edge lies wholly within one of the two residue sets {0,1} or {2,3} modulo 4. Both sets occur among 1,...,n for n≥4, so no spanning path can use only odd gaps. Thus some d is even. The interval factors are q_d=1+2^(d−1) and t_d=2^(d+1)+2(3^(d−1)−1). For even d, q_d=t_d=0 modulo 3. Therefore the products P and R are both divisible by 3 and so are the three sources P+32R, 3P+48R, 5P+64R. The first is at least 33 because P,R≥1; the others are larger. If 7p were divisible by 3 then p=3, making the number 21, excluded by the size bound.


## Exact 3-adic weights

Assume exactly one selected interval for each path gap, with arbitrary singleton padding, and let k be the number of even gaps. The same argument gives 3^k dividing every source and 3^(k+1) dividing the middle source. More precisely, for even d,

  v₃(q_d)=1+v₃(d−1),       v₃(t_d)=1+v₃(d/2).

Odd gaps give 3-adic valuation zero for both factors. For q use the 3-adic lifting-the-exponent formula for 2^(d−1)+1, since d−1 is odd. For t observe t_d/2=(2^d−1)+3^(d−1). When d≥4 is even, the first summand has valuation 1+v₃(d/2), strictly smaller than d−1; d=2 is checked directly. Thus, writing

  a=∑_{d even}(1+v₃(d−1)),    b=∑_{d even}(1+v₃(d/2)),

the exact valuations of P,R are a,b, independent of all singleton padding. If extra non-singleton intervals are included, the exact sums must instead range over all selected interval gaps, counting multiplicity; the divisibility conclusions remain valid. If a≠b, the three source valuations are exactly min(a,b), min(a,b)+1, min(a,b). If a=b cancellation can increase a valuation, so no equality is asserted. For the [15] path in NOTE.md, k=7, a=10 and b=9; its three source valuations are exactly (9,10,9).

There are no square-sum Hamilton paths on [2] or [3]. The singleton graph [1] has an empty gap sequence and is a necessary exception if arbitrary transitive padding is admitted: P=1, R=2^g and the chain is

  1+2^(g+5) → 3+3·2^(g+4) → 5+2^(g+6).

The first two terms are never divisible by 7, because powers of 2 modulo 7 are 1,2,4 and −1 is not one of them. The third is divisible by 7 iff g≡1 mod 3. At g=1 it is 133=7·19, at g=7 it is 8197=7·1171, and both cofactors are prime. This parameter family is exactly the S158 two-bin family, not a new isomorphism class of 7-divisible constructions: Z₂ is X'⊕T₁, and X' (labels 1,2,3,4) maps to S158's four-vertex X (labels 0,1,2,3) under 1→3, 2→0, 3→1, 4→2. Consequently G(Z₂⊕T_g) ≅ G(X⊕T_{g+1}). No infinitude of prime cofactors follows.

## A complementary mod-7 local obstruction

For any gap d, q_d is nonzero modulo 7, while t_d is zero modulo 7 exactly for d≡4 or 5 modulo 6. Hence if a path contains such a gap, R≡0 and all three sources are nonzero modulo 7. Every square-sum Hamilton path on [n] for 10≤n≤166 has such a gap, as seen by deleting all square-sum edges whose positive gap is 4 or 5 modulo 6:

* For 10≤n≤53 vertex 10 is isolated: all its possible neighbors through 53 are 6,15,26,39 and the respective gaps are 4,5,16,29.
* For 52≤n≤143 vertex 52 is isolated: its possible neighbors through 143 are 12,29,48,69,92,117; their gaps are 40,23,4,17,40,65.
* For 144≤n≤166 vertices 52,55,58 all have degree one, with respective sole neighbors 144,141,138. A Hamilton path cannot have three degree-one vertices.

For n>166 this particular certificate says nothing; the separate mod-3 obstruction above remains valid for all n≥4. Thus these sources contain no multiple of 7 in this finite n range, a stronger statement in that range than merely excluding 7p.

When no gap has d≡4,5 mod6, both P and R are invertible mod7. Writing rho=R/P mod7, the three source residues vanish precisely for rho=5,3,2 respectively. These are exact residue tests, not assertions that such rho arise from a Hamilton path or that the quotients are prime.

## Controls

Run prime_obstruction_check.py under Python and Python -O. It checks the gap-factor residue tables over two periods, the small low-degree certificate through n=166, the explicit isomorphism of four-vertex factors, and the 133 and 8197 factorizations with trial-division primality of their cofactors. The structural arguments prove the all-order results; these finite checks are controls.
