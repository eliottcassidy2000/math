# Ten actual sources, minimal selectors, and a composite-source repair

Status: **PROVED scoped**. Every source in this packet is a literal Boolean tournament, not a weighted orientation or a formal profile. Actual rows and independent endpoint/path-cover counts are recorded in `bank.json`.

For a tournament X let p_k(X) count unordered partitions into k nonempty directed paths and put

    F_X(r)=sum_k p_k(X)(r)_k,     (r)_k=r(r-1)...(r-k+1).

This counts placements into r distinguishable ordered path bins, with empty bins permitted. H(X)=F_X(1). The inherited identities are

    F_(X ordinal Y)=F_X F_Y,
    F_(X*)=F_X,
    F_(G(X))(r)=r^2 F_X(r)+r(F_X(r+1)-F_X(r-1)),

where G(X)=C3[X,1,1] has actual arcs X->c->d->X. Therefore H(G(X))=h+q, with h=F_X(1), q=F_X(2). See the scoped homogeneous-context equivalence in `adaptive_state_induction_20261103` and the path-bin calculus in `session159_constructions_20261103`.

## B1. An exact ten-symbol profile alphabet

The bank symbols have H values(3,49,343,63,5,11,13,35,77,91), corresponding to(a,c,d,e,p5,p11,p13,b5,b11,b13) in [NORMAL_FORMS.md](NORMAL_FORMS.md). The H343 source is given explicitly; the others use the inherited S156 two-apex witness bank. With y=r^2:

| H | Order | F(r) | F(2) |
|---|---:|---|---:|
| 3 | 3 | r(y+2) | 12 |
| 49 | 7 | r(y^3+16y^2+28y+4) | 872 |
| 343 | 9 | r(y^4+28y^3+148y^2+152y+14) | 10076 |
| 63 | 8 | r^2(y+2)(y^2+14y+6) | 1872 |
| 5 | 4 | r^2(y+4) | 32 |
| 11 | 5 | r(y^2+8y+2) | 100 |
| 13 | 5 | r(y^2+8y+4) | 104 |
| 35 | 7 | r(y^3+12y^2+20y+2) | 676 |
| 77 | 7 | r(y^3+18y^2+48y+10) | 1108 |
| 91 | 7 | r(y^3+18y^2+58y+14) | 1196 |

These ten polynomials are multiplicatively independent. Take the displayed non-r factor for each, using only y^2+14y+6 as the private factor for63. The ten selected factors are pairwise coprime; `bank.json` supplies45 integer Bezout identities. None vanishes at y=0. Every factor except y+2 occurs in only its own source; it determines that source's multiplicity in any product. Once the63 multiplicity is removed, y+2 determines the3 multiplicity. This proves the claim without requiring each selected factor to be irreducible.

Thus on ordinal words from this bank, the full profile remembers exactly the ten multiplicities. It forgets their order. Evaluation at r=1 forgets more, with fibres described by the arithmetic rewrite system. Ten is a number of multiplicities in this particular source bank, not a universal carrier dimension; triangle attachment produces new profiles.

## B2. Primitive attachments: minimal seven states, then thirteen

Call an ordinal word primitive when gcd(h,q)=1. Its forbidden single letters are{3,13,63,91}. On the six remaining letters, the only forbidden unordered pairs are

    (5,11), (11,35), (11,343), (343,77).

This follows directly from gcd(product H_i, product q_i)=1 iff gcd(H_i,q_j)=1 for every chosen ordered pair i,j. Letter49 is isolated. The maximal permitted alphabets are

    {11,77,49}, {5,35,343,49}, {5,35,77,49}.

By inclusion-exclusion the number of primitive words of length m is exactly

    2*4^m - 2^m.                                      (B2.1)

For a prefix retain the union of the forbidden neighbors of its letters. Precisely six such live sets are reachable; add one rejecting state. Distinct live sets have different forbidden next letters, so one letter separates them; the rejecting state differs by empty-word acceptance. This is a **minimal seven-state DFA** for the word predicate. `bank.json` gives its entire transition table. The seven-leaf native hub compiles every transition in two native moves, retaining the actual source word externally.

On a live word q/h modulo3 equals(-1)^j, where j is the number of letters other than5. Therefore H(G(X)) avoids3 exactly when j is even. The exact accepted-word count becomes

    (2*4^m - 2^m + (-3)^m + (-2)^m - (-1)^m)/2.       (B2.2)

The signed term is inclusion-exclusion on the same three alphabets, assigning weight+1 to5 and-1 to all other letters. The refined DFA has exactly13 states: six live forbidden sets times two parities, plus reject. All12 live combinations are reachable because49 toggles parity without changing the set. States of different parity are separated by the empty word; states of equal parity but different sets are separated by a distinguishing letter, followed if necessary by49. Thus13 is minimal for this precise acceptance predicate. This thirteen is distinct from the13 feasible divisibility snapshots in C2; no equivalence of their operation languages is claimed.

## B3. Rational class generating functions and a hidden register

Let

    B(z)=z^3+z^4+2z^5+4z^7+z^8+z^9.

All bank blocks are strong and have at least three vertices. In G(X), c is the unique vertex of outdegree1 and d the unique vertex of indegree1. Deleting them recovers X, whose ordered strong components recover the entire bank word. Consequently nonempty words give distinct ordinary strong tournament classes, enumerated by

    z^2 B(z)/(1-B(z)).                               (B3.1)

Equivalently the word counts have initial w_0=1, w_n=0 for n<0 and recurrence

    w_n=w_(n-3)+w_(n-4)+2w_(n-5)+4w_(n-7)+w_(n-8)+w_(n-9)

for n>=1. The wrapped count at order n is w_(n-2), excluding the empty word. For example this family has six classes at order10. These are an injectively decoded subfamily of A000568, not its total values.

Likewise G(X ordinal X*) and G(X ordinal 1 ordinal X*) are self-converse: swap X with its converse copy and swap c,d, fixing the middle vertex when present. Their exact class series are respectively

    z^2 B(z^2)/(1-B(z^2)),  z^3 B(z^2)/(1-B(z^2)).     (B3.2)

Their counts are h^2+q^2 and h^2+2q^2. Recovery is again intrinsic: identify c,d, take the first half of the ordered components. These formulas preserve the size clock and actual block sources.

There is a particularly explicit equal-profile bit register. Let A be the two-component word[63,63] and B the five-component word[49,3,3,3,3]. Both have scalar value3969, but different profiles. Encode bit0 by AB and bit1 by BA. Each bit packet has seven components and35 vertices. For every m-bit word w:

    |X_w|=35m,
    H(X_w)=3969^(2m),
    F_(X_w)=(F63^2 F49 F3^4)^m.

All profiles agree, but the component sequence decodes every bit. Wrapping gives2^m distinct strong classes of order35m+2 with the same entire profile. The two reflected wrappers give2^m distinct self-converse classes of orders70m+2 and70m+3, again with a common profile at each order.

The general reason for these collisions is inherited: module reversals preserve labelled cyclic triples and all labelled-subset Hamilton/path-cover data, as established in `session81_counterexample_families_20261028/LIFTING.md`; the noncommutative/cocenter distinction also appears in `session36_information_curl_contexts_20261023`. The present contribution is the concrete arithmetic-relation codec, its intrinsic decoding, and exact size rates. A cursor reading ordered components distinguishes what homogeneous module contexts cannot.

## B4. One composite source opens a missing7p lane

First classify the original bank's limitation for this exact operation G. No nonempty word X over the ten symbols has H(G(X))=7p for an odd prime p!=3.

If a7-divisible letter occurs, h is divisible by7 but q is not, by the table; hence h+q is not divisible by7. Otherwise only{3,5,11,13} remain. A3 or13 letter, or a mix of5 and11, supplies a common divisor3,13 or5 of h and q. Such a divisor precludes7p, except possibly values91 or35; words of length>=2 have h+q>=153, and the four length-one values15,37,111,117 are not7p. The remaining primitive possibilities are pure5 and pure11 words. For pure11, q/h=4 mod7 has odd order3, so it never equals-1. For pure5, q/h=5 mod7 equals-1 exactly at length m=3 mod6. Then5^m+32^m is divisible also by37, and exceeds7*37. This finishes the complete family classification.

Now add the literal strong H25 source

    rows=(14,28,40,32,45,3),
    F25=r^2(r^4+12r^2+12),   F25(2)=304.

The actual eight-vertex G(X25) has H=25+304=329=7*47. Trial division certifies47 prime. One additional sourced block suffices, and zero cannot, by the preceding classification. This is minimum number of extra blocks for this bank/operation; no global minimum order claim is made.

The structural point is that25=5^2 in the arithmetic monoid, but F25 is not F5^2. Replacing the actual strong25 source by the ordinal pair of5 sources changes q from304 to1024 and the source order from6 to8. Arithmetic irreducibility is therefore insufficient for choosing a complete attachment library.

## B5. An all-order recurrence with compatible prime-power clocks

For every t>=0 form the actual word X_t=X25 ordinal (X5)^(ordinal 6t), and wrap it by G. Its order is24t+8 and

    H_t=25 A^t+304 B^t,  A=5^6=15625, B=32^6=1073741824.

Thus

    H_(t+2)=(A+B)H_(t+1)-AB H_t,
    sum H_t z^t = (329-(25B+304A)z)/((1-Az)(1-Bz)).

Every H_t is divisible by7. For Q_t=H_t/7 and all e>=0,

    Q_(t+7^e)-Q_t = 5*7^e mod 7^(e+1).             (B5.1)

Indeed A,B=1 mod7, and their binomial first-order increments give coefficient
`[25(A-1)+304(B-1)]/7=5 mod7`. Hence there is a unique phase modulo7^e with Q_t=0 mod7^e, lifting one digit at a time. In particular49|H_t iff t=6 mod7.

Modulo11, H_t=3*5^t+7 and5 has order5; hence11|H_t iff t=1 mod5. Modulo13 both A and B equal-1, so H_t=4(-1)^t and13 never divides. The first two exclusion charts combine compatibly with period35. These statements preserve the actual sequence index; an unordered set of forbidden residues is not used as an orbit substitute.

The first terms are329,326417905121,350488137406584996329. Only the first quotient47 is asserted prime. The second quotient has factor11. Producing infinitely many prime quotients, or covering every requested prime/7p value, remains open.
