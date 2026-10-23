# The proposed 13=4 cycle: exact exponent graph and quotient obstruction

Status: PROVED elementary statements. Group-theoretic provenance of the input integer is treated in the companion Monster note. No tournament-spectrum explanation follows from this arithmetic graph alone.

Let N=2^5·3^7·5^3·7^4·11·13²=39,046,339,332,000. The prime-to-exponent arrows are

13→2→5→3→7→4, and 11→1.

Therefore the proposed cycle is exactly the chain in the valuation map p↦v_p(N), closed by identifying the output4 with the input prime13. This is a precise interpretation of the supplied integer, not an arbitrary match of unrelated digits.

## 1. A set quotient is legitimate, but its arithmetic meaning must be stated

On the directed graph of the above arrows, identifying13 with4 produces the five-cycle

[13=4]→2→5→3→7→[13=4],

and a separate arrow11→1. There are five vertices on this cycle, not six. The graph quotient is well-defined because4 was an output-only vertex, so no second outgoing arrow at4 was prescribed.

This construction is not an equality of integers or primes. It preserves the selected incidence chain but destroys the prime/composite distinction and changes the endpoint into a recurrent state. No map to tournament isomorphism, Hamiltonian paths or a group action has yet been specified.

## 2. Modulo9 does not give a well-defined quotient of the full exponent map

If the identification is intended as13≡4 modulo9, the full prime-support map fails to descend: 11≡2 modulo9, but

v_11(N)=1, while v_2(N)=5,

and1 is not congruent to5 modulo9. The omitted11 arrow is therefore a decisive obstruction, not a harmless decoration.

More strongly, there is no nontrivial quotient ring Z/mZ identifying13 and4 on which the whole prime-support exponent map descends. Such a modulus must divide9. For m=9 the preceding witness applies. For m=3, inputs2 and5 agree but their outputs5 and3 disagree. Only the one-point quotient m=1 remains.

Even the selected five-state cycle is not the restriction of an integer-polynomial map modulo9. An integer polynomial respects congruence modulo3. Here inputs2 and5 differ by3, but desired outputs5 and3 differ by−2. Thus no integer polynomial can realize those five arrows modulo9. An arbitrary finite-state permutation can realize them; that is a weaker, explicitly combinatorial interpretation.

## 3. The cycle is sensitive to the group-theoretic normalization

The actual Monster2A centralizer is2.B, rather than B. When passing from the ratio using |B| to the involution-class index using |2.B|, the core integer changes from N to N/2. Its exponent at2 changes from5 to4. Under the same set identification13=4, the five-cycle becomes

[13=4]→2→[13=4],

with5→3→7→[13=4] feeding into it, and11→1 separate. Thus it is a two-cycle, not a three-cycle: the repeated representative4 is the same vertex as13. Counting printed symbols rather than quotient states would give a false period.

This does not invalidate the exact ratio based on B. It shows that the cycle cannot be called an invariant of the involution conjugacy class without choosing and justifying the normalization. The prime valuation strata remain meaningful, but the exponent-arrow dynamics are additional structure imposed on their integer labels.

## 4. Why this is useful

The correct positive connection is the exact prime-valuation comparison between two group orders. The cycle is a concrete way to display part of that valuation vector. The correct obstruction is that this display does not descend to the obvious arithmetic quotient, and it changes under the centralizer's factor2.

A proposed further connection must specify the source group object, target tournament or square-sum object, the map, and a preserved invariant. Shared labels7 or21 alone do not supply that map.
