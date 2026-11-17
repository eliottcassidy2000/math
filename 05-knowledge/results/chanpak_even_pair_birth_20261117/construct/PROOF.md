# Two independently read quotients with physical sevenfolds

The pair-envelope operator raises the slot index, but on odd-order inputs makes the ordinary count divisible by three. This packet changes the pair contacts instead. It gives three primitive induced two-vertex extensions of arbitrary modules, and one common-source pair of quotient reads with an exact integer inverse. The inverse is a test for this bank, not a solution of prescribed-profile existence or Chan–Pak.

## Three all-order induced ×7 templates

In each table row use outgoing adjacency masks with vertex zero the low bit. Uniformly replace the last vertex by any nonempty tournament B (its H and optional two-slot count are h and P). Append q beating precisely the tag mask in the old quotient (so no vertex of B is tagged), and append d beating q and losing to every other vertex. The old tournament is induced in the new one.

| old quotient rows | mask | old count |
|---|---:|---:|
| (118,56,18,37,8,20,62) | 52 | 18P+13h |
| (16,181,33,71,44,9,55,125) | 12 | 67P+9h |
| (350,68,16,6,42,207,28,95,254) | 14 | 109P+7h |

In every row the new count is seven times the old count. Both quotients are strong, and strongness survives replacing a quotient vertex by any nonempty tournament: a directed cycle through that vertex gives the necessary routes between any two substituted vertices. The slot has a unique exterior in-neighbor, both before and after extending it. Consequently every full path uses one or two slot runs. For fixed quotient the read is exactly αP+βh, since the number of ordered covers by two nonempty paths is P−2h. Two literal tests, B=TT1 and TT2, determine the coefficients uniquely; these are checked in the replay, along with TT3 and a directed triangle as independent controls. No search cutoff is asserted: a targeted finite search supplied the rows, and the two-run proof lifts each to arbitrary B and order.

The same architecture also makes the connection to the earlier marked pair criterion: if W counts ordered old two-path covers whose nonempty second path starts tagged and U counts those whose first path ends tagged too, the two appended vertices contribute H_new=H_old+2W−U. The two-slot checks above realize the exact constraint 2W−U=6H_old as a linear identity on the arbitrary input, rather than assuming a numerical ratio coming from one particular module.

## A pair of independent odd quotients from one module

In the last induced row, put an ordinal transitive triple before a common arbitrary B *inside the slot*. This multiplies P by 8 and leaves h fixed, yielding p=872P+7h and its induced count 7p. On another copy of the same B use the previously published 128P+h sevenfold scaffold (see the dual-flank packet), with lower and upper masks

    (782,780,776,768,751,975,783,847,16,256),
    (782,652,776,784,679,975,543,845,82,256),

and slot at vertex 9. Its sole exterior exit is unchanged, so its two coefficients are again certified by TT1 and TT2. The upper quotient changes exterior edges but keeps the common slot untouched. This independently gives q=128P+h and 7q. Thus all **four** p,q,7p,7q have strong physical witnesses, even though only the p lift is asserted to be a two-vertex induced extension in this packet.

The two rows have determinant −24, yielding the exact joint inverse

    P=(7q−p)/24,        h=(16p−109q)/3.

The factor 24 is a threefold interaction with the dyadic slot scaling of the three neutral vertices, not an identification of Hamilton counts with the arithmetical 3,4,5 labels. It provides immediate absolute obstructions: both inverse quantities must be integers, h must be positive odd, and P≥2h. When the source order n is specified the stronger P≥(n+1)h follows by cutting each full path at its n+1 positions. Writing X=P/h, the possible ratio p/q=(872X+7)/(128X+1) lies in the narrow interval (109/16,1751/257]. Thus an invertible residue read emphatically does not place all independently prescribed absolute pairs.

Modulo every m coprime to 6 the linear map is invertible. The inherited joint one-apex profile compiler (the mixed-pair packet's LOCAL.md) independently programs h,P modulo every such m, so this actual four-count bank can hit any prescribed residues of p,q there. For this modular statement one may also prescribe odd order: initialize the compiler with the even-length spine word 10, augment its finite-state action by the parity of length, and note that each group command used in the vertical program has net even length (W=B A^(-2), its conjugates, and their products). A formal inverse is implemented by a positive power of the augmented permutation, preserving its length parity. This does not constrain the absolute magnitudes.

There is an essential three-adic gate for prime quotients. P₂≡(-1)^n H modulo 3 by parity of the odd-pack polynomial. In even order both p,q vanish modulo 3, while in odd order both equal −h modulo 3. Thus a simultaneous-prime source must have odd order and h nonzero modulo 3. This complements, rather than removes, the coordinate obstruction of the repeated pair envelope.

An actual simultaneous example uses the strong seven-vertex B=(90,40,59,64,42,9,54), with (h,P)=(71,1012). The quotients are the two certified primes p=882961 and q=129607, giving 7p=6180727 and 7q=907249. Primality is elementary trial division. This is one finite witness, not a statement about prime infinitude or arbitrary p.

Replay `python replay.py` checks the literal inputs, strongness and sole entry/exit, two-run coefficients, induced containment of the three extensions, the common-source inverse and the primality receipts. Normal and optimized output coincide.
