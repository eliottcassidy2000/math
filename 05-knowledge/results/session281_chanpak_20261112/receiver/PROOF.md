# SCC re-factorization and receiver completion: renewable six, then twelve

**Session281. PROVED; independently reviewed before recovery.** Chan–Pak remains OPEN.
This recovered text preserves the final theorem statements and proof from the
session. The all-order claims use the inherited rooted-fusion law; the literal
seed certificate is reproduced by the self-contained recovery verifier.

## 1. Inheritance and change of construction language

[Session280 escape](../../session280_arbitrary_switching_20261112/escape/PROOF.md)
provides rooted fusion, a renewable three-bank, and a paid doubling gate.
Its six-output union has unequal receivers; common fusion sandwiches cannot
repair that union. We instead change the physical deletion core. Equal
Hamilton-count factors have different endpoint-weight systems. Their root
cuts can align receiving coordinates while their scalar counts interlace.

The receiver is `(h,V,e)`, where h is the root-deletion count and V,e are
unsigned and signed root-ending path-cover masses. The distant
[adaptive receiver work](../../session48_conditioned_native_sections_20261024/receiver-repairs/THEOREM.md)
is a method analogy: move to another resource fiber rather than continue
an excluded common repair. No theorem about that distinct carrier is assumed.

## 2. Exact ordinal endpoint-weight compiler

For nonempty D, with H(empty)=1, define

    h(D)=H(D), Q(D)=sum_S H(S)H(D\S),
    K(D)=sum_S (-1)^|S| H(S)H(D\S),
    E_i(D)=sum_(S contains i) end_i(S)H(D\S),
    F_i(D)=sum_(S contains i) (-1)^|S| end_i(S)H(D\S).

Here end_i(S) counts directed Hamilton paths in the induced tournament on S
ending at i. For `D=D1 -> ... -> Dm`, and i belonging to Dj,

    E_i(D)=(product_(a<j) Q(Da)) E_i(Dj) (product_(a>j) h(Da)),
    F_i(D)=(product_(a<j) K(Da)) F_i(Dj) (product_(a>j) h(Da)).       (1)

**Proof.** A path ending at i cannot visit any later ordinal factor. Every
later factor must therefore be wholly in the second path. Each earlier
factor splits independently between two labelled path slots, contributing
Q. The marked factor contributes E_i. Concatenation in ordinal order is a
bijection, including empty slot pieces. Keeping the sign of the first slot
replaces preceding Q by K and E_i by F_i. This proves both identities.

If a new root is beaten precisely by I, its receiver satisfies

    h=H(D), V=h+sum_(i in I) E_i(D), e=h+sum_(i in I) F_i(D).       (2)

The singleton root-ending path supplies h. Each longer root-ending path
enters the root from exactly one i in I, giving (2), including its sign.
Thus matching receivers is a concrete digit-weight problem, but does not
alone match Hamilton counts.

Let A=(6,12,24,17,3), B=(2,4,1), C=(2,4,9,3), using outgoing row masks.
Their respective `(h,Q,E,F)` are

    A: (15,120,(21,21,21,21,21),(-3,-3,-3,-3,-3)),
    B: (3,12,(3,3,3),(-1,-1,-1)),
    C: (5,32,(9,7,5,6),(-1,-1,-1,-2)).

A singleton has `(1,2,(1),(-1))`; all four K values are0. No universal
K-vanishing theorem is needed. For `TT_p -> A -> TT_q`, let a,b encode the
incoming singleton cuts in increasing-position least-significant-bit order,
and r be the incoming core cardinality. Then

    V=15(1+a)+2^p(21r+120b).                              (3)

For `TT_p -> B -> TT_b -> C -> TT_q`, with singleton cut integers a,u,c,
incoming triangle cardinality r, and `j=9x0+7x1+5x2+6x3`,

    V=15+15a+15·2^p r+60·2^p u
         +12·2^(p+b)j+384·2^(p+b)c.                    (4)

These are all-size receiver digit systems. The other H=15 cores used below
are D=(36,1,2,7,15,30) and E=(44,1,2,6,15,30). Their exact endpoint weights
are printed by the verifier. They came from a finite order-six search, not
an assumed universal classification of every deletion count15 tournament.

## 3. A small renewable six-bank

Use these deletion graphs on vertices0,...,9. TT denotes an ordinal
transitive block. Append a root as new vertex0, shifting deletion labels
by1. The decimal bitmask I specifies deletion vertices **beating** the root.

|H|Root-deleted graph|I|
|---:|---|---:|
|979|TT3 -> E -> TT1|236|
|981|B -> TT1 -> C -> TT2|408|
|983|TT4 -> E|76|
|985|E -> TT4|326|
|987|TT3 -> C -> B|204|
|989|TT1 -> C -> B -> TT2|368|

These six eleven-vertex tournaments share `(h,V)=(15,1203)`. Their e values
are `(15,15,15,13,15,15)`. The mismatch is retained, not ignored.

The inherited context C1 (root word10, transitive deletion) obeys for every
pointed input

    H'=H+2V, h'=h, V'=4V−h, e'=h.                       (5)

**Receiver-completion lemma.** Any actual count bank sharing only `(h,V)`,
with no assumption on e, becomes fully renewable under one common C1. Its
count gaps stay unchanged, V becomes4V−h and e becomes h. This relaxes the
uniform seed-supply obligation from three common coordinates to two. The
action identity is inherited; the new use completes a mismatched receiver.

Indeed C1 has fusion packet `(h,Q,H,U,V,d,e)=(1,4,3,3,3,1,1)`; its coefficient
of the old e in the count readout is zero. Applying (5) to the six seeds gives

    H=3385,3387,3389,3391,3393,3395;
    (h,V,e)=(15,4797,15).                               (6)

These thirteen-vertex tournaments form a **renewable SIX-bank**. Every
common left unit-deletion context translates all six counts by the same
amount and retains a common receiver. Explicitly, inherited fusion gives

    H(L*T)=H_T+(H_L−1)h+(U_L−1)(V−h)+(d_L−1)(e−h),
    h'=h, V'=Q_L V+(V_L−Q_L)h,
    e'=K_L e+(e_L−K_L)h,

when h_L=1. This can be iterated for arbitrarily long finite programs.
The repair is at a new absolute location, not a transformation of the old
six outputs to those same prescribed counts. The verifier counts every raw
and completed literal graph; no minimality or exhaustiveness is asserted.

## 4. Twelve consecutive odd counts at infinitely many scales

After k further C1 operations on (6),

    V_k=4792·4^k+5,
    A_k=3385+(9584/3)(4^k−1)+10k.                        (7)

The bank is `A_k,A_k+2,...,A_k+10`, with h=e=15. For every even k≥0,
V_k≡27 mod30, so `z_k=(V_k−27)/15` is a positive even integer. Put P_k=2z_k.
For the physical one-apex context L_P whose root word is the ordinary binary
expansion of P, the inherited count recursion is

    g(0)=1, g(2v+1)=g(v),
    g(2v)=g(v)+v+(v mod2), v≥1.

Its shift on this bank is `T_k(P)=15(g(P)−1)+(V_k−15)P`. The exact carry gives

    T_k(P_k+1)−T_k(P_k)=V_k−15−15z_k=12.                (8)

The two translated banks therefore realize TWELVE consecutive odd counts
`N_k,N_k+2,...,N_k+22`, with `N_k=A_k+T_k(P_k)`. At k=0,

    z=318, P=636, g(P)=483, N_0=3,051,967.

The twelve first witnesses have23 vertices. Larger even k have order
`13+2k+bitlength(P_k)` and N_k tends to infinity. The proof, not finite
extrapolation, gives this all-scale conclusion. For these larger witnesses
the verifier uses exact fusion, not a claimed literal2^23 subset-DP replay.

**Boundary.** The second half incurs receiver defect
`−15·2^(bitlength(P_k)−1)`, as in the inherited paid doubling theorem.
The twelve-output union is not rearmed. This is a renewable six-bank and
infinite twelve-wide output family, not generic repeated doubling, an
eventual odd tail, or Chan–Pak completeness.

## 4.1 General modular preparation of shared(h,V) banks

**Preparation lemma.** Suppose r actual pointed tournaments have counts
`H0,H0+2,...,H0+2r−2` and share positive h,V; their e values may differ.
Since h is an actual tournament path count it is odd. If

    r ≡ 2^j V (mod h) for some integer j≥0,              (9)

then physical contexts produce2r consecutive odd counts at infinitely many
absolute locations. They do not generally produce a renewable2r-bank.

**Proof.** Choose n≥2 with `2^(n−1)V≡r mod h`. Periodicity of powers of2
modulo odd h allows n to be increased arbitrarily; h=1 is vacuous. Apply the
context whose root word is1 followed by n−1 zeros. Its inherited action is

    H'=H+2^(n−1)V, V'=2^n V−h, e'=h, h'=h.

This completes the receiver and prepares `V'≡h+2r mod2h`. Choose n large
enough that `V'>h+2r`. The inherited paid doubling criterion yields2r
consecutive odd counts. For infinitely many scales, choose fixed s≥2 with
`2^s≡1 mod h` and repeatedly apply the root word1 followed by s−1 zeros.
For receiver W, the update is `W'=2^s W−h`. The congruence W≡h+2r mod2h
is preserved because `(2^s−1)r≡0 mod h`. W and output baselines grow without
bound. The phase condition is sufficient, not a characterization of every
possible physical preparation. QED.

The new step is orbit preparation and its use with two-coordinate banks;
the context identities and paid doubling gate are inherited. A gcd test alone
cannot replace (9): h=7,V=1 has doubling orbit{1,2,4}, so r=3 passes the gcd
condition but fails this preparation. No abstract(h,V) pair is asserted
physically realizable. Our raw bank meets (9) with h15,V1203,r6,j1; n2 is
exactly its C1 completion.

## 5. Discovery scope and remaining structural question

The optional search uses deletion skeletons consisting of A with singleton
padding, B and C in either order with singleton gaps, and D or E with singleton
padding. At fixed order it groups exact root cuts by(h,V), then tests e or
uses completion(5). This is a named bounded architecture, not an all-order
tournament classification. The positive theorem requires only the displayed
literal witnesses.

The earlier A/B/C-only search found a renewable six-bank at22145..22155 with
(h,V,e)=(15,27459,15). D/E supply a much smaller seed. The open question is
whether equal-product SCC re-factorizations yield unbounded consecutive
receiver fibers, possibly with changing h and additional factorization types.
Equations(1)–(4) retain the exact cut weights; equal h alone is insufficient.

Reproduce with `python3 verify.py` and `python3 -O verify.py`. The recovered
verifier is self-contained, includes explicit exception checks, and prints
all literal graph rows, profiles and exact all-scale formula controls.
`search.cpp` is optional bounded discovery, not needed for the proof.
