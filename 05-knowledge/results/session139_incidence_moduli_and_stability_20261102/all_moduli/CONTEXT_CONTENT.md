# Context content: every odd-modulus zero exists, but content is not multiplicative

**Status: PROVED; exact small graph receipts supplied.** This is a corollary/addendum to `PROOF.md`; its universal free-action input is explicitly inherited, not claimed new.

## 1. An exact integer invariant of all contexts

For a nonempty tournament T define its **context content** by

    g(T)=gcd {H(Q[T]) : Q is any homogeneous one-port context}.

The identity context is included, so g(T) is a positive divisor of the odd integer H(T). The sharp coefficient-detector theorem of Session137 immediately gives the finite formula

    g(T)=gcd_(1<=c<=|T|) [c!(c-1)! p_T(c)].                    (C1)

Indeed every response is an integer combination of the displayed weighted coordinates. Conversely, the actual triangular detector E_(c-1) has leading coefficient d_c and lower coefficients divisible by d_j. Inductively each weighted coordinate is an integer combination of detector responses. Thus the two ideals in Z agree, which proves(C1).

Because `p_T(n)=1` at n=|T|, the same formula gives the arithmetic order restriction

    g(T) divides gcd(H(T), n!(n-1)!).                          (C1a)

In particular no prime greater than n divides g(T). A zero module modulo M has order at least `kappa(M)=min{n>=1 : M divides n!(n-1)!}`. This is a necessary order bound, not a sufficient existence criterion at that order.

Consequently T is contextually zero modulo M, meaning **every** homogeneous-context response is divisible by M, exactly when M divides g(T). This invariant has an honest all-context quantifier and a finite computable certificate.

## 2. The zero contextual state exists at every odd modulus

The classical free action of Aut(T) on ordered Hamiltonian paths is recorded as LEM-003 in the inherited

`01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md`.

Every automorphism of T extends to Q[T] by fixing the exterior vertices. An automorphism fixing an ordered Hamiltonian path fixes each vertex at its position, hence is the identity. Therefore

    |Aut(T)| divides H(Q[T]) for every Q,
    |Aut(T)| divides g(T).                                  (C2)

This is a contextual lift of the inherited free-action lemma, not an assertion about automorphisms being recoverable from g.

For every odd M>=3 take the cyclic tournament R_M on Z/MZ, with

    i->j iff j-i mod M lies in {1,...,(M-1)/2}.

Its translation group has order M and is contained in Aut(R_M). Thus R_M is contextually zero modulo M. It is strong because the step-one arcs form a directed Hamiltonian cycle. At M=1 zero is the unique residue and K1 suffices. At every even M zero is impossible, since the identity-context response H(T) is odd.

Hence:

> The zero contextual state is represented by an actual nonempty tournament modulo M **if and only if M is odd**. For odd M>1 it has a strong witness of order M.

For an odd prime modulus p, the lower bound is p and R_p attains it: the minimum zero-module order is exactly p. More generally, any tournament on p vertices with an automorphism subgroup of order p has `v_p(g)=1` exactly, since `v_p(p!(p-1)!)=1`.

For arbitrary odd `M=product p^a>1`, ordinal-join a copies of R_p for each prime p. The product of the rotation groups, of order M, acts on the join, producing a zero module on `sum a*p` vertices. A final `C3[T,K1,K1]` wrap is strong and adds only two vertices while retaining that group. This supplies an explicit order bound, not the exact minimum at composites.

No bound-three substitution-width assertion is made for this cyclic construction. In particular it does not render the width-three p=5/p=11 constructions redundant.

### Symmetry is sufficient but not exact

For the regular five-vertex cyclic tournament R5, direct counting gives

    p_R5=(15,45,35,10,1),       g(R5)=15,       |Aut(R5)|=5.

The automorphism group has order5 because fixing one vertex fixes its two outneighbors in their directed order and its two inneighbors in their directed order; rotations supply all five possibilities. Thus(C2) can be strict. Because every zero module modulo15 needs at least5 vertices by(C1a), R5 also proves that the minimum zero-module order modulo15 is exactly5. R5 is already a strong contextual zero modulo15, although it has no order3 automorphism. This hostile blocks an automorphism-only interpretation of context content.

### The exact 11-vertex boundary: 33 versus 55

One bounded path-cover computation on the standard cyclic R11 with connection set `{1,2,3,4,5}` gives

    p_R11=(93027,787479,2034813,2376110,1481546,537075,
           117975,15840,1265,55,1),
    g(R11)=33.

It is therefore NOT a zero module modulo55. It does prove that the minimum zero-module order modulo33 is exactly11, because the prime11 lower bound is sharp there. The receipt is `cyclic11_probe.json`; no broad circulant census was performed.

The actual Paley tournament P11 supplies the different55 construction. Its explicit connection set is the five nonzero squares `{1,3,4,5,9}` modulo11. Each affine map `x->a*x+b`, with a in that set and b arbitrary modulo11, preserves its arcs. These55 distinct maps form an automorphism subgroup. Hence g(P11) is divisible55 and the minimum zero-module order modulo55 is exactly11. P11 is strong since step1 is an arc. This recovers the inherited THM-4747 affine5-by11 action as an exact homogeneous-context result. No claim about the full value of g(P11) or an all-circulant classification is needed.

## 3. Exact failure of multiplicativity: primitive H5/H11 factors

Use the actual W5 and W11 from `PROOF.md`, Section4, and Z5=W5->W11. Independent subset/partition counting in `verify_content.py` gives their profiles (entry c starts at1):

    W5:  (5,11,6,1),
    W11: (11,39,33,10,1),
    Z5:  (55,1545,7495,12310,8665,2898,474,36,1).

The factorial restriction(C1a) and the mod5 projector construction already prove that both W5 and W11 have content1 and their join has content5: primes5 and11 exceed their respective factor orders, while `g(Z5)` divides `gcd(55,9!8!)=5`. The literal profiles independently verify the same conclusion. In particular

    g(W5)=1,       g(W11)=1,       g(W5->W11)=5.               (C3)

The last equality also follows without taking a long gcd: Section4 proves all responses divisible by5, whereas the first two weighted coordinates are55 and3090, whose gcd is5. The first two weighted coordinates of W5 are5 and22, and of W11 are11 and78, proving both primitive claims.

Thus two modules that are individually context-primitive can develop a nontrivial common divisor of **every** contextual response after ordinal join. H and the full polynomial F are multiplicative under join, but the extracted integer invariant g is not.

The mechanism is the complementary support decomposition on F5*, not a scalar-factor coincidence: W5 vanishes on +/-1, W11 on +/-2, and their product vanishes everywhere. This is the exact arithmetic operation that connects the actual counts5 and11 here. It neither identifies the old arithmetic cubic's exceptional parameters nor proves any RH claim.

## 4. Verification and limits

Run:

    python3 -B verify_content.py > content.normal.json
    python3 -B -O verify_content.py > content.optimized.json
    cmp content.normal.json content.optimized.json

The standalone standard-library script constructs the literal graphs, computes unordered path-cover profiles, checks(C3), checks the rotation automorphisms for every odd order3..33, and verifies the contextual zero signatures directly for cyclic orders3,5,7,9. The all-M conclusion is the group-action proof, not extrapolation from those orders.

(C2) is only a sufficient symmetry mechanism for divisibility. No converse, general minimum-witness-order formula, general contextual-state realization theorem, exact Hamiltonian spectrum theorem, or generic formula for g(T->U) is claimed.
