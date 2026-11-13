# Exact addresses, growth, and a sharp scope boundary for the new grammar

**PROVED corollaries of the physical chain-return theorem.** Fix one actual
primitive seed on `C(q,s)=(q+3s,s,q,4q+2s)`, with q odd positive and s
even positive. Only the macros L_k, k≥2, are permitted in this note.
Other neutral repairs, balanced exits, substitutions and arbitrary graphs
are outside this language.

## 1. Every word has a unique exact address

The [physical theorem](constructors/CHAIN_REUSE.md) gives

    S=2^(k−1)(3q+2s),        Q=q+3S.

Since3q+2s is odd, the output(Q,S) uniquely recovers

    q=Q−3S,   k=v_2(S)+1,
    s=(oddpart(S)−3q)/2.

The previous state must have positive q and positive even s. Its H is
strictly smaller. Thus from a fixed seed no two distinct words can have
the same final pair: peel their forced last letters until the seed. A
nonempty word cannot return to that seed because H strictly increases.
Checking equality with the seed before peeling gives a terminating
membership test for the constructed source language. A positive formal
point failing this test is not declared unrealizable by other methods.

For an exact-prefix seed from [EXACT_REPAIR.md](EXACT_REPAIR.md), the
formal preceding q is `c=35h−21f≤−7h<0`, since f≥2h. These seeds lie
outside the positive image, so different such seed pairs generate disjoint
address trees. Different left-factor graphs with the same(h,f) can give
the same seed pair; this statement does not identify their full graphs.

## 2. Exact constructed-family count by added vertices

L_k adds exactly k+2 vertices. Hence the generating function for finite
words, including the empty word, counted by added vertices is

    A(z)=1/(1−sum_(k>=2) z^(k+2))
        =(1−z)/(1−z−z^4).                          (1)

Writing A(z)=sum a_n z^n gives a_0=1, a_1=a_2=a_3=0 and
`a_n=a_(n−1)+a_(n−4)` for n≥4 (negative subscripts mean0).
Different words have different rooted profiles, so their actual graphs
are nonisomorphic as rooted tournaments. Thus(1) is an exact count of
this constructed rooted family, not just a count of formal instructions.
It is not the number of unmarked classes or of all tournaments realizing
a particular Hamilton total.

For a fixed L_2 loop the physical matrix is `[[19,12],[6,4]]`. Therefore
each of q,s,H obeys the exact recurrence

    u_(j+2)=23u_(j+1)−4u_j.

Its characteristic discriminant is513=27·19. From the smaller H9-derived
seed C(87,20), the first totals are147,3699,84489. This is a recurrence
for actual tournaments, obtained from the two retained coordinates. The
repeated513 is also a useful curiosity beside the earlier asymmetric
eight-core count, but no map identifying those two objects is presently
supplied; that numerical coincidence is not a dependency of a result here.

## 3. Indefinite reuse alone is still sparse in the target count

Let H=q+3s. The exact forward law is

    H'=4v H+(14v+1)q,       v=2^(k−1).

Every positive input has H'>4vH. Every output has Q>3S, hence Q/H'>1/2.
Once this latter region has been entered, each further step satisfies
`H'>11vH`. Define

    A=sum_(k>=2) (11*2^(k−1))^(−1/2)
      =(1+sqrt(2))/sqrt(11) <1,
    B=sum_(k>=2) (4*2^(k−1))^(−1/2)
      =(1+sqrt(2))/2.

For the fixed seed of count H0, sum H(word)^(−1/2) over all finite words.
The empty word contributes H0^(−1/2); the first-letter sum is bounded by
B, and every subsequent-letter sum by A. Therefore

    number{words with H<=X}
       <= (1+B/(1−A))*sqrt(X/H0).                  (2)

If the seed already satisfies q>3s, the constant improves to1/(1−A).
To justify the counting bound, each qualifying word contributes at least
one to sqrt(X)*H(word)^(−1/2); the geometric sum bounds the larger sum
over all words. The same bound applies to distinct scalar totals, which
can only be fewer than words.

Thus the new grammar is physically renewable and exactly decodable, but
its fixed-seed target set has density zero. This does not limit enlarged
languages with fresh exact-target repairs, balanced exits or new seed
families. It pinpoints the next missing ingredient: source supply or
additional target-changing moves whose combined image exceeds this
specific language. Primitivity and infinite branching alone do not
establish the Chan–Pak spectrum.

Run `python3 -B 05-knowledge/results/h17_exact_repairs_20261112/chain_language.py`.
It enumerates every word with at most24 added vertices for three named
seeds, tests exact inversion and profile uniqueness, and compares the
per-order word counts with(1). These are recurrence-level controls of
an already proved physical grammar, not fresh graph enumeration.
