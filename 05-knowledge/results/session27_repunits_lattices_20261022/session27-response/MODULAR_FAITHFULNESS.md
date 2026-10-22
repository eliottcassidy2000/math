# Fixed-endpoint faithfulness modulo every modulus at least three

Status: PROVED by exact cancellation. This is a separate extension to the frozen `faithfulness-proof.md`, not a replacement of its valid averaged-endpoint identity. The fixed-endpoint improvement and modulus consequence were proposed by root and independently checked here. No novelty/priority claim is made.

## Exact fixed-endpoint identity

Fix a labelled k-core tournament Q, k≥2, distinct core labels u,v, and m=k−1 ordered transitive base vertices. As in the preceding faithfulness proof, choose the alternating path

R=(v,b_1,w_1,b_2,...,w_{k−2},b_{k−1}),

where the w_i are the remaining core vertices in any predetermined order. Its d=2k−3 cross edges receive independent signs x_e=+1 along R and −1 against R. Fix the extra edge b_{k−1}→u, rather than averaging it. Freeze all other attachment orientations arbitrarily. Let T(x) be the resulting tournament and χ(x)=Π_e x_e.

Then

Σ_{x∈{−1,+1}^d} H(T(x))χ(x)=2Q_uv.

If instead the extra edge is fixed u→b_{k−1}, the same sum equals 2Q_uv−2. Summing the two fixed-endpoint identities recovers the earlier averaged-endpoint value 2(2Q_uv−1).

Proof: as in the earlier proof, every Hamiltonian vertex word omitting an edge of R cancels, leaving only uR, R^rev u, Ru, uR^rev. At the core endpoint exactly one of uR and R^rev u is permitted by Q_uv; its signed contribution is s=2Q_uv−1. There is no factor two, since the extra edge is now fixed. At the base endpoint b_m→u permits only Ru, contributing +1; uR^rev is forbidden. The sum is s+1=2Q_uv. With the extra edge reversed, only uR^rev is allowed instead, contributing −1 because d is odd, and the sum is s−1=2Q_uv−2. This reasoning includes k=2, where R=(v,b_1) and d=1.

## Half-size query battery

The fixed-endpoint identity uses 2^d=2^(2k−3) response queries per unordered core pair, all at exactly m=k−1 base vertices. The entire labelled core is reconstructed from at most

C(k,2)·2^(2k−3)

queries, half the earlier averaged-endpoint bound. The query words depend only on k, the chosen pair, and predetermined fixed orientations, not on the unknown Q. This is an upper bound and does not claim optimality. For k≤1 no queries are required.

## Sharp modulus threshold for raw-count residues

Fix any integer modulus M≥3. Reduce the response counts modulo M before taking the signed sum. The resulting residue is 2Q_uv modulo M. Since 0 and 2 are distinct residues for every such M, one recovers Q_uv exactly by distinguishing those two residues. No inverse of 2 is used or required; the proof applies unchanged to even M, including M=4 and M=6.

Consequently, for every M≥3, equality of the fixed-labelled residue response functions w↦F_Q(w) mod M at exactly length k−1 implies equality of the cores. Consistent port relabelling and the word transform reverse/complement give the same faithful isomorphism and converse-isomorphism orbit statements as for exact counts.

For M=2, Rédei's parity theorem says every tournament has an odd number of directed Hamiltonian paths. Every response value is therefore 1 mod 2, for every core and every attachment word. Thus modulo-two raw responses cannot distinguish any two different labelled cores, even with unbounded base length. M=1 is trivial as well. For k≥2 the threshold M≥3 is therefore exact among positive integer moduli of the raw counts. For k=0,1 there is only one core and distinguishability is vacuous.

The parity input is classical, not proved or claimed new here. A verified primary exposition is Darij Grinberg and Richard P. Stanley, “The Redei–Berge symmetric function of a directed graph,” arXiv:2307.05569v1, Theorem 6.1 and its proof in Section 6: https://arxiv.org/html/2307.05569v1 . The paper's Theorem 7.1 gives a modulo-four refinement in terms of nontrivial odd directed cycles. Only the classical oddness conclusion is needed for the obstruction above; fixed-endpoint reconstruction uses the independent cancellation proof.

## One nonconstant response bit is enough

Do not interpret the modulo-two obstruction as an impossibility for every binary-valued observable. For a tournament T, define

q(T)=((H(T)−1)/2) mod 2.

Because H(T) is odd, q is precisely its second least significant binary bit and is equivalent to knowing whether H(T) is 1 or 3 modulo 4. For the fixed-endpoint experiment b_m→u, the stronger binary statement is

XOR_{x∈{−1,+1}^d} q(T(x)) = Q_uv.

Indeed χ has zero total sum because d≥1. Substitute H(T)=1+2r(T) in the exact identity and cancel the constant term, obtaining Σ_xχ(x)r(T(x))=Q_uv. Reduce this integer equality modulo two: each sign χ becomes 1, yielding the XOR identity. Thus one binary response coordinate, the second bit rather than the constant parity bit, already forms a faithful marked invariant through k−1 base vertices. For the opposite fixed extra orientation, the XOR is 1−Q_uv.

This is still an exponential-size structured query battery and not a cheap graph-enumeration algorithm. Knowing a residue or one bit of H may itself be difficult; this result counts oracle queries, not the time required to evaluate them. It neither proves universal two-port reconstruction nor produces metagraph wiring from scalar degree/count data.

## Audit and boundaries

`modular_faithfulness.py` independently constructs full tournaments and counts Hamiltonian paths by endpoint subset DP. It checks both choices of the fixed extra orientation, arbitrary deterministic frozen attachment masks, and every unordered pair in every labelled core of orders two through four. Additional order-five examples exercise both endpoint-label orders. For each query it checks classical parity, reduces the response individually modulo 2,3,4,5,6,8,11 before accumulating, and verifies the second-bit XOR identity. It never infers modular responses from the final exact sum alone.

Run normal and optimized Python and compare outputs. The expected complete audit comprises 840 fixed-endpoint arc tests and 27,528 direct DP queries. `modular_faithfulness.out` records completion, and `MODULAR-FAITHFULNESS-SHA256SUMS` freezes this addendum and its controls. Prior frozen proof/code files are unchanged.
