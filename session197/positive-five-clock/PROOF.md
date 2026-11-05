# A physical positive-five ternary family with prime and 7p witnesses

## Result and scope

For every integer m>=1, there is an explicitly constructed strong tournament G_m of order m+11 with

    H(G_m) = 1152*3^m+5.

The source G_1 has H=3461, a prime. The fifteen-vertex G_4 has H=93317=7*13331, with13331 prime. A literal operation on the same four physical ports returns H_next=3H-10 indefinitely.

These are exact graph constructions and two certified prime/7p examples. No assertion is made that infinitely many family members are prime or7 times a prime, that the family covers arbitrary prime atoms, or that it is globally minimal in amplitude, order or Hamilton count. The parameter m is always a positive integer.

## 1. Physical graph

Use the quotient tournament with outgoing row masks

    Q=[14,32,2,22,103,13,47]

in the order(T,P,U,a,R,b,V). Substitute the following nonempty modules:

- T is the transitive tournament of order m;
- P and R are transitive tournaments of order2;
- U has rows[8,1,3,6];
- a,b,V are singletons.

Every cross-module edge follows Q. Within each transitive block, earlier vertices beat later vertices. Concatenate the modules in the displayed order to obtain the canonical labels. A bit j in row i means i->j.

The quotient contains the directed spanning cycle(0,2,1,5,3,4,6,0). Because every module is nonempty, quotient paths connect any two different modules. A quotient cycle also connects vertices within the same module when their internal edge has the wrong direction. Thus every G_m is strong, regardless of the nontransitive internal structure of U.

The twelve-vertex source G_1 is

    [254,1028,1024,70,14,30,54,894,3711,3199,249,1279].

Its marked clone root is vertex0. The fifteen-vertex G_4 is

    [2046,2044,2040,2032,8224,8192,560,112,240,432,7152,29695,25599,1999,10239].

The verifier constructs these rows directly from the quotient rather than treating them as opaque asserted witnesses.

## 2. Exact all-order identity

There are two consistent proofs.

### Whole-module formula

The previously audited `session197-odd-cycle-path-polynomials/UNIT_OFFSET_RIGIDITY.md` gives, for arbitrary nonempty modules in this quotient,

    H = A B_T(3)+B B_T(2)+C H(T),

where

    A=2 b_U b_V (b_P-h_P)(b_R-h_R),
    B=b_V(b_R-h_R)(4h_P-b_P)h_U
      +b_U(b_P-h_P)(4h_R-b_R)h_V,
    C=h_U h_V(b_P-3h_P)(b_R-3h_R).

Here h_X=H(X), and b_X=B_X(2) counts covers in two labelled bins, allowing empty bins. The four-vertex U has bin polynomial x^4+4x^2, hence(h_U,b_U)=(5,32). P,R each have(h,b)=(1,4), and V has(h,b)=(1,2). Substitution gives

    (A,B,C)=(1152,0,5).

Since B_(Tr_m)(d)=d^m, this proves the claimed identity for every m>=1. It uses the all-order whole-module theorem, not extrapolation from a few values.

### Exact root-run certificate

For the marked source root0, let a_j count words that use every other vertex once and j indistinguishable root symbols, with no adjacent root symbols and with every consecutive pair a directed arc. The complete source run profile at positive indices is

    (a_1,a_2,a_3)=(3461,3456,1152),
    a_j=0 for j>=4.

The direct repeated-root DP in the verifier computes these values without extracting them from clone differences.

A uniform mathematical bound explains the vanishing higher indices. The source root's out-neighbors are the P and U vertices and vertex a=7. Its in-neighbors are R's vertices8,9, b=10 and V=11. The complete crossing-edge list from out-neighbors to in-neighbors is

    (1,10),(2,10),(7,8),(7,9).

It is covered by the two fixed vertices7 and10. Each internal root-to-root fragment must use a crossing edge, and different fragments are vertex-disjoint. Thus there can be at most two internal fragments, or three root runs. The disjoint edges(1,10),(7,8) also confirm matching number two.

Replacing the root by m transitive clones lifts a j-run word in j!S(m,j) ways, by distributing the ordered clone vertices surjectively among the j runs. Therefore, for every positive m,

    H_m=a_1+a_2(2^m-2)+a_3(3^m-3*2^m+3)
       =1152*3^m+5.

This supplies an independently checkable finite run certificate together with an all-order degree bound. No deletion identity at m=0 is needed or claimed.

## 3. Literal four-port return

In the source labels retain the same four vertices

    ports=[8,9,10,11].

At every step add a left-born vertex: it loses to those four ports and beats every other old vertex. Keep all four old ports and do not expose the newborn as a retained port. This is L(mask15,drop4) in the usual four-port convention.

Each newborn has the root's former incidences to the exterior and beats every earlier root clone. Thus it is another transitive root clone. Starting at the twelve-vertex source, after m-1 births the canonical root block maps to the physical labels

    [12+m-2,12+m-3,...,12,0],

with the list interpreted as[0] when m=1. Every other canonical module vertex, in order, maps to the unchanged original labels1,...,11. This is an explicit all-edge isomorphism at every m.

Consequently the same physical interface returns

    H_(m+1)=3H_m-10.

The four ports are never changed. All graphs remain strong. This is a single-sign nonunit-offset return; it is not a common-state choice between two ternary digits.

## 4. Prime certificates and the 7-divisibility condition

Direct endpoint counts for m=1,2,3,4 are

    3461,10373,31109,93317.

For3461, the verifier tests every possible divisor from2 through floor(sqrt(3461))=58 and records every nonzero remainder. For13331 it tests every divisor from2 through115. These complete trial-division certificates prove that both numbers are prime: every composite positive integer has a divisor in the tested range. The full certificate arrays are in `witness.json`.

The second requested graph therefore has the exact factorization

    H(G_4)=93317=7*13331.

More generally,1152=4 modulo7, and3 has multiplicative order six modulo7. Hence

    7 divides H(G_m)  if and only if  m=4 modulo6.

The six residues for m=0,...,5 are2,3,6,1,0,4. This classifies divisibility by7, not primality of the quotient. In particular it does not assert that later such members are7 times a prime.

## 5. Reproduction and provenance

Run `python verify.py` and `python -O verify.py`. They emit the same deterministic receipt and identical `witness.json` bytes. The verifier uses only the standard library and explicit exception gates. It validates tournament rows, graph strongness, all four physical counts, the complete root-run vector, the exact crossing support, the four-port isomorphism, and the primality certificates.

The arbitrary-module identity is the frozen `UNIT_OFFSET_RIGIDITY.md`, SHA2563c116086bfd697eaeaf5fcbf46308f1e31cdf653eb2a651053f01bce46ce374b, independently reviewed in `session197-unit-offset-audit/REPORT.md`. Its inherited transitive quotient is `session197-hybrid-constructions/positive-clock-family/PROOF.md`, SHA256be3c701188a4a354cee62f21b6ee76bf172f67351028d88a872f4c95c9b865c3. The general root-run lifting argument is also proved in the frozen `MINIMUM_TERNARY.md`.

The new result here is the explicit nonunit +5 ray and its certified prime/7p graphs, with a literal four-port return. The companion `PRIME_OFFSET_BRIDGE.md` separately explains the exact limitations of prime-offset targeting in this quotient.
