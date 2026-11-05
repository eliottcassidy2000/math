# Square-factor reduction and genuine omissions in the self-converse spectrum

Separate Session197 result. **The complete self-converse Hamilton-path spectrum through41, at arbitrary tournament order, is**

    {1,3,5,9,11,13,15,17,25,27,29,31,33,37,41}.       (1)

In particular19,23,35 and39 are attained by ordinary tournaments but never by a self-converse tournament. These are restrictions on a construction class, not new gaps in the unrestricted spectrum. No literature-priority claim is made.

## 1. Exact all-order SCC reduction

Let A be the unrestricted Hamilton-count spectrum, let A_sc be its self-converse subclass, and let B_sc contain1 and the counts of nontrivial strong self-converse tournaments.

**Proposition1.**

    A_sc = { d^2 h : d in A, h in B_sc }.             (2)

**Proof.** Every tournament has a uniquely ordered list of strong components C1->...->Ct, and its Hamilton count is their product. An isomorphism to its converse reverses that list, so Ci is isomorphic to the converse of C_(t+1-i). Their Hamilton counts agree. If t is odd, the central component is self-converse; if t is even, use h=1. Let d be the product of the counts in the first half of the list. Those components form an actual tournament U with H(U)=d. Thus H=d^2 h.

Conversely, take any tournament U with H(U)=d and a strong self-converse S with H(S)=h, or omit S when h=1. The actual ordered join

    U -> S -> converse(U)

is self-converse and has H=d^2 h. For h=1 a central singleton could equally be inserted. This proves both inclusions. QED.

In particular every attainable square d^2 has a self-converse witness, even if the original d-witness is not self-converse. Conversely every squarefree self-converse count must already occur in a strong self-converse component. Therefore a restriction of a search to self-converse graphs can lose prime targets and7p targets with p!=7; paired weak components cannot repair those lanes.

The SCC mechanism is elementary inherited tournament structure. The value of (2) here is its exact use in distinguishing a complete self-converse enumeration universe from unrestricted arithmetic coverage.

## 2. A finite complete all-order bound through41

The externally checked Busch strong-order theorem, cited in `PRIME_CUBES.md` and the original Session197 entropy packet, gives for every nontrivial strong component S

    5^(|S|-1) <= H(S)^3.

Since5^7=78125>41^3=68921, any strong component of count at most41 has at most7 vertices. This theorem is never applied to the full weak tournament or to singleton padding.

The complete self-converse normalized universe from `CONVERSE_MASS.md` has2^floor((n-1)^2/4) graphs at each positive order n. Exhausting it through n=7 and retaining H<=41 gives exactly the set in (1). Each graph is physically generated with an anti-involution and a Hamilton path. The verifier additionally exhausts the larger universe of all tournaments preserved by fixed reversal as an anti-involution, without forcing a Hamilton path, and obtains the same support. Every self-converse tournament belongs to this larger universe after relabeling: an anti-involution has at most one fixed vertex, hence the same permutation cycle type as reversal. The anti-involution existence proof is in `CONVERSE_MASS.md`.

Now let T be self-converse with H(T)<=41. By (2), H(T)=d^2 h. Since d is positive odd and d^2<=41, it is one of1,3,5. The central strong component (if any) has h in the finite set (1). Multiplication by9 or25 gives no additional value at most41: the only products are9,27 and25, already in (1). Thus the finite support is complete at arbitrary order. Every listed value is already physically attained in the enumerated universe, proving equality.

This also supplies an independent self-converse exclusion for7 and21; no unrestricted omission theorem is required for the reduction.

## 3. Explicit unrestricted escapes

Outgoing rows for actual tournaments with the four extra omitted values are:

    H19: [2,20,57,51,33,3],
    H23: [38,8,58,49,35,2],
    H35: [54,48,122,115,96,64,3],
    H39: [250,252,249,48,64,208,136,24].

The first three are strong. The last is the ordered join of a directed triangle and a five-vertex H13 tournament, so its product count is39. Direct subset DP and literal Hamilton-word enumeration independently verify all four counts. Equation(1) proves that none of these graphs can be self-converse, regardless of presentation.

The H35 witness has a short physical two-apex description: apexes a->b, increasing chain with orientation letters(1,0,3,3,0), where bit0 records a beating the chain vertex and bit1 records b beating it. The displayed rows use a,b followed by the chain. Thus the7p target35 is lost by the self-converse restriction even though it has a small concrete unrestricted construction.

The H19 and H23 examples are inherited physical two-apex banks from earlier Session197/S195 work; they are not claimed as new numerical constructions. The role of these rows is to make the model limitation explicit rather than merely record an unfilled census entry.

## 4. Reproduction

Run `python3 verify_converse_spectrum.py` and `python3 -O verify_converse_spectrum.py`, comparing receipts and witness files. It checks both complete generation universes through7, independently enumerates literal Hamilton words on every path-normalized graph, checks the exact integer bound and square closure, and checks all four unrestricted witnesses by both Hamilton methods. `converse-spectrum-witness.json` includes one actual self-converse witness and its anti-involution for every value in (1).

The bound and SCC reduction, not finite extrapolation, make the result all-order. Neither the self-converse mass theorem nor positive fixed-word counts can imply unrestricted prime/7p coverage, as these explicit omissions demonstrate.
