# State-dependent two-ended growth still loses the lower branch after H=13

**PROVED all-order obstruction, with complete finite certificates.** In the
fixed-two-apex word model, an end extension means

    w -> u w v,                                       (1)

where u,v are arbitrary finite words. Their lengths and letters may depend
on the entire current state. Old apex names, their mutual arc, old chain
names, and all old incidences stay fixed. Any sequence of end extensions is
again an end extension.

There is **no** end-growth chain

    H(T(w_0))=13,    H(T(w_1))=41,
    H(T(w_2)) in {125,127},                            (2)

with w_1 an end extension of w_0 and w_2 an end extension of w_1.
This permits arbitrary word lengths, arbitrary source/sink padding, arbitrary
state-dependent end blocks, and nonstrong states. The forbidden second-stage
counts are exactly 3*41+2 and 3*41+4.

Thus changing the fixed prepend/append macros in `END_COMMUTATION.md` cannot
repair this particular lower second-generation branch while staying in the
end-extension architecture. The theorem does not restrict arbitrary insertion
inside the old chain, changing apex marks, changing old arcs, or rebuilding
a different H=41 witness.

## 1. Why an unbounded question has an exact finite cutoff

Three elementary facts retain the needed physical information.

**SCC product.** H(T) is the product of the Hamilton counts of its strong
components, because Hamilton paths concatenate in the condensation order.

**Induced monotonicity.** H(A)<=H(B) whenever A is an induced subtournament of
B. It suffices to add one vertex z. Every old Hamilton word has a legal
insertion of z: insert first if z beats its first name, last if its last name
beats z, or at a transition from a name beating z to a name beaten by z.
Choose the first legal place. Deleting z recovers the old word, proving an
injection. Iteration proves the claim.

**Contiguous strong core.** A strong component containing both apexes has a
contiguous interval of chain names. A chain name between two names of the
component is reachable from the component and reaches it. All external chain
names are then earlier source letters 0 or later sink letters 3.

The H=13 normal form in `APPEND_OBSTRUCTION.md` is

    0^p c 3^q,   c in {120,212,230,302,310,321}.         (3)

It follows from the exact H=13 strong-core cutoff, the one-apex exclusion,
and a complete 64-word check. Both physical apexes belong to this core.

Now suppose an end extension has H=41. The strong component C_1 containing
the old H=13 core has H(C_1)>=13 by induced monotonicity, and H(C_1) divides
41. Hence H(C_1)=41. It contains both apexes, so stripping its exterior source
and sink chains leaves a word that is an end extension of c in (3).
Busch's bound gives |C_1|<=7 because

    5^7 > 41^3.                                       (4)

Thus at most two new chain names can occur inside C_1, regardless of how much
neutral padding was present in either original word.

Next suppose an end extension has H=125 or 127. Its strong component C_2
containing C_1 has H(C_2)>=41 and divides that target. For 127, primality
forces H(C_2)=127. For 125, its only proper divisors are 1,5,25, all below
41, so H(C_2)=125. Therefore C_2 again contains both apexes and has a
contiguous chain word extending that of C_1. The cutoff is |C_2|<=10 since

    5^10 > 125^3,    5^10 > 127^3.                     (5)

The first-stage enumeration below shows that every possible C_1 has exactly
seven vertices. Hence at most three further chain names occur inside C_2.
This proves the all-order reduction; it does not apply a strong size bound
to an arbitrarily padded nonstrong parent.

## 2. The complete two stages

For each of the six three-letter cores c, enumerate all end extensions adding
zero, one, or two names. At b births there are (b+1)4^b presentations, because
the left/right split is chosen and every letter is free. The complete first
universe therefore has

    6*(1+8+48) = 342 presentations.

Exactly nine distinct strong words have H=41:

    10310, 10321, 11310, 11321, 12011,
    12031, 13101, 31011, 31031.                        (6)

Every one has five chain names, hence seven vertices. This list is of physical
words, not isomorphism classes. Their provenance from the six H=13 words is
retained in `ends_receipt.json`.

From each of the nine words enumerate all end extensions adding zero through
three names. The complete second universe has

    9*(1+8+48+256) = 2817 presentations.                (7)

**None has H=125 or H=127.** No strongness filter is applied at this stage;
all weak and strong outputs are checked. Literal duplicate words are cached,
but their presentation multiplicities remain in the receipt. Every distinct
second-stage count is independently recomputed by full vertex-subset DP as
well as the 36-state transfer. Equations (3)--(7) prove (2) for all orders.

## 3. The failure is branch-specific, not failure of every end fan

The strong parent word 120 has H=13 and strong end-grown children

    12011 with H=41,          3120 with H=43.

The H=43 child does have a further two-ended ternary fan:

    31200 with H=131=3*43+2,
    33120 with H=133=3*43+4.

All are strong, with exact counts checked independently. Therefore the
all-order obstruction identifies the lower H=41 branch specifically. A few
successful fans, even at two successive levels along the upper branch, do
not establish a common-state induction that preserves both branches.

A further positive control is the strong word `321`, with H=13, and its strong end extension `32120`, with H=125. Thus the obstruction is specifically the intermediate H=41 ancestry, not growth from H=13 to H=125 in general.

A potential proof may treat 41 as a fresh seed represented by a different
physical word, or use operations outside (1). The theorem does not show that
no sufficiently large closed grammar exists after finitely many exceptions,
or that either 125 or 127 is globally absent. It proves that this complete
end-growth ancestry through 13 and 41 is impossible.

## 4. Reproduction and source scope

    python3 verify_ends.py
    python3 -O verify_ends.py

The standard-library verifier checks the exact cutoff inequalities, complete
presentation universes, every first-stage count and every distinct second-stage
count by independent subset DP, and the positive upper-branch controls. Its
receipt retains the complete count histogram and middle-word provenance.

The only external bound is the same cited Busch theorem used in the packet:
Arthur H. Busch, EJC 13 (2006), N3, https://doi.org/10.37236/1141. Its full proof
has not been independently re-verified here. No literature-priority claim is
made. The mathematical contribution is the typed all-order reduction and
complete physical-word ancestry obstruction, rather than a larger finite
spectrum prefix.
