# Return banks: a bulk valuation obstruction and an actual channel switch

**Status: PROVED scoped return-bank obstruction and actual infinite switched
replacement family; FINITE-EXACT small matrix/source audits; unrestricted
spectrum and a closed ternary constructor remain OPEN.** No external novelty
or literature-priority claim is made. Input commit is
`ff94a6d74130612bf3e34f2ac278c1f848168d88`.

The new point is that a constant-error eigenspace permits a transported affine
law but can prevent the return required for induction. For the inherited
`{b,q2}` and `{r21,r38}` channels, a ternary replacement raises an exact
3-adic bulk level; every allowed return word preserves that level. Thus a
finite bank of fixed seed interfaces of the explicit form below cannot close.
This is an all-length return theorem, not an inference from a modular scan.

An actual escape from that particular invariant also exists: the q5 channel
can raise the bulk level, and a singleton changes its error line into the
b/q2 line. The switch supplies a new explicit infinite family, but does not
provide a return cycle or an all-target construction.

## 1. Inheritance and the concept board

[INHERITANCE.md](INHERITANCE.md) gives the source maps and scope distinctions.
Closest positive mechanism: the actual seed replacements in
[golden inverse TRANSPORT T1-T6](../golden_inverse_defects_20261103/TRANSPORT.md)
and [inverse65 INCOMING N2-N3](../inverse65_zeckendorf_20261103/INCOMING.md).
Canonical hostile: the unrealizable next states `(287,122,225)`,
`(289,121,225)`, `(293,119,225)`, `(295,118,225)`, although their scalar counts
are attainable elsewhere. Corrected near miss: the finite affine certificate
needs an actually closed parent language and the corrected negative-offset
descent bound. The least-used sidecar here is the integer lattice in the
two-dimensional quotient by the error line.

| Live concept | Object and operation | Lost coordinate / cheap hostile |
|---|---|---|
| Actual ternary induction, anchor | Replace a named source seed, retain its continuation | Test whether the output returns to an input interface |
| Error versus charge | Primitive left error tau and dual charge ell | A rational projector need not preserve the source lattice |
| Bulk quotient, niche | Remove the error line, retain a saturated integer quotient | Test prime divisibility of quotient determinants |
| Channel switching, wildcard | Apply a bank symbol to an error direction | Compute every one-letter map between the five known lines |
| Longer return words | Products of actual transfer matrices | A two-letter negative test has no all-length consequence |
| Scalar coverage | Repaired multi-interface affine certificate | Each inverse must land in an independently proved input language |

The perspective generator was `actual sources × error/bulk coordinates ×
valuation × seed replacement/return × charge-line quotient × unbounded word
length`. It produced the following three rounds rather than a new large
support scan.

## 2. Round one: quotient the correction before looking for returns

Write states as row vectors `v=(h,x,y)` in **processing order**. A letter acts
by `v -> vP`; literal source-component order is reversed. Suppose primitive
integer row tau obeys `tau P = sigma tau`, with sigma=+1 or-1. Choose an
integer 3-by-2 matrix B whose columns form a basis of the lattice annihilating
tau. Then

`P B = B Q`

for an integer 2-by-2 matrix Q. The bulk of an actual state is `beta=vB`.
Different choices of this saturated basis differ by GL(2,Z), so

`level_p(beta)=min(v_p(beta_1),v_p(beta_2))`

is basis-independent for nonzero beta. Zeros among the coordinates are
allowed. In these channels beta cannot vanish for an actual state: the kernel
is the rational line through tau, whose mixed signs exclude `h>=1,x,y>=0`.

The exact quotient matrices are:

| Error line | Integer bulk beta | Letter | Q | det Q |
|---|---|---|---|---:|
| `(2,-1,0)` | `(y, h+2x-y)` | b | `[[5,8],[2,3]]` | -1 |
| same | same | q2 | `[[11,18],[4,7]]` | 5 |
| `(3,-1,-1)` | `(x-y, h+3x)` | r21,r38 | `[[3,-49],[-2,37]]` | 13 |
| `(2,-2,1)` | `(h+x, x+2y)` | q5 | `[[9,14],[6,7]]` | -21 |

Each displayed B has a two-row minor of determinant +/-1, proving
saturation without a rational-coordinate assumption. The probe checks PB=BQ
entry by entry. The first quotient is the inherited `(Y,X2)` bulk; the
norm13 quotient is another integral presentation of that inherited channel.

**Valuation lemma.** If det Q is prime to p, then
`level_p(beta Q)=level_p(beta)`. Divide beta by its maximal common power of
p. Its reduction is a nonzero row, and multiplication by the invertible
matrix Q modulo p leaves it nonzero. This proves equality; it does not merely
bound the level.

Thus every word over `{b,q2}` and every word over `{r21,r38}` preserves the
3-adic bulk level, irrespective of length or sign phase.

## 3. A precise finite return-bank theorem

This theorem concerns the following deliberately explicit source-replacement
model. It does not quantify over every possible finite interface or algorithm.

Fix one correction line tau and one transport alphabet with the quotient
matrices above. Interface i has a fixed actual seed word U_i, state v_i and
nonzero bulk beta_i. Its witnesses are words `U_i w`, with w in a declared
compatible continuation language. An edge i->j realizes a ternary source
replacement and returns to interface j by inserting a **fixed transport word**
r before the retained continuation, with exact full-state certificate

`v_j P_r = 3 v_i + m_e tau`.                         (R1)

Language compatibility must also be supplied: the new word `r w` must belong
to interface j's language. The integer m_e, terminal observation and sign
guard determine the scalar offset. This model includes the natural attempted
repair of the inherited seed replacement: reexpress its output as another
fixed seed followed by a transport reset, while retaining the old continuation.
It permits changing the actual source graph; it is broader in that direction
than an induced-vertex extension or the old append-only scalar-chart model.

**Theorem R.** For either transport alphabet `{b,q2}` or `{r21,r38}`, every
edge satisfying (R1) obeys

`level_3(beta_j) = level_3(beta_i) + 1`.              (R2)

Consequently no finite such bank has a directed cycle. From a seed at level
l, at most `max_i level_3(beta_i)-l` ternary replacements are possible. In
particular such a finite bank cannot implement arbitrarily many ternary
digits, even when its resets have arbitrarily large but fixed word lengths.

**Proof.** Multiply (R1) by B. The error disappears and
`beta_j Q_r=3 beta_i`. Every factor of Q_r has determinant prime to3, so the
valuation lemma identifies the left-hand level with `level_3(beta_j)`.
Multiplication by3 raises the nonzero right-hand level by exactly1. A cycle
would make an integer level equal to itself plus its positive length. The
depth bound follows by summing (R2). No computational word cutoff is used.

The same proof works at any prime p with nonzero integer slopes a_e in place of3: edges obey
`level_p(beta_j)=level_p(beta_i)+v_p(a_e)`. A cycle containing an edge whose
slope is divisible by p is impossible when all quotient returns are p-units.
Zero-increment edges need not be excluded by this argument.

For the inherited common seed and its two outputs, the bulk values are

| State | beta | level_3 |
|---|---|---:|
| `(31,14,25)` | `(25,34)` | 0 |
| `(95,41,75)` | `(75,102)` | 1 |
| `(97,40,75)` | `(75,102)` | 1 |

Hence no b/q2 reset word can turn either output into the original level-zero
seed interface. Adding finitely many higher-level seeds only postpones this
failure. The identical argument applies before the terminal b in the
inherited norm13 `3H-10` and `3H-20` constructions: their preterminal
replacement errors lie in tau=`(3,-1,-1)`.

**Algorithmic consequence.** Given a proposed bank in this model, compute its
seed bulk levels and check (R1) exactly on every advertised edge. The levels
are a finite, auditable obstruction certificate. A necessary way out is a
return whose quotient determinant is divisible by3, a change of correction
line, an unbounded seed bank, or a source-dependent reencoding outside (R1).
These are alternatives to investigate, not sufficient conditions for success.

## 4. Round two: a channel that can pay the level cost

The q5 quotient has determinant -21, so the preceding unit argument does not
apply. This is not just a missing proof hypothesis. The actual seed q2 has
state `(5,3,4)` and q5-bulk `(8,11)` of level0. One q5 gives
state `(87,51,69)` and bulk `(138,189)` of level1. A second q5 gives bulk
`(2376,3255)`, still level1. Thus q5 really can absorb one factor of3 in an
actual bulk, but it does not simply implement repeated scalar tripling.

The procedural next probe was the complete graph of one-letter maps between
the five primitive signed-error lines recovered from the seven charged
symbols. There are exactly ten maps, all with multiplier +/-1:

| Source error line | Letter(s) | Target line | Sign |
|---|---|---|---:|
| bq2 | b | bq2 | -1 |
| bq2 | q2 | bq2 | +1 |
| q1 | q1 | q1 | -1 |
| q4 | q4 | q4 | -1 |
| q5 | q5 | q5 | -1 |
| q5 | a | bq2 | -1 |
| q5 | r21,r38 | bq2 | +1 |
| r21 | r21,r38 | r21 | +1 |

Here row labels mean the error vectors in section2 and
q1=`(3,-2,1)`, q4=`(3,-2,0)`. The classification is **FINITE-EXACT** on
these five lines and65 supplied matrices. It says nothing about other lines
appearing inside a longer word. The q5->bq2 bridge is one-way in this graph;
there is no reverse one-letter edge.

For perspective, a separate exact audit of all4225 two-letter words finds
only eleven products with eigenvalue+1 or-1: bb, bq2, q2b, q2q2, q1q1,
q4q4, q5q5 and the four words in `{r21,r38}^2`. Thus this short probe finds no
new mixed return channel. Longer products remain OPEN; this finite result
is not used to prove Theorem R.

## 5. Round three: make the channel switch an actual construction

Take processing seeds

`U=q2`, `v(U)=(5,3,4)`;
`V=a a a a a a`, `v(V)=(17,7,13)`.

They satisfy `v(V)=3v(U)+tau5`, where tau5=`(2,-2,1)`. Exact inherited
matrices give

`tau5 P_q5=-tau5`, `tau5 P_a=-tau2`,
`tau2 P_b=-tau2`, `tau2 P_q2=tau2`,

where tau2=`(2,-1,0)`. Therefore for every k>=0 and every word z over
`{b,q2}` the actual words

`L = q2 q5^k a z`,
`R = a^6 q5^k a z`

satisfy the all-parameter identity

`H(R)=3 H(L)+2(-1)^(k+1+#b(z))`.                    (R3)

This proof multiplies a difference of two actual states; it never declares
the mixed-sign error vector itself to be a tournament. Both words have
multiple source components and at least five vertices, so their retained
Hamilton spine closes through the complemented endpoint arc to a directed
Hamilton cycle. They are strong. Their orders are respectively
`5+4k+cost(z)` and `7+4k+cost(z)`, where b costs3 and q2 costs4.

Independent recursive path counting on the reconstructed Boolean graphs
checks the following controls; the first right-hand graph is also counted
by all7! literal permutations:

| k | z | H(L) | H(R) | Offset |
|---:|---|---:|---:|---:|
| 0 | empty | 11 | 31 | -2 |
| 0 | b | 59 | 179 | +2 |
| 1 | empty | 189 | 569 | +2 |

The sign bit matters. Dropping it would turn the first row into the false
claim31=3*11+2. The switch has not returned to a q5 input interface, supplied
the third ternary digit, or closed the two positive-offset branches. Its
value is an actual source-level operation connecting the singular quotient
to the older constant-error channel.

A small discovery saturation of all reachable states H<=500 has3592 states.
Requiring both source and output states inside this bound finds25 bq2,
5 q1,4 q4,8 q5 and3 r21 source states with at least one of the tested lifts.
The tested m values are1,2 for even first coordinate of tau and2,4 otherwise.
Only the bq2 channel has a jointly lifted source in this bounded test. These
counts classify neither all sources H<=500 with unrestricted outputs nor
all future channel words. They were used to select (R3), not as coverage
evidence. All source and target words are in the receipt.

## 6. What a genuine next induction must provide

The positive signal and the obstruction now meet at a small proof interface:

1. Construct an actual return into a q5-type or another 3-singular quotient
   channel after a b/q2-type branch. A scalar equality does not identify the
   output source. Retain the full state, source word, sign and quotient level.
2. Supply all needed residue branches and prove that every output is again
   in the input language. A finite list of seed pairs with formulas valid
   under transport is insufficient.
3. Only then apply the repaired multi-interface affine certificate, with
   all seeds up to its threshold and all inverse descendants covered. If
   reducing only to the prime and7p lanes, arbitrary affine inverses still
   require a larger closed parent language.

A different operation, including re-marking, both-ended joins, or the
adaptive cap3 matrix-state family, can leave (R1). Neither Theorem R nor the
finite scalar-affine spectral-radius test rules out such a family. Conversely,
observed modular surjectivity, eigenvalue1, or a rational projector does not
certify its source realization or closure.

## Reproduction and audit scope

Run `python3 return_bank_probe.py` and `python3 -O return_bank_probe.py` from
this directory. The normal and optimized receipts agree byte for byte:

`return_bank.json` SHA256
`65ed181a49907876a3bd84fbbc96f5ca8193e8215fc15e18416e4cdbc9cef803`.

The script loads only the pinned atlas, independently reconstructs literal
graphs, and uses explicit failure checks that survive -O. Its finite universe
is stated in the receipt. The general return obstruction and switched-family
formula are algebraic proofs above. No external literature was needed or
consulted, and none of the inherited large censuses was rerun.
