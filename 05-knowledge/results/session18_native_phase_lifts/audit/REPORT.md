# Independent audit: native affine phase lifts

**Verdict: PASS in the stated scope. No blocking proof or implementation defect found.**

This review covers the frozen packet with manifest SHA-256
`f20cb439a11829cd1ba0529d06f5e598f5fa496faf19ac930697b9ad88850ffa`
and receipt SHA-256
`d1c9ce40b69d4952dee72b2f10b8b39ba7d960eaa6df8de5ba8e84d18a059de6`.
The producer's files were read without modification. All audit work is separate.
No theorem identifier or priority claim is assigned.

## 1. Scope and inheritance

The proved content is the exact native affine lift, its sparse unit-defect
phase law, guarded synchronous zero packets, the source-disjoint same-scale
compiler, and the residual external cross-term identity. The computational
receipts are finite exact controls, not a substitute for the general proofs.

The two inverse chains are correctly credited to THM-4690, statement 2,
`01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md`.
That theorem itself credits THM-4667(8). The period-three affine notation
does not establish a new inverse classification. The fixed-source 64-scale
identity and 153N source-norm budget are inherited from
`05-knowledge/results/native_cancellation_105_20261018g.md`, section 3.

All four manuscript source hashes in the producer manifest matched the
actual repository bytes. Its frozen annulus manifest hash also matched.
The frontier and ranks agree exactly with the full inherited JSON, whose
SHA-256 is `d2d4d152593a03eb653d1fd6fecfbc13f573a6eb74e1310f59b7a6cb78b01721`.
The correction boundary in MISTAKE-631 was read: these are native strip
clocks, not shortcut stopping-time differences at ROOT, and inherited
preimages are not claimed as novel. MISTAKE-634's unrelated geometric
correction does not alter the inverse formula used here.

## 2. Algebra and boundary review

- For the actual branch a in {1,2}, the constants c=21+7*2^a give
  `(3(64s+c)+1)/2^a = B^3((3s+1)/2^a)` exactly. The lifted residue
  remains in the same branch. Three additional B-strips therefore give
  the same primitive P-target, strip count r+3, and integer clock t+6.
- Admission is an integer iff, not an approximate cutoff. If s<4N then
  s<=4N-1 and 64s+c<256N; if an odd s is not below the even boundary
  4N, then s>=4N+1 and its lift exceeds the enlarged cutoff.
- The phase rotation follows from B^j(lift(s)) modulo 3. Both palettes
  rotate under the same Eisenstein unit, and the difference of the two
  threshold indicator representations gives precisely the stated sign
  and coefficients in the defect formula.
- At a fixed j and K, the source crossing interval has width
  `(c-21)/64+21/(64*4^j)<1`. The two branch intervals are nested, so
  there is at most one integer source across both branches, not one
  per branch. The j=0 endpoint is integral and cannot cross.
- If B^j(s)<4N, integrality yields 4^(j+1)<=12N, hence j<h_N.
  Thus the candidate-event bound is exactly 2(h_N-1), including N=1.
- The strict comparison `m < B^j(lift(s))/64 < 2m` confines every
  crossing to N<m<4N. The next B-chain point is already greater than
  4N. At the unique possible point, crossing the lower cutoff requires
  m<2N while crossing the upper requires m>2N. There is at most one
  event per source, so no hidden sum of multiple units occurs.
- The literal N=2019, s=1009 hostile has old charge (-1,-2), rotated
  charge (2,1), actual lifted charge (1,0), and defect (-1,-1), norm 1.
  Omitting the phase guard leaves the triple residual (0,-1), norm 1.

## 3. Stopping sets, horizons, and no reuse

The exact frontier has minimum excluded source 1009 and maximum 63665.
The bound `64*1009+35>63665` is sufficient to keep every lift of every
outside source outside V. It is not inferred merely from V being finite.
Since the sources then all take their first edge into the same target,
their actual terminal and P-depth agree whenever they arrive, and the
clock differences are 6 and 12. Non-entry is shared too; no universal
entry premise is introduced.

The `{1,113,7281}` hostile was independently reconstructed: first arrivals
are `(terminal,clock,depth)=(1,0,0),(113,0,0),(1,14,1)`.
Their unlabelled phase sum is zero but their stopped terminal energy is 6.
Deleting this false cross-terminal packet would be invalid.

For valid packets, every finite P-depth arrival cut sees all three sources
or none. This is not true of integer-clock cuts. The two named packets
have nonzero norm-1 partial vectors at their first clock thresholds.
The polynomial factorization is correct as a polynomial over the
Eisenstein coefficient ring; it vanishes at z=1 and z=-1 but is not the
zero polynomial for a nonzero charge. At native depth zero, the sources
still occupy different physical coordinates. The one-edge cancellation
does not remove their initial full-state energy.

At fixed M=64N, each admitted adjacent lift has its predecessor in S_N.
The fixed-source transport theorem therefore converts the cross-scale
phase law into the required same-scale relation. Splitting chains at V
members and bad edges, then grouping consecutive triples, is valid.
The unique branch and strip count prevent chain overlap; block grouping
prevents source reuse. Short segments and all residual coherence remain.

### Additional hostile: the second adjacent guard is independently needed

At N=129251 and M=8272064, take the same sources
`1009,64625,4136049`. The first edge is good; the second is bad:

`B(64625)=258501<258502=2N`,

`B(4136049)=16544197>=16544128=128N`.

The three actual same-scale charges are
`(-1,-2),(2,1),(0,1)`, whose sum is `(1,0)`, not zero.
The compiler correctly leaves all three sources unspent. This complements
the original N=2019 first-edge hostile: neither adjacent guard is optional.

## 4. External interactions and the final scope fence

With beta(v,w)=Q(v+w)-Q(v)-Q(w), rotation gives
Q(v+w)=Q(v) and beta(v,w)=-Q(v). Consequently the pair's added diagonal
is exactly paid by its internal cross term. Expansion with an arbitrary
remaining terminal vector Z leaves exactly beta(Z,w), as stated.
That term can have either sign. Between-pair interactions are not removed
by independently paying each internal pair term.

Complete source-disjoint triples cancel their cross terms with all other
sources in the same terminal. But a triple mixing old and new sources
still changes the old baseline in an annulus difference. The correction
in the frozen annulus packet's equation (D) remains necessary.

This work proves neither a uniform 153N annulus-increment inequality nor
universal entry. It gives no lower bound on the fraction of the annulus
budget cancelled. It does not bound the full moving-endpoint Gram energy
or erase partial integer-clock packets. These distinctions are explicit
in the frozen manuscript and are preserved by this verdict.

## 5. Independent exact checks and reproduction

The producer replay was run from a copied standalone directory, normally
and under Python -O. Both matched its saved receipt at 945113 checks.

The separate `independent_audit.py` imports no producer or repository code.
It uses all three original real readouts for charges, a valuation formula
instead of a strip loop for primitive edges, closed-form affine iteration
for inverse reconstruction, and separate integer lower/upper bounds for
each branch's crossing interval. It also implements and tests the compiler,
which the producer manuscript states but does not itself enumerate.

Normal and optimized independent runs match at **6801996 explicit checks**:

- all 197120 primitive source rows for every N=1,...,512, including zero
  charges, fixed-source transport, both admission directions, exact lift
  clocks, and complete sparse exception equality;
- independent reproduction of both saved per-source digests;
- all 98304 primitive sources below 262144 reconstructed by inverse address;
- all primitive targets below 16384, both branches, and k=0,...,4 checked
  against the target-side inverse formulas;
- all saved sparse controls, including N=10^40+7;
- actual named routes, full clocks, registers, both original hostiles,
  the added second-edge hostile, and integer-clock-prefix hostiles;
- compiler checks at eleven declared N values, with 5 nonzero disjoint
  triples at M=1048576, 982 at M=4194304 and M=4194752, and 2444 at
  M=8272064; these are finite controls, not density estimates;
- pair payment, clock-polynomial coefficients and residual beta(Z,w)
  over an integer vector box larger than the source alphabet.

From this audit directory, with the frozen producer files in PACKET:

```text
python3 -B PACKET/check_phase_lifts.py
python3 -B -O PACKET/check_phase_lifts.py
python3 -B independent_audit.py PACKET
python3 -B -O independent_audit.py PACKET
```

The `portable/` directory is only the audit's local immutable replay copy.
The publication deliverables are this report, the independent checker,
the two independent JSON receipts, two producer replay logs, provenance
receipt, and `AUDIT_MANIFEST.json`. No repository or producer bytes were
changed by the audit.
