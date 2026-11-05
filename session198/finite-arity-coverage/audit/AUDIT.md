# Independent audit: finite-arity physical coverage reduction

Verdict: PASS. The cofinite range lemma remains an explicit open hypothesis.

- PROOF.md SHA256: a30dee0615c0332236e0f31ddb15cac5522fa1d2fdd921e67f96cda3d5e86657
- verify.py SHA256: efe82f416656a5f7405309839d9961d367f5569d8049fd8627a0635cbd0b2814
- Both receipts SHA256: 8a182bd774a806aa2b84b75ee3a64348c4b8d2445b4bd3dfa7749882828397cf
- REPAIRS.json SHA256: efb96a9d833112f7ec2897a00926a4451b96af1f4d6bc245dc510a898eeb1505

The absorption identity holds in all four physical coordinates. For r>=1,
O^r resets E and scales U by2^r; the even index2^r L resets exactly the same
endpoint term, and the g/Q/B recurrences give the remaining entries, including
the new deletion count D. The singleton and L=1 boundaries are correct. This
normalizes every source with at most four indexed gadgets and arbitrary sink
runs. It does not assert a bounded-gadget grammar for all sources.

The two arithmetic families are genuinely physical for every allowed tuple.
Projected T_L updates do not depend on the omitted deletion count. The final
power-of-two map gives exactly the C1 closure readout with the specified sink
run. The sole value1 boundary is excluded; every remaining output is strong
C1. The one-apex converse profiles w(L), L>=2, are actual nontrivial sources,
so both bilinear forms in the two-argument family are strong C2 witnesses.

The logarithmic graph-order proof is valid. Each indexed gadget multiplies U
by at least its index, and the terminal target bounds2^r times the resulting U.
Thus2^r product(L_i)<=n. The exact vertex count is
k+2+r+sum ceil(log2 L_i), at most floor(log2 n)+10 for k<=4.
The C2 product bound and the six literal repairs also obey the same looser
bound. This is unconditional for represented tuples; existence of a tuple for
every large n remains conditional.

The bounded enumerator is complete for exactly these families. F>=U and strict
F growth bound all relevant prefixes and indices. Breadth-first insertion makes
the retained depth minimal, so merging projected triples cannot lose any
continuation with the four-gadget budget. All terminal exponents are tested.
Symmetry and3LK+1 bound the complete C2 pair loop. The independently evaluated
gadget constants and every selected positive recipe are checked, not just the
presence of a target bit.

Isolated normal and optimized runs reproduce the frozen receipts and all review
manifest hashes. They certify49,992 positive tuples through100,001, including
49,956 from the restricted C1 family, with the stated eight missing values.
All six repair adjacency records are unchanged from the previously audited
C1-monotone packet and independently count correctly, leaving only7 and21.

Additional controls:
- I inspected, copied and recompiled the independent C++ enumerator from the
  separate adversarial lane. It marks terminal fourth-gadget states without
  the Python verifier's final deduplication/parent-record machinery. Its exact
  full100,001 support, C1 count and142,907 pair tests agree. The normalized
  C1 gap list through50,001 also agrees with the earlier unnormalized probe.
- A separate one-vertex physical compiler and forward subset counter verify
  96 literal absorption pairs, including singleton and nonunit deletion inputs,
  through graph order12. All four profiles, starting counts and C1 counts agree.
- The same independent counter verifies all six repair graphs and order bounds.

Reproduction artifacts are ../finite-arity-control.cpp,
../FINITE_ARITY_RECOMPILED.json, ../finite-arity-physical.py and
../FINITE_ARITY_PHYSICAL.json. The copied C++ control is attributable to the
independent adversarial worker, not newly authored in this review.

No extrapolation to the unproved threshold lemma is justified. If that lemma
is later proved, the conditional all-target and logarithmic-order conclusions
follow; this packet alone proves only its all-order reductions and exact bounded
certificate. No unrestricted nonexistence theorem for7 or21 is claimed.
