# Odd eighths: uniform tail banks and conditional source energy

2026-10-30. **PROVED scoped / FINITE-EXACT; same-agent proof review.**
No independent formal certification or literature-priority claim.
LRC(14), Collatz convergence and fixed-seed H1 remain OPEN.

The four phases1/8,3/8,5/8,7/8 now connect two exactly specified families.
They are the complete maximum-time set of the LRC core{1,...,7,9}, and
the complete source-address set of positive Collatz valuation words of
halving cost three. Character evaluation gives their common phase map;
the directional owners and word clocks remain different required data.

**LRC advance.** A first-arrival inequality handles EVERY finite resonant
tail bank over a core with a proved top-window description. For the entire
family C_q={1,...,q-1,q+1}, a finite residue test completely classifies ALL
additional speeds preserving a strict-branch single-tail maximum. It also
gives the exact terminal mass coefficient using the fastest tied owners.
For example, adjoining16 to C_8 gives maximum2/19; any finite bank of new
speeds preserves this maximum iff every speed avoids0,8,11 modulo19.
Any four distinct new allowed speeds give an infinite thirteen-speed family.

**Collatz advance.** An exact energy factorization works at EVERY dyadic
source depth h and EVERY odd multiplier. The conditional covariance has
bandwidth exactly h, with its outermost entry recovering the entire source
residue. The formula gives exact periodic coefficients for finite remainders,
including, for the3n+1 source-character sums S_M and M>=3,

    E[|S_M+1/2|^2 | R=c mod8]
      =3^-M [1/4+(-1)^M chi_8(c) 3 sqrt(2)/64],

where chi_8 is+1 on1,7 and-1 on3,5. This is an exact conditioned Haar
identity, not a deterministic convergence claim for the integer c.

**Incoming transport classification, extended to odd offsets.** Matching the same words between odd multipliers
q and p preserves every character product through depth D iff
p=q modulo2^(D-1). The first failure occurs at v_2(p-q)+2 and is exactly
an antipodal sign defect. Thus the3-versus5 obstruction appears first at
the odd eighths. Session115 published the same multiplier classification
concurrently at65cec73488 and receives credit; our proof is an independent
replay, with an explicit odd-offset extension.
Changing only the odd offset preserves products, including the plus/minus
reflection, which exchanges the positive and negative integer cones.

**Exact synthesis.** Incoming115's harmonic LRC family has endpoint weights
88/105+(32/105)chi_8(c). This and our energy are two affine readings of the
same reflection-even phase bit. A proved decoder recovers the actual LRC
mass slope from the phase-conditioned energy ensemble, retaining clock parity
and physical multiplicities. The full oriented phase still requires another bit.

- [LRC proofs L1--L3](LRC.md): joint arrivals, all-bank residue classification,
  exact mass, and a vanishing-mass hostile family.
- [Collatz proofs E1--E4](COLLATZ.md): all-depth energy, sharp memory,
  periodic remainders and complete product-transport classification.
- [Recovery and synthesis](RECOVERY.md): exact maps, rejected coincidences,
  current incoming work, and scoped next questions.
- [Phase bridge B1--B2](PHASE_BRIDGE.md): the exact cross-response decoder,
  its multiplicity and clock guards, and its remaining orientation bit.
- [Audit](AUDIT.md), [checker](check.py), [receipt](receipt.json),
  [integrity manifest](manifest.json).

Normal and optimized runs agree on4778 local gates plus seven inherited
81/27-control gates. Finite controls include60 independently constructed
conditional Gram universes,320 joint-tail decisions,1432 individual attachment
decisions,84 literal mass slopes, and all-depth transport checks. General
claims rest on the displayed proofs, not extrapolation of these universes.
