# Physical gap-one cores and a positive-density obstruction

**PROVED in the stated scope; independently audited.** The physical compiler
and its positive-density obstruction have a separate audit in `audit/REPORT.md`.
This does not enlarge the independent deficit-core audit or prove universal
coverage.

The audited sibling packet is `../canonical-deficit-cores/`. It supplies the
physical grammar, core parameterization, and distinction between physical
source profiles and formal integer tuples. The whole spectrum remains open.

## 1. Identify what the one-return construction retains

For any physical source v=(D,E,U,d), write G(v)=D+2U. Then

    v.OJ=(G(v),D,G(v),D).

The OJ construction therefore supplies the balanced-core subfamily
`(D_core,u_core,d_core)=(G,G,D)`.
The synchronized one-return word `C=Z^j OJ` has dagger `JZO^j`.
The chronological identity `JZ=OJ` makes its source `v.OJ.O^j`, so its
joint readouts are

    m=(2^(j+1)+1)G,
    h=2^(j+1)G-D.

This identifies the physical restriction behind the return: its core has
u_core=D_core. Repeating a scalar-ray return is not the same as exploring
general physical cores. The separate return investigation found a prime
obstruction for its repeated-return family; that result is not assumed here.

## 2. A genuine nonbalanced core with an additive offset

In general q=v.J has

    q=(D-E+U,D-E,D+U,D),
    u_core-D_core=E.

Take a unit-deletion source v=(g(t),1,t,1), t odd and t>=3. Such a source
exists for every odd t: begin with cover 2 at the seed, use O to double the
cover and Z to change u to 2u-1. Recursively, the unique zero-J word at cover
t strips O when t is even and Z when t is odd. Its path-count function is

    g(1)=1, g(2t)=g(t),
    g(2t-1)=g(t)+t-(t mod 2), t>=2.

For t odd the last letter is Z and E=d=1. Applying J gives the core

    (D_core,u_core,d_core)=(g(t)+t-1,g(t)+t,g(t)).

It has gap u_core-D_core=1 and is therefore outside the balanced-core slice.
Put R=2^(a+1), M=2^(k+1)-1 for a,k>=0. The exact target compiler becomes

    n=[M(R-1)+2]g(t)+[MR+1]t-1.                    (1)

The physical recipe is the zero-J word at cover t, then J, then O^(a+1),
then Z^k, and finally the closing J. No relaxed-cone tuple is introduced.
The supplied precursor uses one J; the final C1 closing adds another.

At k=0 the formula simplifies to `(R+1)(g(t)+t)-1`. The additive -1 means
these targets are not constrained to be products of the return scalars.
That observation is a mechanism change, not proof of arbitrary target landing.

## 3. Immediate hostile: the family is still too sparse by itself

The separate continuation calculation found an all-order counting barrier.
Its detailed rational bound and finite verification are recorded in
`DENSITY-CALCULATION.md`. In particular, the entire family in (1), with t,a,k
all free, has an upper counting coefficient strictly below 1/2. It therefore
still misses a positive fraction of odd targets. The gap-one family can be
a component of a larger construction, but cannot alone close universal
coverage.

This negative control has its own audit and is not part of the earlier audit
of the sibling core theorem. It prevents treating the new offset as an
automatic positive solution.

## 4. Reproduction and scope

    python3 verify.py > VALIDATION_NORMAL.json
    python3 -O verify.py > VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json

The checker independently compiles the unit-deletion words and verifies (1)
on the stated finite universe and an affine-coordinate basis. Finite formula
checks are not an all-target theorem. The audited deficit-core proof remains
byte-identical apart from its recorded status-only promotion.

## Independent audit and frozen dependencies

`audit/REPORT.md` reviews the all-order physical construction and density
argument. It independently checks 122,500 physical formula choices and 32
literal C1 graphs, and enumerates the complete family through 20,001. These
finite checks corroborate the implementation, not universal coverage.

The original source files are archived in `audit/` and identified by
`audit/INPUTS.json`. `audit/PROMOTION.md` checks the exact status, wording,
and navigation changes for this release. Author verifier and receipts are
unchanged; their pre-audit UNREVIEWED strings are historical run metadata,
not this packet's current review status.

The separately audited core proof, its audit/promotion records, and the
physical adjoint and binary-word dependencies are frozen in `dependencies/`.
The published core source is also available at
`../canonical-deficit-cores/PROOF.md`. This packet does not change its claims.

Reproduce the independent audit from `audit/`:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

From this packet's root, `sha256sum -c MANIFEST.sha256` verifies the complete
release. Both verifiers use only the Python standard library.
