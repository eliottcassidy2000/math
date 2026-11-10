---
id: THM-5114
title: "Scalar countermodels with constant-error anchoring and the inherited compact gap-image loss"
status: "PROVED, independently audited. Inherited image bound; new abstract logical nonimplication. Not spectrum exclusions."
session: Session244-even-gap-anchoring-20261110
source: 05-knowledge/results/session244_even_gap_anchoring_20261110/selector/proof.md
audit: 05-knowledge/results/session244_even_gap_anchoring_20261110/audit/REPORT.md
---

# The missing absolute-count coordinate is substantive

The density estimate below is INHERITED from Session199
`session199/hidden-structure/PROOF.md`, section4, for the same balanced image.
Its application identifies precisely the loss in the Session243 gap compiler;
no new density theorem is claimed.

The compact THM-5112 pair has lower L(N)=g(N)+2N and upper R(N)=g(N)+4N,
with g(1)=1. Since g(2N)=g(N), R(N)=L(2N): all upper endpoints are redundant
as full-count values. The image S={L(N):N>=1} satisfies

    limsup |S intersect[1,X]|/X <= c,
    c=sum_{j>=0} 1/(2^(j+2)+1)=0.431166447015110... .

Therefore the lower limiting fraction of odd targets omitted by this
particular compiler is at least1-2c=0.137667105969779... . Partition
N=2^j m with m odd and use g(m)>=(m+1)/2; summing the corresponding input
bounds proves the estimate, with an O(log X) rounding error. This does not
exclude any target from the full tournament spectrum.

There is a stronger abstract obstruction to scalar-only reasoning. The set
of positive odds excluding7,21 and every prime>=151 is multiplicatively
closed, has arbitrarily long consecutive odd intervals, and has infinitely
many pairs for every even gap. For each N>=1 it has a pair(a,a+2N) with
2N+1<=a and a+2N<=5N. For every odd anchor X>=51 and every N>=1 it has such
a gap pair with X<=a<=X+28, without the compact bound being asserted for
that separately chosen pair.

For N>=30, choose odd a in[2N+1,3N] with 3|a and5|(a+2N), one residue
modulo30. Both endpoints are composite and outside7,21. Small N use
verified compact witnesses below151. The same CRT progression gives the
constant-error and infinite-fibre claims. Long intervals follow from CRT
forcing consecutive odd entries to have specified small divisors.

Moreover, for any prescribed finite prefix and any prime p>=151 beyond it,
the set of positive odds excluding only7,21,p has all these properties and
contains that allowed prefix. Thus finite evidence, multiplicativity,
compact every-gap existence and even uniform bounded-error placement do not
logically force exact coverage. These are abstract scalar sets, not claimed
tournament spectra; they do not satisfy every additional physical identity.

[Image proof and exact inverses](../../05-knowledge/results/session244_even_gap_anchoring_20261110/selector/proof.md) ·
[Countermodels](../../05-knowledge/results/session244_even_gap_anchoring_20261110/logic/PROOF.md) ·
[Finite-prefix strengthening](../../05-knowledge/results/session244_even_gap_anchoring_20261110/logic/ROUND2.md).
