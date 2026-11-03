# Incoming complementary guards: exact identification and new consequences

Received before publication: Session164 checkpoints `88eb65624c` through
`89416555c4`, and Session165 through `d0bda46cae`. The input worktree had been
at `6ef21267ad` when this session began. These are mathematical inputs, not
just merge traffic.

## I1. The user's two-parameter chain is now identified literally

[Session164 COMPLEMENTARY_OFFSETS](../session164_emergent_20261103/complementary_guard/COMPLEMENTARY_OFFSETS.md)
already proves the protected two-sided attachment theorem, its commuting
physical ports, and the exact family giving5->17->53->161. W5 in this packet
is credited to that result, with a canonical prime-quotient representation;
it is not a claim of a new scalar formula.

The incoming seed is `(8,1,19,6,43,15)`. Relabelling by old-vertex order
`(1,0,3,2,4,5)` gives our Q=(2,4,9,19,39,15). Its protected left guard4 stays4
and grows block5 by prepending a vertex; its protected right guard0 becomes1
and grows block0 by appending a vertex. Thus its T(k_left,l_right) is exactly
our T_(l_right+1,k_left+1), with both physical ports transported. Its old
observer port3 becomes our vertex2, not our default vertex1 observer.
`incoming_checks.py` verifies the full relations and these orders, not only H.

What is added here is the four-block enlargement, its dot-product law and
realization guards; the complete fixed-boundary C4/C6 two-ear classification;
the all-size induced-extension obstruction at the resulting161 source;
the actual native485/487 escapes; and the family class/marked recurrence.
The incoming31-seed gives31->97 but its two-sided corner is307, not295.
That distinct +4 obstruction remains intact.

## I2. The native escapes enter a larger infinite continuation bank

[Session164 CANONICAL_RETURNS](../session164_emergent_20261103/profile_algebras/CANONICAL_RETURNS.md)
uses (s,u,d), H=s+u. Temporal CD returns a canonical state t(1,2,0) with
multiplier3, and temporal DCD has multiplier5. The two native escape targets
here have a literal zero-defect guard (after the specified converse for485).
Applying CD once gives the canonical state h(1,2,0).

Consequently for h in {485,487}, every a,b>=0 has a strong literal continuation

    H=3^(a+1)5^b h,    order=13+2a+3b.                    (I1)

This extends the earlier zero-only suffix. The C/D words, moving guard and
vertex clock are retained. It does not turn those multipliers into +2/+4
returns, and it does not assert that arbitrary canonical replacements keep
the original F49 observer unchanged.

The universal rooted guard lift in the same incoming packet agrees with the
prior unrestricted reset H'=H+2F. Its correction J=F-H is essential: access
to a canonical bank changes H by a known slack, not by a free source reset.

## I3. Two actual probes and a finite polynomial test bank for a wing

[Session165 guard geometry](../session165_geometry_representation_20261103/guard_geometry/NOTE.md)
separates positive physical words from their full matrix algebra and proves
a bounded polynomial identity test using actual guarded sources. The wing
family supplies a lower-dimensional instance of the same proof mechanism.

Pair a wing (x,b) first with the actual wing (3,0), then with (6,1), and let
their full tournament counts be h0,h1. Equation W2 gives

    x=(h0+1)/6,   b=(h1-2h0-1)/2.                         (I2)

Both probes are actual strong substitutions. They recover the integer joint
response, including its realization guard; modulo7 the divisions remain
valid. They do not recover the hidden clock at the ramified wing (5,3): the
extra legal growth query in W4 is still needed.

For an **integral linear** observer valid modulo every modulus, add the actual
partner (5,0). Normalize each odd count by K=(H+1)/2, and order the three
partners as (3,0),(5,0),(6,1). Their observation rows are

    [[3,0],[5,0],[6,1]],
    x=2K0-K1,   b=K2-6(2K0-K1).                         (I2a)

The maximal minors include3 and5, so their gcd is1 and the Smith factors are
(1,1). Three is minimal among such partner probes for an integer left inverse
on arbitrary two-coordinate response weights: a wing with l=1 has b=0 and
odd x>=3; with l>=2 it has even x and odd b. Two rows of the latter type have
even determinant, while a mixed pair has determinant of magnitude at least3;
two rows of the former type have determinant zero. None has determinant1.
This minimality concerns integral linear recovery on the ambient module,
not nonlinear recovery on the physically realizable dyadic subset.

Halving is an exact normalization of known odd counts. Raw H modulo2M suffices
to form K modulo M; if M is odd, H modulo M already suffices because2 is a
unit. No precision is silently discarded. For example the first
two partners of I2, with determinant12 before normalization, fail modulo3
even on actual wings (3,0) and(5,0); I2a removes that odd-prime defect.

If a rational polynomial P(x,b) of total degree at most delta is proposed as
an identity on all wings, it suffices to test the actual grid

    1<=m,l<=delta+1.                                     (I3)

Proof. Substitution x=M(L-1)-L/2+2, b=L/2-1 gives a polynomial of degree at
most delta separately in M,L. This substitution is dominant over Q: outside
L=1 it has the displayed rational inverse. Therefore a nonzero P remains
nonzero. Two univariate root bounds on the distinct dyadic grid establish
I3. The largest h1 probe has at most2delta+7 vertices. This is a certified
test bank only after polynomial dependence and its degree are proved; it
does not make arbitrary response fitting reliable or solve integer coverage.

This is the useful synthesis: old quotient walks establish the finite response
space, actual probes read its coordinates, and realization guards determine
which algebraic operations can be performed on sources. Dropping any one of
those three steps would misidentify scalar success with reusable closure.

## I4. Final incoming scope: regular-five obstruction and common CRT witnesses

Incoming through `2af5e7549d` adds the
[six-root return forest](../triple_reentry_interfaces_20261103/RETURN_FOREST.md)
and the [strong induced-R5 readiness obstruction](../lambda_emergent_20261103/READY_REPAIR.md).
Their scope is preserved: describing all responses is different from finding
one physical interface meeting all requirements; an R5-containing source
cannot acquire a one/two-ready root in a **strong** induced extension. Neither
statement proves universal exact ternary closure.

Our whole four-module family is induced-R5-free. Indeed Q has two vertices
of outdegree1, so every five-vertex subtournament contains a vertex of
outdegree at most1 and cannot be regular. The regular five-vertex tournament
is prime (directly checking its proper subsets suffices); any induced copy
in a transitive substitution would therefore use at most one vertex per
block and project to a forbidden copy in Q. The two literal escape targets
485 and487 are also R5-free: all462 five-subsets of each have been checked.
Thus the new escapes avoid this particular inherited obstruction. Its absence
is not a readiness certificate, as the incoming work explicitly warns.

Reproduce the latter finite check from the repository root:

```python
import sys, json
from pathlib import Path
from itertools import combinations
sys.path.insert(0, '05-knowledge/results/defect_offset_dynamics_20261103')
from wing import induced
p = Path('05-knowledge/results/defect_offset_dynamics_20261103/native_escape_161.json')
for h, row in json.loads(p.read_text()).items():
    T = tuple(row['target'])
    print(h, sum(all(r.bit_count() == 2 for r in induced(T, S))
                 for S in combinations(range(len(T)), 5)))
# 487 0; 485 0
```

The new common-CRT-witness obstruction also clarifies I2a. Coprime minors
license reconstruction by combining several actual observations; they do
not assert that one of those probes is a common unit for every prime factor
of a modulus. Joint interface identity and observation-module invertibility
remain separate requirements. The incoming P4 path-index correction is
accepted; none of our arguments uses the historical incomplete P4 proof.
