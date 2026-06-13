# HYP-2446 - Grothendieck-Katz p-curvature gives the LRC14 denominator-curvature ledger

**Status:** OPEN synthesis; toy operator/scalar mismatch verified.
**Source:** codex-2026-06-12.
**Companions:** HYP-2443, HYP-2444, HYP-2445, HYP-2438, THM-492,
HYP-2241, HYP-2256.
**Computation:** `04-computation/p_curvature_lrc14_support_gate_codex.py`;
stored output `05-knowledge/results/p_curvature_lrc14_support_gate_codex.out`.
**External anchors:** Shankar arXiv:1610.05674; Tang arXiv:1412.7875;
Lam-Litt arXiv:2501.13175; Katz 1972.

## Statement

The Grothendieck-Katz p-curvature conjecture should be imported into the repo
as the following proof schema:

```text
local tests at almost all primes are not scalar residues;
they are operator/carry tests.

vanishing of the operator obstruction at almost all primes
forces finite global monodromy.
```

For LRC14 the proposed analogue is:

```text
raw q-blocking is not the proof object;
the proof object is a denominator-curvature ledger

K_q(S) = (blocked unit twists,
          marked blocker support tau_q,
          shell-27/Pisano class,
          13-clock debt,
          divisor fiber,
          Bprime/owner-private target).
```

If a primitive LRC14 speed set blocks a long fibered ladder, then these
local ledgers must be compatible under denominator lifts such as

```text
q -> p*q,        q -> d*m with d|14,
shell 27 -> shell 9 -> shell 3,
plain shell -> Q27 fiber.
```

The conjectural support-pressure theorem becomes:

```text
Either a finite denominator witness appears,
or the compatible curvature ledger descends to a finite wall type
AP, Vstar, 2AP,
or the incompatibility opens Bprime(any runner) / owner-private deletion.
```

In short: LRC14 needs a p-curvature-style "operator" for denominator covers.
HYP-2443 supplies the marked support part; HYP-2444 supplies the Pisano class
part; HYP-2445 supplies the Frobenius support-gate model from geometry.

## Toy Evidence

The computation works on local p-jets

```text
R = F_p[z]/(z^p)
```

for rank-one connections `D + M_a`, comparing the true local operator

```text
(D + M_a)^p
```

against the deliberately lossy scalar multiplication by `a^p`.

For primes `3,5,7,11,13`:

```text
a = 1/(1-z):
  p-curvature ranks = [0,0,0,0,0]
  naive a^p ranks   = [3,5,7,11,13]

a = z/(1-z):
  p-curvature ranks = [3,5,7,11,13]
  naive a^p ranks   = [0,0,0,0,0]
```

Thus the scalar shadow can be wrong in both directions:

- it can see a full-rank obstruction where the operator cancels it;
- it can see no scalar obstruction where the derivative/carry term creates one.

This is the exact warning for LRC14.  The scalar statement "denominator `q`
is blocked" has to be replaced by a marked, functorial ledger of how it is
blocked.

## LRC14 Evidence

The same run records the current one-stranger denominator table:

```text
row       plain<=27  plain<=41  Q27       Q41
AP        none       none       none      none
Vstar     none       none       none      none
S611      none       9/41       1/91      9/41
S702      none       9/41       1/91      9/41
S793      none       11/40      11/40     11/40
S1053     none       11/40      11/40     11/40
```

The interpretation matches HYP-2444:

- AP and Vstar are finite-monodromy-like wall atoms, not failures of the
  finite search.
- `S611` and `S702` look blocked to the plain shell `q<=27`, but the fibered
  operator ledger catches them at `q=91=7*13`.
- `S793` and `S1053` are caught at `q=40`, already inside `Q27` because
  `40=2*20`.

This is precisely the p-curvature pattern: the naive local face can miss the
actual operator witness until the retained lift/fiber coordinate is attached.

## Proposed Lemmas

1. **Curvature-ledger definition.**  Make `K_q(S)` functorial under the
   denominator maps used in THM-492 and HYP-2444.  The maps should track not
   only whether a twist is blocked, but which runners pay for it.
2. **Compatibility implies wall.**  If the ledgers stay compatible through
   the full LRC14 fibered ladder without producing a witness or Bprime route,
   then the speed set lies in the finite wall family AP/Vstar/2AP after
   normalization.
3. **Curvature creates a witness.**  Any incompatibility between adjacent
   denominator ledgers should produce a finite witness denominator, just as
   nonzero p-curvature detects non-finite monodromy in the differential model.
4. **Multi-stranger resource bound.**  Generalize HYP-2444 from one stranger:
   each stranger can spend shell-27 class, 13-clock, or low-clock coverage, but
   the p-curvature ledger should prevent these resources from being spent
   independently at every lift.

## Tournament Analysis

Vertices are proof routes, not primes, denominators, or runners:

```text
p_curvature_operator
> partial_frobenius_support_gate
> lrc14_Q27_curvature_ledger
> marked_blocker_hypergraph
> raw_q_blocking
> naive_scalar_modp
```

The pairwise observable is

```text
(operator_retained, exactness, LRC leverage, computability, -risk).
```

The stored route tournament is transitive:

```text
score histogram: {0:1,1:1,2:1,3:1,4:1,5:1}
directed 3-cycles: 0
Hamiltonian paths: 1
```

Assumption challenge: alternate vertex sets considered were primes,
denominators, unit twists, runners, local p-jets, blocker supports, Pisano
classes, Bprime targets, Frobenius twists, monodromy generators, and proof
obligations.  The chosen proof-route quotient preserves whether a scalar local
test retains its operator/carry side channel.  It destroys actual differential
geometry and arbitrary multi-stranger LRC interactions; those remain the open
work.
