# Exact four-register physical construction behind the three-port near-match

Status: proved all-word physical construction and candidate-specific obstruction. Not a universal compiler, not an F(8,137) proof, and not an obstruction to other finite-control machines.

## Physical construction

Start with the six-vertex tournament whose outgoing bit rows are (4,5,48,7,11,27), retaining vertices (0,1,2) as labeled ports. A newly introduced vertex loses to every bulk vertex.

- Digit 0: introduce a vertex beating ports 0 and 1, losing to port 2, then forget the newborn. Keep the same labeled ports.
- Digit 1: introduce a vertex beating port 2, losing to ports 0 and 1. Forget old port 1. Relabel (newborn,old port 0,old port 2) as (0,1,2).

After any finite word, close by adding two vertices successively, each beating ports 0 and 1 and losing to every other existing vertex. Mark old port 2 as root. Each digit adds one actual vertex; closing adds two. No deletion of a vertex occurs: forgetting only removes its port label.

For the closed tournament let h be its Hamilton path count, s the count starting at the marked root, and L the number of unordered spanning covers by one or two directed paths with one component ending at the root. Put B=h+L-s and D=3B+2s.

## All-word register theorem

Write M0,M1 for the open path-forest transfers, and D(x),s(x) for the fixed closing readouts. Define the additional integer readouts

u(x)=D(M0 x)-2D(x), v(x)=D(M1 x)-2D(x).

The following exact registers hold after EVERY word:

Initial (D,u,v,s)=(563,-23,-17,1).

0: (D,u,v,s) -> (2D+u, u, (D+7u)/6-36, 1).

1: (D,u,v,s) -> (2D+v, 2u-v-36, (D+13u)/6-v-78, 1).

Although fractions occur in the presentation, all these readouts are integers on every actual reachable state. In particular s=1 holds at every length, without assuming a finite scan generalizes. The construction therefore supplies actual rooted tournaments with a permanently unit start count and a four-dimensional exact linear representation.

## Exact finite proof

The self-contained verifier regenerates the 210-coordinate unordered path-forest presentation with at most two bulk-start components, both digit matrices, both closing rows, and the seed vector. The seed vector comes from an independent subset-path and canonical partition DP.

The full reachable linear subspace of the 210-coordinate state has dimension EIGHT, not four. Its image under the four readouts (D,u,v,s) has dimension four. The distinction matters: this is an observable quotient, not a claim that four numbers recover every forest count.

Let V span the full reachable space, R be the four-readout matrix, and A0,A1 be the displayed 4 by 4 transfer matrices. Exact rational elimination verifies

(R Mi - Ai R)V=0, i=0,1.

The span algorithm terminates only after closure under both Mi has been checked. The identities therefore prove the register theorem by induction for all words. Seven separate literal graph checks, for all words of lengths 0,1,2, compare closed B,s with the transfer. These are regression controls, not the all-order proof.

As usual, discarding patterns with more than two bulk starts is safe for these readouts: a future newborn cannot enter a bulk-start component, so bulk starts never decrease. No vector-realizability assumption is used.

## Why a bounded carry cannot repair this candidate

With nominal binary integer N initialized at 90, the earlier candidate sought D=6N+23. Its homogeneous runs match, but exact cross-run formulas are

D(0^n1)=1125*2^n-13 for n>=1;
D(1^n0)=1050*2^n+53 for n>=0.

Thus the discrepancies from 6N+23 are respectively

45*2^n-42 for n>=1;
42*(1-2^n) for n>=0.

The verifier proves the scalar recurrence f(n+2)-3f(n+1)+2f(n)=0 on the entire invariant span of each repeated-digit orbit, then checks the two initial values. This is an exact proof of the formulas, not extrapolation.

Consequently no finite-state, bounded additive correction to this fixed closing readout makes these same two physical operations encode the nominal binary integer. Indeed the first discrepancy is unbounded. The claim permits arbitrary finite control for choosing the additive correction, but NOT different physical digit operations or different physical closing contexts. It also does not rule out nonlinear recoding of words.

## Reproduce

In this directory, run `python verify.py` and `python -O verify.py`. Dependencies: SymPy. The proof packet consists of RESULT.md, verify.py, literal_forest.py, REVIEW.md and independent.py. Run `python independent.py` for the additional 15-case independent endpoint-DP audit. Exploratory compression JSON and scripts are not needed by the certificate.

## Remaining constructive gap

These identities make the actual small-width positive architecture explicit and reduce its inverse problem from 210 forest coordinates to three nonconstant registers. They do not show that its scalar outputs, or outputs of a larger controlled family, contain any tail. In particular they do not solve the all-order physical digit compiler or the eventual odd Hamiltonian-path spectrum threshold. A future compiler must use the auxiliary registers constructively, change the operations/contexts, or use a different architecture; bounded affine carry repair of this pair is insufficient.
