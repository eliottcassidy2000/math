# Finite synchronous witness selectors cannot cover the fixed join bank

Status: PROVED by the elementary argument below; awaiting independent review. Not a Chan–Pak proof or a refutation of F2 completeness.

## Setting and exact scope

Let g(0)=1, g(2u+1)=g(u), g(2u)=g(u)+u+(u mod 2), and e(v)=1−(v mod 2), including v=0. For a fixed finite bank Y of triples (a,b,c) with a,c positive integers and b nonnegative, let

    F_y(v)=a g(v)−b e(v)+c(v+1).

The order-at-most-seven rooted bank of THM-4875(D) has these properties. Its finite coverage through 2^35 remains FINITE-EXACT. The first moment of THM-4885 does not supply pointwise coverage.

A synchronous binary witness selector is a finite automaton recognizing some set of valid tuples (y,v,N) with N=F_y(v), reading the binary expansions of v and N together from their most significant ends, with leading zero padding to equal length. The finite label y can be a separate finite track or selected by separate start states. The automaton need not recognize every witness, need not choose a unique witness, and can omit arbitrary valid witnesses.

**Theorem.** For every such selector, there is a fixed R such that every accepted v has at most R binary runs. Consequently the number of distinct selected outputs 1<=N<=X is O((log(X+2))^(R+1)). In particular no such selector witnesses all sufficiently large odd integers, or any positive-density set of integers.

This is stronger than nonregularity of the full graph, but applies only to synchronous finite-state relations retaining the literal binary witness v. It does NOT rule out a finite automaton for output support alone, an automaton using an additional unbounded coordinate, a different witness encoding, or a non-finite-state construction. In particular an automaton recognizing all odd N without computing any v is outside its scope.

## 1. Positive drift for every mixed binary period

Let w be a binary word of length k containing both digits, B=2^k, and r its integer value. Appending w maps v to Bv+r. Once the incoming parity is fixed to the final bit of w, the g increment is affine:

    g(Bv+r)−g(v)=C v+K.                         (1)

Indeed each appended zero adds the current prefix plus its parity, and an appended one adds zero. Every current prefix is affine in v; the first parity is fixed, and later parities are literal preceding bits.

Write v_m for the integer represented by u w^m, for an arbitrary fixed prefix u. For m>=1, incoming parity is stable. Summing (1) and v_(m+1)=Bv_m+r gives constants A,E with

    g(v_m)=A B^m+D m+E,
    D=K−Cr/(B−1).                              (2)

Crucially D>0 for every mixed w. To prove this, evaluate the affine increment at the virtual periodic prefix x=−r/(B−1). Under x↦2x+bit, the successive x values are negatives of cyclic rotations of r divided by B−1. Thus each lies strictly between −1 and 0. Assign parity equal to the preceding literal bit, including across the cyclic boundary. The sum of contributions x+parity at the zero positions is exactly D.

Partition the cyclic word into its maximal zero runs. If a run has length l, write −q for x just before its first zero. It follows a one, so the parity contributions are 1,0,...,0. Its contribution is

    1−q(1+2+...+2^(l−1))=1−(2^l−1)q.

After those l zeroes, the virtual prefix is −2^l q in (−1,0). Therefore 0<q<2^(−l), and the run contributes strictly more than 2^(−l)>0. Summing proves D>0.

Now append any fixed suffix z. Its effect on g is affine in v_m with coefficient of g(v_m) exactly one; the parity dependence is constant for m>=1. Hence the positive linear coefficient D in (2) survives every fixed suffix. The integer represented by u w^m z itself has the form A_0 B^m+E_0, and its parity is fixed for m>=1. Therefore

    F_y(value(u w^m z))=A_1 B^m+aD m+E_1,       (3)

with aD>0.

## 2. Pumping contradicts every mixed loop

Fix one bank label y by taking the corresponding regular-language slice. Suppose an accepting automaton path contains a loop whose input v-track word w contains both bits. Pumping that loop produces accepted pairs for all repetition counts m. The output-track word repeated on the same loop has the same length k. Consequently its numerical output is necessarily

    N_m=A_2 B^m+E_2.                            (4)

Leading zeroes do not affect this identity. But validity of the relation requires (3)=(4) for all m>=1, which is impossible: the second-order difference operator t_(m+2)−(B+1)t_(m+1)+B t_m annihilates (4), while its value on (3) is (1−B)aD, nonzero. Thus no loop occurring on an accepting path can have a mixed v-track label.

Discard states not lying on an accepting path. Every strongly connected component containing edges must consequently have all its internal v-track labels equal: if internal edges of both labels existed, strong connectivity would produce a closed walk containing both, contradicting the previous paragraph. Along any accepting path, components occur without repetition in the condensation DAG. Internal edges in each component contribute at most one run; each intercomponent edge contributes at most one additional run. If the trimmed automaton has q states, R=2q+1 is a safe uniform bound. For finitely many y slices take the maximum of their bounds.

## 3. Counting the outputs

Since g(v)>=1 and e(v)<=1,

    F_y(v)>=c(v+1)+a−b,

so any accepted 1<=N<=X has v<=max_y((X+b)/c), which is O_Y(X). Leading zero padding adds no new integer values. The number of binary strings of length l>=1 having at most R runs is

    2 sum_(j=1)^R binom(l−1,j−1).

Summing through l=O(log(X+2)), adding v=0, and multiplying by the fixed number of labels proves the stated polynomial logarithmic bound. Distinct output counts cannot exceed witness counts.

## 4. Explicit hostile and reproducibility

For v_m=(10)^m=2(4^m−1)/3,

    9g(v_m)=4^(m+1)+6m+5,
    9F_y(v_m)=(4a+6c)4^m+6am+5a+3c−9b.

This makes the obstruction visible without asymptotics: the linear drift cannot be encoded by pumping fixed input and output blocks in lockstep.

Run `python3 verify.py`. It checks all 2,026 mixed words of lengths 2 through 10 using exact rational arithmetic, independently computes virtual-prefix drift, and checks 24,312 recurrence identities for four starting prefixes. It also checks the displayed alternating identity for m=1..100. These computations are regressions; the proof covers every length.

## Inheritance and connection ledger

- Closest proved mechanism: THM-4885, `01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md`, binary reversal with logarithmic correction; THM-4875, `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, exact F2 carrier.
- Related earlier barrier: `05-knowledge/results/session217_exact_inverse_and_gate_limits_20261107/run-entropy/PROOF.md`, logarithmic run defect and low-run zero-density. That proof concerns a different fixed-core padding family and a preimposed run budget; the present result derives a bounded run budget from synchronous regularity itself.
- Hostile: every mixed pumped period, explicitly (10)^m. Corrected near miss: positive finite coverage and first moment are not a covering proof.
- Source/target/map: a finite-state witness relation maps a repeated paired word to actual integers (v,N). Preserved: exact equality N=F_y(v). The tempting quotient discards the linear drift in g. Sidecar required to restore it: an unbounded drift/carry coordinate or richer encoding, not a finite label.
- Positive frontier: the F2 output support might still be cofinite; the theorem rules out this particular constructive finite-state certificate, not the desired support statement.
- No literature novelty claim. No repository files modified and no theorem ID reserved.

## 5. Inherited curvature obstruction, specialized to F2

This is a specialization, not a new obstruction mechanism, of Session216 numeration `05-knowledge/results/session216_overnight_constructions_20261107/numeration/PROOF.md` Section 3. That result compares the unchanged quadratic f coefficient under universal free-block rewrites and credits the broader Session203 uniform-incidence curvature theorem. The identity g(v)+f(v)=v+e(v) identifies its h with the present a. The present calculation covers only fixed suffix insertion with unchanged u, a narrower operation class, but applies to every physical rooted bank because a is odd.

The exact parity split is

    F_(a,b,c)(2u+1)=a g(u)+2c(u+1),
    F_(a,b,c)(2u)=a g(u)+(a+2c)u+a(u mod2)+c−b.

The increasing c coefficient does give a complete target-dependent search, since every step halves the v argument. It does not close within a fixed physical bank after dividing the target by two.

More generally fix a suffix r of length k>=1. The inherited append law gives

    F_y(2^k u+r)=a_y g(u)+L u+P(u mod2)+Q.

There is no identity

    F_y(2^k u+r)=2^j F_z(u)+constant,

valid for all u, for any physical rooted triples y,z and j>=1. Their a coordinates are odd Hamilton-path counts; equality would require a_y=2^j a_z, impossible.

For completeness the coefficient comparison is justified, not merely formal: g(u) is not in the rational span of u, parity(u), and 1. If g(u)=A u+B parity(u)+C, then u=2^n−1 forces A=0 because g(2^n−1)=1, while u=2^n forces an unbounded value g(2^n)=2^n+1, contradiction. The same comparison holds if the additive constant is replaced by a parity-dependent constant. This does not exclude a piecewise target-dependent descent or a change of witness coordinate.

## 6. Stronger corollary: no rational finite-state witness selector, even unlabelled

**Corollary.** The bounded-run and polylogarithmic-output conclusions hold for any finite-state transducer recognizing a relation of ordinary binary string pairs (v,N) contained in the union of the graphs F_y, with y in a fixed finite bank. The transducer may have epsilon transitions and need not report the bank label. Leading zero padding is permitted. The integers, rather than the lengths of padded strings, determine validity. No uniqueness or totality is assumed.

First, uniformly for a fixed finite bank,

    F_y(v)=Theta(v) as v tends to infinity.     (5)

The lower bound was already proved. For the upper bound, g(v)<=v+1. One quick independent induction uses f(0)=0, f(2u)=f(u), f(2u+1)=f(u)+u+(u mod2), for which f>=0 and

    g(v)+f(v)=v+1−(v mod2).

Now take any loop on a successful transducer path whose input word w is mixed. Let k>=2 be its input length and l>=0 its output length. Pumping the loop m times gives an input whose value grows as a positive constant times 2^(km). Its output value is either constant, or grows as a positive constant times 2^(lm). This follows directly from the integer-value formula for a repeated word, and includes empty or all-zero loop words. The constant case contradicts (5). Comparing growth rates forces l=k.

With lengths equal, the output is exponential plus constant as in (4). For each possible bank label y the required value has the form (3), with its own positive linear coefficient a_yD. Equality for that y is possible for only finitely many integer m: if the exponential coefficient difference is zero the remaining affine function has nonzero slope; otherwise its exponential term dominates the linear term in absolute value eventually. A finite union of finite sets cannot account for all pump counts. This rules out every mixed input loop, even when the bank label is allowed to change as m changes.

The SCC argument survives epsilon transitions verbatim. A component cannot contain both a 0-labelled and a 1-labelled input edge; otherwise connecting paths give a mixed closed walk. Components with only epsilon input edges add no run. Thus input integers still have a uniform finite run bound, and the same counting argument proves the corollary.

This direct proof avoids assuming an external synchronization theorem. It covers ordinary asynchronous finite-state transducers and finite labelled bank choices. It does not cover automata with an unbounded arithmetic register, a non-binary witness representation, or a certificate that forgets the witness entirely.

## Independent review

The independent audit worker reviewed Sections 1–3 and reported PASS for the positive-drift calculation, suffix preservation, SCC bound and counting argument. Its independent alternating-word derivation agrees. It also supplied the unlabelled-bank finite-pigeonhole extension and direct asynchronous growth comparison recorded in Section 6. A separate full review artifact is being maintained by that worker. Claims in this packet remain scoped to the specified witness relations; Chan–Pak and F2 completeness remain OPEN.

## 7. What a different encoding must change

Fixed block radix changes, such as binary to base 4 or base 8, do not evade Section 6: finite-state decoding of each fixed-size digit block and composition of finite-state transducers would give the forbidden binary relation. More generally, the same conclusion applies to any proposed witness encoding with rational finite-state decoders for both v and N. This is the standard closure of rational relations under composition, also seen directly by taking product machines and matching the intermediate finite alphabet.

A variable base or an extra unbounded register is a genuinely different model and is not excluded. For example, g already has an exact streaming compiler with integer registers (v,g) and a parity bit p, initially (0,1,0):

    append 1: (v,g,p) -> (2v+1,g,1),
    append 0: (v,g,p) -> (2v,g+v+p,0).

The register g retains precisely the drift a finite state loses; the terminal readout is a g−b(1−p)+c(v+1). This compiler evaluates any chosen witness but supplies no target-dependent inverse-selection or all-target coverage theorem. No alternate encoding yielding such a theorem was established in this lane.
