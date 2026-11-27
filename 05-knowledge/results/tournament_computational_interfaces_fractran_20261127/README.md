# The tournament port as a computation: exact memory, finite quotients, and proof obligations

This links the [guarded machines](../guarded_moore_pollock_collatz_cullen_machines_20261127/README.md), [simultaneous one-apex inverse](../chanpak_decimal_kernel_20261122/LOCAL-WORDS.md), and [odd-modulus tournament selection](../collatz_modular_multiplication_tournament_classes_20261127/README.md). A useful test of a representation is whether it preserves **the response to the next legal edit**, not merely today's output.

## 1. Equal heights need not be the same computational state

Fix a transitive spine and an apex. Read its contact word as a binary integer \(x\), with 1 meaning that the apex beats the corresponding spine vertex. Set \(c=(H-1)/2\). The established one-apex append rule is
\[
(x,c)\overset{b}{\longmapsto}
\big(2x+b,\ c+(1-b)\lceil x/2\rceil\big).
\]
The physical edit is to append a last spine vertex and orient its apex contact by \(b\). The height is a readout, not the full state. Indeed, **two strong order-five examples** have words 1000 and 1110 and both have height 9. Appending exactly the same word 10 gives strong order-seven examples of heights 27 and 39. The replay counts paths independently by subset DP.

The complete continuation behavior over the integers distinguishes every different address \(x\): compare the outputs after 10 and 110, which are respectively \(c+x+1\) and \(c+2x+2\) in the height register. These determine both \(x\) and \(c\). The point is not that the exact count \(c\) and address need independent storage as integers (on reachable words \(c\) is a function of \(x\)); it is that the scalar count cannot replace the address.

### A sharp finite-state theorem

For \(M\ge2\), an automaton reading arbitrary apex words beginning in 1 and outputting \(c\bmod M\) needs **exactly \(2M^2\) states after the first bit**. More precisely, all of these states remain distinguishable even if output is requested only after a suffix ending with 0, when the host is strong.

The pair \((x\bmod 2M,c\bmod M)\) is sufficient by the update formula. All pairs are reachable by the [existing binary-controller theorem](../chanpak_decimal_kernel_20261122/LOCAL-WORDS.md), starting in particular from the state of the first bit 1. To prove necessity, the two probes 10 and 110 recover \(x\bmod M\) and \(c\bmod M\). The only remaining ambiguity is which lift of \(x\) to mod \(2M\) was used. Probe 0 distinguishes those lifts: its increment changes by \(M/2\) for even \(M\), and by one of \((M-1)/2,(M+1)/2\) for odd \(M\), always nonzero mod \(M\). Thus the three probes \(0,10,110\) distinguish every pair of states. All end in 0, so they are actual strong tournaments rather than formal terminal queries. At \(M=1\) the output is constant and there is only one observable state. The replay independently checks reachability and probe injectivity for every \(M=2,\ldots,128\).

This is a minimal-state statement for the specified append interface. It is not a lower bound on all algorithms for Hamiltonian path counts or for Chan–Pak.

## 2. The same edit in FRACTRAN, tape language, and a local lattice

In [Conway's FRACTRAN](https://gwern.net/doc/cs/computable/1987-conway.pdf), scan an ordered list of fractions and multiply by the first one for which the answer is integral. Registers can be prime exponents; an extra prime encodes control phase. Here is a concrete program for **one** append bit \(b\) on the unbounded integer state. Input \(7\cdot2^x3^c\); put \(P=3^{1-b}\); use the ordered list
\[
\frac{11\cdot5^4P}{7\cdot2^2},\quad
\frac7{11},\quad
\frac{13\cdot5^2P}{7\cdot2},\quad
\frac{13}7,\quad
\frac{17\cdot5^b}{13},\quad
\frac{19\cdot2}{17\cdot5},\quad
\frac{17}{19},\quad
\frac{23}{17}.
\]
The first two rules consume pairs from the exponent of 2 and put four units per pair in the exponent of 5. The next two cases consume a leftover bit if present, and the half-ceiling contribution goes to the exponent of 3 exactly for \(b=0\). The second phase writes \(b\) and moves all of the temporary exponent of 5 back to 2. The terminal factor is 23, and there are no remaining applicable rules. The exact output is \(23\cdot2^{2x+b}3^{\,c+(1-b)\lceil x/2\rceil}\). This is an actual finite fraction list per bit, not a claim of universality of the tournament family. Iterating arbitrary words additionally requires a bit stream and a control program to issue these calls in order.

* In a Turing implementation, the same registers can be unary marks between delimiters and a head sweeps through them. Each control phase, loop test and register update has the same effect as above; time and representation size may differ.
* To encode this regular tape machine in a radius-one cellular automaton, use cells holding a tape symbol and an optional head-state marker. At each tick the old head writes its new symbol; its neighbor copies the arriving head-state according to the transition; all other cells hold. Starting with one head gives exactly the tape evolution. This does not require choosing a universal cellular automaton rule or adding simultaneous second heads.
* The one-apex update and halved Collatz both use parity-sensitive affine arithmetic, but branch selection differs. Here the tape supplies \(b\), the next externally chosen arc. Collatz computes its branch from the current parity. Given a legal parity word, it likewise has an affine update, but that does not transfer a global termination argument from one machine to the other.

In particular a machine representation is faithful only after specifying its input, legal scheduling, observable output, and what counts as halting. Forgetting the address here destroys fidelity even at order five. Conversely nothing about the existence of universal [generalized Collatz](https://people.cs.uchicago.edu/~simon/RES/collatz.pdf) or FRACTRAN programs makes this specific one-apex transducer or the \(3x+1\) problem undecidable.

## 3. A resource bound makes the Chan–Pak target test genuinely finite

A tournament's strong components form an ordered join and its path count is the product of their counts. Delete singleton components, which have count one. For a strong component of order \(s\ge4\), [Busch's exponential bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/) gives \(H_s\ge 5^{(s-1)/3}> 3^{s/3}\), since \(5^{s-1}>3^s\) for \(s\ge4\). For \(s=3\) the strong component is a directed triangle of count 3. Consequently every realizable \(H>1\) has a witness of order at most **\(3\log_3 H\)**, after deleting singletons. The constant is sharp: an ordered join of \(r\) directed triangles has count \(3^r\), and every singleton-free witness at that height has order at most \(3r\). More exactly the bound is attained by those joins, not necessarily by the *minimum-order* witness for that count. Target \(H=1\) has a singleton witness.

As a result, membership in the exact spectrum is uniformly decidable in time \(2^{O((\log H)^2)}\) by enumerating all tournaments through this order and counting paths. This is a finite resource guarantee even when the answer is **no**. The alternative branch of a search for one target is not an unbounded waiting room.

| Problem | What a particular input gives you | What the conjecture still asks |
|---|---|---|
| Chan–Pak | A bounded exact witness-or-rejection for each \(H\) | No failed target outside the prescribed exceptions |
| Pollock | A finite five-track positive-summand reachability test for each target | No failed target anywhere |
| Degree-57 Moore | One finite rooted permutation constraint search | A witness for that fixed size, or an exhaustion proof |
| Collatz | A finite certificate if the orbit from \(x\) reaches 1 | Such a certificate exists for every positive start |
| Cullen or Woodall infinitude | A finite prime/composite decision for each index | A prime index beyond every given cutoff |

Thus Chan–Pak and Pollock admit a program that halts upon locating a counterexample; their universal statements have a decidable per-target predicate (a \(\Pi^0_1\) presentation). Collatz and the two infinitude statements have natural \(\forall x\,\exists t\) presentations (\(\Pi^0_2\)); this is an upper-level classification in the arithmetical syntax, not a completeness or independence claim. A finite Moore instance is a different kind of obstacle: search size, not unbounded quantifier alternation. The same hardware can encode all five problems, but the proof obligations do not become equivalent.

Replay: run the supplied JavaScript with Node, for the literal tournament split, minimum-state finite checks, and the 1290 exact FRACTRAN input transitions.
