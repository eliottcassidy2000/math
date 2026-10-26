# Session65: the supplied sequence, weak Collatz, and a guarded 3-adic bridge

Status: source-pinned research and elementary derived results; no universal Collatz or all-input halting claim. No external novelty claim. No repository files edited.

## 1. Exact identification, not sequence fitting

The user's fourteen values are the nonzero macro states in the blank-tape run of the Marxen–Buntrock five-state binary busy beaver. Put

    g(3k)=5k+6,  g(3k+1)=5k+9,  g(3k+2)=HALT.

Then the complete macro orbit is 0, followed by the supplied list, followed by HALT. The successful branch word is 00111011010001; the next residue is 2. In particular there is no numerical continuation after12284 in this specified system.

Primary sources:
- Pascal Michel, [Behavior of busy beavers](https://bbchallenge.org/~pascal.michel/beh), updated September2026, winner subsection and Collatz-like-problems section. The native configuration is C(n)=...0(A0)1^n0... . Costs for inputs3k,3k+1,3k+2 are respectively5k²+19k+15,5k²+25k+27,6k+12. Their sum on this orbit is47,176,870; the final tape has4098 ones.
- Michel, [historical survey v7](https://www.cs.cmu.edu/~cdm/projects/BusyBeaver/Michel2022-busy-beaver.pdf), arXiv:0906.3749v7, §5.2.1,p26, provides the same formulas.
- Aaronson, [The Busy Beaver Frontier](https://scottaaronson.com/bb.pdf), §4 identifies the orbit; §5.6,Conjecture20 asks whether every natural seed halts. Its equivalence to halting on every finite input is credited to Michel. Michel's September2026 page still lists the corresponding all-input question as open.

These are three distinct predicates: this one seed halts (verified); every nonnegative seed halts (open); ordinary Collatz has no nontrivial positive cycle (a separate open statement).

For every nonnegative live x, g(x)>x, since g(x)-x is(2x+18)/3 or(2x+22)/3. Thus g has no nonnegative cycles. A nonhalting nonnegative orbit necessarily diverges to infinity. This is precisely the failure type a weak-Collatz cycle detector does not detect.

## 2. Source and finite audit of tinycryptids

Pin: [calr3/tinycryptids commit3c43f97768d4ea06d70e0eb154ba88060bc36dec](https://github.com/calr3/tinycryptids/tree/3c43f97768d4ea06d70e0eb154ba88060bc36dec),2026-10-24. Exact source blob hashes are in tinycryptids/source-manifest.json. The local files are source snapshots, not publication payloads.

[WeakCollatz/Main.lean](https://github.com/calr3/tinycryptids/blob/3c43f97768d4ea06d70e0eb154ba88060bc36dec/lean/TinyCryptids/WeakCollatz/Main.lean) states TM.WeakCollatz.halts_iff:

    Halts(delta,start) iff exists n,k>0, C^k(n)=n and n notin{1,2,4}.

The14-state binary machine and10-state ternary version use shortcut T(n)=n/2 for even n,(3n+1)/2 for odd n. They enumerate a round maximum M, start v=M, advance T until1 or a value≥M, continue to the next M at1 or>M, and halt on return M. To prove completeness choose the least maximum of a nontrivial cycle. An earlier round cannot run forever trapped below its threshold: finite-state arithmetic confinement would yield a smaller cycle. No claim about divergent ordinary Collatz orbits follows.

The downloaded14-state and10-state main theorem import closures contain no sorry,admit,axiom,unsafe,native_decide tokens. Root TinyCryptids.lean imports both main files. Check.lean asks Lean to print theorem axioms; PROOF.md reports propext,Classical.choice,Quot.sound only. **Neither lean nor lake is installed here; no local Lean build or actual #print-axioms evaluation was performed.** These remain source-level audit facts and an external formalization claim, not a locally rebuilt theorem.

The [Cryptids page revision9619](https://wiki.bbchallenge.org/w/index.php?title=Cryptids&oldid=9619) lists14-state binary and10-state ternary weak Collatz,23-state ternary odd-multiply-perfect(unverified), and37-state binary odd-perfect(with formalization,not independently verified). Do not substitute odd perfect for odd multiply-perfect.

The wiki's independent-verification link points to [sligocki/busy-beaver validator.rs](https://github.com/sligocki/busy-beaver/blob/266ff1328852a87a008cc6bca0c22e215f050338/rust/src/validator.rs), pinned266ff1328852a87a008cc6bca0c22e215f050338. test_weak_collatz_14 contains symbolic chain and inductive rules and the exact machine string; that test has no Admit, although unrelated tests in the file do. We did not compile/run the Rust verifier or audit its whole trusted implementation. Its local rule certificate should not be inflated into our independent global proof.

Actual executed finite checks:
-14-state: all1829 C(M,v),2≤M≤60,1≤v≤M, including synthetic halting configurations;1839 blank-run macro visits for rounds2..300;100,000,000 native steps without halting.
-10-state ternary: the analogous1829 configurations;1837 blank-run macro steps for rounds3..300;100,000,000 native steps without halting.
-The native five-state MB table, separately simulated by the C interpreter, halts at47,176,870 steps with4098 ones, matching the independently summed macro costs.

Reproduction: cd tinycryptids; gcc -O3 -o tools/sim tools/sim.c; python3 weakcollatz/test_e2e.py; python3 weakcollatz/sym3/test_e2e.py. See verify_bridge.py and bridge_results.json for arithmetic controls.

## 3. Elementary derived theorem: exact binary survivor cylinders

Extend g to the two live residue classes of Z_3. For r∈{0,1}, set h_r(y)=(3y-18-4r)/5. Because5 is a3-adic unit, h_r maps Z_3 bijectively onto r+3Z_3, scales3-adic distances by1/3, and inverts the r branch of g.

For a word w=(r_0,...,r_(m-1)), put h_w=h_(r_0)∘...∘h_(r_(m-1)). Its image is one ball of radius3^(-m), hence one residue class modulo3^m. Different words have disjoint images: at the first different digit the corresponding forward iterates have different residues. Therefore exactly2^m residue classes survive m successful steps.

For every infinite binary word r, nested balls h_(r|m)(Z_3) have one point

    pi(r)=-sum_(j>=0) (18+4r_j)3^j/5^(j+1),  convergence in Z_3.

This gives a continuous bijection from compact binary sequence space to the closed infinite-survivor set K, with continuous inverse obtained by residue itineraries; g∘pi=pi∘shift. Thus K is a compact binary full-shift Cantor repeller. It is NOT asserted to equal the standard base3 digit Cantor set: its branches scale by3/5,not3. The displayed pi is the explicit coding conjugacy.

The mth survival set has Haar measure(2/3)^m. Hence K has Haar measure0. This alone would not control integers. But the exact residue count also gives natural density(2/3)^m for integer seeds surviving m steps. The upper natural density of nonhalting nonnegative seeds is at most(2/3)^m for every m, therefore0. More quantitatively, for N≥1 and m=ceil(log_3 N), each safe residue class contains at most one integer in[0,N), so

    #{0≤x<N : x never halts} ≤2^ceil(log_3 N) ≤2 N^(log_3 2).

These are elementary derived density-one-halting and sparse-exception results, NOT all-input halting. No exchange of density and expectation limits is justified here. Constant codes give the hostile fixed points pi(000...)=-9 and pi(111...)=-11; K is certainly nonempty. Every finite word has a negative rational periodic-point extension, so all finite safe-prefix tests have nonhalting3-adic completions.

This is a useful hostile to ternary-digit arguments around Erdős/Sloane: a zero-density survivor carrier need not exclude a specified sparse arithmetic orbit. One must retain the particular arithmetic embedding and prove intersection emptiness, rather than use ambient density.

## 4. Exact arithmetic-tournament bridge and its limitation

Use the existing THM-4754,01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md: x→y when the3-adic unit part of y-x is1 modulo3. This is an intrinsic binary relation, with no ties between distinct3-adic points.

On any common admissible m-step cylinder,

    g^m(y)-g^m(x)=5^m(y-x)/3^m.

Removing powers of3 does not change unit part, while5≡-1mod3. Hence the orientation reverses iff m is odd. The map is exact only when the common branch address is retained; crossing different branches introduces different offsets.

More strongly, if r,s first differ at j, then

    v_3(pi(s)-pi(r))=j,
    unit_3(pi(s)-pi(r))=(-1)^j(s_j-r_j) mod3.

Proof: the common inverse prefix scales by(3/5)^j, and before that prefix, the difference between the branch1 and branch0 images is1mod3. Thus the induced tournament on K is transitive: order codes lexicographically after complementing every odd-position digit. Every finite safe-cylinder leaf tournament is likewise transitive. This fact is about the pairwise carrier, not absence of dynamical cycles: K contains the fixed points-9,-11.

Loss ledger: source=guarded macro system on nonnegative integer seeds; target=3-adic survivor/cylinder tournament; map=integer embedding and explicit itinerary pi; preserved=admissible finite words, exact affine differences, orientation; lost if only tournament retained=integer membership, branch cylinders, native runtime, designated seed, halting exit. Required sidecars are branch address, arithmetic embedding, and native clock. Cheapest hostile=-9,-11 and the known terminating seed0. No multiply-perfect or LRC implication follows.

## 5. Multiply-perfect hostile checks on the actual list

Direct divisor summation shows only6 in the supplied list has integer abundancy σ(n)/n. Already16 has31/16, so the sequence is not a list of multiply-perfect numbers. There is one successor coincidence, σ(946)=1584=g(946), because946=2·11·43 and3·12·44=1584. It is not an orbit-wide divisor-sum rule; all other displayed live positive transitions fail σ(x)=g(x).

The general bridge to odd-multiply-perfect search is logical, not a hidden identity: a bounded finite certificate (n,k,σ(n)=kn) can be enumerated and cause a blank-tape machine to halt, whereas an infinite divergent orbit cannot be certified by merely waiting. Any stronger arithmetic connection needs a new explicit map.

## 6. Finite universe and controls

verify_bridge.py enumerates all residues modulo3^m for1≤m≤9, independently constructs inverse-word representatives, compares the sets and all word labels, and checks shifted representatives. It checks43,435 pair orientations on all word sets of lengths1..8 against alternating lex order. It also verifies the supplied orbit, all exact macro costs, negative fixed-point controls, and divisor-sum hostiles. Output is bridge_results.json. These finite checks support the elementary proofs above; they do not replace them or establish all-input halting.

## 7. Sealed replay and authorship boundary

`verify_bridge.py` is independently authored here and uses only Python's standard library. It has explicit optimization-safe `require` gates, reads no input files, and has no dependency on the source checkout or ambient PYTHONPATH. The commands `python3 verify_bridge.py > bridge_results.json` and `python3 -O verify_bridge.py > bridge_results_optimized.json` both completed successfully with byte-identical JSON. Hashes and observed external outputs are recorded in `verification_receipt.json`.

The fourteen-state and ternary test suites, Python TM wrapper, and C native interpreter are **external tinycryptids code run here**, not independently authored simulators. The five-state machine table was independently entered from Michel's source, but its native execution used that same external C interpreter. Its comparison with the separately authored arithmetic macro-cost sum is an independent arithmetic path, not an independent second native simulator. External stdout in the receipt is explicitly identified as transcription of observed tool output rather than a newly captured raw log. No Lean or Rust build is claimed.
