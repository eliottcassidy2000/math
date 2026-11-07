# Session210: a local cycle-spectrum bridge and closure observations

The degree-defined class of tournaments with an arc a->b satisfying d+(a)=d-(b)=1 has, at order n+2, exactly the Hamiltonian-cycle count spectrum of n-vertex tournament Hamiltonian paths. From order4 the gate is unique, giving an unmarked class bijection and preserving automorphism groups. This is an elementary structural transport, not a proof of the unrestricted Chan–Pak spectrum conjecture.

The stronger result is a universal operator on ordered path-cover polynomials:

    P_Gamma(T)(q)=q[P_T(q+1)+qP_T(q)-P_T(q-1)].

For known core order n, the first floor((n-1)/2) Hamiltonian counts in the nested-closure sequence recover P_T exactly. The same observations recover the complete permutation forward-adjacency distribution by an explicit binomial transform. Scalar H first fails to suffice at n=5: two H=15 cores close to H=123 and135. A separate order5 hostile shows that a sole entry and sole exit without their connecting gate can have an even cycle count.

See PROOF.md for definitions, elementary proofs, size gradings, inherited sources, hostile controls and limitations. RESULT.json contains deterministic exact checks; verify.py uses standard Python only. An independent audit accompanies publication. No broad novelty claim is made.
