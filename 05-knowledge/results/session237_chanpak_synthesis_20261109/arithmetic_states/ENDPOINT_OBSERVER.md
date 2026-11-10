# Minimal modular observer for the physical apex endpoint

Status: **PROVED**, with FINITE-EXACT corroboration. This is a synthesis of two Session237 mechanisms, not a coverage theorem. Chan–Pak remains OPEN.

## Source, map, and exact scope

The audited physical family in `../endpoint_slack/PROOF.md`, section2, consists of the transitive chain x1->...->x_l->a and an apex b, with a->b and b->x_i exactly when bit l-i of m is one. Its endpoint count at a and full Hamilton-path count are

    E(m)=m-Q(m),       H(m)=1+2m-Q(m).

The digit-pair energy Q is the inherited Session205/THM-4888 statistic, not a new definition. `PROOF.md` in this directory proves reachability of every pair (m mod2M,Q(m)/2 modM) by the older Session206 block-pumping mechanism. This note combines those two exact maps.

The object counted is spanning directed paths ending at the designated vertex a. E(m) is not a complete Hamilton-path count, and its residue behavior does not imply complete-H coverage. Empty binary input m=0 is allowed; it gives E(0)=0 in a directed edge. Initial zeros only add vertices with the same scalar count and leave the automata below at their initial state.

## Theorem

For M>=1, read the binary digits of m most significant first and output E(m) modulo M after each complete input word. The exact minimum number of states of a deterministic finite automaton with output is

    1          if M=1,
    2 M^2      if M>=3 is odd,
    M^2 / 2    if M is even.

The automata are completely explicit. All their states are reachable, and any two distinct states are distinguished by an input suffix of length at most one.

### Recurrences

Since Q(2a)=Q(a) and Q(2a+1)=Q(a)+a+(a mod2),

    E(2a)=a+E(a),
    E(2a+1)=E(a)+1-(a mod2).

The second formula shows exactly why retaining the parity lift matters at odd moduli.

### Odd M>=3: 2M^2 states

States are pairs (r,q) in Z/(2M) x Z/M, initially (0,0), representing m modulo2M and Q(m) moduloM. The output is r-q moduloM. For b in {0,1}, the transition is

    (r,q) -> (2r+b mod2M, q+b(r+(r mod2)) modM).

This is well-defined: r modulo2M retains both r moduloM and its parity. Every state is reachable, because the simultaneous pumping theorem supplies any m modulo2M and any Q(m)/2 moduloM, and multiplication by2 is invertible modulo odd M.

For distinction, first compare the present outputs. If equal, the difference after a zero is r-s moduloM. Unless this vanishes, a zero distinguishes. If r=s moduloM and r=s as representatives in [0,2M), equal present outputs imply q=t, so the states coincide. Otherwise r and s differ by M, hence have opposite parity. After a one the output is present_output+1-(r mod2), so the two states differ by1 moduloM. Thus every pair of distinct states is distinguishable by the empty word, zero, or one. Reachability plus distinction proves minimality.

### Even M: M^2/2 states

States are (r,q) with r in Z/M and q an **even residue** moduloM, initially (0,0). The output is again r-q. Transitions are

    (r,q) -> (2r+b modM, q+b(r+(r mod2)) modM).

Because M is even, parity is determined by r moduloM, so these transitions are well-defined. The increment r+(r mod2) is even, preserving the allowed second-coordinate residues.

Every state is reachable: lift r to a chosen residue modulo2M, choose z moduloM with 2z=q moduloM, and apply the simultaneous reachability theorem to (m modulo2M,Q(m)/2 moduloM). Projecting gives the desired (r,q). There are M*(M/2) states.

Suppose two states have identical present outputs. Their outputs after zero differ by r-s moduloM. If that also vanishes, r=s and therefore q=t moduloM, so the states are identical. Thus the empty word and zero distinguish every distinct pair, proving minimality. In particular M=2 gives exactly two states, not eight.

### M=1

Every output is zero. A single looping state suffices and is minimal.

## Relation to the full-count observer

`PROOF.md` establishes that the **normalized full count** (H(m)-1)/2 moduloM needs exactly2M^2 states for M>=2. The endpoint E(m) moduloM needs the same number at odd moduli and exactly one quarter as many at even moduli.

This fourfold saving is an exact finite-state phenomenon, not an asymptotic estimate. The endpoint output uses Q rather than Q/2. At even modulus this retains only M/2 values of Q and permits reducing the parameter coordinate from modulo2M to moduloM. At odd modulus multiplication by2 is invertible, but reducing m moduloM would forget parity and corrupt the digit-one transition.

A useful check: for odd u>=3, E modulo u and E modulo2u both require2u^2 states. Adding the parity of the endpoint to the output requires no extra states because the odd-modulus observer already retains the parameter parity and Q is even.

These state lower bounds concern deterministic MSB-first residue-output automata only. They do not constrain arbitrary computational methods or prove a lower bound on tournament order. The exact bridge retains the input digit word and designated endpoint; it loses the general rooted response state used by port operations.

## Verification and hostile tests

Run:

    python3 04-computation/experiments/session237_arithmetic_states/verify_endpoint.py

The standalone verifier imports no endpoint producer or modular producer code. It enumerates reachable automata for M=2,...,50 and verifies the stated distinguishing signatures. For all m<65536 it checks both recurrences and eight observers. For m<256 it constructs literal apex tournaments with the correct bit order and counts all endpoint totals using backward subset DP, checking both E and H independently of the recurrences.

Hostile checks show that (r,q)=(0,0),(M,0) for odd M have equal present output but different digit-one continuations, so the parity lift cannot be discarded. M=2 and M=1 test the exceptional boundaries. Explicit checks use exceptions rather than assertions.

Output: `04-computation/experiments/session237_arithmetic_states/endpoint_results.json`; optimized-Python replay is byte-identical. No completeness, asymptotic distribution, or arbitrary finite-bank conclusion is inferred.
