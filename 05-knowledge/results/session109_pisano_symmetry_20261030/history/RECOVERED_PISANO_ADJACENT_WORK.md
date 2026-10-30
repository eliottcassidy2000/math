# Session 109: recovered Pisano-adjacent work and a ramified extension

Read-only source recovery, 2026-10-30. Source checkout `math`, HEAD
`1b91742ce55a202b73c80d56561ba56d00b503a6`. This checkout is dirty; no source
files were edited. Requested `session103_words_necklaces/source` was absent.
The historical claims below are reports of inspected repository material,
not a fresh external literature or current record audit.

## 1. Main inheritance and precise credits

### Metallic units and invariant tournaments

Primary source: `01-canon/theorems/THM-4656-wall-sun-sun-golden-orbit-explosion.md`,
source **mac-mini-2026-10-14-tourn**. Depends on opus S33's THM-4644.
For mu^2-a mu-1=0 and odd unramified p, the odd-order part of multiplication
by mu organizes Cayley-tournament connection sets into antipodally paired
orbits. The unusual prime-square lift delay is the metallic Wall–Sun–Sun
condition. The source credits the unit criterion to Bragman–Rowland (2022,
Theorem 10), McIntosh–Roettger (2007), and classical order lifting to Wall
(1960). Its new reading is the invariant-tournament orbit explosion, not
the arithmetic criterion. Golden rigidity is equivalent to a power-of-two
Pisano period; golden Paley primes connect to Fibonacci primitive roots.

Critical later extension: `05-knowledge/results/unit_lift_clock_20261015a.md`,
**2026-10-15**, Theorem U1. It ALREADY proves all depths in BOTH split and
inert full rank-two rings. Put pi=ord(mu mod p), m=oddpart(pi),
g=mu^(pi/m), h=v_p(mu^pi-1). Then

    ord(g mod p^k)=m p^max(0,k-h),
    N_k=sum_{l=1}^k (p^(2l)-p^(2l-2))/(m p^max(0,l-h)),
    number of invariant labelled Cayley tournaments=2^(N_k/2).

The split zero-divisor issue is already settled: norm one gives reciprocal
components, hence equal valuations of g^d-1. A broad arbitrary split unit
still needs two component valuations. Explicitly excluded: ramified primes
and p=2. Verification shipped: 61 full-ring/Burnside cases, 2,442,949 elements.
At a=2,p=31,k=2, 64 orbits refers to one scalar factor, whereas 61,568
refers to the full ring. Do not conflate the vertex sets.

**Duplicate credit correction:** opus-2026-10-14-S36's THM-4657 was ceded to
mac-mini THM-4656. Sources `HYP-9281-fibonacci-wieferich-wall-sun-sun-lucas-binomial-eleven-fifty-five.md`
and `strong_balance_pisano_golden_eleven_20261014.md` say so explicitly.
The eleven-square/icosahedron link is only a C5 orbit-pattern analogy, not
an isomorphism: cell adjacency has degree six and does not match the punctured
icosahedron. The 55=F_10 coincidence does not supply a general mechanism.

### Cat map, modular representations, Lucas traces

Sources THM-4762, THM-4763, HYP-9372, **opus-2026-10-22-S55**.
Q=[[1,1],[1,0]], A=Q^2=[[2,1],[1,1]]=[t,s]. For N>=3,
ord(A mod N)=pi(N)/2, with N=2 exceptional (order 3).
THM-4762: C=TST^-1S^-1 is scalar-normalization independent. Congruence
representation theorems bound its order. The D(Z_N) representation is
conjugate to the actual cat-map permutation on (Z/N)^2 and saturates the
bound. The 135-case modular-data census is numerical, not an all-category
classification. HYP-9372's deficiency in {1,2,4} remains OPEN. Central powers
need not act as scalars; the source explicitly corrects that converse.

THM-4763 uses the KNOWN FGST 2006 decomposition
R_(p+1) + C^2 tensor R_(p-1). Its cat-map consequence is spectral radius
phi^2 and trace R_(p+1)(A^n)+L_(2n) trace R_(p-1)(A^n). Do not claim the
decomposition itself as new. Likewise matrix/Pisano identification and
finite/quantum cat-map orders are classical, explicitly credited there.

### Zeckendorf: exact arithmetic versus lossy color projections

Source `05-knowledge/results/zeckendorf_guard_automaton_20261003.md`,
**2026-10-03**, with inherited 2026-09-25/26 color/carry notes cited inside.
For weights f_i=F_(i+2), retain X=sum e_i f_i and Y=sum e_i f_(i+1).
Appending digit d sends (X,Y) to (Y+d,X+Y+2d). Previous-digit memory rejects
11. Modulo m this gives at most 2m^2 states, exactly recognizing numerical
residue guards. Modulo 2 the linear part has order three. This already
provides the rigorous Fibonacci/Pisano finite-state connection.

The second register is indispensable: n=2 and n=8 have the same X mod3,
previous digit, and charge, but differing Y mod3; appending zero separates
outputs. The user's 35-symbol stripe is not the four-state XOR charge and
not an established infinite Wythoff word: the latter misses position 35.
Colored markers are typed representation data, not numerical residues.
Any symmetry quotient must explicitly state which arithmetic predicate and
which marker data it preserves. Existing audits cover 70,007 register cases,
524 typed supports and 14,760 compiler/automaton/decoder interface cases.

### Older Pisano/LRC bridges: retain scoped content, correct global shortcuts

THM-486, **claudebox-2026-06-11-S5**, credits the classical fixed points
pi(n)=n iff n in {1} union {24*5^k:k>=0}. Its Cayley–Dickson/Golay dimensional
connections are interpretations, not deductions from numerical equality.
Its divisibility pi(p)|(p^2-1) for p!=5 is valid also at p=2,
since pi(2)=3. The odd-prime proof via 2(p+1)|(p^2-1) must handle
p=2 separately.

THM-491, **claudebox-2026-06-11-S6**, connects shell valuations to Pisano
towers. Its formula pi(p^k)=p^(k-1)pi(p) must NOT be imported universally:
Wall–Sun–Sun delays are precisely the missing hypothesis. The displayed
3-adic example pi(3^k)=8*3^(k-1) is valid; the metaphor of matching towers
is not a general nonexistence proof for exceptional primes.

`07-reflections/lrc14-pisano-quotient-and-the-fibered-band-ladder.md`,
**codex-2026-06-11-P6**, says the quotient (Z/27)^*/± has nine classes,
numerically pi(27)/pi(3)=9. Its exact single-stranger obstruction includes
both a missing ±10 class mod27 and a killed 13-clock. The fibered Q27
lattice already covers the 936-row family; plain-shell band two adds no
coverage. Equality of nine clocks is not itself a canonical identification
of the doubling action with the Fibonacci matrix action.

## 2. Concrete new extension proposed: complete the ramified golden prime

This was NOT located in the inspected chain; novelty is relative to that
chain, not a claim against all literature. It fills U1's explicit p=5 gap.

Let O=Z_5[phi], varpi=sqrt(5)=2phi-1, R_k=O/(5^k), and g=phi^4=3phi+2.
Then g-1=3phi+1=(5+3varpi)/2 has varpi-valuation one. For every positive
integer j,

    v_varpi(g^j-1)=1+2v_5(j).

Proof: a prime-to-five power preserves valuation one; in the fifth-power
binomial expansion of 1+u with v_varpi(u)>=1, 5u is the unique least-valued
term (valuation v(u)+2), because u^5 has valuation >=v(u)+4 and intermediate
terms are larger. Iterate. The norm-one unit g therefore has order 5^k
in R_k, consistent with pi(5^k)=4*5^k.

For nonzero x, let r=2k-v_varpi(x). Its multiplication orbit has size
5^floor(r/2); there are 4*5^(r-1) elements in that valuation layer. Thus

    N_k=sum_{r=1}^{2k} 4*5^(r-1-floor(r/2))
       =2*(5^k-1).

All orbits have odd length, so none equals its negative. Consequently:

    # g-invariant labelled Cayley tournaments on (R_k,+)
       =2^(5^k-1).

This count is on 5^(2k) vertices, not 5^k vertices, and is not an
isomorphism-class count. The staircase in ramified valuation is the new
coordinate absent from U1's content-only unramified shell formula.

Independent quick integer-pair orbit enumeration performed during recovery:

| k | nonzero elements | orbit lengths : orbit counts | total orbits |
|---|---:|---|---:|
|1|24|1:4, 5:4|8|
|2|624|1:4, 5:24, 25:20|48|
|3|15624|1:4, 5:24, 25:120, 125:100|248|

Action in basis (1,phi): (a,b) -> (2a+3b,3a+5b) mod 5^k.
The proof and independent saved audit should be reviewed before publication.

## 3. Other scoped extension candidates

1. For general split units (u,v), classify orbits by TWO coordinate contents;
   lengths are lcm of the corresponding component orders. This precisely
   repairs U1's explicit counterexample (2,1) mod7, whose primitive-element
   orbit lengths differ. Standard arithmetic, potentially useful new synthesis.
2. For canonical Fibonacci digit words, investigate the affine block map
   x -> M^ell x+b_word modulo m together with legality memory. Fixed points
   solve (I-M^ell)x=b_word; kernels alone do not establish solvability.
   This can merge word/necklace symmetry with the existing two-register
   automaton, but cyclic boundary legality and labels require separate checks.
3. For cat-map permutation representations, refine order to cycle index,
   then derive invariant coloring counts. Existing modular-data order census
   alone does not give a cycle index. Do not extend permutation conclusions
   directly to arbitrary non-permutation modular representations.

## 4. Publication controls

- Preserve prior credits and KNOWN / PROVED / FINITE / OPEN / ANALOGY labels.
- Distinguish full quadratic ring, a split scalar factor, and projective quotient.
- Preserve exceptional p=2 and ramified p=5 cases explicitly.
- No newly claimed Wall–Sun–Sun existence/nonexistence result follows here.
- No current search-record assertion was externally checked in this recovery.

## 5. Current-source update and general ramified metallic theorem

Parent supplied isolated `session109_source` at
`0db004b11ca0ef8a52986a04abac89f23b12b5da`. Git-object inspection confirms
Session103/107/108 artifacts there (the initially requested scratch checkout
was absent, not the published research). Read current Session108 README at
`05-knowledge/results/session108_sigma_symmetry_20261030/README.md`.
Its relevant inheritance is the distinction between a kernel-preserving
symmetry and an observer-preserving symmetry. It explicitly distinguishes
sigma-squareclass, numerical sigma, and aliquot dynamics, just as this
Pisano extension must distinguish a unit's action, ring reduction, and
isomorphism classes. Session108's exact recursion machinery is handled by
other Session109 contributors; no wholesale rediscovery is claimed here.

**Parent-proposed generalization, independently checked here.** Let
mu²−a mu−1=0, D=a²+4, and p be an odd prime with v_p(D)=1. Put
O=Z_p[mu], delta=2mu−a, and R_k=O/(p^k). Then delta²=D, so O is a
ramified quadratic DVR with v_delta(p)=2 and residue field F_p.
Furthermore p∤a and a²≡−4 modp, hence p≡1 mod4 and p>=5.

The exact identity

    mu^4−1=a mu² delta

shows v_delta(g−1)=1 for g=mu^4. The pth-power binomial expansion gives
v_delta((1+u)^p−1)=v_delta(u)+2 for v_delta(u)>=1: pu has that valuation,
the intermediate terms have larger valuation, and u^p does as well since
p>=5. Thus v_delta(g^j−1)=1+2v_p(j). The nonzero layer of remaining depth r
has (p−1)p^(r−1) elements and orbit length p^floor(r/2). Therefore

    N_k=2(p^k−1),   tournament count=2^(p^k−1).

Equivalently the orientation-bit count B_k=p^k−1 satisfies
B_k=p B_(k−1)+(p−1), with B_0=0. All orbit lengths are powers of odd p,
so antipodal pairing is free. The additive vertex set has p^(2k) elements.
The hypothesis v_p(D)=1 is essential to this proof: nonmaximal orders at
higher discriminant valuation are not silently included. The next clean
example after a=1,p=5 is a=3,p=13.

## 6. THM-491 correction lineage checked at current HEAD

The Git object for THM-491 at `0db004b11` retains the original unqualified
prime-power formula and has no warning. Path history reports commit
`881a9903` (`codex-tiling: prove finite-horizon scale bounds and accelerate
tournament counts`). A targeted current `01-canon/MISTAKES.md` search for
THM-491, Pisano lifting, and Wall–Sun–Sun found the later MISTAKE-614/615
credit/prior-art corrections, but no THM-491 lifting warning. U1 already
contains the correct delay coordinate but does not cite or repair THM-491.
No duplicate correction was located by these targeted searches.

The proposed narrowly scoped correction is: for odd prime p and
h=min entrywise p-adic valuations of Q^pi(p)−I,

    pi(p^k)=pi(p) p^max(0,k−h).

The original h=1 specialization is not a proved universal assertion;
Wall–Sun–Sun exceptions remain an open issue. Say **unproved without the
no-delay hypothesis**, not **refuted**. At p=3, h=1, so the concrete tower
8,24,72 remains correct. This review does not audit unrelated LRC claims
in THM-491. Its old LRC status assertions should not be promoted as current
facts in the new session.

**Search limitation.** A broad current-source Git grep over the Session10*
artifact paths triggered a partial-clone pack fetch, hung without results,
and was interrupted. Git reported RPC/SSL unexpected EOF and `expected
'packfile'` (exit130 after interruption). Therefore the current-source
novelty check is targeted, not exhaustive. Locally available sources and the
specific current THM-491 and Session108 README Git objects were read
successfully. No absence claim relies on the failed broad grep.
