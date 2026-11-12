# All-order receiver inversion and a renewable eight-count bank

**Session282. PROVED from the inherited rooted-fusion carrier, subject to final independent review. Chan–Pak remains OPEN.**

## 1. Inheritance and new operation universe

The starting mechanism is [Session281 receiver completion](../../session281_chanpak_20261112/receiver/PROOF.md): actual tournaments sharing `(h,V)` can be normalized to a fully renewable bank. Its hostile is Session280's widened output interval with unequal V, which common fusion sandwiches cannot repair. The least-used sidecar is the exact endpoint-weight system of a root-deleted core.

We instead prescribe V and solve the unbounded padding variables backwards. This is a signed digit-carry inverse which supplies literal witnesses, not a modular existence test. The carrier `(H,h,Q,U,V,d,e,k)` and physical fusion (identify roots; all nonroot left vertices beat all nonroot right vertices) are inherited from [Session264](../../session264_constructive_20261111/fusion/PROOF.md). A singleton is the identity with all entries1.

Fix any nonempty finite deletion tournament D. Let h=H(D)>0 be odd and let Q be its ordered two-path-cover count, including empty slots. For each incoming root cut x, let C_x be its rooted packet and s_x=V(C_x)−h. An incoming bit means the corresponding vertex beats the root. Endpoint-cover positivity gives s_x≥0, with equality exactly at the empty cut; this is Session280's path-splitting injection.

Our universe initially consists of every rooted tournament with deletion

    TT_a -> D -> TT_b,       a,b≥0,

with arbitrary root incidence. It is not the universe of every deletion-count-h tournament unless a separate core classification is supplied.

## 2. Exact finite inverse of the all-order fixed-V fiber

Let p∈[0,2^a), q∈[0,2^b) encode incoming singleton cuts in increasing-position least-significant-bit order. The ordinal endpoint-weight law, equivalently rooted fusion, gives

    V=h(1+p)+2^a(s_x+Qq).                              (1)

Fix V≥h. Start z_0=V−h,p_0=0. At state a accept every core cut x with

    z_a≥s_x,   z_a≡s_x modQ,

and output the literal witness with q=(z_a−s_x)/Q and b=bitlength(q), b=0 for q=0. Continue by

    epsilon_a=z_a mod2,
    p_(a+1)=p_a+epsilon_a 2^a,
    z_(a+1)=(z_a−h epsilon_a)/2.                       (2)

Stop after processing zero or when z becomes negative.

**Theorem1.** This returns exactly the scalar Hamilton-count support at fixed V in the stated all-order universe, with actual tournament witnesses.

**Proof.** Inductively z_a=(V−h(1+p_a))/2^a and p_a is the unique residue in[0,2^a) solving hp_a≡V−h mod2^a. Oddness of h makes the next bit unique. Thus every solution of(1) is tested, while each acceptance solves it with nonnegative q and a physical core cut. Sufficiency is literal construction, not merely arithmetic consistency.

Extra zero bits after the minimal suffix are global sinks in ordinal order and do not change H. For positive z, z_next≤z/2, so at most floor(log2(V−h))+2 states are processed when V>h. Negative z never returns. At zero, only s_x=q=0 can be accepted: the core root is a source and the count is h times the prefix count. Later states only append incoming-zero prefix vertices; that prefix H does not change. Stopping at zero therefore loses no scalar count, including the singleton/empty-padding boundary. QED.

There are O(2^|D| logV) structural candidates after the core packets are precomputed. Rebuilding all binary packets directly costs O(2^|D|(logV)^2) scalar arithmetic operations; integer bit-complexity is not declared constant. The discarded padding can still change order,Q,U,d and isomorphism type. This is not an enumeration of all tournament classes.

A finite union of core lists with the same h is handled exactly by taking their union of outputs. It remains a named architecture unless the core list is proved complete.

## 3. A sharp GLOBAL obstruction at deletion count3

For D=C3, h=3,Q=12 and every incoming endpoint weight is3. Hence V=3(r+1) for an integer r≥0. Let b_i be the bits of r, zero beyond its end, and set b_(−1)=1. Let B(r) be the Hamilton count of the rooted transitive word with incoming integer r; trailing incoming-zero bits leave B unchanged.

At triangle position a, the incoming triangle cardinality is k=b_a+2b_(a+1). Put t=r>>(a+2), p=b_(a−1), f=3t+2p−(t mod2). The EXACT all-order support is

    {3B(r)+delta_a:a≥0},
    delta_a=0       for k0 or3,
    delta_a=f       for k1,
    delta_a=−f−2    for k2.                           (3)

To derive this, replace two adjacent binary singleton components by a triangle. Subtract three times the two-singleton matrix. In coordinates `(H,h,Q,U,V,d,e,k)` the difference is zero for k0,3; it is `(2,0,0,3,0,1,2,0)` for k1 and `(−4,0,0,−3,0,−1,−2,0)` for k2. Multiplication by actual prefix/suffix matrices gives(3). Prefix d is its last incoming bit, except d=1 for the singleton prefix; suffix V=1+t,e=1−(t mod2), including empty suffix. These boundary entries are essential.

**Theorem2.** Every fixed(h,V) fiber with h=3 contains at most TWO consecutive odd Hamilton counts. The bound is attained by H3,5 at V6.

**Binary no-three proof.** Positive defects arise from10 transitions (least-significant digit first); negative defects arise from01. For two10 transitions a<b, let d=b−a≥2 and write

    t_a=2^d t_b+2^(d−2)+L, 0≤L<2^(d−2).

Substitution into f shows that t_b≥1 forces a difference greater than2. When t_b=0, difference2 is possible only at d2,L0, yielding the pair0,2; for d≥3 the difference is at least4. Thus positive defects have no step-two adjacent pair except0,2. For two01 transitions the corresponding decomposition has2^(d−1) in place of2^(d−2), and the f difference is at least4. Distinct negative defects are never step-two adjacent.

Finally delta=−2 requires a top01 transition with t=p=0: the highest1 has two zeros immediately below it. A positive defect2 requires either the highest10 transition with predecessor1 or a10 transition two positions below the highest1. Both contradict those zeros. Therefore−2,+2 cannot coexist. These facts exclude every step-two triple, including those crossing0. The base defect0 always exists by placing the triangle beyond all nonzero bits. QED.

**Why this is global.** Tournament SCCs are totally ordered and Hamilton counts multiply across their ordinal sum. Classical Camion gives a Hamilton cycle in every nontrivial strong tournament of order n; breaking its n arcs gives n distinct Hamilton paths. Thus the only strong component with count3 is C3, while count1 components are singletons. Every deletion with H3 is exactly TT_a→C3→TT_b. Camion is a cited dependency, not a new result; see [THM-4866](../../../../01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md).

## 4. Global inverse at h5 and an independent renewable four-bank

For the strong four-vertex D=(2,4,9,3), outgoing row masks, h=5,Q=32 and the incoming endpoint weights are(9,7,5,6). A finite exhaustive check through order5 shows this is the unique strong count5 core up to isomorphism. Camion excludes higher orders; SCC multiplication and primality5 permit only one nontrivial factor. Thus Theorem1 also computes the GLOBAL scalar support for every prescribed(h,V)=(5,V).

At V1906 the complete support is

    {195,741,819,821,823,825}.

The four consecutive values have minimal-padding witnesses `(H;a,p,x,q,b)`:

    (819;7,73,6,0,0), (821;3,1,10,7,3),
    (823;0,0,10,59,6), (825;8,73,8,0,0).

They share h5,V1906 but have e0,0,2,0. Applying C2 (root word100) gives

    H'=H+4V, V'=8V−h, e'=h,

hence renewable8443,8445,8447,8449 with(h,V,e)=(5,15243,5). The inherited paid carry uses z=(15243−13)/5=3046,P6092 and produces eight consecutive odds beginning92,862,629. Repeating C3 (root word1000, exponent4) preserves the carry congruence and gives infinitely many scales. The widened eight-output union has the known receiving defect; it is not recursively renewable. This is not a larger-width record than Section5. A universal maximum fiber-run length4 at h5 is NOT proved.

## 5. Inversion across distinct cores produces a renewable EIGHT-bank

The compiler was applied to distinct deletion-count15 cores, not just one core with a larger size cutoff. Define

    A=(6,12,24,17,3), E=(44,1,2,6,15,30),
    B=(2,4,1), C=(2,4,9,3), BC=B->C.

All these deletion choices A,E,BC have h15. The following explicit minimal-padding witnesses have V139878:

|H|core|a|p|incoming core mask x|q|
|---:|:---|---:|---:|---:|---:|
|54365|A|8|57|7|4|
|54367|A|8|57|11|4|
|54369|BC|10|313|96|0|
|54371|A|2|1|3|291|
|54373|A|2|1|5|291|
|54375|E|4|9|57|40|
|54377|A|0|0|7|1165|
|54379|A|0|0|11|1165|

Take b=bitlength(q). Every graph has17 vertices except54369, which has18. The verifier constructs their row masks and independently checks all receiver coordinates by literal subset DP. Neutral sink padding puts them all at18 if desired.

The common C1 context (root word10) has

    H'=H+2V, V'=4V−h, e'=h.

It therefore produces the renewable EIGHT-bank

    334121,334123,...,334135; (h,V,e)=(15,559497,15).    (4)

Any common left unit-deletion context translates all eight counts by the same amount and retains their common receiver, indefinitely. This is a proof from literal seeds and the inherited action law, not extrapolation of an observed run.

**Load-bearing boundary.** The eight-bank does NOT satisfy the existing width-doubling phase criterion. V139878 mod15 is3; its power-of-two orbit is{3,6,12,9}, excluding width8. Completion does not repair that incompatibility. Selecting six members recovers a permitted six-to-twelve route, but not sixteen. The new bank is wider and genuinely renewable, while eventual target coverage and unbounded renewable width remain OPEN.

The discovery search included A, two six-vertex cores, and B/C in both orders with at most8 internal singleton vertices, followed by the exact all-order outer-padding inverse. It found the displayed eight-bank at V139878. The process searching for nine was interrupted by an environment reset; NO negative conclusion or completed range from that process is reported. Only the displayed positive witnesses are used.

## 6. What changes for the constructive problem

The inverse replaces an unbounded padding search by a finite signed carry descent at a prescribed ACTUAL receiver. It keeps core identity, binary incidence and absolute H. Every accepted state owns a literal witness. It still prescribes V rather than H and does not give a uniform core supplier.

The global h3 obstruction says the most natural low-count architecture cannot supply a three-bank. Unequal endpoint weights at h5 escape that obstruction. Crossing core types at h15 then supplies a new renewable eight-bank. The missing step is not another congruence check: it is a construction controlling both the absolute H target and the continuing receiver while increasing renewable width.

Run `python3 verify.py` and `python3 -O verify.py`. Literal controls and the positive all-scale constructions are separated in RESULTS.json. The mathematical all-order guarantees come from the proofs, not finite census extrapolation.
