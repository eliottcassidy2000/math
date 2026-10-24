# Primorial offsets, exact finite membership, and the Ogg comparison

24 October 2026. **PROVED** elementary wheel, cutoff, and logical-composition statements; **FINITE-EXACT** certified arithmetic through the value cutoff 71; **CITED** supersingular counting formulas and Ogg's classification; **OPEN** general Fortune and lesser-Fortune conjectures. No new supersingular correspondence, prime-generation theorem, native hole closure, or LRC result is claimed.

## 1. Inheritance and the question being answered

The closest proved mechanism is the elementary finite-wheel rejection argument in `05-knowledge/results/abc_difference_families_20261011c.md`, Section 2. Its canonical hostile is an arithmetic progression that survives a finite wheel while both intended endpoints have proper factors. That note already warns that a finite sieve is not a prime oracle.

The supersingular comparison inherits `05-knowledge/results/monster_levels_square_sum_7_21_20261022.md`, Sections 1–2, and `supersingular_divisor_shapes_qseries_20261022.md`, Section 2. The corrected level filtration is genuine: levels 1, 2, and 6 recover the stated Monster/Baby Monster/Fi22 prime sets. The present offset rules are different predicates and do not inherit this geometric interpretation. MISTAKE-666 and MISTAKE-668 are preserved boundaries, not superseded by numerical overlap.

The least-used relevant sidecar is the actual branch-conditioned support in Session49, `session49_relational_hole_repairs_20261024/hole-theory/THEOREM.md`, Sections 2–3. Its exact hole intersection requires all compatible full records, not a necessary residue screen.

Live concepts: primorial-indexed offsets; offset value sets; prime certificates; Frobenius rationality; conditional support versus residue survival. Their carriers and quantifiers are compared below.

## 2. The four supplied sequences have different indexing rules

Let p_n be the nth prime, starting p_1=2, and let P_n=product_{i=1}^n p_i. Thus P_6=13#=30030. Here 13# means the product of primes at most 13, not the product of the first 13 primes.

Define

    F_n = min{d in Z : d>1 and P_n+d is prime}, n>=1;
    L_n = min{d in Z : d>1 and P_n-d is prime}, n>=2.

Equivalently, F_n is the difference between P_n and the first prime strictly above P_n+1; L_n is the difference between P_n and the last prime strictly below P_n-1. The strict inequalities are essential. L_1 is undefined because 2-d is not prime for d>1.

- `3,5,7,13,23,17,19,23,37,61,67,61,71,47,...` is F_n, in **primorial-index order**, OEIS [A005235](https://oeis.org/A005235).
- `7,11,13,17,29,23,43,41,73,...` is L_n starting at **n=3**, omitting L_2=3. The complete sequence is [A055211](https://oeis.org/A055211), whose indexing starts at 2.
- `3,5,7,13,17,19,23,37,47,59,61,67,...` is the **increasing distinct-value list** of Fortunate primes, [A046066](https://oeis.org/A046066). This sorting discards n and repeated occurrences.
- `11,29,31,41,43,...` is the list of odd primes that never occur as F_n. Including 2 gives [A045493](https://oeis.org/A045493): `2,11,29,31,41,43,53,...`.

The name “Fortunate primes” is conventional; asserting that **every** F_n is prime is Fortune's conjecture, not the definition or a proved theorem. General primality of L_n is likewise conjectural. The primary article [Čejchan–Křížek–Somer, JIS 25 (2022), Section 3](https://cs.uwaterloo.ca/journals/JIS/VOL25/Krizek/krizek3.pdf) explicitly distinguishes these conjectures from short-window sufficient conditions. No claimed large verification limit is used here.

### Exact initial table

Every entry below is certified by `check.py` and `certificates.json`, including all smaller forbidden offsets, not merely copied from a sequence database.

| n | p_n | F_n | L_n |
|---:|---:|---:|---:|
|1|2|3|undefined|
|2|3|5|3|
|3|5|7|7|
|4|7|13|11|
|5|11|23|13|
|6|13|17|17|
|7|17|19|29|
|8|19|23|23|
|9|23|37|43|
|10|29|61|41|
|11|31|67|73|
|12|37|61|59|
|13|41|71|47|
|14|43|47|89|
|15|47|107|67|
|16|53|59|73|
|17|59|61|107|
|18|61|109|89|
|19|67|89|101|

## 3. The elementary mechanism and its sharp boundary

### Proposition 1. Wheel transfer

Fix n and an integer d>1. If P_n+d or P_n-d is a prime greater than p_n, then

    gcd(d,P_n)=1.

**Proof.** Every prime r dividing gcd(d,P_n) also divides P_n+d and P_n-d. Since r<=p_n, it would be a proper factor of the specified prime endpoint. QED.

Consequently d>=p_{n+1}. If additionally d<p_{n+1}^2, then d is prime: a composite d in this interval would have a prime factor at most sqrt(d)<p_{n+1}, hence dividing P_n. This is a one-way implication from an **already prime endpoint** to the small offset's primality, not a way to make the endpoint prime.

The endpoint condition greater than p_n is automatic for positive offsets. It also holds for the endpoint defining L_n when n>=3: Bertrand's postulate gives a prime q with p_n<q<2p_n, and P_n>=6p_n, so q<P_n-1. The greatest prime below P_n-1 is therefore greater than p_n. The exception L_2=3 has gcd(3,6)=3 and must not be silently fed into the coprimality lemma.

Thus F_n>=p_{n+1} for all n>=1, and L_n>=p_{n+1} for n>=3. A composite value in either regular range must be at least p_{n+1}^2. This does not prove that the first prime endpoint is always reached before that square threshold.

### Corollary 2. Exact finite membership, without either conjecture

For any fixed integer B>=3, every occurrence F_n<=B has p_n<B. The same holds for L_n<=B except that L_2=3 is handled separately. Hence checking all p_n<B settles membership **for all indices**, not just for an observed prefix.

For B=71, only n=1,...,19 can contribute. For a fixed prime q, the positive-offset predicate is the finite statement

    q is prime AND there exists n with p_n<q such that
    P_n+q is prime AND every P_n+d, 2<=d<q, is nonprime.

This is an explicit arithmetic predicate. It includes a minimality condition and is not equivalent to “q is prime” or “q passes the wheel.” The negative predicate replaces + by -, with its initial exception retained. The existence of a finite test for each proposed value does not settle universal primality of all F_n.

### The user's 30030 and 17 claim: verified, and strengthened

Every d=2,...,16 has a prime divisor among 2,3,5,7,11,13. Both 30030-d and 30030+d exceed 13, so both are composite. At d=17,

    30030-17=30013 and 30030+17=30047

are prime. Thus F_6=L_6=17, and 17 is the smallest offset greater than 1 giving a prime on either side, hence also the smallest giving primes on both sides. Their gap is 34; they are not a twin-prime pair.

There is a sharp, larger regularity window: if 1<d<289 and either endpoint 30030±d is prime, then d itself is prime. The threshold cannot be enlarged to include 289:

    30030-17^2=29741 and 30030+17^2=30319

are both prime, although 17^2 is composite. Independent trial division certifies all four displayed endpoints. The paper cited above also records the endpoint examples; the bidirectional wheel proof here states its exact size condition explicitly.

For a specified sign, the reverse implication fails well inside the window:

    d=19 is prime, gcd(19,30030)=1,
    30030-19=30011 is prime,
    30030+19=30049=151*199 is composite.

Even the weaker reverse implication “prime offset implies at least one prime endpoint” fails: d=23 is prime and coprime to 30030, but 30007=37*811 and 30053=41*733 are both composite.

The offset-1 boundary also matters: 30029 is prime, while 30031=59*509 is composite. The lesser Fortunate definition deliberately skips 30029.

### Proposition 3. The two signs share the same small-prime wheel

For every r dividing P_n, the congruences P_n+d=0 mod r and P_n-d=0 mod r are both equivalent to d=0 mod r. Thus demanding small-prime survival on both sides does **not** square the single-side survival fraction. The survivor residues are exactly the units modulo P_n, with cardinality product_{r|P_n}(r-1). For P_6 they number 5760, with fraction 192/1001.

This is an exact finite residue count, not an asymptotic prime-pair density. Beyond the included primes the signs need not have the same obstruction, as d=19 shows. More generally, around an arbitrary center a, a prime r forbids offset residues {a,-a}; these coincide exactly when r divides 2a. This explains precisely why a primorial center synchronizes the two local screens.

## 4. What the supersingular comparison actually says

Let O denote Ogg's fifteen primes:

    {2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}.

Their defining geometric predicate is that all supersingular j-invariants in characteristic p are F_p-rational, equivalently genus(X_0^+(p))=0. This is Ogg's theorem, recalled and generalized in [Aricheta, SIGMA 15 (2019), Theorem 1.1](https://arxiv.org/html/1809.07421v2). “Supersingular prime” here means this moonshine condition, not a prime of supersingular reduction for one fixed elliptic curve.

We do not identify O by overlap with a numerical list. For p>3, define

    s(p)=floor(p/12)+epsilon(p),
    epsilon(p)=0,1,1,2 for p=1,5,7,11 mod 12;
    r(p)=(h(-4p)+[p=3 mod 4] h(-p))/2;
    gplus(p)=(s(p)-r(p))/2.

Here h(D) counts proper classes of primitive positive definite binary quadratic forms of discriminant D, including nonfundamental discriminants. The standard supersingular count and class-number formulas identify s with total j-invariants, r with rational ones, and gplus with the number of nonrational Frobenius pairs. See [Delfs–Galbraith, Section 1](https://arxiv.org/pdf/1310.7789) and the inherited audited supersingular note. The class-number formula can equivalently be written r=h(-4p)/2 for p=1 mod4, r=h(-p) for p=7 mod8, and r=2h(-p) for p=3 mod8.

The exact script enumerates reduced forms for every prime through 71. As an independent check it enumerates all j in F_p, counts points on one nonsingular curve of each j, and checks that exactly r(p) have trace zero. At p=2,3 the unique supersingular j is 0, treated separately. Thus the test compares two explicit arithmetic predicates: **primorial first-prime offset** and **vanishing Frobenius-pair count**.

Writing F={F_n} and L={L_n}, the complete value decisions are

    F ∩ [2,71] = {3,5,7,13,17,19,23,37,47,59,61,67,71};
    L ∩ [2,71] = {3,7,11,13,17,23,29,41,43,47,59,67}.

Consequently

    O ∩ F = {3,5,7,13,17,19,23,47,59,71};
    O \ F = {2,11,29,31,41};
    (F ∩ [2,71]) \ O = {37,61,67};
    O \ L = {2,5,19,31,71};
    O \ (F union L) = {2,31}.

All these are global set statements because O is contained in [2,71] and Corollary 2 closes future indices. They do not depend on Fortune's conjecture. The union has additional non-Ogg primes {37,43,61,67} below 72, so even combining signs does not characterize O.

Decisive hostile examples retain arithmetic witnesses:

- 11 is in O but never equals F_n: only n<=4 could contribute, and their F-values are 3,5,7,13.
- 31 is in O but never occurs in either offset sequence; the complete cutoff test excludes it.
- 37=F_9, but h(-148)=2, s(37)=3, r(37)=1, and gplus(37)=1, so 37 is not in O.
- 61=F_10=F_12=F_17 has gplus=1; 67=F_11=L_15 has gplus=2.
- 17=F_6=L_6 has two supersingular j-values, 0 and 8 modulo 17, both rational. The two facts coexist but neither calculation supplies a map between their witnesses.

The inherited level differences {29,41,59,71}, {17,19,23,31,47}, and {2,3,5,7,11,13} are not reproduced by Fortunate membership. Their intersections with F are respectively {59,71}, {17,19,23,47}, and {3,5,7,13}. The original level filtration remains geometric; this cross-tabulation is only an exact comparison.

## 5. Euclid–Mullin is not either Fortunate construction

Euclid–Mullin starts e_1=2 and takes e_{k+1} to be the **least prime factor** of 1+product_{i<=k}e_i. Its first seven terms are 2,3,7,43,13,53,5. Every term is prime by definition; no term repeats because an earlier prime divides the old product and cannot divide that product plus 1. The sequence is therefore infinite, but the proposition that it contains every prime is a different open problem. See [Booker–Irvine, Introduction](https://arxiv.org/html/1508.03039v2).

Fortunate construction multiplies the **first n ordinary primes**, then searches additively for a prime endpoint and reports the offset. Euclid–Mullin multiplies **its own earlier outputs**, factors one number, and reports a factor. Equality of a few small primes does not identify these operations. For example, 2*3*7*43+1=1807=13*139 gives the fifth Euclid–Mullin term 13; this is not the fourth primorial 210 or its Fortunate offset 13.

## 6. The legitimate connection to conditional hole filters

Let I be a finite offset domain on which both endpoints exceed the sieving primes, and define the **actual prime supports**

    A_plus={d in I:P_n+d is prime}, A_minus={d in I:P_n-d is prime}.

For holes H_plus=I\A_plus and H_minus=I\A_minus, elementary complement identities give

    no prime on either side: H_plus intersection H_minus;
    failure of simultaneous primality: H_plus union H_minus.

The existential-versus-conjunctive quantifier determines whether holes intersect or unite. In particular, “prime on at least one side” and “prime on both sides” are not interchangeable tasks.

Replacing actual prime supports by the local wheel W={d:gcd(d,P_n)=1} gives only a **superset** of each A. Failure of W is a sound rejection once endpoint-equals-divisor exceptions are excluded; survival is not a full witness. It forgets all factors outside P_n, endpoint primality, and first-success minimality. The d=19 example supplies an explicit false positive.

Session49's branch identity instead uses A_z={q(x):(x,z) is a fully compatible actual assignment} and therefore exactly gives

    H = intersection_z (c(z)+H_z).

For two branches this is H'=H_0 intersection (lambda+H_1), with branch-specific domains and grades retained. It is exact **because an element of A_z already has a full witness**. It is not made exact by analogy to a prime wheel. Neither a residue survivor nor an unsupported charge label provides a native route. No map from offsets, primality certificates, or supersingular points to those native records has been supplied here.

The transferable lesson is precise: keep the source domain, the AND/OR quantifier, the excluded endpoint exceptions, and the witness sidecar. The identities are elementary; no new native existence theorem is asserted.

## 7. Reproduction and validity ledger

Run from this directory:

    python3 check.py
    python3 -O check.py

Both write the same `results.json` and produce identical saved outputs. The normal verifier uses only the Python standard library. Optional `python3 check.py --generate` rebuilds the proposal data with SymPy, then immediately runs the exact verifier; SymPy probable-primality outputs are never accepted as proof.

`certificates.json` includes 200 recursive Lucas prime certificates, reaching 7858321551080267055879179, plus proper-divisor or strong-Miller–Rabin **failure** witnesses for every composite endpoint before each first prime. A failed strong-Miller–Rabin test is an exact compositeness witness; passing such a test is never treated as proof of primality.

For a prime certificate at N, the verifier checks the full factorization N-1=product q^e with recursively certified prime q, a^(N-1)=1 modulo N, and gcd(a^((N-1)/q)-1,N)=1 for each q. For any prime divisor r of N these conditions force a to have order N-1 modulo r, hence N-1 divides r-1 and r>=N. Thus N is prime. This explains why the check is exact at sizes beyond deterministic 64-bit primality bounds.

The finite universe is all 19 possible primorial indices contributing values at most 71, every prime through 71 for the Ogg comparison, all offsets 2..289 at P_6, and the first seven Euclid–Mullin terms. All endpoint tests around 30030 have a separate trial-division path. Form counting is cross-checked by elliptic-curve point counting. Explicit `require` checks, not Python assertions, remain active under `-O`.

Open boundaries: neither Fortune conjecture is proved; no future prime-gap bound is assumed; no prime density or independence law is inferred; Ogg's all-prime classification is cited rather than re-proved by the bounded computation; no supersingular/Fortunate or native-routing bijection is asserted. No theorem ID is reserved and no literature-priority claim is made.
