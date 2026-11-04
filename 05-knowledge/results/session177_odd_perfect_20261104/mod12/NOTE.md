# Mod-12 sigma transfers realize every finite tournament

Status: **PROVED; independent mathematical and arithmetic audit PASS.** Classical inputs are quadratic reciprocity, the Chinese remainder theorem, and Dirichlet's theorem. No literature-priority claim. Odd perfect number existence remains **OPEN**.

## Main theorem: common-exponent universality

For every finite labelled tournament T, there exist distinct primes p_1,...,p_n, all congruent to 11 modulo 12, and a single positive even integer a such that, for all i != j,

    i -> j in T  iff  p_j divides sigma(p_i^a).

The exponent a can additionally be required to be 0 modulo 4, or instead 2 modulo 4. For any prescribed positive integer H, all displayed arc valuations can simultaneously be required to be at least H.

These are induced divisibility graphs on the selected primes. They are NOT claims that the selected set contains all prime divisors of the sigma values, or that the resulting product is perfect.

### Proof

First realize T as a quadratic-symbol tournament. Suppose p_1,...,p_k have been selected. For each i<=k choose a nonzero residue c_i modulo p_i whose Legendre symbol (c_i/p_i) is -1 if i->k+1 is desired, and +1 otherwise. Both signs exist. CRT gives a reduced residue class c modulo 12 p_1...p_k with c=11 modulo 12 and c=c_i modulo p_i. Dirichlet gives a new prime p_{k+1} in that class, larger than the preceding labels. Quadratic reciprocity, since all labels are 3 modulo 4, gives

    (p_i/p_{k+1}) = -(p_{k+1}/p_i).

Thus (p_i/p_j)=+1 precisely along the chosen arcs.

For an arc p->r, r=3 modulo 4, so its quadratic-residue subgroup has odd order (r-1)/2. Consequently d=ord_r(p) is odd. Define the odd edge requirement

    t(p,r) = r if p=1 modulo r;
             d otherwise.

If an odd m is divisible by t(p,r), then r divides 1+p+...+p^(m-1). The first case needs r|m; merely taking d=1 would be incorrect because the geometric-series denominator vanishes modulo r. In the second case the denominator is invertible and p^m=1 modulo r.

Choose a common odd m>=3 divisible by every t(p,r), and set a=m-1. Every prescribed arc appears. A reversed arc is impossible: if r divides sigma(p^a), then r divides p^(a+1)-1, so ord_r(p) divides the odd a+1 and (p/r)=+1. Reciprocity forbids the opposite quadratic symbol. This proves the iff.

Multiplying m by an odd factor 1 or3 makes m=1 or3 modulo 4 as desired, yielding the two exponent classes. To force valuations at least H, choose m divisible by every d*r^H for d>1, and r^H when d=1. The elementary odd-prime lifting formula yields the claimed lower bounds. Again an odd multiplier handles the modulo 4 condition. QED.

## Exact local gate and valuations

For distinct odd primes p,r and positive even a, put m=a+1. Then:

- If p=1 modulo r, v_r(sigma(p^a))=v_r(m).
- If p!=1 modulo r and d=ord_r(p) does not divide m, the valuation is zero.
- Otherwise v_r(sigma(p^a))=v_r(p^d-1)+v_r(m/d).

The last two identities follow directly from the geometric sum and odd-prime LTE. In particular, possible edges into an r=3 modulo 4 target are exactly its quadratic residues, but actual edges retain the exponent-divisibility condition. The valuation sidecar cannot be replaced by an unweighted tournament. A hostile denominator control is p=23, r=11: ord_11(23)=1, but sigma(23^2)=553 is not divisible by 11. The first positive even exponent producing that edge is 10, because 11 must divide a+1.

## Why 12 matters, and why 11 and 13 differ

For every prime p>3 and e>=1, p^2=1 modulo 12, hence

    sigma(p^(2e)) = 1+e(p+1) modulo 12.

The four residue rows p=1,5,7,11 give respectively 1+2e, 1+6e, 1+8e, and 1. Thus every p=11 modulo 12 source outputs a number equal to1 modulo 12, regardless of e. Its prime factors are all coprime to 6. If c_r is the valuation-weighted number of its factors in residue class r modulo 12, then

    c_5+c_11 = 0 modulo2,
    c_7+c_11 = 0 modulo2.

Equivalently c_5,c_7,c_11 have identical parity. These two character checks are exact, but do not constrain tournament topology by themselves.

For a hypothetical 3-free odd perfect number, Euler's special prime is 1 modulo 12. Every nonspecial prime p=1 modulo 3 has even exponent a=0 or4 modulo 6, since a=2 modulo 6 makes 3 divide sigma(p^a). The special exponent is 1 or 9 modulo 12, not 5. Concrete contrast:

    sigma(11^2)=133=7*19;
    sigma(13^2)=183=3*61.

Consequently the exact component 11^2 in a 3-free candidate forces both 7 and 19, each to exponent at least 4. The component 13^2 is forbidden in that branch. These are local lower costs, not exclusions of 11 or 13 as prime divisors. The special component 13^1 is permitted locally and forces 7 through sigma(13)=14.

## Hostile cycle and the missing information

The exact directed cycle

    11^40 ->83, 83^88 ->179, 179^4 ->11

has multiplicative orders 41,89,5. Each displayed target occurs to valuation exactly 1; every reverse arc is forbidden for every even exponent. Each sigma value is 1 modulo 12. Dividing by its displayed target leaves a cofactor 11 modulo 12, so each source necessarily has additional prime factors outside this three-vertex set, including at least one 3 modulo 4 factor to odd valuation. This is a genuine cycle in local transfers, not a closed sigma component and not an odd-perfect candidate.

For odd perfectness one must restore an Euler prime q=1 modulo 4, the exceptional odd exponent at q, closure under every sigma factor, exact incoming valuation budgets, and sigma(N)=2N. The companion odd-order reachability result requires every prime to reach q; arbitrary tournament universality does not invalidate that necessary global condition. It instead proves that an induced nonspecial tournament cannot be excluded merely by its finite shape, even with all its primes in the single 11 modulo 12 class and a common even exponent.

## Repository bridge and open next target

Source: exponent-labelled divisor-transfer networks. Target: intrinsic tournaments on 3 modulo 4 prime vertices. Map: p->r iff the quadratic symbol is+1, followed by exponent requirements. Preserved: allowed divisibility direction. Lost: all valuations, exponent congruences, exterior prime factors, Euler exception, and abundancy. Sidecar needed: exact sigma factorization plus incoming/outgoing valuation accounting. This is not the usual Paley tournament on differences in one finite field.

The cheapest next decisive test is to combine 3-free exponent costs with the companion unique-Euler-sink filter and a valuation budget. For each sink candidate q, can one certify that every realization of the required return paths introduces more valuation demand than its suppliers can pay? This is **OPEN** here. A theorem using only finite induced tournament exclusions cannot accomplish it, by the universality theorem. Arbitrary realizable arc lower bounds also do not supply exact budget balance.

Inheritance boundary: the repository reflection `07-reflections/two-twelves-the-universal-bernoulli-and-the-arithmetic-runner-and-why-they-coincide-at-n14-kps.md` distinguishes an analytic 1/12 from an arithmetic 12; neither identifies those objects structurally. The present 12 is specifically the unit group of Z/12Z and its two quadratic characters.

## Exact computation and controls

Run `python3 check.py` in this directory. Standard library only. It constructs and checks all 1099 labelled tournaments on 1 through 5 vertices, with prime labels at most 3023, and 5580 direct mod12 geometric-sum checks. The chosen per-source exponents in the finite checker suffice to test the basic theorem; their common odd-length least common multiple independently tests the stronger common-exponent version. No integer factorization of enormous sigma values is needed: geometric sums are tested exactly modulo target primes. Negative controls check every reversed pair. The explicit cycle and its valuations are separately checked.

These finite checks support implementation correctness; the proof above handles all finite tournaments. They are not an exhaustive search for odd perfect numbers.

## Sources checked 2026-11-04

- Andrew V. Sutherland, MIT 18.785 Number Theory I, Fall 2021, Dirichlet theorem: https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_full_lec.pdf
- Abhinav Kumar, MIT 18.781, Quadratic Residues and Quadratic Reciprocity: https://ocw.mit.edu/courses/18-781-theory-of-numbers-spring-2012/resources/mit18_781s12_lec9/
- S. Adam Fletcher, Pace P. Nielsen, Pascal Ochem, *Sieve Methods for Odd Perfect Numbers*, factor-chain context and Euler form: https://mathdept.byu.edu/~pace/OPNSieves_web.pdf

The last paper is a historical primary reference, not a source for current record bounds. The universality theorem is derived here; searches did not establish literature novelty.

## Independent audit receipt

Independent proof review checked the CRT reciprocity signs, positive even common exponent, both modulo 4 exponent classes, arbitrary valuation lower bounds, local 3-free OPN gates, and the three-cycle exterior-factor claim. No mathematical defect found. `independent_audit.py` uses quotient residues modulo r(p−1), avoiding the primary checker's modular-inverse calculation. Across all 1099 labelled tournaments it passes 21300 per-source incidence checks, 42600 common-exponent/modulo 4 checks, and 42600 H=2 valuation/modulo 4 checks. Both checkers use explicit checks active under Python -O; normal and optimized runs pass. General H and arbitrary tournament size are proved in the note rather than inferred from finite tests.
