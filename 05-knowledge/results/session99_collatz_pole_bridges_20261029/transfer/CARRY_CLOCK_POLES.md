# Fixed-size clock poles and the terminal carry defect

Status: PROVED below; finite exact checks in `check_carry_clock_poles.py`.
This is a refinement of, not a replacement for, the inherited affine moment recurrence.

## Source and conventions

Pinned source: `d9d7e33349653b6b65dbad477841c9ed518c081a`.
Inspected: `05-knowledge/results/affine_tree_recurrences_20261028/{RECURRENCES,MOMENTS}.md`,
`session95_niche_recurrences_20261029/valuation/RESIDUE_CLOCK_BRANCHING.md`, and
`session97_overnight_arithmetic_carriers_20261029/controllers/{README,INTEGER_TERMINATION,SOURCE_PRESERVING_RETRIES,FAILURE_TRANSDUCER}.md`
(the last two directory names are under `05-knowledge/results/`).
We use the cheaper (3,3) rule, not the critical parity rule: a valuation divisible by 3 gives queue increment +2, every other valuation gives -1. Start at 1 and stop first at 0.

For k branches put J=3k+1 and c_k=binom(3k+1,k)/(3k+1).
Let W_k be all chronological first-return valuation words with k branches.
For w=(a_1,...,a_J), D=sum a_i and

    C(w)=sum_{i=1}^J 3^(J-i) 2^(a_1+...+a_(i-1)).

Define the *counting* clock series H_(k,r)(z)=sum_(w in W_k) C(w)^r z^D,
with integer moment r>=0. Haar-weighted clock series is H_(k,r)(z/2).
These are raw integer carries, not normalized intercepts C/2^D.

## 1. Exact fixed-size clock law

Set L(z)=(z+z^2)/(1-z^3), B(z)=z^3/(1-z^3). Then

    H_(k,0)(z)=c_k B(z)^k L(z)^(2k+1)
             =c_k z^(5k+1)(1+z)^(2k+1)/(1-z^3)^(3k+1).

There are c_k legal ordered full ternary shapes; each contributes k branch slots and 2k+1 leaf slots. Consequently all three cube roots of unity are poles of exact order 3k+1: none makes the numerator zero. This is the inherited residue-clock coefficient law specialized to the subcritical cheaper controller, with its poles made explicit.

In particular a 3k+1 pole order arises here because there are 3k+1 independently summable valuation slots. This is a structural parallel to other rational enumeration series, not identification with a tournament transfer operator or with the arithmetic expression 3n+1.

## 2. Mark profiles give a finite rational carry algorithm

For each legal shape sigma and each weak composition m_1+...+m_J=r, let

    h_j=sum_(i>j) m_i,     A(m)=r!/prod(m_i!) * 3^(sum_i (J-i)m_i).

Let V_j be B at a branch slot and L at a leaf slot. Then the exact finite formula is

    H_(k,r)(z)=sum_sigma sum_m A(m) prod_(j=1)^J V_j(2^(h_j)z).    (1)

Proof: expand the rth power of the displayed carry sum multinomially. A valuation a_j occurs in precisely those summands marked at i>j, hence exponent h_j. Sum each valuation in its allowed residue class independently. This is finite over shapes and marks, even though valuations are unbounded. It retains the chronology before the final summation.

The heights satisfy r>=h_1>=...>=h_J=0. Thus a denominator bound is

    (1-z^3)^J prod_(e=1)^r (1-2^(3e)z^3)^(J-1).                 (2)

This is an annihilating denominator, not a claim of exact multiplicity at every nonreal pole.

## 3. Sharp positive-real pole theorem

For k>=1 and integer r>=0:

* z=1 has exact pole order J=3k+1;
* for each 1<=e<=r, z=2^(-e) has exact pole order J-1=3k;
* no other positive-real pole occurs.

At z=1, order J requires h_j=0 for all j, so all r marks are at i=1. The unique maximal-order contribution is

    3^(r(J-1)) H_(k,0)(z).

Its leading coefficient is nonzero. Other mark profiles have fewer zero-height slots.

At z=2^(-e), e>0, a profile has at most J-1 slots at height e because h_J=0. Exactly J-1 requires m_1=r-e, m_J=e and every intermediate m_i=0. Its contribution summed over shapes is

    c_k binom(r,e) 3^((J-1)(r-e))
      * B(2^e z)^k L(2^e z)^(2k) L(z).                       (3)

Every complete shape ends in a leaf, and the preceding J-1 positions have k branches and 2k leaves. Multiplying (3) by (1-2^(3e)z^3)^(J-1) and setting z=2^(-e) gives

    c_k binom(r,e) 3^((J-1)(r-e)) 2^(2k) L(2^(-e)) > 0.

All remaining profiles have lower pole order there. Cancellation of the top term is therefore impossible. At z=1 the analogous normalized leading coefficient is

    3^(r(J-1)) c_k 2^(2k+1) > 0.

For k=0 there is just one leaf, C=1, and H_(0,r)=L for all r; no shifted carry pole occurs. Haar weighting moves the positive poles to z=2^(1-e), with unchanged orders. Normalized-intercept moments instead replace z by z/2^(r+1), and must not be confused with either convention.

The missing one in 3k is arithmetic: the terminal valuation divides the endpoint but cannot affect the numerator carry accumulated before that final division. This distinguishes the full clock from the affine extension coordinate.

### Moment-boundary corollary

For every fixed k>=1 and every integer r>=1,

    E[C^r ; K=k] = infinity,
    E[(C/2^D)^r ; K=k] < infinity.

Indeed the positive-coefficient counting series has radius 2^(-r): its smallest positive pole is there, and all candidate poles have at least that modulus. Haar evaluation for raw carries is at z=1/2. For r=1 this is the genuine pole itself; for r>1 it is beyond the convergence radius. In contrast the normalized-intercept moment evaluates at z=2^(-(r+1)), strictly inside the radius. The event K=k has strictly positive finite mass, so the same distinction holds for conditional moments given K=k. There is also the elementary check C/2^D <= sum_(i=1)^J 3^(J-i)/2^(J-i+1) = (3/2)^J-1, because all terminal suffix valuations are at least one. At k=0 the carry is identically one. Thus summing a finite formal raw-carry jet and then evaluating at Haar mass is not automatically a finite numerical statistic. This does not contradict the inherited finite unconditional mean of the normalized carry for the cheaper rule.

## 4. Source guard and nonidentification

This computation is an aggregate clock law, not a source-complete quotient. Each underlying word still requires

    n ≡ (2^D-C) 3^(-J) (mod 2^(D+1)),
    2^D n_out=3^J n+C.

Concatenation must use C_uv=3^(J_v)C_u+2^(D_u)C_v and pull back the second source guard. The words (3,1,1,2) and (3,1,2,1) are both cheaper first returns with J=4,D=7, but carries 179 and 211; their least positive odd source residues are 189 and 157 modulo 256. Clock poles alone cannot distinguish them. The checker verifies guards by native execution.

The shifted poles come from scalar substitutions attached to carry marks. A pole's multiplicity in a reduced scalar generating function records a reachable and observable generalized eigenvalue length in a finite linear realization, not merely the number of matrix eigenvalues. No operator conjugacy, shared Jordan blocks, or equality with tournament poles 4^k follows from a numerical resemblance. The present theorem neither proves nor needs native Collatz convergence; Session97's universal termination applies to each auxiliary first return, not to unbounded source-preserving retries.

## 5. Reproduction receipt

Run `python3 -B check_carry_clock_poles.py` and `python3 -O -B check_carry_clock_poles.py`.
Both completed successfully and produced byte-identical JSON in this session (`checks.json`, `checks_optimized.json`).
The checker requires SymPy (run with version 1.14.0), uses explicit exceptions rather than assertions, verifies eight (k,r) cases with 64 independently enumerated clock coefficients, all stated positive-real leading coefficients in these cases, the k=0 exception, and the two native source guards. The all-k, all-r statements are supplied by the proof, not inferred from those bounded tests.
