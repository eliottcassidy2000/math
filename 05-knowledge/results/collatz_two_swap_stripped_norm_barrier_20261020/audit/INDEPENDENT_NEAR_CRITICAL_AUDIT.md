# Independent audit: an infinite near-critical stripped-norm barrier

Date: 20 October 2026. Verdict: **PASS**.

Scope: independent read-only audit of the proposed upgrade. The frozen two-edit packet and its symbolic files were not changed. This note audits the family construction, actual swap legality, reduction at 2 and 3, limiting argument, and the scope of the resulting obstruction. It does not prove a two-edit cycle exclusion or an actual-content lower bound.

## 1. Exact family and the lower boundary

Let alpha=log_2 3. For each integer a>=2, define

    p=2*3^a,
    q_0=ceil(alpha*p),
    q=q_0+((1-q_0) mod 6),
    r=q-p,
    n=p-r=2p-q.

The residue in the definition of q is the member of {0,1,2,3,4,5}. The following claims all hold:

- 0<q-alpha*p<6. Irrationality of alpha makes 0<q_0-alpha*p<1, and q-q_0<=5.
- gcd(p,q)=1, since the only prime factors of p are 2 and 3, whereas q=1 mod 6.
- p<q<2p. The lower inequality follows from alpha>1. For the upper inequality, 3^3=27<32=2^5 gives alpha<5/3. Since p>=18,

      q<alpha*p+6<(5/3)*p+6<=2p.

- r>p/2, because 3^2=9>8=2^3 gives alpha>3/2 and q>alpha*p.
- p=0 mod 6, r=1 mod 6, and n=5 mod 6. The proved positivity of n therefore implies n>=5.

The a>=2 boundary is essential for this construction. At a=0, p=2 gives q_0=4,q=7,n=-3; at a=1, p=6 gives q_0=10,q=13,n=-1. Neither is in p<q<2p. The first valid member is p=18,q=31,r=13,n=5.

All members have positive clock D=2^q-3^p. They are near-critical in the explicitly bounded additive sense above, not necessarily exact critical clocks q=q_0.

## 2. The endpoint ++ chain is an actual two-swap path

Use the inherited mechanical labels a_j=(-jr mod p), and consider the unique one whose residue is A=p-1. Its following one has residue

    B=A-r=p-r-1=n-1.

Because A>=r, the gap from the A-labelled one to that next one is 1: the first one cannot initially move right. Because B=n-1<r, the following one has a zero immediately on its right. Hence the source contains the labelled local pattern 110, and the legal ordered path is

    110 -> 101 -> 011.

The B-labelled leading one moves first; the A-labelled following one moves second. This argument keeps the intermediate word and includes original-cut crossings. The source is primitive because gcd(p,q)=1. No same-one displacement or illegal formal sign is used.

Its exact integral plus-plus norm representative is

    Q_p(Z)=3-Z^(p-1)+Z^n-Z^(n-1).

The general one/two-edit reduction therefore supplies the ordinary content divisibility gate g | M, where M is the positive 2,3-free part of the nonzero resultant Res(Z^p-2,Q_p). That gate does not identify g with M.

## 3. The resultant is a 2-adic unit

Here n>=5, and p-1>0, so Q_p has constant coefficient 3 and all its other displayed exponents are positive. Over F_2,

    Z^p-2=Z^p,
    Res(Z^p,Q_p)=Q_p(0)^p=1.

Consequently the integer resultant is odd. There is no suppressed dyadic denominator in this plus-plus representative.

## 4. The resultant is a 3-adic unit

Over F_3, Frobenius gives the exact identity

    Z^p-2=Z^(2*3^a)+1=(Z^2+1)^(3^a).

The quadratic Z^2+1 is irreducible over F_3. Let zeta be either of its two roots in F_9. Then zeta^2=-1, so zeta has order 4. Moreover

    (zeta-1)^2=zeta,
    (zeta-1)^4=-1,
    (zeta-1)^8=1.

Thus zeta-1 has order exactly 8. Every zeta^r has order dividing 4, so zeta^r cannot equal zeta-1, for any integer r.

Reduction of the edit polynomial gives

    Q_p(Z)=Z^(n-1)(-Z^r+Z-1)   over F_3,

because p=n+r. Its nonzero-root factor cannot vanish at either zeta. Therefore Q_p and Z^p-2 have no common root over the algebraic closure of F_3, hence are coprime. The repeated multiplicity 3^a in the first polynomial does not change that criterion. Their resultant is nonzero modulo 3.

Combining Sections 3 and 4 proves the exact equality

    M=|Res(Z^p-2,Q_p)|

for every a>=2. This resolves the stripping issue on this selected infinite family only; it does not resolve it for all critical or near-critical shapes.

## 5. The norm limit applies without exact criticality

The earlier torus proof requires only p->infinity and

    n/p -> beta=2-log_2 3,

not exact equality q=ceil(p log_2 3). The family satisfies that limit because 0<q-alpha*p<6 and p->infinity. The number beta is irrational.

For completeness, modulo Z^p-2,

    Z Q_p(Z)=H_p(Z)=3Z-2+(Z-1)Z^n.

The absolute norms of H_p and Q_p differ by exactly 2. Write rho_p=2^(1/p). The p conjugates are rho_p exp(2pi i j/p). Their logarithmic integrands are samples on the two-torus of

    f_p(theta,phi)=log|3rho_p exp(i theta)-2
                      +rho_p^n(rho_p exp(i theta)-1)exp(i phi)|.

They converge uniformly to

    f(theta,phi)=log|3 exp(i theta)-2
                    +(4/3)(exp(i theta)-1)exp(i phi)|,

because rho_p->1, rho_p^n->2^beta=4/3, and the perturbation is uniformly dominated. Indeed, for rho>=1,

    max_theta |rho exp(i theta)-1|/|3rho exp(i theta)-2|
      =(rho+1)/(3rho+2)<=2/5,

while rho_p^n<2. Thus the perturbation ratio is <4/5 and no logarithmic singularity occurs.

The sample points (2pi j/p,2pi n j/p) equidistribute on the two-torus. For a fixed nonzero integer character (h,l), its sample mean is zero unless p divides h+ln. This divisibility eventually fails: if l!=0, h/p+l(n/p) tends to the noninteger l beta; if l=0, a fixed nonzero h is eventually smaller than p. Approximation by trigonometric polynomials proves the asserted equidistribution.

Averaging f first over phi gives log|3 exp(i theta)-2| because the linear term dominates. Its theta-average is log 3. Consequently

    (1/p) log |Norm H_p| -> log 3,
    M^(1/p)=|Norm Q_p|^(1/p) -> 3.

The finite 2 factor between the norms disappears in the limit. The unit results in Sections 3 and 4 make this a stripped-norm limit, not merely a raw-norm limit.

## 6. Exact obstruction and limits of its scope

For every fixed C>0 and 0<c<3, the constructed family eventually satisfies

    M>C*c^p.

For example, choose d with c<d<3. The root limit gives M>d^p eventually, and (d/c)^p eventually exceeds C. Thus no eventual bound M<=C*c^p with c<3 can cover all actual positive near-critical two-edit shapes, even when q-p log_2 3 is restricted to (0,6), the source is primitive, and only these overlapping ++ paths are used.

This does NOT establish any of the following:

- failure of such a bound restricted to exact critical clocks q=ceil(p log_2 3);
- a corresponding lower bound on actual content g;
- an integer cycle or failure of two-swap cycle exclusion;
- that exact critical clocks exhaust the residual two-edit problem.

The missing clock/common-root information remains relevant to g even though no 2 or 3 factor is lost from this family's norm.

## 7. Independent finite spot checks

These are exact integer checks of the first three valid members, not a general census. The q_0 values were computed as bit_length(3^p), the ordered swaps were checked directly on the mechanically generated position set, and norms were computed by an independent SymPy integer resultant call.

| a | p | q_0 | q | r | n | ordered zero-based edges | absolute norm |
|---|---:|---:|---:|---:|---:|---|---:|
| 2 | 18 | 29 | 31 | 13 | 5 | 13, 12 | 261006559 |
| 3 | 54 | 86 | 91 | 37 | 17 | 33, 32 | 57871986057582534182735809 |
| 4 | 162 | 257 | 259 | 97 | 65 | 252, 251 | 79822855858369572242004884889147943552874279970694998545424737296828349812177 |

All three displayed norms are 1 modulo 6. That residue is a finite observation, not an additional claim for every a. Their positive p-th roots are approximately 2.93489, 2.99973, and 2.98335; non-monotonic finite values are consistent with the proved limit.

Verdict remains **PASS**, with no required mathematical correction. The a>=2, near-critical-only, norm-versus-content, and fixed-prefactor qualifications are essential.
