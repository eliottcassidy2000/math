# When ternary Fourier fails: a lossless Newton guard compiler

Status: PROVED for supplied finite guard tables; exact tests are supplementary. This does not change the original electrical guard's modulus or authenticate a new native source network.

## Inheritance and the meaningful other side

The incoming `native_guard_moves_farey_20261023/fourier/THEOREM.md` proves a Fourier compiler over `(Z/8)[w]/(w^2+w+1)`. Its mod-nine hostile, the nonzero constant table `(3,3,3)` with zero unnormalized Fourier transform, is correct. The lost information is caused by the nonunit averaging factor three. It is not a reason that mod-nine guard tables cannot be encoded exactly.

The older `01-canon/theorems/THM-2206-quadratic-deck-augmentation-and-hamiltonian-characteristic-incompatibility.md` gives the binary law `(T-I)^j=(-2)^(j-1)(T-I)` for an involution. The ternary replacement below has a genuinely different finite-difference filtration. The local-versus-split distinction also occurs in `THM-2845-local-residue-versus-split-trace-unit-observability.md`; no Jacobian or unit-detection theorem is inferred from the present guard encoding.

## 1. Integral interpolation without averaging

Let R be ANY commutative ring and let f:{0,1,2}^r -> R be a supplied complete table. For alpha in {0,1,2}^r define

    c_alpha = sum_(beta<=alpha) (-1)^(|alpha|-|beta|)
                  product_i binom(alpha_i,beta_i) f(beta).

Then exactly

    f(a) = sum_alpha c_alpha product_i binom(a_i,alpha_i).

Binomial coefficients here are ordinary integers evaluated on the representatives 0,1,2 BEFORE mapping into R. No division by two or three in R is used. In one coordinate the evaluation matrix is

    [[1,0,0], [1,1,0], [1,2,1]],

which has determinant one over the integers. Its inverse gives c_0=f(0), c_1=f(1)-f(0), c_2=f(2)-2f(1)+f(0). Tensoring proves existence and uniqueness for r coordinates. Thus all powers of two, all powers of three and mixed moduli are equally lossless for this representation.

An axis i is essential exactly when some nonzero c_alpha has alpha_i>0. Independence of that axis forces those coefficients to vanish by uniqueness; their vanishing visibly removes the axis from the reconstruction.

This is an interpolation on a chosen three-element representative set, not an ordinary polynomial function obtained by substituting the relation 3=0 in the target ring. The wrap of the cyclic input must be retained separately.

## 2. The cyclic wrap and the sharp ternary filtration

Let T act by (Tf)(a)=f(a+1 modulo3), and write D=T-I. On the three Newton coefficients,

    T(c0,c1,c2) = (c0+c1, c1+c2, -3c1-2c2).

This is an integral order-three matrix. Since T^3=I,

    D^3=-3TD,
    D^(2k+1)=(-3)^k T^k D,
    D^(2k)=(-3)^(k-1) T^(k-1) D^2       (k>=1).

The latter two identities follow by induction because T and D commute. On ALL three-point tables over Z/(3^k), the operator D therefore has nilpotency index EXACTLY 2k+1. The upper bound follows from the displayed identity. For the lower bound, D^2 applied to a delta table has an entry equal to one or minus two, hence not divisible by three; multiplying by 3^(k-1) and permuting entries does not kill it modulo3^k.

In particular D^5=0 but D^4!=0 modulo9. The same ternary permutation splits into eigensummands after adjoining roots at modulus8, but has a nontrivial finite-difference filtration at modulus9. These are alternative useful representations of the same input action, not competing claims about which modulus is correct.

When three is a unit, constants split off from the zero-sum submodule. D is invertible on that zero-sum submodule: its restriction has determinant three in an integral basis. Thus no nonzero power of D kills that whole submodule. This proves the contrasting behavior without numerical roots.

The sharp index stated here is for ONE cyclic coordinate. It is not a claimed nilpotency index for the total multivariable augmentation ideal.

## 3. Exact linear observers still exist without Fourier inversion

For ANY output set, define the literal translation stabilizer

    K = {v in F3^r : f(a+v)=f(a) for every a}.

It is a subgroup and therefore an F3-linear subspace: translations compose and invert, regardless of the output ring. The table factors through F3^r/K. If B is any linear map with kernel K, then f(a)=g(Ba), with exactly 3^dim(K) lifts for each B-value. This has minimum rank among linear observations through which f factors, since every sufficient C must satisfy ker(C) subset K.

A lossless compiler can therefore:

1. retain the complete Newton table (or the original values);
2. compare actual translated tables to find K;
3. row-reduce K over F3 and retain each actual fibre coordinate.

This may take exponential time. It does not inherit the incoming Fourier method's character-support shortcut. In particular the span of indices of nonzero Newton coefficients is NOT a substitute for the Fourier frequency span. For f(a,b)=(a+b modulo3), embedded as values0,1,2 in Z/9, the literal kernel is span(1,2), so the minimal observer has rank one; the nonzero Newton indices include (1,0) and (0,1), whose span has rank two. The meaning of a coefficient index changed.

## 4. Native guard scope and shared frames

For an authenticated native frame parameterization, apply this to the ENTIRE guard record, with validity separately retained and fixed owner/arc pins substituted first. Equal source/target owner labels, converse transport and exact clocks remain required external data. Encoding mod-nine tables does not turn an original mod-eight electrical character into a mod-nine character. The two moduli can be kept as separate record components if both have independently authenticated meaning.

A ternary decision tree enumerates these tables by prefixes. The Newton representation retains local finite differences; the translation quotient retains global periods. Neither representation guarantees that equal-valued branches are connected by allowed native moves. That last obligation belongs to the actual coupled-move and holonomy certificates.

## Reproduction

Run `python3 -B check.py` and `python3 -O -B check.py` from this directory. The checker exhausts all729 one-coordinate mod-nine functions, all512 two-coordinate Boolean tables, the sharp operator filtration for k=1,...,6, and the contrasting mod-eight zero-sum action. It separately checks the false Newton-support-rank shortcut. These finite checks complement the all-ring proofs above.
