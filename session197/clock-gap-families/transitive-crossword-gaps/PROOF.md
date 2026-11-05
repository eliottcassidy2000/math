# Exact cancellation criterion for an arbitrary-gap transitive cross word

Session 197. All-order four-parameter formula and an exact pure-ternary criterion within one explicit word skeleton. This is separate from the previously frozen four-seed constructive packet. It does not classify arbitrary four-letter cross words or arbitrary unit clocks.

## 1. Physical family and closed form

Take a,u,b in the transitive order a->u->b, a->b. A transitive chain I has word

    1 0^i 1 0^j 2 0^k 2,

now with arbitrary integers i,j,k>=0. Type 1 loses to a and beats u,b; type 2 loses to b and beats a,u; type 0 beats a,u,b. Earlier I vertices beat later ones. Add a transitive module R of order m>=0, with I->R->{a,u,b}. Empty modules are omitted. Write H(i,k,j,m) for this actual tournament's Hamilton-path count.

Define

    U_t = 1 if t=0, and U_t=2^t-1 if t>=1,
    f_i = 2^(i+1)+1,
    A = 4 f_i (3^(k+1)+2^(k+1)) 2^j,
    D = -2 [ f_i (3^k-2^(k+1)+1) 2^j
              +(2^(k+1)+1)(U_i-1) ],
    C = U_i U_k - 2[j=0].

**Theorem 1.** For all i,k,j,m>=0,

    H(i,k,j,m)=A 3^m+D 2^m+C.                        (1)

The zero-module cases are included through explicit boundary checks below. They are not inferred by substituting zero into a positive-multiplicity Stirling identity.

**Theorem 2.** The R-clone family is a pure nonconstant ternary affine clock, even eventually, if and only if

    i,k in {0,1}.                                   (2)

In that case D=0, A=4(3+2i)(5+8k)2^j, and C=-1 for j=0 or +1 for j>=1. Thus every pure ternary clock in this arbitrary-gap skeleton automatically has unit offset. The four primitive amplitudes are 60,100,156,260.

For the proof of (2), put g_k=3^k-2^(k+1)+1. Then g_0=g_1=0 and g_2=2. The recurrence

    g_(k+1)=3g_k+2^(k+1)-2

shows g_k>0 for every k>=2. Also U_i-1 is zero exactly for i=0,1 and positive thereafter. All other factors in -D are positive, so D<=0 and D=0 exactly in (2). The coefficient A is always positive. Finally

    H(i,k,j,m+1)-3H(i,k,j,m)=-D 2^m-2C

is eventually constant exactly when D=0. This proves the eventual as well as exact version without assuming bounded parameters.

## 2. Uniform matching bounds for all four modules

For a transitive module replaced by n>=1 clones, collapse each maximal clone run in a Hamilton word. If the module has t runs, there are t-1 disjoint exterior O-to-I fragments, where O and I are the out- and in-neighbors of a representative. Selecting one O-to-I edge per fragment gives a matching. Thus matching rank at most nu bounds the number of runs by nu+1.

A fixed t-run word lifts in t! S(n,t) ways, assigning the distinct clone labels to ordered nonempty runs. Since this is a linear combination of q^n for q=1,...,t, rank at most nu puts the positive-n sequence in the span of 1,2^n,...,(nu+1)^n.

The required bounds are uniform in every other module size, including their zero boundaries:

- The i-gap representative has only the first type-1 anchor before it in I. Its crossing graph is a nonempty star centered at a, hence rank one.
- The j-gap representative has the two type-1 anchors and the i-gap before it. Only a points from its out-neighborhood into that earlier chain segment, hence rank one.
- The k-gap representative has both type-1 anchors and the first type-2 anchor before it. All crossing edges meet {a,b}; a and b point to distinct earlier anchors, hence rank two.
- An R representative has out-neighbors a,u,b and the entire chain as in-neighbors. All crossing edges meet {a,b}, with disjoint edges into the type-1 and type-2 anchors, hence rank two.

In these arguments every later chain vertex loses to the earlier chain segment. R vertices also lose to all I. The isolated outside vertex u points to no chain vertex. Thus adding any number of vertices to any other module introduces no unaccounted crossing edge. When the module under consideration is positive, all but one of its vertices are first collapsed; the argument concerns the source of that module substitution.

## 3. A complete finite tensor certificate, including zero boundaries

Let delta_n be 1 at n=0 and 0 otherwise. A function on nonnegative integers that is in the span of 1,2^n on positive n, with an arbitrary value at zero, lies in

    span{delta_n,1,2^n}.

Likewise a positive-n sequence spanned by 1,2^n,3^n with arbitrary deletion value lies in span{delta_n,1,2^n,3^n}. The first space is determined uniquely by values at n=0,1,2; the second by n=0,1,2,3. Positive Vandermonde interpolation determines the exponential coefficients, and the zero value determines the remaining delta coefficient.

By the uniform matching bounds, the actual Hamilton count belongs separately to these spaces in i,j and in k,m. Therefore it is uniquely determined on the entire nonnegative parameter lattice by the rectangular grid

    i in {0,1,2},  k in {0,1,2,3},
    j in {0,1,2},  m in {0,1,2,3}.                    (3)

This is 144 actual tournaments. Uniqueness follows by extending successively in each coordinate, not by assuming a total-degree bound. All zero/nonzero regimes are covered by the delta basis.

The right side of (1) belongs to the same four separate spaces: U_t=2^t-1+delta_t; all other displayed factors use only the listed exponential bases and delta_j. The standalone verifier counts every graph in (3) by two direct methods and checks (1). Consequently tensor uniqueness proves (1) for all parameters. No bound on i,j,k,m and no finite-word stabilization hypothesis remains in the theorem.

## 4. Certificate and additional consequences

`verify.cpp` constructs every arc explicitly and performs exact endpoint subset DP and a separate literal Hamilton-word DFS on every grid graph. The largest order is 17. The receipt lists all 144 parameter tuples, orders, A,D,C and independently counted H values. The implementation also checks each relevant one-representative crossing rank on the grid and the exact finite D=0 controls.

All path counts fit unsigned 64-bit integers: even 17! is below 2^64. Formula arithmetic is explicitly signed 64-bit, with safe magnitudes on this finite grid. All checks throw exceptions rather than depend on assertions. Optimized and unoptimized C++17 runs give identical receipts. The verifier is standalone and imports neither the earlier four-seed packet nor the exploratory cross-word transfer.

For an R-singleton source, equation (1) and its separately proved m=0 value give the entire root-run polynomial

    R(z)=(z+1)[A(z+1)^2+D(z+1)+C].                   (4)

Thus the alternating root boundary vanishes throughout this arbitrary-gap family, including cases with D<0 that are not pure ternary clocks. Zero boundary alone therefore does not impose pure exponential cancellation.

Leading zero letters before the first type-1 anchor are universal-source padding and do not alter H. Trailing zero letters after the last type-2 anchor are transitive preclones of R and replace m by m+d. They do not change the criterion D=0. This extends the exact pure-clock characterization to those two arbitrary outer gaps as well.

For i,k in {0,1}, the theorem recovers the independently packaged four-seed unit families, with physical prime and 7p examples in `session197-prime-atoms/two-module-unit-clocks/PROOF.md`. The primitive-amplitude-60 binary/ternary bridge is the earlier construction in `session197-hybrid-constructions/binary-ternary-bridge/PROOF.md`, SHA256 `65c57c4f232dc77713b98048765ac901e4039ee761d518066817c84395ef3fc7`. Neither earlier packet is needed to run or prove this certificate.

The restriction to the nonzero skeleton1122 is essential to the stated scope. No theorem here excludes words containing type3, more nonzero letters, different orders of the nonzero letters, or other orientations of a,u,b. A finite scan suggested a larger regular-language conjecture, but that conjecture is not used or promoted by this packet. No full prime/7p coverage or unbounded primality claim follows from (1).
