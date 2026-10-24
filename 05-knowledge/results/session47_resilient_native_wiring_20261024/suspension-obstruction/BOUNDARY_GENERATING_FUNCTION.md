# A two-port generating function for suspended target cycles

Status: PROVED formal identities; FINITE-EXACT controls. These count closed
walks and, at lengths3 and4, cyclic subgraphs. They do NOT count Hamiltonian
paths or supply an electrical harmonic transfer.

## The retained information is two boundary ports

Use the incoming suspension notation and the accompanying THEOREM.md:
K=K_P(T) is the inner target, f and l are its path's first and last named
owners, and H=K_(a,P,b)(S(T)). Treat K,H as adjacency matrices. Put

    u=1-e_f, v=e_l, w=1-e_l.

With old owners first, then a,b, the full actual matrix is

    H = [ K       u   v ]
        [ e_f^t   0   0 ]
        [ w^t     1   0 ].

Thus the needed boundary data is explicitly the pair of named endpoint
vectors, not a choice of matrix eigenvalues or an unmarked isomorphism class.

## Exact determinant and walk series

Over formal power series, R(z)=(I-zK)^(-1) exists because its constant term
is I. Define the four scalar boundary series

    a=e_f^t R u, b=e_f^t R v, c=w^t R u, d=w^t R v.

Then the exact identity is

    det(I-zH)/det(I-zK)
       = 1-z^2(a+d)-z^3 b+z^4(ad-bc).                 (1)

Proof. Apply the block determinant identity to I-zH, with old block I-zK.
Its two-by-two Schur block has entries

    [1-z^2 a,       -z^2 b]
    [-z-z^2 c,      1-z^2 d].

Expanding gives(1). The inversion is formal, so no assumption about numerical
convergence, invertibility at a particular temperature, or self-adjointness
is involved. The final product with det(I-zK) is a polynomial.

Writing the right side as Q(z), formal differentiation gives

    -z Q'(z)/Q(z) = sum_(r>=1) (tr(H^r)-tr(K^r)) z^r.  (2)

This follows from log det(I-zA)=-sum_(r>=1) tr(A^r) z^r/r,
equivalently from formal differentiation of a determinant. Since K is an
induced subgraph, the trace difference counts rooted closed directed walks
that touch at least one outer port. At higher orders these walks can repeat
vertices; interpreting them as simple cycles would be incorrect.

## Triangle conservation and a fourth-order refinement

There are no loops or directed two-cycles. At length3 every closed walk is a
cyclic triangle counted three times; at length4 it is a simple cyclic
four-cycle counted four times. Therefore

    boundary cyclic triangles = n+1-d_T(f)+d_T(l),

as proved independently by the degree-square identity, and

    boundary cyclic four-cycles
       = e_f^t K^2(1-e_f)
         +(1-e_l)^t K^2 e_l
         +K_(f,l).                                    (3)

For(3), the first two terms count cycles using a alone and b alone; the final
term counts the unique possible pattern using both, a->f->l->b->a. Equally,
extract the fourth coefficient in(1)-(2), noting a(0)=d(0)=b(0)=0.

This makes the distinction between third- and fourth-order data concrete.
The total enlarged triangle count depends only on the source score sequence
(and n), while the fourth-order boundary term retains two-step incidence at
the old endpoints. Equal global counts do not erase that endpoint memory.

## The connection contract

Source: a full marked source/path record.
Target: an inner adjacency matrix plus two named boundary response ports.
Preserved: the entire enlarged adjacency matrix and all rooted closed-walk
counts through(1)-(2).
Lost if the ports are removed: which rows/columns define the response series.
Required extra data: e_f,e_l and the actual oriented inner matrix.
Decisive test: the accompanying same-old-class/different-enlarged-class
counterexample and direct exact matrix-series checks.

The block calculation is standard linear algebra applied to this native
suspension. Its use here is an explicit boundary-accounting tool, not a claim
that tournaments are physical scattering systems or a solution of an open
spectral problem.

## Reproduction

`check_generating.py` checks(3) on all1,099 forward-path tournaments of
orders2..6. On all75 such tournaments of orders2..5 it independently builds
both determinant polynomials from exact Newton trace recurrences and checks
(1)-(2) through degree n+6, including all coefficients beyond the determinant
degrees. No floating point is used. Normal and optimized runs agree.
