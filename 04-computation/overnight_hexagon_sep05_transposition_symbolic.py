"""Exact polynomial controls; not a finite proxy for the all-order proof."""
import sympy as s

checks = 0


def check(value, label):
    global checks
    if not value:
        raise RuntimeError(label)
    checks += 1


def identity(a, b, label):
    check(s.expand(a - b) == 0, label)


n, r, t, k, q = s.symbols("n r t k q")
v = n - 2 - r
position = 4*r*v*(n-4) + 2*r*v*((r-1)*(r-2)+(v-1)*(v-2))
closed = 2*r*v*((n-2)*(n-3)-2*r*v)
identity(position, closed, "three cycle-position cases")
m = n-2
EX = 2*r*v/(m-1)
EX2 = 4*r*v*(r*v-1)/((m-1)*(m-2))
identity(s.cancel((m-1)*(m-2)*((m+1)*EX-EX2)), closed,
         "independent crossing-gap moments")
F = q*((n-2)*(n-3)-2*q)
C = (n-2)*(n-3)*(n-4)
identity(F.subs(q, 2*(n-4))-C, (n-4)*(n*n-13*n+38), "lower endpoint")
identity(F.subs(q, (n-2)**2/4)-C,
         (n-2)*(n-4)*(n*n-12*n+28)/8, "upper real endpoint")
identity(s.diff(F,q,2), -4, "strict concavity")
for poly in (n*n-13*n+38, n*n-12*n+28):
    shifted=s.Poly(poly.subs(n,t+9),t)
    check(all(c>0 for c in shifted.all_coeffs()), "all-height endpoint positivity")
identity((n*(n-1)-2*(k-1)*(n-k+2)).subs(n,k+t),
         (k-1)*(k-4)+t+t*t, "antibalanced odd/even pairing")
identity((n-1)*(n-k)/(n-2)/(n-k), (n-1)/(n-2), "deletion multiplicity ratio")
check(s.cancel((n-1)/(n-2)-1-1/(n-2))==0, "strict deletion excess")
# Quantitative strengthening: minimum over the exact integer r range.
for order in range(9,201):
    vals=[s.Integer(split)*(order-2-split)*
          ((order-2)*(order-3)-2*split*(order-2-split))
          for split in range(2,order-3)]
    lo=2*(order-4)
    hi=(order-2)**2//4
    formula=min(lo*((order-2)*(order-3)-2*lo),
                hi*((order-2)*(order-3)-2*hi))
    check(min(vals)==formula, "concavity on all integer splits")
    check(formula>(order-2)*(order-3)*(order-4), "strict transposition threshold")
print(f"symbolic_identities_and_exact_controls={checks}")
print("all_height_proof=exact_concavity_and_positive_shifted_coefficients")
print("finite_auxiliary_integer_endpoint_universe=9..200")
print("RESULT=PASS")
