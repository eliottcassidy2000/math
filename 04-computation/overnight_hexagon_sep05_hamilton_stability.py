"""Exact controls for analytic Hamilton spectral isolation, no large census."""
from itertools import combinations, permutations
from math import factorial
import sympy as s

gates = 0


def check(value, label):
    global gates
    if not value:
        raise RuntimeError(label)
    gates += 1


def identity(left, right, label):
    check(s.cancel(left-right) == 0, label)


n, t, q, k = s.symbols("n t q k")
F = lambda q: q*((n-2)*(n-3)-2*q)
low = F(2*(n-4))
even_difference = (n-10)*(n-6)**2*(n-4)/8
odd_difference = (n-7)*(n-5)*(n*n-14*n+41)/8
identity(F((n-2)**2/4)-low, even_difference, "even integer endpoint difference")
identity(F(((n-2)**2-1)/4)-low, odd_difference, "odd integer endpoint difference")
check(all(c>=0 for c in s.Poly(even_difference.subs(n,10+t),t).all_coeffs()),
      "even endpoint positive for all n>=10")
check(all(c>0 for c in s.Poly(odd_difference.subs(n,11+t),t).all_coeffs()),
      "odd endpoint positive for all n>=11")
lam = 2*(n*n-9*n+22)/((n-2)*(n-3))
upper = 2*(n-4)/(n-2)
identity(low/((n-2)*(n-3)*(n-4)), lam, "lower isolation ratio")
identity(2-lam, 8*(n-4)/((n-2)*(n-3)), "asymptotic lower deficit")
identity(upper-lam, 4*(n-5)/((n-2)*(n-3)), "explicit bracket width")
identity(s.Rational(min(10*(42-20),12*(42-24))*factorial(4),factorial(7)),
         s.Rational(36,35), "n9 exact integer endpoint")
identity((2-2/(n-2))-(2-4*(k-3)/((n-2)*(n-3))),
         2*(2*k-n-3)/((n-2)*(n-3)), "adjacent/disjoint scale transition")
check(s.limit(lam,n,s.oo)==2 and s.limit(upper,n,s.oo)==2, "second minimum squeeze")

cycle_total=0
for order in range(5,10):
    for length in range(3,order+1):
        total=both_adj=both_disj=odd_adj=odd_disj=0
        for chosen in combinations(range(order),length):
            root=chosen[0]
            for rest in permutations(chosen[1:]):
                if rest[0]>rest[-1]:
                    continue
                word=(root,)+rest
                edges={tuple(sorted((word[j],word[(j+1)%length])))
                       for j in range(length)}
                one=(0,1) in edges
                adj=(0,2) in edges
                disj=(2,3) in edges
                total+=1
                both_adj+=one and adj
                both_disj+=one and disj
                odd_adj+=one != adj
                odd_disj+=one != disj
        e=factorial(order-2)//factorial(order-length)
        aj=factorial(order-3)//factorial(order-length)
        dj=(2*(length-3)*factorial(order-4)//factorial(order-length)
            if length>=4 else 0)
        check(total==factorial(order)//(2*length*factorial(order-length)), "full cycle universe")
        check(both_adj==aj and both_disj==dj, "two-edge intersection counts")
        check(odd_adj==2*e-2*aj and odd_disj==2*e-2*dj, "literal pair parity")
        cycle_total+=total
    print(f"n={order} all_layers=3..{order} adjacent_disjoint_counts_and_scale_transition=PASS")
print(f"exact_failure_gates={gates} directly_enumerated_cycles={cycle_total}")
print("second_minimum_lower_n9=36/35")
print("second_minimum_lower_n_ge10=2*(n*n-9*n+22)/((n-2)*(n-3))")
print("second_minimum_upper_n_ge9=2*(n-4)/(n-2)")
print("RESULT=PASS")
