"""The two-slot invariant Q_w and the threshold function theta_w, for many weights."""
from fractions import Fraction
from math import factorial, comb

def Q(w,a,b):
    return w(2*a)*w(b)**2 - 2*w(a)*w(b)*w(a+b) + w(a)**2*w(2*b)

SEQS = {
 "j!":                 lambda j: factorial(j),
 "geometric 3^j":      lambda j: 3**j,
 "central binom C(2j,j)": lambda j: comb(2*j,j),
 "Catalan":            lambda j: comb(2*j,j)//(j+1),
 "ballot C(2j,j-1)":   lambda j: comb(2*j,j-1) if j>=1 else 0,
 "(j!)^2":             lambda j: factorial(j)**2,
 "double fact (2j-1)!!": lambda j: factorial(2*j)//(2**j*factorial(j)),
 "Fibonacci":          lambda j: [0,1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765][j],
 "2^(j^2)":            lambda j: 2**(j*j),
 "1/j! (as j!^-1 scaled)": lambda j: factorial(12)//factorial(j) if j<=12 else 0,
}
print("Q_w(a,b) = w_2a w_b^2 - 2 w_a w_b w_(a+b) + w_a^2 w_2b")
print("Q != 0 for all a!=b  =>  two-slot threshold is 2 (two moments force f=0)")
print("Q == 0 identically   <=> w geometric <=> L_w is an evaluation, FC FAILS")
print()
print(f"{'weight w_j':26s} {'log-convex?':12s} {'Q>0?':8s} {'Q=0 anywhere?':14s} verdict")
for name,w in SEQS.items():
    lc=True; pos=True; zero=False
    for a in range(0,7):
        for b in range(a+1,8):
            try:
                if w(a)==0 or w(b)==0: continue
                q=Q(w,a,b)
                if q==0: zero=True
                if q<0: pos=False
                if w(a+b)**2 >= w(2*a)*w(2*b): lc=False
            except Exception: pass
    verdict = "FC FAILS (evaluation)" if zero and not pos else ("threshold 2" if not zero else "degenerate")
    print(f"{name:26s} {str(lc):12s} {str(pos):8s} {str(zero):14s} {verdict}")

# --- sharpness: Fibonacci-recurrence weights and f = s(1-s) ---
print()
print("SHARPNESS (THM-3022 sec 4): w_{n+1}=w_n+w_{n-1} => Delta w_n = w_{n-1},")
print("so L_w(f^m) = (-1)^m w_0 for f = s(1-s), every m.")
from math import comb as _c
for name, seq in (("Fibonacci (w_0=0)", [0,1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765]),
                  ("Lucas (w_0=2)",     [2,1,3,4,7,11,18,29,47,76,123,199,322,521,843,1364,2207,3571,5778,9349,15127])):
    vals=[sum(_c(m,j)*(-1)**j*seq[m+j] for j in range(m+1)) for m in range(1,9)]
    print(f"   {name:20s} L(f^m), m=1..8: {vals}")
print("   -> Fibonacci is an exact two-slot counterexample; Lucas never vanishes.")
print("   The factorial escapes because Delta(n!) = n*n!, never a shift.")
