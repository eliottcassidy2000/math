"""MOSEK: rigorous UPPER bound on max_E L_y (THM-534 sector crux). Boolean-moment relaxation over J(A)=meas{all
runners avoid sectors A}, A subset {1..6}: moment matrix M[a,b]=J(a|b)>>0, single-runner bound J(A)<=(7-|A|)/7,
monotone. max L_y=sum_r y_r sum_{|A|=r} J(A). If <= cap_k => proves L_y<=cap for ALL E (closes the crux)."""
import numpy as np, cvxpy as cp
def popcount(a): return bin(a).count("1")
def solve_Ly_bound(y,cap,Lconsec,extra_bounds=None):
    n=64; J=cp.Variable(n)
    M=cp.bmat([[J[a|b] for b in range(n)] for a in range(n)])
    cons=[J[0]==1, M>>0, J>=0]
    for a in range(n): cons.append(J[a] <= (7-popcount(a))/7.0)         # single-runner achievability
    for a in range(n):                                                   # monotone
        for j in range(6):
            if not (a>>j)&1: cons.append(J[a] >= J[a|(1<<j)])
    if extra_bounds:
        for a,ub in extra_bounds.items(): cons.append(J[a] <= ub)
    Ly = cp.sum([y[popcount(a)]*J[a] for a in range(n)])
    prob=cp.Problem(cp.Maximize(Ly),cons); prob.solve(solver=cp.MOSEK)
    return prob.value, prob.status
# k=8: L_y = 1 - S1 + S2 - 0.9 S3 + 0.6 S4
y8=[1,-1,1,-0.9,0.6,0,0]; cap8=0.3815; Lc8=0.35823
val,st=solve_Ly_bound(y8,cap8,Lc8)
print(f"k=8: MOSEK max_E L_y (boolean-moment + single-runner bound) = {val:.5f} ({st})")
print(f"     L_y(consec)={Lc8:.5f}, cap_8={cap8:.5f}")
print(f"     => {'<= cap: CLOSES the k=8 sector crux!' if val<=cap8+1e-4 else 'LOOSE (>cap): need more achievability constraints; gap='+format(val-cap8,'+.4f')}")
