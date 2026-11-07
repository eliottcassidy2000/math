"""Exact finite totient prefix and Camion k=50 hostile, stdlib only."""
from math import comb,isqrt
B=2*76**2
# phi(n)>=sqrt(n/2): each odd prime-power contributes phi(p^a)/sqrt(p^a)>=1;
# the 2-power contributes at least 1/sqrt(2). Hence this bound is complete.
phi=list(range(B+1))
for p in range(2,B+1):
    if phi[p]==p:
        for n in range(p,B+1,p): phi[n]-=phi[n]//p
image=set(phi[1:])
missing=[n for n in range(2,77,2) if n not in image]
assert missing==[14,26,34,38,50,62,68,74,76]
p=251
assert all(p%d for d in range(2,isqrt(p)+1))
n=2**50-1
assert n%p==0
r=(n-1)%(p-1)+1
T=[pow(2,k*(k-1)//2,p) for k in range(r+1)]
s=[0]*(r+1)
for k in range(1,r+1):
    s[k]=(T[k]-sum(comb(k,j)*s[j]*T[k-j] for j in range(1,k)))%p
D=(s[r]-pow(2,(r-1)//2,p))%p
assert (r,D)==(123,229)
print({'complete_preimage_bound':B,'even_nontotients_through76':missing,'k50_obstruction':{'prime':p,'reduced_index':r,'D_mod_p':D}})
