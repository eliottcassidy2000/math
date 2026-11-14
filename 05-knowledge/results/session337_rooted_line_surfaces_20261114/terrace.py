from itertools import combinations
def pathcount(gaps,z=2):
    k=len(gaps)+1
    G=[0]
    for r in gaps:G.append(G[-1]+r)
    P=[1]
    for j in range(k):
        total=P[-1]
        for i in range(j):
            R=G[j]-G[i]
            c=2**(R-1) if R else int((j-i)%2==0)
            total+=z*c*P[i]
        P.append(total)
    return P[-1]
def comp(n):
    for mask in range(1<<(n-2)):
        end=(0,)+tuple(i+1 for i in range(n-2) if (mask>>i)&1)+(n-1,)
        yield tuple(b-a-1 for a,b in zip(end,end[1:]))
def canonical(n):
    for mask in range(1<<(n-4)):
        end=(0,)+tuple(i+2 for i in range(n-4) if (mask>>i)&1)+(n-1,)
        yield tuple(b-a-1 for a,b in zip(end,end[1:]))
if __name__=='__main__':
    for n in range(4,15):
        d=sorted({pathcount(g) for g in canonical(n)})
        print(n,len(d),d[:25],d[-1],flush=True)

def four_state(gaps):
    # E=closed, A=open contains positive gap, O/V=pure-zero open odd/even.
    E,A,O,V=1,0,0,0
    for r in gaps:
        C=E+A+2*V
        if r:E,A,O,V=C,(E+A+O+V)*(1<<r),0,0
        else:E,A,O,V=C,A,E+V,O
    return E+A+2*V

def pair_inverse(N):
    "Solve N=(1+2^a)(1+2^b), a,b>=1, unordered."
    if N==9:return (1,1)
    if N<15:return None
    s=N.bit_length()-1
    D=N-1-(1<<s)
    if D<2:return None
    ex=[i for i in range(1,D.bit_length()) if D>>i&1]
    if len(ex)==2 and sum(ex)==s and D==sum(1<<i for i in ex):return tuple(ex)
    if s%2==0 and s>=2 and D==1<<(s//2+1):return (s//2,s//2)
    return None
