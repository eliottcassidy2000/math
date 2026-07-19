import itertools
def bad(g):
    cuts=[]
    for x in g:
        h=7/6*x; a=h-1/6; a=a if a>0 else 0.0; b=h if h<1 else 1.0
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=0.0; L=0.0
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1-cur>L: L=1-cur
    return L<=1/6+1e-9
def soj(DS,N):
    c=0
    for t in range(N):
        u=t/N; c+= bad(tuple((-d*u)%1.0 for d in DS))
    return c/N
N=1000000
perms=list(itertools.permutations([1,2,3]))
print("2/21 = %.6f"%(2/21))
print("the 6 permutation-rays of (1,2,3), sojourn at N=1e6:")
for p in perms:
    print("  %s: %.6f"%(str(p),soj(p,N)))
print("(3,9,6)=3*(1,3,2): %.6f"%soj((3,9,6),N))
# max over non-proportional-to-(1,2,3) directions in [1,8]^3 at moderate N, then refine top
best=[]
for DS in itertools.product(range(1,9),repeat=3):
    if DS==(1,2,3) or (DS[1]==2*DS[0] and DS[2]==3*DS[0]): continue
    best.append((soj(DS,25200),DS))
best.sort(reverse=True)
print("top non-(1,2,3)-proportional sojourns [1,8]^3 (N=25200):")
for s,d in best[:6]: print("  %s: %.5f  (perm-ray: %s)"%(d,s, tuple(sorted(d))==(1,2,3) or (d[1]%d[0]==0 and sorted(x//d[0] for x in d)==[1,2,3]) ))
