import collections, sys
def odd(x):
    x=abs(x)
    while x and x%2==0: x//=2
    return x
def ispow2(x): return x>0 and x&(x-1)==0
fam=collections.defaultdict(lambda: None)   # (al, C, D) -> first witness (order, code, x, y); ray variable is the module with zero single coefficient
tot=0
for m in range(3,9):
    for line in open("mf_%d.txt"%m):
        p=line.split(); code=int(p[0]); x,y=int(p[1]),int(p[2]); al,be,ga,de,ex=map(int,p[3:8])
        if not ex: continue
        tot+=1
        # N = al + be 2^a + ga 2^b + de 2^(a+b).  ray in b if be == 0:  N = al + 2^b (ga + de 2^a)
        if be==0 and abs(al)==1 and de!=0:
            k=(al,ga,de)
            if fam[k] is None: fam[k]=(m,code,x,y)
print("exact two-module families read:",tot)
for sign in (1,-1):
    ks=sorted(k for k in fam if k[0]==sign)
    print("\noffset %+d: %d distinct slope forms  ga + de 2^a  (ga != 0: %d)"%(sign,len(ks),sum(1 for k in ks if k[1]!=0)))
    nz=[k for k in ks if k[1]!=0]
    for k in nz:
        al,ga,de=k
        tag=""
        if ispow2(abs(ga)): tag+=" [ga=%s2^%d]"%("-" if ga<0 else "",abs(ga).bit_length()-1)
        if ispow2(abs(de)): tag+=" [de=2^%d]"%(abs(de).bit_length()-1)
        print("  ga=%5d de=%5d  wit=%s  c(a=1..6)=%s%s"%(ga,de,fam[k],[odd(ga+de*2**a) for a in range(1,7)],tag))
    z=sorted(set(odd(k[2]) for k in ks if k[1]==0))
    print("  constant-odd-part families (ga = 0): odd parts of de:",z[:80],"..." if len(z)>80 else "")
