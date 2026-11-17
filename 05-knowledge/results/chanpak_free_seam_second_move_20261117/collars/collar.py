import sys
from pathlib import Path
root=Path(__file__).resolve().parents[2]
for base in ('wheel_square_selector','chanpak_wheel_square_selector_20261117'):
    if (root/base).is_dir(): sys.path.insert(0,str(root/base))
from ray.replay import ends,check_graph
from check_free import D,mask,apex,TT,F

def sub(D,q,idx):
    return tuple(sum((D[i]>>j&1)<<k for k,j in enumerate(idx)) for i in idx),sum((q>>j&1)<<k for k,j in enumerate(idx))

def sig(d,v):
    a=ends(d)[1]
    b=ends(tuple(sum(((x>>i)&1)<<j for j,x in enumerate(d)) for i in range(len(d))))[1]
    M=len(a)-1
    H=ends(d)[0];L=[0,0]; R=[0,0];P=0
    for s in range(M+1):
        j=s.bit_count()%2;h=H[M^s]
        P+=H[s]*h
        L[j]+=sum(z for i,z in enumerate(b[s]) if v>>i&1)*h
        R[j]+=sum(z for i,z in enumerate(a[s]) if not v>>i&1)*h
    return H[-1],P,sum(L),sum(R),L[0]-L[1],R[0]-R[1],ends(apex(d,v))[0][-1]

def join(X,Y):
    h,p,l,r,s,e,k=X
    H,P,L,R,S,E,K=Y
    return h*H,p*P,l*P+h*L,r*H+p*R,h*S,e*H,k*H+h*K-h*H+l*R+s*E

def RC(x,n):
    return ((1<<n)-1)^int(format(x,'0'+str(n)+'b')[::-1],2)

prefix=[3,6,0]
rest=[1,2,4,5]

G={
 (0,0):((4,0,(-2,5)),(4,4,(0,0))),
 (-2,5):((7,3,(3,3)),(7,7,(5,-2))),
 (3,3):((3,3,(3,3)),(3,7,(5,-2))),
 (5,-2):((0,0,(-2,5)),(0,4,(0,0)))
}
mid={(0,0):4,(-2,5):7,(3,3):3,(5,-2):0}
def words(k,state=(0,0)):
    if k==0:
        if state[0]==state[1]:yield ''
    elif k==1:
        yield str(mid[state])
    else:
        for a,b,st in G[state]:
            for w in words(k-2,st):
                yield str(a)+w+str(b)

def require(b):
    if not b: raise ArithmeticError('control failed')
def check_all():
    from itertools import product
    def rho(x):return 7-RC(x,3)
    for l,h in G:
        generated=[]
        for b,a in product(range(8),repeat=2):
            z=4*a+5*rho(b)+l-21
            y=21+8*h-4*b-5*rho(a)
            if z%8==0 and -3<=z//8<=6 and -3<=y<=6:
                generated.append((b,a,(z//8,y)))
        require(set(generated)==set(G[l,h]))
        require([x for x in range(8) if 4*x+5*rho(x)+l==21+8*h]==[mid[l,h]])
    X=sig(*sub(D,mask,prefix));Y=sig(*sub(D,mask,rest))
    require(X==(1,8,4,6,0,0,5) and Y==(5,32,14,14,-2,-2,15))
    for k in range(1,15):
        z=sorted(int(w,8) for w in words(k))
        require(len(z)==1<<((k-1)//2))
        for u in z:
            n=3*k;pw=1<<n;v=RC(u,n)
            require(4*u-5*v==-2*pw+2)
            if k<=2:
                midmask=int(format(u,'0'+str(n)+'b')[::-1],2)
                M=sig(TT(n),midmask)
                h,p,l,r,s,e,K=join(join(X,M),Y)
                require((h,p,l,r,s,e,K)==(5,pw*256,128*pw+32*u+14,128*pw+32*u+14,-2,0,5*F(u)+30*u+64*pw+22+2*(u%2)))
        for a,b in zip(z,z[1:]):
            require(b-a>=2*8**(k//2))
            require(5*(F(b)-F(a))+126*(b-a)+2*((b%2)-(a%2))>=26*(b-a)-7)
    # h=5, P2=256: h5 strong component (8,1,3,6), and
    # exactly three singleton SCCs; scan the four ordinal placements.
    a1=(8,1,3,6)
    group=[]
    for i in range(4):
        a=TT(i);b=TT(3-i);left=len(a);last=left+4
        d=tuple(r|(127^((1<<left)-1)) for r in a)+tuple((r<<left)|(127^((1<<last)-1)) for r in a1)+tuple(r<<last for r in b)
        for q in range(128):
            s=sig(d,q)
            if s[2]==s[3] and abs(s[4]-s[5])==2: group.append((i,q,s))
    require(group==[(0,53,(5,256,142,142,0,-2,91)),(0,74,(5,256,109,109,-5,-3,55)),(3,41,(5,256,142,142,-2,0,91)),(3,86,(5,256,109,109,-3,-5,55))])
    print('Four-state grammar, 14 depths, literal collar profiles, sharp fixed-scale bank pass.')
if __name__=='__main__':check_all()
