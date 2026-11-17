import sys
from pathlib import Path

root=Path(__file__).resolve().parents[2]
for name in ('wheel_square_selector','chanpak_wheel_square_selector_20261117'):
    if (root/name).is_dir(): sys.path.insert(0,str(root/name))
from ray.replay import check, check_graph, ends, strong
from check_free import F, TT, apex

D=(52,365,176,69,266,280,437,315,13)
M=(113,43,39)
def profile(d,q):
    f=(1<<len(d))-1
    h,fin=ends(d)
    back=tuple(sum(((x>>i)&1)<<j for j,x in enumerate(d)) for i in range(len(d)))
    _,beg=ends(back)
    w=wmin=r=rmin=P=0
    for s in range(f+1):
        h0=h[f^s]
        sign=(-1)**s.bit_count()
        a=h0*sum(c for i,c in enumerate(beg[s]) if q>>i&1)
        b=h0*sum(c for i,c in enumerate(fin[s]) if not(q>>i&1))
        P+=h0*h[s]
        w+=a;wmin+=sign*a;r+=b;rmin+=sign*b
    return h[-1],P,w,r,wmin,rmin,ends(apex(d,q))[0][-1]

def ordinal(a,b):
    n=len(a); return tuple(x|(((1<<len(b))-1)<<n) for x in a)+tuple(x<<n for x in b)

def calc(p,u,w):
    H,P,L,R,Lm,Rm,k=p
    return k+H*(F(u)+F(w)-2)+R*u+L*w+P*u*w-Rm*(u&1)-Lm*(w&1)

def literal(d,q,u,w):
    s=bin(u)[2:]
    t=''.join('1' if z=='0' else '0' for z in reversed(bin(w)[2:]))
    a=ordinal(ordinal(TT(len(s)),d),TT(len(t)))
    m=sum(int(b)<<i for i,b in enumerate(s))+(q<<len(s))+sum(int(b)<<(len(s)+len(d)+i) for i,b in enumerate(t))
    return apex(a,m)

def main():
    check_graph(D)
    for i,m in enumerate(M):
        p=(1521,26412,10422,12894,-750,-890,7285+2*i)
        check(profile(D,m),p)
        for u,w in [(1,1),(1,2),(2,1),(2,2),(3,5),(4,3),(5,4),(4,4)]:
            t=literal(D,m,u,w)
            check(strong(t),True)
            check(ends(t)[0][-1],calc(p,u,w))
        print(m,p)
    print("same literal centers; all four parity orientations; 24 raw graph counts")
if __name__=='__main__':main()
