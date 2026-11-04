"""Independent geometric enumeration and boundary checks; runs under -O."""
import importlib.util, random
spec=importlib.util.spec_from_file_location('detector','/workspace/shared/session178/detector.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def check(c,msg):
    if not c: raise RuntimeError(msg)
# Direct shell walker: no inverse formula, and adjacency tested using cube norm.
coords={1:(0,0)}
index=2
for r in range(1,201):
    x,y=r,0
    for dx,dy in [(0,-1),(-1,0),(-1,1),(0,1),(1,0),(1,-1)]:
        for _ in range(r):
            coords[index]=(x,y);index+=1;x+=dx;y+=dy
check(len(set(coords.values()))==len(coords),'walker collision')
for n,p in coords.items():
    check(m.xy(n)==p,('xy mismatch',n))
    check(m.label(*p)==n,('label mismatch',n))
edges=[]
for n in range(1,40001):
    target=n//2 if n%2==0 else 3*n+1
    p,q=coords[n],coords[target]
    dx,dy=p[0]-q[0],p[1]-q[1]
    if max(abs(dx),abs(dy),abs(dx+dy))==1: edges.append((n,target))
expected=[(1,4),(2,1),(3,10),(28,14),(30,15),(32,16),(34,17),(36,18)]
check(edges==expected,edges)
for r in [1,2,6,100,10**10,10**100]:
    base=2+3*r*(r-1)
    check(m.xy(base)==(r,0),('start',r))
    check(m.xy(base+6*r-1)==(r-1,1),('end',r))
    for j in [0,r-1,2*r,5*r,6*r-1]:check(m.label(*m.xy(base+j))==base+j,('large',r,j))
rng=random.Random(178)
for bits in [1,2,17,100,1024,4096]:
    for _ in range(10):
        n=rng.getrandbits(bits) or 1
        for A,B in [(2,0),(3,1)]:
            y=A*n+B;c=m.certify(n,y,A,B)
            check(c is not None,'missing certificate')
            check(c[0]==B and c[-1]==0,'boundary')
            for i in range(len(c)-1):check(A*((n>>i)&1)+c[i]==((y>>i)&1)+2*c[i+1],'cell')
            check(m.certify(n,y+1,A,B) is None,'false acceptance')
print('PASS independent: 120601 spiral sites; 40000 geometric edges; huge shells; 120 long carry certificates and false cases.')
