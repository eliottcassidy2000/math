#!/usr/bin/env python3
"""Independent exact moving-centre Rouché prototype on one frozen braid path.

Only rational witness data and a literal, hash-pinned coefficient table are
inherited. No old root-continuation or braid-extraction code is imported.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import ast,gzip,hashlib,json,time,sys

class C:
    __slots__=('r','i')
    def __init__(self,r=0,i=0):self.r,self.i=Q(r),Q(i)
    def __add__(self,o):
        if not isinstance(o,C):o=C(o)
        return C(self.r+o.r,self.i+o.i)
    __radd__=__add__
    def __neg__(self):return C(-self.r,-self.i)
    def __sub__(self,o):return self+-o
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        if not isinstance(o,C):return C(self.r*o,self.i*o)
        return C(self.r*o.r-self.i*o.i,self.r*o.i+self.i*o.r)
    __rmul__=__mul__
    def __eq__(self,o):return isinstance(o,C) and self.r==o.r and self.i==o.i
    def A(self):return abs(self.r)+abs(self.i)
    def B(self):return max(abs(self.r),abs(self.i))
    def wire(self):return [str(self.r),str(self.i)]

GATES=0
def check(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
def powers(z,n):
    out=[C(1)]
    for _ in range(n):out.append(out[-1]*z)
    return out

COEFS={(6,0):Q(1),(5,0):Q(3,2),(4,1):Q(6),(4,0):Q(91,48),
 (3,2):Q(-2),(3,1):Q(145,18),(3,0):Q(925,162),(2,2):Q(15,2),
 (2,1):Q(-1,3),(2,0):Q(337,54),(1,3):Q(-6),(1,2):Q(38,3),
 (1,1):Q(-2),(1,0):Q(17,9),(0,4):Q(1),(0,3):Q(-4),(0,2):Q(34,9)}

def taylor(coefs,u,z):
    pu=powers(u,max(i for i,k in coefs));pz=powers(z,max(k for i,k in coefs))
    out={}
    for (i,k),c in coefs.items():
        for j in range(i+1):
            for l in range(k+1):
                key=(j,l)
                out[key]=out.get(key,C())+pu[i-j]*pz[k-l]*(c*comb(i,j)*comb(k,l))
    return out

def moving(a,du,dz):
    pu=powers(du,max(j for j,k in a));pz=powers(dz,max(k for j,k in a))
    out={}
    for (j,k),c in a.items():
        if c==C():continue
        for l in range(k+1):
            key=(j+k-l,l)
            out[key]=out.get(key,C())+c*pu[j]*pz[k-l]*comb(k,l)
    return out

def margin(coeffs,radius):
    left=coeffs.get((0,1),C()).B()*radius
    right=sum((v.A()*radius**k for (j,k),v in coeffs.items() if (j,k)!=(0,1)),Q(0))
    return left-right

def affine_distance(d0,delta):
    # The L-infinity norm is the maximum of four affine functions. All
    # changes of active affine function occur at the candidates below.
    cand={Q(0),Q(1)}
    for v0,dv in [(d0.r,delta.r),(d0.i,delta.i),
                  (d0.r-d0.i,delta.r-delta.i),(d0.r+d0.i,delta.r+delta.i)]:
        if dv:
            h=-v0/dv
            if 0<=h<=1:cand.add(h)
    return min((d0+delta*h).B() for h in cand)

def segment(a,z0,du,z1):
    distances={}
    for i in range(4):
        for j in range(i):
            distances[i,j]=affine_distance(z0[i]-z0[j],z1[i]-z0[i]-z1[j]+z0[j])
    if not all(distances.values()):return None
    rr=[min(d for (i,j),d in distances.items() if h in(i,j))/16 for h in range(4)]
    for (i,j),distance in distances.items():
        if not rr[i]+rr[j]<distance:return None
    mm=[]
    for i in range(4):
        m=margin(moving(a[i],du,z1[i]-z0[i]),rr[i])
        if m<=0:return None
        mm.append(m)
    return rr,mm

def word(points):
    rotate=C(1,Q(1,4));out=[]
    current=sorted(range(4),key=lambda i:(points[0][i]*rotate).r)
    check(current==list(range(4)),'declared initial projected order')
    for row0,row1 in zip(points,points[1:]):
        row0=[z*rotate for z in row0];row1=[z*rotate for z in row1]
        events=[]
        for i in range(4):
            for j in range(i):
                d0=row0[i]-row0[j];d1=row1[i]-row1[j]
                check(d0.r!=0 and d1.r!=0,'generic projected endpoints')
                if d0.r*d1.r<0:events.append((d0.r/(d0.r-d1.r),i,j))
        check(len(set(h for h,i,j in events))==len(events),'distinct exact crossing times')
        for h,i,j in sorted(events):
            p=current.index(i);q=current.index(j)
            check(abs(p-q)==1,'adjacent projected strands')
            p=min(p,q);left,right=current[p:p+2]
            gap=(row0[left]-row0[right])*(1-h)+(row1[left]-row1[right])*h
            check(gap.i!=0,'crossing strands remain distinct')
            letter=(p+1)*(1 if gap.i<0 else -1)
            if out and out[-1]==-letter:out.pop()
            else:out.append(letter)
            current[p:p+2]=[right,left]
    return out

def multiply(p,q):
    out={}
    for (i,j),v in p.items():
        for (k,l),w in q.items():out[i+k,j+l]=out.get((i+k,j+l),Q(0))+v*w
    return {k:v for k,v in out.items() if v}

def controls():
    f={(0,0):Q(1)}
    for a in range(4):f=multiply(f,{(0,1):Q(1),(1,0):Q(-1),(0,0):Q(-a)})
    delta=C(10**6,10**6);zs=[C(k) for k in range(4)]
    aa=[taylor(f,C(),z) for z in zs]
    for a in aa:
        m=moving(a,delta,delta)
        check(all(v==C() for(j,k),v in m.items() if j>0),'all common-translation terms cancel')
        check(margin(m,Q(1,16))>0,'arbitrarily long translation exact tube')
    check(segment(aa,zs,delta,[z+delta for z in zs]) is not None,'full labelled translation segment')
    check(affine_distance(C(1),C(-2))==0,'endpoint-separated swapped strands collide')
    check(segment(aa,zs,C(),[zs[1],zs[0],zs[2],zs[3]]) is None,'swapped endpoint labels rejected')
    check(affine_distance(C(1,0),C(-1,1))==Q(1,2),'interior Re=Im minimum')
    check(affine_distance(C(2,-1),C(0,2))==2,'flat maximum interval')
    # Vanishing first-order residual alone cannot control the whole path.
    f={(0,1):Q(1),(2,0):Q(-1),(3,0):Q(1)}
    for k in [2,3,4]:f=multiply(f,{(0,1):Q(1),(0,0):Q(-k)})
    a=taylor(f,C(),C());m=moving(a,C(1),C())
    check(m.get((1,0),C())==C(),'curvature hostile has C10 zero')
    check(margin(m,Q(1,16))<0,'higher terms reject curvature hostile')
    check(Q(4,27)>Q(1,16),'actual intermediate root exits endpoint disk')

def main():
    started=time.perf_counter();controls()
    root=Path(__file__).resolve().parents[1]
    old=root/'04-computation/planar_jc48_sep08_two_cusp_braid.py'
    check(hashlib.sha256(old.read_bytes()).hexdigest()=='d5f0d1ce16b439c13d936f35ce38f558b1c5210f892a11218ab7f4f4695198e1','frozen polynomial producer pin')
    def literal(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int):return n.value
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -literal(n.operand)
        if isinstance(n,ast.Tuple):return tuple(literal(v) for v in n.elts)
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Q':
            return Q(*(literal(v) for v in n.args))
        raise RuntimeError('unexpected coefficient literal')
    tables=[n.value for n in ast.parse(old.read_text()).body if isinstance(n,ast.Assign)
            and any(isinstance(t,ast.Name) and t.id=='COEFS' for t in n.targets)
            and isinstance(n.value,ast.Dict) and n.value.keys]
    inherited={literal(k):literal(v) for k,v in zip(tables[-1].keys,tables[-1].values)}
    check(COEFS==inherited,'complete literal coefficient table independently checked')
    check(max(k for i,k in COEFS)==4 and {(i,k):v for(i,k),v in COEFS.items() if k==4}=={(0,4):1},'exact monic quartet degree')
    cert=root/'05-knowledge/results/planar_jc48_sep08_two_cusp_braid_smooth_certificate.json.gz'
    compressed=cert.read_bytes()
    check(hashlib.sha256(compressed).hexdigest()=='3ae1dfc0ada5d2fe148cb4da1b4643fcfd10280e66ac8a73ad783b6a007b22b9','frozen input gzip pin')
    raw=gzip.decompress(compressed)
    check(hashlib.sha256(raw).hexdigest()=='06d050671046d9145ce23c82fa7fbbbeaa2d8a73ba62fd664b74cdc8046820be','frozen input raw pin')
    data=json.loads(raw)
    check(all(isinstance(v,list) and len(v)==2 for v in data['u']),'every base centre has two coordinates')
    check(all(isinstance(row,list) and len(row)==4 and
              all(isinstance(v,list) and len(v)==2 for v in row) for row in data['z']),
          'every fibre has exactly four two-coordinate centres')
    us=[C(*v) for v in data['u']];zs=[[C(*v) for v in row] for row in data['z']]
    check(len(us)==len(zs)==9554 and data['word']==[3,1,-3],'named input universe')
    base=C(4,3);rad=Q(1,32)
    corners=[base,C(rad),C(0,rad),C(-rad),C(0,-rad),C(rad),base]
    edge_indices=[0];previous=0
    for endpoint in corners[1:]:
        hit=next(i for i in range(previous+1,len(us)) if us[i]==endpoint)
        edge_indices.append(hit);previous=hit
    check(edge_indices[-1]==len(us)-1,'complete six-edge input')
    for left,right in zip(edge_indices,edge_indices[1:]):
        delta=us[right]-us[left];last=Q(0)
        for i in range(left+1,right+1):
            h=(us[i].r-us[left].r)/delta.r if delta.r else(us[i].i-us[left].i)/delta.i
            check(last<h<=1 and us[i]==us[left]+delta*h,'preserved original base-edge parameter')
            last=h
    selected=[0];accepted=[];attempts=0;hint=32
    for edge,(left,right) in enumerate(zip(edge_indices,edge_indices[1:])):
        i=left
        while i<right:
            aa=[taylor(COEFS,us[i],z) for z in zs[i]];cache={}
            def attempt(j):
                nonlocal attempts
                if j not in cache:
                    attempts+=1;cache[j]=segment(aa,zs[i],us[j]-us[i],zs[j])
                return cache[j]
            lo=i;hi=min(right,i+hint)
            if attempt(hi) is not None:
                lo=hi
                while lo<right:
                    hi=min(right,i+2*(lo-i))
                    if attempt(hi) is None:break
                    lo=hi
            while hi-lo>1:
                mid=(lo+hi)//2
                if attempt(mid) is None:hi=mid
                else:lo=mid
            if lo==i:
                check(attempt(i+1) is not None,'even the adjacent original step must certify')
                lo=i+1
            certseg=attempt(lo)
            check(certseg is not None,'accepted segment strict inequalities')
            selected.append(lo)
            accepted.append([i,lo,[[str(r),str(m)] for r,m in zip(*certseg)]])
            hint=max(2,lo-i);i=lo
        if '--progress' in sys.argv:print('edge',edge,'moving segments',len(selected)-1,'attempts',attempts,flush=True)
    result=word([zs[i] for i in selected]);check(result==data['word'],'identical complete exact braid word')
    check(us[selected[0]]==us[selected[-1]]==base,'closed original base loop')
    check(sorted(z.wire() for z in zs[0])==sorted(z.wire() for z in zs[-1]),'closed unordered fibre')
    semantic=hashlib.sha256(json.dumps(accepted,separators=(',',':')).encode()).hexdigest()
    print('Exact moving-centre tubes PASS; one complete inherited base path')
    print('Input segments:',len(us)-1,'Moving segments:',len(selected)-1,'Exact segment attempts:',attempts)
    print('Preserved word:',result,'Always-active gates:',GATES)
    print('Controls: unbounded common translation; swapped-label collision; zero-linear-residual curvature hostile')
    print('Semantic SHA256:',semantic)
    if '--progress' in sys.argv:print('Wall seconds:',round(time.perf_counter()-started,3))

if __name__=='__main__':main()
