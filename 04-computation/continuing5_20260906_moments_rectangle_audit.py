"""Independent rectangle referee: row polarization, value-to-Bernstein tensor map.
No producer imports. The producer certificate is compared, never used as oracle.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
from pathlib import Path
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
N=0

def need(ok,label):
    global N
    N+=1
    if not ok:raise ArithmeticError(label)

def clean(a):return {k:F(v) for k,v in a.items() if v}
def plus(a,b):
    c=a.copy()
    for k,v in b.items():c[k]=c.get(k,F(0))+v
    return clean(c)
def scale(a,q):return clean({k:v*q for k,v in a.items()})
def times(a,b):
    c={}
    for i,u in a.items():
        for j,v in b.items():c[i+j]=c.get(i+j,F(0))+u*v
    return clean(c)
def power(a,n):
    b={0:F(1)}
    for _ in range(n):b=times(b,a)
    return b
def val(a,s):return sum(c*s**j for j,c in a.items())

def rows(x,y,z):
    beta=dict(zip(range(-1,5),map(F,(1,13,55,x,y,z))))
    c=dict(zip(range(-1,4),map(F,(1,12,45,2*x/3,3*y/7))))
    d=dict(zip(range(-1,4),map(F,(1,11,36,5*x/12,y/7))))
    odd={j:F(comb(14,2*j+1)) for j in range(7)}
    even={j:F(comb(14,2*j)) for j in range(8)}
    carrier=plus(times(odd,odd),times({-1:F(1)},times(even,even)))
    raw=plus(times(beta,beta),scale(times({1:F(1)},times(c,d)),2))
    return clean({j:v*odd.get(j,0) for j,v in beta.items()}),clean({j:v*carrier.get(j,0) for j,v in raw.items()})

def eliminated(x,y):
    # Polarize the quadratic z dependence using three numerical rows.
    raw={z:rows(x,y,F(z))[1] for z in (0,1,-1)}
    zcoeff=[raw[0],scale(plus(raw[1],scale(raw[-1],-1)),F(1,2)),scale(plus(plus(raw[1],raw[-1]),scale(raw[0],-2)),F(1,2))]
    zroot={-1:F(12,7)*y,-2:-x,-3:F(10),-4:F(-1,11)}
    out={}
    for i,coeff in enumerate(zcoeff):
        tpoly={j+1:v*(-1 if j%2 else 1) for j,v in coeff.items()}
        out=plus(out,times(tpoly,power(zroot,i)))
    out=scale(out,F(-11,14))
    need(all(0<=k<=8 for k in out),'all lower Laurent powers cancel, degree at most8')
    return out

def transform(h,a,b):
    if b is None:
        # Horner translation h(a+u).
        out={}
        for e in range(8,-1,-1):out=plus(times(out,{0:a,1:F(1)}),{0:h.get(e,F(0))})
        return out
    A={0:a,1:b};B={0:F(1),1:F(1)}
    out={}
    for e,c in h.items():out=plus(out,scale(times(power(A,e),power(B,8-e)),c))
    return out

def bernstein_from_values(v):return [v[0],2*v[1]-(v[0]+v[2])/2,v[2]]
def power_from_values(v):
    q=(v[2]-2*v[1]+v[0])/2
    return [v[0],v[1]-v[0]-q,q]

def main():
    stem=Path(__file__).with_name('continuing5_20260906_moments_rectangle')
    pins={'.py':'da6236b1828b05bf5d30fe09d97e95005db36fe1baaa390c1b422aba85f13a25','.out':'4aab430109b707da57cd2654ee0138d1b141002a1d67d33164cbe66ba0bb4a28','_certificate.json':'1e86411634a8f62f98a512d741d4b8df8ba5a069fb262d92bafe64e40569a4f8'}
    for suffix,pin in pins.items():
        path=stem.with_suffix(suffix) if suffix[0]=='.' else stem.with_name(stem.name+suffix)
        if not path.exists():path=stem.parent.parent/'05-knowledge'/'results'/path.name
        need(hashlib.sha256(path.read_bytes()).hexdigest()==pin,'frozen producer input pin')
    certpath=stem.with_name(stem.name+'_certificate.json')
    if not certpath.exists():certpath=stem.parent.parent/'05-knowledge'/'results'/certpath.name
    cert=json.loads(certpath.read_text())
    # Reconstruct every raw polynomial coefficient by exact tensor interpolation.
    rawgrid={(i,j,k):rows(F(i),F(j),F(k))[1] for i,j,k in product(range(3),repeat=3)}
    mono={}
    for e in range(-1,9):
        grid={ijk:row.get(e,F(0)) for ijk,row in rawgrid.items()}
        for axis in range(3):
            nxt={}
            for rest in product(range(3),repeat=2):
                keys=[]
                for u in range(3):
                    k=list(rest);k.insert(axis,u);keys.append(tuple(k))
                for k,v in zip(keys,power_from_values([grid[k] for k in keys])):nxt[k]=v
            grid=nxt
        for k,v in grid.items():
            if v:mono[k+(e,)]=v
    expected={tuple(r['powers']):F(r['coefficient']) for r in cert['Q_raw']}
    need(mono==expected,'complete raw Q reconstructed through value polarization')
    need(mono[(0,0,0,-1)]==28 and all(c>0 for c in mono.values()),'full lower carry and raw positive monomials')
    bound=F(0)
    for (i,j,k,e),v in mono.items():
        c=v*(-1 if e%2 else 1)
        x,y,z,s=(F(85),F(36),F(5),F(1,100)) if c>0 else (F(83),F(34),F(0),F(1,102))
        bound+=c*x**i*y**j*z**k*s**(e+1)
    need(bound==F(cert['raw_first_upper']) and bound<-427,'raw first phase upper bound')
    # Every polynomial has separate x/y degrees <=2 by the original bilinear rows.
    hgrid={(i,j):eliminated(F(83+i),F(34+j)) for i,j in product(range(3),repeat=2)}
    for (i,j),h in hgrid.items():
        expected={}
        for r in cert['eliminated_h']:
            a,b,z,e=r['powers'];need(z==0,'certificate has no eliminated z')
            expected[e]=expected.get(e,F(0))+F(r['coefficient'])*(83+i)**a*(34+j)**b
        need(h==clean(expected),'nine-grid complete bidegree identity for h')
        for s in (F(2,17),F(7,5),F(11)):
            zr=F(12,7)*(34+j)/s-(83+i)/s**2+10/s**3-F(1,11)/s**4
            pp,qq=rows(F(83+i),F(34+j),zr)
            need(val(pp,-s)==0 and s*val(qq,-s)==F(-14,11)*val(h,s),'literal original-phase consumer, all carries')
    windows=[(F(1,9),F(13,100)),(F(1),F(8,5)),(F(19,2),None)]
    mins=[]
    for (a,b),record in zip(windows,cert['positive_transforms']):
        grid={ij:transform(h,a,b) for ij,h in hgrid.items()}
        vals=[]
        for k in range(9):
            v=[[grid[i,j].get(k,F(0)) for j in range(3)] for i in range(3)]
            alongx=[bernstein_from_values([v[i][j] for i in range(3)]) for j in range(3)]
            block=[bernstein_from_values([alongx[j][i] for j in range(3)]) for i in range(3)]
            for i,j in product(range(3),repeat=2):
                c=block[i][j];need(c==F(record['coefficients'][k][i][j]) and c>0,'complete independently reconstructed Bernstein coefficient');vals.append(c)
        need(min(vals)==F(record['minimum']),'exact transform minimum');mins.append(str(min(vals)))
    signs=[1,-1,-1,1,1,-1,-1]
    endpoints=list(map(F,('1/102','1/100','1/9','13/100','1','8/5','19/2')))
    for s,sgn,record in zip(endpoints,signs,cert['phase_endpoints']):
        v=[]
        for x,y,z in product((83,85),(34,36),(0,5)):
            pp,_=rows(F(x),F(y),F(z));value=val(pp,-s)/2002
            need(value*sgn>0,'literal first-row brackets at every corner');v.append(value)
        need(min(v)==F(record['min']) and max(v)==F(record['max']),'complete bracket extrema')
    eta=F(8,25)
    derivative=5*eta**4-52*eta**3+165*eta**2-166*eta+36
    need(derivative==F(-146524,78125)<0,'first critical point uniformly before8/25')
    need(F(36)**2/(4*(83-55*eta))==F(540,109)<5,'strict beta fifth-coefficient cap')
    for x,y in product((83,85),(34,36)):
        for v,sgn in zip(map(F,('0','1/10','1','3','5','7')),(-1,1,-1,1,-1,1)):
            need(sgn*(v**5-13*v**4+55*v**3-x*v*v+y*v-1)>0,'whole rectangle beta nonvacuity')
    # Independent charge enumeration: choose number of +15 terms; solve -27 count.
    model=rows(F(84),F(35),F(1))
    sizes=[]
    for mass,expected in zip((14,28),model):
        literal={}
        for c in range(mass+1):
            numerator=mass+14*c
            if numerator%28:continue
            a=numerator//28;b=mass-a-c
            if min(a,b,c)<0:continue
            e=a-mass//14
            literal[e]=literal.get(e,0)+comb(mass,a)*comb(mass-a,c)
        need(literal==expected,'literal charge-solved factorial channels');sizes.append(len(literal))
    need(sizes==[5,10],'all actual channels including doubled lower carry')
    # Distinct formal series inversion by powers of 1-denominator.
    def moments(x,y,z,cflag):
        inv={};base={1:F(13),2:F(-55),3:x,4:-y,5:z}
        for k in range(9):inv=plus(inv,{e:c for e,c in power(base,k).items() if e<=8})
        num={0:F(1),1:F(-12 if cflag else -11),2:F(45 if cflag else 36),3:-x*(F(2,3) if cflag else F(5,12)),4:y*(F(3,7) if cflag else F(1,7))}
        a=times(inv,num);return [a.get(i,F(0)) for i in range(9)]
    def det(a):
        # Permutation determinant: small5x5, independent of producer elimination.
        from itertools import permutations
        out=F(0)
        for p in permutations(range(len(a))):
            sign=(-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
            value=F(sign)
            for i,j in enumerate(p):value*=a[i][j]
            out+=value
        return out
    x,y,s=F(78071,1000),F(601,50),F(57,2)
    z=12*y/(7*s)-x/s**2+10/s**3-F(1,11)/s**4
    pp,qq=rows(x,y,z)
    need(val(pp,-s)==0 and val(qq,-s)>0,'original positive surrogate control')
    for flag,expected in zip((True,False),cert['degree7_hostile_H5_determinants']):
        m=moments(x,y,z,flag);v=det([[m[i+j] for j in range(5)] for i in range(5)])
        need(v==F(expected)<0,'exact degree8 domain rejects surrogate')
    lo,hi=F(16693,2000),F(41733,5000)
    pp,qq=rows(F(84),F(35),F(6))
    need(val(pp,-lo)<0<val(pp,-hi),'z6 positive-response original root exists')
    t=transform(hgrid[1,1],lo,hi)
    need(len(t)==9 and all(v<0 for v in t.values()),'negative h on entire z6 bracket')
    m=moments(F(84),F(35),F(6),False)
    need(det([[m[i+j] for j in range(5)] for i in range(5)])==-25668,'z6 exact-domain hostile')
    print('AUDIT PASS full closed coefficient prism x83..85 y34..36 z0..5; every positive original phase Q<0')
    print('PATH raw row polarization;9-grid same-root Laurent cancellation; value-to-Bernstein dual map;243 exact coefficient matches')
    print('BERNSTEIN_MINIMA',mins,'FIRST_RAW_BOUND',bound)
    print('SCOPE beta nonnegative-root cap540/109; all escaping and weak boundaries; no interlacer needed; outside rectangle remains open')
    print('CONTROLS all15 charge-solved channels, exact phase corners, z1 full rectangle nonvacuity, both domain hostiles')
    print('PASS',N,'always-active exact gates; raw LF normal/-O')

if __name__=='__main__':main()
