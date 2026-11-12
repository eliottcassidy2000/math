"""Exact switching-pressure certificates. Standard library only; no census input."""
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
from pathlib import Path
import json, hashlib
RAYS=((1,1,0),(1,1,1),(0,1,0))
C63=(899022,845140,951686,817725,910828,866012,977578,803895,897084,846064,954747,826899,925006,875228,990417,796949,890369,840096,947131,817052,917862,864493,978358,808178,903903,852118,962488,830751,932264,879572,996810,793469,887049,837144,943370,814261,914331,861366,974341,803287,900481,849080,958602,825542,928620,874016,990614,799016,893714,843082,950924,819862,921420,867658,982409,809933,907353,855178,966406,832527,935936,881683,1000000)
def require(x,msg):
    if not x: raise ValueError(msg)
@lru_cache(None)
def g(p):
    if p==0:return 1
    v=p//2
    return g(v) if p%2 else g(v)+v+v%2
def block(w):
    p=int(w,2);q=sum(1<<i for i,b in enumerate(w) if b=='0');t=1<<len(w)
    return ((g(p),q+1,-1),(p+1,t,0),(1,0,0))
def mv(m,v):return tuple(sum(a*b for a,b in zip(row,v)) for row in m)
def words(L):return ['1'+''.join(mid)+'0' for n in range(2,L+1) for mid in product('01',repeat=n-2)]
def ratio_bounds(m):
    values=[F(mv(m,r)[0],mv(m,r)[1]) for r in RAYS]
    return min(values),max(values)
def inverse_growth(ms,lower):
    qs=[ratio_bounds(m)[0 if lower else 1] for m in ms]
    return [[1/(v[1][1]+v[1][0]*q) for v in ms] for q in qs]
def ratios(K,c):return [sum(a*b for a,b in zip(row,c))/c[i] for i,row in enumerate(K)]
def discover_integer_vector(K):
    # Floating arithmetic proposes only a witness. Every conclusion below rechecks
    # the resulting integer vector against the exact rational matrix.
    x=[1.0]*len(K);kf=[[float(z) for z in row] for row in K]
    for _ in range(100):
        y=[sum(a*b for a,b in zip(row,x)) for row in kf];s=max(y);x=[v/s for v in y]
    return [round(v*1000000) for v in x]
def pisot(m):
    h,q,_=m[0];p,t,_=m[1];c=h*t-p*q+1
    return 1-h-t+c-t<0 and -1-h-t-c-t<0

def main():
    ws63=words(7);ms63=[block(w) for w in ws63]
    require(len(ws63)==63,'alphabet size')
    K=inverse_growth(ms63,True);rs=ratios(K,C63)
    require(all(r<F(973,1000) for r in rs),'subcritical exact witness')
    require(256*64**400<65**400,'fractional exponent bound')
    rho=F(65,64)*F(973,1000)
    chunk=1+F(195,128)*F(max(C63),min(C63))/(1-rho)
    require(chunk<164,'chunk constant')
    # The physical projective bound is checked independently on every cone ray.
    require(all(F(1,4)<ratio_bounds(m)[0]<=ratio_bounds(m)[1]<=1 for m in ms63),'active cone bound')
    extra=['1'+''.join(mid)+'0' for mid in product('01',repeat=6)][:20]
    ws83=ws63+extra;ms83=[block(w) for w in ws83]
    require(all(pisot(m) for m in ms83),'every83 letter individually Pisot')
    L=inverse_growth(ms83,False);d=discover_integer_vector(L);ls=ratios(L,d)
    require(all(r>F(1003,1000) for r in ls),'supercritical exact witness')
    low=json.loads((Path(__file__).with_name('low_oscillation_pressure.json')).read_text())
    wslow=[w for w in words(8) if sum(a!=b for a,b in zip(w,w[1:]))<=3]
    require(low['alphabet']==wslow and len(wslow)==98,'low oscillation alphabet')
    mslow=[block(w) for w in wslow]
    require(all(pisot(m) for m in mslow),'low oscillation individual Pisot')
    lowratios=ratios(inverse_growth(mslow,False),low['vector'])
    require(all(r>F(101,100) for r in lowratios),'low oscillation pressure')
    require(512*1000**4096<1003**4096,'S83 collision exponent')
    require(512*100**1024<101**1024,'S98 collision exponent')
    require([sum(len(w)==n for w in ws83) for n in range(2,9)]==[1,2,4,8,16,32,20],'S83 order polynomial')
    require([sum(len(w)==n for w in wslow) for n in range(2,9)]==[1,2,4,8,15,26,42],'S98 order polynomial')
    sync=json.loads(Path(__file__).with_name('synchronized_pressure151.json').read_text())
    metricpath=Path(__file__).resolve().parents[1]/'switching/metric151.json'
    metric=json.loads(metricpath.read_text())
    require(sync['alphabet']==metric['words'] and len(sync['alphabet'])==151,'synchronized151 shared alphabet')
    require(sync['metric_file_sha256']==hashlib.sha256(metricpath.read_bytes()).hexdigest(),'metric151 binding')
    mssync=[block(w) for w in sync['alphabet']]
    require(all(pisot(m) for m in mssync),'synchronized individual Pisot')
    syncratios=ratios(inverse_growth(mssync,False),sync['vector'])
    require(all(r>F(1069,1000) for r in syncratios),'synchronized151 pressure')
    require(1024*1000**128<1069**128,'synchronized151 collision exponent')
    require([sum(len(w)==n for w in sync['alphabet']) for n in range(2,10)]==[1,2,4,8,15,26,37,58],'synchronized151 order polynomial')
    # Literal graph controls use inherited DP independently of our bounds.
    import importlib.util
    old=Path(__file__).resolve().parents[2]/'session279_interval_tightening_20261112/physical/verify.py'
    spec=importlib.util.spec_from_file_location('literal279',old);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    controls=0
    for word in ws63[:7]:
        for seed in (('10',),('01','0')):
            oldx=module.literal_packet(seed);newx=module.literal_packet(seed+(word,))
            require(mv(block(word),oldx)==newx,'literal physical control')
            q=F(oldx[0],oldx[1]);fac=F(newx[1],oldx[1]);m=block(word)
            require(fac==m[1][1]+m[1][0]*q,'exact W factor');controls+=1
    # Finite endpoint tests: every consecutive factor lies inside its declared
    # last-letter lower/upper growth interval. Not used as all-word proof.
    checks=0;seeds=((1,1,1),(9,9,1),(9,9,3))
    for x in seeds:
        for i,u in enumerate(ms83):
            a=mv(u,x);lo,hi=ratio_bounds(u)
            for v in ms83:
                b=mv(v,a);factor=F(b[1],a[1])
                require(v[1][1]+v[1][0]*lo<=factor<=v[1][1]+v[1][0]*hi,'factor bound');checks+=1
    # Scalar collisions survive the capacity transition. Output examples are
    # physical one/two-letter descendants of a literal rooted singleton.
    seen={};collision=None
    for i,u in enumerate(ms83):
        for j,v in enumerate(ms83):
            x=mv(v,mv(u,(1,1,1)));label=(ws83[i],ws83[j])
            if x[0] in seen and seen[x[0]][0]!=label:
                collision={'count':x[0],'first':seen[x[0]],'second':(label,x)};break
            seen[x[0]]=(label,x)
        if collision:break
    out={'status':'exact rational pressure certificates; no claim of scalar density or Chan-Pak completeness',
         'subcritical':{'alphabet':ws63,'vector':C63,'max_ratio_exact':str(max(rs)),'max_ratio_display':float(max(rs)),'certificate_upper':'973/1000','s':'399/400','rho':str(rho),'chunk_bound_exact':str(chunk),'usable_chunk_bound':164},
         'supercritical':{'alphabet':ws83,'vector':d,'min_ratio_exact':str(min(ls)),'min_ratio_display':float(min(ls)),'certificate_lower':'1003/1000','all_letters_Pisot':True},
         'low_oscillation':{'letters':98,'certificate_lower':'101/100','min_ratio_exact':str(min(lowratios))},
         'synchronized151':{'letters':151,'certificate_lower':'1069/1000','min_ratio_exact':str(min(syncratios)),'collision_exponent':'1/128','metric_file_sha256':sync['metric_file_sha256']},
         'literal_controls':controls,'exact_transition_controls':checks,'scalar_collision':collision,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
