"""Frozen native coefficients and new weighted-edge clock/pin audits.

Independent tuple Gray enumeration audits dynamic-programming coefficients.
The heterogeneous integer weights are an explicit control, not ROOT weights.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from math import comb
import hashlib,json,subprocess,sys,random

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=HERE.parent/'three_target570_20261024'
REV='34b0a4cbce72b6986fabd6d71041448a17a77004'
DEPS={}
for name in ('common.py','bank.py','route.py','twelve.py','menus.json','bank.json'):
    p='05-knowledge/results/three_target570_20261024/'+name
    blob=subprocess.check_output(['git','show',REV+':'+p],cwd=ROOT)
    if not (OLD/name).is_file() or (OLD/name).read_bytes()!=blob:
        raise RuntimeError('Materialize the frozen inherited dependency without changing it: '+p)
    DEPS[p]=hashlib.sha256(blob).hexdigest()
sys.path.insert(0,str(OLD))
from route import Router,add,sub,charge
from twelve import TwoRoots
from common import data

CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(msg)

def dp(labels,weights):
    D=[defaultdict(int) for _ in range(len(labels)+1)];D[0][0]=1
    for i,(g,w) in enumerate(zip(labels,weights)):
        for r in range(i+1,0,-1):
            for h,c in list(D[r-1].items()):D[r][add(g,h)]+=w*c
    return D

def extrema(D,a,b):
    f=[Q(v,comb(b-a,r-a)) for r in range(a,b+1) for v in D[r].values()]
    need(all(len(D[r])==2048 for r in range(a,b+1)),'all group debts in every full band grade')
    return min(f),max(f),max(f)/min(f)

def coefficients():
    R=Router();labels=[R.labels[i] for i in R.R];ids=[R.ids[i] for i in R.R]
    weights=[1+s%3 for s in ids];pin=ids.index(176)
    D=dp(labels,[1]*20);W=dp(labels,weights)
    DD=dp([g for i,g in enumerate(labels) if i!=pin],[1]*19)
    conditional=[]
    for z in (0,1):
        C=[defaultdict(int) for _ in range(21)]
        for r,d in enumerate(DD):
            for g,v in d.items():C[r+z][add(g,labels[pin]) if z else g]=v
        conditional.append(C)
    bins=[defaultdict(int) for _ in range(21)];wbins=[defaultdict(int) for _ in range(21)]
    cbins=[[defaultdict(int) for _ in range(21)] for _ in range(2)]
    y=x=k=prev=0;wt=1
    for j in range(1<<20):
        m=j^(j>>1)
        if j:
            bit=(m^prev).bit_length()-1;sgn=1 if m>>bit&1 else -1
            yy,xx=R.rows[R.R[bit]]['label'];y^=yy;x=(x+sgn*xx)%128;k+=sgn
            if sgn==1:wt*=weights[bit]
            else:
                need(wt%weights[bit]==0,'integer Gray weight division');wt//=weights[bit]
        g=(y,x);bins[k][g]+=1;wbins[k][g]+=wt;cbins[(m>>pin)&1][k][g]+=1;prev=m
    for r in range(21):
        need({(g&15,g>>4):v for g,v in D[r].items()}==bins[r],'full tuple count census')
        need({(g&15,g>>4):v for g,v in W[r].items()}==wbins[r],'independent nonuniform weighted census')
        for z in (0,1):need({(g&15,g>>4):v for g,v in conditional[z][r].items()}==cbins[z][r],'conditioned tuple count census')
    lo,hi,kap=extrema(D,5,15);wlo,whi,wkap=extrema(W,5,15)
    F2=[Q(v,comb(10,r-5)*2**(r-5)) for r in range(5,16) for v in W[r].values()]
    kap2=max(F2)/min(F2)
    need(kap2==Q(13048,115) and wkap/kap2>240,'matched auxiliary activity improves bound without changing native weights')
    need((lo,hi,kap)==(Q(67,252),Q(16),Q(4032,67)),'exact unweighted comparison constant')
    cs=[extrema(C,7,13) for C in conditional];clo=min(x[0] for x in cs);chi=max(x[1] for x in cs)
    native=[]
    for m in (60,61):
        H=1+9*m+4*comb(m,2)+2*comb(m,3)+comb(m,4)+comb(m,5)
        gamma=1/(2*kap*comb(m+10,2))
        native.append({'outside_sources':m,'H_uniform':H,'interior_gap':str(gamma),'full_gap':str(gamma/(8*H))})
    return {'enumerated_subsets':1<<20,'base':{'min':str(lo),'max':str(hi),'kappa':str(kap)},
            'moving_pin':{'min':str(clo),'max':str(chi),'kappa':str(chi/clo)},
            'explicit_weight_fixture':{'reservoir_V_odds':dict(zip(map(str,ids),weights)),
                'min':str(wlo),'max':str(whi),'kappa':str(wkap),'all_U_weights':1,
                'auxiliary_odds_2_kappa':str(kap2)},
            'native_uniform_bounds':native}

def events():
    rng=random.Random(57054101);_,_,_,X,Y=data();chars=[X]+[[64*y for y in row] for row in Y[:4]]
    reports=[]
    for remove in (False,True):
        for mode in ('base','one_pin','twelve_pins'):
            C=TwoRoots(remove) if mode=='twelve_pins' else None
            R=C.R if C else Router(remove,mode=='one_pin');outside=C.free if C else R.E
            a,b=(5,15) if mode=='base' else (7,13)
            paths={s:{o['target']:o['path'] for o in row['options']} for s,row in zip(R.ids,R.rows)}
            stats=defaultdict(lambda:{'events':0,'max_clocks':0,'max_support':0})
            def fill(J,r,q,z):return C.fill(J,r,q,z) if C else R.replace(J,r,q)
            def check_pair(x,y,kind):
                need(x.bit_count()==y.bit_count() and R.q(x)==R.q(y),'full integer margin and exact group debt')
                for v in (x,y):
                    if C:
                        need(C.valid(v),'all twelve bit constraints')
                        for g in C.guards:
                            s,h=g['root'],g['helper']
                            pp=paths[s][2234 if v>>R.ids.index(s)&1 else 1965]
                            qq=paths[h][2234 if v>>R.ids.index(h)&1 else 1965]
                            need(C.F['embeddings'][str(s)][pp[g['root_clock']]]==C.F['embeddings'][str(h)][qq[g['helper_clock']]],'literal shared owner, not only phase code')
                    elif mode=='one_pin':
                        p=paths[176][2234 if v>>R.ids.index(176)&1 else 1965]
                        q=paths[368][2234 if v>>R.ids.index(368)&1 else 1965]
                        need(p[2]==[8,9,10,2,6,11,12,13][q[1]],'single literal moving-owner equation')
                clocks=0
                for i,s in enumerate(R.ids):
                    if not (x^y)>>i&1:continue
                    p=paths[s][2234 if x>>i&1 else 1965];q=paths[s][2234 if y>>i&1 else 1965]
                    clocks+=sum(u!=v for u,v in zip(p,q))
                    need(all(p[j]==q[j]==v for j,v in R.rows[i]['pins']),'inherited fixed indexed pins')
                    need(sum(u!=v for u,v in zip(p[::-1],q[::-1]))==sum(u!=v for u,v in zip(p,q)),'converse has identical integer clock work')
                support=(x^y).bit_count();budget=(30,101) if C else (22,54)
                need(support<=budget[0] and clocks<=budget[1],'weighted interior edge budget')
                for chi in chars:
                    debt=defaultdict(int)
                    for i,s in enumerate(R.ids):
                        if not (x^y)>>i&1:continue
                        u=2234 if x>>i&1 else 1965;v=2234 if y>>i&1 else 1965
                        old=2*(chi[s]-chi[u]);new=2*(chi[s]-chi[v])
                        debt[s]+=new-old;debt[u]+=old;debt[v]-=new
                    need(all(v%256==0 for v in debt.values()),'authentic physical current at every source and target owner')
                st=stats[kind];st['events']+=1;st['max_clocks']=max(st['max_clocks'],clocks);st['max_support']=max(st['max_support'],support)
            for r in range(a,b+1):
                for trial in range(48):
                    chosen=rng.sample(outside,len(outside)//2);J=sum(1<<i for i in chosen)
                    q=rng.randrange(2048);z=rng.randrange(2);x=fill(J,r,q,z)
                    i=rng.choice(chosen);j=rng.choice([i for i in outside if not J>>i&1]);K=J^(1<<i)^(1<<j)
                    check_pair(x,fill(K,r,q,z),'outside_swap')
                    if r>a:check_pair(x,fill(J|(1<<j),r-1,q,z),'outside_add')
                    if r<b:check_pair(x,fill(J^(1<<i),r+1,q,z),'outside_remove')
                    if C:check_pair(x,fill(J,r,q,1-z),'root176_refresh')
            reports.append({'binary_sources':R.n,'mode':mode,'edge_probes':dict(stats)})
    return reports

def main():
    co=coefficients();ev=events()
    return {'status':'PASS exact native coefficient census and sampled full weighted-edge replay',
            'checks':CHECKS,'frozen_revision':REV,'dependency_sha256':DEPS,'coefficients':co,'native_events':ev,
            'scope':'entire 2^20 coefficient census; edge probes sample all update types and band grades; proof pays all representatives'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
