"""Uniform shallow-horizon budget from guarded inverse-clock counting.

The proof is in README.md. Counts overcount invalid words; no arrival is
assumed, and no negative interaction is spent by this certificate.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from math import comb
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('shallow_native',HERE.parent/'session18_native_annulus/verify_annulus.py')
native=module_from_spec(spec);sys.modules[spec.name]=native;spec.loader.exec_module(native)
CHECKS=0

def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def integer(x,lo=0):
    if type(x) is not int or x<lo:raise ValueError('literal nonnegative integer')
    return x

def word_bound(K,d):
    integer(K);integer(d)
    return sum((1<<r)*comb((K+8+r)//6+r,r) for r in range(d+1))

def energy_bound(K,d,ports=627):
    integer(ports,1)
    return 3*ports*word_bound(K,d)**2

def tail_certificate(K,d,ports=627):
    integer(K);integer(d);integer(ports,1)
    # The six residue classes of K each advance by six. Every summand
    # increases by at most 1+d/(floor((K+9)/6)+1), which is at most eight.
    h=(K+9)//6
    monotone=d<=7*(h+1)
    rows=[[k,word_bound(k,d),energy_bound(k,d,ports),153*(1<<k)]
          for k in range(K,K+6)]
    return monotone and all(row[2]<=row[3] for row in rows),rows

def inverse(u,a,j):
    integer(u,1);integer(a,1);integer(j)
    if a not in (1,2) or not native.source(u):raise ValueError('primitive parent and branch')
    num=(1<<a)*(3*u+1)*4**j-((1<<a)+3)
    if num%9:raise ValueError('inverse integrality')
    s=num//9
    if not native.source(s):raise ValueError('primitive child')
    need(native.edge(s)==(u,a+2*j),'literal inverse edge and clock')
    return s

def rejects(fn):
    try:fn()
    except ValueError:return True
    return False

def finite_arrivals(N,D,V):
    """Independent forward prefixes, never require an unarrived route to finish."""
    W=[{} for _ in range(D+1)];count=[[0]*len(V) for _ in range(D+1)]
    index={u:i for i,u in enumerate(sorted(V))}
    largest_clock=[0]*(D+1)
    for s in range(1,256*N,2):
        if not native.source(s):continue
        u=s;t=0
        for d in range(D+1):
            if u in V:
                q=native.charge(64*N,s)
                need(native.Q(q)<=3,'universal source charge alphabet')
                signed=tuple(x if t%2==0 else -x for x in q)
                for h in range(d,D+1):
                    old=W[h].get(u,(0,0));W[h][u]=tuple(a+b for a,b in zip(old,signed))
                    count[h][index[u]]+=1;largest_clock[h]=max(largest_clock[h],t)
                need((1<<t)*u<=4**d*s,'actual multiplicative clock bound')
                need((s%15)%5==(s%35)%5,'compatible source registers')
                break
            if d<D:
                u,a=native.edge(u);t+=a
    K=N.bit_length()-1
    rows=[]
    for d in range(D+1):
        E=sum(native.Q(z) for z in W[d].values())
        need(max(count[d])<=word_bound(K,d),'per-terminal first-arrival word bound')
        need(E<=energy_bound(K,d),'forward-prefix energy bound')
        rows.append([d,sum(count[d]),max(count[d]),largest_clock[d],E,energy_bound(K,d)])
    return rows

def run():
    fb=(HERE.parent/'session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n',b'\n')
    need(sha256(fb).hexdigest()=='391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f','fixed grounded frontier')
    obj=json.loads(fb);routes=native.Routes(obj['frontier'],dict(obj['frontier_ranks']));V=routes.V
    for r in range(3):
        for length in range(12):
            q=tuple(sum(native.L[(r+i)%3][c] for i in range(length)) for c in range(2))
            for upper in (False,True):
                z=tuple(q[c]+(native.U[(r+length)%3][c] if upper else 0) for c in range(2))
                need(native.Q(z)<=3,'all reduced lower/upper palette cases')
    for u in range(1,1000,2):
        if not native.source(u):continue
        for a in (1,2):
            legal=[]
            for j in range(12):
                num=(1<<a)*(3*u+1)*4**j-((1<<a)+3)
                if num%9==0:
                    inverse(u,a,j);legal.append(j)
            need(len(legal)==4 and all(j%3==legal[0]%3 for j in legal),'one inverse residue class')
    for K in range(8,80):
        for d in range(13):
            h=(K+9)//6
            need(Fraction(word_bound(K+6,d),word_bound(K,d))<=1+Fraction(d,h+1),'six-step growth majorant')
    thresholds=[]
    for d in range(1,21):
        K=0
        while not tail_certificate(K,d)[0]:K+=1
        need(K==0 or not tail_certificate(K-1,d)[0],'least threshold for declared six-band certificate')
        thresholds.append({'depth':d,'K':K,'N_lower':1<<K,'six_bands':tail_certificate(K,d)[1]})
    need(rejects(lambda:word_bound(True,1)),'Boolean scale rejected')
    need(rejects(lambda:tail_certificate(12,1.0)),'floating horizon rejected')
    need(not tail_certificate(2,8)[0],'N7 deep hostile is not certified')
    need(Fraction(1881*9,153*4)**3<2**16,'exact cubed constant for infinite logarithmic cone')
    for d in range(65):
        need(4*(d+1)**2<=9*(1<<d),'quadratic prefactor majorant controls')
        need(energy_bound(8*d+12,d)<=153*(1<<(8*d+12)),'closed-form cone boundary control')
    controls={str(N):finite_arrivals(N,d,V) for N,d in ((7,8),(64,3),(4096,1))}
    return {'status':'PROVED uniform shallow-horizon budget; all-horizon bound OPEN',
            'frontier_size':len(V),'closed_cone':'all d>=0, all N>=2^(8d+12)',
            'thresholds':thresholds,'finite_forward_controls':controls,'checks':CHECKS}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run();data=(json.dumps(result,indent=2)+'\n').encode();path=HERE/'receipt.json'
    if a.write:path.write_bytes(data)
    else:need(path.read_bytes().replace(b'\r\n',b'\n')==data,'immutable receipt comparison')
    print(json.dumps({'checks':CHECKS,'sha256_lf':sha256(data).hexdigest(),
                      'thresholds':[[x['depth'],x['K']] for x in result['thresholds']]}))
