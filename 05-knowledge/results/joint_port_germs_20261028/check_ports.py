"""Exact port-register compiler, selector costs, spectra, and hostile controls."""
from pathlib import Path
from itertools import product
from math import comb
import json
import port_core as c
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,m):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(m)
def check_decode(t,ports,shift):
    encoded=c.register(t,ports);N=len(encoded)
    order=tuple(range(shift%N,N))+tuple(range(shift%N))
    shuffled=c.relabel(encoded,order)
    decoded,marks,blocks=c.decode_register(shuffled,len(ports))
    owners={order[v]:i for i,block in enumerate(blocks) for v in block}
    need(set(range(len(t)))<=set(owners),'all old vertices represented')
    need(len(decoded)==len(t),'decoded size')
    need(tuple(owners[p] for p in ports)==marks,'ordered joint port data')
    need(all(c.edge(t,i,j)==c.edge(decoded,owners[i],owners[j]) for i in range(len(t)) for j in range(len(t))),'literal relation recovered')
    need(len(set(owners[i] for i in range(len(t))))==len(t),'no old vertices identified')
def polynomial_unsigned(group,n,ports):
    # Named attachment vertices: one copy of the vertex-permutation action per port.
    out=[0]*(n*ports+1)
    for g in group:
        lengths=[];seen=set()
        for i in range(n):
            if i in seen:continue
            j=i;k=0
            while j not in seen:seen.add(j);k+=1;j=g[j]
            lengths.append(k)
        poly=[1]
        for _ in range(ports):
            for k in lengths:
                q=poly+[0]*k
                for i,a in enumerate(poly):q[i+k]+=a
                poly=q
        for i,a in enumerate(poly):out[i]+=a
    need(all(x%len(group)==0 for x in out),'spectrum integrality')
    return [x//len(group) for x in out]
def automorphism_audit(t,ports):
    rooted=sum(all(g[p]==p for p in ports) for g in c.oracle.isomorphisms(t,t,True))
    factor=3 if t[ports[0]]==0 and ports[0] not in ports[1:] else 1
    actual=len(c.oracle.isomorphisms(c.register(t,ports),c.register(t,ports),True))
    need(actual==factor*rooted,'exact ternary kernel of the unlabelled register')
def path_counts(t):
    n=len(t);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not dp[mask][v]:continue
            for w in range(n):
                if not mask>>w&1 and c.edge(t,v,w):dp[mask|1<<w][w]+=dp[mask][v]
    return dp
def dent_identity(t,p):
    n=len(t);full=(1<<n)-1
    dp=path_counts(t);rev=tuple(sum(c.edge(t,j,i)<<j for j in range(n)) for i in range(n))
    backward=path_counts(rev)
    # When the new vertex is last, the previous endpoint cannot be p.
    total=sum(dp[full])-dp[full][p]
    for S in range(1<<n):
        if S>>p&1:continue
        total+=(sum(dp[S]) if S else 1)*backward[full^S][p]
    return total
def main():
    census=[]
    for n in range(1,8):
        count=1<<((n-1)*(n-2)//2);registers=0
        for mask in range(count):
            t=c.presentation(n,mask)
            for p in range(n):
                check_decode(t,(p,),mask+p);registers+=1
                if n<=4:automorphism_audit(t,(p,))
            if n<=6:
                for p,q in product(range(n),repeat=2):
                    check_decode(t,(p,q),mask+p+q);registers+=1
                    if n<=4:automorphism_audit(t,(p,q))
            if n<=5:
                for p in range(n):
                    need(dent_identity(t,p)==c.oracle.hamilton(c.dent(t,p)),'two-path dent identity')
        census.append({'order':n,'Hamilton_presentations':count,'registered_tuples':registers,
                       'arity':2 if n<=6 else 1})
        print(json.dumps(census[-1]),flush=True)
    pair=[]
    for turn in (0,1):
        a,b=tuple(map(c.tau,c.PAIR)) if turn else c.PAIR
        groups=[c.oracle.isomorphisms(t,t,True) for t in (a,b)]
        need(list(map(len,groups))==[1,1],'rigid seven-pair')
        spectra=[]
        for q in (1,2,3):
            f=[polynomial_unsigned(g,7,q) for g in groups]
            need(f[0]==f[1]==[comb(7*q,k) for k in range(7*q+1)],'all named joint-attachment spectra agree')
            spectra.append({'ports':q,'coefficients':f[0]})
        samples=[]
        for p in range(7):
            A,B=c.blow(a,p),c.blow(b,p)
            samples.append({'port':p,'C3_isomorphic':bool(c.oracle.isomorphisms(A,B)),
                            'C3_H':[c.oracle.hamilton(A),c.oracle.hamilton(B)]})
        pair.append({'native_complement':bool(turn),'spectra':spectra,'C3':samples,
                     'successor_dent_H':[c.oracle.hamilton(c.dent(t,1)) for t in (a,b)]})
        for t in (a,b):
            for ports in product(range(7),repeat=2):check_decode(t,ports,sum(ports))
    need(pair[1]['successor_dent_H']==[307,293],'one-vertex successor separator')
    need(pair[1]['C3'][1]['C3_H']==[777,981],'smallest strong-atom successor separator')
    # One dent without the temporary sink is not an unmarked pointed-core decoder.
    trans=(6,4,0);cycle=(2,4,1)
    need(bool(c.oracle.isomorphisms(c.dent(trans,0),c.dent(cycle,0))),'bare dent decoder hostile')
    witness=c.oracle.hostile()
    hostile=[]
    for rows in (witness['rows_AAA'],witness['rows_BAA']):
        group=c.oracle.isomorphisms(rows,rows,True)
        f=polynomial_unsigned(group,12,1)
        hostile.append({'group_order':len(group),'binary_attachment_spectrum':f,'orbits':sum(f)})
        for ports in ((0,1),(1,5),(4,9),(2,2)):
            check_decode(rows,ports,3)
    need([r['orbits'] for r in hostile]==[176,512],'81/27 attachment orbit control')
    # Repeated ports and terminal-C3 collapse, with arbitrarily long tested registers.
    tails=[]
    for n in range(1,6):
        t=c.presentation(n,(1<<((n-1)*(n-2)//2))-1)
        for r in range(1,7):
            for ports in ((n-1,)*r,tuple(i%n for i in range(r))):
                check_decode(t,ports,r)
                tails.append({'order':n,'ports':ports,'added_vertices':r+1})
    out={'status':'FINITE-EXACT','census':census,'seven_pair':pair,'hostile_81_27':hostile,
         'long_register_controls':tails,'checks':CHECKS}
    (HERE/'ports.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','checks':CHECKS,'seven_pair_dent':[r['successor_dent_H'] for r in pair],
                      'hostile_orbits':[r['orbits'] for r in hostile]}))
if __name__=='__main__':main()
