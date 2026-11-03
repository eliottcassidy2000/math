"""Independent exact checks of the return catalog, cones and compatible gates."""
from collections import Counter
from fractions import Fraction as F
from math import gcd
import json,sympy as sp
from support import *

def nu3(n):
    n=abs(int(n));need(n>0,'nonsingular quotient');e=0
    while n%3==0:n//=3;e+=1
    return e

def clock(beta,m,c):
    start=tuple(x%m for x in beta);z=start;seen={};hits=[]
    while z not in seen:
        k=len(seen);seen[z]=k
        if (z[0]+c*z[1])%m==0:hits.append(k)
        z=((-3*z[0]+2*z[1])%m,(8*z[0]-5*z[1])%m)
    need(z==start,'invertible pure clock')
    return {'modulus':m,'period':len(seen),'accepted_indices':hits}

def main():
    data=json.loads((HERE/'return_census.json').read_text());counts={T:1}
    # This independent census merges rows, retaining multiplicities instead of words.
    for declared in data['census']:
        nxt=Counter();r=0
        for v,c in counts.items():
            for p in OLD_MATS.values():
                z=step(v,p)
                if min(z)>=0 or max(z)<=0:continue
                nxt[z]+=c
                if z[2]==0 and z[0]==-2*z[1]:r+=c
        need((len(nxt),sum(nxt.values()),r)==(declared['mixed_rows'],declared['mixed_words'],declared['return_words']),'independent complete return census')
        counts=nxt
    records=data['returns']+[{'word':['s78'],'matrix':MATS['s78'],'bulk_matrix':quotient(matrix(('s78',))),'gain':1}]
    B=sp.Matrix([[0,1],[0,2],[1,-1]]);cert=[];intervals=[]
    for r in records:
        p=sp.eye(3)
        for l in r['word']:p=p*sp.Matrix(ATLAS[l]['matrix'])
        need(p==sp.Matrix(r['matrix']) and sp.Matrix([T])*p==r['gain']*sp.Matrix([T]),'independent full product and error')
        q=sp.Matrix(r['bulk_matrix']);need(p*B==B*q,'full quotient identity')
        need(all(x>0 for x in q),'positive quotient')
        e=nu3(q.det());mu=min(sum(q.row(i)) for i in range(2));need(mu>=5*3**e,'strict growth versus determinant valuation')
        lo,hi=sorted(F(int(q[0,1]+i*q[1,1]),int(q[0,0]+i*q[1,0])) for i in (0,3))
        intervals.append((lo,hi));cert.append({'word':r['word'],'v3_det':e,'minimum_row_sum':int(mu),'margin':str(F(int(mu),3**e))})
    intervals.sort();merged=[]
    for lo,hi in intervals:
        if not merged or lo>merged[-1][1]:merged.append([lo,hi])
        else:merged[-1][1]=max(hi,merged[-1][1])
    gaps=[];clocks=[]
    for v in data['seed_states']:
        y,x=bulk(v);ratio=F(x,y)
        need(not any(lo<=ratio<=hi for lo,hi in merged),'all-depth terminal cone exclusion')
        lower=max((hi for lo,hi in merged if hi<ratio),default=F(0))
        upper=min((lo for lo,hi in merged if lo>ratio),default=F(3))
        gaps.append({'state':v,'ratio':str(ratio),'excluded_gap':[str(lower),str(upper)]})
        clocks.append({'state':v,'gate49':clock((y,x),49,-21),'gate100':clock((y,x),100,50)})
    qB=((5,8),(2,3));R49=((3395,6688),(1232,2427));R100=((238350,336547),(89350,126161))
    v=(295,216,260);beta=bulk(v);joint=[];z=tuple(x%4900 for x in beta);start=z;k=0
    while True:
        if (z[0]-21*z[1])%49==0 and (z[0]+50*z[1])%100==0:joint.append(k)
        z=((-3*z[0]+2*z[1])%4900,(8*z[0]-5*z[1])%4900);k+=1
        if z==start:break
    need(k==1680 and joint==[9,849],'compatible joint clock with shorter acceptance period')
    # True mixtures of the two raisers retain an additional noncommuting unit phase.
    U=sp.Matrix(R49);V=sp.Matrix(R100)
    need(nu3(U.det())==nu3(V.det())==0,'each normalized raiser is a 3-unit')
    def modmat(p):return tuple(int(x)%3 for x in p)
    ident=modmat(sp.eye(2));group={ident};queue=[ident]
    for g in queue:
        for gen in (U,V):
            z=modmat(sp.Matrix(2,2,g)*gen)
            if z not in group:group.add(z);queue.append(z)
    need(len(group)==48,'full GL2(F3) normalized mixed action')
    mixed=[]
    for w in (('s78','s78')+D,D+('s78','s78')):
        vactual=vector(SEEDS[0]+w);rows=actual_target(SEEDS[0]+w)
        need(literal_frontier(rows)[0]==vactual[0],'actual mixed-raiser control')
        beta_actual=bulk(vactual);need(nu3(gcd(*beta_actual))==2,'two exact level raises')
        mixed.append({'word':w,'state':vactual,'vertices':len(rows),'port':len(rows)-1,
            'normalized_phase_mod3':[x//9%3 for x in beta_actual]})
    need(mixed[0]['normalized_phase_mod3']!=mixed[1]['normalized_phase_mod3'],'same gain, level and cost do not determine the unit phase')
    z=beta
    for _ in range(9):z=inverse_row(z,qB)
    witnesses=[]
    for name,R,word in [('49',R49,('s78','s78')),('100',R100,D)]:
        pre=inverse_row(z,R);need(all(t.denominator==1 for t in pre) and min(pre)<0,'integral gate escape still outside positive cone')
        full=word+('b',)*9;p=sp.Matrix(matrix(full));rows=[]
        for j in (1,2):
            u=sp.Matrix([[3*v[i]+j*T[i] for i in range(3)]])*p.inv()
            rows.append({'digit':j,'full_predecessor':[str(a) for a in u]})
        witnesses.append({'gate':name,'return_word':full,'cost':sum(COSTS[l] for l in full),
            'bulk_predecessor':[int(t) for t in pre],'full_predecessors':rows})
    # Equality case: strict growth in the generic cycle theorem is essential.
    equality=sp.Matrix([[2,1],[1,2]]);need(sp.Matrix([[1,1]])*equality==sp.Matrix([[3,3]]) and equality.det()==3,'equality-cycle hostile')
    result={'status':'PASS','catalog_size':len(records),'distinct_matrices':len({tuple(map(tuple,r['matrix'])) for r in records}),
        'growth_certificates':cert,'projective_union':[[str(a),str(b)] for a,b in merged],
        'seed_gaps':gaps,'clocks':clocks,'joint_gate':{'pair_period':1680,'acceptance_period':840,'hits':joint,'witnesses':witnesses},
        'mixed_raisers':{'normalized_action_group_size':48,'group_matrices':sorted(group),'controls':mixed},
        'equality_hostile':[[2,1],[1,2]],'scope':'All compositions of the finite complete short-return catalog plus s78; not all arbitrary longer 65-atom returns.'}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':'PASS','catalog':len(records),'intervals':len(merged),'all_seven_excluded':True,'joint_clock':[1680,840,joint],'minimum_growth_margin':min(F(c['margin']) for c in cert)},default=str))
if __name__=='__main__':main()
