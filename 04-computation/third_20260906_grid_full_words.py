#!/usr/bin/env python3
"""Exact full-complement-word marginal costs on the audited finite clock set."""
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
PROFILE_SHA='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
BASELINE_SHA='a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58'
GATES=0

def need(ok,why):
    global GATES
    GATES+=1
    if not ok: raise ArithmeticError(why)

def canonical(value):
    return json.dumps(value,separators=(',',':'))

def load():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    need(sha256(raw).hexdigest()==PROFILE_SHA,'pinned complete inherited profile table')
    P={int(k):{(c,tuple(w)) for c,w in level['profiles']} for k,level in json.loads(raw)['levels'].items()}
    baseline=[]
    for line in (ROOT/'05-knowledge/results/third_20260906_grid_bootstrap.out').read_text().splitlines():
        if line.startswith('SCALE_SET '):
            obj=json.loads(line[10:])
            if obj['name']=='profile6_coupled': baseline=obj['survivors']
    need(sha256(canonical(baseline).encode()).hexdigest()==BASELINE_SHA,'pinned audited baseline clock set')
    need(len(baseline)==8301 and max(baseline)==23760,'complete declared clock universe')
    sub={}
    for k,rows in P.items():
        for c,w in rows:
            for n in range(k+1):
                sub.setdefault((k,c,n),set()).update(tuple(w[i] for i in I) for I in combinations(range(k),n))
    return P,baseline,sub

class Optimizer:
    def __init__(self,P,sub):
        self.P,self.sub=P,sub
        self.values=sorted({c for c,w in P[6]})
    def valid_prefix(self,S):
        n=len(S)
        for k in range(1,min(n,6)+1):
            for I in combinations(range(n),k):
                c=gcd(*(S[i] for i in I))
                w=tuple(sorted(gcd(c,S[j]) for j in range(n) if j not in I))
                if w not in self.sub.get((7-k,c,n-k),()): return False
        return n<7 or gcd(*S)==1
    def literal_full(self,S,check=False):
        if len(S)!=7 or gcd(*S)!=1: return False
        for k in range(1,7):
            for I in combinations(range(7),k):
                c=gcd(*(S[i] for i in I))
                w=tuple(sorted(gcd(c,S[j]) for j in range(7) if j not in I))
                ok=(c,w) in self.P[7-k]
                if check: need(ok,'literal maximizing full complement word')
                if not ok: return False
        return True
    def solve(self,t,forced=()):
        values=sorted((d for d in self.values if t%d==0),key=lambda d:(-d*((-(t//d))%7),d))
        costs=[d*((-(t//d))%7) for d in values]
        best=-1;owner=None;visited=0
        def visit(S,start,score):
            nonlocal best,owner,visited
            visited+=1
            if len(S)==7:
                if score>best: best,owner=score,S
                return
            for j in range(start,len(values)):
                # Every completion here uses only values at indices >=j.
                # Their costs cannot exceed costs[j]; repetitions are allowed.
                if score+(7-len(S))*costs[j]<=best: break
                candidate=S+(values[j],)
                if self.valid_prefix(candidate): visit(candidate,j,score+costs[j])
        need(all(d in values for d in forced) and self.valid_prefix(forced),'admissible fixed pair prefix')
        visit(tuple(forced),0,sum(d*((-(t//d))%7) for d in forced))
        need(owner is not None and best%7==0,'attained integer seven-tail excess')
        need(all(t%d==0 for d in owner),'all maximizing sheet states divide their clock')
        need(self.literal_full(owner,True),'complete maximizing owner passes direct membership')
        direct=sum(d*((t+7*d-1)//(7*d)) for d in owner)-t
        need(direct==best//7,'literal marginal count equals cost arithmetic')
        return best//7,list(owner),visited


def main():
    P,baseline,sub=load();solver=Optimizer(P,sub)
    brute_cases=0
    for t in range(1,25):
        value,owner,_=solver.solve(t)
        allowed=[d for d in solver.values if t%d==0]
        exact=-1
        for S in combinations_with_replacement(allowed,7):
            brute_cases+=1
            if solver.literal_full(S):
                exact=max(exact,sum(d*((t+7*d-1)//(7*d)) for d in S)-t)
        need(value==exact,'small-clock complete unpruned multiset enumeration')
    for t,E in ((23760,168),(27360,164),(28080,112)):
        need(solver.solve(t)[0]==E,'large relaxed-endpoint hostile repaired by full words')
    need(solver.solve(7)[0]==0 and solver.solve(49)[0]==0,'zero-cost septimal controls')
    print('FINITE-EXACT: exact maximum E over full seven-state word profiles, not physical realizability.')
    L=lcm(*solver.values)
    need(L==5388292800 and all(gcd(d,7)==1 for d in solver.values),'full core projection excludes septimal sheet gcds')
    for tau,expected in enumerate((210,270,192,239,197,224),1):
        t=L*((tau*pow(L,-1,7))%7)
        E,owner,visited=solver.solve(t)
        need(t%7==tau and E==expected,'sharp full-profile residue envelope')
        print('GLOBAL',canonical([tau,t,E,owner,visited]))
    low,low_owner,visits=solver.solve(16704,(6,12))
    need(low==133,'first prospective boundary pair is excluded by full words')
    print('CONDITIONED',canonical([16704,[6,12],low,low_owner,visits]))
    high,high_owner,visits=solver.solve(16704,(12,16))
    need(high==188,'another boundary pair survives the full-word refinement')
    witness=(12,16,72,58,64,9,9)
    need(solver.literal_full(witness,True),'explicit compatible alternative boundary owner')
    need(sum(d*((16704+7*d-1)//(7*d)) for d in witness)-16704==188,'boundary excess attained')
    t,e,p,q=16704,4,3,308
    need((e*gcd(t//e,p),e*gcd(t//e,q))==(12,16),'exact pair marginals at boundary')
    need(p+q==311 and all(311%d for d in range(2,18)),'literal inert prime atlas sum')
    den=14*p*q
    numerators=[2*p]+[min(2*p,p+q-14*r) for r in range(1,(p+q-1)//14+1) for _ in range(2)]
    credit=e*sum(((t//e)*n+den-1)//den-1 for n in numerators)
    need(credit==172<188,'finite component credit genuinely leaves this relaxation open')
    print('BOUNDARY',canonical([t,e,p,q,[12,16],188,credit,list(witness)]))
    print('PROFILE_SHA256',PROFILE_SHA)
    print('BASELINE_SET_SHA256',BASELINE_SHA)
    summary=[];total_nodes=0
    for t in baseline:
        E,owner,visited=solver.solve(t)
        summary.append([t,E]);total_nodes+=visited
        print('MAXIMUM',canonical([t,E,owner,visited]))
    print('MAXIMA_SHA256',sha256(canonical(summary).encode()).hexdigest())
    print('CLOCKS',len(baseline),'PREFIX_NODES',total_nodes,'UNPRUNED_SMALL_MULTISETS',brute_cases)
    print('PASS gates='+str(GATES))

if __name__=='__main__': main()
