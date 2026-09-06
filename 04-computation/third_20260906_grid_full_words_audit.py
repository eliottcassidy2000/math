#!/usr/bin/env python3
"""Independent full-word E maxima: numerical order, scalar/singleton relaxation,
bounded-knapsack pruning, and literal final words. No producer imports.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,combinations_with_replacement
from math import gcd,lcm
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
PROFILE_SHA='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
BASELINE_SHA='a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58'
MAXIMA_SHA='ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca'
GATES=0
TRACE=sha256()

def canonical(x):return json.dumps(x,separators=(',',':'))
def digest(x):return sha256(canonical(x).encode()).hexdigest()
def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data!r}')
    GATES+=1;TRACE.update((label+':'+repr(data)+'\n').encode())

def load():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    gate(sha256(raw).hexdigest()==PROFILE_SHA,'pinned-complete-word-supplier')
    P={int(k):{(c,tuple(w)) for c,w in L['profiles']} for k,L in json.loads(raw)['levels'].items()}
    baseline=[]
    for line in (ROOT/'05-knowledge/results/third_20260906_grid_bootstrap.out').read_text().splitlines():
        if line.startswith('SCALE_SET '):
            row=json.loads(line.split(' ',1)[1])
            if row['name']=='profile6_coupled':baseline=row['survivors']
    gate(len(baseline)==8301 and max(baseline)==23760 and digest(baseline)==BASELINE_SHA,'exact-audited-clock-universe')
    rows=[];global_rows=[];boundaries=[]
    for line in (ROOT/'05-knowledge/results/third_20260906_grid_full_words.out').read_text().splitlines():
        if line.startswith('MAXIMUM '):rows.append(json.loads(line.split(' ',1)[1]))
        if line.startswith('GLOBAL '):global_rows.append(json.loads(line.split(' ',1)[1]))
        if line.startswith('BOUNDARY '):boundaries.append(json.loads(line.split(' ',1)[1]))
    gate([row[0] for row in rows]==baseline,'every-clock-candidate-record')
    gate(digest([r[:2] for r in rows])==MAXIMA_SHA,'pinned-claimed-maxima-table')
    gate([r[0] for r in global_rows]==list(range(1,7)),'six-global-residue-candidates')
    return P,rows,global_rows,boundaries

class IndependentOptimizer:
    def __init__(self,P):
        self.P=P;self.cores={k:{c for c,w in rows} for k,rows in P.items()}
        self.values=tuple(sorted(self.cores[6]));self.subwords={}
        for d,w in P[6]:
            for n in range(7):
                self.subwords.setdefault((d,n),set()).update(tuple(w[i] for i in I) for I in combinations(range(6),n))
        self.caps={d:(7 if d==1 else max(k for k in range(1,7) if d in self.cores[7-k])) for d in self.values}
        gate(len(self.values)==42 and all(gcd(d,7)==1 for d in self.values),'all-allowed-core-values-coprime-seven')
    @lru_cache(maxsize=350000)
    def relaxed_prefix(self,S):
        n=len(S)
        # Unlike the producer, only singleton full-word projections are used
        # here; larger subsets retain just their scalar allowed gcd sets.
        for i,d in enumerate(S):
            w=tuple(sorted(gcd(d,v) for j,v in enumerate(S) if i!=j))
            if w not in self.subwords[d,n-1]:return False
        for k in range(2,min(6,n)+1):
            if any(gcd(*I) not in self.cores[7-k] for I in combinations(S,k)):return False
        return n<7 or gcd(*S)==1
    def full(self,S,check=False):
        if len(S)!=7 or gcd(*S)!=1:return False
        for k in range(1,7):
            for I in combinations(range(7),k):
                d=gcd(*(S[i] for i in I));w=tuple(sorted(gcd(d,S[j]) for j in range(7) if j not in I))
                ok=(d,w) in self.P[7-k]
                if check:gate(ok,'literal-candidate-full-word',(S,I))
                if not ok:return False
        return True
    def score(self,t,S):return sum(d*((-(t//d))%7) for d in S)
    def candidate(self,t,E,owner):
        S=tuple(owner)
        gate(all(t%d==0 for d in S) and self.full(S,True),'literal-attaining-candidate',(t,E,S))
        direct=sum(d*((t+7*d-1)//(7*d)) for d in S)-t
        gate(self.score(t,S)==7*direct==7*E,'candidate-cost-and-direct-ceiling-identity',(t,E))
        return S
    def solve(self,t,owner=None,forced=()):
        ds=[d for d in self.values if t%d==0];ws=[d*((-(t//d))%7) for d in ds]
        fixed=Counter(forced)
        gate(all(d in ds and n<=self.caps[d] for d,n in fixed.items()),'conditioned-state-capacity',forced)
        capacities=[self.caps[d]-fixed[d] for d in ds]
        n=len(ds);dp=[[-10**9]*8 for _ in range(n+1)];dp[n][0]=0
        # Exact upper bound for the relaxed suffix multiset knapsack.
        for j in range(n-1,-1,-1):
            for r in range(8):
                dp[j][r]=max(k*ws[j]+dp[j+1][r-k] for k in range(min(capacities[j],r)+1))
        best=-1;best_owner=None
        if owner is not None:
            gate(self.full(owner) and all(t%d==0 for d in owner) and not (fixed-Counter(owner)),'verified-incumbent-lower-bound',(t,forced))
            best=self.score(t,owner);best_owner=tuple(owner)
        nodes=leaves=0
        def visit(selected,start,score):
            nonlocal best,best_owner,nodes,leaves
            nodes+=1;r=7-len(forced)-len(selected)
            if not r:
                leaves+=1;S=tuple(sorted(tuple(forced)+selected))
                if self.full(S) and score>best:best,best_owner=score,S
                return
            for j in range(start,n):
                d=ds[j];available=capacities[j]-selected.count(d)
                if available<=0:continue
                bound=max(k*ws[j]+dp[j+1][r-k] for k in range(1,min(available,r)+1))
                if score+bound<=best:continue
                extension=selected+(d,);S=tuple(sorted(tuple(forced)+extension))
                if self.relaxed_prefix(S):visit(extension,j,score+ws[j])
        gate(not forced or self.relaxed_prefix(tuple(sorted(forced))),'conditioned-prefix-feasible',forced)
        visit((),0,self.score(t,forced))
        gate(best>=0 and best%7==0 and self.full(best_owner),'independent-attained-optimum',(t,forced))
        return best//7,best_owner,nodes,leaves

def literal_component_credit(t,e,p,q):
    cuts={F(0),F(1)}
    for speed in (p,q):
        for j in range(speed+1):
            for wall in (F(14*j-1,14*speed),F(14*j+1,14*speed)):
                if 0<wall<1:cuts.add(wall)
    cuts=sorted(cuts)
    def bad(x,speed):
        z=(x*speed)%1
        return min(z,1-z)<F(1,14)
    flags=[bad((a+b)/2,p) and bad((a+b)/2,q) for a,b in zip(cuts,cuts[1:])]
    lengths=[b-a for a,b,yes in zip(cuts,cuts[1:],flags) if yes]
    if flags[0] and flags[-1]:lengths=[lengths[0]+lengths[-1]]+lengths[1:-1]
    n=t//e
    ceil=lambda x:-((-x.numerator)//x.denominator)
    return e*sum(ceil(n*ell)-1 for ell in lengths),len(lengths)

def main():
    P,rows,global_rows,boundaries=load();opt=IndependentOptimizer(P)
    patterns={}
    for t,E,owner,_ in rows:
        opt.candidate(t,E,owner)
        key=(tuple(d for d in opt.values if t%d==0),t%7)
        if key in patterns:gate(patterns[key][1]==E,'same-divisor-residue-same-candidate-value',(t,E))
        else:patterns[key]=(t,E,tuple(owner))
    gate(len(patterns)==1440,'complete-divisor-residue-pattern-quotient')
    total_nodes=total_leaves=0;pattern_records=[]
    for key,(t,E,owner) in sorted(patterns.items(),key=lambda item:item[1][0]):
        actual,witness,nodes,leaves=opt.solve(t,owner)
        gate(actual==E,'independent-optimum-for-entire-clock-pattern',(t,E))
        total_nodes+=nodes;total_leaves+=leaves;pattern_records.append([t,E,nodes,leaves])
        print('AUDITED_PATTERN',canonical(pattern_records[-1]))
    gate(lcm(*opt.values)==5388292800,'all-value-lcm')
    global_results=[]
    for tau,t,E,owner,_ in global_rows:
        gate(t%7==tau and all(t%d==0 for d in opt.values),'global-residue-realizes-full-value-domain',(tau,t))
        opt.candidate(t,E,owner)
        actual,witness,nodes,leaves=opt.solve(t,tuple(owner))
        gate(actual==E==[210,270,192,239,197,224][tau-1],'independent-global-word-optimum',(tau,E))
        global_results.append(E)
        print('AUDITED_GLOBAL',canonical([tau,t,E,witness,nodes,leaves]))
    # Complete multiset enumeration is deliberately unpruned and independent
    # of either optimizer for this finite small-clock control bank.
    brute_count=0
    for t in range(1,33):
        allowed=[d for d in opt.values if t%d==0];best=-1;owner=None
        for S in combinations_with_replacement(allowed,7):
            brute_count+=1
            if opt.full(S):
                score=opt.score(t,S)
                if score>best:best,owner=score,S
        actual,_,_,_=opt.solve(t,(1,)*7)
        gate(actual*7==best,'full-unpruned-small-clock-enumeration',t)
    for t,E,owner in ((23760,168,(72,44,30,16,16,9,9)),(27360,164,(1,9,16,18,32,60,72)),(28080,112,(1,4,9,24,36,40,54))):
        opt.candidate(t,E,owner);actual,_,_,_=opt.solve(t,owner)
        gate(actual==E,'independent-large-control-optimum',t)
    for t in (7,49,343,7*5388292800):
        E,_,_,_=opt.solve(t,(1,)*7)
        gate(E==0,'septimal-zero-cost-control',t)
    controls=[((6,12),133,(6,12,58,64,36,9,1),6,22,205,174),((12,16),188,(12,16,72,58,64,9,9),4,3,308,172)]
    for margins,E,owner,e,p,q,credit in controls:
        opt.candidate(16704,E,owner)
        actual,witness,nodes,leaves=opt.solve(16704,owner,margins)
        gate(actual==E,'conditioned-full-word-optimum',margins)
        gate(tuple(sorted((e*gcd(16704//e,p),e*gcd(16704//e,q))))==margins,'literal-forced-pair-sheet-marginals',margins)
        literal,components=literal_component_credit(16704,e,p,q)
        gate(literal==credit,'independent-literal-component-credit',(e,p,q,credit))
        gate((credit>E)==(margins==(6,12)),'first-closure-versus-stopping-boundary',margins)
        print('AUDITED_CONDITIONED',canonical([16704,list(margins),E,witness,e,p,q,credit,components,nodes,leaves]))
    gate(any(row[:7]==[16704,4,3,308,[12,16],188,172] for row in boundaries),'frozen-positive-stopping-boundary')
    print('ALL_CLOCKS',len(rows),'PATTERNS',len(patterns),'INDEPENDENT_NODES',total_nodes,'RELAXED_LEAVES',total_leaves,'UNPRUNED_SMALL_MULTISET_CASES',brute_count)
    print('AUDITED_MAXIMA_SHA256',MAXIMA_SHA,'PATTERN_TRACE_SHA256',digest(pattern_records))
    print('SCOPE exact maxima in full seven-state profile quotient; physical realizability and strict failure are not inferred')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())

if __name__=='__main__':main()
