"""Independent literal checks of the 32-fraction apex and 20-fraction core programs."""
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'fractran'))
sys.path.insert(0,str(Path(__file__).parent))
from audit_path_counter import count_dp

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def raw(rows,N,limit=2_000_000):
    steps=0
    while steps<limit:
        for a,b in rows:
            if N%b==0:
                N=N//b*a;steps+=1;break
        else:return N,steps
    raise RuntimeError('raw limit')

def valuation(N,p):
    k=0
    while N%p==0:N//=p;k+=1
    return k,N

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    apex=load('apex',ROOT/'path_counter/one_apex.py')
    emitted=[tuple(map(int,line.split('/'))) for line in (ROOT/'path_counter/one_apex.fractions.txt').read_text().splitlines()]
    need(tuple(emitted)==apex.FRACTIONS,'apex emitted fractions differ')
    apex_cases=0;maxsteps=0;outputs=set()
    for m in range(9):
        for word in itertools.product((0,1),repeat=m):
            # Vertex0 is apex, vertexj+1 is the chain vertex indexed by j.
            bits=0
            for j in range(1,m+1):
                for i in range(j):
                    edge=word[j-1] if i==0 else 1
                    bits|=edge<<(j*(j-1)//2+i)
            H=count_dp(m+1,bits)
            cycles=0
            for subset in range(1<<m):
                vertices=[j for j in range(m) if subset>>j&1]
                if len(vertices)>=2 and len(vertices)%2==0:
                    cycles+=word[vertices[0]]==1 and word[vertices[-1]]==0
            E=(1<<m)+sum(b<<j for j,b in enumerate(word))
            N,steps=raw(emitted,38*5**E)
            Q=0
            for b in word:Q=2*Q+b
            need(steps==4*E+3*Q-6*Q.bit_count()-(Q%2)+4*m-2,
                 ('apex exact runtime',word))
            actual_H,N=valuation(N,2);actual_cycles,N=valuation(N,3)
            _,N=valuation(N,11);_,N=valuation(N,13)
            need(N==107,('apex terminal',word,N))
            need((actual_H,actual_cycles)==(H,cycles),('apex count',word,H,cycles))
            need(H==1+2*cycles,('independent OCF check',word))
            maxsteps=max(maxsteps,steps);apex_cases+=1;outputs.add(H)
    N,_=raw(emitted,38)
    need(N==2*109,'zero encoding must error')
    core=load('core',ROOT/'fractran/core_counter.py')
    # Independent min-product over all SCC-size compositions, using Busch f.
    dp=[None]*151;dp[0]=1
    def strong(n):
        q,r=divmod(n,3)
        return [3*5**(q-1),5**q,9*5**(q-1)][r]
    for n in range(3,151):
        dp[n]=min(strong(k)*dp[n-k] for k in range(3,n+1) if dp[n-k] is not None)
        need(core.expected(n)==dp[n],('core min-product',n))
    core_cases=[]
    for tail in [False,True]:
        for n in range(6 if tail else 3,22):
            initial=(11 if tail else 43)*5**(n-(6 if tail else 3))
            N,steps=raw(core.TAIL14 if tail else core.FULL20,initial)
            output,N=valuation(N,2)
            need(N==41 and output==dp[n],('core literal',tail,n))
            need(steps==core.exact_steps(n,tail),('core exact cost',tail,n))
            core_cases.append({'tail':tail,'n':n,'output':output,'steps':steps})
    print(json.dumps({'status':'PASS','apex_all_words_max_length':8,
          'apex_literal_cases':apex_cases,'apex_max_steps':maxsteps,
          'apex_attained_H':sorted(outputs),'core_partition_DP_through':150,
          'core_literal_cases':core_cases,
          'scope':'All-order program semantics are established by loop invariants; finite checks independently compare actual raw fractions to tournament DP/cycle enumeration and component-product DP.'},indent=2))

if __name__=='__main__':main()
