"""Reuse one exact defect DAG for tournaments and self-converse tournaments.

A_n(2^beta,2^gamma)=sum_(odd lambda) 2^(beta*e+gamma*length)/z.
Only the scalar edge costs change; no further active-state quotient is claimed.
"""
import argparse,importlib.util,json,sys
from fractions import Fraction
from math import gcd,factorial
from itertools import combinations
from collections import Counter
from pathlib import Path
from time import perf_counter

HERE=Path(__file__).resolve().parent
ENGINE=HERE.parent/'session64_lca_reverse_valuation_20261026/session64-counting'
sys.path.insert(0,str(ENGINE))
import contracted_ordinary as engine
base=engine.base
if hasattr(sys,'set_int_max_str_digits'): sys.set_int_max_str_digits(0)

def evaluate(n):
    started=perf_counter(); dag=engine.build_contracted_dag(n)
    build=perf_counter()-started
    channels=((1,0),(2,0),(2,1))
    values=[[0]*len(dag.vertices) for _ in channels]
    for row in values: row[0]=1
    odd=base.odd_factorials(n); factors=0
    for source,outgoing in enumerate(dag.edges):
        for edge in outgoing:
            T=base.edge_factor(edge,odd); factors+=1
            for (beta,gamma),row in zip(channels,values):
                J=beta*edge.shift+(beta-1)*base.vfact(edge.k)+gamma*edge.k
                base.need(J>=0,'nonnegative deformed transfer')
                row[edge.target]+=(row[source]*T)<<J
        # Each source is fully consumed and never queried again.
        if source+1<len(dag.vertices):
            for row in values: row[source]=0
    answers=[]
    for row in values:
        a,rem=divmod(row[-1],odd[-1])
        base.need(rem==0,'exact terminal odd division')
        answers.append(a)
    return answers,{'vertices':len(dag.vertices),'edges':factors,'shared_odd_factor_evaluations':factors,'separate_odd_factor_evaluations':3*factors,'build_seconds':round(build,3),'total_seconds':round(perf_counter()-started,3)}

def odd_parts(n,low=1):
    if not n: yield ()
    for d in range(low,n+1,2):
        for tail in odd_parts(n-d,d): yield (d,)+tail

def direct(n,beta,gamma):
    total=Fraction()
    for parts in odd_parts(n):
        e=sum((d-1)//2 for d in parts)+sum(gcd(a,b) for a,b in combinations(parts,2))
        z=1
        for d,m in Counter(parts).items(): z*=d**m*factorial(m)
        total+=Fraction(1<<(beta*e+gamma*len(parts)),z)
    base.need(total.denominator==1,'direct integrality')
    return int(total)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--large',type=int,default=150); args=ap.parse_args()
    if not 30<=args.large<=200: raise ValueError('explicit bounded large order30..200')
    reference=json.loads((HERE.parent/'tournament_future_boundary_20261020/receipt.json').read_text(encoding='utf-8'))['large']['coefficients']
    rows=[]
    for n in list(range(31))+sorted(set((50,100,args.large))-{*range(31)}):
        answers,work=evaluate(n)
        base.need(answers[0]==int(reference[n]),'independent future-response A000568')
        if n<=30:
            base.need(answers==[direct(n,b,g) for b,g in ((1,0),(2,0),(2,1))],'three independent partition sums')
        rows.append({'m':n,'A000568_m':str(answers[0]),'self_converse_2m':str(answers[1]),'self_converse_2m_plus1':str(answers[2]),'work':work})
        if n in (30,50,100,args.large): print(json.dumps({'m':n,**work}),flush=True)
    native=json.loads((HERE.parent/'session67_metagraph_enumeration_information_20261026/session67-root/output/census.json').read_text(encoding='utf-8'))
    for row in native:
        n=row['n']; actual=int(rows[n//2]['self_converse_2m_plus1' if n%2 else 'self_converse_2m'])
        base.need(actual==row['self_converse'],'native metagraph self-converse census')
    report={'status':'FINITE-EXACT','channels':[[1,0],[2,0],[2,1]],'direct_partition_control_through_m':30,'native_self_converse_control_through_n':8,'independent_A000568_control':'future-response exact sequence through200','rows':rows}
    (HERE/'sequence_channels.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('PASS: common structural carrier, three exact counting channels.')

if __name__=='__main__':main()
