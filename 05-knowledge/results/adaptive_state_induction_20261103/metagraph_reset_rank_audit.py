#!/usr/bin/env python3
"""Independent actual-graph audit of the concurrent cap-two query minors.
Reads the producer's certificate as JSON; imports none of its implementation.
Counts Hamilton words by recursive completion on the actual named graph.
"""
from pathlib import Path
from functools import lru_cache
from fractions import Fraction
import subprocess,json,hashlib,sys


def need(p,m):
    if not p:raise RuntimeError(m)


def birth(T,r):
    n=len(T)
    return tuple(row|((1<<n) if i==r else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<r),)


def active(T):return tuple(v for v,row in enumerate(T) if len(T)-1-row.bit_count()<=2)


@lru_cache(None)
def graph(word):
    T=birth((2,4,1),0)
    for letter in word:
        A=active(T);need(0<=letter<len(A),'actual active-rank guard');T=birth(T,A[letter])
    return T


@lru_cache(None)
def H(T):
    full=(1<<len(T))-1
    @lru_cache(None)
    def finish(v,remaining):
        if not remaining:return 1
        total=0;choices=T[v]&remaining
        while choices:
            b=choices&-choices;choices-=b;total+=finish(b.bit_length()-1,remaining^b)
        return total
    return sum(finish(v,full^(1<<v)) for v in range(len(T)))


def query(word):return H(graph(tuple(word)))


def guarded_query(word):
    T=birth((2,4,1),0)
    for a in word:
        A=active(T)
        if not 0<=a<len(A):return 0
        T=birth(T,A[a])
    return H(T)


def det(A):
    B=[r[:] for r in A];n=len(B);sgn=1;prev=1
    if n==1:return B[0][0]
    for k in range(n-1):
        i=next((i for i in range(k,n) if B[i][k]),None)
        if i is None:return 0
        if i!=k:B[i],B[k]=B[k],B[i];sgn=-sgn
        pivot=B[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=B[i][j]*pivot-B[i][k]*B[k][j]
                need(numerator%prev==0,'fraction-free determinant exact division')
                B[i][j]=numerator//prev
            B[i][k]=0
        prev=pivot
    return sgn*B[-1][-1]


def main():
    source=Path(__file__).with_name('active_rank_exact.py')
    receipt=json.loads(subprocess.check_output([sys.executable,'-I','-B',str(source)]))
    charts={tuple(c['control']):c for c in receipt['charts']};results=[];entries=0;transfer_queries=0
    for control,c in charts.items():
        for prefix in c['reachable_prefixes']:
            T=graph(tuple(prefix));actual=tuple(len(T)-1-T[v].bit_count() for v in active(T))
            need(actual==control,'actual reachable-prefix control')
        prefixes=c['hankel_prefixes'];suffixes=c['query_suffixes']
        actual=[[query(tuple(p)+tuple(s)) for s in suffixes] for p in prefixes]
        need(actual==c['hankel_matrix'],'all actual Hamilton-query minor entries')
        determinant=det(actual);need(determinant==Fraction(c['hankel_determinant']) and determinant!=0,'independent exact Hankel determinant')
        entries+=len(prefixes)*len(suffixes)
        results.append({'control':control,'dimension':len(suffixes),'determinant':determinant,
                        'largest_actual_order':max(len(graph(tuple(p)+tuple(s))) for p in prefixes for s in suffixes)})
    # A genuine actual prefix in each chart, then every available root and
    # every target query. This is a named-graph check of the rational transport.
    for t in receipt['transfers']:
        c=charts[tuple(t['source'])];dest=charts[tuple(t['target'])]
        p=tuple(c['hankel_prefixes'][0]);a=t['root_rank']
        y=[query(p+tuple(s)) for s in c['query_suffixes']]
        M=[[Fraction(x) for x in row] for row in t['matrix']]
        actual=[query(p+(a,)+tuple(s)) for s in dest['query_suffixes']]
        predicted=[sum(y[i]*M[i][j] for i in range(len(y))) for j in range(len(actual))]
        need(predicted==actual,'actual named query transport')
        transfer_queries+=len(actual)
    stable=receipt['stable_return_control'];prefix=tuple(stable['actual_seed_active_rank_prefix'])
    initial_tail=[query(prefix+(1,)*j) for j in range(5)]
    need(initial_tail==stable['root1_tail_first12'][:5],'actual stable return scalar prefix')
    need(any(not p for p in charts[tuple(receipt['initial_control'])]['reachable_prefixes']),'seed included in exact invariant span')
    c=charts[(2,2,2,1)];y=[query((1,)+tuple(s)) for s in c['query_suffixes']]
    need(y==[13,23,29,25,29,49,43],'actual mod75 source queries')
    coefficients=[-236,0,299,-40,-12,-29,37]
    value=Fraction(sum(a*b for a,b in zip(coefficients,y)),75)
    fake=y[:];fake[0]=15;fake_value=Fraction(sum(a*b for a,b in zip(coefficients,fake)),75)
    need(value==query((1,0,0))==59 and fake_value==Fraction(3953,75),'literal lattice guard')
    global_data=receipt['global_guarded_hankel']
    global_matrix=[[guarded_query(tuple(p)+tuple(s)) for s in global_data['suffixes']] for p in global_data['prefixes']]
    need(global_matrix==global_data['matrix'],'all guarded global actual query entries')
    global_det=det(global_matrix)
    need(global_det==global_data['determinant']==12470867040665600,'independent global guarded determinant')
    return {'status':'PASS independent actual Hamilton-query minors and transports',
            'producer_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'charts':results,
            'literal_query_entries':entries,'actual_transition_queries':transfer_queries,
            'actual_return_prefix':initial_tail,'distinct_actual_graphs_counted':H.cache_info().currsize,
            'mod75_guard':{'actual_queries':y,'actual_H00':int(value),'fake_H15_H00':str(fake_value)},
            'global_guarded_minor':{'dimension':len(global_matrix),'entries':sum(map(len,global_matrix)),
                                    'illegal_zero_entries':sum(v==0 for row in global_matrix for v in row),'determinant':global_det},
            'scope':'All198 local and784 guarded global lower-bound entries,19 transitions on one actual prefix per source chart, and actual scalar prefix; full rational closure proof separately code-audited, not reimplemented here. Global minimum is scoped to illegal words assigned zero.'}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
