"""Big-integer implementation cross-check and short independent literal controls."""
import json,gzip
from functools import lru_cache
from math import gcd
from repair23 import literal_profile,join,birth
from context_transport import glue
from fibre_probe import HERE,save,need

def load(name):return json.loads(gzip.decompress((HERE/(name+'.json.gz')).read_bytes()))

def reference(T,r):
    n=len(T);full=(1<<n)-1
    @lru_cache(None)
    def ending(S,v):
        R=S^(1<<v)
        if not R:return 1
        return sum(ending(R,u) for u in range(n) if R>>u&1 and T[u]>>v&1)
    @lru_cache(None)
    def h(S):return sum(ending(S,v) for v in range(n) if S>>v&1) if S else 1
    @lru_cache(None)
    def starting(S,v):
        R=S^(1<<v)
        if not R:return 1
        return sum(starting(R,u) for u in range(n) if R>>u&1 and T[v]>>u&1)
    return h(full),starting(full,r),h(full^(1<<r)),sum(ending(S,r)*h(full^S) for S in range(full+1) if S>>r&1)

def run():
    selector=load('prefix_selector');selected=0;maxn=0
    for rows in selector['outputs'].values():
        for row in rows:
            T=tuple(row['target_rows']);need(reference(T,row['root'])==tuple(row['profile']),'Python integer target recount');selected+=1;maxn=max(maxn,len(T))
    words=load('transport_words');wordgraphs=0;literal=0
    for obj in words['outputs']:
        for T,z in zip(obj['output_rows'],obj['output_profiles']):
            need(reference(tuple(T),obj['root'])==tuple(z),'Python integer transport recount');wordgraphs+=1
            if len(T)<=7:need(literal_profile(T,obj['root'])==tuple(z),'permutation transport recount');literal+=1
    # The H9 one-arc generator has a physical slope conversion, with a hostile gcd.
    A=(6,4,0);root=1;pair=((22,8,2,5,14),(22,8,2,21,6));states=[]
    for B,expected in zip(pair,((15,0,3,30),(15,0,5,30))):
        T=glue(A,root,B,0);z=reference(T,root);need(z==expected==literal_profile(T,root),'H9 slope-half conversion')
        h,s,q,F=z;need(gcd(h-s+3*q,F-q)==q,'conversion does not repair chain gcd')
        states.append(dict(rows=T,root=root,profile=z))
    # Exact integer inverse of the formal fibre lattice; includes signed coordinates.
    lattice=0
    for h in (1,9,23):
        for M in (7,21,81,143):
            for c in range(-20,21):
                for d in range(-20,21):
                    if (c+2*d-M-2)%4:continue
                    q=(c+2*d-M)//2;s=(c+4*d+2*h-3*M)//2;F=s+M-h
                    need(q%2==1 and (h-s+4*q-F,F-q)==(c,d),'parity lattice inverse');lattice+=1
    save('audit',dict(status='PASS',independent_agent_review=False,target_recounts=selected,max_target_order=maxn,transport_bigint_recounts=wordgraphs,transport_permutation_recounts=literal,H9_conversion=states,signed_lattice_checks=lattice,limitations='Self-audit, with independent literal/permutation readers and Python integers; normal/-O replay is a regression check, not independent mathematical evidence.'))
    print('PASS',selected,'targets;',wordgraphs,'transport graphs;',literal,'permutation controls;',lattice,'signed lattice checks')

if __name__=='__main__':run()
