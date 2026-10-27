#!/usr/bin/env python3
"""Exact 15-move certificate check using adjacency matrices, not finder routines."""
import hashlib,json
from pathlib import Path
from itertools import combinations,permutations

def need(x,msg):
    if not x:raise RuntimeError(msg)

def scores(A):return [sum(row) for row in A]
def L(A):return sum(abs(v-i) for i,v in enumerate(sorted(scores(A))))//2

def native(A,p):
    n=len(A);need(sorted(p)==list(range(n)),'nonpermutation');need(all(A[a][b] for a,b in zip(p,p[1:])),'illegal current path')
    keep={frozenset((a,b)) for a,b in zip(p,p[1:])};B=[[False]*n for _ in range(n)]
    for a,b in combinations(range(n),2):
        ab=A[a][b] if frozenset((a,b)) in keep else not A[a][b]
        B[a][b]=ab;B[b][a]=not ab
    return B

def fas(A):
    n=len(A);inc=[sum(1<<u for u in range(n) if A[u][v]) for v in range(n)];dp=[0]*(1<<n)
    for S in range(1,1<<n):
        todo=S;best=n*n
        while todo:
            bit=todo&-todo;todo-=bit;v=bit.bit_length()-1;R=S^bit;best=min(best,dp[R]+(inc[v]&R).bit_count())
        dp[S]=best
    return dp[-1]

def main():
    p=Path(__file__).resolve().parent;c=json.loads((p/'certificate15.json').read_text());n=c['n'];qr={x*x%n for x in range(1,n)}
    A=[[(b-a)%n in qr for b in range(n)] for a in range(n)];need(L(A)==15,'initial score');table=[]
    for k,s in enumerate(c['steps'],1):
        need(s['depth']==k,'step numbering');before=A;A=native(A,s['path']);rows=[sum(1<<b for b in range(n) if A[a][b]) for a in range(n)]
        need(rows==s['rows'],'row mismatch');need(native(A,s['path'])==before,'inverse failure');need(L(A)==15-k,'not score-geodesic')
        table.append({'depth':k,'L':L(A),'f':fas(A),'scores':scores(A)})
    need(len(table)==15,'wrong length');order=c['final_forward_order'];need(sorted(order)==list(range(n)),'final order invalid')
    need(all(A[order[i]][order[j]] for i in range(n) for j in range(i+1,n)),'final tournament not transitive in claimed order')
    need(sorted(scores(A))==list(range(n)) and fas(A)==0,'final control')
    # Independent finite hostile check for the no-strict-descent source/sink lemma.
    controls=[]
    for nn in range(2,6):
        es=list(combinations(range(nn),2));pairs=0
        for mask in range(1<<len(es)):
            B=[[False]*nn for _ in range(nn)]
            for j,(a,b) in enumerate(es):B[a][b]=bool(mask>>j&1);B[b][a]=not B[a][b]
            if not any(v in (0,nn-1) for v in scores(B)):continue
            for path in permutations(range(nn)):
                if all(B[a][b] for a,b in zip(path,path[1:])):
                    pairs+=1;need(L(native(B,path))>=L(B),'source/sink descent hostile')
        controls.append({'n':nn,'all_legal_moves_from_source_or_sink_states':pairs})
    stalled=[4,2037,512,1015,709,725,517,1093,757,129,893]
    B=[[bool(stalled[a]>>b&1) for b in range(n)] for a in range(n)];need(L(B)==1 and fas(B)==4,'stalled control')
    legal_count=0;descents=0
    def visit(path,unused):
        nonlocal legal_count,descents
        if not unused:
            legal_count+=1;descents+=L(native(B,path))==0;return
        for v in sorted(unused):
            if not path or B[path[-1]][v]:visit(path+[v],unused-{v})
    visit([],set(range(n)));need(legal_count==243 and descents==0,'stalled exact scan mismatch')
    # Search files are reproducible but not imported or trusted as theorem evidence.
    digest=hashlib.sha256((p.parent/'certificate.json').read_bytes()).hexdigest()
    print(json.dumps({'status':'PASS','initial':'P11','final_forward_order':order,'moves':table,'exact_distance_to_transitive_class':15,'distance_justification':'15 explicit legal moves plus inherited |delta L|<=1 and L(P11)=15; not exhaustive distance search','parent_six_move_certificate_sha256':digest,'source_sink_hostile_checks':controls,'stalled_greedy_control':{'L':1,'f':4,'legal_paths':legal_count,'one_step_transitive_descents':descents},'caveat':'No distance to every fixed labelled transitive tournament is claimed.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
