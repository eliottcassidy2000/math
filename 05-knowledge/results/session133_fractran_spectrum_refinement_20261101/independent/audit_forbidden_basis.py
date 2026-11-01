"""Independent full n6 orbit-cover and 56*64 order7 augmentation audit."""
from collections import Counter
from itertools import combinations,permutations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def h_all(a):
    n=len(a);end=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
    for v in range(n):end[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            value=end[mask][v];h[mask]+=value
            if not value:continue
            for w in range(n):
                if a[v]>>w&1 and not(mask>>w&1):end[mask|1<<w][w]+=value
    return h

def code(a,p=None):
    if p is None:p=tuple(range(len(a)))
    return sum(((a[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)))

def canon(a):return min(code(a,p) for p in permutations(range(len(a))))

def strong(a):
    n=len(a)
    if n<3:return False
    for v in range(n):
        reached=1<<v
        while True:
            new=reached
            for u in range(n):
                if reached>>u&1:new|=a[u]
            if new==reached:break
            reached=new
        if reached!=(1<<n)-1:return False
    return True

def main():
    atlas=json.loads((ROOT/'spectrum/atlas.json').read_text())
    basis=json.loads((ROOT/'spectrum/basis21.json').read_text())
    rows6=[r['rows'] for r in atlas['classes'] if r['n']==6]
    need(len(rows6)==56,'need56 reps')
    coverage=set()
    for a in rows6:
        orbit={code(a,p) for p in permutations(range(6))}
        need(not coverage.intersection(orbit),'duplicate n6 class')
        coverage.update(orbit)
    need(coverage==set(range(1<<15)),'incomplete n6 cover')
    expected6={canon(a) for a in rows6 if h_all(a)[-1]>21}
    got6={canon(r['rows']) for r in basis['strong_obstructions'] if r['n']==6}
    need(expected6==got6 and len(got6)==29,'n6 obstruction mismatch')
    candidates=0;minimal_presentations=0;expected7=set()
    for parent in rows6:
        for attachment in range(64):
            a=[row for row in parent]+[0]
            for v in range(6):
                if attachment>>v&1:a[6]|=1<<v
                else:a[v]|=1<<6
            counts=h_all(a);candidates+=1
            if counts[-1]>21 and max(counts[127^(1<<v)] for v in range(7))<=21:
                need(strong(a),'unexpected nonstrong order7 obstruction')
                expected7.add(canon(a));minimal_presentations+=1
    got7={canon(r['rows']) for r in basis['strong_obstructions'] if r['n']==7}
    need(expected7==got7 and len(got7)==13,'n7 obstruction mismatch')
    need(len(basis['strong_obstructions'])==42,'extra strong obstruction')
    # Every listed witness is independently checked on all one-vertex deletions.
    for r in basis['strong_obstructions']+basis['nonstrong_obstructions']:
        a=r['rows'];n=len(a);counts=h_all(a);full=(1<<n)-1
        need(counts[full]==r['H']>21,'bad witness H')
        cards=[counts[full^(1<<v)] for v in range(n)]
        need(cards==r['deletion_H'] and max(cards)<=21,'bad deletion witness')
    types=Counter(tuple(x['n'] for x in r['atoms']) for r in basis['nonstrong_obstructions'])
    need(types=={(3,5):6,(5,3):6,(4,4):1,(3,3,3):1},'SCC word types')
    # Six order5 strong atoms are all possible, not a selected list of counts.
    need(sum(r['n']==5 and r['strong'] for r in atlas['classes'])==6,'order5 inventory')
    for shape in [(3,5),(5,3)]:
        words=[r for r in basis['nonstrong_obstructions'] if tuple(x['n'] for x in r['atoms'])==shape]
        identities={tuple((x['n'],x['degree_code']) for x in r['atoms']) for r in words}
        need(len(identities)==6,'duplicate nonstrong SCC word')
    print(json.dumps({'status':'PASS','n6_independent_orbit_cover':len(coverage),
           'n7_augmentation_presentations':candidates,
           'n7_minimal_presentations':minimal_presentations,
           'n6_strong_obstructions':sorted(got6),'n7_strong_obstructions':sorted(got7),
           'nonstrong_shapes':{str(k):v for k,v in types.items()},
           'all_obstruction_vertex_deletions_checked':56,
           'all_order_basis_logic':'Moon excludes strong obstructions>=8 because strong order7 has H>=25; all nonstrong minimal words use good atoms and <=3 components.'},indent=2))

if __name__=='__main__':main()
