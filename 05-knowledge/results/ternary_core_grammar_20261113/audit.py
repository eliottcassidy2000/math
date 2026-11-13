"""Independent literal receipts, complete source-fibre summaries and scope checks."""
from pathlib import Path
from itertools import permutations
import json,gzip,hashlib
from triangle_probe import HERE,need,save
from atomic_ears import pure_ends,profile
from peeling import induced,strong,near,parse,cc

def read(name):return json.loads(gzip.decompress((HERE/(name+'.json.gz')).read_bytes()))

def direct_profile(T,r):
    h=s=F=0
    for p in permutations(range(len(T))):
        good=all(T[a]>>b&1 for a,b in zip(p,p[1:]))
        h+=good;s+=good and p[0]==r
        k=p.index(r)+1
        if all(T[a]>>b&1 for a,b in zip(p[:k],p[1:k])) and all(T[a]>>b&1 for a,b in zip(p[k:],p[k+1:])):F+=1
    q=sum(pure_ends(induced(T,[v for v in range(len(T)) if v!=r]))[-1])
    return h,s,q,F

def validate_tree(t):
    T=tuple(t['rows']);r=t['root'];need(strong(T),'tree strong')
    if t['kind']=='atom':need(near(T) is not None,'atom recognition');return
    keep=t['keep'];removed=t['removed'];v=t['shared']
    need(r in keep and v in keep and set(keep).isdisjoint(removed),'protected names')
    need(set(keep)|set(removed)==set(range(len(T))),'whole parent/atom union')
    need(induced(T,keep)==tuple(t['parent']['rows']),'literal parent restriction')
    need(induced(T,[v]+removed)==tuple(t['atom_rows']),'literal atom restriction')
    need(strong(tuple(t['atom_rows'])),'strong atom guard')
    need(keep[t['parent']['root']]==r,'persistent root identity')
    validate_tree(t['parent'])

def run():
    d=read('triangle_probe');p=read('peeling');a=read('atomic_ears');n=read('neutral_entry')
    need(d['rooted_inputs']==2618 and d['cross_matrices']==10014992,'complete declared universe')
    controls=0
    for key,w in p['witnesses'].items():
        T=tuple(w['child_rows']);A=tuple(w['parent_rows']);r=len(A)-1
        need(induced(T,list(range(len(A))))==A,'parent retained')
        u,v=len(A),len(A)+1
        need(T[r]>>u&1 and T[u]>>v&1 and T[v]>>r&1,'literal attached triangle')
        need(strong(T) and sum(pure_ends(T)[-1])==int(key),'independent full target recount')
        need(parse(len(A),w['parent_code']) is not None,'recursive source is in grammar')
        controls+=1
    expected={str(h) for h in range(9,2882,2) if h!=21}
    need(set(p['witnesses'])==expected,'entire claimed finite interval')
    for rec in a['records']:validate_tree(rec['tree'])
    literal=[]
    for rows,r in ((n['old_rows'],n['old_root']),(n['new_rows'],5),(n['H79_rows'],n['H79_root'])):
        z=direct_profile(tuple(rows),r);need(z==profile(tuple(rows),r),'permutation/cover profile check');literal.append(z)
    need(cc.canonical(tuple(n['old_rows']))==131076,'the blocked count25 class')
    fibres=[]
    for h in sorted({r['core_H'] for r in d['records'] if r['core_H']<45}):
        rr=[r for r in d['records'] if r['core_H']==h];vals={v for r in rr for v in map(int,r['spectrum'])}
        fibres.append(dict(H=h,rooted_classes=len(rr),ordinary_classes=len({(r['n'],r['ordinary_code']) for r in rr}),missing_offsets=[e for e in (0,2,4) if 3*h+e not in vals]))
    inherited=HERE.parent/'strong_H_spectrum_m9_values_kps_S134.out'
    full=set(map(int,inherited.read_text().split()))
    alltri={v for r in d['records'] if r['n']==7 for v in map(int,r['spectrum'])}
    recursive={v for r in d['records'] if r['n']==7 and parse(7,r['ordinary_code']) is not None for v in map(int,r['spectrum'])}
    need(alltri<=full and recursive<=alltri,'inherited complete order9 spectrum control')
    result=dict(status='PASS',independent_target_recounts=controls,rooted_parse_certificates=len(a['records']),literal_profiles=literal,complete_strong_source_fibres_below45=fibres,order9_spectrum_comparison=dict(inherited_values=len(full),all_triangle_values=len(alltri),all_triangle_missing=sorted(full-alltri),recursive_triangle_values=len(recursive),recursive_triangle_missing=sorted(full-recursive)),arithmetic_scope='All fixed-target negatives have their all-order source bound. Finite positive range is not a new coverage record.')
    save('audit',result)
    print('PASS',controls,'literal target recounts;',len(a['records']),'rooted certificates')
    print('ORDER9',result['order9_spectrum_comparison'])

if __name__=='__main__':run()
