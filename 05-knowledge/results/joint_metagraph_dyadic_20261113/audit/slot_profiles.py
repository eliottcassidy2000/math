#!/usr/bin/env python3
"""Independent full-slot profiles on the frozen root n3..5 class/line census.
No producer module is imported. Two exact methods count every ordered cover:
permutations with break subsets; subset Hamilton DP with ordered partitions.
"""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import comb
from pathlib import Path
import json
import subprocess


def need(ok,msg):
    if not ok: raise RuntimeError(msg)


def profiles(n,code):
    pairs=list(combinations(range(n),2))
    pos={e:i for i,e in enumerate(pairs)}
    def arc(a,b):
        i,j=sorted((a,b))
        return bool((code>>pos[i,j])&1)==(a<b)
    literal=[0]*(n+1)
    for p in permutations(range(n)):
        bad=sum((not arc(p[i],p[i+1]))<<i for i in range(n-1))
        for cuts in range(1<<(n-1)):
            if cuts&bad==bad:
                literal[1+cuts.bit_count()]+=1
    paths=[[0]*n for _ in range(1<<n)]
    hp=[0]*(1<<n)
    hp[0]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v)&1: continue
            smaller=mask^(1<<v)
            paths[mask][v]=1 if not smaller else sum(paths[smaller][u] for u in range(n)
                                                        if (smaller>>u)&1 and arc(u,v))
        hp[mask]=sum(paths[mask])
    covers=[[0]*(1<<n) for _ in range(n+1)]
    covers[0][0]=1
    for j in range(1,n+1):
        for mask in range(1,1<<n):
            sub=mask
            while sub:
                covers[j][mask]+=hp[sub]*covers[j-1][mask^sub]
                sub=(sub-1)&mask
    dp=[covers[j][-1] for j in range(n+1)]
    need(literal==dp,('ordered-cover methods',n,code))
    # Exact conversion from binomial to monomial basis.
    ordinary=[Fraction(0)]*(n+1)
    falling=[Fraction(1)]
    for j in range(1,n+1):
        new=[Fraction(0)]*(len(falling)+1)
        for i,a in enumerate(falling):
            new[i]-=a*(j-1)/j
            new[i+1]+=a/j
        falling=new
        for i,a in enumerate(falling): ordinary[i]+=literal[j]*a
    need(all(a.denominator==1 for a in ordinary),('integral monomial control',n,code))
    ordinary=list(map(int,ordinary))
    for t in range(n+2):
        need(sum(literal[j]*comb(t,j) for j in range(min(t,n)+1))==
             sum(a*t**j for j,a in enumerate(ordinary)),('basis conversion',n,code,t))
    score=sorted(sum(arc(u,v) for v in range(n) if u!=v) for u in range(n))
    return literal,ordinary,score


def main():
    folder=Path(__file__).resolve().parent
    data=json.loads((folder/'small.normal.json').read_text())
    rows=[]
    for old in data['rows']:
        n=old['n']
        if n>5:continue
        nodes=[]; groups=defaultdict(list)
        for meta in old['node_metadata']:
            profile,Q,score=profiles(n,meta['code'])
            need(profile[1]==meta['H'],('H metadata',n,meta['code']))
            node=dict(meta)
            node.update(ordered_nonempty_covers=profile[1:],Q_monomial_coefficients=Q,score=score)
            nodes.append(node)
            groups[tuple(Q)].append(meta['code'])
        quotient=sorted(groups.values(),key=lambda g:g[0])
        lookup={c:i for i,g in enumerate(quotient) for c in g}
        bycode={x['code']:x for x in nodes}
        lumped={}
        for c,node in bycode.items():
            lumped[c]=[sum(Fraction(node['neighbors'].get(str(t),0),node['degree']) for t in g)
                       for g in quotient]
        failures=[]
        for g in quotient:
            for a,b in combinations(g,2):
                if lumped[a]!=lumped[b]:
                    failures.append({'codes':[a,b],'rows':[list(map(str,lumped[a])),list(map(str,lumped[b]))]})
        item={'n':n,'merged_classes':len(nodes),'Q_classes':len(groups),
              'Q_blocks':quotient,'node_metadata_and_profiles':nodes,
              'Q_partition_lumpable':not failures,'lumpability_failures':failures}
        if n==5:
            a,b=bycode[8],bycode[10]
            need(a['Q_monomial_coefficients']==b['Q_monomial_coefficients'],'H9 same full Q')
            need(a['score']==b['score'] and a['H']==b['H']==9 and a['aut']==b['aut']==1
                 and a['degree']==b['degree']==9 and a['self_converse'] and b['self_converse'],
                 'strengthened same-observer control')
            need(a['neighbors'].get('0',0)==0 and b['neighbors'].get('0',0)==1,'TT incidence difference')
            actual=sum(Fraction(bycode[0]['neighbors'].get(str(j),0),bycode[0]['degree'])*
                       Fraction(bycode[j]['neighbors'].get('0',0),bycode[j]['degree']) for j in bycode)
            source_group=quotient[lookup[0]]
            masses=[sum(bycode[c]['degree'] for c in g) for g in quotient]
            P=[[sum(bycode[c]['neighbors'].get(str(t),0) for c in g for t in h)/Fraction(masses[i])
                for h in quotient] for i,g in enumerate(quotient)]
            qret=sum(P[lookup[0]][j]*P[j][lookup[0]] for j in range(len(quotient)))
            need(actual==Fraction(1,9) and qret==Fraction(1,18),'Q quotient return loss')
            item['transitive_return2']=str(actual)
            item['Q_stationary_quotient_return2']=str(qret)
            need(len(groups)==8 and sorted(g for g in quotient if len(g)>1)==[[2,5],[8,10]],'exact n5 Q collisions')
        rows.append(item)
    need([r['Q_classes'] for r in rows]==[2,3,8],'small exact Q boundary')
    commit='5479c530c7907b5442dd885cbd6b47d76b41b8ba'
    root=folder.parents[3]
    source_paths=[
      '01-canon/theorems/THM-5153-exact-slot-pairings-of-merged-tournament-classes-and-self-converse-structure-to-order-nine.md',
      '05-knowledge/results/kps_s28_exact_pairing_sc_census_20261113/README.md',
      '05-knowledge/results/session310_positive_clocks_20261113/NEUTRAL_TAILS.md']
    pins={}
    for path in source_paths:
        raw=subprocess.check_output(['git','show',commit+':'+path],cwd=root)
        pins[path]={'commit':commit,'sha256':sha256(raw).hexdigest()}
    return {'status':'FINITE-EXACT: full Q independently counted for every merged class n3..5',
            'methods':['all vertex permutations and all break subsets','subset Hamilton path DP and ordered-partition convolution'],
            'root_census_sha256':sha256((folder/'small.normal.json').read_bytes()).hexdigest(),
            'source_pins':pins,'rows':rows,'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
