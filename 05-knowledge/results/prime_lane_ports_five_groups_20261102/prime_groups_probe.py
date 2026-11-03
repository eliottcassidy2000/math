#!/usr/bin/env python3
"""Exact recovery of the five prime groups, physical graph ports, square carry.
No third-party mathematical libraries. An inherited order-seven mask list is
used only as a candidate bank; each selected tournament is independently
verified by subset DP, and the H=49 anchor also by all vertex permutations.
"""
from collections import Counter
from itertools import combinations,permutations
from math import isqrt
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
GROUPS=[[2,3,5,7,11,13],[17,19,23,31,47],[29,41,59,71],
        [37,43,53,61,67],[73,79,83,89,97]]

def require(c,m):
    if not c:raise RuntimeError(m)

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def square(n):return isqrt(n)**2==n

def graph_edges(vertices):
    return [(p,q,isqrt(p+q)) for p,q in combinations(vertices,2) if square(p+q)]

def components(vertices,edges):
    adj={v:set() for v in vertices}
    for p,q,*_ in edges:adj[p].add(q);adj[q].add(p)
    result=[];seen=set()
    for v in vertices:
        if v in seen:continue
        todo=[v];comp=set()
        while todo:
            x=todo.pop()
            if x in comp:continue
            comp.add(x);todo.extend(adj[x]-comp)
        result.append(sorted(comp));seen|=comp
    return sorted(result,key=lambda x:(-len(x),x)),adj

def matrix(n,mask):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(mask>>k)&1;a[j][i]=1-a[i][j]
    return a

def hamilton(a):
    n=len(a)
    if not n:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            c=dp[mask][v]
            if not c:continue
            for w in range(n):
                if a[v][w] and not(mask>>w&1):dp[mask|1<<w][w]+=c
    return sum(dp[-1])

def sccs(a):
    n=len(a);reach=[{i}|{j for j in range(n) if a[i][j]} for i in range(n)]
    for k in range(n):
        for i in range(n):
            if k in reach[i]:reach[i]|=reach[k]
    seen=set();cs=[]
    for i in range(n):
        if i in seen:continue
        c={j for j in range(n) if j in reach[i] and i in reach[j]};seen|=c;cs.append(sorted(c))
    return cs

def join(a,b):
    m,n=len(a),len(b)
    return [row+[1]*n for row in a]+[[0]*m+row for row in b]

def rowstrings(a):return [''.join(map(str,row)) for row in a]

def small_witness(h,n):
    for mask in range(1<<(n*(n-1)//2)):
        a=matrix(n,mask)
        if hamilton(a)==h:return a
    raise RuntimeError('no small witness')

def diagram(primes,edges,adj):
    from html import escape
    colors=['#2764a5','#a14ca1','#cf8425','#168a75','#d45162']
    color={p:colors[i] for i,g in enumerate(GROUPS) for p in g}
    left=[p for p in primes if p==2 or p%4==1];right=[p for p in primes if p%4==3]
    xy={p:(260,145+43*i) for i,p in enumerate(left)}|{p:(690,145+43*i) for i,p in enumerate(right)}
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="1260" height="920" viewBox="0 0 1260 920">',
           '<rect width="1260" height="920" fill="#f8fafc"/>']
    try:
        from PIL import Image,ImageDraw,ImageFont
        image=Image.new('RGB',(1260,920),'#f8fafc');draw=ImageDraw.Draw(image)
        def font(s):
            try:return ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',s)
            except OSError:return ImageFont.load_default()
    except ImportError:image=draw=None
    def line(x,y,x2,y2,col='#b3bdcc',w=2):
        parts.append(f'<line x1="{x}" y1="{y}" x2="{x2}" y2="{y2}" stroke="{col}" stroke-width="{w}"/>')
        if draw:draw.line((x,y,x2,y2),fill=col,width=w)
    def text(x,y,t,s=19,col='#14213d'):
        parts.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{s}" fill="{col}">{escape(t)}</text>')
        if draw:draw.text((x,y-s),t,font=font(s),fill=col)
    def dot(x,y,r,col):
        parts.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{col}"/>')
        if draw:draw.ellipse((x-r,y-r,x+r,y+r),fill=col)
    text(30,43,'The actual five prime groups: names and square-sum ports survive',27)
    text(30,80,'First 25 primes = all primes below 100. An edge means p + q is an integer square.',21)
    text(85,116,'A: 2 and primes 1 mod 4 (12)',18);text(580,116,'B: primes 3 mod 4 (13)',18)
    for p,q,s in edges:
        x,y=xy[p];x2,y2=xy[q];line(x,y,x2,y2)
    for p in primes:
        x,y=xy[p];dot(x,y,9,color[p]);text(x-48 if p in left else x+18,y+7,str(p),20,color[p])
        if not adj[p]:text(x-143 if p in left else x+64,y+7,'isolated',16,'#b91c1c')
    text(825,145,'Recovered groups',23)
    for i,g in enumerate(GROUPS):
        text(825,193+i*86,f'G{i+1}: '+','.join(map(str,g)),18,colors[i])
        if i==0:label='prime divisors of Fi22'
        elif i==1:label='new prime divisors in Baby Monster'
        elif i==2:label='new prime divisors in Monster'
        elif i==3:label='first five primes outside Monster support'
        else:label='remaining five below 100'
        text(825,220+i*86,label,15)
    text(825,665,'24 edges; components 22 + 1 + 1 + 1',18)
    text(825,696,'Cycle rank 3; no spanning path.',18)
    line(30,740,1230,740,'#cbd5e1')
    text(30,778,'Exact tournament bridge: choose K with H(K)=49; send T to the ordered join K -> T.',21)
    text(30,816,'H becomes 49H, tournament order increases by 7, and each square root s becomes 7s.',20)
    text(30,854,'5 + 11 = 4^2  ->  245 + 539 = 28^2.  Multiplication by 7 alone gives 35 + 77 = 112.',20)
    text(30,891,'The 24 odd primes lift to count labels 49p; an odd label for prime 2 cannot keep its edges to 23 and 47.',17)
    parts.append('</svg>');(HERE/'prime_groups.svg').write_text('\n'.join(parts)+'\n')
    if image:image.save(HERE/'prime_groups.png')

def main():
    primes=[n for n in range(2,100) if prime(n)]
    require(len(primes)==25 and sorted(sum(GROUPS,[]))==primes,'not the exact first25 partition')
    cited_orders={'Fi22':64561751654400,
                  'B_from_ATLAS_2B':8309562962452852382355161088000000//2,
                  'M':808017424794512875886459904961710757005754368000000000}
    support_orders={}
    for key,value in cited_orders.items():
        remaining=value;support=[]
        for p in primes:
            if remaining%p==0:support.append(p)
            while remaining%p==0:remaining//=p
        require(remaining==1,'group order factor outside listed primes')
        support_orders[key]=support
    require(support_orders['Fi22']==GROUPS[0],'Fi22 shell mismatch')
    require(support_orders['B_from_ATLAS_2B']==sorted(sum(GROUPS[:2],[])),'B shell mismatch')
    require(support_orders['M']==sorted(sum(GROUPS[:3],[])),'Monster shell mismatch')
    edges=graph_edges(primes);comps,adj=components(primes,edges)
    # Independent root-first enumeration, including the no-loop guard.
    byroot=sorted((p,s*s-p,s) for s in range(1,isqrt(194)+1) for p in primes
                  if p<s*s-p and s*s-p in primes)
    require(edges==byroot,'root-first and pair-first edge lists disagree')
    require(len(edges)==24 and list(map(len,comps))==[22,1,1,1],'prime census drift')
    require([p for p in primes if not adj[p]]==[37,43,67],'wrong isolated primes')
    group={p:i for i,g in enumerate(GROUPS) for p in g}
    incidence=[[0]*5 for _ in range(5)]
    for p,q,s in edges:incidence[group[p]][group[q]]+=1;incidence[group[q]][group[p]]+=1
    profile={str(p):[sum(group[q]==i for q in adj[p]) for i in range(5)] for p in primes}
    require(profile['37']==[0]*5 and profile['53']==[1,1,0,0,0],'physical port hostile drift')
    require(all(incidence[i][j] for i,j in [(0,1),(1,2),(2,0)]),'quotient triangle missing')
    require(not any(q in adj[p] and r in adj[p] and r in adj[q] for p,q,r in combinations(primes,3)),
            'prime graph triangle contradicts bipartition')
    firstpartners={str(p):next(q for q in range(2,500) if q!=p and prime(q) and square(p+q)) for p in [37,43,67]}
    # The existing complete order-seven class list is a candidate source only.
    source=REPO/'05-knowledge/results/n7_class_reps_boxeph_S152.txt'
    targets=set(primes)-{2,7};bank={};anchor=None;anchor_mask=None
    for mask in map(int,source.read_text().split()):
        a=matrix(7,mask);h=hamilton(a)
        if h==49 and anchor is None:anchor=a;anchor_mask=mask
        if h in targets and h not in bank:
            pieces=sccs(a);nontrivial=[c for c in pieces if len(c)>1]
            require(len(nontrivial)==1,'prime H has multiple nontrivial SCCs')
            c=nontrivial[0];strong=[[a[i][j] for j in c] for i in c]
            require(hamilton(strong)==h and len(sccs(strong))==1,'strong extraction failed')
            bank[h]={'source_order7_mask':mask,'transport_new_to_old':c,'n':len(c),'rows':rowstrings(strong)}
    require(set(bank)==targets,'finite prime lane missing candidate')
    require(anchor is not None and len(sccs(anchor))==1,'49 anchor not strong')
    literal=sum(all(anchor[p[i]][p[i+1]] for i in range(6)) for p in permutations(range(7)))
    require(literal==49,'independent literal anchor verification')
    carry343=matrix(8,136003753)
    require(hamilton(carry343)==343 and len(sccs(carry343))==1,'343 carry failed')
    literal343=sum(all(carry343[p[i]][p[i+1]] for i in range(7)) for p in permutations(range(8)))
    require(literal343==343,'literal 343 carry failed')
    # Odd H assigned to prime2, with all other H labels fixed at49p, cannot
    # retain its three edges to primes7,23,47. Enumerate exact factorizations.
    two_edge_candidates=set()
    for u in range(1,isqrt(784)+1):
        if 784%u:continue
        v=784//u
        if (v-u)%2:continue
        a,b=(v-u)//2,(v+u)//2
        if a%2 or b%2:continue
        h=a*a-343
        if h>0 and h%2:two_edge_candidates.add(h)
    require(two_edge_candidates=={8873},'extension candidate enumeration changed')
    require(all(not square(h+49*47) for h in two_edge_candidates),'unexpected prime2 extension')
    controls=[]
    for h,n in [(5,4),(11,5)]:
        a=small_witness(h,n);j=join(anchor,a);v=hamilton(j)
        require(v==49*h,'joined DP not equal product')
        controls.append({'source_H':h,'source_n':n,'source_rows':rowstrings(a),
                         'target_H':v,'target_n':len(j),'target_rows':rowstrings(j)})
    require(square(5+11) and square(49*5+49*11) and not square(7*5+7*11),'square dilation hostile')
    require(all(square(49*p+49*q) for p,q,s in edges),'49 root map failed')
    require(all(not square(7*p+7*q) for p,q,s in edges),'7 should destroy every edge')
    even_square_residues={x*x%16 for x in range(0,16,2)}
    differences={(a-b)%16 for a in even_square_residues for b in even_square_residues}
    require(even_square_residues=={0,4} and 8 not in differences,'mod16 obstruction failed')
    require(all((24*c*c)%16==8 for c in range(1,16,2)),'odd-square scaling obstruction failed')
    birth_roots=lambda p:[s for s in range(1,isqrt(2*p-1)+1) if p<s*s<2*p]
    prime_birth=lambda p:[s*s-p for s in birth_roots(p) if prime(s*s-p)]
    require(birth_roots(41)==birth_roots(43)==[7,8,9],'root-chart hostile drift')
    require(prime_birth(41)==[23] and prime_birth(43)==[],'prime-port hostile drift')
    groupdetails=[]
    for i,g in enumerate(GROUPS):
        es=graph_edges(g)
        groupdetails.append({'group':i+1,'primes':g,'internal_edges':es,
            'bipartition':[sum(p==2 or p%4==1 for p in g),sum(p%4==3 for p in g)],
            'equitable':len({tuple(profile[str(p)]) for p in g})==1})
    ogg=sorted(sum(GROUPS[:3],[]));oe=graph_edges(ogg);oc,oa=components(ogg,oe)
    require(len(oe)==14 and len(oc)==1,'Ogg tree control')
    result={'status':'FINITE-EXACT','cited_group_orders':cited_orders,'independently_factored_prime_supports':support_orders,'groups':groupdetails,'edges_root_labelled':edges,
            'components':comps,'cycle_rank':len(edges)-len(primes)+len(comps),
            'oriented_group_incidence_matrix':incidence,'per_prime_group_ports':profile,
            'first_partners_of_finite_isolates':firstpartners,
            'integer_graph_comparison':{'Q25_edges':len(graph_edges(list(range(1,26)))),
                                        'Q97_edges':len(graph_edges(list(range(1,98)))),
                                        'prime25_edges':len(edges)},
            'ogg_control':{'vertices':ogg,'edges':oe,'leaves':[p for p in ogg if len(oa[p])==1],
                           'bipartition_sizes':[6,9]},
            'prime_H_witnesses':{str(h):bank[h] for h in sorted(bank)},
            'candidate_bank_source':str(source.relative_to(REPO)),
            'candidate_bank_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'anchor49':{'n':7,'mask':anchor_mask,'rows':rowstrings(anchor),'H_subset_DP':49,'H_literal_orders':literal},
            'carry343':{'source':'05-knowledge/results/session150_hamiltonian_spectrum_20261103/status/verify_witnesses.py',
                        'n':8,'mask':136003753,'rows':rowstrings(carry343),'H_subset_DP':343,'H_literal_orders':literal343},
            'odd_prime49_graph_embedding':{'prime_vertices':primes[1:],'H_labels':[49*p for p in primes[1:]],
                       'exception7_uses_carry343':True,'prime2_extension_candidates_for_neighbors7_23':sorted(two_edge_candidates),
                       'prime2_extension_third_sum':8873+49*47,'prime2_extension_impossible':True},
            'square_carry_controls':controls,'all25_arithmetic_edges_scaled49':True,
            'all25_arithmetic_edges_destroyed_scaled7':True,
            'unrealizable_H_prime_labels':[2,7],'unrealizable_H_7p_prime_labels':[2,3],
            'odd_neighbor_obstruction':{'even_squares_mod16':sorted(even_square_residues),
                       'differences_mod16':sorted(differences),'hostile_prime_neighbors':[23,47],
                       'difference_after_any_odd_square_scale_mod16':8},
            'prime_birth_clock_hostile':{'prime_values':[41,43],'prime_ranks':[13,14],
                       'root_chart':[7,8,9],'prime_parent_lists':[[23],[]],
                       'root8_parent_transport':[23,21]},
            'scope_guard':'finite graph and selected actual witnesses only; no global spectrum completeness claim'}
    out=json.dumps(result,indent=2,sort_keys=True)+'\n';(HERE/'prime_groups.json').write_text(out);print(out,end='')
    diagram(primes,edges,adj)

if __name__=='__main__':main()
