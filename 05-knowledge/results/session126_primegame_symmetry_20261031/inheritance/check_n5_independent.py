#!/usr/bin/env python3
"""Portable independent direct permutation/HP audit of the inherited n5 atlas.
Source expectations: 8f76b84d2960f9e14ad62a154b8533c77593702f,
05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/census_n5.json.
No data files or repository checkout are required.
No imports from producer, integer-only standard library, optimized-mode-safe gates.
"""
import itertools,json
from collections import defaultdict,Counter

def need(ok,msg):
 if not ok: raise RuntimeError(msg)
n=5
perms=list(itertools.permutations(range(n)))
pairs=list(itertools.combinations(range(n),2)); pos={p:i for i,p in enumerate(pairs)}
free=[(4,0),(3,0),(2,0),(4,1),(3,1),(4,2)]

def arc(code,i,j):
 return ((code>>pos[min(i,j),max(i,j)])&1) == (i<j)
def relabel(code,p):
 return sum(1<<k for k,(i,j) in enumerate(pairs) if arc(code,p[i],p[j]))
def canon(code): return min(relabel(code,p) for p in perms)
def make(mask):
 dirs={(i,i-1) for i in range(1,n)}
 dirs.update((y,x) if mask>>k&1 else (x,y) for k,(x,y) in enumerate(free))
 return sum(1<<k for k,(i,j) in enumerate(pairs) if (i,j) in dirs)
def reflect(mask):
 return sum(((mask>>free.index((n-1-y,n-1-x)))&1)<<k for k,(x,y) in enumerate(free))
def hp(code): return sum(all(arc(code,p[i],p[i+1]) for i in range(n-1)) for p in perms)
def aut(code): return sum(relabel(code,p)==code for p in perms)

classes=defaultdict(list); nodes=defaultdict(list)
for mask in range(64):
 c=canon(make(mask));m=min(c,canon(1023^make(mask)))
 classes[c].append(mask);nodes[m].append(mask)
need(len({canon(c) for c in range(1024)})==12,'full labelled tournament census')
need(len(classes)==12 and len(nodes)==10,'class/node counts')
sc=sum(canon(1023^c)==c for c in classes)
need(sc==8,'self-converse count')
ordered=sorted(nodes.values(),key=min); lookup={t:i for i,f in enumerate(ordered) for t in f}
pure=[i for i,f in enumerate(ordered) if all(reflect(t)==t for t in f)]
blue=Counter();black=Counter()
for t in range(32):
 u,v=sorted((lookup[t],lookup[t^63])); (blue if reflect(t)==t else black)[u,v]+=1
# Pinned compact expectations from the inherited atlas; computation is independent.
inherited={'pi': [0, 1, 2, 3, 4, 5, 2, 1, 2, 3, 6, 7, 5, 6, 3, 3, 8, 9, 2, 6, 2, 6, 2, 5, 2, 6, 5, 3, 1, 9, 1, 1, 4, 5, 5, 6, 5, 6, 6, 9, 2, 1, 3, 3, 6, 9, 7, 6, 2, 6, 1, 9, 3, 7, 3, 6, 2, 5, 1, 1, 3, 6, 3, 5], 'blue_edges': [[0, 5, 1], [1, 3, 1], [3, 5, 1], [3, 7, 1], [5, 6, 1], [6, 7, 1], [6, 8, 1], [7, 9, 1]], 'black_edges': [[1, 2, 2], [1, 4, 2], [1, 5, 2], [1, 6, 2], [2, 3, 2], [2, 5, 2], [2, 6, 2], [2, 9, 2], [3, 3, 2], [3, 6, 2], [5, 9, 2], [6, 6, 2]], 'classes': [{'representative': 0, 'H': 1, 'aut': 1}, {'representative': 1, 'H': 9, 'aut': 1}, {'representative': 2, 'H': 5, 'aut': 1}, {'representative': 3, 'H': 11, 'aut': 1}, {'representative': 4, 'H': 3, 'aut': 3}, {'representative': 5, 'H': 9, 'aut': 1}, {'representative': 8, 'H': 5, 'aut': 1}, {'representative': 10, 'H': 13, 'aut': 1}, {'representative': 11, 'H': 15, 'aut': 5}, {'representative': 16, 'H': 3, 'aut': 3}, {'representative': 17, 'H': 15, 'aut': 3}, {'representative': 32, 'H': 3, 'aut': 3}]}
need(lookup==dict(enumerate(inherited['pi'])),'all64 inherited node fibres')
need(dict(blue)=={(u,v):m for u,v,m in inherited['blue_edges']},'blue adjacency')
need(dict(black)=={(u,v):m for u,v,m in inherited['black_edges']},'black adjacency')
need(pure==[0,7,8],'pureblue nodes')
for row in inherited['classes']:
 code=make(row['representative'])
 need(hp(code)==row['H'] and aut(code)==row['aut'],'Hamilton/automorphism counts')
source_singletons=[f[0] for i,f in enumerate(ordered) if i in pure and len(f)==1]
need(len(source_singletons)==2 and (source_singletons[0]^source_singletons[1]).bit_count()==1,'PB singleton distance1')
requested=[7,17,55]; distances=[(a^b).bit_count() for a,b in itertools.combinations(requested,2)]
need(1 not in distances,'requested mask-distance hostile')
print(json.dumps({'status':'PASS','universe':{'labelled_tournaments':1024,'fixed_path_tilings':64,'permutations':120},'classes':12,'SC_classes':sc,'merged_nodes':10,'pureblue_nodes':pure,'pureblue_fibres':[ordered[i] for i in pure],'pureblue_H_Aut':[[hp(make(ordered[i][0])),aut(make(ordered[i][0]))] for i in pure],'line_counts':{'blue':sum(blue.values()),'black':sum(black.values()),'black_loops':sum(v for (a,b),v in black.items() if a==b)},'requested_masks':requested,'requested_nodes':[lookup[t] for t in requested],'requested_pairwise_Hamming_distances':distances,'signed_coordinate_gauge_can_select_three_PB_representatives':False},sort_keys=True,indent=2))
