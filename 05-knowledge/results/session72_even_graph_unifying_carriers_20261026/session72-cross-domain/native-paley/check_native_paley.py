#!/usr/bin/env python3
"""Native three-port alignment hostile: exact P7 adjacency, no repo imports."""
from itertools import permutations, combinations
from pathlib import Path
import hashlib,json,sys

def need(ok,label):
    if not ok:raise RuntimeError(label)
if '--inject-failure' in sys.argv:need(False,'injected native Paley failure')
r=Path(__file__).resolve().parent
pins=json.loads((r/'SOURCE_PINS.json').read_text())
for p in pins['files']:
    b=(r/p['snapshot']).read_bytes()
    need(hashlib.sha256(b).hexdigest()==p['sha256'],'source SHA256')
    need(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==p['git_blob'],'source Git blob')
V=tuple(range(7));R={1,2,4}
def arc(i,j):return (j-i)%7 in R
all_perms=list(permutations(V))
auts=[p for p in all_perms if all(arc(i,j)==arc(p[i],p[j]) for i,j in permutations(V,2))]
affine={tuple((a*x+b)%7 for x in V) for a in R for b in V}
need(len(auts)==21 and set(auts)==affine,'full exact affine Aut(P7)')

def orbit(x):return {tuple(p[i] for i in x) for p in auts}
def canon(x):return min(orbit(x))
def out_common(x):return [v for v in V if all(arc(i,v) for i in x)]
def in_common(x):return [v for v in V if all(arc(v,i) for i in x)]
x=(0,1,3);y=(0,1,5)
need(y not in orbit(x),'different triple orbit')
witnesses=[]
for i,j in combinations(range(3),2):
    w=next(p for p in auts if (p[x[i]],p[x[j]])==(y[i],y[j]))
    witnesses.append({'positions':[i,j],'permutation':w})
need(out_common(x)==[] and out_common(y)==[2],'native common-out guard')
need(in_common(x)==[6] and in_common(y)==[],'native common-in guard')
need(sum(all(p[i]==i for i in [0,1]) for p in auts)==1,'two-point stabilizer trivial')
triple_orbits={canon(t) for t in permutations(V,3)}
need(len(triple_orbits)==10 and all(len(orbit(t))==21 for t in triple_orbits),'ten ordered triple orbits')
pair_records={}
for t in permutations(V,3):
    record=tuple(canon((t[i],t[j])) for i,j in combinations(range(3),2))
    pair_records.setdefault(record,set()).add(canon(t))
need(sorted(map(len,pair_records.values()))==[1]*6+[2]*2,'exact pair-record ambiguity')
cycles=[]
for S in combinations(V,3):
    if all(sum(arc(i,j) for j in S)==1 for i in S):cycles.append(S)
classes={0:[],1:[]}
for S in cycles:
    q=len(out_common(S));need(q in [0,1],'cycle common-out count');classes[q].append(S)
need(len(cycles)==14 and all(len(c)==7 for c in classes.values()),'two seven-cycle families')
for q,C in classes.items():
    orb={tuple(sorted(p[i] for i in C[0])) for p in auts}
    need(set(C)==orb,'each common-out class is one actual Aut orbit')
    for pair in combinations(V,2):need(sum(set(pair)<=set(S) for S in C)==1,'native Fano pair incidence')
# Prove the point-stabilizer upper bound from its actual adjacency data.
O=tuple(j for j in V if arc(0,j));I=tuple(j for j in V if arc(j,0))
need(O==(1,2,4) and I==(3,5,6),'root neighbour sets')
need(all(sum(arc(i,j) for j in O)==1 for i in O),'root outneighbours induce C3')
profiles={i:tuple(int(arc(i,j)) for j in O) for i in I}
need(len(set(profiles.values()))==3,'in-neighbour profiles separate')
# Exact F8 coordinate interpretation, not merely an equal-cardinality label.
def mul(a,b):
    out=0
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return out
powers=[1]
for j in range(1,7):powers.append(mul(powers[-1],2))
need(powers==[1,2,4,3,6,7,5],'F8 exponent coordinates')
def trace(a):return a^mul(a,a)^mul(mul(a,a),mul(a,a))
for i,j in permutations(V,2):need(arc(i,j)==(trace(powers[(j-i)%7])==0),'literal trace-ratio arc predicate')
def rank(vectors):
    piv={}
    for v in vectors:
        while v:
            lead=v.bit_length()-1
            if lead in piv:v^=piv[lead]
            else:piv[lead]=v;break
    return len(piv)
field_x=[powers[i] for i in x];field_y=[powers[i] for i in y]
need(rank(field_x)==2 and rank(field_y)==3,'native triple F2 dependence')
need(field_x[0]^field_x[1]^field_x[2]==0 and field_y[0]^field_y[1]^field_y[2]==4,'field triple sums')
for t in [field_x,field_y]:
    for pair in combinations(t,2):need(rank(pair)==2,'pairwise F2 ranks agree')
for q,C in classes.items():
    for S in C:need(rank([powers[i] for i in S])==2+q,'cyclic family common-out equals rank-minus2')
print(json.dumps({'status':'PASS','input_pin':pins['pin'],'permutations_exhausted':len(all_perms),'automorphisms':len(auts),'F8':{'exponent_masks':powers,'triple_masks':[field_x,field_y],'triple_ranks':[rank(field_x),rank(field_y)],'triple_sums':[field_x[0]^field_x[1]^field_x[2],field_y[0]^field_y[1]^field_y[2]],'arc_trace_checks':42},'ordered_triples':210,'ordered_triple_orbits':len(triple_orbits),'pair_record_fibre_sizes':sorted(map(len,pair_records.values())),'unordered_cyclic_triples':len(cycles),'cyclic_orbit_sizes':[len(classes[0]),len(classes[1])],'hostile':{'first':x,'second':y,'pairwise_witnesses':witnesses,'common_out':[out_common(x),out_common(y)],'common_in':[in_common(x),in_common(y)]},'root_in_profiles':profiles,'source_files_verified':len(pins['files'])},sort_keys=True,indent=2))
