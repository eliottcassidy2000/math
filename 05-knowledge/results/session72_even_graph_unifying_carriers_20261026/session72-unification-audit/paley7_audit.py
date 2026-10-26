#!/usr/bin/env python3
"""Independent exact native noncyclic pair-orbit boundary, no producer imports."""
from itertools import permutations,combinations
import argparse,json

def need(ok,why):
    if not ok:raise RuntimeError(why)

def mul(x,y):
    z=0
    while y:
        if y&1:z^=x
        y>>=1;x<<=1
        if x&8:x^=11
    return z

def bitrank(xs):
    pivots={}
    for x in xs:
        while x:
            p=x.bit_length()-1
            if p in pivots:x^=pivots[p]
            else:pivots[p]=x;break
    return len(pivots)

def main():
    p=argparse.ArgumentParser();p.add_argument('--inject',action='store_true');args=p.parse_args()
    V=tuple(range(7));R={1,2,4}
    arc=lambda x,y:(y-x)%7 in R
    automorphisms=[p for p in permutations(V) if all(arc(x,y)==arc(p[x],p[y]) for x,y in combinations(V,2))]
    affine={tuple((a*x+b)%7 for x in V) for a in R for b in V}
    need(set(automorphisms)==affine and len(automorphisms)==21,'full native automorphism group exactly affine21')
    x=(0,1,3);y=(0,1,5)
    need(y not in {tuple(p[t] for t in x) for p in automorphisms},'native triples inequivalent')
    pair_matches=[]
    for i,j in combinations(range(3),2):
        matches=[p for p in automorphisms if (p[x[i]],p[x[j]])==(y[i],y[j])]
        need(len(matches)==1,'named pair orbit agrees with unique matching automorphism')
        pair_matches.append({'slots':[i,j],'permutation':matches[0]})
    need([p for p in automorphisms if p[0]==0 and p[1]==1]==[V],'ordered pair stabilizer trivial')
    out=lambda t:sorted(v for v in V if all(arc(u,v) for u in t))
    inc=lambda t:sorted(v for v in V if all(arc(v,u) for u in t))
    need(out(x)==[] and out(y)==[2] and inc(x)==[6] and inc(y)==[],'native common-neighbor observer separates')
    out0=[v for v in V if arc(0,v)];in0=[v for v in V if arc(v,0)]
    profiles={v:tuple(int(arc(v,w)) for w in out0) for v in in0}
    need(len(set(profiles.values()))==3,'stabilizer restriction injective by distinct incoming profiles')
    powers=[1]
    for i in range(1,7):powers.append(mul(powers[-1],2))
    need(powers==[1,2,4,3,6,7,5] and mul(powers[-1],2)==1,'F8 exponent dictionary')
    trace=lambda z:z^mul(z,z)^mul(mul(z,z),mul(z,z))
    need(all(arc(i,j)==(trace(powers[(j-i)%7])==0) for i in V for j in V if i!=j),'all42 native arcs are trace-ratio predicates')
    for a in R:
        for b in V:
            for i in V:
                z=powers[i]
                if a==2:z=mul(z,z)
                elif a==4:z=mul(mul(z,z),mul(z,z))
                need(powers[(a*i+b)%7]==mul(powers[b],z),'full native group is F2-linear semilinear1 action')
    cyclic=[t for t in combinations(V,3) if all(sum(arc(u,v) for v in t)==1 for u in t)]
    need(len(cyclic)==14,'native cyclic-triangle census')
    families={r:[t for t in cyclic if len(out(t))==r] for r in (0,1)}
    for t in cyclic:need(len(out(t))==bitrank([powers[i] for i in t])-2,'cyclic common-out rank-minus-two detector')
    for family in families.values():
        need(len(family)==7,'Fano family size')
        need(all(sum(set(e)<=set(t) for t in family)==1 for e in combinations(V,2)),'every pair once in Fano family')
    triple_orbits={min(tuple(p[t] for t in z) for p in automorphisms) for z in permutations(V,3)}
    need(len(triple_orbits)==10,'all ordered triple orbits')
    ranks=[bitrank([powers[i] for i in t]) for t in (x,y)]
    sums=[powers[t[0]]^powers[t[1]]^powers[t[2]] for t in (x,y)]
    need(ranks==[2,3] and sums==[0,4],'explicit triple dependence detector')
    if args.inject:need(out(x)==out(y),'deliberate pairwise-to-joint conflation')
    print(json.dumps({'status':'PASS','permutations_examined':5040,'automorphisms':21,'triples':[x,y],'pair_matching_automorphisms':pair_matches,'common_out':[out(x),out(y)],'common_in':[inc(x),inc(y)],'in0_profiles_to_out0':profiles,'F8':{'exponent_masks':powers,'trace_arcs_checked':42,'triple_ranks':ranks,'triple_sums':sums,'cyclic_triangles':len(cyclic),'Fano_family_sizes':[len(families[r]) for r in (0,1)],'ordered_triple_orbits':len(triple_orbits)},'scope':'native Paley7 automorphism action, no cyclic all-arity counterexample; rank-minus-two only on cyclic triangle locus'},sort_keys=True,indent=2))
if __name__=='__main__':main()
