"""Shared named-port recurrences; independent group/partition controls at12."""
import importlib.util
from collections import Counter
from support import ROOT,need,save
from hostile_embedding import tournament,order
def partition(p):
    seen=set();blocks=[]
    for i in range(len(p)):
        if i in seen:continue
        b=[];j=i
        while j not in seen:seen.add(j);b.append(j);j=p[j]
        blocks.append(tuple(sorted(b)))
    return tuple(sorted(blocks))
def join(A,B):
    parent=list(range(sum(map(len,A))))
    def find(x):
        while parent[x]!=x:x=parent[x]
        return x
    for block in A+B:
        for x in block:parent[find(x)]=find(block[0])
    groups={}
    for x in range(len(parent)):groups.setdefault(find(x),[]).append(x)
    return tuple(sorted(tuple(v) for v in groups.values()))
def regular_polynomial(G):
    # Every nontrivial stabilizer in a 3-group contains an element of order3.
    events=sorted({partition(p) for p in G if order(p)==3});coeff={tuple((i,) for i in range(12)):1}
    for E in events:
        old=list(coeff.items())
        for P,c in old:
            J=join(P,E);coeff[J]=coeff.get(J,0)-c
        coeff={P:c for P,c in coeff.items() if c}
    poly=Counter()
    for P,c in coeff.items():poly[len(P)]+=c
    return {k:v for k,v in sorted(poly.items()) if v},len(events),len(coeff)
def main():
    path=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec=importlib.util.spec_from_file_location('ports_hostile',path);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    result=[]
    for flipped in (False,True):
        T=tournament(flipped);G=h.iso_maps(T,T);poly,events,terms=regular_polynomial(G);rows=[]
        for p in (1,2,3):
            q=2**p;a=q**4-q**2;o=(q**4+2*q**2)//3
            predicted=a**3 if flipped else a**3-9*a
            regular=sum(c*q**d for d,c in poly.items())
            orbits=sum(q**len(partition(g)) for g in G)//len(G)
            expected_orbits=o**3 if flipped else (o**3+2*o)//3
            need(regular==predicted and orbits==expected_orbits,'independent partition and Burnside checks')
            rows.append({'named_ports':p,'colors':q,'regular_attachments':regular,'attachment_orbits':orbits})
        stabilizers=Counter(sum(h.permute_mask(mask,g)==mask for g in G) for mask in range(1<<12))
        need(set(stabilizers)=={3**j for j in range(4 if flipped else 5)},'all and only ternary stabilizer orders at12')
        result.append({'flipped':flipped,'automorphisms':len(G),'order3_events':events,'partition_terms':terms,
                       'regular_count_polynomial':poly,'ports':rows,'binary_stabilizer_order_counts':dict(sorted(stabilizers.items())),
                       'two_port_orbit_excess':rows[1]['attachment_orbits']-rows[0]['attachment_orbits']**2})
    levels=[]
    for q in (2,4,8):
        A=B=3;R=F=q**4-q**2;O=E=(q**4+2*q**2)//3
        for k in range(1,5):
            levels.append({'k':k,'colors':q,'aut_uniform':A,'aut_one_flip':B,'regular_uniform':R,'regular_one_flip':F,
                           'orbits_uniform':O,'orbits_one_flip':E})
            A,B,R,F,O,E=3*A**3,B*A**2,R**3-R*A**2,F*R**2,(O**3+2*O)//3,E*O**2
    save('attachments.json',{'status':'EXACT_ALL_LEVEL_RECURRENCE_WITH_INDEPENDENT12_VERTEX_CHECKS','actual_controls':result,
      'recurrence_levels':levels,'scope':'unrestricted orientations to named fixed ports; no Hamilton-count or strong-output filter'})
    print('PASS shared-port recurrence: independent Burnside and partition checks for1,2,3 ports',flush=True)
if __name__=='__main__':main()
