"""All-size disjointness test for inflation cones of the98 small births.

Compress maximal unselected transitive twin runs, retaining their positive
integer lengths. Quotient automorphisms turn cone intersection into a finite
test of fixed-coordinate constraints. No length cutoff is involved.
"""
from pathlib import Path
from itertools import combinations,permutations,product
from collections import defaultdict,Counter
import json,gzip
import sympy as s
from contraction_dag import key
from joint_port_extension import cones
HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def canonical_orders(rows,mask):
    cells=defaultdict(list)
    for i,r in enumerate(rows):cells[(mask>>i&1,r.bit_count())].append(i)
    groups=[cells[c] for c in sorted(cells)];best=None;orders=[]
    for chunks in product(*(permutations(g) for g in groups)):
        order=sum((tuple(g) for g in chunks),())
        rr=tuple(sum(1<<j for j,v in enumerate(order) if rows[u]>>v&1) for u in order)
        mm=sum(1<<j for j,v in enumerate(order) if mask>>v&1);v=(len(rows),mm,*rr)
        if best is None or v<best:best=v;orders=[order]
        elif v==best:orders.append(order)
    return best,orders
def quotient(rows,mask):
    n=len(rows);adj={i:set() for i in range(n) if not mask>>i&1}
    for a,b in combinations(adj,2):
        if all(bool(rows[a]>>v&1)==bool(rows[b]>>v&1) for v in range(n) if v not in [a,b]):adj[a].add(b);adj[b].add(a)
    components=[];seen=set()
    for a in range(n):
        if a in seen:continue
        comp={a};todo=[a]
        while todo:
            for b in adj.get(todo.pop(),[]):
                if b not in comp:comp.add(b);todo.append(b)
        seen.update(comp);components.append(sorted(comp))
    for c in components:
        need(len({sum(bool(rows[u]>>v&1) for v in c if u!=v) for u in c})==len(c),'transitive twin component')
        for v in range(n):
            if v not in c:need(len({bool(rows[u]>>v&1) for u in c})==1,'component module')
    rr=[sum(1<<j for j,d in enumerate(components) if i!=j and rows[c[0]]>>d[0]&1) for i,c in enumerate(components)]
    mm=sum(1<<i for i,c in enumerate(components) if mask>>c[0]&1)
    return rr,mm,components
def run():
    bank=json.loads((HERE/'cyclic_port_receipt.json').read_text())['positive_screen_certificates'];st={key(r['rows'],r['selected']):r['individually_stable'] for r in bank}
    dag=json.loads((HERE/'contraction_dag.json').read_text())['records'];births=[r for r in dag if not r['children']]
    classes=defaultdict(list);records=[];literal=defaultdict(set);presented=Counter();aut_orders=Counter();hist=Counter();z=s.Symbol('z');gf=0
    for r in births:
        rows=r['rows'];mask=r['selected'];ports=st[key(rows,mask)];O=[v for v in range(len(rows)) if mask>>v&1]
        d,u=cones.kernel(rows,ports,O);need(d=={(1,)*len(ports):1},'complete birth joint kernel')
        rr,mm,comps=quotient(rows,mask);q,orders=canonical_orders(rr,mm)
        variables={i for i,c in enumerate(comps) if any(v in ports for v in c)}
        need(all(len(comps[i])==1 for i in variables),'birth variable must be singleton')
        charts=[]
        for order in orders:
            fixed=tuple(0 if i in variables else len(comps[i]) for i in order)
            charts.append(fixed);classes[q].append({'birth_record':r['record'],'fixed_coordinates':fixed})
        need(len(charts)==len(set(charts)),'marked birth stabilizer')
        # Separate literal graph construction and canonicalizer: all cone
        # members of total order <=10, including intersections at small sizes.
        for lengths in product(range(1,12-r['n']),repeat=len(ports)):
            total=r['n']+sum(lengths)-len(ports)
            if total>10:continue
            expanded,selected,_=cones.expand(rows,O,ports,lengths)
            sm=sum(1<<v for v in selected);literal[total].add(key(expanded,sm));presented[total]+=1
        hist[(r['n'],len(ports))]+=1;aut_orders[len(orders)]+=1;gf+=z**r['n']/(1-z)**len(ports)
        records.append({'birth_record':r['record'],'order':r['n'],'stable_ports':ports,'quotient':q,'components':comps,'quotient_automorphisms':len(orders),'charts':charts,'joint_defect_terms':[[*k,v] for k,v in sorted(d.items())]})
    overlaps=[];disjoint_certificates=[];pair_tests=0
    for q,charts in classes.items():
        for a,b in combinations(charts,2):
            pair_tests+=1;aa=a['fixed_coordinates'];bb=b['fixed_coordinates']
            conflicts=[{'coordinate':i,'first_fixed_length':x,'second_fixed_length':y} for i,(x,y) in enumerate(zip(aa,bb)) if x and y and x!=y]
            if conflicts:
                disjoint_certificates.append({'quotient':q,'first':a,'second':b,'conflicts':conflicts})
                continue
            w=[max(x,y,1) for x,y in zip(aa,bb)]
            overlaps.append({'quotient':q,'first':a,'second':b,'minimal_intersection_weights':w,'minimal_order':sum(w)})
    out={'scope':'Finite exact certificate for all-length cone intersections after canonical unselected-twin compression, for every stable-clone birth through8. All original stable ports inflated jointly and independently; no newly emergent ports are added to the allowed operation set.','births':len(records),'birth_order_port_histogram':[[*k,v] for k,v in sorted(hist.items())],'quotient_aut_order_histogram':dict(aut_orders),'distinct_quotients':len(classes),'chart_pair_tests':pair_tests,'overlaps':overlaps,'disjoint_certificates':disjoint_certificates,'provisional_sum_GF':str(s.factor(gf)),'records':records}
    if not overlaps:
        vals=s.series(gf,z,0,65).removeO().expand();out['disjoint_family_GF']=str(s.factor(gf));out['counts_n1_to64']=[int(vals.coeff(z,n)) for n in range(1,65)]
        for lev in json.loads((HERE/'cyclic_port_receipt.json').read_text())['levels']:need(vals.coeff(z,lev['n'])==sum(int(v) for k,v in lev['all_length_stable_port_histogram'].items() if int(k)>0),'all bank counts through8')
        n9=json.loads(gzip.decompress((HERE/'descendants_nine.json.gz').read_bytes()));need(vals.coeff(z,9)==n9['marked_descendant_types'],'complete descendant count9')
        for n,keys in literal.items():
            need(len(keys)==presented[n]==vals.coeff(z,n),'literal cone count')
            if n<=8:need(keys=={key(r['rows'],r['selected']) for r in dag if r['n']==n},'complete bank key equality')
            if n==9:need(keys=={key(r['rows'],r['selected']) for r in n9['records']},'complete descendant key equality')
        out['literal_controls']={'all_cone_members_through_order':10,'marked_graph_counts':{n:len(keys) for n,keys in literal.items()},'complete_bank_key_equality_through8':True,'complete_order9_descendant_key_equality':True}
        edges=s.factor(sum(count*k*z**(b+1)/(1-z)**k for (b,k),count in hist.items()))
        diamonds=s.factor(sum(count*z**(b+2)/(1-z)**2 for (b,k),count in hist.items() if k==2))
        ee=s.series(edges,z,0,65).removeO().expand();dd=s.series(diamonds,z,0,65).removeO().expand()
        for n in range(1,9):need(ee.coeff(z,n)==sum(len(r['children']) for r in dag if r['n']==n),'full contraction graph edge control')
        need(ee.coeff(z,9)==n9['edges_to_order8_types'],'order9 edge control')
        need(dd.coeff(z,9)==sum(len(r['children_order8_records'])==2 for r in n9['records']),'order9 diamond control')
        need(s.cancel(edges-(sum(ee.coeff(z,n)*z**n for n in range(1,9))+38*z*s.diff(z**9/(1-z),z)-215*z**9/(1-z)))==0,'all-order edge tail')
        need(s.cancel(diamonds-(19*z*s.diff(z**9/(1-z),z)-166*z**9/(1-z)))==0,'all-order diamond tail')
        out['original_port_growth_graph']={'components':{'rays':79,'integer_quadrants':19},'contraction_edge_GF':str(edges),'two_child_state_GF':str(diamonds),'edges_n1_to64':[int(ee.coeff(z,n)) for n in range(1,65)],'two_child_states_n1_to64':[int(dd.coeff(z,n)) for n in range(1,65)],'for_every_n_ge9':{'vertices':'19*n-49','edges_to_previous_level':'38*n-215','two_child_states':'19*n-166','one_child_states':117},'scope':'Only the original stable-port operations, with twin clones merged to one child type. No assertion that later emergent ports cannot add edges.'}
        need(all(vals.coeff(z,n)==19*n-49 for n in range(8,65)),'linear tail identity coefficients')
        need(s.cancel(gf-(sum(vals.coeff(z,n)*z**n for n in range(1,8))+sum([19*(z* s.diff(z**8/(1-z),z))-49*z**8/(1-z)])))==0,'exact linear-tail GF identity')
    tag='.optimized' if not __debug__ else ''
    (HERE/('birth_cones'+tag+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['records','disjoint_certificates']}),flush=True)
if __name__=='__main__':run()
