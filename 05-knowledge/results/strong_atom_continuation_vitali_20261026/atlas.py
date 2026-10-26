"""Complete strong k=1 responses through seven; all run arities through five.

Default prints a compact receipt. --write also creates the complete frozen atlas.
No finite census is an all-order minimality theorem.
"""
from collections import Counter, defaultdict
from itertools import permutations
from pathlib import Path
import json, sys
from core import (require, edge, converse, relabel, sigma, complement, presentation,
                  free_mask, canonical, is_strong, is_prime, paths, hamilton,
                  automorphisms, substitute)

def dump(obj): return json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n'

def main():
    atoms=[]; summary=[]; first_lift=None; first_tau=None
    for n in range(1,8):
        total=1<<max(0,(n-1)*(n-2)//2)
        records=defaultdict(lambda: {'F':Counter(),'B':Counter(),'masks':[]})
        state_cells=defaultdict(list); guard_cells=set(); coarse_tau=defaultdict(set)
        for mask in range(total):
            x=presentation(n,mask);y=complement(x);a=canonical(x)[0]
            h=canonical(y,(0,n-1));g=canonical(x,(0,n-1));blue=sigma(x)==x
            key=(g,h,blue);state_cells[key].append(mask)
            coarse_tau[(a,h,blue)].add(g)
            # Complete labelled legal re-marking guards; no floating hashes.
            guard=tuple(sum(v*(n**i) for i,v in enumerate(p)) for p in paths(x))
            guard_cells.add(guard)
            if is_strong(x):
                records[a]['F'][h]+=1;records[a]['B'][h]+=int(blue);records[a]['masks'].append(mask)
        if first_lift is None:
            for vals in state_cells.values():
                if len(vals)>1:
                    i,j=vals[:2];x=presentation(n,i);y=presentation(n,j)
                    gx=set(paths(x));gy=set(paths(y));p=min(gx^gy)
                    first_lift={'n':n,'masks':[i,j],'rows':[list(x),list(y)],'permutation':list(p),
                                'legal':[p in gx,p in gy],'same_both_joint_rooted_endpoints':True,
                                'blue':sigma(x)==x}
                    break
        if first_tau is None:
            collision=next((k for k,v in coarse_tau.items() if len(v)>1),None)
            if collision is not None:
                ids=[]
                for mask in range(total):
                    x=presentation(n,mask)
                    if (canonical(x)[0],canonical(complement(x),(0,n-1)),sigma(x)==x)==collision:ids.append(mask)
                first_tau={'n':n,'masks':ids,'rooted_source_types':len(coarse_tau[collision])}
        rhist=Counter();primes=0;maxfiber=0;blue_atoms=0
        for a,rec in sorted(records.items()):
            aut=len(automorphisms(a));H=hamilton(a)
            require(H==aut*len(rec['masks']),'H/Aut source mass')
            prime=is_prime(a);primes+=prime;blue_atoms+=any(rec['B'].values())
            hist=Counter(rec['F'].values());rhist.update(hist);maxfiber=max(maxfiber,max(hist))
            response=[{'rows':list(h[0]),'roots':list(h[1]),'total':v,'blue':rec['B'][h]}
                      for h,v in sorted(rec['F'].items())]
            atoms.append({'n':n,'source':list(a),'prime':prime,'H':H,'Aut':aut,
                          'presentations':len(rec['masks']),'response':response})
        summary.append({'n':n,'states':total,'strong_classes':len(records),'prime_classes':primes,
                        'strong_presentations':sum(len(r['masks']) for r in records.values()),
                        'blue_strong_classes':blue_atoms,'joint_response_fiber_hist':dict(sorted(rhist.items())),
                        'max_response_fiber':maxfiber,'both_rooted_endpoint_cells':len(state_cells),
                        'guard_cells':len(guard_cells),'guard_injective':len(guard_cells)==total})
    # All ordered path covers of each strong atom through five. Cuts are mandatory
    # only at missing forward arcs; optional cuts may split existing directed arcs.
    covers=[]
    for atom in atoms:
        n=atom['n']
        if n>5:continue
        t=tuple(atom['source']);byk=defaultdict(Counter);raw=Counter()
        for p in permutations(range(n)):
            missing=sum((not edge(t,p[i],p[i+1]))<<i for i in range(n-1))
            for cuts in range(1<<max(0,n-1)):
                if cuts&missing!=missing:continue
                runs=[];start=0
                for i in range(n-1):
                    if cuts>>i&1:runs.append(p[start:i+1]);start=i+1
                runs.append(p[start:]);u=list(converse(t))
                for run in runs:
                    for a,b in zip(run,run[1:]):u[a]^=1<<b;u[b]^=1<<a
                roots=tuple(v for run in runs for v in (run[0],run[-1]));h=canonical(tuple(u),roots)
                byk[len(runs)][h]+=1;raw[len(runs)]+=1
        response=[]
        for k,hist in sorted(byk.items()):
            require(all(v%atom['Aut']==0 for v in hist.values()),'ordered-cover free action')
            response.append({'k':k,'raw':raw[k],'orbits':raw[k]//atom['Aut'],
                             'states':[{'rows':list(h[0]),'roots':list(h[1]),'multiplicity':v//atom['Aut']}
                                       for h,v in sorted(hist.items())]})
        require({(tuple(h['rows']),tuple(h['roots'])):h['multiplicity'] for h in response[0]['states']}==
                {(tuple(h['rows']),tuple(h['roots'])):h['total'] for h in atom['response']},'k=1 matches atom atlas')
        covers.append({'source':atom['source'],'n':n,'responses':response})
    # A source with a prime quotient and repeated child visits; global endpoints
    # do not determine its target, even before iterating the operation.
    t,owner=substitute((2,4,1),[(2,0)]*3);bins=defaultdict(list);run_hostile=None
    for p in paths(t):
        u=complement(t,p);runs=[]
        for v in p:
            if not runs or owner[runs[-1][-1]]!=owner[v]:runs.append([v])
            else:runs[-1].append(v)
        bins[(p[0],p[-1])].append((p,u,runs))
    for vals in bins.values():
        if run_hostile:break
        for i,a in enumerate(vals):
            for b in vals[i+1:]:
                if canonical(a[1])[0]==canonical(b[1])[0]:continue
                if not any(Counter(owner[r[0]] for r in z[2]).most_common(1)[0][1]>1 for z in [a,b]):continue
                run_hostile={'source':list(t),'paths':[list(a[0]),list(b[0])],'runs':[a[2],b[2]],
                             'targets':[list(a[1]),list(b[1])],'target_H':[hamilton(a[1]),hamilton(b[1])]};break
            if run_hostile:break
    require(run_hostile is not None,'six-vertex run hostile')
    receipt={'universe':summary,'first_paired_port_guard_failure':first_lift,
             'first_one_target_response_tau_failure':first_tau,'six_vertex_run_hostile':run_hostile,
             'all_run_cover_universe':{'max_n':5,'strong_atoms':len(covers),
               'raw_covers':sum(r['raw'] for a in covers for r in a['responses']),
               'response_states':sum(len(r['states']) for a in covers for r in a['responses'])}}
    here=Path(__file__).resolve().parent
    if '--verify-atlases' in sys.argv:
        require((here/'strong-atlas.json').read_bytes()==dump(atoms).encode('utf-8'),'frozen strong atlas')
        require((here/'run-atlas.json').read_bytes()==dump(covers).encode('utf-8'),'frozen all-arity atlas')
    if '--write' in sys.argv:
        (here/'strong-atlas.json').write_text(dump(atoms),encoding='utf-8',newline='\n')
        (here/'run-atlas.json').write_text(dump(covers),encoding='utf-8',newline='\n')
    print(json.dumps(receipt,sort_keys=True,indent=2))

if __name__=='__main__':main()
