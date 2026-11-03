"""Complete minimum carries for the supplied endpoints and return catalog."""
from fractions import Fraction as F
import json
from support import *
from inverse66 import decode

def main():
    census=json.loads((HERE/'return_census.json').read_text())
    records=census['returns']+[{'word':['s78'],'matrix':MATS['s78']}]
    unique={}
    for r in sorted(records,key=lambda r:(sum(COSTS[l] for l in r['word']),tuple(r['word']))):
        p=tuple(map(tuple,r['matrix']))
        lower=p[0][0]+min(F(p[1][0]+p[2][0],2),p[2][0])
        need(lower>=5,'all return blocks multiply H by at least five on the physical cone')
        unique.setdefault(p,(tuple(r['word']),lower))
    sources=[vector(w) for w in SEEDS]+list(REPAIRS)
    targets={3*v[0]+2*j:(v,j) for v in sources for j in (1,2)}
    cap=max(targets)//5;states=[(1,0,0)];words={states[0]:()}
    for v in states:
        for l,p in MATS.items():
            z=step(v,p)
            if z[0]>cap or z in words:continue
            words[z]=words[v]+(l,);states.append(z)
    need(cap==177 and len(states)==939,'complete bounded parent universe')
    for v in states:
        decoded=decode(v);need(decoded is not None and vector(decoded)==v,'independent inverse versus forward parent closure')
    hits={h:[] for h in targets};tested=0
    for w in states[1:]:
        for p,(r,lower) in unique.items():
            if lower*w[0]>max(targets):continue
            tested+=1;u=step(w,p)
            if u[0] not in targets:continue
            v,j=targets[u[0]];dx=u[1]-(3*v[1]-j);dy=u[2]-3*v[2]
            if dx%7 or dy%7:continue
            cx,cy=dx//7,dy//7;full=words[w]+r
            hits[u[0]].append((abs(cx)+abs(cy),sum(COSTS[l] for l in full),full,cx,cy,w,r,words[w],u))
    results=[]
    for h,(v,j) in targets.items():
        options=sorted(hits[h]);answer=None
        if options:
            radius,n,full,cx,cy,w,r,parent,u=options[0]
            rows=actual_target(full);count,width,peak=literal_frontier(rows)
            need(count==h and vector(full)==u,'literal minimum-carry source')
            answer={'radius':radius,'carry':[0,cx,cy],'output_state':u,'parent_state':w,
                'parent_word':parent,'return_word':r,'output_word':full,'vertices':n,'port':n-1,
                'target_rows':rows,'literal_H':count,'width':width,'peak':peak,
                'candidates_in_minimal_layer':sum(o[0]==radius for o in options)}
        results.append({'input_state':v,'digit':j,'H_target':h,'total_candidates':len(options),'answer':answer})
        print(v,j,answer and (answer['radius'],answer['return_word'],answer['output_state']),flush=True)
    result={'status':'PASS','parent_H_bound':cap,'parent_states_including_empty':len(states),
        'return_words':len(records),'distinct_return_matrices':len(unique),'tested_state_matrix_pairs':tested,
        'scope':'Complete all-word, all-carry test: parent and intermediate words in the 66 bank; nonempty return in the catalog monoid. H is fixed and carry is 7(0,cx,cy). No carry cutoff.',
        'results':results,'parent_certificates':[{'state':v,'word':words[v]} for v in states]}
    (HERE/'carry_complete.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
if __name__=='__main__':main()
