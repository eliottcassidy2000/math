"""Exhaust frozen nonselfconverse trivial-aut sources for fresh 570 pins."""
from common import *
from collections import defaultdict,Counter
def main():
    nodes,owner,merged,X,Y=data();old={r['source'] for r in inventory()};T=(570,1965,2234)
    counts=Counter();records=[];triples=[]
    for s,row in enumerate(nodes):
        if row['class_a']==row['class_b'] or s in T:continue
        B=matrix(row['canon_a']);P=list(paths(B))
        if 2*len(P)!=int(row['fibre']):continue
        counts['eligible_original_sources']+=1
        menus=defaultdict(lambda:defaultdict(list));caps=Counter()
        for p in P:
            w=word(B,p);t=merged[owner[FULL^w]]
            if t not in T:continue
            need(merged[owner[w]]==s,'actual source')
            caps[t]+=1
            for j in range(1,6):menus[(p[0],p[-1],j,p[j],p[j+1])][t].append({'path':p,'word':w,'target':t})
        best={}
        for pins,D in menus.items():
            if 570 not in D:continue
            if all(t in D for t in T):triples.append({'source':s,'pins':pins,'domains':D})
            for t in T[1:]:
                for p in D[570]:
                    for q in D[t]:
                        key=(sum(a!=b for a,b in zip(p['path'],q['path'])),pins,p['path'],q['path'])
                        if t not in best or key<best[t][0]:best[t]=(key,p,q)
        for t,(key,p,q) in best.items():
            records.append({'source':s,'existing_source':s in old,'source_label':label(X,Y,s),'cost':key[0],'pins':key[1],
                            'before':p,'after':q,'full_target_path_capacities':caps})
    return {'scope':'all nonselfconverse original canonical sources, trivial automorphisms, outside the three targets; fixed endpoints and at least one indexed internal arc',
            'target_labels':{str(t):label(X,Y,t) for t in T},'target_X':{str(t):X[t] for t in T},
            'counts':dict(counts),'pair_records':records,'three_target_pin_contexts':triples,'input_revision':REV,'inputs':INPUTS}
if __name__=='__main__':
    out=main();(HERE/'seek570.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'pair_records':len(out['pair_records']),'fresh_sources':sorted({r['source'] for r in out['pair_records'] if not r['existing_source']}),'three_target_pin_contexts':len(out['three_target_pin_contexts']),'target_X':out['target_X']},indent=2))
