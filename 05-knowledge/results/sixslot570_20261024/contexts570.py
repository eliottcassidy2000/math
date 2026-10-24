"""Compare actual two-target context predicates on ALL original source classes."""
from common import *
from collections import Counter,defaultdict
def main():
    nodes,owner,merged,X,Y=data();old={r['source'] for r in inventory()};T=(570,1965,2234)
    totals=Counter();best={};examples={};all_pairs=[]
    for s,row in enumerate(nodes):
        if s in T:continue
        B=matrix(row['canon_a']);D=defaultdict(list)
        for p in paths(B):
            w=word(B,p);t=merged[owner[FULL^w]]
            if t in T:D[t].append((p,w))
        for t in T[1:]:
            if not D[570] or not D[t]:continue
            totals['common_source_pairs']+=1
            for p,w in D[570]:
                for q,z in D[t]:
                    equal=[i for i in range(8) if p[i]==q[i]]
                    cost=8-len(equal);ends=(p[0]==q[0] and p[-1]==q[-1]);arcs=[j for j in range(1,6) if p[j:j+2]==q[j:j+2]]
                    kinds=['unrestricted']
                    if ends:kinds.append('both_endpoints')
                    if arcs:kinds.append('indexed_internal_arc')
                    if ends and arcs:kinds.append('endpoints_and_indexed_arc')
                    if p[0]==q[0]:kinds.append('first_endpoint')
                    if p[-1]==q[-1]:kinds.append('last_endpoint')
                    for kind in kinds:
                        totals[kind]+=1;key=(cost,s,t,p,q)
                        if kind not in best or key<best[kind]:
                            best[kind]=key;examples[kind]={'source':s,'target':t,'cost':cost,'before':{'path':p,'word':w,'target':570},'after':{'path':q,'word':z,'target':t},'fixed_clocks':equal,'indexed_arcs':arcs,'existing_source':s in old,'source_label':label(X,Y,s)}
                    all_pairs.append({'source':s,'target':t,'cost':cost,'before_path':p,'after_path':q,'fixed':equal})
    return {'scope':'all 3528 canonical merged sources outside three target rows, including selfconverse and nontrivial automorphism sources; all literal Hamiltonian paths in canonical ordinary representative','counts':dict(totals),'best':examples,'pairs':all_pairs,'input_revision':REV,'inputs':INPUTS}
if __name__=='__main__':
    out=main();(HERE/'contexts570.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'counts':out['counts'],'best':out['best']},indent=2))
