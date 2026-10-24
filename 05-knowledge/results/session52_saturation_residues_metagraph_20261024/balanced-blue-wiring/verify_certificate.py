#!/usr/bin/env python3
"""Replay literal records and adversarial controls without regenerating the atlas."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from collections import Counter
import copy,importlib.util,json,hashlib,argparse
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('producer',HERE/'check.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def validate(d):
    nodes={r['atlas_index']:r for r in d['nodes']};targets={r['boundary_class_id']:r for r in d['boundary_target_classes']}
    need(len(nodes)==229 and len(targets)==786,'node count')
    ports={r[0]:r for r in d['ports']}
    need(len(ports)==len(d['ports'])==41699,'unique complete port records')
    adj=Counter();badj=Counter();degree=Counter();boundary=Counter();blue=Counter()
    for word,idx,path,cword,kind,tid,order in d['ports']:
        need(sorted(path)==list(range(10)),'source path permutation')
        T=tuple(nodes[idx]['canonical_rows']);U=m.normalize(T,path)
        need(all(U[i]>>(i+1)&1 for i in range(9)),'Hamiltonian source path')
        need(m.decode(word)==U,'source literal matrix')
        need(all((U[i]>>j&1)==(U[9-j]>>(9-i)&1) for i in range(10) for j in range(10)),'blue reflection')
        need(cword==m.FULL^word,'literal native C')
        need(sorted(order)==list(range(10)),'target permutation')
        V=m.decode(cword);blue[idx]+=1
        if kind=='internal':
            need(cword in ports and ports[cword][1]==tid and ports[cword][3]==word,'reciprocal internal record')
            need(V==m.normalize(tuple(nodes[tid]['canonical_rows']),order),'internal witness')
            degree[idx]+=1
            if word<cword:adj[tuple(sorted((idx,tid)))]+=1
        else:
            need(kind=='boundary' and cword not in ports,'boundary classification')
            need(m.normalize(V,order)==tuple(targets[tid]['canonical_rows']),'boundary canonical witness')
            need(sorted(r.bit_count() for r in V)==[3]+[4]*4+[5]*4+[6],'boundary scores')
            badj[idx,tid]+=1;boundary[idx]+=1
    need([[a,b,v] for (a,b),v in sorted(adj.items())]==d['adjacency_multiplicities'],'exact internal multiplicities')
    need([[a,b,v] for (a,b),v in sorted(badj.items())]==d['boundary_adjacency'],'exact boundary multiplicities')
    for i,r in nodes.items():need((degree[i],boundary[i],blue[i])==(r['internal_degree'],r['boundary_edges'],r['blue_presentation_words']),'node conservation')
    for key,value in m.graph_statistics(nodes,adj).items():need(d[key]==json.loads(json.dumps(value)),'graph statistic '+key)
    need(sum(adj.values())==d['internal_blue_edges']==11452,'edge count')
    need(sum(v for (a,b),v in adj.items() if a==b)==d['internal_blue_loops']==24,'loop count')
    need(sum(badj.values())==d['boundary_blue_edges']==18795,'boundary edge count')
    return {'ports_replayed':len(ports),'internal_edges_replayed':sum(adj.values()),'boundary_edges_replayed':sum(badj.values())}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default=str(HERE/'normal.json'));ap.add_argument('--output',default=str(HERE/'certificate-check.json'));args=ap.parse_args()
    d=json.loads(Path(args.input).read_text());result=validate(d)
    # Targeted corruptions stop early; no dependence on Python assert.
    rejected=[]
    for name,mutate in [
      ('wrong_native_C',lambda x:x['ports'][0].__setitem__(3,x['ports'][0][3]^1)),
      ('nonpermutation_source_path',lambda x:x['ports'][0][2].__setitem__(0,x['ports'][0][2][1])),
      ('duplicated_port',lambda x:x['ports'].__setitem__(1,x['ports'][0])),
      ('wrong_target_witness',lambda x:x['ports'][0][6].__setitem__(0,x['ports'][0][6][1]))]:
        damaged=copy.deepcopy(d);mutate(damaged)
        try:validate(damaged)
        except RuntimeError:rejected.append(name)
        else:raise RuntimeError('mutation accepted: '+name)
    fixtures=[('single_bridge',[0,1],{(0,1):1},1,[(0,1)],[]),('parallel_not_bridge',[0,1],{(0,1):2},1,[],[]),('loop_and_isolate',[0,1],{(0,0):1},2,[],[]),('path_cut',[0,1,2],{(0,1):1,(1,2):1},1,[(0,1),(1,2)],[1]),('triangle',[0,1,2],{(0,1):1,(1,2):1,(0,2):1},1,[],[])]
    for name,vertices,edges,cc,bridges,cuts in fixtures:
        g=m.graph_statistics(vertices,Counter(edges));need((g['connected_component_count'],g['bridge_pairs'],g['articulation_nodes'])==(cc,bridges,cuts),'hostile topology '+name)
    result.update(status='PASS',rejected_mutations=rejected,topology_fixtures=[r[0] for r in fixtures],certificate_sha256=hashlib.sha256(Path(args.input).read_bytes()).hexdigest())
    Path(args.output).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
