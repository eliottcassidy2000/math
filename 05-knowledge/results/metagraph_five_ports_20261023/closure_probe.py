"""Reflection-implied converse margins: forests and actual singleton hostiles.

Uses the incoming Session38 forest theorem. It does not confuse the joint
pair of signature labels with the two separate margin families.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,gzip,hashlib,json,struct

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];checks=0
def need(x,why):
    global checks
    checks+=1
    if not x:raise RuntimeError(why)
def inventory(n):
    if n==9:
        p=ROOT/'05-knowledge/results/blue_quotient_20261020b/n9_classes.tsv.gz'
        rows=list(csv.DictReader(gzip.open(p,'rt'),delimiter='\t'))
        dual={int(r['code']):int(r['converse_code']) for r in rows}
    else:
        p=ROOT/'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv'
        rows=list(csv.DictReader(p.open(),delimiter='\t'))
        dual={int(r['canon']):int(rows[int(r['sigmaClass'])]['canon']) for r in rows}
    return dual,hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    probe=json.loads((HERE/'probe.json').read_text());types=probe['rooted_type_representatives'];out={'status':'FINITE-EXACT singleton raw-rooted-five coordinate menus under native converse margins','orders':{}}
    n9=None
    for n in (8,9):
        menu=json.loads((HERE/f'menu_n{n}.json').read_text());base=next(x for x in menu['menus'] if x['name']=='v4' and x['round']==1);groups=base['collisions'];gid={c:i for i,g in enumerate(groups) for c in g}
        rec={r['code']:r['rooted5'] for r in probe['orders'][str(n)]['records']};dual,h=inventory(n);need(set(gid)==set(rec) and all(dual[c] in gid for c in gid),'closed complete residual bank')
        menus=[];allkeys=[]
        for k,word in enumerate(types):
            key={c:(gid[c],tuple(sorted(row[k] for row in rec[c]))) for c in gid};allkeys.append(key)
            parent={};byedge=defaultdict(list);cycles=0
            def find(v):
                parent.setdefault(v,v)
                while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
                return v
            for c in sorted(gid):
                left=(0,key[c]);right=(1,key[dual[c]]);a=find(left);b=find(right)
                if a==b:cycles+=1
                else:parent[b]=a
                byedge[key[c],key[dual[c]]].append(c)
            parallel=sorted(sorted(g) for g in byedge.values() if len(g)>1);components=len({find(v) for v in parent});need(cycles==len(gid)-len(parent)+components,'multigraph cycle rank')
            row={'index':k,'rooted_word':word,'class_collisions':len(gid)-len(set(key.values())),'incidence_vertices':len(parent),'incidence_edges':len(gid),'components':components,'cycle_rank':cycles,'forest':cycles==0,'parallel_class_groups':parallel}
            if parallel:
                c,d=parallel[0][:2];need(key[c]==key[d] and key[dual[c]]==key[dual[d]],'parallel witness')
                row['two_cycle_witness']={'classes':[c,d],'converse_classes':[dual[c],dual[d]],'left_signature':[key[c][0],list(key[c][1])],'right_signature':[key[dual[c]][0],list(key[dual[c]][1])]}
            menus.append(row)
        forests=[r['rooted_word'] for r in menus if r['forest']];repairs=[r['rooted_word'] for r in menus if r['forest'] and r['class_collisions']]
        out['orders'][str(n)]={'class_universe_sha256':h,'complete_classes':len(dual),'residual_classes':len(gid),'forest_words':forests,'noninjective_forest_words':repairs,'menus':menus}
        print('order',n,'forests',len(forests),'noninjective repairs',repairs,flush=True)
        if n==9:n9=gid,rec,dual,allkeys,menus
    need(out['orders']['8']['noninjective_forest_words']==[65,73,198,203,206,219,459,475],'n8 reflection repair changed')
    need(not out['orders']['9']['forest_words'] and all(r['parallel_class_groups'] for r in out['orders']['9']['menus']),'n9 singleton forest boundary')
    gid,rec,dual,keys,menus=n9;raw=gzip.decompress((HERE/'focused_ports.bin.gz').read_bytes());need(hashlib.sha256(raw).hexdigest()==json.loads((HERE/'verification.json').read_text())['port_bytes_sha256'],'native input hash')
    cells=[(i,j) for i in range(9) for j in range(i+2,9)];ix={e:k for k,e in enumerate(cells)};leftmask=sum(1<<ix[0,j] for j in range(2,9));rightmask=sum(1<<ix[i,8] for i in range(7));flat=Counter();examples={};bytarget=defaultdict(lambda:defaultdict(set))
    for mask,s,t,order,path in struct.iter_unpack('<IQQQQ',raw):
        ref=sum(((mask>>k)&1)<<ix[8-j,8-i] for k,(i,j) in enumerate(cells));color=int(mask==ref);a=1+(mask&leftmask).bit_count();b=7-(mask&rightmask).bit_count();e=s,a,b,color,t;flat[e]+=1;examples.setdefault(e,mask)
        if not color:bytarget[t][a,b,color].add(s)
    for (s,a,b,color,t),amount in flat.items():need(flat[dual[s],8-b,8-a,color,dual[t]]==amount,'native reflected capacity')
    two={c:tuple(sorted((r[6],r[24]) for r in rec[c])) for c in gid};switches=[]
    for k,menu in enumerate(menus):
        signature=keys[k];found=None
        for group in menu['parallel_class_groups']:
            if found:break
            for i,c in enumerate(group):
                if found:break
                for d in group[i+1:]:
                    if found:break
                    for typ in sorted(set(bytarget[c])&set(bytarget[d])):
                        if found:break
                        a,b,color=typ
                        for s in sorted(bytarget[c][typ]):
                            if found:break
                            for t in sorted(bytarget[d][typ]):
                                if min(s,dual[s])==min(t,dual[t]):continue
                                delta=Counter({(s,a,b,color,c):-1,(s,a,b,color,d):1,(t,a,b,color,c):1,(t,a,b,color,d):-1})
                                for (u,x,y,z,v),value in list(delta.items()):delta[dual[u],8-y,8-x,z,dual[v]]+=value
                                if any(flat[e]+value<0 for e,value in delta.items()):continue
                                rd=Counter();cd=Counter();merged=Counter();new=Counter()
                                for (u,x,y,z,v),value in delta.items():
                                    rd[u,x,y,z,signature[v]]+=value;cd[v,9-x,7-y,z]+=value;merged[tuple(sorted((min(u,dual[u]),min(v,dual[v]))))]+=value;new[u,x,y,z,gid[v],two[v]]+=value
                                need(not any(rd.values()) and not any(cd.values()),'singleton counterfeit margins')
                                need(all(value%2==0 for value in merged.values()),'black-pair parity')
                                if not any(merged.values()):continue
                                need(any(new.values()),'two-coordinate refinement did not reject switch')
                                for (u,x,y,z,v),value in delta.items():need(delta[dual[u],8-y,8-x,z,dual[v]]==value,'switch reflection equality')
                                found={'coordinate_index':k,'rooted_word':types[k],'parallel_targets':[c,d],'converse_targets':[dual[c],dual[d]],
                                    'edits':[{'source':u,'scores':[x,y],'color':z,'target':v,'delta':value,'native_count':flat[u,x,y,z,v],'native_mask':examples.get((u,x,y,z,v))} for (u,x,y,z,v),value in sorted(delta.items()) if value],
                                    'merged_black_delta':[[list(e),value] for e,value in sorted(merged.items()) if value],
                                    'two_coordinate_rejection_buckets':sum(bool(value) for value in new.values())};break
        need(found is not None,('no grounded singleton counterfeit',k));switches.append(found)
    out['n9_grounded_singleton_counterfeits']=switches;out['native_port_bytes_sha256']=hashlib.sha256(raw).hexdigest();out['scope']='class-signature/source-endpoint/color margins plus explicit reflection equality and target endpoint capacities; not an alternative native realization';out['checks']=checks
    if args.write:(HERE/'closure_probe.json').write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads((HERE/'closure_probe.json').read_text())==out,'frozen closure receipt')
    print('grounded singleton hostiles',len(switches),'checks',checks,flush=True)
if __name__=='__main__':main()
