"""Complete old pins at three targets; frozen native matrices and literal paths."""
from common import *
from collections import Counter, defaultdict

SUPPLY_REV='31c8e07bab98ee047edf8f3537b36ef3d7bca11f'
def supplied(path):
    raw=subprocess.check_output(['git','show',SUPPLY_REV+':05-knowledge/results/'+path],cwd=ROOT)
    INPUTS[path]=hashlib.sha256(raw).hexdigest()
    return json.loads(raw)

def main():
    nodes,owner,merged,X,Y=data()
    supply=supplied('compatible_charged_helpers_20261023/supply/selected_supply.json')
    fresh=supplied('sixslot570_20261024/native570.normal.json')
    selected=[]
    for r in inventory():selected.append((r['source'],'core',r['before']['path'],r['after']['path']))
    for r in supply['helpers']:selected.append((r['source'],'helper',r['selected']['before']['path'],r['selected']['after']['path']))
    for r in fresh['literal_sources']:selected.append((r['source'],'fresh',r['domain'][0]['path'],r['domain'][1]['path']))
    need(len(set(s for s,_,_,_ in selected))==85,'actual source identities')
    rows=[];resources=defaultdict(list)
    for s,role,p,q in sorted(selected):
        B=matrix(nodes[s]['canon_a']);pins=[(i,p[i]) for i in range(8) if p[i]==q[i]]
        raw=[];options=[]
        for h in paths(B):
            if not all(h[i]==v for i,v in pins):continue
            w=word(B,h);t=merged[owner[FULL^w]];raw.append(t)
            if t not in (570,1965,2234):continue
            need(merged[owner[w]]==s,'literal source')
            O=orbit(w);need(len(O)==4,'closed resource size')
            option={'path':h,'target':t,'word':w,'resource':min(O),'closed_words':O}
            options.append(option);resources[min(O)].append((s,h,t))
        need(tuple(p) in [o['path'] for o in options] and tuple(q) in [o['path'] for o in options],'inherited records retained')
        rows.append({'source':s,'role':role,'pins':pins,'label':label(X,Y,s),'X':X[s],
                     'raw_domain_size':len(raw),'raw_targets':dict(Counter(raw)), 'options':options})
    collisions=[v for v in resources.values() if len(v)>1]
    need(not collisions,'all literal eligible records occupy distinct reserved units')
    return {'status':'FINITE-EXACT complete 85-source pin-defined three-target menus',
            'targets':[570,1965,2234],'source_count':len(rows),'option_count':sum(len(r['options']) for r in rows),
            'three_way_sources':[r['source'] for r in rows if len(set(o['target'] for o in r['options']))==3],
            'target570_sources':[r['source'] for r in rows if any(o['target']==570 for o in r['options'])],
            'resource_collisions':collisions,'rows':rows,'inputs':INPUTS,'source_revision':REV,'supply_revision':SUPPLY_REV}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
