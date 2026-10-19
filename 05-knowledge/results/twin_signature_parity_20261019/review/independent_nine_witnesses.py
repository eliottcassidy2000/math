from independent_embedding_witnesses import count
from independent_tree_capacity import canon
from independent_cherry_extensions import linear_extensions,arcs
from pathlib import Path
import json,hashlib

def require(x,m):
    if not x:raise RuntimeError(m)

def main():
    here=Path(__file__).parent
    source=here/'nine_vertex_bad_witnesses.json'
    data=json.loads(source.read_text());prior=json.loads((here/'nine_criterion_receipt.json').read_text())
    expected={x['canonical'] for x in prior['excluded_classes']};observed=set();rows=[]
    for w in data['witnesses']:
        e=tuple(map(tuple,w['edges']));m=w['orientation_mask'];h=w['host_out_bits']
        require(len(h)==9 and all(0<=x<512 for x in h),'host shape')
        require(all(not(h[i]>>i&1) for i in range(9)),'host loops')
        require(all((h[i]>>j&1)+(h[j]>>i&1)==1 for i in range(9) for j in range(i+1,9)),'host tournament')
        a=linear_extensions(9,arcs(e,m));b=count(e,m,h)
        require(a==w['transitive_count'] and b==w['host_count'] and (a^b)&1,'witness exact count')
        observed.add(canon(9,e));rows.append([a,b])
    require(observed==expected,'complete excluded class coverage')
    print(json.dumps({'status':'PASS','all_11_witness_counts':rows,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'scope':'Every supplied arbitrary tournament count independently recomputed by all 9! permutations; transitive count independently recomputed by subset linear-extension DP. Complete source-class coverage checked against independent criterion census.'},indent=2))

if __name__=='__main__':main()
