#!/usr/bin/env python3
"""Joint shared-dictionary and physical-port selection at the declared row."""
from hashlib import sha256
import importlib.util
from itertools import product
from pathlib import Path
import json
from selector import make_tree,nodes,build_table,query,select,cover,valid,leaves,require

def main():
    base=Path(__file__).resolve().parent
    source=base.parent/'history/check_shared_dictionary.py'
    spec=importlib.util.spec_from_file_location('owned_dictionary',source)
    dictionary=importlib.util.module_from_spec(spec); spec.loader.exec_module(dictionary)
    tags=dictionary.tagged_blocks(7)[7][:9]
    t=make_tree(list(enumerate((0,0,0,1,1,1,0,0,1))))
    E=sum(v.active for v in nodes(t)); entries=[]
    for deleted in (False,True):
        free=sorted({owner for _,owner in tags if owner and (not deleted or owner%2==0)})
        results=[]; physical_sets=set()
        for bits in product((0,1),repeat=len(free)):
            assignment=dict(zip(free,bits)); prefix=[]
            for color,owner in tags:
                if assignment.get(owner,0): color='B' if owner%2 else 'A'
                prefix.append(color)
            row=prefix[::-1]; roles={i:'CAB'.index(c) for i,c in enumerate(row)}
            table=build_table(t,511,3,roles,(1,1,1),True)
            total=sum(n for _,n in query(table,3,E,(1,1,1)))
            actual={select(table,3,E,r,(1,1,1)) for r in range(total)}
            direct={S for S in range(512) if S.bit_count()==3 and valid(S,511,True)
                    and cover(t,S)==E and all(sum(roles[i]==c for i in leaves(S))==1 for c in range(3))}
            require(actual==direct,'joint dictionary/marking lift mismatch')
            physical_sets.update(actual)
            results.append({'owner_bits':assignment,'row':''.join(row),'solutions':total,
                            'rank0':list(leaves(select(table,3,E,0,(1,1,1)))) if total else None})
        entries.append({'grammar':'D' if deleted else 'R','visible_free_owners':free,
                        'dictionary_assignments':len(results),
                        'dictionary_marking_pairs':sum(x['solutions'] for x in results),
                        'distinct_physical_markings':len(physical_sets),'assignments':results})
    receipt={'status':'FINITE-EXACT control of the joint finite selector theorem',
             'row_length':9,'canonical_tags':tags,'physical_palette':'BBB,DDD,BBD',
             'target':{'marks':3,'nonadjacent':True,'role_quota_CAB':[1,1,1],'residual_group_order':1},
             'owner_source':'../history/check_shared_dictionary.py',
             'owner_source_sha256':sha256(source.read_bytes()).hexdigest(),'entries':entries,
             'scope':'Pairs of consistent dictionaries and actual marking sets; not independent role choices or Hamiltonian-path targets.',
             'reproduce':'python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/charts/joint_selector.py'}
    (base/'joint_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in x.items() if k!='assignments'} for x in entries],indent=2))

if __name__=='__main__': main()
