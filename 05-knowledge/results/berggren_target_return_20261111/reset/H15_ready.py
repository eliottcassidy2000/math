"""Exact rooted retention optima onto the four strong two-ready 5993 classes."""
from pathlib import Path
from itertools import combinations
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'uniform_return_metagraph_20261111/reset'))
from classify_modules import H,strong,need,inherited as a
from probe_ready import embeddings
from finite_bridge import replace

def main():
    data=json.loads((HERE/'H15_inverse.json').read_text())
    S=next(tuple(h['rows']) for row in data['cores'] if row['name']=='S4'
           for h in row['hits'] if (h['cut'],h['A'],h['B'])==(1,3,57))
    candidates=json.loads((HERE/'ready_two_inverse.json').read_text())['hits']['5993']
    need(len(candidates)==4 and len(S)==13,'exact declared finite search universe')
    out=[]
    for i,rec in enumerate(candidates):
        T=a.birth_beaten_only_by(rec['parent'],rec['root']);r=len(T)-1;p=rec['root']
        need(H(T)==5993 and strong(T) and a.root_profiles(T)[r]['gap']==2,'literal two-ready target')
        failures=[];answer=None
        for loss in range(max(0,len(S)-len(T)),len(S)-1):
            tested=0
            for drop in combinations(range(1,12),loss):
                kept=[v for v in range(len(S)) if v not in drop];tested+=1
                mapping=next(embeddings(a.induced(S,kept),T,{0:r,len(kept)-1:p}),None)
                if mapping is not None:
                    answer=(kept,list(mapping));break
            if answer:break
            failures.append(dict(loss=loss,all_source_subsets=tested))
        need(answer is not None,'root and sole incoming neighbor always yield at least a two-vertex map')
        kept,mapping=answer
        U,names,rewrite=replace(S,[f'old{v}' for v in range(13)],T,kept,mapping,
                              [f'new{v}' for v in range(len(T)) if v not in mapping])
        need(names[r]=='old0' and names[p]=='old12','same physical root and sole in-neighbor survive')
        nextT=a.birth_beaten_only_by(a.birth_beaten_only_by(T,r),len(T))
        need(H(nextT)==17983,'next lawful plus-four target')
        row=dict(class_index=i,order=len(T),root=r,port=p,rows=T,
                 retained=len(kept),source_loss=loss,failed_higher_retention_levels=failures,
                 witness=rewrite,next_target=17983)
        out.append(row);print('CLASS',i,'ORDER',len(T),'RETAIN',len(kept),'TRACE',rewrite['trace'],flush=True)
    result=dict(status='PASS',scope='Four strong two-ready 5993 rooted classes, fixed source S4 chart(1,3,57), fixed root0 and port12; exhaustive failed source subsets and exact injective adjacency maps at every higher retention level',source_rows=S,targets=out)
    (HERE/'H15_ready.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: exact four-class rooted retention optima',flush=True)

if __name__=='__main__':main()
