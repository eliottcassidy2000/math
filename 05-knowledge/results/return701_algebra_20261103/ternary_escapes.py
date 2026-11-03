"""Changed-relation ternary outputs from the actual supplied701 source."""
from collections import Counter
import sys
from common import HERE,ROOT,need,junction,save,source701,values
from return_native import paths,phases
sys.path.insert(0,str(ROOT/'05-knowledge/results/defect_offset_dynamics_20261103'))
from ear_cuts import ear_values,tables

def main():
    U,r,b,names=source701();targets=(2103,2105,2107);counts=Counter();selected={};all_hits=[];tested=0
    for number,P in enumerate(paths(U)):
        V=junction.native(U,P);h=junction.H(V)
        if h>max(targets):continue
        vals=ear_values(V);tested+=len(vals)
        for mask,hh in enumerate(vals):
            if hh not in targets:continue
            W=junction.ear(V,mask)
            if not junction.strong(W):continue
            need(junction.H(W)==hh,'independent literal ternary target')
            prof=junction.profile(W,r);counts[hh]+=1
            row={'target':hh,'path':P,'intermediate_H':h,'birth_mask':mask,'retained_root_profile':prof}
            all_hits.append(row)
            score=(prof[3]-hh,h,number,mask)
            if hh not in selected or score<selected[hh]['score']:
                selected[hh]={'score':score,**row,'rows':W,'source_names':names+[15],'root':r,'spare':b}
        tables.cache_clear()
    for row in selected.values():
        T=tuple(row['rows']);row['root_phases']=phases(T,r);row['spare_phases']=phases(T,b)
    save('ternary_escapes.json',{'status':'FINITE_EXACT_NATIVE_TERNARY_BRANCHES','source_rows':U,'source_names':names,
       'root':r,'spare':b,'evaluated_birth_masks':tested,'counts':dict(counts),'selected':selected,'all_hits':all_hits})
    print('ternary escape counts',dict(counts),'selected profiles',{h:x['retained_root_profile'] for h,x in selected.items()},flush=True)
if __name__=='__main__':main()
