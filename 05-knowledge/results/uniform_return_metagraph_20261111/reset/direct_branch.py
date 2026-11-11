"""All positive dyadic module templates extending the retained six-vertex core.

The source is T_a(1) with port b deleted. A complete strong-ear search is
pruned only by universal induced-count monotonicity and positive growth.
"""
from pathlib import Path
import itertools, json, sys, time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'chanpak_five_core_bridges_20261110/alternatives'))
import check as a
import two_birth as b
from classify_modules import H,strong,need

def key(rows,base):
    # Old named vertices fixed; only the newly born names are quotiented.
    return min(a.induced(rows,tuple(range(base))+p)
               for p in itertools.permutations(range(base,len(rows))))

def main():
    started=time.perf_counter()
    T=a.ready_source((2,3,1),1)
    W=a.induced(T,[v for v in range(len(T)) if v!=5])
    bound=221
    base=len(W)
    first={W:[]}
    levels=[]
    witnesses=[]
    while first:
        second={}
        examined=0
        for rows,path in first.items():
            for out,h in enumerate(b.all_birth_counts(rows)):
                examined+=1
                if h>bound:continue
                child=a.attach(rows,out)
                if not strong(child):continue
                co=b.module_runs(child,0)
                if any(co[3:]) or co[2]>102:continue
                canonical=key(child,base)
                if canonical in second:continue
                second[canonical]=path+[out]
                if h==221 and co[2]==102:
                    witnesses.append({'rows':canonical,'co':co})
        levels.append({'order':len(next(iter(first)))+1,'parents':len(first),
                       'cuts':examined,'children':len(second),
                       'targets':len(witnesses)})
        print(levels[-1],flush=True)
        first=second
    report={'source':W,'source_h':H(W),'bound':bound,'module':0,'levels':levels,
            'witnesses':witnesses,'seconds':time.perf_counter()-started}
    (HERE/'direct_branch.json').write_text(json.dumps(report,indent=2)+'\n')
    print('DONE',len(witnesses),flush=True)
if __name__=='__main__':main()
