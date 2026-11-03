"""Guard-prioritized actual native/birth searches in the retained old frame."""
import sys,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'sheet_return_compiler_20261103'))
from next_escapes import profile
from bulk import tables,native,birth,strong,need
from escape_lifts import ear_values
from return_program import opposite


def paths(T,priority):
    starts=tables(opposite(T))[1]
    def go(prefix,left):
        if not left:yield prefix;return
        for v in priority:
            if left>>v&1 and starts[left][v] and (not prefix or T[prefix[-1]]>>v&1):
                yield from go(prefix+(v,),left^(1<<v))
    return go((),(1<<len(T))-1)


def bank():
    old=json.loads((HERE.parent/'sheet_return_compiler_20261103/next_escapes.json').read_text())
    return {int(h):tuple(w['target']) for row in old for h,w in row['witnesses'].items()}


def search(T):
    initial=profile(T);h=initial['observer'][0];wanted={3*h+2,3*h+4}
    # Old names survive; only enumeration priority changes with the exact state.
    priority=tuple(sorted(range(len(T)),key=lambda v:(initial['root_defects'][v]!=0,initial['root_defects'][v],v)))
    found={};tested=0
    print('begin',h,'priority',priority,flush=True)
    def programs():
        yield None,T
        for p in paths(T,priority):yield p,native.native(T,p)
    for p,U in programs():
        tested+=1;hi=tables(U)[0][-1]
        if hi<=max(wanted):
            for mask,value in enumerate(ear_values(U)):
                if value in wanted-found.keys():
                    V=birth(U,mask)
                    if strong(V):
                        found[value]=dict(path=p,intermediate=U,intermediate_H=hi,mask=mask,target=V,profile=profile(V))
                        print(h,'->',value,'program',tested,'intermediate',hi,flush=True)
        tables.cache_clear()
        if set(found)==wanted:break
        if tested%250==0:print('searched',h,tested,'of',h+1,flush=True)
    return dict(source=T,H=h,profile=initial,path_priority=priority,programs_examined=tested,witnesses=found)


def main():
    rows=[]
    for h,T in sorted(bank().items()):
        rows.append(search(T));(HERE/'branches.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


if __name__=='__main__':main()
