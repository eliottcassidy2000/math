"""An actual reusable1-ready chart, with absolute-count change explicit."""
from pathlib import Path
import json
from classify_modules import H,strong,need,inherited as a,two_birth as b,converse
from ready_rebase import rebase,full_state
HERE=Path(__file__).resolve().parent
BASE=(54,20,48,55,32,578,159,319,127,479)


def family(k):return b.inflate(BASE,0,k)


def acquire(k):
    R=rebase(k)[-1];T=a.two_birth_reset(R,len(R)-1)
    P=a.induced(T,[v for v in range(len(T)) if v!=k+2])
    return a.birth_beaten_only_by(P,k+4)


def cycle(k):
    Q=family(k);first=a.birth_beaten_only_by(Q,k+8)
    second=a.birth_beaten_only_by(first,len(first)-1)
    kept=[v for v in range(len(second)) if v!=k+5]
    P=a.induced(second,kept);last=a.birth_beaten_only_by(P,k+4)
    mapping=list(range(2,k+2))+[k+2,k+3,1,k+5,k+6,k+4,0,k+7,k+8,k+9,k+10]
    return first,second,kept,P,last,mapping


def formula(k):
    D=1<<k
    return [210*D+11,3,42*D+14,966*3**k+720*D-45,210*D+12,
            42*D+15,966*3**k+768*D-30]


def main():
    P=a.induced(BASE,range(9));z=5;n=len(P)
    separator=tuple(row|(1<<n) for row in P)+((1<<n)-1,)
    forced_start=tuple(row&~(1<<z) for row in separator)
    forced_end=tuple(0 if i==z else row for i,row in enumerate(P))
    runs={key:b.module_runs(T)[1:] for key,T in
          [('parent',P),('child',BASE),('two_covers',separator),
           ('parent_start_covers',forced_start),('parent_end_paths',forced_end)]}
    expected={'parent':[99,42],'child':[431,210],'two_covers':[4404,3666,966],
              'parent_start_covers':[111,48],'parent_end_paths':[3]}
    for key,head in expected.items():
        need(runs[key][:len(head)]==head and not any(runs[key][len(head):]),
             'complete positive marked coefficients: '+key)
    literal=[]
    for k in range(1,6):
        Q=family(k);need(Q==acquire(k),'actual acquisition from earlier425 rebase family')
        need(strong(Q) and full_state(Q,len(Q)-1)==formula(k),'independent full rooted response')
        P=a.induced(Q,range(len(Q)-1));path=[k+4,k+5,k+6,k+7,k+2,*range(k),k,k+1,k+3]
        need(len(set(path))==len(P) and all(P[x]>>y&1 for x,y in zip(path,path[1:]))
             and P[path[-1]]>>path[0]&1,'spanning unique-start cycle')
        need(a.root_profiles(P)[k+4]['s']==1,'literal unique-start parent')
        first,second,kept,card,last,mapping=cycle(k);D=1<<k
        trace=[H(T) for T in(Q,first,second,card,last)]
        need(trace==[210*D+11,378*D+9,630*D+35,168*D+15,840*D+11],
             'independent physical cycle count trace')
        target=family(k+2)
        need(a.induced(target,mapping)==last,'actual all-block pointed isomorphism')
        need(mapping[-1]==len(target)-1,'new marked root is preserved by isomorphism')
        need(a.induced(last,range(len(card)))==a.induced(second,kept),'all persistent arcs')
        need((H(Q)-1)%4==2,'square-plus-one exclusion')
        need((3*H(Q)+2)==35*(18*D+1),'fixed composite ternary excursion')
        for side,T in [('forward',Q),('converse',converse(Q))]:
            ready=[(p['root'],p['gap']) for p in a.root_profiles(T) if p['gap'] in(1,2)]
            need(ready==([(k+8,1)] if side=='forward' else []),'all marked/converse ready controls')
        literal.append({'k':k,'rows':Q,'state':formula(k),'trace':trace,
                        'kept_labels_after_two_births':kept,'return_rows':last,'mapping_to_Q_k_plus_2':mapping})
    iteration=[]
    for start in(1,2):
        for j in range(10):
            k=start+2*j;result=cycle(k);last=result[-2];mapping=result[-1]
            need(a.induced(family(k+2),mapping)==last,'arbitrarily composed fixed graph identity')
            need(formula(k)[0]==210*(1<<start)*4**j+11,'explicit iteration count')
            iteration.append({'start':start,'step':j,'k':k,'order':len(family(k)),'state':formula(k)})
    report={'status':'PASS','scope':'Fixed pointed-template identity for every k; independent literal DP k1..5',
            'complete_marked_run_coefficients':runs,'literal_controls':literal,'iteration_controls':iteration}
    (HERE/'iterable_return.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS: literal reusable ready-return chart, full state, block map and ten-step compositions')


if __name__=='__main__':main()
