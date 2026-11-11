"""Recover the growing-width cap as repeatedone-letter funnels; pay E(2)."""
from pathlib import Path
import json,sys
from funnel_six_state import advance,SINGLE,a
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'uniform_return_metagraph_20261111'/'reset'))
from classify_modules import H,strong,need,inherited,two_birth,converse

def cap(j):
    n=2*j+1;T=[sum(1<<v for v in range(u+1,n)) for u in range(n)]
    for i in range(1,j+1):
        parent=2*(i-1);child=2*i
        T[parent]^=1<<child;T[child]|=1<<parent
    return tuple(T)

def birth2(state):
    h,F,s=state;return(h+2*F,2*h+2*F-s,F)

def main():
    c=SINGLE;core=(0,);rows=[];e2=[]
    for j in range(13):
        if j:
            c=advance(c,(1,));G,r=a.funnel(core,{0},(1,));core=a.induced(G,[r]+list(range(r)))
        T=converse(cap(j));path=a.unique_start_path(T,0)
        need(path is not None and a.induced(T,path)==core,'literal all-depth cap/funnel graph identity controls')
        need(c[0]==c[2]==3**j,'exact endpoint-cap scalar coordinates')
        ready_h=c[0]+c[4]-1;e2.append(c[4])
        if j<=6:
            p,e,s=a.bin_values(core,0,3)
            need(c==(p[1],e[1],s[2],p[2],e[2],s[3]),'independent full six-state literal cap controls')
            a.tables.cache_clear()
            child=inherited.birth_beaten_only_by(core,0)
            profile=inherited.root_profiles(child)[-1]
            need(profile['h']==ready_h and profile['gap']==1,'literal cap ready newborn')
        if j>=3:need(e2[j]==6*e2[j-1]+3*e2[j-2]+12*e2[j-3],'exact three-coordinate recurrence')
        rows.append(dict(depth=j,state=c,ready_H=ready_h,feedback_parent=(j+1)//2))
    # The separate growing-defect birth orbit starting at the newQ family.
    growing=[]
    from target_ready_return import family
    for k in range(1,4):
        T=family(k,1);p=inherited.root_profiles(T)[-1];state=(p['h'],p['F'],p['s'])
        for step in range(3):
            expected=birth2(state)
            U=inherited.birth_beaten_only_by(T,len(T)-1)
            V=inherited.birth_beaten_only_by(U,len(U)-1)
            vp=inherited.root_profiles(V)[-1]
            need(expected==(vp['h'],vp['F'],vp['s']),'independent physical growing-defect transition')
            need(vp['h']-(3*state[0]+2)==2*(state[1]-state[0]-1),'exact compensation invoice')
            growing.append(dict(k=k,step=step,source_state=state,target_state=expected,
                                required_compensation=2*(state[1]-state[0]-1)))
            T,state=V,expected
    result=dict(status='PASS',scope='Fixed recursive cap graph identity, exact inherited compiler and independent bounded literal counts',
                cap_transfer=[[4,-2,6],[2,-1,0],[2,0,3]],
                cap_characteristic=[1,-6,-3,-12],cap_rows=rows,
                birth2_transfer=[[1,2,0],[2,2,-1],[0,1,0]],
                birth2_characteristic=[1,-3,-1,-1],growing_defect_controls=growing)
    (HERE/'cap_cover_state.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: endpoint cap needs growingcover state; exact growingdefect birth orbit')

if __name__=='__main__':main()
