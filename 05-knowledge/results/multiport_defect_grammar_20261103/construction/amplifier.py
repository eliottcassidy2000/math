"""Finite hostile probes for the all-order C2-image tripling theorem."""
from pathlib import Path
import json
from algebra import operator,response,need,strong,induced,birth,tables
HERE=Path(__file__).resolve().parent

def main():
    initial=operator((8,1,3,6),0,1)
    rows=tuple(initial['C2']);a=4;w=0;states=[]
    for level in range(4):
        h,f=response(rows,w);active=operator(rows,w,a)
        need(strong(rows),'strong current core')
        need(active['S_a_2']==h and active['L']+active['K']==f,'closed amplifier state')
        candidates=[];gaps=[]
        for aux in range(len(rows)):
            item=operator(rows,w,aux);gaps.append(item['S_a_2']-h)
            for recipe in ('predicted2','predicted4'):
                hh,ff=item[recipe]
                if ff==3*f and hh-3*h in (2,4):
                    candidates.append(dict(aux=aux,recipe=recipe,offset=hh-3*h))
        states.append(dict(level=level,rows=rows,w=w,aux=a,H=h,F=f,
                           C2_start_defects=gaps,two_recipe_digit_guards=candidates))
        if level==3:break
        next_rows=tuple(active['C2'])
        need(response(next_rows,w)==(3*h,3*f),'exact tripling')
        need(induced(next_rows,list(range(len(rows))))==rows,'all old arcs preserved')
        need(strong(active['X']) and strong(next_rows),'strong intermediate and target')
        a=len(rows);rows=next_rows
    clone_rows=tuple(initial['X']);clones=[]
    for m in range(1,7):
        cut=15|sum(1<<v for v in range(5,len(clone_rows)))
        clone_rows=birth(clone_rows,cut)
        hf=response(clone_rows,0)
        need(hf==(5+6*(1<<m),7+7*(1<<m)),'clone-block formula')
        need(strong(clone_rows),'strong clone block')
        clones.append(dict(m=m,H=hf[0],F=hf[1],rows=clone_rows))
    receipt=json.loads((HERE/'probe.json').read_text())
    clone_hits=receipt['two_birth_menus']['child2']['hits']
    need(len(clone_hits)==6,'six fixed-mark clone witnesses')
    for hit in clone_hits:
        t=hit['rows'];yy=(5,6,7)
        need(all(t[b]>>y&1 for b in range(4) for y in yy),'base dominates clone block')
        need(all(t[y]>>4&1 for y in yy),'clone block dominates x')
        need(tables(induced(t,yy))[0][-1]==1,'three-clone block transitive')
    result=dict(status='FINITE_EXACT_PROBES_OF_PROVED_GENERAL_OPERATOR',states=states,clone_blocks=clones,
        scope='The all-order tripling statement follows algebraically. The absence or presence of digit guards is asserted only for these listed states and two recipes.')
    (HERE/'amplifier.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{'level':v['level'],'H':v['H'],'F':v['F'],
        'C2_start_defects':v['C2_start_defects'],'digit_guards':v['two_recipe_digit_guards']}for v in states]))

if __name__=='__main__':main()
