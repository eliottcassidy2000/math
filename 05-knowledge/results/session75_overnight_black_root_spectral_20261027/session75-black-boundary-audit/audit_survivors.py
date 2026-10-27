#!/usr/bin/env python3
"""Separate survivor audit; not a modification of initial obstruction packet."""
import json,pathlib
from audit_checks import *
from audit_families import build
H=pathlib.Path(__file__).resolve().parent

def main():
    report=[]
    for a in list(range(1,26))+[100]:
        for odd in (False,True):
            t,u=build(a,odd);m=a+4;r=m-1;w=2*m;v=a+2;n=len(t)
            p=tuple(range(a+1,-1,-1))+tuple(range(m,2*m))+((w,) if odd else ())+(r,v)
            q=(r,a+1)+((w,) if odd else ())+tuple(range(2*m-1,m-1,-1))+tuple(range(a))+(a,v)
            terminal,moves=replay(t,[p,q],require_core=True,detail=True)
            need(core(terminal)==[v] and degrees(terminal)[v]==n-1,'wrong named source/core')
            TU=tuple(tuple(terminal[i][j] for j in range(m,2*m)) for i in range(m,2*m));need(TU==u,'U not restored')
            Vcore=[i for i in range(m) if i!=v];need(strong(terminal,Vcore),'V-minus-source not strong')
            need(all(terminal[i][r] for i in range(a)) and terminal[r][a] and terminal[a][a+1] and terminal[a+1][r],'V-core directed triangle links')
            need(all(terminal[b][i] for b in (a,a+1) for i in range(a)),'B to A core links')
            U=range(m,2*m)
            if odd:
                need(all(terminal[j][i] for i in U for j in Vcore),'odd direct blocks not uniform')
                need(terminal[m+r][w] and terminal[w][a+1],'odd exceptional core bridge')
                need(all(terminal[w][i] for i in U if i!=m+r) and all(terminal[j][w] for j in Vcore if j!=a+1),'odd other middle arcs')
            else:
                need([(i,j) for i in U for j in Vcore if terminal[i][j]]==[(m+r,a+1)],'even unique core bridge')
            mid=reverse(t,p);need(degrees(mid)[v]==n-2 and degrees(mid)[r]==a+6+int(odd),'second core step endpoint scores')
            report.append({'a':a,'n':n,'odd':odd,'named_source':v,'terminal_rows':rows(terminal),'moves':moves,'core_distance':2})
    # Bind all forty producer survivor cases to the independently built named edges.
    src=json.loads((H/'producer/survivors/normal.json').read_text())
    for theirs in src['cases']:
        mine=next(z for z in report if (z['a'],z['odd'])==(theirs['a'],theirs['odd']))
        for x,y in zip(mine['moves'],theirs['route']):
            need(list(x['R_path'])==y['R_word'] and list(x['K_path'])==y['K_word'],'survivor named words differ')
            need(x['K_source_rows']==y['K_source'] and x['K_target_rows']==y['K_target'],'survivor native endpoints differ')
    parent10=json.loads((H/'inputs/parent-counterexample10.json').read_text());a10=from_rows(parent10['rows'])
    end10,moves10=replay(a10,[z['R_path'] for z in parent10['moves']],require_core=True,detail=True)
    for x,y in zip(moves10,parent10['moves']):
        need(x['K_source_rows']==y['K_source_rows'] and x['K_target_rows']==y['K_target_rows'],'parent n10 native endpoints differ')
    # Producer parent's n11 certificate is a relabeling, not same raw source.
    c=json.loads((H/'inputs/candidate11.json').read_text());a=from_rows(c['rows']);standard,_=build(1,True);phi=list(range(5))+list(range(6,11))+[5]
    need(all(standard[i][j]==a[phi[i]][phi[j]] for i in range(11) for j in range(11)),'parent n11 isomorphism mismatch')
    cert=next(w for w in c['result']['B'] if w['verdict']=='WITNESS');words=[p[::-1] if cert['dual'] else p for p in cert['route']];end,moves=replay(a,words,require_core=True,detail=True)
    print(json.dumps({'status':'PASS','scope':'separate all-family core survivor; independent exact controls a1..25 and100; structural proof assessed separately','cases':report,'parent_n10_distinct_core_route':{'moves':moves10,'terminal_core_ports':core(end10)},'parent_n11_named_label_binding':{'family_to_parent_permutation':phi,'parent_rows':rows(a),'moves':moves,'terminal_core_ports':core(end)},'unique_unbalanced_minimum_shore_repair':'structural proof accepted; no finite-census inference'},sort_keys=True,indent=2))
if __name__=='__main__':main()
