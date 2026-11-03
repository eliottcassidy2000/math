#!/usr/bin/env python3
"""Arbitrary-child controls, checking both literal words and cover coefficients."""
from pathlib import Path
import json
from scout import hp
from audit import need,perform,add,skeleton,ordered_covers_literal

HERE=Path(__file__).resolve().parent


def substitute(t,v,child):
    addresses=[(i,j)for i in range(len(t))for j in range(len(child)if i==v else 1)]
    rows=[sum(1<<k for k,(b,y)in enumerate(addresses)
              if (a==b==v and child[x]>>y&1)or(a!=b and t[a]>>b&1))for a,x in addresses]
    return rows,addresses


def expand_word(p,addresses,v,selected):
    return [addresses.index((i,j))for i in p for j in (selected if i==v else [0])]


def lift_birth(base,actual,mask,v,child):
    _,addresses=substitute(base,v,child)
    lifted=sum(1<<j for j,(i,k)in enumerate(addresses)if mask>>i&1)
    actual=add(actual,lifted);base=add(base,mask)
    need(actual==substitute(base,v,child)[0],'child birth naturality')
    return base,actual


def run(r,child,selected,repair=False):
    source=r['source'];target=r['target'];v=r['clone'];actual,addresses=substitute(source,v,child)
    if repair:
        base=source
        for mask in r['permanent_birth_masks']:base,actual=lift_birth(base,actual,mask,v,child)
        for word in r['exterior_flips']:
            old_count=len(actual);ext=base
            for mask in word['masks']:ext,actual=lift_birth(ext,actual,mask,v,child)
            _,addresses=substitute(ext,v,child)
            for key in ('first_path','second_path'):actual=perform(actual,expand_word(word[key],addresses,v,selected))
            actual=[row&((1<<old_count)-1)for row in actual[:old_count]]
            base=list(base);a,b=word['a'],word['b'];base[a]^=1<<b;base[b]^=1<<a
            need(actual==substitute(base,v,child)[0],'child exterior repair')
    else:
        for key in ('first_path','second_path'):actual=perform(actual,expand_word(r[key],addresses,v,selected))
        # The expanded two-complement endpoint has the old target restriction.
        base=[row&((1<<len(source))-1)for row in target[:len(source)]]
        need(actual==substitute(base,v,child)[0],'child common-frame endpoint')
        _,actual=lift_birth(base,actual,r['mask'],v,child)
    want=substitute(target,v,child)[0];need(actual==want,'full named child target')
    before=hp(tuple(substitute(source,v,child)[0]));after=hp(tuple(want))
    child_h=ordered_covers_literal(child,1)
    for core,h in ((source,before),(target,after)):
        coeff=skeleton(core,v)
        by_covers=sum(coeff[k]*ordered_covers_literal(child,k)for k in range(1,min(len(coeff),len(child)+1)))
        need(h==by_covers,'independent arbitrary-child cover evaluation')
    need(after==3*before+r['offset']*child_h,'exact child-scaled offset')
    return {'base_source_H':hp(tuple(source)),'offset':r['offset'],'child':child,'selected_path':selected,
            'child_H':child_h,'source_H':before,'target_H':after,'constant_offset_excess':r['offset']*(child_h-1)}


if __name__=='__main__':
    first=json.loads((HERE/'scout31.json').read_text())['programmes']
    short=[first[1],first[0],json.loads((HERE/'next_0.json').read_text())['programmes'][0]]
    repairs=json.loads((HERE/'returns.json').read_text())
    report=[]
    for child,selected in (([2,4,1],[0,1,2]),([2,4,9,3],[0,1,2,3])):
        report.extend(run(r,child,selected)for r in short)
        report.extend(run(r,child,selected,True)for r in repairs)
    (HERE/'child_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    for r in report:print(r['base_source_H'],r['child_H'],r['source_H'],r['target_H'],r['constant_offset_excess'])
