"""Companion finite-menu quotient and literal physical repair verification."""
import hashlib
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent/'dependencies'))
import core_verify as c

def main():
    root=Path(__file__).parent
    menus=json.loads((root/'MENUS.json').read_text())
    kept,_=c.normal_menu();ks=set(kept);counts={};dagger_checks=0
    for limit in (5,6,8):
        seen=set();expected=[]
        for w in kept:
            if len(w)>limit or w in seen:continue
            p=min(c.local_component(c.dagger(w)))
            c.need(p in ks and len(p)==len(w),'dagger preserves retained components')
            c.need(min(c.local_component(c.dagger(p)))==w,'dagger is an involution')
            seen|={w,p}
            expected.append({'word':w,'dagger_partner':p,'matrix':[list(row) for row in c.mm(c.Q,c.matrix(w))]})
            for u,v in ((2,3),(3,2),(3,5),(11,23),(55,3)):
                c.need(c.K(c.tau(u),c.state(w,c.tau(v)))==c.K(c.tau(v),c.state(p,c.tau(u))),
                       'actual index interchange')
                dagger_checks+=1
        c.need(expected==menus[str(limit)],'exact frozen menu')
        counts[limit]=len(expected)
    c.need(counts=={5:12,6:26,8:137},'finite quotient counts')
    raw=json.loads((root/'RAW-REPAIRS.json').read_text())
    repairs=json.loads((root/'REPAIRS.json').read_text());forms=c.form_data()
    c.need(raw['cap']==2134117 and raw['gaps5']==repairs['gaps5'],'exact finite prefix')
    c.need(len(raw['gaps5'])==65 and max(raw['gaps5'])==2134117,'length-five gap count')
    c.need(sorted(set(raw['gaps5'])-set(c.GAPS))==[r['target'] for r in repairs['repairs']],
           'every new finite gap has an actual repair')
    c.need(len(repairs['repairs'])==36,'repair count')
    for r,raw_record in zip(repairs['repairs'],raw['repairs']):
        n,i,u,v,k=raw_record;w=forms[i][0]
        c.need((n,u,v,w,k)==(r['target'],r['left_index'],r['right_index'],r['center'],len(w)),
               'C++ witness matches named certificate')
        word=c.graph_word(w,u,v);s=c.state(word[:-1])
        c.need(word==r['graph_word'] and list(s)==r['source_profile'] and c.F(s)==n,
               'full physical profile replay')
        rows=c.graph(word)
        c.need(len(rows)==r['vertices']<=n.bit_length()-1+6,'literal graph and integer size bound')
    old=json.loads((root/'dependencies/REPAIR_GRAPHS.json').read_text())
    c.need(len(old)==27 and all(r['order']<=r['H'].bit_length()-1+6 for r in old),
           'old literal repair order bound')
    # Direct counting is an extra control, not a replacement for the all-order compiler.
    literal=[]
    for r in repairs['repairs'][:6]:
        c.need(c.graph_counts(c.graph(r['graph_word']))==(r['target'],1),'literal C1 repair')
        literal.append(r['target'])
    return {'status':'PASS','menus':counts,'dagger_index_checks':dagger_checks,
            'length_five_finite_gaps':len(raw['gaps5']),'C1_repair_certificates':36,
            'inherited_repair_order_checks':27,'extra_literal_C1_checks':literal,
            'universal_occupancy':'OPEN'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
