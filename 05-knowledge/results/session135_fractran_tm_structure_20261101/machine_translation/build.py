"""Rebuild finite lists/tables from the frozen actual programs."""
from pathlib import Path
import json
import inherited_linear as linear
from machines import load_fractions,to_ir,to_tm,ir_to_fractran,save_json,save_tm
HERE=Path(__file__).resolve().parent

def main():
    records={}
    def put(name,rows,data,controls,start,initial,kind,extra=None):
        p=to_ir(rows,data,controls);p.update(name=name,initial_state='C'+str(start),initial_vector=initial)
        if extra:p.update(extra)
        save_json(HERE/(name+'.ir.json'),p)
        tm=to_tm(p,kind);tm.update(initial_state=p['initial_state'],initial_vector=initial)
        save_tm(HERE/(name+'.tm.json'),tm)
        back=ir_to_fractran(p)
        save_json(HERE/(name+'.roundtrip.json'),back)
        records[name]={'source_fractions':len(rows),'registers':len(data),'control_primes':len(controls),
                       'IR_instructions':len(p['instructions']),'TM_model':kind,'TM_states':len(tm['states']),
                       'TM_transitions':len(tm['rows']),'TM_tapes':tm['tapes'],'TM_alphabet_size':len(tm['alphabet']),
                       'roundtrip_fractions':len(back['fractions'])}
    put('pure_blue',load_fractions(HERE/'pure_blue.fractran'),[2,3,5,7,11,13,17,19],
        [23,29,31,37,41,43,47,53,59,61,67,71,73,79,83],23,[1,1,1,2,1,3,0,0],'packed_single')
    a=linear.compact_tail_matrix();p,_=linear.compile_linear(a);seed=[linear.ordinary_table(6)[6][h] for h in linear.COMPACT_H]
    rows=load_fractions(HERE/'compact_tail21.fractran')
    if rows!=p['fractions']:raise ValueError('S133 serialized input mismatch')
    put('compact_tail21',rows,p['x']+p['y']+[p['clock']],p['control']+[p['checkpoint']],p['checkpoint'],seed+[0]*7+[6],'unary_multi',{'matrix':a,'initial_order':6})
    full=[[64,44,46,18],[0,0,0,0],[0,0,0,0],[56,38,41,17]]
    for name,a,seed in [('migration_full',full,[0,0,0,1]),('migration_face',[[64,18],[56,17]],[0,1])]:
        p,_=linear.compile_linear(a)
        (HERE/(name+'.fractran')).write_text(''.join(str(n)+'/'+str(d)+'\n' for n,d in p['fractions']))
        put(name,p['fractions'],p['x']+p['y']+[p['clock']],p['control']+[p['checkpoint']],p['checkpoint'],seed+[0]*len(seed)+[0],'unary_multi',{'matrix':a})
    save_json(HERE/'BUILD.json',records)
    print(json.dumps(records,indent=2))
if __name__=='__main__':main()
