"""Literal all-k one-deletion/three-birth branch and its finite core proof."""
from pathlib import Path
import json
from classify_modules import H,strong,need,inherited as a,two_birth as b,converse
HERE=Path(__file__).resolve().parent

def main():
    traces=[]
    canonical=(4,1,10,51,71,151,175,287,127)
    W=a.induced(canonical,range(7));ports=(5,6)
    hp,end=a.dp(W);_,start=a.dp(converse(W));full=(1<<7)-1
    E=sum(sum(end[mask][r] for r in ports)*hp[full^mask] for mask in range(1,full+1))
    C=sum(end[mask][i]*start[full^mask][j] for mask in range(1,full)
          for i in ports for j in ports if i!=j)
    need((H(W),E,C)==(41,102,24),'fixed seven-vertex response')
    residual=a.induced(canonical,range(5))
    coremap=a.iso(residual,a.CORES['S5a'])
    need(H(residual)==9 and coremap is not None,'five-core residual is S5a')
    for k in range(1,7):
        D=1<<k;source=a.ready_source((2,3,1),k)
        retained=[v for v in range(len(source)) if v!=k+1]
        base=a.induced(source,retained)
        need(base==a.ready_source((2,1),k),'middle deletion returns actual BA pentagon')
        sink1=a.attach(base,0);sink2=a.attach(sink1,0)
        final=a.attach(sink2,33*D)
        target=b.inflate(canonical,8,k)
        mapping=tuple(range(8,8+k))+(3,4,6,5,7,1,0,2)
        need(final==a.induced(target,mapping),'literal all-k target map')
        values=[H(x) for x in(source,base,sink1,sink2,final)]
        need(values==[34*D+5,14*D+3,14*D+3,14*D+3,102*D+17],'entire count trace')
        need(values[-1]==3*values[0]+2,'desired uniform ternary branch')
        need(strong(source) and strong(base) and strong(final),'strong endpoints and returned BA source')
        need(not strong(sink1) and not strong(sink2),'neutral padding retained honestly')
        need(a.induced(final,range(len(base)))==base,'every persistent arc retained')
        prof=a.root_profiles(final)[k+4]
        need((prof['s'],prof['q'],prof['gap'],prof['delta'])==(41,41,24,24),'neutral root response')
        traces.append({'k':k,'source_rows':source,'deleted_source_vertex':k+1,
                       'retained_source_vertices':retained,'returned_BA_rows':base,
                       'birth_outmasks':[0,0,33*D],'final_rows':final,
                       'target_map':mapping,'count_trace':values,'root_profile':prof})
    base=a.ready_source((2,1),1);padded=a.attach(a.attach(base,0),0)
    candidates=[]
    for out,h in enumerate(b.all_birth_counts(padded)):
        if h not in (221,223):continue
        rows=a.attach(padded,out);co=b.module_runs(rows,0)
        candidates.append({'outmask':out,'h':h,'run_coefficients':co,
                           'uniform':co[2]==102 and not any(co[3:])})
    need([c['outmask'] for c in candidates if c['h']==221 and c['uniform']]==[66],
         'unique uniform+2 cut in this exact tail schedule')
    need(not any(c['h']==223 and c['uniform'] for c in candidates),
         'no+4 repair by changing only final cut of this schedule')
    report={'status':'PASS','fixed_core_rows':W,'H_core':41,'E_ports_2':E,'C_ports_2':C,
            'residual_rows':residual,'residual_H':9,'residual_iso_to_S5a':coremap,
            'canonical_target_base_rows':canonical,'module':8,
            'run_coefficients':b.module_runs(canonical,8),'traces':traces,
            'final_cut_universe':1<<len(padded),'final_cut_candidates':candidates}
    (HERE/'uniform_return.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print('PASS: uniform141->425 analogue, fixed core(41,102,24), six literal scales')
if __name__=='__main__':main()
