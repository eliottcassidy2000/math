"""Targeted residual-gap transport followed by every old-apex recut."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'uniform_return_metagraph_20261111/reset'))
from classify_modules import H,strong,need,inherited as a,two_birth as b

BASE=(2046,2044,2024,2032,4068,4032,3968,3840,3584,1024,2048,527)

def flip(T,u,v):
    out=list(T);out[u]^=1<<v;out[v]^=1<<u;return tuple(out)

def main():
    source=BASE;standard=a.birth_beaten_only_by(a.birth_beaten_only_by(source,0),len(source))
    need(H(source)==1997 and H(standard)==5997,'actual defect3 positive control')
    recs=[]
    for direction,triangle in [('right',(2,3,4)),('left',(0,12,13))]:
      for mask in range(1,8):
        name=direction+str(mask)
        moved=standard
        for bit,v in enumerate(triangle):
            if mask>>bit&1:moved=flip(moved,1,v)
        ids=[v for v in range(len(moved)) if v!=11];R=a.induced(moved,ids)
        need(H(R)==(9 if mask==7 else 15),'exact partial/full gap transport residual H')
        vals=b.all_birth_counts(R)
        row=dict(direction=name,transport_arcs=[(1,v) for bit,v in enumerate(triangle) if mask>>bit&1],
                 residual=R,apex_recut_universe=len(vals),hits=[])
        for cut,h in enumerate(vals):
            if h not in (5993,5995):continue
            T=a.attach(R,cut);need(H(T)==h,'independent actual recut count')
            need(not strong(T),'all successful local partial-transport outputs are weak')
            runs=b.module_runs(T,ids.index(1))
            uniform=(runs[2]==3402 and not any(runs[3:]))
            hit=dict(cut=cut,H=h,rows=T,strong=False,module=ids.index(1),runs=runs,uniform=uniform)
            row['hits'].append(hit);print('HIT',name,cut,h,'uniform',uniform,'runs',runs,flush=True)
        print('DIRECTION',name,'hits',len(row['hits']),flush=True);recs.append(row)
    (HERE/'gap_return_probe.json').write_text(json.dumps(dict(status='FINITE-EXACT',scope='All14 nonempty partial/full transports across either neighboringtriangle and every apex recut; complete module coefficients for all scalar hits',records=recs),indent=2)+'\n')

if __name__=='__main__':main()
