"""Replay the proved-claim engines and independent audits in isolated processes."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time

HERE=Path(__file__).resolve().parent
JOBS=[
 ('arithmetic/check.py',[]),
 ('arithmetic/TRIANGLE_local.py',[]),
 ('arithmetic/TRIANGLE_check.py',[]),
 ('fibres/h9_fibres.py',[]),
 ('fibres/extension_content.py',[]),
 ('fibres/triangle_filtration.py',[]),
 ('fibres/gap_transport.py',[]),
 ('fibres/triangle_enumeration.py',[]),
 ('reset/funnel_six_state.py',['--order','15','--controls','9']),
 ('reset/ready_two_inverse.py',[]),
 ('reset/target_ready_return.py',[]),
 ('reset/finite_bridge.py',[]),
 ('reset/cap_cover_state.py',[]),
 ('reset/low_defect_escape.py',[]),
 ('reset/gap_return_probe.py',[]),
 ('reset/two_triangle_inverse.py',[]),
 ('reset/H15_inverse.py',[]),
 ('reset/H15_ready.py',[]),
 ('reset/audit_berggren.py',[]),
 ('berggren/check.py',[]),
 ('berggren/wide_return.py',[]),
 ('berggren/audit_reset.py',[]),
 ('berggren/audit_triangle.py',[]),
 ('berggren/audit_gap.py',[]),
 ('berggren/audit_all_k.py',[]),
 ('berggren/audit_defect_escape.py',[]),
 ('berggren/audit_H15.py',[]),
]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--optimized',action='store_true')
    args=parser.parse_args();records=[];failed=False
    for name,extra in JOBS:
        path=HERE/name;command=[sys.executable]+(['-O'] if args.optimized else [])+['-B',str(path)]+extra
        start=time.monotonic()
        try:
            process=subprocess.run(command,cwd=path.parent,capture_output=True,text=True,
                                   timeout=600,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
            ok=process.returncode==0
            record={'engine':name,'arguments':extra,'passed':ok,'returncode':process.returncode,
                    'stdout_tail':process.stdout[-2200:],'stderr_tail':process.stderr[-2200:]}
        except subprocess.TimeoutExpired:
            ok=False;record={'engine':name,'arguments':extra,'passed':False,'error':'600 second timeout'}
        records.append(record);failed|=not ok
        print(('PASS ' if ok else 'FAIL ')+name+' '+str(round(time.monotonic()-start,2))+'s',flush=True)
        if not ok:print(json.dumps(record),flush=True)
    report={'status':'PASS' if not failed else 'FAIL','optimized':args.optimized,
            'engine_count':len(records),'records':records,
            'scope':'New packet claims and inherited dependencies invoked by their declared engines; discovery-only probe_ready.py excluded.'}
    filename='REPLAY_OPTIMIZED.json' if args.optimized else 'REPLAY_NORMAL.json'
    (HERE/filename).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if failed:raise SystemExit(1)

if __name__=='__main__':main()
