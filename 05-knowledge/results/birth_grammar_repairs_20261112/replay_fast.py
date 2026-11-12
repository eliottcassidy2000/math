"""Replay remaining producers under-O; long extension layers are separate."""
from pathlib import Path
import subprocess,sys,json,time
HERE=Path(__file__).resolve().parent
JOBS=[
('even_valuations.py',[]),('synchrony_probe.py',[]),('synchrony_kernels.py',[]),('synchrony_grammar.py',[]),
('repair_obstructions.py',[]),('lifting_roots.py',[]),('hamilton_clocks.py',['9']),
('direct_head_profiles.py',[]),('odd_axis_bank.py',[]),('new_axis_bank.py',['12']),('new_axis_bank.py',['13']),
('extension_heads.py',[]),('targeted_heads.py',[]),('rank_probe.py',[]),('rank_kernel.py',[]),
('free_word_growth.py',[]),('mixed_birth_sources.py',[]),('birth_extensions.py',['12','0','synchrony']),
('mixed_family_probe.py',[]),('mixed_grammar.py',[]),('repaired_rays.py',[]),('birth_inventory.py',[]),
('small_head_barriers.py',[]),('ozj_periodic.py',['9']),('ozj_mod2.py',[]),('audit_ozj.py',[]),('audit_math.py',[])]
def run():
    start=time.time();done=[]
    for script,args in JOBS:
        p=subprocess.run([sys.executable,'-O','-X','utf8',str(HERE/script),*args],capture_output=True,text=True)
        if p.returncode:
            print(p.stdout);print(p.stderr);raise RuntimeError(script)
        done.append([script,args]);print(json.dumps({'completed':script,'args':args,'jobs':len(done),'seconds':round(time.time()-start,1)}),flush=True)
    (HERE/'replay_jobs.json').write_text(json.dumps({'scope':'Exact listed-O producer runs completed successfully. Byte equality is checked separately.','jobs':done},indent=2)+'\n',encoding='utf8')
if __name__=='__main__':run()
