"""Normal/-O reproduction with exact-output comparison and source pins."""
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
JOBS=[('entry.py','entry.json'),('repair.py','repair.json'),('probe.py','probe.json'),
      ('next_returns.py','next_returns.json'),('response_recurrence.py','response_recurrence.json'),('star.py','star.json')]
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
    receipts=[]
    for script,result in JOBS:
        commands=[];normal=None
        for flags in ([],['-O']):
            cmd=[sys.executable,'-B',*flags,'-X','utf8',str(HERE/script)]
            output=subprocess.run(cmd,cwd=ROOT,text=True,encoding='utf-8',capture_output=True,check=True)
            content=(HERE/result).read_bytes()
            if normal is None:normal=content
            elif content!=normal:raise RuntimeError('Normal/-O receipt mismatch: '+result)
            commands.append(dict(optimized=bool(flags),stdout=output.stdout.strip()))
        receipts.append(dict(script=script,result=result,sha256=digest(HERE/result),commands=commands))
        print(script,'normal/-O exact match',flush=True)
    pins=list(HERE.glob('*.py'))+list(HERE.glob('*.md'))
    pins += [HERE.parent/p for p in [
        'strong_atom_continuation_vitali_20261026/strong-atlas.json',
        'polarized_path_calculus_20261103/FUNNELS.md',
        'zero_ready_component_returns_20261103/recovery/NOTE.md',
        'recursive_seed_switching_20261103/selector/bank.json',
        'recursive_seed_switching_20261103/selector/followup_named.json',
        'recursive_seed_switching_20261103/selector/inverse.py',
        'recursive_seed_switching_20261103/transport/followup.py',
        'recursive_seed_switching_20261103/transport/fast_gap.cpp',
        'zero_ready_component_returns_20261103/interface/gap.py',
        'zero_ready_component_returns_20261103/interface/check.py',
        'zero_ready_component_returns_20261103/interface/phase.py',
        'adaptive_metagraph_reformulation_20261103/metagraph_probe.py',
        'return56941_two_feedback_20261103/endpoint_checked.cpp',
        'session103_words_necklaces_20261029/bridge/BRIDGE_SINK_STAR_LEMMA.md',
        'session188_structural_bridges_20261104/skeleton/PROOF.md']]
    missing=[str(p) for p in pins if not p.exists()]
    if missing:raise RuntimeError('Missing inherited source pins: '+str(missing))
    report=dict(status='PASS_SAME_AGENT_WITH_INDEPENDENT_ALGORITHMS',normal_optimized_equal=True,jobs=receipts,
                hash_convention='SHA256 of UTF8 source bytes after CRLF-to-LF normalization',
                source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in sorted(pins)})
    (HERE/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
