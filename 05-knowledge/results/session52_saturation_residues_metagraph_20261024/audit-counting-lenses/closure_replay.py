#!/usr/bin/env python3
"""Copy only manifest-listed bytes, replay in isolation, and test warning preservation."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def require(condition,description):
    if not condition:
        raise RuntimeError(description)

def main():
    manifest=json.loads((HERE/'manifest.json').read_text())
    records=[]
    with tempfile.TemporaryDirectory(prefix='session52-counting-audit-') as temp:
        target=Path(temp)/'packet';target.mkdir()
        for name,expected in manifest.items():
            rel=Path(name)
            require(not rel.is_absolute() and '..' not in rel.parts,('unsafe manifest name',name))
            source=HERE/rel
            require(digest(source)==expected,('manifest source',name))
            destination=target/rel;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source,destination)
            require(digest(destination)==expected,('manifest copy',name))
        require(set(manifest)=={str(p.relative_to(target)) for p in target.rglob('*') if p.is_file()},'copy closure')
        # All runtime sources and inputs are inside the manifest-only temporary tree.
        for label,script,expected in (
            ('producer','inputs/producer/check.py','inputs/producer/receipt.json'),
            ('independent','independent_audit.py','audit-normal.json')):
            for optimized in (False,True):
                receipt=Path(temp)/f'{label}-{optimized}.json'
                command=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(target/script),'--receipt',str(receipt)]
                result=subprocess.run(command,cwd=temp,env={'PATH':os.environ['PATH'],'PYTHONHASHSEED':'12345'},capture_output=True,text=True,check=False)
                require(result.returncode==0,('relocated run failed',label,optimized,result.stdout,result.stderr))
                require(receipt.read_bytes()==(target/expected).read_bytes(),('relocated output changed',label,optimized))
                data=json.loads(receipt.read_text())
                records.append({'lane':label,'optimized':optimized,'exit_code':result.returncode,'checks':data['checks'],'receipt_sha256':digest(receipt),'isolated_python':True})
        source=target/'inputs/source-snapshots'
        provenance=json.loads((target/'inputs/producer/provenance.json').read_text())
        for item in provenance['sources']:
            require(digest(source/item['repository_path'])==item['sha256'],('source snapshot hash',item['repository_path']))
        original=(source/'05-knowledge/results/sc_partition_formula.md').read_bytes()
        patch=target/'inputs/producer/historical-warning.patch'
        result=subprocess.run(['patch','-p1','--batch','--forward','-i',str(patch)],cwd=source,capture_output=True,text=True,check=False)
        require(result.returncode==0,('historical patch',result.stdout,result.stderr))
        updated=(source/'05-knowledge/results/sc_partition_formula.md').read_bytes()
        require(updated.endswith(original),'historical body not preserved')
        added=updated[:-len(original)]
        require(added.startswith(b'> **HISTORICAL / SUPERSEDED') and added.endswith(b'\n\n'),'warning boundary')
        require(len(updated)==len(original)+len(added),'only prepend changed')
        require((source/'01-canon/theorems/THM-283-sc-burnside-formula.md').is_file(),'warning theorem link')
        require((source/'05-knowledge/results/converse_tournament_count_channels_20261021/README.md').is_file(),'warning lineage link')
        # A changed expected external baseline is detected in both Python modes.
        external=target/'inputs/external-oeis-prefixes.json'
        baseline=json.loads(external.read_text());baseline['A000568'][4]+=1
        external.write_text(json.dumps(baseline))
        for optimized in (False,True):
            command=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(target/'independent_audit.py'),'--receipt',str(Path(temp)/'must-not-exist.json')]
            result=subprocess.run(command,cwd=temp,capture_output=True,text=True,check=False)
            require(result.returncode!=0 and "external A" in result.stderr,('mutated baseline not rejected',optimized,result.stderr))
        result={'status':'PASS','copied_manifest_files':len(manifest),'manifest_sha256':digest(HERE/'manifest.json'),'runs':records,
                'all_provenance_source_snapshots_verified':len(provenance['sources']),
                'historical_warning':'Applies cleanly in temporary source tree, preserves entire original body byte-for-byte, and both redirect targets exist.',
                'hostile_baseline_mutation_rejected_normal_and_optimized':True,
                'closure':'Four executions use isolated Python with an unrelated working directory and only manifest-listed producer/auditor source, inputs, and frozen outputs.'}
        (HERE/'closure-replay.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,sort_keys=True))
if __name__=='__main__':
    main()
