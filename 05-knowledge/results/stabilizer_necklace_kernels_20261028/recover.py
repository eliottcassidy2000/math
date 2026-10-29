"""Targeted niche recovery; reference counts are proxies, never priority claims."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
WAVES={
 'stabilizers_and_necklaces':r'Lyndon|Witt transform|necklace polynomial|distinguishing colo[u]?r|stabilizer polynomial',
 'automatic_and_hidden_coefficients':r'Mahler|Cartier|Moser|fibbinary|prefix.residue|weighted boundary',
 'exact_kernel_probe':r'\b4864\b|1\s*\+\s*[txz]\s*-\s*[txz](\^2|²)|two rotation.rigid|regular subset',
}
SOURCES={
 '01-canon/theorems/THM-2172-frobenius-collapse-of-mahler-and-twisted-cyclic-packets.md':'Dormant operation: use characteristic3 on a new ternary signed kernel; retain the integer lift.',
 '01-canon/theorems/THM-2005-support-dirichlet-automatic-tournament-atlas.md':'Under-titled Mahler block laws inside a reciprocal-profile atlas.',
 '01-canon/theorems/THM-2352-q-adic-prefix-residue-collision-spectrum.md':'Boundary: finite digit observations do not settle infinite tail predicates.',
 '05-knowledge/results/session68_tournament_composition_classification_20261026/session68-classification/CLASSIFICATION.md':'Under-titled stabilizer polynomial and two rigid palettes; recovered before any new-result claim.',
 '05-knowledge/results/signed_tournament_palette_rigidity_20261021/session23b-stabilizers/REPORT.md':'Classical regular-subset existence and setwise versus pointwise stabilizers.',
 '05-knowledge/results/context_leaves_prime_weights_and_union_20261011d.md':'Source-conditioned extension weights, as opposed to independent marginal counts.',
 '05-knowledge/results/collatz_circulant_20260930_circulants_lucas_cubic_monotile.md':'Adjacent lane: a fixed circulant spectrum does not control a changing-level cocycle.',
}
def rg(args):
    p=subprocess.run(['rg',*args],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    if p.returncode not in [0,1]: raise RuntimeError(p.stderr)
    return p.stdout.splitlines()

def main():
    waves={name:{'query':pat,'matched_files':rg(['-l','--glob','*.md',pat,'01-canon/theorems','05-knowledge/results'])}
           for name,pat in WAVES.items()}
    sources=[]
    for name,reason in SOURCES.items():
        path=ROOT/name
        needle='/'.join(path.parts[-3:]) if path.name in {'REPORT.md','CLASSIFICATION.md','README.md','RESULT.md','THEORY.md'} else path.name
        refs=rg(['-l','-F','--glob','*.md',needle,'00-navigation','01-canon','05-knowledge/results','05-knowledge/hypotheses','07-reflections'])
        # Exclude the new packet, whose own references would make the census drift.
        refs=[r for r in refs if 'stabilizer_necklace_kernels_20261028' not in r]
        text=path.read_text(encoding='utf-8')
        sources.append({'path':name,'selection':reason,'citation_needle':needle,'referring_document_count':len(refs),
                        'referring_documents':refs,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                        'opening_excerpt':'\n'.join(text.splitlines()[:12])})
    result={'status':'DISCOVERY/provenance, not mathematical evidence or a historical rediscovery census',
            'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'waves':waves,'selected_sources':sources}
    (HERE/'recovery.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'waves':{k:len(v['matched_files']) for k,v in waves.items()},
                      'sources':[{k:v for k,v in s.items() if k in ['path','selection','referring_document_count']} for s in sources]},indent=2))

if __name__=='__main__': main()
