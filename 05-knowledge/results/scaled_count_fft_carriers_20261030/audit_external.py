"""Record a pinned external Lean build and exact source bytes; do not run Lean."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess

HERE = Path(__file__).resolve().parent
PIN = 'c91cc3f4947caca8977e81a5d14244f73b5ae2b1'


def sha(data): return hashlib.sha256(data).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--build-log',type=Path,required=True)
    parser.add_argument('--pdf',type=Path,required=True)
    args=parser.parse_args()
    def git(*cmd): return subprocess.check_output(['git',*cmd],cwd=args.repo)
    if git('rev-parse','HEAD').decode().strip()!=PIN:
        raise RuntimeError('External checkout is not the pinned revision')
    if git('status','--porcelain'):
        raise RuntimeError('External checkout has unrecorded modifications')
    files=git('ls-tree','-r','--name-only','HEAD').decode().splitlines()
    source_files=[p for p in files if p.endswith('.lean') or p in
                  ['lean-toolchain','lakefile.toml','lake-manifest.json','comparator.json']]
    blobs={p:git('show',f'HEAD:{p}') for p in source_files}
    def code(path):
        # Only for a small syntactic receipt; the axiom audit is supplied by Lean.
        text=blobs[path].decode('utf-8')
        text=re.sub(r'/\*.*?\*/','',text,flags=re.S)
        # Lean nested block comments are removed by depth, not a flat regex.
        out=[]; i=0; depth=0
        while i<len(text):
            if text[i:i+2]=='/-': depth+=1; i+=2
            elif depth and text[i:i+2]=='-/': depth-=1; i+=2
            elif depth: i+=1
            elif text[i:i+2]=='--':
                j=text.find('\n',i); i=len(text) if j<0 else j
            else: out.append(text[i]); i+=1
        if depth: raise RuntimeError('Unclosed comment in syntactic audit')
        return ''.join(out)
    challenge=code('Challenge.lean'); definitions=code('FactoringAux/Definitions.lean')
    solution=code('Solution.lean')
    compact=lambda x: re.sub(r'\s+','',x)
    cdefs=challenge[challenge.index('def factor'):challenge.index('theorem subexponentialRate')]
    pdefs=definitions[definitions.index('def factor'):definitions.rindex('end Factoring')]
    if compact(cdefs)!=compact(pdefs): raise RuntimeError('Definition surface mismatch')
    sig=lambda s: compact(s[s.index('theorem subexponentialRate'):s.index(':=',s.index('theorem subexponentialRate'))])
    if sig(challenge)!=sig(solution): raise RuntimeError('Theorem signature mismatch')
    closure=set(); todo=['Solution.lean']
    while todo:
        path=todo.pop()
        if path in closure: continue
        closure.add(path)
        for module in re.findall(r'^\s*(?:public\s+)?import\s+([\w.]+)',code(path),re.M):
            dep=module.replace('.','/')+'.lean'
            if dep in blobs: todo.append(dep)
    if 'Challenge.lean' in closure: raise RuntimeError('Statement placeholder reached by proof')
    forbidden={p:re.findall(r'\b(?:sorry|admit|axiom|native_decide|unsafe|implemented_by)\b|\bpartial\s+def\b',code(p))
               for p in closure}
    forbidden={p:v for p,v in forbidden.items() if v}
    if forbidden: raise RuntimeError(forbidden)
    log=args.build_log.read_text(encoding='utf-8-sig')
    ax="'Factoring.subexponentialRate' depends on axioms: [propext, Classical.choice, Quot.sound]"
    if ax not in log or 'Build completed successfully (1949 jobs).' not in log:
        raise RuntimeError('Missing expected successful build and axiom output')
    pdf=args.pdf.read_bytes()
    if sha(pdf)!='e57f24481d2474a6f2c342d86312572f024710835ca583beb31e2bb1f9ceae40':
        raise RuntimeError('Wrong PDF version')
    sources=HERE/'sources'; sources.mkdir(exist_ok=True)
    (sources/'hexagon2610.00600v1.pdf').write_bytes(pdf)
    lines=['Pinned local build receipt; complete original log hash is in source_receipt.json.',
           'Commands: lake exe cache get; lake build',f'External commit: {PIN}',
           'Toolchain: leanprover/lean4:v4.35.0',
           'Challenge is the separate statement target and contains an intentional sorry.',
           'It is not imported by Solution; the proof axiom output below excludes sorryAx.']
    lines += [line for line in log.splitlines() if 'Built Solution' in line or ax in line or
              'Build completed successfully' in line or 'warning: Challenge.lean' in line]
    (HERE/'lean_build_receipt.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    manifest=json.loads(blobs['lake-manifest.json'])
    mathlib=next(p for p in manifest['packages'] if p['name']=='mathlib')
    receipt={'paper':{'title':'Deterministic Factoring in Time N^{o(1)}',
                     'authors':['Isaac Rudich','Louis-Martin Rousseau'],
                     'url':'https://hexagonmath.org/pdf/2610.00600v1',
                     'namespace':'Hexagon, not arXiv', 'license':'CC BY 4.0 (landing page)',
                     'bytes':len(pdf),'sha256':sha(pdf)},
             'repository':'https://github.com/IsaacRudich/SubexponentialDeterministicFactoring',
             'commit':PIN,'mathlib_commit':mathlib['rev'],
             'toolchain':blobs['lean-toolchain'].decode().strip(),
             'verification_date':'2026-10-30',
             'local_build':'SUCCESS (1949 jobs), exit 0',
             'axioms':['propext','Classical.choice','Quot.sound'],
             'full_build_log_sha256':sha(args.build_log.read_bytes()),
             'definitions_syntactically_equal_after_comments_and_whitespace':True,
             'main_theorem_signatures_equal':True,
             'local_solution_import_closure':sorted(closure),
             'challenge_in_solution_import_closure':False,
             'forbidden_token_scan_in_local_closure':forbidden,
             'source_blob_sha256':{p:sha(b) for p,b in sorted(blobs.items())},
             'scope':'Local Lean build, root reach, axiom output and stated definitions; not independent human review of every proof or a performance benchmark.'}
    (HERE/'source_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'pinned_source_files':len(blobs),'local_imports':len(closure),
                      'definition_match':True,'build':'SUCCESS','axioms':receipt['axioms']}))


if __name__=='__main__': main()
