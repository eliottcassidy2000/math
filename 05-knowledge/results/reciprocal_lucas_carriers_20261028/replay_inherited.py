"""Replay, without rewriting, the preceding packet's actual transporter controls."""
from pathlib import Path
import sys,json,contextlib,io,types
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'moving_selector_atoms_20261028'))
import frames

def main(write=True):
    with contextlib.redirect_stdout(io.StringIO()):data=frames.controls()
    runtime=[];pending=list(sys.modules.values());seen=set()
    # Some inherited loaders do not register their modules in sys.modules.
    # Follow module-valued references too, rather than silently omitting them.
    while pending:
        module=pending.pop()
        if id(module) in seen:continue
        seen.add(id(module))
        source=getattr(module,'__file__',None)
        if not source:continue
        path=Path(source).resolve()
        if path.is_relative_to(ROOT) and not path.is_relative_to(HERE):
            runtime.append(path.relative_to(ROOT).as_posix())
            pending.extend(v for v in vars(module).values() if isinstance(v,types.ModuleType))
    out={'status':'FINITE-EXACT inherited implementation replay',
         'controls':data,'runtime_sources':sorted(set(runtime))}
    if write:(HERE/'inherited.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'groups':[r['group'] for r in data['atoms'][-2:]],
                      'hostile':data['hostile'],'F64':data['F64_shifts'],'F256':data['F256_shifts']}))
    return out
if __name__=='__main__':main()
