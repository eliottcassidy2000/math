"""Same-language, same-target isolation of projective gauge and {1,3} contraction."""
import json,time,statistics,hashlib
from pathlib import Path
from counter import count
rows=[]
for n in (120,200):
    modes={'raw':dict(gauge=False),'projective_direct':dict(terminal3=False),'projective_terminal':dict(terminal3=True)}
    times={k:[] for k in modes};ref=None
    for rep in range(3):
        names=list(modes);names=names[rep:]+names[:rep]
        for name in names:
            t=time.perf_counter();a,stats=count(n,**modes[name]);times[name].append(time.perf_counter()-t)
            if ref is None:ref=a
            if ref!=a:raise AssertionError((n,name))
    rows.append({'n':n,'seconds':times,'median_seconds':{k:statistics.median(v) for k,v in times.items()}})
report={'counter_sha256':hashlib.sha256(Path(__file__).with_name('counter.py').read_bytes()).hexdigest(),'rows':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
