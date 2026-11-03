"""Rebuild literal G25/G26 sources and expose the two-port birth polynomial."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import json
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'eleven_wedge_honeycomb_20261102'/'extension26.py'
spec=importlib.util.spec_from_file_location('native_birth',SOURCE)
native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)

def active(patterns,dimension):
    result={}
    for m,v in patterns:
        free=[e for e in range(dimension) if not m>>e&1]
        for bits in range(1<<len(free)):
            x=v|sum(((bits>>j)&1)<<e for j,e in enumerate(free))
            result[x]=result.get(x,0)+1
    return result

def main():
    E,P=native.bank(25);EE,PP=native.bank(26);new=[e for e in EE if e not in E]
    old=[native.encode(p,E) for p in P];fresh=[native.encode(p,E+new) for p in PP]
    A=active(old,len(E));B=active(fresh,len(EE));mask=(1<<len(E))-1
    hist=Counter((int((x&mask) in A),x>>len(E),h) for x,h in B.items())
    dead=Counter(x&mask for x in B if (x&mask) not in A)
    rows=[]
    for a in range(4):
        for was_active in (0,1):
            values={h:c for (flag,b,h),c in hist.items() if (flag,b)==(was_active,a)}
            rows.append({'old_active':bool(was_active),'new_bits':a,'positive_count_histogram':values,
                         'active_sources':sum(values.values()),'H_sum':sum(h*c for h,c in values.items()),
                         'H_square_sum':sum(h*h*c for h,c in values.items())})
    decoders=[]
    for attachment in (1,2):
        alive={x&mask:h for x,h in B.items() if x>>len(E)==attachment and (x&mask) not in A}
        fixed=value=0;trace=[]
        def energy(fixed,value):
            population=(1<<(32-fixed.bit_count()))-sum(x&fixed==value for x in A)
            if not population:return None,0
            correction=sum(h*h-2*h for x,h in alive.items() if x&fixed==value)
            return F(1)+F(correction,population),population
        initial,_=energy(fixed,value)
        if not initial<1:raise RuntimeError('missing guarded integer rounding threshold')
        for e in range(32):
            before,mass=energy(fixed,value);newfixed=fixed|1<<e
            options=[energy(newfixed,value|(bit<<e)) for bit in (0,1)]
            bit=min((b for b in (0,1) if options[b][1]),key=lambda b:options[b][0])
            if options[bit][0]>before:raise RuntimeError('conditional clamping increased energy')
            trace.append({'edge':E[e],'old_energy':str(before),'conditional_energies':[str(z[0]) for z in options],
                          'conditional_populations':[z[1] for z in options],'chosen':bit})
            fixed=newfixed;value|=bit<<e
        if value in A or alive.get(value)!=1:raise RuntimeError('rounding failed actual dead-source birth')
        fullsource=value|(attachment<<32)
        legalpaths=[p for p,(m,v) in zip(PP,fresh) if fullsource&m==v]
        if len(legalpaths)!=1:raise RuntimeError('unique actual target path missing')
        decoders.append({'attachment':attachment,'initial_energy':str(initial),'source25_mask':value,
                         'source26_mask':fullsource,'square_sum_path':legalpaths[0],'trace':trace})
    result={'edges25':len(E),'edges26':len(EE),'paths25':len(P),'paths26':len(PP),
            'new_edges':new,'old_active':len(A),'new_active':len(B),'old_dead_parents':len(dead),
            'old_dead_extension_multiplicities':dict(Counter(dead.values())),
            'rows':rows,'dead_source_decoders':decoders}
    (HERE/'birth_scout.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({**{k:v for k,v in result.items() if k!='dead_source_decoders'},
                      'dead_source_decoders':[{k:v for k,v in z.items() if k!='trace'} for z in decoders]},indent=2))

if __name__=='__main__':main()
