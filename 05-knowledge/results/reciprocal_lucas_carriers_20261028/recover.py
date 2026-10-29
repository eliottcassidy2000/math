"""Four reproducible, bounded waves of corpus discovery; hits are not truth."""
from pathlib import Path
import subprocess,json,random,base64
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
WAVES={
 'quadratic_modular':r'Ramanujan|Heegner|Weber|idoneal|singular moduli|class number',
 'lucas_transport':r'Lucas prime|rank of apparition|golden prime|primitive divisor|Fricke|Chebyshev|Pell',
 'quotient_geometry':r'equidecompos|equinumeros|Dehn invariant|scissors|Vitali|holonomy',
 'specific_coordinates':r'2508951|4372|4096/m|4096b|g_?\{?130|G_?\{?37|23, ?29, ?31|6/[bu]|276e|13, ?37',
}
def main():
    result={}
    for label,pattern in WAVES.items():
        args=['rg','--json','-i','-m','2','-g','*.md','-g','*.out',
              '-g','!**/sources/**','-g','!**/inputs/**','-g','!**/source/**',
              '-g','!**/audit*/**','-g','!**/replay/**','-g','!**/INDEX.md',
              '-g','!**/SESSION-LOG*','-g','!**/MISTAKES.md',
              '-g','!**/reciprocal_lucas_carriers_20261028/**',pattern,
              '01-canon/theorems','05-knowledge/results','07-reflections']
        run=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,check=False)
        if run.returncode not in (0,1):raise RuntimeError(run.returncode)
        hits={}
        for line in run.stdout.decode('utf-8').splitlines():
            v=json.loads(line)
            if v['type']!='match':continue
            d=v['data'];p=d['path']['text'].replace('\\','/')
            raw=d['lines'];excerpt=raw.get('text')
            if excerpt is None:excerpt=base64.b64decode(raw['bytes']).decode('utf-8',errors='replace')
            hits.setdefault(p,[]).append({'line':d['line_number'],'excerpt':excerpt.strip()[:320]})
        paths=sorted(hits);old=[p for p in paths if p.endswith('.out') or p.startswith('07-')]
        sample=random.Random('reciprocal-lucas-20261028-'+label).sample(old,min(10,len(old)))
        result[label]={'pattern':pattern,'matched_files':len(paths),'paths':paths,
                       'seeded_provenance_sample':{p:hits[p] for p in sample}}
        print(json.dumps({'wave':label,'matched_files':len(paths),'sample':sample}))
    (HERE/'recovery.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
if __name__=='__main__':main()
