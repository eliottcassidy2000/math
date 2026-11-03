"""Cross-check the independently published seam bank in the physical basis."""
from pathlib import Path
import argparse,hashlib,json
import sympy as sp
from verify import stats,join,reverse_labels,actual_complement,converse,literal_target,h
from transfer import SMALL,apply3,vector3
from extension import strong
HERE=Path(__file__).resolve().parent
OTHER=HERE.parent/'golden_seam_recurrence_20261103'

def need(ok,why):
    if not ok:raise RuntimeError(why)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='integration.json');args=parser.parse_args()
    atlas=json.loads((OTHER/'atom_atlas.json').read_text())
    L=sp.Matrix([[1,0,0],[1,-sp.Rational(1,2),sp.Rational(1,2)],[0,sp.Rational(1,2),sp.Rational(1,2)]])
    counts={};derived={}
    for label,entry in atlas.items():
        rows=tuple(entry['source_rows']);n=len(rows)
        need(strong(rows),('source strong',label))
        need(all(rows[i]>>(i+1)&1 for i in range(n-1)),('source spine',label))
        counts[n]=counts.get(n,0)+1
        block=reverse_labels(actual_complement(rows));gs=[]
        for j in range(3):
            inverse=stats(reverse_labels(converse(block)))
            hh=sum(inverse[:2]);ss=inverse[1];dd=inverse[2]
            gs.append(sp.Matrix([[hh,hh-ss,dd]]))
            block=join(block,(0,))
        cc=(gs[2]-gs[0])/2;aa=gs[1]-cc;bb=gs[0]-aa
        mat=sp.Matrix.vstack(aa,bb,cc).T
        need(mat*L==L*sp.Matrix(entry['matrix']),('independent coordinate transport',label))
        need(all(x.is_Integer and x>=0 for x in mat),('physical integer positive matrix',label))
        derived[label]=tuple(tuple(int(mat[i,j]) for j in range(3)) for i in range(3))
    need(counts=={1:1,3:1,4:5,5:50,6:903},('all supplied charts',counts))
    for label,mat in SMALL.items():
        other={'s':'a','c':'b'}.get(label,label)
        need(mat==derived[other],('golden matrices agree',label))
    patches={47:('a','s14'),83:('a','s268'),103:('a','s330'),491:('a','a','s248','a'),
             709:('a','s261','a','a','a'),833:('a','s85','b'),3227:('a','r35','s262'),4619:('a','s453','r29')}
    extra={k:tuple(v['source_rows']) for k,v in atlas.items()}
    witnesses={}
    for target,source_word in patches.items():
        processing=tuple(reversed(source_word));rows=literal_target(processing,extra)
        v=(1,0,0)
        for letter in processing:v=apply3(v,derived[letter])
        need(sum(v[:2])==target and h(rows)==target and strong(rows),('literal patch',target))
        witnesses[str(target)]={'source_order_word':source_word,'processing_word':processing,'order':len(rows),'rows':rows,'state':v}
    need(vector3(('s','c'))==vector3(('q5',))==(4,1,9),'one-end identical states')
    need(h(literal_target(('s','s','c')))==9 and h(literal_target(('s','q5')))==13,'other-end hostile')
    coverage=json.loads((OTHER/'coverage.json').read_text())
    pins={name:hashlib.sha256((OTHER/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
          for name in ('README.md','THEORY.md','COVERAGE.md','atom_atlas.json','coverage.json')}
    receipt={'status':'PASS','independently_transported_matrices':len(atlas),'strong_chart_counts':counts,
             'basis_L':[[str(x) for x in row] for row in L.tolist()],
             'literal_six_vertex_patches':witnesses,
             'inherited_holdout_coverage':{k:v for k,v in coverage['selected_six'].items() if k!='values'},
             'inherited_coverage_scope':'30000 holdout census cited from concurrent packet, not independently rerun here',
             'source_hash_basis':'LF normalized SHA256','source_pins':pins}
    (HERE/args.output).write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','matrices':len(atlas),'chart_counts':counts,'literal_patch_values':list(patches)}))
if __name__=='__main__':main()
