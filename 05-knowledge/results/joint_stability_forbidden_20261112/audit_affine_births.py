"""Independent complete skeleton replay of every certified birth-clock slope."""
import hashlib,json
from pathlib import Path
from census_nine import read,save,need
from joint_nine import cones
HERE=Path(__file__).resolve().parent

def run():
    data=read('affine_births.json.gz');count=0
    for r in data['records']:
        O=[i for i in range(len(r['rows'])) if r['selected']>>i&1]
        D,U=cones.kernel(r['rows'],[r['port']],O)
        need(D=={(1,):1},'complete independent unit kernel')
        expected={(j,):c for j,c in enumerate(r['Stirling_coefficients'],1) if c}
        need(U==expected,'complete independent slope skeleton including absence of higher modes')
        count+=1
        if count%200==0:print('full independent birth-clock skeletons',count,flush=True)
    out={'scope':'Every one of1379 stable birth clocks through10, recounted by the independent full split-cover skeleton compiler. Entire defect and U tensors agree with one-sided-bound interpolation, not only sampled lengths.','full_kernels':count,'input_sha256':hashlib.sha256((HERE/'affine_births.json.gz').read_bytes()).hexdigest(),'all_equal':True}
    save('audit_affine_births.json.gz',out);print(json.dumps(out),flush=True)
if __name__=='__main__':run()
