"""Verify sealed receipts and replay the small exhaustive mathematics."""
from hashlib import sha256
from pathlib import Path
import contextlib,io,json,sys
import reconstruct as r

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(path):
    data=path.read_bytes()
    if path.suffix!='.gz':data=data.replace(b'\r\n',b'\n')
    return sha256(data).hexdigest()

def seals():
    for name,base in [('sources.json',ROOT),('manifest.json',HERE)]:
        pins=json.loads((HERE/name).read_text())
        for path,wanted in pins.items():r.need(digest(base/path)==wanted,('changed pinned source or artifact',path))

def main():
    sys.set_int_max_str_digits(0);seals()
    import check_theory,boundary_profile,pair_square,cross_domain,incoming_check
    with contextlib.redirect_stdout(io.StringIO()):
        theory=check_theory.main();boundary=boundary_profile.main();pair_square.main();cross_domain.main()
        incoming=incoming_check.main()
    seals()
    a=json.loads((HERE/'counts-bucket-350.json').read_text())
    b=json.loads((HERE/'counts-reference-350.json').read_text())
    r.need(a['values']==b['values'] and len(a['values'])==351,'complete 350-prefix agreement')
    r.need(a['nonempty_strong']==b['nonempty_strong'],'strong prefix agreement')
    for packet in (a,b):
        h=sha256(json.dumps(packet['values'],separators=(',',':')).encode()).hexdigest()
        r.need(h==packet['values_sha256'],'exact integer prefix digest')
        A=packet['values'];S=packet['nonempty_strong']
        r.need(A[0]==1 and S[0]==0,'empty conventions')
        for n in range(1,351):
            r.need(A[n]==sum(S[k]*A[n-k] for k in range(1,n+1)),'ordered strong-component convolution')
    old=json.loads((HERE.parent/'strong_attachment_tail_compiler_20261028/ordinary-300.json').read_text())['values']
    r.need(a['values'][:301]==old and len(str(a['values'][-1]))==17646,'saved prefix and endpoint digits')
    root=json.loads((HERE/'lazy-root-receipt.json').read_text())
    for row in root['results']:
        N=row['order'];r.need(row['classes']==a['values'][N],'constructive/scalar comparison')
        r.need(row['classification_requests']==N*row['classes']==row['audit']['trace'],'root trace')
    r.need(root['independent_bank']['classes']==191536,'complete bank receipt')
    saved_incoming=json.loads((HERE/'incoming-receipt.json').read_text())
    r.need(all(saved_incoming[k]==v for k,v in incoming.items()),'incoming cheap rechecks')
    r.need(saved_incoming['full_bank_audit']['classes']==9733056 and
           saved_incoming['full_bank_audit']['full_class_uniqueness'] and
           saved_incoming['full_bank_audit']['labelled_mass']==1<<45,'incoming complete ten-bank receipt')
    result=dict(status='PASS',python_optimized=not __debug__,
                replayed_contexts=sum(theory['witness']['context_levels']),
                replayed_profile_updates=sum(x['updates'] for x in boundary['profile_updates']),
                closed_cells_at_seven=boundary['seven_closed_carrier_cells'],
                scalar_prefix=350,scalar_digits=17646,lazy_root_complete_order=9,incoming_full_bank_order=10,
                all_size_scope='Sufficient K; reconstruction refuted first at seven; minimality open')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
