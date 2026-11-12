"""Make the24 synchronized order12 births available to the extension grammar."""
from common import *
from synchrony_probe import expand

def run():
    out=[]
    for idx,r in enumerate(read('synchrony_kernels.json.gz')['records']):
        rr,O,labs=expand(r['rows'],r['selected'],{r['ports'][0]:2,r['ports'][1]:2})
        P,cert,_,_=certify(rr,[O]);P=P[O];h,u,w,_=tables(np.array(rr,dtype=np.int64))
        need(P and strong(rr) and defect(h,w,[v for v in range(len(rr)) if O>>v&1])==1,'physical synchronized source')
        _,_,C=quotient(rr,O);need(not any(len(c)>1 and set(c)&set(P) for c in C),'birth source')
        out.append({'birth_id':f'S12:{idx}','rows':rr,'selected':O,'ports':P,'stable_ports':P,'H':int(h),'U':sum(int(u[v]) for v in range(len(rr)) if O>>v&1),'source_family':idx,'labels':labs})
    save('birth_sources_synchrony_12.json.gz',{'scope':'All24 synchronized birth families at m=2 and free-port length1, ordinary marked classes.','records':out});print('sources',len(out),flush=True)
if __name__=='__main__':run()
