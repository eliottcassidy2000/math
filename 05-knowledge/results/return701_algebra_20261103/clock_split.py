"""Coprime 12/13 clocks on a formal joint module; no legal-word claim."""
import json
import sympy as s
from common import HERE,need,save
from radical_clocks import modular,power

def main():
    d=json.loads((HERE/'radical_clocks.json').read_text())
    B=s.Matrix([[0,1],[1,5]])
    U=modular(s.eye(27)+s.Matrix(d['nilpotent_matrix']),13)
    T=s.kronecker_product(B,U).applyfunc(lambda x:x%13);I=s.eye(54)
    need(power(B,12,13)==s.eye(2) and all(power(B,12//p,13)!=s.eye(2) for p in (2,3)),'primitive quadratic clock12')
    need(power(U,13,13)==s.eye(27) and U!=s.eye(27),'radical clock13')
    T12=power(T,12,13);T13=power(T,13,13)
    need(T12==s.kronecker_product(s.eye(2),power(U,12,13)),'twelfth sampling retains only nilpotent clock')
    need(T13==s.kronecker_product(B,s.eye(27)),'thirteenth sampling retains only quadratic clock')
    need((T13*power(T12,12,13)).applyfunc(lambda x:x%13)==T,'two clocks reconstruct full operation')
    need(power(T,156,13)==I and all(power(T,156//p,13)!=I for p in (2,3,13)),'exact joint order156')
    observation=s.Matrix([[2,5],[5,27]])
    need(int(observation.det())%13!=0,'primitive two-trace selector separates169 quadratic source states')
    save('clock_split.json',{'status':'EXACT_FORMAL_TENSOR_CLOCK','modulus':13,'quadratic_matrix':[[0,1],[1,5]],
        'quadratic_order':12,'unipotent_order':13,'tensor_dimension':54,'tensor_order':156,
        'sample12_order':13,'sample13_order':12,'primitive_trace_determinant':int(observation.det()),
        'scope':'formal joint linear module; not a native tournament attachment or a Collatz transition'})
    print('PASS formal tensor clock156: samples12/13 separate commuting nilpotent/quadratic factors',flush=True)
if __name__=='__main__':main()
