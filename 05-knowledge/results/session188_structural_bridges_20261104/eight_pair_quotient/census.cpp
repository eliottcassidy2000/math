#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
vector<unsigned> graph(unsigned n,unsigned code){vector<unsigned>a(n);unsigned k=0;for(unsigned i=0;i<n;i++)for(unsigned j=i+1;j<n;j++,k++){if(code>>k&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
unsigned seven_endpoint_bit(const vector<unsigned>&a){
 unsigned d=0,out=0;for(unsigned i=0;i<7;i++)d|=(__builtin_popcount(a[i])&1)<<i;
 for(unsigned i=0;i<7;i++){unsigned square=0;for(unsigned j=0;j<7;j++)if(a[i]>>j&1)square^=a[j];unsigned incoming=127^a[i]^(1u<<i);unsigned bit=((d>>i)&1)^(__builtin_popcount(square&(d^incoming))&1);out|=bit<<i;}return out;
}
unsigned delete_code(const vector<unsigned>&a,unsigned v){unsigned code=0,k=0;for(unsigned i=0;i<8;i++)if(i!=v)for(unsigned j=i+1;j<8;j++)if(j!=v){code|=((a[i]>>j)&1)<<k;k++;}return code;}
bool prime(const vector<unsigned>&a){for(unsigned m=1;m<255;m++){if(__builtin_popcount(m)<2)continue;bool module=true;for(unsigned i=0;i<8;i++)if(!(m>>i&1)&&((a[i]&m)!=0)&&((a[i]&m)!=m)){module=false;break;}if(module)return false;}return true;}
int main(){
 vector<unsigned char>table(1u<<21);vector<unsigned>zero;
 for(unsigned c=0;c<(1u<<21);c++){table[c]=seven_endpoint_bit(graph(7,c));if(!table[c])zero.push_back(c);}
 ofstream tablefile("t7.bin",ios::binary);tablefile.write(reinterpret_cast<const char*>(table.data()),table.size());tablefile.close();if(!tablefile)throw runtime_error("seven-vector write failure");
 uint64_t tried=0,blind=0,primes=0;vector<unsigned>first,firstprime;ofstream bank("blind8.bin",ios::binary);
 for(unsigned c:zero)for(unsigned mask=0;mask<128;mask++){
  unsigned code=(c<<7)|mask;auto a=graph(8,code);tried++;bool pass=true;
  for(unsigned j=1;j<8;j++)if(table[delete_code(a,j)]){pass=false;break;}
  if(pass){blind++;if(first.size()<16)first.push_back(code);if(prime(a)){primes++;if(firstprime.size()<16)firstprime.push_back(code);}for(unsigned k=0;k<4;k++){char b=char((code>>(8*k))&255);bank.write(&b,1);}}
 }
 bank.close();if(!bank)throw runtime_error("residual write failure");
 cout<<"{\"seven_zero_vectors\":"<<zero.size()<<",\"extension_candidates\":"<<tried<<",\"eight_shared_cubic_blind\":"<<blind<<",\"prime_blind\":"<<primes<<",\"first_codes\":[";for(unsigned i=0;i<first.size();i++){if(i)cout<<",";cout<<first[i];}cout<<"],\"first_prime_codes\":[";for(unsigned i=0;i<firstprime.size();i++){if(i)cout<<",";cout<<firstprime[i];}cout<<"]}\n";
}
