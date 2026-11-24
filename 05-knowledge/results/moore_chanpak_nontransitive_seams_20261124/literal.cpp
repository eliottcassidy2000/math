// Compile with c++ -O3 -std=c++17 literal.cpp -o literal.
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;
using Rows=vector<unsigned>;
Rows sub(const Rows& r,const Rows& x,const Rows& y){
 vector<unsigned> tag,idx;vector<const Rows*> child;
 for(unsigned k=0;k<7;k++){
  const Rows* p=k==3?&x:k==4?&y:nullptr;
  for(unsigned i=0;i<(p?p->size():1);i++){tag.push_back(k);idx.push_back(i);child.push_back(p);}
 }
 int n=tag.size();Rows a(n);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){
  bool ok=tag[i]==tag[j]?((*child[i])[idx[i]]>>idx[j]&1):(r[tag[i]]>>tag[j]&1);
  if(ok)a[i]|=1u<<j;
 }
 return a;
}
U count(const Rows& r){
 unsigned n=r.size(),full=(1u<<n)-1;
 vector<U> dp((size_t)(full+1)*n);
 for(unsigned i=0;i<n;i++)dp[(size_t)(1u<<i)*n+i]=1;
 for(unsigned s=1;s<full;s++){
   U *F=&dp[(size_t)s*n];
   unsigned cur=s;
   while(cur){
     unsigned u=__builtin_ctz(cur);cur&=cur-1;
     if(!F[u])continue;
     unsigned q=r[u]&~s;
     while(q){
       unsigned bit=q&-q;q^=bit;
       dp[(size_t)(s|bit)*n+__builtin_ctz(bit)]+=F[u];
     }
   }
 }
 U ans=0;for(unsigned i=0;i<n;i++)ans+=dp[(size_t)full*n+i];return ans;
}
int main(int argc,char**){
 Rows r={122,100,25,82,34,76,20};
 Rows x={66,284,409,209,33,463,150,275,89};
 Rows y={62,124,56,112,32,64,21};
 if(argc>1){x={258,388,481,247,263,211,19,81,232};y={2,132,97,119,7,83,19,125};}
 U A=count(sub(r,x,y));r[5]^=1u<<6;r[6]^=1u<<5;U B=count(sub(r,x,y));
 if(argc>1){assert(A==85587188881ULL&&B==85587188879ULL);}
 else{assert(A==19116233339ULL&&B==19116233337ULL);}
 cout<<x.size()+y.size()+5<<"-vertex literal counts "<<A<<" -> "<<B<<endl;
}
