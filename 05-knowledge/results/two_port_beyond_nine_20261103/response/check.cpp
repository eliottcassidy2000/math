// Literal permutation/cut audit of full ordered endpoint tensors, independently
// compared with squarefree subset multiplication and the census implementation.
#define main census_entry_unused
#include "../census/census.cpp"
#undef main
#include <numeric>

struct Tensor {
 int n;std::vector<U> s,e,m,c,u,v,d;
 explicit Tensor(int size):n(size),s(n),e(n),m(n*n),c(n*n),u(n*n*n),v(n*n*n),d(n*n*n*n){}
 int i2(int a,int b)const{return a*n+b;}
 int i3(int a,int b,int c)const{return (a*n+b)*n+c;}
 int i4(int a,int b,int c,int d)const{return ((a*n+b)*n+c)*n+d;}
 bool operator==(const Tensor&t)const{return s==t.s&&e==t.e&&m==t.m&&c==t.c&&u==t.u&&v==t.v&&d==t.d;}
};
Rows from_code(int n,U code){Rows a{};int b=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++b)
 if((code>>b)&1)a[i]|=1u<<j;else a[j]|=1u<<i;return a;}
U encode_rows(int n,const Rows&a){U c=0;int b=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++b)if((a[i]>>j)&1)c|=U(1)<<b;return c;}
Tensor literal(int n,const Rows&a){
 Tensor t(n);std::vector<int> p(n);std::iota(p.begin(),p.end(),0);
 do{
  unsigned bad=0;for(int k=0;k<n-1;++k)if(!((a[p[k]]>>p[k+1])&1))bad|=1u<<k;
  if(!bad){++t.s[p[0]];++t.e[p[n-1]];++t.m[t.i2(p[0],p[n-1])];}
  for(int k=0;k<n-1;++k)if(!(bad&~(1u<<k))){
   ++t.c[t.i2(p[k],p[k+1])];++t.u[t.i3(p[k],p[k+1],p[n-1])];++t.v[t.i3(p[0],p[k],p[k+1])];
  }
  for(int k=0;k<n-1;++k)for(int l=k+1;l<n-1;++l)if(!(bad&~((1u<<k)|(1u<<l))))
   ++t.d[t.i4(p[k],p[k+1],p[l],p[l+1])];
 }while(std::next_permutation(p.begin(),p.end()));return t;
}
Tensor squarefree(int n,const Rows&a){
 Tensor t(n);unsigned full=(1u<<n)-1;
 std::vector<U> g((full+1)*n*n),e((full+1)*n),s((full+1)*n);
 auto at=[&](unsigned mask,int start,int end)->U&{return g[(mask*n+start)*n+end];};
 for(int r=0;r<n;++r){at(1u<<r,r,r)=1;
  for(unsigned mask=1;mask<=full;++mask)if(mask&(1u<<r))for(int last=0;last<n;++last){
   U value=at(mask,r,last);if(!value)continue;unsigned next=a[last]&~mask&full;
   while(next){int v=__builtin_ctz(next);next&=next-1;at(mask|(1u<<v),r,v)+=value;}
  }
 }
 for(unsigned mask=1;mask<=full;++mask)for(int i=0;i<n;++i)for(int j=0;j<n;++j){U value=at(mask,i,j);s[mask*n+i]+=value;e[mask*n+j]+=value;}
 for(int i=0;i<n;++i){t.s[i]=s[full*n+i];t.e[i]=e[full*n+i];for(int j=0;j<n;++j)t.m[t.i2(i,j)]=at(full,i,j);}
 for(unsigned A=1;A<full;++A){unsigned B=full^A;
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){
   t.c[t.i2(i,j)]+=e[A*n+i]*s[B*n+j];
   for(int k=0;k<n;++k){t.u[t.i3(i,j,k)]+=e[A*n+i]*at(B,j,k);t.v[t.i3(i,j,k)]+=at(A,i,j)*s[B*n+k];}
  }
  unsigned rest=full^A;
  for(unsigned B=(rest-1)&rest;B;B=(B-1)&rest){unsigned C=rest^B;if(!C)continue;
   for(int i=0;i<n;++i)if(e[A*n+i])for(int j=0;j<n;++j)for(int k=0;k<n;++k)if(at(B,j,k))
    for(int l=0;l<n;++l)if(s[C*n+l])t.d[t.i4(i,j,k,l)]+=e[A*n+i]*at(B,j,k)*s[C*n+l];
  }
 }
 return t;
}
U response(const Tensor&t,unsigned a,unsigned c){
 // a_i=old i->x; c_i=old i->y; x->y. All four tensors retain named indices.
 int n=t.n;auto bit=[](unsigned m,int i){return (m>>i)&1;};U z=0;
 for(int i=0;i<n;++i)z+=t.e[i]*bit(a,i)+t.s[i]*(1-bit(c,i));
 for(int i=0;i<n;++i)for(int j=0;j<n;++j){
  z+=t.c[t.i2(i,j)]*bit(a,i)*(1-bit(c,j));
  z+=t.m[t.i2(i,j)]*((1-bit(a,i))*bit(c,j)+(1-bit(c,i))*bit(a,j));
  for(int k=0;k<n;++k){
   z+=t.u[t.i3(i,j,k)]*(bit(a,i)*(1-bit(a,j))*bit(c,k)+bit(c,i)*(1-bit(c,j))*bit(a,k));
   z+=t.v[t.i3(i,j,k)]*((1-bit(a,i))*bit(c,j)*(1-bit(c,k))+(1-bit(c,i))*bit(a,j)*(1-bit(a,k)));
   for(int l=0;l<n;++l)z+=t.d[t.i4(i,j,k,l)]*(bit(a,i)*(1-bit(a,j))*bit(c,k)*(1-bit(c,l))
                         +bit(c,i)*(1-bit(c,j))*bit(a,k)*(1-bit(a,l)));
  }
 }
 return z;
}
U hamilton(int n,const std::vector<unsigned>&a){
 unsigned full=(1u<<n)-1;std::vector<U> dp((full+1)*n);std::vector<unsigned> pred(n);
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)if((a[i]>>j)&1)pred[j]|=1u<<i;
 for(unsigned mask=1;mask<=full;++mask){unsigned vs=mask;
  while(vs){int v=__builtin_ctz(vs);vs&=vs-1;unsigned rest=mask^(1u<<v),bs=rest&pred[v];U val=!rest;
   while(bs){int u=__builtin_ctz(bs);bs&=bs-1;val+=dp[rest*n+u];}dp[mask*n+v]=val;
  }
 }
 U z=0;for(int i=0;i<n;++i)z+=dp[full*n+i];return z;
}
U direct(int n,const Rows&old,unsigned a,unsigned c){
 std::vector<unsigned> rows(n+2);for(int i=0;i<n;++i){rows[i]=old[i];
  if((a>>i)&1)rows[i]|=1u<<n;else rows[n]|=1u<<i;
  if((c>>i)&1)rows[i]|=1u<<(n+1);else rows[n+1]|=1u<<i;
 }rows[n]|=1u<<(n+1);return hamilton(n+2,rows);
}
Canon brute_canon(const Data&d){
 std::vector<int> p(N);std::iota(p.begin(),p.end(),0);Canon out;
 do{Key key;for(int u:p){key.push_back(d.s[u]);key.push_back(d.e[u]);}
  // Reject by the full vertex prefix before constructing the edge suffix.
  if(out.aut&&std::lexicographical_compare(out.key.begin(),out.key.begin()+2*N,key.begin(),key.end()))continue;
  for(int u:p)for(int v:p)key.push_back(d.c[u][v]);++out.leaves;
  if(out.aut&&key>out.key)continue;U code=code_at(d.a,p);
  if(!out.aut||key<out.key){out.key=key;out.aut=1;out.mincode=out.maxcode=code;}
  else{++out.aut;out.mincode=std::min(out.mincode,code);out.maxcode=std::max(out.maxcode,code);}
 }while(std::next_permutation(p.begin(),p.end()));return out;
}
void array_json(const std::vector<U>&v){std::cout<<"[";for(std::size_t i=0;i<v.size();++i){if(i)std::cout<<",";std::cout<<v[i];}std::cout<<"]";}
int main(int argc,char**argv){try{
 need(argc==2,"usage: response-check controls.txt");std::ifstream input(argv[1]);need(bool(input),"controls open");
 int count;input>>count;std::cout<<"{\"controls\":[";
 for(int index=0;index<count;++index){int n;U code;input>>n>>code;need(bool(input)&&n>=3&&n<=10,"valid control");auto rows=from_code(n,code);
  auto a=literal(n,rows),b=squarefree(n,rows);need(a==b,"literal/squarefree full tensor disagreement");
  U aut=0;if(n==10){auto d=make(rows,true);for(int i=0;i<n;++i){need(d.s[i]==a.s[i]&&d.e[i]==a.e[i],"census endpoints");for(int j=0;j<n;++j)need(d.c[i][j]==a.c[a.i2(i,j)],"census covers");}
   auto x=canon(d),y=brute_canon(d);need(x.aut==y.aut,"refined/brute symmetry order");need((x.mincode==x.maxcode)==(y.mincode==y.maxcode),"refined/brute extra symmetries");aut=x.aut;
  }
  std::vector<unsigned> masks;for(unsigned m=0;m<(1u<<n);++m)if(__builtin_popcount(m)<=2)masks.push_back(m);
  std::vector<U> table;for(unsigned x:masks)for(unsigned y:masks){U actual=direct(n,rows,x,y);need(response(a,x,y)==actual,"full response decomposition");table.push_back(actual);}
  // Explicit high-weight queries exercise interpolation outside the battery.
  std::vector<std::int64_t> coefficients(table.size());
  for(std::size_t i=0;i<masks.size();++i)for(std::size_t j=0;j<masks.size();++j){std::int64_t val=0;
   for(std::size_t k=0;k<masks.size();++k)if(!(masks[k]&~masks[i]))for(std::size_t l=0;l<masks.size();++l)if(!(masks[l]&~masks[j])){
    int sign=(__builtin_popcount(masks[i]^masks[k])+__builtin_popcount(masks[j]^masks[l]))%2?-1:1;
    val+=sign*std::int64_t(table[k*masks.size()+l]);
   }coefficients[i*masks.size()+j]=val;
  }
  for(unsigned x:std::vector<unsigned>{(1u<<n)-1,((1u<<n)-1)^1,((1u<<n)-1)^3})for(unsigned y:std::vector<unsigned>{(1u<<n)-1,((1u<<n)-1)^2,((1u<<n)-1)^5}){
   std::int64_t val=0;for(std::size_t i=0;i<masks.size();++i)if(!(masks[i]&~x))for(std::size_t j=0;j<masks.size();++j)if(!(masks[j]&~y))val+=coefficients[i*masks.size()+j];
   need(val>=0&&U(val)==direct(n,rows,x,y)&&U(val)==response(a,x,y),"high-weight interpolation");
  }
  if(index)std::cout<<",";std::cout<<"{\"n\":"<<n<<",\"code\":"<<code<<",\"sidecar_aut\":"<<aut<<",\"queries\":"<<table.size()+9<<",\"S\":";array_json(a.s);
  std::cout<<",\"E\":";array_json(a.e);std::cout<<",\"M\":";array_json(a.m);std::cout<<",\"C\":";array_json(a.c);std::cout<<",\"U\":";array_json(a.u);std::cout<<",\"V\":";array_json(a.v);std::cout<<",\"D\":";array_json(a.d);
  std::cout<<",\"battery\":";array_json(table);std::cout<<"}";std::cerr<<"tensor control complete n="<<n<<" code="<<code<<"\n";
 }
 std::cout<<"],\"status\":\"PASS\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
