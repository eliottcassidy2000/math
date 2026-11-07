#include <bits/stdc++.h>
using namespace std; using I=long long; using U=uint64_t;
struct S { I h,q,a,g,t,s,e; bool operator==(const S&o)const{return tie(h,q,a,g,t,s,e)==tie(o.h,o.q,o.a,o.g,o.t,o.s,o.e);} };
vector<vector<int>> cores={{8,1,3,6},{2,4,9,19,7},{2,8,3,21,7}};
S literal(const vector<int>& rows,int cut){
 int n=rows.size(), M=1<<n; vector<I> end(M*n),start(M*n),h(M),e(M),s(M);h[0]=1;
 for(int v=0;v<n;v++)end[(1<<v)*n+v]=start[(1<<v)*n+v]=1;
 for(int mask=1;mask<M;mask++)for(int v=0;v<n;v++)if(mask>>v&1){
  I a=end[mask*n+v],b=start[mask*n+v];h[mask]+=a;if(cut>>v&1){e[mask]+=a;s[mask]+=b;}
  for(int w=0;w<n;w++)if(!(mask>>w&1)){
   if(rows[v]>>w&1)end[(mask|(1<<w))*n+w]+=a;
   if(rows[w]>>v&1)start[(mask|(1<<w))*n+w]+=b;
  }
 }
 S z{h[M-1],0,0,0,0,s[M-1],e[M-1]};
 for(int A=0;A<M;A++){int B=(M-1)^A;z.q+=h[A]*h[B];z.a+=h[A]*s[B];z.g+=e[A]*h[B];z.t+=e[A]*s[B];}
 return z;
}
S prepend(S z,int bit){return {z.h,2*z.q,z.a+bit*z.q,2*z.g+bit*z.h,z.t+bit*(z.g+z.s),bit*z.h,z.e};}
S append(S z,int bit){return {z.h,2*z.q,2*z.a+bit*z.h,z.g+bit*z.q,z.t+bit*(z.a+z.e),z.s,bit*z.h};}
I fword(int z){I out=0;for(int i=0;i<16;i++)if(z>>i&1)for(int j=i+1;j<16;j++)if(z>>j&1)out+=(1LL<<(j-i-1))+(j==i+1);return out;}
I closed(S z,int p,int r,int b){int rev=0;for(int i=0;i<b;i++)rev=2*rev+(r>>i&1);return (z.q*p+z.a)*((1LL<<(b+1))-rev)+z.h*(2*r+1-fword(r))-p*z.g-z.h*fword(p)-z.t-(p%2)*z.s-z.e*(b?(r>>(b-1)&1):0);}
vector<int> padded(int ci,int a,int b){auto c=cores[ci];int k=c.size(),n=a+k+b;vector<int>r(n);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(i>=a&&j<a+k){if(c[i-a]>>(j-a)&1)r[i]|=1<<j;else r[j]|=1<<i;}else r[i]|=1<<j;return r;}
int fullcut(int ci,int cut,int a,int b,int p,int r){int k=cores[ci].size();int out=cut<<a;for(int i=0;i<a;i++)if(p>>i&1)out|=1<<(a-1-i);for(int i=0;i<b;i++)if(r>>i&1)out|=1<<(a+k+b-1-i);return out;}
vector<int> close(vector<int>r,int cut){int n=r.size();r.resize(n+2);for(int i=0;i<n;i++){r[i]|=1<<n;if(cut>>i&1)r[n+1]|=1<<i;else r[i]|=1<<(n+1);}r[n]|=1<<(n+1);return r;}
I paths(const vector<int>&r){int n=r.size(),M=1<<n;vector<I>d((size_t)M*n);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;for(int mask=1;mask<M;mask++)for(int i=0;i<n;i++)if(mask>>i&1){int next=r[i]&~mask;while(next){int bit=next&-next,j=__builtin_ctz((unsigned)bit);next-=bit;d[(mask|bit)*n+j]+=d[mask*n+i];}}I h=0;for(int i=0;i<n;i++)h+=d[(M-1)*n+i];return h;}
bool strong(const vector<int>&r){int n=r.size();for(int s=0;s<n;s++){int seen=1<<s,last=0;while(seen!=last){last=seen;for(int i=0;i<n;i++)if(seen>>i&1)seen|=r[i];}if(seen!=(1<<n)-1)return false;}return true;}
int main(int argc,char**argv){if(argc!=3)return 2;ifstream st(argv[1]);vector<S>C;vector<pair<int,int>>names;int statechecks=0,dpchecks=0,hostile=0;set<I> smallstrong;
 for(int ci=0;ci<3;ci++){if(!strong(cores[ci]))return 3;for(int cut=0;cut<(1<<(int)cores[ci].size());cut++){S z=literal(cores[ci],cut),w;st>>w.h>>w.q>>w.a>>w.g>>w.t>>w.s>>w.e;if(!(z==w))return 4;C.push_back(z);names.emplace_back(ci,cut);
 for(int a=0;a<=3;a++)for(int b=0;a+b<=3;b++)for(int p=0;p<(1<<a);p++)for(int r=0;r<(1<<b);r++){
 S u=z;for(int i=0;i<a;i++)u=prepend(u,p>>i&1);for(int i=b-1;i>=0;i--)u=append(u,r>>i&1);
 auto L=padded(ci,a,b);int sc=fullcut(ci,cut,a,b,p,r);if(!(u==literal(L,sc)))return 5;statechecks++;
 auto G=close(L,sc);I h=paths(G);if(h!=2*u.a+u.h-u.t||h!=closed(z,p,r,b))return 6;dpchecks++;if(u.t&&h!=2*u.a+u.h)hostile++;
 bool expected=a?(p>>(a-1)&1):cut!=0;if(strong(G)!=expected)return 7;if(expected)smallstrong.insert(h);
 }}}
 if(smallstrong.count(47)||smallstrong.count(63))return 14;
 ifstream in(argv[2],ios::binary);int checked=0,maxorder=0;vector<int>missing;map<int,int>bycore;for(int index=0;index<=500000;index++){U value=0;for(int j=0;j<8;j++){int v=in.get();if(v==EOF)return 8;value|=U(v)<<(8*j);}int target=2*index+1;if(!value){if(target>=23&&target<=1000000)missing.push_back(target);continue;}U pack=value-1;int c=pack&127,a=(pack>>7)&31,b=(pack>>12)&31,p=(pack>>17)&65535,r=(pack>>33)&65535;if(c>=80||a+b>16||p>=(1<<a)||r>=(1<<b)||(pack>>49))return 9;auto[ci,cut]=names[c];if(a?!(p>>(a-1)&1):!cut)return 10;S z=C[c];for(int i=0;i<a;i++)z=prepend(z,p>>i&1);for(int i=b-1;i>=0;i--)z=append(z,r>>i&1);if(2*z.a+z.h-z.t!=target||closed(C[c],p,r,b)!=target)return 11;checked++;maxorder=max(maxorder,(int)cores[ci].size()+a+b+2);bycore[ci]++;}
 if(in.get()!=EOF)return 12;
 vector<bool> hits(1000001);I scanned=0;
 function<void(S,int)> suffix=[&](S z,int remaining){I v=2*z.a+z.h-z.t;scanned++;if(v>=23&&v<=1000000){if(!(v&1))throw runtime_error("even");hits[v]=true;}if(remaining){suffix(append(z,0),remaining-1);suffix(append(z,1),remaining-1);}};
 for(int c=0;c<80;c++){
  function<void(S,int,int)> prefix=[&](S z,int depth,int last){if((depth&&last)||(!depth&&names[c].second))suffix(z,16-depth);if(depth<16){prefix(prepend(z,0),depth+1,0);prefix(prepend(z,1),depth+1,1);}};
  prefix(C[c],0,0);
 }
 vector<int> scanmissing;for(int v=23;v<=1000000;v+=2)if(!hits[v])scanmissing.push_back(v);if(scanmissing!=missing)return 13;
 cout<<"{\"independent_scan_count\":"<<scanned<<",\"scan_agrees_with_all_recipes\":true}"<<endl;
 cout<<"{\"state_controls\":"<<statechecks<<",\"literal_closure_controls\":"<<dpchecks<<",\"omitted_tau_failures\":"<<hostile<<",\"verified_recipes\":"<<checked<<",\"max_order\":"<<maxorder<<",\"missing\":[";for(int i=0;i<(int)missing.size();i++){if(i)cout<<",";cout<<missing[i];}cout<<"],\"core_counts\":["<<bycore[0]<<","<<bycore[1]<<","<<bycore[2]<<"]}"<<endl;
}
