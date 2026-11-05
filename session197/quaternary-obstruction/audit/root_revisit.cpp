#include <bits/stdc++.h>
using namespace std;using Z=long long;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
vector<Z> runs(vector<int>a,int r){
 vector<int>v;for(int i=0;i<(int)a.size();i++)if(i!=r)v.push_back(i);
 int n=v.size(),F=(1<<n)-1,K=n/2+1;vector<array<Z,7>>d((1<<n)*(n+1));
 auto at=[&](int S,int last)->array<Z,7>&{return d[S*(n+1)+last];};
 at(0,n)[1]=1;for(int i=0;i<n;i++)at(1<<i,i)[0]=1;
 for(int S=0;S<=F;S++){
  if(S)for(int i=0;i<n;i++)if(S>>i&1&&a[v[i]]>>r&1)for(int j=1;j<=K;j++)at(S,n)[j]+=at(S,i)[j-1];
  for(int next=0;next<n;next++)if(!(S>>next&1)){
   auto&dest=at(S|1<<next,next);
   for(int j=0;j<=K;j++){
    if(a[r]>>v[next]&1)dest[j]+=at(S,n)[j];
    for(int i=0;i<n;i++)if(S>>i&1&&a[v[i]]>>v[next]&1)dest[j]+=at(S,i)[j];
   }
  }
 }
 vector<Z>out(K+1);for(int j=0;j<=K;j++)for(int i=0;i<=n;i++)out[j]+=at(F,i)[j];return out;
}
int main(){string line;long long graphs=0,signs=0,roots=0,second=0,third=0,fact=0;vector<string>survivors;
 while(getline(cin,line)){replace(line.begin(),line.end(),',',' ');istringstream ss(line);int n,h;ss>>n>>h;vector<int>a(n);for(int&v:a)ss>>v;graphs++;
 for(int eps:{-1,1}){if((h-eps)%4)continue;int A=(h-eps)/4;if(A<1||A>54||A%6)continue;signs++;
 for(int r=0;r<n;r++){roots++;auto c=runs(a,r);need(c[1]==h,"a1 original count");Z f=1;for(int j=1;j<(int)c.size();j++){if(j>1)f*=j-1;need(c[j]%f==0,"free internal-fragment permutations");fact++;}
 Z h2=c[1]+(c.size()>2?2*c[2]:0);Z h3=c[1]+(c.size()>2?6*c[2]:0)+(c.size()>3?6*c[3]:0);
 if(h2==16*A+eps){second++;ostringstream z;z<<"{\"n\":"<<n<<",\"H1\":"<<h<<",\"root\":"<<r<<",\"A\":"<<A<<",\"offset\":"<<eps<<",\"H2\":"<<h2<<",\"H3\":"<<h3<<",\"expected_H3\":"<<64*A+eps<<",\"runs\":[";for(int j=1;j<(int)c.size();j++){if(j>1)z<<',';z<<c[j];}z<<"]}";survivors.push_back(z.str());if(h3==64*A+eps)third++;}
 }}}
 need(graphs==2456&&signs==410&&roots==3426&&second==2&&third==0,"complete exclusion counts");
 cout<<"{\"status\":\"PASS\",\"graphs\":"<<graphs<<",\"eligible_graph_signs\":"<<signs<<",\"roots\":"<<roots<<",\"factorial_controls\":"<<fact<<",\"H2_survivors\":[";for(int i=0;i<(int)survivors.size();i++){if(i)cout<<',';cout<<survivors[i];}cout<<"],\"H3_survivors\":"<<third<<"}\n";
}
