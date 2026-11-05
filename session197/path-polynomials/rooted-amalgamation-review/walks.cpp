#include <bits/stdc++.h>
using namespace std;
using Z=__int128_t;
string str(Z v){if(!v)return "0";string s;bool neg=v<0;if(neg)v=-v;while(v){s+=char('0'+v%10);v/=10;}if(neg)s+='-';reverse(s.begin(),s.end());return s;}
int main(int argc,char**argv){if(argc!=2)throw runtime_error("graphs input required");ifstream in(argv[1]);int n,checked=0;long long wantH,wantC;while(in>>n>>wantH>>wantC){vector<unsigned>a(n);for(auto &x:a)in>>x;Z h=0,c=0;unsigned full=(1u<<n)-1;
for(unsigned mask=1;mask<=full;mask++){vector<Z>v(n),u(n);for(int i=0;i<n;i++)v[i]=(mask>>i)&1;u[0]=(mask&1);for(int step=0;step<n-1;step++){vector<Z>w(n),t(n);for(int i=0;i<n;i++){unsigned out=a[i]&mask;while(out){int j=__builtin_ctz(out);out&=out-1;w[j]+=v[i];t[j]+=u[i];}}v.swap(w);u.swap(t);}Z hp=0,cy=0;for(int i=0;i<n;i++){hp+=v[i];if(a[i]&1)cy+=u[i];}int sign=((n-__builtin_popcount(mask))&1)?-1:1;h+=sign*hp;c+=sign*cy;}
if(h!=wantH||c!=wantC)throw runtime_error("walk inclusion-exclusion mismatch "+str(h)+","+str(c));cout<<n<<" "<<str(h)<<" "<<str(c)<<"\n";checked++;}
cout<<"PASS "<<checked<<" independent walk-inclusion-exclusion graphs\n";}
