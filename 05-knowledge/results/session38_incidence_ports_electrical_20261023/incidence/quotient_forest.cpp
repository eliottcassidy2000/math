// Complete downward forest search through quotients of the fixed nine E33 symbols.
// Distinct class edges are never merged; signature equality is exact string equality.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
struct Data {int conv; array<vector<int>,8> rows;};
struct Result {bool ok; int signatures; int first_cycle;};
struct UF {vector<int> p,s; void add(){int v=p.size();p.push_back(v);s.push_back(1);} int find(int a){while(p[a]!=a){p[a]=p[p[a]];a=p[a];}return a;} bool join(int a,int b){a=find(a);b=find(b);if(a==b)return false;if(s[a]<s[b])swap(a,b);p[b]=a;s[a]+=s[b];return true;}};
string signature(const Data& d,const string& part){array<string,8> rows;for(int u=0;u<8;u++){string r;for(int x:d.rows[u])r+=part[x];sort(r.begin(),r.end());rows[u]=char('0'+r.size())+r;}sort(rows.begin(),rows.end());string r;for(auto &x:rows)r+=x;return r;}
Result test(const vector<Data>& D,const string& part){unordered_map<string,int> ids;ids.reserve(14000);vector<int> colors(D.size(),-1);UF uf;auto get=[&](int i){if(colors[i]>=0)return colors[i];string s=signature(D[i],part);auto it=ids.find(s);if(it!=ids.end())return colors[i]=it->second;int x=ids.size();ids.emplace(move(s),x);uf.add();uf.add();return colors[i]=x;};for(int i=0;i<(int)D.size();i++){int a=get(i),b=get(D[i].conv);if(!uf.join(2*a,2*b+1))return {false,(int)ids.size(),i};}return {true,(int)ids.size(),-1};}
string mergepart(string p,char a,char b){for(char &x:p)if(x==b)x=a;array<int,9> m;fill(m.begin(),m.end(),-1);int z=0;for(char &x:p){int q=x-'0';if(m[q]<0)m[q]=z++;x='0'+m[q];}return p;}
int main(int argc,char**argv){if(argc!=2)throw runtime_error("usage: quotient_forest prepared_input.txt");ifstream f(argv[1]);int n;f>>n;if(n!=6880)throw runtime_error("wrong class count");vector<Data>D(n);for(int i=0;i<n;i++){string row;f>>D[i].conv>>row;if(row.size()!=36)throw runtime_error("bad row bytes");int at=0;for(int u=0;u<8;u++){int k=row.at(at++)-'0';for(int j=0;j<k;j++){int x=row.at(at++)-'0';if(x<0||x>8)throw runtime_error("bad symbol");D[i].rows[u].push_back(x);}}if(at!=36)throw runtime_error("bad parse");}for(int i=0;i<n;i++)if(D.at(D[i].conv).conv!=i)throw runtime_error("bad converse");set<string> parents={"012345678"};cout<<"level\tpartition\tforest\tsignatures_seen\tfirst_cycle_class\n";for(int k=8;k>=1;k--){set<string> parts,good;for(auto p:parents)for(int a=0;a<=k;a++)for(int b=a+1;b<=k;b++)parts.insert(mergepart(p,'0'+a,'0'+b));for(auto p:parts){auto r=test(D,p);cout<<k<<'\t'<<p<<'\t'<<r.ok<<'\t'<<r.signatures<<'\t'<<r.first_cycle<<'\n';if(r.ok)good.insert(p);}cerr<<"level "<<k<<" tested "<<parts.size()<<" passed "<<good.size()<<endl;cout.flush();if(good.empty())break;parents=move(good);}}
