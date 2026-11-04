#define BIRTH_LIBRARY
#include "two_births.cpp"
int main(){int n;while(cin>>n){if(n<1||n>18)return 2;T t(n);for(auto&x:t)cin>>x;auto s=stats(t);auto b=births(t);I best=INT64_MAX;int h2=0,h4=0;for(I x:b){if(x>s.h)best=min(best,x-s.h);h2+=x==s.h+2;h4+=x==s.h+4;}cout<<"{\"H\":"<<s.h<<",\"minimum_increment\":"<<best<<",\"plus2\":"<<h2<<",\"plus4\":"<<h4<<",\"E2\":[";for(int i=0;i<n;i++){if(i)cout<<',';cout<<s.E[i];}cout<<"]}\n"<<flush;}return 0;}
