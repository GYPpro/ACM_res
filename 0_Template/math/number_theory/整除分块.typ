#import "@preview/codelst:2.0.1": sourcecode
// Display inline code in a box
#set text(font:("Times New Roman","Source Han Serif SC"))
#show raw.where(block: false): box.with(
  fill: luma(230),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
)
#show raw.where(block: true): block.with(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
)
#show raw: set text(
    font: ("consolas", "Source Han Serif SC")
  )
#set page(
  // flipped: true,
  // background: [#image("background.png")]
  paper: "a4",
)
#set text(
    font:("Times New Roman","Source Han Serif SC"),
    style:"normal",
    weight: "regular",
    size: 13pt,
)
#show math.equation:set text(font:("New Computer Modern Math","Source Han Serif SC"))
#let nxtIdx(name) = box[ #counter(name).step()#counter(name).display()]
#set math.equation(numbering: "(1)")

#set page(
  paper:"a4",
  number-align: right,
  margin: (x:2cm,y:2.5cm),
  header: [
    #box(baseline:5pt)[#set text(
      size: 11pt,
    )
    #align(
      left+bottom, 
      [
        #smallcaps[ ]
        #h(1fr)#text(" ",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)
#set math.mat(delim: "[")
#set math.vec(delim: "[")

#set page(
  paper:"a4",
  number-align: right,
  margin: (x:2cm,y:2.5cm),
  header: [
    #box(baseline:5pt)[#set text(
      size: 11pt,
    )
    #align(
      left+bottom,
      [
        #smallcaps[Templetes]
        #h(1fr)#text("Github GYPpro/Acm_res",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)
///MAIN---MAIN///

#figure(caption: "整除分块.png")[#image("1.png")]

#sourcecode()[
```cpp
#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define inf 100000000000000
const int N=2e6+7;
//题目，在1<=a<b<=n的条件下，求gcd(a,b)的和

//这里用到了欧拉函数
//欧拉函数也可用欧拉筛求出
int main()
{
    void solve();
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    solve();
    return 0;
}
bool isprime[N];
vector<ll> p;
ll phi[N]={0,1};//边界条件
void eular(int n)
{
    for(int i=2;i<=n;i++)
    {
        if(!isprime[i])
        {
            p.push_back(i);
            phi[i]=i-1;
        }
        for(auto re:p)
        {
            if(i*re>n)break;
            isprime[i*re]=1;
            if(i%re==0){phi[i*re]=phi[i]*re;break;}
            phi[i*re]=phi[i]*(re-1);
        }
    }
    //到此，欧拉函数就求出来了
    for(int i=1;i<=n;i++)phi[i]+=phi[i-1];  //此处在求函数的前缀和，用于整除分块
}


ll cal(ll n,ll m)
{
    ll l=1,r=0,ans=0;
    while(l<=n)
    {
        r=min((n/(n/l)),(m/(m/l))); 
        ans+=(phi[r]-phi[l-1])*(n/l)*(m/l);
        l=r+1;
    }
    return ans;
}
void solve()
{
    ll n;
    cin>>n;
    eular(n);
    ll ans=(cal(n,n)-n*(n+1)/2)/2;
    cout<<ans<<"\n";
}
```
]