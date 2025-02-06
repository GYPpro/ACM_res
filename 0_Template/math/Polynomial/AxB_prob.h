#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define ld double
#define inf 10000000
#define mod 998244353
const ld pi=acos(-1);
signed main()
{
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    void solve();
    solve();
    return 0;
}
//高精乘法
//但是FFT
//注意只能FFT，不能NTT,NTT会乱
//只需将每一个位置上的数视作系数，对应一个10的k次方即可
const int N=(1<<22)|10;
string s1,s2;
int n,m;
int rev[N]={0};
struct comp    //手搓复数
{
    ld r=0,i=0;
}a[N],b[N];
ll c[N];
comp mul(comp x,comp y)
{
    comp res;
    res.r=x.r*y.r-x.i*y.i;
    res.i=x.r*y.i+x.i*y.r;
    return res;
}
comp add(comp x,comp y,int op)
{
    comp res;
    res.r=x.r+op*y.r;
    res.i=x.i+op*y.i;
    return res;
}

void chang(comp *A,int n)
{
    for(int i=0;i<n;i++)rev[i]=(rev[i>>1]>>1)|((i&1)?(n>>1):0);
    for(int i=0;i<n;i++)if(i<rev[i])swap(A[i],A[rev[i]]);
}
void FFT(comp *A,int n,int op)
{
    chang(A,n);
    for(int mid=1;mid<n;mid<<=1)
    {
        comp g1;
        g1.r=cos(pi/(ld)mid);
        g1.i=sin(op*pi/(ld)mid);
        for(int j=0;j<n;j+=(mid<<1))
        {
            comp gk;
            gk.r=1;gk.i=0;
            for(int k=0;k<mid;k++,gk=mul(gk,g1))
            {
                comp x=A[j+k],y=mul(A[j+k+mid],gk);
                A[j+k]=add(x,y,1);
                A[j+k+mid]=add(x,y,-1);
            }
        }
    }
    if(op==1)return ;
    for(int i=0;i<n;i++)c[i]=(ll)(A[i].r/n+0.5);
}

void solve()
{
    stack<int> ans;
    cin>>s1;
    cin>>s2;
    n=s1.length();
    m=s2.length();
    for(int i=0;i<n;i++)a[n-1-i].r=s1[i]-'0',a[n-1-i].i=0;
    for(int i=0;i<m;i++)b[m-1-i].r=s2[i]-'0',b[n-1-i].i=0;
    int mx=1;
    while(mx<n+m)mx<<=1;
    for(int i=n;i<mx;i++)a[i].i=0,a[i].r=0;
    for(int i=m;i<mx;i++)b[i].i=0,b[i].r=0;
    FFT(a,mx,1);FFT(b,mx,1);
    for(int i=0;i<mx;i++)a[i]=mul(a[i],b[i]);
    FFT(a,mx,-1);
    for(int i=0;i<mx;i++)
    {
        if(c[i]>=10)
        {
            c[i+1]+=(ll)c[i]/10;
            c[i]%=10;
        }
        ans.push(c[i]);
    }
    bool flag=true;
    while(!ans.empty())
    {
        int temp=ans.top();
        ans.pop();
        if(temp==0&&flag)continue;
        flag=false;
        cout<<temp;
    }
    cout<<"\n";
}