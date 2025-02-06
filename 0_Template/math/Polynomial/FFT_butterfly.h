#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double
#define inf 0x3f3f3f3f
#define mod 1000000007
const ld PI=acos(-1.0);  //取PI
const ll N=1<<22;  //注意这个N
typedef complex<ld> comp;
int main()
{
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int t =1;
    void solve();
    // cin>>t;
    while(t--)solve();
}
//FFT最基础的板子
//FFT+蝴蝶优化
//规避了动态数组带来的时间复杂度
//O(n*logn)
//FFT需要运用到复数，建议手写复数，会更快
//FFT会有精度问题
int mx=0;
int rev[N]={0};
ll p1[N],p2[N];
comp tmp1[N],tmp2[N];

void chang(comp *tmp,int len)  //FFT每次递归到底都有一定的规律，即下标二进制翻转
{
    for(int i=0;i<len;i++) //二进制翻转
    {
        rev[i]=rev[i>>1]>>1;
        if(i&1)rev[i]|=len>>1;
    }
    for(int i=0;i<len;i++)
    {
        if(i<rev[i])swap(tmp[i],tmp[rev[i]]);  //交换
    }
}

void FFT(comp *f,int n,int op)    
{
    chang(f,n);  //先变换
    //采用非递归方法求解
    for(int mid=1;mid<n;mid<<=1)   //遍历交换中点，亦是半周期
    {
        comp w(cos(PI/mid),sin(PI*op/mid));   //单位根

        for(int R=mid<<1,j=0;j<n;j+=R)   //周期为R,区间起点j
        {
            comp cur(1,0);  //变化的自变量

            for(int k=0;k<mid;k++,cur*=w)   //傅里叶变换
            {

                comp x=f[j+k],y=cur*f[j+k+mid];  

                f[j+k]=x+y;     //直接覆写

                f[j+k+mid]=x-y;  //直接覆写
            }
        }
    }
}

void solve()
{
    int n,m;
    cin>>n>>m;
    for(int i=0;i<=n;i++)
    {
        cin>>p1[i];
        tmp1[i]=p1[i];
    }
    for(int i=0;i<=m;i++)
    {
        cin>>p2[i];
        tmp2[i]=p2[i];
    }
    mx=1;
    while(mx<n+m+1)mx<<=1;

    for(int i=n+1;i<=mx;i++)
    {
        p1[i]=0;  // 高位补上0,保证是2的幂次
        tmp1[i]=p1[i];
    }
    for(int i=m+1;i<=mx;i++)
    {
        p2[i]=0;
        tmp2[i]=p2[i];
    }
    
    FFT(tmp1,mx,1);  //两个fft求出点值式
    FFT(tmp2,mx,1);

    for(int i=0;i<=mx;i++)tmp1[i]=tmp1[i]*tmp2[i];  //值相乘

    FFT(tmp1,mx,-1);  //乘完后用fft转化为系数式

    for(int i=0;i<n+m+1;i++)cout<<(ll)((ld)tmp1[i].real()/(ld)mx+0.5)<<" ";  //注意实部取四舍五入
    
    cout<<"\n";
}