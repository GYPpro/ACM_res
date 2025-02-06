#include<bits/stdc++.h>
using namespace std;
#define ll long long
const ll N=250000+7;
int main()
{
	// ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
	void solve();
	int t=1;
	while(t--)solve();
	return 0;
}
//无向图
//判断欧拉回路和通路
map<string,int> point;
int cnt=0,tot=0;
struct node{int u,v,i;};
vector<node> edge;
vector<int> e[N<<1];
int ind[N<<1],vis[N];
void yes(){cout<<"Possible\n";}
void no(){cout<<"Impossible\n";}
int dfs(int u)  //深搜找欧拉回路/通路
{
    int ans=0;
    for(auto re:e[u])
    {
        int v=edge[re].u^edge[re].v^u;
        if(!vis[re])continue;
        vis[re]--;
        ans++;
        ans+=dfs(v);
    }
    return ans;   //这里在数经过的边数
    //如果要找路，可以在这将点加入栈
} 
bool check(int m)
{
    for(int i=1;i<=m;i++)  //建图、数度数
    {
        ind[edge[i].u]++;    //有向图则要分出度和入度
        ind[edge[i].v]++;

        vis[i]++;  //每条边只走一次

        e[edge[i].u].push_back(i);  //邻接表存边号
        e[edge[i].v].push_back(i);
    }

    int st=1;   //起点

    int flag=3; //

    for(int i=1;i<=tot;i++)
    {
        if(ind[i]&1)   //数奇数度的点，超过2个直接返回false
        {               //有向图这里要判断出度和入度相等
            st=i;
            flag--;
        }
        if(flag==0)break;
    }

    if(!flag)return flag;  //不满足充要条件————恰好两个奇数点，或者没有奇数点
    if(dfs(st)==cnt)return flag;  //满足条件，且能一笔画
    else return 0;
}
/// 以上是核心算法////
void solve()
{
    int op=1;
    char c;
    string s1,s2;
    edge.push_back({-1,-1,0});
    for(int i=1;i<=tot;i++)e[i].clear();
    while(scanf("%c",&c)!=EOF)
    {
        if(c=='0')break;
        if(c==' ')op^=1;
        else if(c=='\n')
        {
            cnt++;
            if(!point[s1])point[s1]=++tot;
            if(!point[s2])point[s2]=++tot;
            edge.push_back({point[s1],point[s2],cnt});
            s1.clear();
            s2.clear();
            op^=1;
        }
        else if(op)
        {
            s1.push_back(c);
        }
        else 
        {
            s2.push_back(c);
        }
    }
    if(check(cnt))yes();
    else no();
}