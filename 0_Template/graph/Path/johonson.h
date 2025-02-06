// Johnson算法

// 带负边的全源最短路

// 这个算法最重要的功能是使得dj能用来处理负边

// 只跑一次的话复杂度和spfa同阶
// 跑多次的话会比spfa优秀

// 先建立一个虚拟源点o，从该点向所有边连一条边权为0，跑一遍spfa()，记录点o到任意点i的最短路h[i]
// 再将原边权w(u,v),改造为w(u,v)+h[u]-h[v]
// 再以每个点作为起点，跑n轮dj即可

// 最终d'(s,t)=d(s,t)+h[s]-h[t]
#define ll long long
#define inf 0x3f3f3f3f
int main()
{
    void solve();
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    solve();
    return 0;
}
//Johnson算法
//带负边的全源最短路
//先建立一个虚拟源点o，从该点向所有边连一条边权为0，跑一遍spfa()，记录点o到任意点i的最短路h[i]
//再将原边权w(u,v),改造为w(u,v)+h[u]-h[v]
//再以每个点作为起点，跑n轮dj即可
int n,m;
struct p5905
{
    int to;
    ll w;
};
vector<p5905> edge[3005];
ll dis[3005]={0},h[3005]={0};
int cnt[3005]={0};
int vis[3005]={0};
queue<int> q;
void spfa()
{
    memset(h,inf,sizeof(h));
    h[0]=0;
    q.push(0);
    vis[0]++;
    while(!q.empty())
    {
        int t=q.front();
        q.pop();
        cnt[t]++;
        if(cnt[t]>n)
        {
            cout<<-1<<"\n";
            exit(0);
        }
        for(auto re:edge[t])
        {
            int v=re.to;
            ll w=re.w;
            if(h[v]>h[t]+w)
            {
                h[v]=h[t]+w;
                if(!vis[v])
                {
                    vis[v]++;
                    q.push(v);
                }
            }
        }
        vis[t]--;
    }
}
struct node
{
    int id;
    ll w;
    bool operator < (const node &x) const {return x.w<w;}
};
void dj(int f)
{
    memset(dis,inf,sizeof(dis));
    int vv[3005]={0};
    priority_queue<node> qq;
    dis[f]=0;
    qq.push((node){f,dis[f]});
    while(!qq.empty())
    {
        node temp=qq.top();
        qq.pop();
        int t=temp.id;
        ll w=temp.w;
        if(vv[t])continue;
        vv[t]++;
        dis[t]=w;
        for(auto re:edge[t])
        {
            int v=re.to;
            ll tw=re.w;
            if(vv[v])continue;
            if(dis[v]>w+tw)
            {
                dis[v]=w+tw;
                qq.push((node){v,dis[v]});
            }
        }
    }
}
void solve()
{
    cin>>n>>m;
    while(m--)
    {
        int u,v;
        ll w;
        cin>>u>>v>>w;
        edge[u].push_back({v,w});
    }
    for(int i=1;i<=n;i++)
    {
        edge[0].push_back({i,0});
    }
    spfa();
    for(int i=1;i<=n;i++)
    {
        for(auto &re:edge[i])
        {
            int v=re.to;
            re.w+=h[i]-h[v];   //修改边权
        }
    }
    ll ans=0;
    for(int i=1;i<=n;i++)
    {
        dj(i);
        ans=0;
        for(int j=1;j<=n;j++)
        {
            if(dis[j]>1e9)ans+=j*(1e9);
            else 
            {
                ans+=j*(dis[j]+h[j]-h[i]);//这里j*是题目要求
                //真正的i~j的距离是括号内的表达式
            }
        }
        cout<<ans<<"\n";
    }
}