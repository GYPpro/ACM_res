//V1
//Dinic算法，求解最大流
//O（M*N^2）
//是对EK算法的优化
//考虑到BFS每次都只能找一条路，思考可不可以一口气找多条路，于是就有了Dinic算法
//先用BFS求分层图，操作和EK算法一样，但不计算，只分层
//然后用DFS在分层图上找增广路，并计算流量，更新残量图
//其中，运用了先前弧剪枝，体现在now数组上，即改变遍历链式前向星的起点，找过的节点不重复寻找
const ll M=1000000000000;
int n,m,s,t,si;
ll sum=0;
struct edge
{
    int to;
    int nex;
    ll w;
}a[500005];
int head[500005]={0},cnt=1;  //链式前向星相关   
void add(int u,int v,ll w)  //链式前向星
{
    //a[++cnt].from=u;
    a[++cnt].to=v;
    a[cnt].w=w;
    a[cnt].nex=head[u];
    head[u]=cnt;

    //a[++cnt].from=v;   //regard edge
    a[++cnt].to=u;
    a[cnt].w=0;
    a[cnt].nex=head[v];
    head[v]=cnt;
}
ll dis[500005]={0};
int now[500005]={0};

int bfs() //在残量网络中构造分层图
{     //事实上，深度标记dis[i]就是分层图，无穷深代表不在图内
    for(int i=1;i<=n;i++)dis[i]=M;//重置分层图
    queue<int> q;
    dis[s]=0;     //标记源点为0深度
    now[s]=head[s];
    q.push(s);
    while(!q.empty())
    {
        int temp=q.front();
        q.pop();
        for(int i=head[temp];i;i=a[i].nex)
        {
            int v=a[i].to;
            if(dis[v]!=M||a[i].w<=0)continue;
            q.push(v);
            dis[v]=dis[temp]+1;  //分层，标记深度


            now[v]=head[v];     //若点在分层图内，就给它的当前弧数组赋值，表示它可以找儿子

            if(v==t)return 1;    //找到汇点就可退出，因为已经标好汇点的深度了
        }
    }
    return 0;
}
ll dfs(int x,ll delta)   //深搜求解 父节点->子节点  的流量
{                        //源点的父节点->源点  就是最大流

    if(x==t)return delta;     //delta最终会是某条路的最小权  
    ll k,res=0;  

    for(int i=now[x];i&&delta;i=a[i].nex)  //从当前弧出发找边，当delta为0时，表示流量到达上限
    {
        now[x]=i;//更新

        int v=a[i].to;//取儿子
        if(a[i].w<=0||dis[v]!=dis[x]+1)continue;//跳过0权边和非法边，即深度不+1的边
        k=dfs(v,min(a[i].w,delta));  //求出x->v的流量

        if(k==0)dis[v]=M;   //发现流量x->v 等于0，删点

        a[i].w-=k;       //更新容量
        a[i^1].w+=k;

        res+=k;       //更新总和
        delta-=k;    //更新剩余容量
    }
    return res;   //返回流入点x的流量
}
void Dinic()
{
    while(bfs())  //bfs构建分层图
    {
        sum+=dfs(s,M);    //dfs找增广路并计算
    }
}
void solve()
{
    cin>>n>>m>>s>>t;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        ll w;
        cin>>u>>v>>w;
        add(u,v,w);
    }
    while(bfs())  //bfs构建分层图
    {
        sum+=dfs(s,M);    //dfs找增广路并计算
    }
    cout<<sum<<"\n";
}

//V2
class maxFlow//AC
{
private:
    class edge
    {
    public:
        ll int nxt, cap, flow;                  
        pair<int, int> revNodeIdx; 
    public:
        edge(int _nxt, int _cap)
        {
            nxt = _nxt,cap = _cap,flow = 0;
        }
        void setRevIdx(int _i, int _j) {  revNodeIdx = {_i,_j}; }
    };
    vector<vector<edge>> edgeList; 
    vector<int> dep,fir;
    ll int maxFlowAns;
    int T, S;
public:
    maxFlow(int _n)
    {
        maxFlowAns = 0;
        S = 1;
        T = _n;
        edgeList.resize(_n + 1);
        dep.resize(_n + 1);
        fir.resize(_n+1);
    }
    void resetTS(int _T, int _S) { T = _T,S = _S; }

    void addedge(int _u, int _v, int _w)
    {
        edgeList[_u].push_back(edge(_v, _w));
        edgeList[_v].push_back(edge(_u, 0)); 
        edgeList[_u][edgeList[_u].size() - 1].setRevIdx(_v, edgeList[_v].size() - 1);
        edgeList[_v][edgeList[_v].size() - 1].setRevIdx(_u, edgeList[_u].size() - 1);
    }

    bool bfs()
    {
        queue<int> que;
        for (auto x = dep.begin(); x != dep.end(); x++) (*x) = 0; 
        dep[S] = 1;
        que.push(S);
        while (que.size())
        {
            ll int at = que.front();
            que.pop();
            for (int i = 0; i < edgeList[at].size(); i++)
            {
                auto tar = edgeList[at][i];
                if ((!dep[tar.nxt]) && (tar.flow < tar.cap))
                {
                    dep[tar.nxt] = dep[at] + 1;
                    que.push(tar.nxt);
                }
            }
        }
        return dep[T];
    }

    ll int dfs(int at, ll int flow)
    {
        if ((at == T) || (!flow))
            return flow; // 到了或者没了
        ll int ret = 0;  // 本节点最大流
        for (int &i = fir[at]; i < edgeList[at].size(); i++)
        {
            auto tar = edgeList[at][i];   // 目前遍历的边
            int tlFlow = 0;   // 目前边的最大流
            if (dep[at] == dep[tar.nxt] - 1) // 遍历到的边为合法边
            {
                tlFlow = dfs(tar.nxt, min((ll)tar.cap - tar.flow, flow - ret));
                if (!tlFlow)
                    continue;   // 若最大流为0，什么都不做
                ret += tlFlow;   // 本节点最大流累加
                edgeList[at][i].flow += tlFlow;  // 本节点实时流量
                edgeList[tar.revNodeIdx.first][tar.revNodeIdx.second].flow -= tlFlow; // 反向边流量
                if (ret == flow)
                    return ret; // 充满了就不继续扫了
            }
        }
        return ret;
    }

    ll int dinic()
    {
        if (maxFlowAns)
            return maxFlowAns;
        while (bfs())
        {
            for(auto x = fir.begin();x != fir.end();x ++) (*x) = 0;
            maxFlowAns += dfs(S, INT_MAX);
        }
        return maxFlowAns;
    }
};
