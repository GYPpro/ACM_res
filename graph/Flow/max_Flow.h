#ifndef _IN_TEMPLATE_
#include <template_overAll.h>
#endif

class maxFlow
{
private:
    class edge
    {
    public:
        ll int nxt,                   // 出度
              cap,                   // 容量
              flow;                  // 流量
        pair<int, int> revNodeIdx; // 反向边
    public:
        edge(int _nxt, int _cap)
        {
            nxt = _nxt;
            cap = _cap;
            flow = 0;
        }
        void setRevIdx(int _i, int _j)
        {
            revNodeIdx.first = _i;
            revNodeIdx.second = _j;
        }
    };
    vector<vector<edge>> edgeList; // 节点列表
    vector<int> dep;               // 深度
    vector<int> fir;               // 节点最近合法边
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

    void resetTS(int _T, int _S)
    {
        T = _T;
        S = _S;
    }

    void addedge(int _u, int _v, int _w)
    {
        edgeList[_u].push_back(edge(_v, _w));
        edgeList[_v].push_back(edge(_u, 0)); // 反向建边
        edgeList[_u][edgeList[_u].size() - 1].setRevIdx(_v, edgeList[_v].size() - 1);
        edgeList[_v][edgeList[_v].size() - 1].setRevIdx(_u, edgeList[_u].size() - 1);
    }

    bool bfs() // 统计深度
    {
        queue<int> que;
        for (auto x = dep.begin(); x != dep.end(); x++)
            (*x) = 0;

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
            auto tar = edgeList[at][i];      // 目前遍历的边
            int tlFlow = 0;                  // 目前边的最大流
            if (dep[at] == dep[tar.nxt] - 1) // 遍历到的边为合法边
            {
                tlFlow = dfs(tar.nxt, min((ll)tar.cap - tar.flow, flow - ret));
                if (!tlFlow)
                    continue;                                                         // 若最大流为0，什么都不做
                ret += tlFlow;                                                        // 本节点最大流累加
                edgeList[at][i].flow += tlFlow;                                       // 本节点实时流量
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
