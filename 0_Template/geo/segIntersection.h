using pii = pair<int, int>

#define fi first
#define se second
    const long double EPS = 1e-9;

template <class T>
int sign(T x)
{
    if (-EPS < x && x < EPS)
        return 0;
    return x < 0 ? -1 : 1;
}

// 叉乘
template <class T>
T cross(pair<T, T> a, pair<T, T> b)
{
    return a.fi * b.se - a.se * b.fi;
}

// 二维快速跨立实验
template <class T>
bool segIntersection(pair<T, T> l1, pair<T, T> l2)
{
    auto [s1, e1] = l1;
    auto [s2, e2] = l2;
    auto A = max(s1.fi, e1.fi), AA = min(s1.fi, e1.fi);
    auto B = max(s1.se, e1.se), BB = min(s1.se, e1.se);
    auto C = max(s2.fi, e2.fi), CC = min(s2.fi, e2.fi);
    auto D = max(s2.se, e2.se), DD = min(s2.se, e2.se);
    return A >= CC && B >= DD && C >= AA && D >= BB &&
           sign(cross(s1, s2, e1) * cross(s1, e1, e2)) == 1 &&
           sign(cross(s2, s1, e2) * cross(s2, e2, e1)) == 1;
}

//三维线段交点，需要P3封装，不相交返回{0,{}}
pair<bool, P3> lineIntersection(L3 l1, L3 l2)
{
    if (!onPlane(l1.a, l1.b, l2.a, l2.b) || lineParallel(l1, l2))
    {
        return {0, {}};
    }
    auto [s1, e1] = l1;
    auto [s2, e2] = l2;
    ld val = 0;
    if (!onPlane(l1.a, l1.b, {0, 0, 0}, {0, 0, 1}))
    {
        val = ((s1.x - s2.x) * (s2.y - e2.y) - (s1.y - s2.y) * (s2.x - e2.x)) /
              ((s1.x - e1.x) * (s2.y - e2.y) - (s1.y - e1.y) * (s2.x - e2.x));
    }
    else if (!onPlane(l1.a, l1.b, {0, 0, 0}, {0, 1, 0}))
    {
        val = ((s1.x - s2.x) * (s2.z - e2.z) - (s1.z - s2.z) * (s2.x - e2.x)) /
              ((s1.x - e1.x) * (s2.z - e2.z) - (s1.z - e1.z) * (s2.x - e2.x));
    }
    else
    {
        val = ((s1.y - s2.y) * (s2.z - e2.z) - (s1.z - s2.z) * (s2.y - e2.y)) /
              ((s1.y - e1.y) * (s2.z - e2.z) - (s1.z - e1.z) * (s2.y - e2.y));
    }
    return {1, s1 + (e1 - s1) * val};
}
