
using dot = pair<int,int>;
using lin = pair<dot,dot>;

#define x first
#define y second
#define fi first
#define se second

const int TRF = 1e11;

int cross(dot a,dot b) {
    return a.x * b.y - a.y * b.x;
}

dot dsc(dot a,dot b) {
    return {a.x - b.x,a.y - b.y};
}

dot add(dot a,dot b) {
    return {a.x + b.x,a.y + b.y};
}

int cross(dot p1,dot p2,dot p0) {
    return cross(dsc(p1,p0),dsc(p2,p0));
}

int sign(int x) {
    if(x == 0) return 0;
    return x < 0 ? -1 : 1;
}

bool onseg(lin l,dot p) {
    return sign( cross(p,l.fi,l.se) == 0 ) && 
    (min(l.fi.x,l.se.x) <= p.x && p.x <= max(l.fi.x,l.se.x)) &&
    (min(l.fi.y,l.se.y) <= p.y && p.y <= max(l.fi.y,l.se.y)) ;
};

bool sic(lin a,lin b) {
    auto [s1,e1] = a;
    auto [s2,e2] = b;
    auto A = max(s1.x,e1.x),AA = min(s1.x,e1.x);
    auto B = max(s1.y,e1.y),BB = min(s1.y,e1.y);
    auto C = max(s2.x,e2.x),CC = min(s2.x,e2.x);
    auto D = max(s2.y,e2.y),DD = min(s2.y,e2.y);

    bool flag_cross = (sign(cross(s1,s2,e1)) * sign(cross(s1,e1,e2))) == 1 &&
                      (sign(cross(s2,s1,e2)) * sign(cross(s2,e2,e1))) == 1;
    bool flag_onseg = onseg(a,s2) || onseg(a,e2) ||
                      onseg(a,s2) || onseg(a,e2);

    return A >= CC && B >= DD && C >= AA && D >= BB && (flag_cross || flag_onseg); 
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
