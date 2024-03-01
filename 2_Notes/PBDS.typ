#import "@preview/tablex:0.0.6": tablex, hlinex, vlinex, colspanx, rowspanx

// Display inline code in a small box
// that retains the correct baseline.
#set text(font:("Times New Roman","Source Han Serif SC"))
#show raw.where(block: false): box.with(
  fill: luma(230),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
)
// #set raw(align: center)
#show raw: set text(
    font: ("consolas", "Source Han Serif SC")
  )
#set page(
//   flipped: true,
  paper: "a4",
//   background: [#image("background.png")]
)
#set text(
    font:("Times New Roman","Source Han Serif SC"),
    style:"normal",
    weight: "regular",
    size: 13pt,
)

#let nxtIdx(name) = box[ #counter(name).step()#counter(name).display()]

// Display block code in a larger block
// with more padding.
#show raw.where(block: true): block.with(
  fill: luma(230),
  inset: 7pt,
  radius: 4pt,
)
#set math.equation(numbering: "(1)")

#let titlePage(titleName,translatedTitle,Discrib) = [
    #set page(footer: [])
    #[
        #text(
            font:("Times New Roman","Source Han Serif SC"),
            style:"normal",
            weight:"regular",
            size: 22pt,
        )[
            #align(
                left+horizon
            )[
                #heading(level: 1,[#strong(titleName)])
                #smallcaps(translatedTitle)
                #line(start: (0pt,11pt),end:(300pt,11pt))
                #[
                    #text(Discrib,size:19pt)
                ]
            ]
        ]
    ]
]

#set heading( 
    numbering: "1.1.1."
 )

= PBDS (Policy-Based Data Structures)中的常用结构

== 基于`tree_policy`
\
*头文件:*

```cpp
#include <ext/pb_ds/assoc_container.hpp> 
#include <ext/pb_ds/tree_policy.hpp> 
using namespace __gnu_pbds;
```

*类型定义：*
```cpp
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;

typedef tree<int, null_type, less_equal<int>, rb_tree_tag,
tree_order_statistics_node_update> ordered_multiset;
```

*使用：*

除正常按`set`、`multiset`使用外，有两种特殊成员函数：

```cpp
    ordered_multiset o_mset;
    ordered_set o_set;
    ...
    o_set.find_by_order(n);//返回第n+1小的元素的迭代器
                           //即序列增序时下标n对应元素
    o_mset.find_by_order(n);//同理

    o_set.order_of_key(n);//返回容器中严格小于n的元素的数量

    o_mset.order_of_key(n);//同理

```