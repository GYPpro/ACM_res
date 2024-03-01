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
// #set par(
//   first-line-indent: 1cm
// )
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

= Lambda表达式在DFS上的应用

Lambda表达式：形如```cpp
[capture](parameters)->return-type{
    ...
}
```
其中`capture`为闭包行为指定标识符：
```cpp
[]      // 沒有定义任何变量。使用未定义变量会引发错误。
[x, &y] // x以传值方式传入（默认），y以引用方式传入。
[&]     // 任何被使用到的外部变量都隐式地以引用方式加以引用。
[=]     // 任何被使用到的外部变量都隐式地以传值方式加以引用。
[&, x]  // x显式地以传值方式加以引用。其余变量以引用方式加以引用。
[=, &z] // z显式地以引用方式加以引用。其余变量以传值方式加以引用。
```
具体的，使用：
```cpp
auto dfs = [&](auto self,int x,int p) -> void {
    ...
}
```
两个`auto`均解释为`class lambda [](auto self, int x, int y)->void`

调用时`dfs(dfs,x,p);`
