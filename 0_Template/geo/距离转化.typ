#import "@preview/codelst:2.0.1": sourcecode
// Display inline code in a box
#set text(font:("Times New Roman","Source Han Serif SC"))
#show raw.where(block: false): box.with(
  fill: luma(230),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
)
#show raw.where(block: true): block.with(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
)
#show raw: set text(
    font: ("consolas", "Source Han Serif SC")
  )
#set page(
  // flipped: true,
  // background: [#image("background.png")]
  paper: "a4",
)
#set text(
    font:("Times New Roman","Source Han Serif SC"),
    style:"normal",
    weight: "regular",
    size: 13pt,
)
#show math.equation:set text(font:("New Computer Modern Math","Source Han Serif SC"))
#let nxtIdx(name) = box[ #counter(name).step()#counter(name).display()]
#set math.equation(numbering: "(1)")

#set page(
  paper:"a4",
  number-align: right,
  margin: (x:2cm,y:2.5cm),
  header: [
    #box(baseline:5pt)[#set text(
      size: 11pt,
    )
    #align(
      left+bottom, 
      [
        #smallcaps[ ]
        #h(1fr)#text(" ",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)
#set math.mat(delim: "[")
#set math.vec(delim: "[")

#set page(
  paper:"a4",
  number-align: right,
  margin: (x:2cm,y:2.5cm),
  header: [
    #box(baseline:5pt)[#set text(
      size: 11pt,
    )
    #align(
      left+bottom,
      [
        #smallcaps[Templetes]
        #h(1fr)#text("Github GYPpro/Acm_res",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)
///MAIN---MAIN///
==== 曼哈顿距离
$ d(A,B) = |x_1 - x_2| + |y_1 - y_2| $

==== 欧几里得距离
$ d(A,B) = sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2) $

==== 切比雪夫距离
$ d(A,B) = max(|x_1 - x_2|, |y_1 - y_2|) $

==== 闵可夫斯基距离
$ d(A,B) = (|x_1 - x_2|^p + |y_1 - y_2|^p)^{1/p} $

==== 曼哈顿转切比雪夫

对于直角坐标中的$A(x_1,y_1),B(x_2,y_2)$

其曼哈顿距离
$ d(A,B) = max(|(x_1+y_1) - (x_2+y_2)|,|(x_1-y_1)-(x_2-y_2|)) $
即为点$A'(x_1+y_1,x_1-y_1),B'(x_2+y_2,x_2-y_2)$的切比雪夫距离。

同理，其切比雪夫距离
$ d(A,B) = max(|(x_1+y_1)/2-(x_2+y_2)/2| + |(x_1-y_1)/2-(x_2-y_2)/2|) $
即为点$A'((x_1+y_1)/2,(x_1-y_1)/2),B'((x_2+y_2)/2, (x_2-y_2)/2)$的曼哈顿距离。

综上：

$
"曼哈顿距离" & =>"切比雪夫距离：" \

(x,y) & => (x+y,x-y) \

"切比雪夫距离"&=>"曼哈顿距离："\

(x,y) &=> ((x+y)/2,(x-y)/2) $