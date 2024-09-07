#import "@preview/tablex:0.0.6": tablex, hlinex, vlinex, colspanx, rowspanx

#import "@preview/codelst:2.0.1": sourcecode
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
  // flipped: true,
  paper: "a4",
//   background: [#image("background.png")]
)
#set text(
    font:("Times New Roman","Source Han Serif SC"),
    style:"normal",
    weight: "regular",
    size: 13pt,
)


#outline()
#set heading(numbering: "1.")

#let nxtIdx(name) = box[ #counter(name).step()#counter(name).display()]

// Display block code in a larger block
// // with more padding.
// #show raw.where(block: true): block.with(
//   fill: luma(230),
//   inset: 7pt,
//   radius: 4pt,
// )
#set par(
  // first-line-indent: 1cm
)
#set math.equation(numbering: "(1)")



// Display block code in a larger block
// with more padding.
#show raw.where(block: true): block.with(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
)

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
        #smallcaps[Templetes]
        #h(1fr)#text("Github GYPpro/Acm_res",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)

