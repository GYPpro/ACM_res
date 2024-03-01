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
#set par(
  first-line-indent: 1cm
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