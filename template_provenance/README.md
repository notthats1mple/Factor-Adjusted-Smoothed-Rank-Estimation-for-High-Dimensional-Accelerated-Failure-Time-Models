# Springer SVJour3 provenance and journal settings

Journal rules checked on 2026-10-14: https://link.springer.com/journal/10985/submission-guidelines

The journal explicitly requests the Springer `smallextended` layout and 12-point new submissions. This build uses SVJour3 v3.2 with `smallextended,nospthms`. The accompanying class and global option file were transported as readable source from the public mirror below, with comments/whitespace condensed. They are not claimed to be a byte-identical download of the publisher's current ZIP.

- https://raw.githubusercontent.com/latextemplates/svjour/master/svjour3.cls
- https://raw.githubusercontent.com/latextemplates/svjour/master/svglov3.clo

`lda-review12.sty` is this project's separate, documented font-size override. SVJour3 v3.2 does not implement a `12pt` class option; merely appending that option would leave the body at the default size. The override sets normal text and references to 12 TeX points while retaining the 119 mm smallextended text width. Captions/footnotes are 10 points; plot labels are 8–10 PDF points. Subscripts retain their normal mathematical sizing. The anonymous manuscript-number box and acceptance dates are not fabricated.

The journal's linked macro-package download can be used to replace the two vendor files on the author's connected computer; keep the local override separate and rebuild all outputs after any replacement. No standalone font files are distributed. PDF/EPS resources have embedded fonts.

The publisher's binary ZIP could not be downloaded from the current container's network. The supplied source compiles as tested; this note records that acquisition limit rather than calling the mirror the latest official release.
