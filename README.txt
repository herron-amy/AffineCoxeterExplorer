╔══════════════════════════════════════════════════════════════╗
║       Affine Coxeter Explorer  —  README                  ║
╚══════════════════════════════════════════════════════════════╝


WHAT IS THIS?
─────────────
A desktop application for computing and visualising:
  • Conjugacy classes
  • Coconjugation sets

in affine Weyl groups of dimensions 2 and 3:

  Dimension 2 (2D plots):   A1xA1,  A2-tilde,  B2-tilde,  C2-tilde,  G2-tilde
  Dimension 3 (3D plots):   A3-tilde,  B3-tilde,  C3-tilde


FILES
─────
  AffineCoxeterExplorer.java     The application
  compute_helper.sage          SageMath computation engine
  Launch (macOS).command       Double-click launcher for macOS
  Launch (Windows).bat         Double-click launcher for Windows
  Launch (Linux).sh            Double-click launcher for Linux
  README.txt                   This file


STEP 1 — INSTALL REQUIREMENTS  (one time only)
───────────────────────────────────────────────
1. SageMath (free)
     https://www.sagemath.org/download.html
     macOS: download the .dmg, drag to Applications

2. Java 11 or newer (free)
     https://adoptium.net
     macOS: download the .pkg for your chip:
       Apple Silicon (M1/M2/M3) → choose aarch64
       Intel Mac               → choose x64
     Not sure which? Run in Terminal:  uname -m
       arm64  = Apple Silicon
       x86_64 = Intel


STEP 2 — LAUNCH THE APP
────────────────────────
macOS
  1. Open Terminal and run (one time only):
       chmod +x "Launch (macOS).command"
  2. Double-click "Launch (macOS).command" in Finder
     First time only: right-click → Open → Open
     (macOS Gatekeeper requires this for downloaded files)

Windows
  Double-click "Launch (Windows).bat"

Linux
  chmod +x "Launch (Linux).sh"
  ./Launch\ \(Linux\).sh

Any platform (from Terminal)
  conda deactivate          ← skip if you don't have Anaconda
  java --source 11 AffineCoxeterExplorer.java


USING THE APP
─────────────
The window has a control panel on the left and a plot area on the right.

  1  Select group       8 affine Weyl groups (dimensions 2 and 3)
  2  Computation        Conjugacy Class  or  Coconjugation Set
  3  Elements           Enter in either of two formats (see below)
                        Leave blank OR enter "e" for the identity
                        Press Enter in any field to run
  4  Bounding box       Slider — controls how many alcoves are shown
                        Keep at 3-5 for dimension-3 groups (A3, B3, C3)

  Quick examples        15 pre-filled examples, click any to auto-fill
  Compute               Runs in background; window stays responsive
  Clear output          Clears the plot and console

  Dimension-2 results        2D image shown directly in the window
  Dimension-3 results        3D interactive viewer opens in your browser
                        (rotate, zoom, pan with the mouse)
                        Requires internet connection to load

Console panel (right side)
  After each run, shows:
    * Input element info:  s_XXX  =  t_(n1, n2, n3) - s_YYY
    * Coroot Basis vectors (columns = simple coroots)
    * A list of every element in the conjugacy class / coconjugation
      set that falls within the bounding box, each shown in BOTH
      formats:
            s_XXX     =  t_(n1, n2, n3) - s_YYY
      where t_(...) is the translation in basis B coordinates and
      s_YYY is the finite Weyl group part.


ELEMENT NOTATION
────────────────
Elements can be entered in EITHER of two formats:

  FORMAT 1 - Reduced word (string of generator indices):
    ""          ->  identity  e
    "e"         ->  identity  e  (same as blank)
    "1"         ->  s_1
    "012"       ->  s_0 s_1 s_2
    "01201210"  ->  s_0 s_1 s_2 s_0 s_1 s_2 s_1 s_0
    (non-reduced words are automatically reduced)

  FORMAT 2 - Translation * finite part:
    Dimension 2:       t_(n1,n2)*s_XX     or  t_(n1,n2)*e
    Dimension 3:       t_(n1,n2,n3)*s_XX  or  t_(n1,n2,n3)*e
    Example:      t_(1,2,1)*s_2312   for A3-tilde
    (n1, n2, n3) are coordinates in the coroot basis;
    s_XX is the finite Weyl part; "e" means identity finite part.

  A1 x A1 (special case - two separate fields):
    Field 1: first A1-tilde factor (accepts both formats, including
             the combined form t_(n1,n2)*(s_X,s_Y) for both factors)
    Field 2: second A1-tilde factor (digit form only, ignored if
             Field 1 uses the combined format)

  Pre-filled examples are available under "Quick examples".


VISUAL FEATURES
───────────────
  * Chosen element (conjugacy class):  bold red outline & label
  * Element 1 (coconjugation set):     red outline & label
  * Element 2 (coconjugation set):     blue outline & label
  * Identity alcove:                   cyan/light fill
  * Striped alcove:                    indicates the identity is
                                       currently in the set being
                                       computed (conjugacy class or
                                       coconjugation set)
  * Wireframe grid:                    all alcoves within the
                                       bounding box shown in gray
  * Coroot axes (dimension 2):              labeled with alpha^v_1 and
                                       alpha^v_2 (simple coroots)


TROUBLESHOOTING
───────────────
"Java not found"
  Install Java from https://adoptium.net and try again.

macOS: "cannot be verified" or "cannot be opened"
  Right-click the .command file -> Open -> Open  (once only).
  Then run:  chmod +x "Launch (macOS).command"

Anaconda/conda interference
  Run:  conda deactivate  before launching.
  The macOS launcher does this automatically.

Computation is slow
  Lower the bounding box. Dimension-3 groups with bbox > 5
  can take several minutes.

3D plot does not appear
  Click "Open 3D in browser" in the toolbar.
  Make sure you have an internet connection (three.js loads from CDN).

Emoji compile error when running Java directly
  Run this once in Terminal:
  python3 -c "
  src = open('AffineCoxeterExplorer.java').read()
  src = src.replace('\\\\U0001F537','').replace('\\\\U0001F310','')
  open('AffineCoxeterExplorer.java','w').write(src)
  print('done')
  "
