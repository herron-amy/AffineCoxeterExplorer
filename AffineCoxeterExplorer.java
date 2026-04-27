// AffineCoxeterExplorer.java
// Requires Java 11+  (run with:  java AffineCoxeterExplorer.java)
// SageMath must be installed:  https://www.sagemath.org

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.imageio.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.*;
import java.net.*;

public class AffineCoxeterExplorer extends JFrame {

    // ── Colours ─────────────────────────────────────────────────────────
    static final Color BG      = new Color(0x0d, 0x11, 0x17);
    static final Color SURF    = new Color(0x16, 0x1b, 0x22);
    static final Color SURF2   = new Color(0x1c, 0x23, 0x33);
    static final Color BORDER  = new Color(0x30, 0x36, 0x3d);
    static final Color ACCENT  = new Color(0x38, 0x8b, 0xfd);
    static final Color ACCENT2 = new Color(0x58, 0xa6, 0xff);
    static final Color TEXT    = new Color(0xe6, 0xed, 0xf3);
    static final Color MUTED   = new Color(0x7d, 0x85, 0x90);
    static final Color OK      = new Color(0x3f, 0xb9, 0x50);
    static final Color WARN    = new Color(0xd2, 0x99, 0x22);
    static final Color ERR     = new Color(0xf8, 0x51, 0x49);
    static final Font  MONO    = new Font(Font.MONOSPACED, Font.PLAIN, 11);

    // ── Group / computation data ────────────────────────────────────────
    static final String[][] GROUPS = {
        {"~A\u2081 \u00d7 ~A\u2081  (dimension 2, product)", "A1xA1"},
        {"~A\u2082  \u2014  Affine A\u2082", "aff2_A"},
        {"~B\u2082  \u2014  Affine B\u2082", "aff2_B"},
        {"~C\u2082  \u2014  Affine C\u2082", "aff2_C"},
        {"~G\u2082  \u2014  Affine G\u2082", "aff2_G"},
        {"~A\u2083  \u2014  Affine A\u2083", "A3"},
        {"~B\u2083  \u2014  Affine B\u2083", "B3"},
        {"~C\u2083  \u2014  Affine C\u2083", "C3"},
    };
    static final int[] DEF_BBOX = {10, 7, 7, 7, 7, 3, 3, 3};

    // ── Examples: {label, group, comp, e1, e2, e3, e4, bbox} ───────────
    static final Object[][] EXAMPLES = {
        {"~A\u2081\u00d7~A\u2081  conj (01,010)", "A1xA1","conj","01","010","","",10},
        {"~A\u2081\u00d7~A\u2081  cocon",          "A1xA1","cocon","01","01","10","01",10},
        {"~A\u2082  conj  s_010",                         "aff2_A","conj","010","","","",10},
        {"~A\u2082  cocon",                               "aff2_A","cocon","01021","0210210201020","","",5},
        {"~C\u2082  conj  s_12012101",                   "aff2_C","conj","12012101","","","",7},
        {"~C\u2082  cocon",                              "aff2_C","cocon","12012101","212101210121","","",7},
        {"~B\u2082  conj  s_0",                          "aff2_B","conj","0","","","",7},
        {"~G\u2082  conj  s_02121",                      "aff2_G","conj","02121","","","",7},
        {"~G\u2082  cocon",                              "aff2_G","cocon","02121","212102121","","",7},
        {"~A\u2083  conj  s_13",                          "A3","conj","13","","","",3},
        {"~A\u2083  cocon  s_13/s_1210\u2026",           "A3","cocon","13","1210301201","","",5},
        {"~B\u2083  conj  s_1",                          "B3","conj","1","","","",3},
        {"~B\u2083  cocon  s_1 / s_1",                   "B3","cocon","1","1","","",3},
        {"~C\u2083  conj  s_13",                         "C3","conj","13","","","",3},
        {"~C\u2083  cocon  s_13 / s_13",                 "C3","cocon","13","13","","",3},
    };

    // ── UI fields ───────────────────────────────────────────────────────
    JComboBox<String> groupCombo, compCombo;
    JPanel  inputsPanel;
    JTextField[] fields = new JTextField[4];
    JSlider bboxSlider;
    JLabel  bboxLabel, statusLabel, toolbarLabel;
    JButton computeBtn, clearBtn, open3dBtn;
    ImagePanel imagePanel;
    JTextArea  consoleArea;
    JPanel     consoleWrapper;

    String sagePath  = "sage";
    File   helperSage;
    File   lastHtmlFile;
    ExecutorService pool = Executors.newSingleThreadExecutor(r -> {
        Thread t = new Thread(r, "sage-worker");
        t.setDaemon(true);
        return t;
    });
    Future<?> currentJob;

    // ════════════════════════════════════════════════════════════════════
    // CONSTRUCTOR
    // ════════════════════════════════════════════════════════════════════
    public AffineCoxeterExplorer() {
        super("  Affine Coxeter Explorer");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setSize(1200, 760);
        setMinimumSize(new Dimension(900, 580));
        applyDarkTheme();
        findSage();
        locateHelper();
        buildUI();
        setLocationRelativeTo(null);
    }

    // ── Dark theme ──────────────────────────────────────────────────────
    void applyDarkTheme() {
        UIManager.put("Panel.background",                BG);
        UIManager.put("ScrollPane.background",           BG);
        UIManager.put("Viewport.background",             BG);
        UIManager.put("SplitPane.background",            BORDER);
        UIManager.put("ComboBox.background",             new Color(0x0d,0x11,0x17));
        UIManager.put("ComboBox.foreground",             TEXT);
        UIManager.put("ComboBox.selectionBackground",    ACCENT);
        UIManager.put("ComboBox.selectionForeground",    Color.WHITE);
        UIManager.put("ComboBox.buttonBackground",       SURF2);
        UIManager.put("TextField.background",            new Color(0x0d,0x11,0x17));
        UIManager.put("TextField.foreground",            TEXT);
        UIManager.put("TextField.caretForeground",       ACCENT2);
        UIManager.put("TextField.border",                BorderFactory.createLineBorder(BORDER));
        UIManager.put("TextField.selectionBackground",   ACCENT);
        UIManager.put("Slider.background",               SURF2);
        UIManager.put("Slider.foreground",               ACCENT);
        UIManager.put("Slider.thumbColor",               ACCENT);
        UIManager.put("ScrollBar.background",            SURF);
        UIManager.put("ScrollBar.thumb",                 BORDER);
        UIManager.put("ScrollBar.thumbDarkShadow",       BORDER);
        UIManager.put("ScrollBar.thumbHighlight",        SURF2);
        UIManager.put("ScrollBar.thumbShadow",           SURF);
        UIManager.put("TextArea.background",             new Color(0x0d,0x11,0x17));
        UIManager.put("TextArea.foreground",             new Color(0x8b,0x94,0x9e));
        UIManager.put("TextArea.caretForeground",        ACCENT2);
        UIManager.put("Label.foreground",                TEXT);
        UIManager.put("Button.background",               SURF2);
        UIManager.put("Button.foreground",               TEXT);
        UIManager.put("OptionPane.background",           SURF);
        UIManager.put("OptionPane.messageForeground",    TEXT);
    }

    // ── Find SageMath ────────────────────────────────────────────────────
    void findSage() {
        // Try PATH first
        if (canRun("sage")) { sagePath = "sage"; return; }
        // macOS: versioned app bundles
        File apps = new File("/Applications");
        if (apps.isDirectory()) {
            File[] ls = apps.listFiles();
            if (ls != null) for (File f : ls) {
                if (f.getName().startsWith("SageMath")) {
                    File s = new File(f, "sage");
                    if (s.canExecute()) { sagePath = s.getAbsolutePath(); return; }
                }
            }
        }
        // macOS Homebrew / common Linux
        for (String p : new String[]{"/usr/local/bin/sage","/opt/homebrew/bin/sage",
                                     "/opt/sage/sage", System.getProperty("user.home")+"/sage/sage"}) {
            if (new File(p).canExecute()) { sagePath = p; return; }
        }
        // Windows: Program Files
        for (int v = 12; v >= 9; v--) {
            String p = "C:\\Program Files\\SageMath "+v+"\\runtime\\bin\\sage.bat";
            if (new File(p).exists()) { sagePath = p; return; }
        }
    }

    boolean canRun(String cmd) {
        try {
            Process p = Runtime.getRuntime().exec(new String[]{cmd, "--version"});
            return p.waitFor() == 0;
        } catch (Exception e) { return false; }
    }

    // ── Locate compute_helper.sage ──────────────────────────────────────
    void locateHelper() {
        // Look in the same directory as this class file
        String here = AffineCoxeterExplorer.class.getProtectionDomain()
                          .getCodeSource().getLocation().getPath();
        File dir = new File(here).getParentFile();
        helperSage = new File(dir, "compute_helper.sage");
        if (!helperSage.exists()) {
            // Fall back to current working directory
            helperSage = new File("compute_helper.sage");
        }
        if (!helperSage.exists()) {
            JOptionPane.showMessageDialog(this,
                "compute_helper.sage not found.\nPlease keep it in the same folder as AffineCoxeterExplorer.java.",
                "File missing", JOptionPane.ERROR_MESSAGE);
        }
    }

    // ════════════════════════════════════════════════════════════════════
    // BUILD UI
    // ════════════════════════════════════════════════════════════════════
    void buildUI() {
        getContentPane().setBackground(BG);

        // Header
        JPanel header = new JPanel(new FlowLayout(FlowLayout.LEFT, 14, 10));
        header.setBackground(new Color(0x1c, 0x23, 0x33));
        header.setBorder(BorderFactory.createMatteBorder(0,0,1,0, BORDER));
        JLabel logo  = label("", 22, Font.PLAIN, TEXT);
        JLabel title = label("Affine Coxeter Explorer", 14, Font.BOLD, TEXT);
        JLabel sub   = label("  Conjugacy classes & coconjugation sets  \u2014  dimensions 2 and 3", 10, Font.PLAIN, MUTED);
        header.add(logo); header.add(title); header.add(sub);
        add(header, BorderLayout.NORTH);

        // Split pane: sidebar | main
        JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                                          buildSidebar(), buildMain());
        split.setDividerLocation(340);
        split.setDividerSize(4);
        split.setBorder(null);
        split.setBackground(BORDER);
        add(split, BorderLayout.CENTER);

        // Status bar
        JPanel statusBar = new JPanel(new FlowLayout(FlowLayout.LEFT, 10, 3));
        statusBar.setBackground(SURF);
        statusBar.setBorder(BorderFactory.createMatteBorder(1,0,0,0, BORDER));
        statusLabel = label("Ready.", 10, Font.PLAIN, MUTED);
        statusBar.add(statusLabel);
        add(statusBar, BorderLayout.SOUTH);
    }

    // ── Sidebar ──────────────────────────────────────────────────────────
    JScrollPane buildSidebar() {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setBackground(SURF);
        panel.setBorder(BorderFactory.createEmptyBorder(10,10,20,10));

        // 1. Group
        panel.add(card("\u2460  Select group", buildGroupCard()));
        panel.add(vgap(8));

        // 2. Computation
        panel.add(card("\u2461  Computation", buildCompCard()));
        panel.add(vgap(8));

        // 3. Elements
        inputsPanel = new JPanel();
        inputsPanel.setLayout(new BoxLayout(inputsPanel, BoxLayout.Y_AXIS));
        inputsPanel.setBackground(SURF2);
        panel.add(card("\u2462  Elements", inputsPanel));
        panel.add(vgap(8));

        // 4. Bounding box
        panel.add(card("\u2463  Bounding box", buildBboxCard()));
        panel.add(vgap(8));

        // Examples
        panel.add(card("\u26a1  Quick examples", buildExamplesCard()));
        panel.add(vgap(8));

        // Compute / Clear
        panel.add(card("  ", buildButtonsCard()));

        JScrollPane sp = new JScrollPane(panel,
            JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        sp.setBorder(null);
        sp.getViewport().setBackground(SURF);
        rebuildInputs();
        return sp;
    }

    JPanel buildGroupCard() {
        String[] labels = new String[GROUPS.length];
        for (int i=0;i<GROUPS.length;i++) labels[i] = GROUPS[i][0];
        groupCombo = new JComboBox<>(labels);
        styleCombo(groupCombo);
        groupCombo.addActionListener(e -> {
            int i = groupCombo.getSelectedIndex();
            bboxSlider.setValue(DEF_BBOX[i]);
            bboxLabel.setText(String.valueOf(DEF_BBOX[i]));
            rebuildInputs();
        });
        JPanel p = darkPanel();
        p.add(groupCombo);
        return p;
    }

    JPanel buildCompCard() {
        compCombo = new JComboBox<>(new String[]{"Conjugacy Class", "Coconjugation Set"});
        styleCombo(compCombo);
        compCombo.addActionListener(e -> rebuildInputs());
        JPanel p = darkPanel();
        p.add(compCombo);
        return p;
    }

    JPanel buildBboxCard() {
        bboxSlider = new JSlider(1, 15, 5);
        bboxSlider.setBackground(SURF2);
        bboxSlider.setForeground(ACCENT);
        bboxLabel = label("5", 12, Font.BOLD, ACCENT2);
        bboxLabel.setPreferredSize(new Dimension(28, 20));
        bboxSlider.addChangeListener(e ->
            bboxLabel.setText(String.valueOf(bboxSlider.getValue())));
        JPanel row = darkPanel();
        row.setLayout(new BorderLayout(6,0));
        row.add(bboxSlider, BorderLayout.CENTER);
        row.add(bboxLabel,  BorderLayout.EAST);
        JPanel p = darkPanel();
        p.add(row);
        p.add(label("Higher = more alcoves but slower.  Keep \u2264 5 for dimension 3.",
                     9, Font.PLAIN, MUTED));
        return p;
    }

    JPanel buildExamplesCard() {
        JPanel grid = new JPanel(new GridLayout(0, 2, 4, 4));
        grid.setBackground(SURF2);
        for (Object[] ex : EXAMPLES) {
            JButton b = new JButton((String)ex[0]);
            b.setBackground(SURF);
            b.setForeground(MUTED);
            b.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
            b.setBorder(BorderFactory.createLineBorder(BORDER));
            b.setFocusPainted(false);
            b.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            b.addMouseListener(new MouseAdapter() {
                public void mouseEntered(MouseEvent e) {
                    b.setBackground(SURF2); b.setForeground(ACCENT2);
                    b.setBorder(BorderFactory.createLineBorder(ACCENT));
                }
                public void mouseExited(MouseEvent e)  {
                    b.setBackground(SURF);  b.setForeground(MUTED);
                    b.setBorder(BorderFactory.createLineBorder(BORDER));
                }
            });
            final Object[] fex = ex;
            b.addActionListener(e -> loadExample(fex));
            grid.add(b);
        }
        JPanel p = darkPanel();
        p.add(grid);
        return p;
    }

    JPanel buildButtonsCard() {
        computeBtn = new JButton("\u25b6   Compute");
        computeBtn.setBackground(ACCENT);
        computeBtn.setForeground(Color.WHITE);
        computeBtn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 13));
        computeBtn.setFocusPainted(false);
        computeBtn.setBorderPainted(false);
        computeBtn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        computeBtn.setPreferredSize(new Dimension(Integer.MAX_VALUE, 38));
        computeBtn.setMaximumSize(new Dimension(Integer.MAX_VALUE, 38));
        computeBtn.addActionListener(e -> doCompute());

        clearBtn = new JButton("\u2715  Clear output");
        clearBtn.setBackground(SURF2);
        clearBtn.setForeground(MUTED);
        clearBtn.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        clearBtn.setFocusPainted(false);
        clearBtn.setBorderPainted(true);
        clearBtn.setBorder(BorderFactory.createLineBorder(BORDER));
        clearBtn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        clearBtn.setPreferredSize(new Dimension(Integer.MAX_VALUE, 30));
        clearBtn.setMaximumSize(new Dimension(Integer.MAX_VALUE, 30));
        clearBtn.addActionListener(e -> doClear());

        JPanel p = darkPanel();
        p.add(computeBtn);
        p.add(vgap(6));
        p.add(clearBtn);
        return p;
    }

    // ── Main area ────────────────────────────────────────────────────────
    JPanel buildMain() {
        JPanel main = new JPanel(new BorderLayout());
        main.setBackground(BG);

        // Toolbar
        JPanel toolbar = new JPanel(new BorderLayout(8,0));
        toolbar.setBackground(SURF);
        toolbar.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createMatteBorder(0,0,1,0,BORDER),
            BorderFactory.createEmptyBorder(5,14,5,10)));
        toolbarLabel = label("Select a group and element, then press Compute.", 10, Font.PLAIN, MUTED);
        open3dBtn = new JButton("Open 3D in browser");
        open3dBtn.setBackground(ACCENT);
        open3dBtn.setForeground(Color.WHITE);
        open3dBtn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
        open3dBtn.setFocusPainted(false);
        open3dBtn.setBorderPainted(false);
        open3dBtn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        open3dBtn.setVisible(false);
        open3dBtn.addActionListener(e -> open3dInBrowser());
        toolbar.add(toolbarLabel, BorderLayout.CENTER);
        toolbar.add(open3dBtn,    BorderLayout.EAST);
        main.add(toolbar, BorderLayout.NORTH);

        // Image panel (top)
        imagePanel = new ImagePanel();

        // Console panel (bottom) \u2014 always visible, resizable
        consoleArea = new JTextArea();
        consoleArea.setFont(MONO);
        consoleArea.setEditable(false);
        consoleArea.setBackground(new Color(0x0d, 0x11, 0x17));
        consoleArea.setForeground(new Color(0x8b, 0x94, 0x9e));
        consoleArea.setLineWrap(true);
        consoleArea.setWrapStyleWord(true);
        consoleArea.setText("SageMath output will appear here after each computation.");

        JScrollPane csp = new JScrollPane(consoleArea,
            JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        csp.setBorder(null);
        csp.getViewport().setBackground(new Color(0x0d,0x11,0x17));
        consoleArea.setLineWrap(false);  // disable wrap so horizontal scroll works

        // Console header bar
        JPanel consoleHeader = new JPanel(new BorderLayout(8, 0));
        consoleHeader.setBackground(new Color(0x16, 0x1b, 0x22));
        consoleHeader.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createMatteBorder(1,0,1,0,BORDER),
            BorderFactory.createEmptyBorder(4,12,4,8)));
        JLabel consoleTitle = label("SAGE OUTPUT", 8, Font.BOLD, MUTED);
        JButton clearConsoleBtn = new JButton("Clear");
        clearConsoleBtn.setBackground(new Color(0x16,0x1b,0x22));
        clearConsoleBtn.setForeground(MUTED);
        clearConsoleBtn.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
        clearConsoleBtn.setFocusPainted(false);
        clearConsoleBtn.setBorder(BorderFactory.createLineBorder(BORDER));
        clearConsoleBtn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        clearConsoleBtn.addActionListener(e -> {
            consoleArea.setText("");
        });
        consoleHeader.add(consoleTitle,    BorderLayout.WEST);
        consoleHeader.add(clearConsoleBtn, BorderLayout.EAST);

        consoleWrapper = new JPanel(new BorderLayout());
        consoleWrapper.setBackground(new Color(0x0d,0x11,0x17));
        consoleWrapper.add(consoleHeader, BorderLayout.NORTH);
        consoleWrapper.add(csp,           BorderLayout.CENTER);

        // Horizontal split: plot on left, console on right
        JSplitPane vertSplit = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                                              imagePanel, consoleWrapper);
        vertSplit.setDividerLocation(680);
        vertSplit.setDividerSize(5);
        vertSplit.setBorder(null);
        vertSplit.setBackground(BORDER);
        vertSplit.setResizeWeight(0.65);   // plot gets 65% of extra space
        main.add(vertSplit, BorderLayout.CENTER);

        return main;
    }

    // ── Dynamic inputs ────────────────────────────────────────────────────
    void rebuildInputs() {
        inputsPanel.removeAll();
        Arrays.fill(fields, null);
        String gkey = groupKey();
        String gkeyVal = gkey;
        boolean isCocon = compCombo.getSelectedIndex() == 1;

        if ("A1xA1".equals(gkey)) {
            if (!isCocon) {
                inputsPanel.add(label("1st ~A\u2081 component", 9, Font.PLAIN, MUTED));
                fields[0] = addField(inputsPanel, "e.g.  01  or  t_(2,0)*(s_1,e)");
                inputsPanel.add(label("2nd ~A\u2081 component  (digit form only)", 9, Font.PLAIN, MUTED));
                fields[1] = addField(inputsPanel, "e.g.  01  or  t_(2,0)*(s_1,e)");
            } else {
                inputsPanel.add(label("Element 1  \u2014  1st ~A\u2081", 9, Font.PLAIN, MUTED));
                fields[0] = addField(inputsPanel, "e.g.  01  or  t_(2,0)*(s_1,e)");
                inputsPanel.add(label("Element 1  \u2014  2nd ~A\u2081  (digit form only)", 9, Font.PLAIN, MUTED));
                fields[1] = addField(inputsPanel, "e.g.  01  or  t_(2,0)*(s_1,e)");
                inputsPanel.add(vgap(6));
                inputsPanel.add(label("Element 2  \u2014  1st ~A\u2081", 9, Font.PLAIN, MUTED));
                fields[2] = addField(inputsPanel, "e.g.  10  or  t_(0,1)*(e,s_1)");
                inputsPanel.add(label("Element 2  \u2014  2nd ~A\u2081  (digit form only)", 9, Font.PLAIN, MUTED));
                fields[3] = addField(inputsPanel, "e.g.  01  or  t_(2,0)*(s_1,e)");
            }
        } else {
            if (!isCocon) {
                inputsPanel.add(label("Element", 9, Font.PLAIN, MUTED));
                fields[0] = addField(inputsPanel, gkeyVal.startsWith("aff2_") ? "e.g.  12  or  t_(n1,n2)*s_xx,  identity means s_xx = e" : "e.g.  13  or  t_(n1,n2,n3)*s_xx,  identity means s_xx = e");
            } else {
                String ph1 = gkeyVal.startsWith("aff2_") ? "e.g.  12  or  t_(n1,n2)*s_xx,  identity means s_xx = e" : "e.g.  13  or  t_(n1,n2,n3)*s_xx,  identity means s_xx = e";
                String ph2 = gkeyVal.startsWith("aff2_") ? "e.g.  12  or  t_(n1,n2)*s_xx,  identity means s_xx = e" : "e.g.  13  or  t_(n1,n2,n3)*s_xx,  identity means s_xx = e";
                inputsPanel.add(label("Element 1", 9, Font.PLAIN, MUTED));
                fields[0] = addField(inputsPanel, ph1);
                inputsPanel.add(label("Element 2", 9, Font.PLAIN, MUTED));
                fields[1] = addField(inputsPanel, ph2);
            }
        }



        inputsPanel.revalidate();
        inputsPanel.repaint();
    }

    JTextField addField(JPanel parent, String hint) {
        JTextField tf = new JTextField();
        tf.setBackground(new Color(0x0d,0x11,0x17));
        tf.setForeground(TEXT);
        tf.setCaretColor(ACCENT2);
        tf.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createLineBorder(BORDER),
            BorderFactory.createEmptyBorder(4,7,4,7)));
        tf.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        tf.putClientProperty("hint", hint);
        tf.setForeground(MUTED);
        tf.setText(hint);
        tf.addFocusListener(new FocusAdapter() {
            public void focusGained(FocusEvent e) {
                if (tf.getText().equals((String)tf.getClientProperty("hint"))) {
                    tf.setText(""); tf.setForeground(TEXT);
                }
            }
            public void focusLost(FocusEvent e) {
                if (tf.getText().isEmpty()) {
                    tf.setText((String)tf.getClientProperty("hint"));
                    tf.setForeground(MUTED);
                }
            }
        });
        tf.addActionListener(e -> doCompute());   // Enter key
        tf.setMaximumSize(new Dimension(Integer.MAX_VALUE, 32));
        parent.add(tf);
        return tf;
    }

    // ── Load example ─────────────────────────────────────────────────────
    void loadExample(Object[] ex) {
        // ex: {label, group, comp, e1, e2, e3, e4, bbox}
        String gkey  = (String) ex[1];
        String ckey  = (String) ex[2];
        int    bbox   = (Integer) ex[7];

        // Set group
        for (int i=0;i<GROUPS.length;i++) {
            if (GROUPS[i][1].equals(gkey)) { groupCombo.setSelectedIndex(i); break; }
        }
        // Set comp
        compCombo.setSelectedIndex("conj".equals(ckey) ? 0 : 1);

        bboxSlider.setValue(bbox);
        bboxLabel.setText(String.valueOf(bbox));
        rebuildInputs();

        // Fill fields
        String[] vals = {(String)ex[3],(String)ex[4],(String)ex[5],(String)ex[6]};
        for (int i=0; i<4; i++) {
            if (fields[i] != null) {
                String v = vals[i];
                if (v != null && !v.isEmpty()) {
                    fields[i].setText(v);
                    fields[i].setForeground(TEXT);
                } else {
                    String hint = (String)fields[i].getClientProperty("hint");
                    fields[i].setText(hint);
                    fields[i].setForeground(MUTED);
                }
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════
    // COMPUTE
    // ════════════════════════════════════════════════════════════════════
    void doCompute() {
        if (currentJob != null && !currentJob.isDone()) {
            currentJob.cancel(true);
        }

        String[] inputs = getInputs();
        for (String v : inputs) {
            if (!v.isEmpty() && !v.equals("e") && !v.matches("[0-9]+") && !v.matches("t_\\(.*\\)\\*s_[0-9]*") && !v.matches("t_\\(.*\\)\\*e") && !v.matches("t_\\(.*\\)\\(.*\\)") && !v.matches("t_\\(.*\\)\\*\\(.*\\)")) {
                setStatus("\u2717  Input error: \""+v+"\" \u2014 use digits (e.g. 013) or t_(n1,n2)*s_XX format.", ERR);
                return;
            }
        }

        String gkey  = groupKey();
        String ckey  = "conj".equals(compKey()) ? "conj" : "cocon";
        int    bbox   = bboxSlider.getValue();

        computeBtn.setEnabled(false);
        setStatus("\u23f3  Computing\u2026 (SageMath is working, please wait)", WARN);
        toolbarLabel.setText("Computing\u2026");
        open3dBtn.setVisible(false);
        consoleArea.setText("SageMath output will appear here after each computation.");

        currentJob = pool.submit(() -> {
            try {
                File outDir = Files.createTempDirectory("weyl_out_").toFile();
                outDir.deleteOnExit();

                // Write params file
                File params = new File(outDir, "params.txt");
                try (PrintWriter pw = new PrintWriter(params)) {
                    pw.println("group="   + gkey);
                    pw.println("comp="    + ckey);
                    for (int i=0;i<4;i++) pw.println("input"+i+"="+inputs[i]);
                    pw.println("bbox="    + bbox);
                    pw.println("outdir="  + outDir.getAbsolutePath());
                }

                // Run sage
                List<String> cmd = new ArrayList<>(Arrays.asList(
                    sagePath, helperSage.getAbsolutePath(), params.getAbsolutePath()));
                ProcessBuilder pb = new ProcessBuilder(cmd);
                pb.redirectErrorStream(true);
                Process proc = pb.start();

                StringBuilder out = new StringBuilder();
                try (BufferedReader br = new BufferedReader(
                        new InputStreamReader(proc.getInputStream()))) {
                    String ln;
                    while ((ln=br.readLine()) != null) out.append(ln).append("\n");
                }
                proc.waitFor();

                // Read result.json
                File rjson = new File(outDir, "result.json");
                if (!rjson.exists()) {
                    final String errOut = out.toString();
                    SwingUtilities.invokeLater(() -> {
                        setStatus("\u2717  Sage error (no result file)", ERR);
                        toolbarLabel.setText("Error.");
                        showConsole(errOut);
                        computeBtn.setEnabled(true);
                    });
                    return;
                }

                String json = new String(Files.readAllBytes(rjson.toPath()));
                String rtype   = jsonStr(json, "type");
                String rpath   = jsonStr(json, "path");
                String error   = jsonStr(json, "error");
                // Read console from separate file to avoid JSON quote issues
                File consolef = new File(outDir, "console.txt");
                String console = consolef.exists()
                    ? new String(Files.readAllBytes(consolef.toPath()))
                    : out.toString();

                if (!error.isEmpty()) {
                    final String ferr = error;
                    final String fcon = console;
                    SwingUtilities.invokeLater(() -> {
                        setStatus("\u2717  " + ferr.lines().reduce("",(a,b)->b.isEmpty()?a:b), ERR);
                        toolbarLabel.setText("Error.");
                        showConsole(ferr + "\n" + fcon);
                        computeBtn.setEnabled(true);
                    });
                    return;
                }

                final String ft=rtype, fp=rpath, fc=console;
                final String desc = buildDesc(gkey, ckey, inputs, bbox);
                SwingUtilities.invokeLater(() -> {
                    showResult(ft, fp, fc, desc);
                    computeBtn.setEnabled(true);
                });

            } catch (Exception ex) {
                SwingUtilities.invokeLater(() -> {
                    setStatus("\u2717  " + ex.getMessage(), ERR);
                    toolbarLabel.setText("Error.");
                    computeBtn.setEnabled(true);
                });
            }
        });
    }

    void showResult(String type, String path, String console, String desc) {
        showConsole(console);
        if ("png".equals(type)) {
            try {
                BufferedImage img = ImageIO.read(new File(path));
                imagePanel.setImage(img);
                open3dBtn.setVisible(false);
            } catch (Exception e) {
                setStatus("\u2717  Could not load image: " + e.getMessage(), ERR);
                return;
            }
        } else if ("html".equals(type)) {
            lastHtmlFile = new File(path);
            imagePanel.set3DPlaceholder();
            open3dBtn.setVisible(true);
            open3dInBrowser();
        }
        toolbarLabel.setText(desc + "  \u00b7  \u2713 Done");
        setStatus("\u2713  Done", OK);
    }

    void open3dInBrowser() {
        if (lastHtmlFile != null && lastHtmlFile.exists()) {
            try {
                Desktop.getDesktop().browse(lastHtmlFile.toURI());
            } catch (Exception e) {
                setStatus("Could not open browser: " + e.getMessage(), ERR);
            }
        }
    }

    void doClear() {
        imagePanel.clear();
        consoleArea.setText("SageMath output will appear here after each computation.");
        toolbarLabel.setText("Select a group and element, then press Compute.");
        open3dBtn.setVisible(false);
        setStatus("", MUTED);
    }

    // ════════════════════════════════════════════════════════════════════
    // HELPERS
    // ════════════════════════════════════════════════════════════════════
    String groupKey() {
        int i = groupCombo.getSelectedIndex();
        return (i>=0 && i<GROUPS.length) ? GROUPS[i][1] : "";
    }

    String compKey() {
        return compCombo.getSelectedIndex()==0 ? "conj" : "cocon";
    }

    String[] getInputs() {
        String[] r = new String[4];
        for (int i=0;i<4;i++) {
            if (fields[i]==null) { r[i]=""; continue; }
            String v = fields[i].getText().trim();
            String hint = (String)fields[i].getClientProperty("hint");
            r[i] = v.equals(hint) ? "" : v;
        }
        return r;
    }

    String buildDesc(String gkey, String ckey, String[] inputs, int bbox) {
        String gl = "";
        for (String[] g : GROUPS) if (g[1].equals(gkey)) { gl=g[0]; break; }
        String cl = "conj".equals(ckey) ? "Conjugacy Class" : "Coconjugation Set";
        StringBuilder sb = new StringBuilder();
        boolean first=true;
        for (String v : inputs) {
            if (v!=null && !v.isEmpty()) {
                if (!first) sb.append(" / ");
                sb.append("s_").append(v); first=false;
            }
        }
        if (first) sb.append("e");
        return cl + "  \u00b7  " + gl + "  \u00b7  " + sb + "  \u00b7  bbox=" + bbox;
    }

    void setStatus(String msg, Color col) {
        statusLabel.setText(msg);
        statusLabel.setForeground(col);
    }

    void showConsole(String text) {
        if (text == null || text.isBlank()) return;
        // Append with a divider so multiple runs are visible
        String current = consoleArea.getText();
        if (!current.isEmpty() &&
            !current.equals("SageMath output will appear here after each computation.")) {
            consoleArea.append("\n" + "─".repeat(60) + "\n");
        } else {
            consoleArea.setText("");
        }
        consoleArea.append(text.stripTrailing());
        consoleArea.setCaretPosition(consoleArea.getDocument().getLength());
    }

    /** Very small JSON string extractor \u2014 avoids external dependency. */
    String jsonStr(String json, String key) {
        String pat = "\""+key+"\"";
        int ki = json.indexOf(pat);
        if (ki<0) return "";
        int q1 = json.indexOf('"', ki+pat.length()+1);
        if (q1<0) return "";
        StringBuilder sb = new StringBuilder();
        int i = q1+1;
        while (i<json.length()) {
            char c = json.charAt(i);
            if (c=='"') break;
            if (c=='\\'  && i+1<json.length()) {
                char n = json.charAt(++i);
                switch(n) {
                    case 'n': sb.append('\n'); break;
                    case 't': sb.append('\t'); break;
                    case 'r': sb.append('\r'); break;
                    default:  sb.append(n);
                }
            } else sb.append(c);
            i++;
        }
        return sb.toString();
    }

    // ── Widget factories ──────────────────────────────────────────────────
    JLabel label(String text, int size, int style, Color col) {
        JLabel l = new JLabel(text);
        l.setFont(new Font(Font.SANS_SERIF, style, size));
        l.setForeground(col);
        l.setBackground(SURF2);
        l.setAlignmentX(Component.LEFT_ALIGNMENT);
        return l;
    }

    JPanel darkPanel() {
        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p, BoxLayout.Y_AXIS));
        p.setBackground(SURF2);
        p.setAlignmentX(Component.LEFT_ALIGNMENT);
        return p;
    }

    Component vgap(int h) {
        return Box.createVerticalStrut(h);
    }

    void styleCombo(JComboBox<String> cb) {
        cb.setBackground(new Color(0x0d,0x11,0x17));
        cb.setForeground(TEXT);
        cb.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        cb.setBorder(BorderFactory.createLineBorder(BORDER));
        cb.setMaximumSize(new Dimension(Integer.MAX_VALUE, 30));
        cb.setAlignmentX(Component.LEFT_ALIGNMENT);
    }

    JPanel card(String title, JComponent content) {
        JPanel outer = new JPanel(new BorderLayout());
        outer.setBackground(SURF2);
        outer.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createLineBorder(BORDER),
            BorderFactory.createEmptyBorder(10,12,10,12)));
        outer.setAlignmentX(Component.LEFT_ALIGNMENT);
        outer.setMaximumSize(new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE));
        JLabel hdr = label(title.toUpperCase(), 8, Font.BOLD, MUTED);
        hdr.setBorder(BorderFactory.createEmptyBorder(0,0,8,0));
        JPanel stripe = new JPanel();
        stripe.setBackground(ACCENT);
        stripe.setPreferredSize(new Dimension(3, 1));
        JPanel inner = new JPanel(new BorderLayout(8,0));
        inner.setBackground(SURF2);
        inner.add(stripe, BorderLayout.WEST);
        JPanel right = new JPanel();
        right.setLayout(new BoxLayout(right, BoxLayout.Y_AXIS));
        right.setBackground(SURF2);
        right.add(hdr);
        right.add(content);
        inner.add(right, BorderLayout.CENTER);
        outer.add(inner, BorderLayout.CENTER);
        return outer;
    }

    // ════════════════════════════════════════════════════════════════════
    // IMAGE PANEL
    // ════════════════════════════════════════════════════════════════════
    static class ImagePanel extends JPanel {
        BufferedImage img;
        String placeholderText;

        ImagePanel() {
            setBackground(new Color(0x0d, 0x11, 0x17));
            placeholderText =
                "\u2B21\n\nNo result yet\n\nChoose a group, enter an element,\nthen press Compute.";
        }

        void setImage(BufferedImage i) {
            img = i; placeholderText = null; repaint();
        }

        void set3DPlaceholder() {
            img = null;
            placeholderText = "\n\n3D interactive plot\n\nOpening in your browser\u2026\n\n"

                            + "(rotate \u00b7 zoom \u00b7 pan with the mouse)";
            repaint();
        }

        void clear() {
            img = null;
            placeholderText =
                "\u2B21\n\nNo result yet\n\nChoose a group, enter an element,\nthen press Compute.";
            repaint();
        }

        @Override protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            int w=getWidth(), h=getHeight();
            if (img != null) {
                double scale = Math.min((double)w/img.getWidth(), (double)h/img.getHeight());
                scale = Math.min(scale, 1.0);
                int dw=(int)(img.getWidth()*scale), dh=(int)(img.getHeight()*scale);
                int x=(w-dw)/2, y=(h-dh)/2;
                g.drawImage(img, x, y, dw, dh, null);
            } else if (placeholderText != null) {
                g.setColor(new Color(0x7d, 0x85, 0x90));
                g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 14));
                FontMetrics fm = g.getFontMetrics();
                String[] lines = placeholderText.split("\\n");
                int totalH = lines.length * (fm.getHeight()+4);
                int y = (h - totalH)/2 + fm.getAscent();
                for (String line : lines) {
                    int x = (w - fm.stringWidth(line))/2;
                    g.drawString(line, x, y);
                    y += fm.getHeight()+4;
                }
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════
    // MAIN
    // ════════════════════════════════════════════════════════════════════
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            try { UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName()); }
            catch (Exception ignored) {}
            new AffineCoxeterExplorer().setVisible(true);
        });
    }
}