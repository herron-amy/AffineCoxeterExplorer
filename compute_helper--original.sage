# compute_helper.sage
# Called by AffineWeylExplorer.java as:
#   sage compute_helper.sage <params_file>
#
# params_file is a key=value text file written by Java.
# Results are written back to the directory specified in params.

import sys, os, contextlib, io, tempfile, shutil, json as _json

# ── Read parameters ─────────────────────────────────────────
_params = {}
with open(sys.argv[1]) as _pf:
    for _line in _pf:
        if "=" in _line:
            _k, _v = _line.strip().split("=", 1)
            _params[_k] = _v

_group   = _params.get("group", "")
_comp    = _params.get("comp", "conj")
_inputs  = [_params.get(f"input{_i}", "") for _i in range(4)]
_bbox    = int(_params.get("bbox", "5"))
_out_dir = _params.get("outdir", "/tmp")

# ── Computation code ────────────────────────────────────────
# ============================================================
# IMPORTS  (run this cell first)
# ============================================================
import re


# ============================================================
# Ã₁ × Ã₁  —  helper functions
# ============================================================

def _a1xa1_add_stripes_to_square(G, corners, vec, zorder=6,
                                  num_stripes=5, stripe_color='darkgray', stripe_width=1.5):
    v0 = vector([corners[0][0]+vec[0], corners[0][1]+vec[1]])
    v1 = vector([corners[1][0]+vec[0], corners[1][1]+vec[1]])
    v2 = vector([corners[2][0]+vec[0], corners[2][1]+vec[1]])
    v3 = vector([corners[3][0]+vec[0], corners[3][1]+vec[1]])
    total_stripes = num_stripes * 2
    for i in range(1, total_stripes + 1):
        t = i / (total_stripes + 1)
        if t <= 0.5:
            p1 = v0 + (2*t)*(v3-v0)
            p2 = v0 + (2*t)*(v1-v0)
        else:
            t_adj = 2*(t-0.5)
            p1 = v3 + t_adj*(v2-v3)
            p2 = v1 + t_adj*(v2-v1)
        G += line2d([p1, p2], color=stripe_color, thickness=stripe_width, zorder=zorder)
    return G


def _a1xa1_get_element_for_interval(start, W, S):
    if start == 0:
        return W.one()
    element = W.one(); current = 0; target = start
    while current != target:
        if current < target:
            element = element * S[0 if current%2==0 else 1]; current += 1
        else:
            element = element * S[1 if current%2==0 else 0]; current -= 1
    return element


def _a1xa1_reduced_word_square(bounding_box, rw_string):
    RS = RootSystem(["A",1,1]); L = RS.root_lattice()
    W = L.weyl_group(); S = W.simple_reflections()
    if rw_string == '':
        return 0
    element = W.one()
    for char in rw_string:
        element = element * S[int(char)]
    for pos in range(-bounding_box, bounding_box):
        if _a1xa1_get_element_for_interval(pos, W, S) == element:
            return pos
    raise ValueError(f"Could not find position for {rw_string}")


def _a1xa1_plot_grid(bounding_box, H0_alcoves_dict):
    RS1 = RootSystem(["A",1,1]); RS2 = RootSystem(["A",1,1])
    L1 = RS1.root_lattice(); L2 = RS2.root_lattice()
    W1 = L1.weyl_group(); W2 = L2.weyl_group()
    S1 = W1.simple_reflections(); S2 = W2.simple_reflections()
    G = Graphics()
    for i in range(-bounding_box, bounding_box+1):
        G += line([(i,-bounding_box),(i,bounding_box)], color='dimgray', thickness=0.5, zorder=20)
        G += line([(-bounding_box,i),(bounding_box,i)], color='dimgray', thickness=0.5, zorder=20)
    label_fontsize = max(8, min(15, int(-bounding_box+20)))
    colors = [(0.7,1.0,0.7),(1.0,0.7,1.0),(1.0,1.0,0.7),(0.7,1.0,1.0)]
    tmp = 0
    for i in range(-1,1):
        for j in range(-1,1):
            corners = [(i,j),(i+1,j),(i+1,j+1),(i,j+1)]
            w1 = _a1xa1_get_element_for_interval(i, W1, S1)
            w2 = _a1xa1_get_element_for_interval(j, W2, S2)
            rw1 = w1.reduced_word(); rw2 = w2.reduced_word()
            color = colors[tmp]; tmp += 1
            cx = i+0.5; cy = j+0.5
            l1 = "e" if not rw1 else f"s_{{{''.join(str(x) for x in rw1)}}}"
            l2 = "e" if not rw2 else f"s_{{{''.join(str(x) for x in rw2)}}}"
            label = f"$({l1},{l2})$"
            G += polygon2d(corners, fill=True, color=color, alpha=1, edgecolor='dimgray', thickness=1, zorder=1)
            G += text(label, (cx,cy), fontsize=label_fontsize, color='black', zorder=2, alpha=1)
            for key_element in H0_alcoves_dict.keys():
                if (rw1,rw2) == (key_element[0].reduced_word(), key_element[1].reduced_word()):
                    H0_alcoves_dict[key_element] = {
                        'corners': corners, 'color': color,
                        'position': (i,j), 'center': (QQ(cx), QQ(cy)),
                        'reduced_words': (rw1,rw2)
                    }
    G += arrow((0,0),(2,0), color='red', width=2, zorder=2, arrowsize=3)
    G += arrow((0,0),(0,2), color='blue', width=2, zorder=2, arrowsize=3)
    G += text(r"First Ã₁ $\alpha^\vee$", (2.2,0.2), fontsize=10, color='red',
              horizontal_alignment='left', vertical_alignment='top')
    G += text(r"Second Ã₁ $\alpha^\vee$", (0,2.25), fontsize=10, color='blue', rotation=90,
              horizontal_alignment='right', vertical_alignment='bottom')
    return G, H0_alcoves_dict


def _a1xa1_chosen_elt(G, elt1, elt2, bounding_box, coroot_vertices, B, H0_alcoves_dict, border_color):
    RS1 = RootSystem(["A",1,1]); RS2 = RootSystem(["A",1,1])
    W1 = RS1.root_lattice().weyl_group(); W2 = RS2.root_lattice().weyl_group()
    S1 = W1.simple_reflections(); S2 = W2.simple_reflections()
    id1 = W1.from_reduced_word([]); id2 = W2.from_reduced_word([])
    # Check for t_(n1,n2)*(w1,w2) format in elt1
    _parsed_a1 = _a1xa1_parse_translation(elt1) if elt1 else None
    if _parsed_a1 is not None:
        _n1,_n2,_w1,_w2 = _parsed_a1
        # s_1 alcove sits at odd grid pos (2*n-1); identity sits at even pos (2*n)
        pos1 = 2*_n1 - (1 if _w1 else 0)
        pos2 = 2*_n2 - (1 if _w2 else 0)
        elt1 = _w1; elt2 = _w2
    ew1 = W1.from_reduced_word([Integer(c) for c in elt1]); rw1 = ew1.reduced_word()
    ew2 = W2.from_reduced_word([Integer(c) for c in elt2]); rw2 = ew2.reduced_word()
    if _parsed_a1 is None:
        pos1 = _a1xa1_reduced_word_square(bounding_box, rw1)
        pos2 = _a1xa1_reduced_word_square(bounding_box, rw2)
    cx = pos1+0.5; cy = pos2+0.5
    l1 = "e" if not rw1 else f"s_{{{''.join(str(x) for x in rw1)}}}"
    l2 = "e" if not rw2 else f"s_{{{''.join(str(x) for x in rw2)}}}"
    label = f"$({l1},{l2})$"
    label_fontsize = max(8, min(15, int(-bounding_box+20)))
    if vector([pos1,pos2]) in coroot_vertices:
        tp_euclid = vector([pos1,pos2]); cp_tmp = (id1,id2)
    elif vector([pos1+1,pos2]) in coroot_vertices:
        tp_euclid = vector([pos1+1,pos2]); cp_tmp = (S1[1],id2)
    elif vector([pos1,pos2+1]) in coroot_vertices:
        tp_euclid = vector([pos1,pos2+1]); cp_tmp = (id1,S2[1])
    else:
        tp_euclid = vector([pos1+1,pos2+1]); cp_tmp = (S1[1],S2[1])
    trans_part = B.inverse()*tp_euclid
    alcove_color = H0_alcoves_dict[cp_tmp]['color']
    corners = [(pos1,pos2),(pos1+1,pos2),(pos1+1,pos2+1),(pos1,pos2+1)]
    G += polygon2d(corners, fill=True, color=alcove_color, alpha=1,
                   edgecolor='dimgray', thickness=1, zorder=1)
    G += polygon2d(corners, fill=False, color=border_color, thickness=3, zorder=30)
    G += text(label, (cx,cy), fontsize=label_fontsize, color='black', zorder=25, alpha=1)
    return G, trans_part, cp_tmp, (cp_tmp[0].reduced_word(), cp_tmp[1].reduced_word())


# ─── Conjugacy class helpers ─────────────────────────────────────────────────

def _a1xa1_span_vectors(vertex1, vertex2, B, H0_alcoves_dict, coroot_vertices_in_B):
    v1B = B.inverse()*vertex1; v2B = B.inverse()*vertex2
    VS = ZZ**2
    # Filter spans against coroot_vertices_in_B to get finite lists (matches original notebook)
    s1_span   = VS.span([2*v1B])
    s2_span   = VS.span([2*v2B])
    s12_span  = VS.span([2*v1B, 2*v2B])
    ms = [
        [vector([0,0])],
        [v for v in coroot_vertices_in_B if v in s1_span],
        [v for v in coroot_vertices_in_B if v in s2_span],
        [v for v in coroot_vertices_in_B if v in s12_span],
    ]
    for i,k in enumerate(H0_alcoves_dict.keys()):
        H0_alcoves_dict[k]['mod_plot'] = ms[i]
    return H0_alcoves_dict


def _a1xa1_reflect_cc(trans_part, cp_tmp, H0_alcoves_dict):
    lpm = [trans_part + vector(H0_alcoves_dict[cp_tmp]['mod_plot'][i])
           for i in range(len(H0_alcoves_dict[cp_tmp]['mod_plot']))]
    results = [
        lpm,
        [vector([-v[0],v[1]]) for v in lpm],
        [vector([v[0],-v[1]]) for v in lpm],
        [vector([-v[0],-v[1]]) for v in lpm],
    ]
    for i,k in enumerate(H0_alcoves_dict.keys()):
        H0_alcoves_dict[k]['reflect'] = results[i]
    return H0_alcoves_dict


def _a1xa1_finite_cc(cp_tmp, H0_alcoves_dict):
    w = cp_tmp
    tmp = []
    for g in H0_alcoves_dict.keys():
        c1 = g[0]*w[0]*g[0].inverse(); c2 = g[1]*w[1]*g[1].inverse()
        tmp.append((c1.reduced_word(), c2.reduced_word()))
    for i,k in enumerate(H0_alcoves_dict.keys()):
        H0_alcoves_dict[k]['finite'] = tmp[i]
    return H0_alcoves_dict


def _a1xa1_conjugacy_class_plot(S1, S2, id1, id2, G, B, H0_alcoves_dict, bounding_box):
    for k in H0_alcoves_dict.keys():
        tmp = []
        for vec in H0_alcoves_dict[k]['reflect']:
            check = B*vec
            if abs(check[0]) <= bounding_box and abs(check[1]) <= bounding_box:
                tmp.append(check)
        H0_alcoves_dict[k]['t_lambda'] = tmp
    for k in H0_alcoves_dict.keys():
        fc = H0_alcoves_dict[k]['finite']
        if fc == ([],[]):           ak = (id1,id2)
        elif fc == ([1],[]):        ak = (S1[1],id2)
        elif fc == ([],[1]):        ak = (id1,S2[1])
        elif fc == ([1],[1]):       ak = (S1[1],S2[1])
        else: continue
        alcove = H0_alcoves_dict[ak]
        corners = alcove['corners']; color = alcove['color']
        for vec in H0_alcoves_dict[k]['t_lambda']:
            ok = all(abs(corners[i][j]+vec[j]) <= bounding_box
                     for i in range(4) for j in range(2))
            if ok:
                G += polygon2d([vector([corners[i][0]+vec[0],corners[i][1]+vec[1]]) for i in range(4)],
                               fill=True, color=color, alpha=0.4)
                if abs(float(vec[0]))<0.001 and abs(float(vec[1]))<0.001:
                    G = _a1xa1_add_stripes_to_square(G, corners, vec, zorder=2)
    G.show(frame=False, axes=False, figsize=12)


# ─── Coconjugation helpers ────────────────────────────────────────────────────

def _a1xa1_cocon_H0(cp1, cp2, H0_alcoves_dict):
    cocon = []
    for g in H0_alcoves_dict.keys():
        c1 = g[0]*cp1[0]*g[0].inverse(); c2 = g[1]*cp1[1]*g[1].inverse()
        if (c1.reduced_word() == cp2[0].reduced_word() and
                c2.reduced_word() == cp2[1].reduced_word()):
            cocon.append(g)
    return cocon


def _a1xa1_mod_set(vertex1, vertex2, B, H0_alcoves_dict, coroot_vertices_in_B, cp_tmp):
    v1B = B.inverse()*vertex1; v2B = B.inverse()*vertex2
    VS = ZZ**2
    spans = [
        [vector([0,0])],
        VS.span([2*v1B]),
        VS.span([2*v2B]),
        VS.span([2*v1B, 2*v2B]),
    ]
    for i,k in enumerate(H0_alcoves_dict.keys()):
        if k == cp_tmp:
            return spans[i]
    return [vector([0,0])]


def _a1xa1_cocon_H0_lambda(tp1, tp2, cocon_H0, mod_set_h0, H0_alcoves_dict):
    keys_list = list(H0_alcoves_dict.keys())
    reflect_list = [tp1, vector([-tp1[0],tp1[1]]), vector([tp1[0],-tp1[1]]), vector([-tp1[0],-tp1[1]])]
    cocon_H0_lambda = []; lpmul = []
    for i,g in enumerate(cocon_H0):
        idx = keys_list.index(g)
        ulambda = reflect_list[idx]
        diff = tp2 - ulambda
        if diff in mod_set_h0:
            cocon_H0_lambda.append(g); lpmul.append(diff)
    return lpmul, cocon_H0_lambda


def _a1xa1_eta(coroot_vertices_in_B, lpmul_i, cocon_H0_lambda_i, cp2):
    l0 = lpmul_i[0]; l1 = lpmul_i[1]
    cp = (cp2[0].reduced_word(), cp2[1].reduced_word())
    if cp == ([],[]):
        eta = coroot_vertices_in_B if (l0==0 and l1==0) else []
    elif cp == ([1],[]):
        eta = [v for v in coroot_vertices_in_B if 2*v[0]==l0 and 0==l1]
    elif cp == ([],[1]):
        eta = [v for v in coroot_vertices_in_B if 0==l0 and 2*v[1]==l1]
    else:
        eta = [v for v in coroot_vertices_in_B if 2*v[0]==l0 and 2*v[1]==l1]
    return eta, [[v, cocon_H0_lambda_i] for v in eta]


def _a1xa1_plot_cocon(G, bounding_box, B, S1, S2, id1, id2, H0_alcoves_dict, eta_cocon):
    plotted = set()
    filtered = []
    for elt in eta_cocon:
        vec = elt[0]; rw = elt[1]
        check = B*vec
        if abs(check[0]) <= bounding_box and abs(check[1]) <= bounding_box:
            filtered.append([check, rw])
    h0_overlaps = {(tuple(e[1][0].reduced_word()), tuple(e[1][1].reduced_word()))
                   for e in filtered if e[0][0]==0 and e[0][1]==0}
    for elt in filtered:
        vec = elt[0]; we = elt[1]; ak = we
        if ak not in H0_alcoves_dict: continue
        alcove = H0_alcoves_dict[ak]; corners = alcove['corners']; color = alcove['color']
        ok = all(abs(corners[i][j]+vec[j]) <= bounding_box for i in range(4) for j in range(2))
        if ok:
            G += polygon2d([vector([corners[i][0]+vec[0],corners[i][1]+vec[1]]) for i in range(4)],
                           fill=True, color=color, alpha=1, zorder=1)
            if vec[0]==0 and vec[1]==0:
                G = _a1xa1_add_stripes_to_square(G, corners, vec, zorder=2)
            rw_key = (tuple(ak[0].reduced_word()), tuple(ak[1].reduced_word()))
            plotted.add((float(vec[0]), float(vec[1]), rw_key))
    for key in H0_alcoves_dict:
        rw = (tuple(key[0].reduced_word()), tuple(key[1].reduced_word()))
        if (0.0, 0.0, rw) not in plotted:
            alcove = H0_alcoves_dict[key]; corners = alcove['corners']; color = alcove['color']
            G += polygon2d([vector([corners[i][0],corners[i][1]]) for i in range(4)],
                           fill=True, color=color, alpha=1, zorder=1)
    return G


# ============================================================
# Ã₁ × Ã₁  —  main computation functions
# ============================================================

def _a1xa1_parse_translation(s):
    """Parse t_(n1,n2)*(w1,w2) format for A1xA1. Returns (n1,n2,w1,w2) or None."""
    import re
    m = re.match(r't_\((-?\d+),(-?\d+)\)\*\((?:s_)?([01]*|e),(?:s_)?([01]*|e)\)$', s.strip())
    if m: return int(m.group(1)),int(m.group(2)),'' if m.group(3)=='e' else m.group(3),'' if m.group(4)=='e' else m.group(4)
    return None


def _Cox_Element_A1xA1(elt_1, elt_2, bounding_box):
    """Conjugacy class in Ã₁ × Ã₁."""
    vertex1 = vector([2,0]); vertex2 = vector([0,2])
    B = matrix([vertex1,vertex2]).transpose()
    coroot_vertices = [vector(i*vertex1+j*vertex2)
                       for i in range(-bounding_box,bounding_box)
                       for j in range(-bounding_box,bounding_box)
                       if abs(i*vertex1[0]+j*vertex2[0])<=bounding_box
                       and abs(i*vertex1[1]+j*vertex2[1])<=bounding_box]
    coroot_vertices_in_B = [vector([ZZ(round((B.inverse()*v)[k].n(),0)) for k in range(2)])
                            for v in coroot_vertices]
    RS1 = RootSystem(["A",1,1]); RS2 = RootSystem(["A",1,1])
    W1 = RS1.root_lattice().weyl_group(); W2 = RS2.root_lattice().weyl_group()
    S1 = W1.simple_reflections(); S2 = W2.simple_reflections()
    id1 = W1.from_reduced_word([]); id2 = W2.from_reduced_word([])
    H0 = {(id1,id2):None, (S1[1],id2):None, (id1,S2[1]):None, (S1[1],S2[1]):None}
    G, H0 = _a1xa1_plot_grid(bounding_box, H0)
    G += point2d(coroot_vertices, size=70, color='black')
    G, tp, cp_tmp, cp = _a1xa1_chosen_elt(G, elt_1, elt_2, bounding_box,
                                           coroot_vertices, B, H0, 'red')
    H0 = _a1xa1_span_vectors(vertex1, vertex2, B, H0, coroot_vertices_in_B)
    H0 = _a1xa1_reflect_cc(tp, cp_tmp, H0)
    H0 = _a1xa1_finite_cc(cp_tmp, H0)
    # Print conjugacy class elements as (s_XXX, s_YYY)
    _cc_a1xa1 = []; _cc_a1xa1_seen = set()
    for k in H0.keys():
        fc1, fc2 = H0[k]['finite']  # reduced word tuples for each A1 factor
        for vec in H0[k]['reflect']:
            check = B * vec
            if abs(check[0]) > bounding_box or abs(check[1]) > bounding_box:
                continue
            # vec[0] = coroot translation units in factor 1, vec[1] in factor 2
            try:
                n1 = ZZ(round(float(vec[0]), 0)); n2 = ZZ(round(float(vec[1]), 0))
            except Exception:
                continue  # skip infinite/invalid translations
            # Compute reduced word of each A1 factor: t_{n*alpha} * w
            def _a1_rw(n, fc):
                if n > 0:
                    elem = (S1[1]*S1[0])**n * W1.from_reduced_word(list(fc))
                elif n < 0:
                    elem = (S1[0]*S1[1])**(-n) * W1.from_reduced_word(list(fc))
                else:
                    elem = W1.from_reduced_word(list(fc))
                rw = elem.reduced_word()
                return 's_' + ''.join(str(i) for i in rw) if rw else 'e'
            e1_str = _a1_rw(n1, fc1)
            e2_str = _a1_rw(n2, fc2)
            fc1_str = 's_' + ''.join(str(x) for x in fc1) if fc1 else 'e'
            fc2_str = 's_' + ''.join(str(x) for x in fc2) if fc2 else 'e'
            _key = f'({e1_str}, {e2_str})'
            if _key not in _cc_a1xa1_seen:
                _cc_a1xa1_seen.add(_key)
                _cc_a1xa1.append(((n1, n2), _key, f'({fc1_str}, {fc2_str})'))
    def _fmt_v(v):
        def _fv(x):
            s = str(x)
            if any(c in s for c in ['sqrt', 'pi', '/', 'e^']):
                return s
            try:
                f = float(x)
                return f"{f:.10g}"
            except Exception:
                return s
        return [_fv(x) for x in v]
    print(f"Coroot Basis:")
    print(f"  b1 = {_fmt_v(vertex1)}  (1st Ã₁ simple coroot)")
    print(f"  b2 = {_fmt_v(vertex2)}  (2nd Ã₁ simple coroot)")
    print(f"Translation vectors are in basis B coordinates (n1, n2) where t = n1*b1 + n2*b2")
    _cc_a1xa1.sort(key=lambda x: x[0])
    _w = max(len(e) for _,e,_ in _cc_a1xa1) if _cc_a1xa1 else 4
    print(f'Conjugacy class contains {len(_cc_a1xa1)} element(s) in bounding box:')
    for (n1,n2), elem_str, fc_str in _cc_a1xa1:
        _t_str = f"t_({int(n1)}, {int(n2)}) · {fc_str}"
        print(f'  {elem_str:<{_w}}  =  {_t_str}')
    _a1xa1_conjugacy_class_plot(S1, S2, id1, id2, G, B, H0, bounding_box)


def _Cocon_Elements_A1xA1(e1, e2, e3, e4, bounding_box):
    """Coconjugation in Ã₁ × Ã₁.
    Element 1 = (e1, e2), Element 2 = (e3, e4).
    """
    vertex1 = vector([2,0]); vertex2 = vector([0,2])
    B = matrix([vertex1,vertex2]).transpose()
    coroot_vertices = [vector(i*vertex1+j*vertex2)
                       for i in range(-bounding_box,bounding_box)
                       for j in range(-bounding_box,bounding_box)
                       if abs(i*vertex1[0]+j*vertex2[0])<=bounding_box
                       and abs(i*vertex1[1]+j*vertex2[1])<=bounding_box]
    coroot_vertices_in_B = [vector([ZZ(round((B.inverse()*v)[k].n(),0)) for k in range(2)])
                            for v in coroot_vertices]
    RS1 = RootSystem(["A",1,1]); RS2 = RootSystem(["A",1,1])
    W1 = RS1.root_lattice().weyl_group(); W2 = RS2.root_lattice().weyl_group()
    S1 = W1.simple_reflections(); S2 = W2.simple_reflections()
    id1 = W1.from_reduced_word([]); id2 = W2.from_reduced_word([])
    H0 = {(id1,id2):None, (S1[1],id2):None, (id1,S2[1]):None, (S1[1],S2[1]):None}
    G, H0 = _a1xa1_plot_grid(bounding_box, H0)
    G += point2d(coroot_vertices, size=70, color='black')
    G, tp1, cp_tmp1, cp1 = _a1xa1_chosen_elt(G, e1, e2, bounding_box,
                                              coroot_vertices, B, H0, 'red')
    # Compute element 2 position first to check if same as element 1
    _G_tmp, tp2, cp_tmp2, cp2 = _a1xa1_chosen_elt(Graphics(), e3, e4, bounding_box,
                                                   coroot_vertices, B, H0, 'blue')
    # Only draw element 2 if it is at a different position
    if tp1 != tp2 or cp_tmp1 != cp_tmp2:
        G = G + _G_tmp
    cocon_H0 = _a1xa1_cocon_H0(cp_tmp1, cp_tmp2, H0)
    mod_set_h0p = _a1xa1_mod_set(vertex1, vertex2, B, H0, coroot_vertices_in_B, cp_tmp2)
    lpmul, cocon_H0_lambda = _a1xa1_cocon_H0_lambda(tp1, tp2, cocon_H0, mod_set_h0p, H0)
    if not cocon_H0_lambda:
        print("Coconjugacy set is empty.")
        G.show(frame=False, axes=False, figsize=15, aspect_ratio=1); return
    eta_cocon = []
    for i in range(len(lpmul)):
        _, tmp = _a1xa1_eta(coroot_vertices_in_B, lpmul[i], cocon_H0_lambda[i], cp_tmp2)
        eta_cocon.extend(tmp)
    G = _a1xa1_plot_cocon(G, bounding_box, B, S1, S2, id1, id2, H0, eta_cocon)
    G.show(frame=False, axes=False, figsize=15, aspect_ratio=1)


# ============================================================
# Affine dimension 2  —  shared helpers
# ============================================================

def _aff2_add_stripes_to_triangle(G, vertices, zorder=6, num_stripes=5,
                                   stripe_color='lightgray', stripe_width=1.5):
    v0,v1,v2 = vertices
    for i in range(1, num_stripes+1):
        t = i/(num_stripes+1)
        p1 = v0 + t*(v2-v0); p2 = v1 + t*(v2-v1)
        G += line2d([p1,p2], color=stripe_color, thickness=stripe_width, zorder=zorder)
    return G


def _aff2_get_keyname(classical_part_in_W):
    d = {(): 'identity', (1,):'s1', (2,):'s2', (1,2):'s12', (2,1):'s21',
         (2,1,2):'s212', (1,2,1):'s121', (1,2,1,2):'s1212', (2,1,2,1):'s2121',
         (1,2,1,2,1):'s12121', (2,1,2,1,2):'s21212', (2,1,2,1,2,1):'s212121'}
    return d.get(tuple(classical_part_in_W))


def _aff2_add_H0_labels(G, H0_alcoves_dict, label_fontsize):
    for key in H0_alcoves_dict:
        alcove = H0_alcoves_dict[key]
        verts = alcove['triangle_vertices']
        centroid = (verts[0]+verts[1]+verts[2])/3
        rw = key.reduced_word()
        label = "$e$" if not rw else f"$s_{{{''.join(str(i) for i in rw)}}}$"
        G += text(label,(centroid[0],centroid[1]),fontsize=label_fontsize,color='black',zorder=45)
    return G


def _aff2_build_keys(key_sequences, S_weyl, identity):
    keys = {}
    for name, indices in key_sequences.items():
        if not indices:
            keys[name] = identity
        else:
            result = S_weyl[indices[0]]
            for i in indices[1:]:
                result = result * S_weyl[i]
            keys[name] = result
    return keys


def _aff2_element_weyl_color(S_weyl, cp, tp_euclid, identity, keys, H0dict, G, border_color):
    kn = _aff2_get_keyname(cp)
    alcove = H0dict[keys[kn]]
    verts = alcove['triangle_vertices']; color = alcove['color']
    tv = [vector([verts[i][0]+tp_euclid[0], verts[i][1]+tp_euclid[1]]) for i in range(3)]
    G += polygon2d(tv, fill=True, color=color, alpha=1)
    G += polygon2d(tv, fill=False, color=border_color, thickness=3, zorder=30)
    return G


# ============================================================
# Affine dimension 2  —  geometry data (alcoves, vertices, keys)
# ============================================================

def _aff2_get_geometry(affine_letter, S_weyl, identity):
    """Returns vertex1, vertex2, vertexa, vertexb, keys, alcove_data,
    alcoves list, label_fontsize_fn, epsilon/alignment params."""

    key_sequences = {'identity':[], 's1':[1], 's2':[2],
                     's12':[1,2], 's21':[2,1], 's121':[1,2,1], 's212':[2,1,2]}
    keys = _aff2_build_keys(key_sequences, S_weyl, identity)

    if affine_letter == 'G':
        vertex1 = vector([-4.5,(3*sqrt(3))/2]); vertex2 = vector([3,0])
        vertexa = vector([0,3*sqrt(3)]); vertexb = vector([3/2,3*sqrt(3)/2])
        eps = (-0.1,0.1,0.5,-0.15,'right','bottom','right','bottom')
        keys.update(_aff2_build_keys(
            {'s121':[1,2,1],'s1212':[1,2,1,2],'s2121':[2,1,2,1],
             's12121':[1,2,1,2,1],'s21212':[2,1,2,1,2],'s212121':[2,1,2,1,2,1]}, S_weyl, identity))
        q = sqrt(3)
        alcove_data = [
            (keys['identity'], [vector([0,0]),vector([3/4,3/4*q]),vector([0,q])], (0.7,1.0,1.0)),
            (keys['s1'],   [vector([0,0]),vector([3/2,q/2]),vector([3/4,3/4*q])], (1.0,0.7,1.0)),
            (keys['s2'],   [vector([0,0]),vector([-3/4,3/4*q]),vector([0,q])], (1.0,1.0,0.7)),
            (keys['s12'],  [vector([0,0]),vector([3/2,0]),vector([3/2,q/2])], (0.7,1.0,0.7)),
            (keys['s21'],  [vector([0,0]),vector([-3/2,q/2]),vector([-3/4,3/4*q])], (1.0,0.8,0.5)),
            (keys['s121'], [vector([0,0]),vector([3/2,0]),vector([3/2,-q/2])], (0.8,0.7,1.0)),
            (keys['s212'], [vector([0,0]),vector([-3/2,0]),vector([-3/2,q/2])], (0.7,0.9,0.7)),
            (keys['s1212'],[vector([0,0]),vector([3/2,-q/2]),vector([3/4,-3/4*q])], (0.7,0.8,1.0)),
            (keys['s2121'],[vector([0,0]),vector([-3/2,0]),vector([-3/2,-q/2])], (0.8,0.8,0.8)),
            (keys['s12121'],[vector([0,0]),vector([3/4,-3/4*q]),vector([0,-q])], (1.0,0.85,0.9)),
            (keys['s21212'],[vector([0,0]),vector([-3/2,-q/2]),vector([-3/4,-3/4*q])], (1.0,0.7,0.7)),
            (keys['s212121'],[vector([0,0]),vector([-3/4,-3/4*q]),vector([0,-q])], (0.7,0.7,1.0)),
        ]
        all_keys = ['identity','s1','s2','s12','s21','s121','s212','s1212','s2121','s12121','s21212','s212121']
        alcoves_list = [keys[k] for k in all_keys]
        lf = lambda bb: max(8,min(15,int(-1.2*bb+18)))

    elif affine_letter == 'A':
        vertex1 = vector([1.5,(1/2)*sqrt(3)]); vertex2 = vector([-1.5,(1/2)*sqrt(3)])
        vertexa = vector([0.5,(1/2)*sqrt(3)]); vertexb = vector([-0.5,(1/2)*sqrt(3)])
        eps = (0.2,0.2,-0.2,0.2,'left','center','right','center')
        q = sqrt(3)
        alcove_data = [
            (keys['identity'],[vector([0,0]),vector([1/2,(1/2)*q]),vector([-1/2,(1/2)*q])], (0.7,1.0,1.0)),
            (keys['s1'],   [vector([0,0]),vector([-1/2,(1/2)*q]),vector([-1,0])], (1.0,0.7,1.0)),
            (keys['s2'],   [vector([0,0]),vector([1,0]),vector([1/2,(1/2)*q])], (1.0,1.0,0.7)),
            (keys['s12'],  [vector([0,0]),vector([-1/2,(-1/2)*q]),vector([-1,0])], (0.7,1.0,0.7)),
            (keys['s21'],  [vector([0,0]),vector([1,0]),vector([1/2,(-1/2)*q])], (1.0,0.8,0.5)),
            (keys['s212'], [vector([0,0]),vector([1/2,(-1/2)*q]),vector([-1/2,(-1/2)*q])], (0.8,0.7,1.0)),
        ]
        alcoves_list = [keys[k] for k in ['identity','s1','s2','s12','s21','s212']]
        lf = lambda bb: max(8,min(15,int(-1*bb+20)))

    elif affine_letter == 'B':
        vertex1 = vector([1,-1]); vertex2 = vector([0,2])
        vertexa = vector([1,0]); vertexb = vector([1,1])
        eps = (0.1,-0.1,0.1,0.2,'left','top','right','bottom')
        keys.update(_aff2_build_keys({'s121':[1,2,1],'s2121':[2,1,2,1]}, S_weyl, identity))
        alcove_data = [
            (keys['identity'],[vector([0,0]),vector([1,0]),vector([1/2,1/2])], (0.7,1.0,1.0)),
            (keys['s1'],   [vector([0,0]),vector([0,1]),vector([1/2,1/2])], (1.0,0.7,1.0)),
            (keys['s2'],   [vector([0,0]),vector([1,0]),vector([1/2,-1/2])], (1.0,1.0,0.7)),
            (keys['s12'],  [vector([0,0]),vector([0,1]),vector([-1/2,1/2])], (0.7,1.0,0.7)),
            (keys['s21'],  [vector([0,0]),vector([0,-1]),vector([1/2,-1/2])], (1.0,0.8,0.5)),
            (keys['s121'], [vector([0,0]),vector([-1,0]),vector([-1/2,1/2])], (0.8,0.7,1.0)),
            (keys['s212'], [vector([0,0]),vector([0,-1]),vector([-1/2,-1/2])], (0.7,0.9,0.7)),
            (keys['s2121'],[vector([0,0]),vector([-1,0]),vector([-1/2,-1/2])], (0.8,0.8,0.8)),
        ]
        alcoves_list = [keys[k] for k in ['identity','s1','s2','s12','s21','s121','s212','s2121']]
        lf = lambda bb: max(8,min(15,int(-1.2*bb+18)))

    elif affine_letter == 'C':
        vertex1 = vector([2,-2]); vertex2 = vector([0,2])
        vertexa = vector([2,0]); vertexb = vector([1,1])
        eps = (0.1,-0.1,0.2,0.2,'left','top','right','bottom')
        keys.update(_aff2_build_keys({'s121':[1,2,1],'s2121':[2,1,2,1]}, S_weyl, identity))
        alcove_data = [
            (keys['identity'],[vector([0,0]),vector([1,0]),vector([1,1])], (0.7,1.0,1.0)),
            (keys['s1'],   [vector([0,0]),vector([0,1]),vector([1,1])], (1.0,0.7,1.0)),
            (keys['s2'],   [vector([0,0]),vector([1,0]),vector([1,-1])], (1.0,1.0,0.7)),
            (keys['s12'],  [vector([0,0]),vector([0,1]),vector([-1,1])], (0.7,1.0,0.7)),
            (keys['s21'],  [vector([0,0]),vector([0,-1]),vector([1,-1])], (1.0,0.8,0.5)),
            (keys['s121'], [vector([0,0]),vector([-1,0]),vector([-1,1])], (0.8,0.7,1.0)),
            (keys['s212'], [vector([0,0]),vector([0,-1]),vector([-1,-1])], (0.7,0.9,0.7)),
            (keys['s2121'],[vector([0,0]),vector([-1,0]),vector([-1,-1])], (0.8,0.8,0.8)),
        ]
        alcoves_list = [keys[k] for k in ['identity','s1','s2','s12','s21','s121','s212','s2121']]
        lf = lambda bb: max(8,min(15,int(-1.2*bb+18)))

    else:
        raise ValueError(f"Unknown affine letter '{affine_letter}'. Use 'A','B','C', or 'G'.")

    return vertex1, vertex2, vertexa, vertexb, keys, alcove_data, alcoves_list, lf, eps


# ============================================================
# Affine dimension 2  —  mod set, reflect, finite, eta functions
# ============================================================

def _aff2_mod_set(letter, v1, v2, B, H0dict, crv_in_B, cp_tmp, keys):
    v1B = B.inverse()*v1; v2B = B.inverse()*v2
    VS = ZZ**2
    if letter == 'C':
        spans = {keys['identity']:[vector([0,0])],
                 keys['s1']:VS.span([v1B]), keys['s2']:VS.span([2*v2B]),
                 keys['s12']:VS.span([v1B,2*v2B]), keys['s21']:VS.span([v1B,2*v2B]),
                 keys['s121']:VS.span([2*v1B+2*v2B]),
                 keys['s212']:VS.span([v1B+2*v2B]),
                 keys['s2121']:VS.span([2*v1B,2*v2B])}
    elif letter == 'A':
        spans = {keys['identity']:[vector([0,0])],
                 keys['s1']:VS.span([v1B]), keys['s2']:VS.span([v2B]),
                 keys['s12']:VS.span([v1B+2*v2B,3*v2B]),
                 keys['s21']:VS.span([v1B+2*v2B,3*v2B]),
                 keys['s212']:VS.span([v1B+v2B]),
                 keys['s121']:VS.span([v1B+v2B])}
    elif letter == 'B':
        spans = {keys['identity']:[vector([0,0])],
                 keys['s1']:VS.span([2*v1B]), keys['s2']:VS.span([v2B]),
                 keys['s12']:VS.span([v1B,v2B]), keys['s21']:VS.span([2*v1B,v2B]),
                 keys['s121']:VS.span([2*v1B+v2B]),
                 keys['s212']:VS.span([2*v1B+2*v2B]),
                 keys['s2121']:VS.span([2*v1B,2*v2B])}
    elif letter == 'G':
        spans = {keys['identity']:[vector([0,0])],
                 keys['s1']:VS.span([v1B]), keys['s2']:VS.span([v2B]),
                 keys['s12']:VS.span([v1B,v2B]), keys['s21']:VS.span([v1B,v2B]),
                 keys['s121']:VS.span([v1B+v2B]),
                 keys['s212']:VS.span([v1B+3*v2B]),
                 keys['s1212']:VS.span([v1B,3*v2B]),
                 keys['s2121']:VS.span([v1B,3*v2B]),
                 keys['s12121']:VS.span([2*v1B+3*v2B]),
                 keys['s21212']:VS.span([v1B+2*v2B]),
                 keys['s212121']:VS.span([2*v1B,2*v2B])}
    mod_plots = {}
    for k,sp in spans.items():
        if k == keys['identity']:
            mod_plots[k] = [vector([0,0])]
        else:
            mod_plots[k] = [v for v in crv_in_B if v in sp]
    for k in H0dict:
        if k in mod_plots:
            H0dict[k]['mod_set'] = mod_plots[k]
    return spans[cp_tmp]


def _aff2_reflect_cc(letter, tp, cp_tmp, H0dict, keys):
    lpm = [tp + vector(H0dict[cp_tmp]['mod_set'][i])
           for i in range(len(H0dict[cp_tmp]['mod_set']))]
    if letter == 'G':
        rd = {keys['identity']:lpm,
              keys['s1']:[vector([-v[0]+v[1],v[1]]) for v in lpm],
              keys['s2']:[vector([v[0],3*v[0]-v[1]]) for v in lpm],
              keys['s12']:[vector([2*v[0]-v[1],3*v[0]-v[1]]) for v in lpm],
              keys['s21']:[vector([-v[0]+v[1],-3*v[0]+2*v[1]]) for v in lpm],
              keys['s121']:[vector([-2*v[0]+v[1],-3*v[0]+2*v[1]]) for v in lpm],
              keys['s212']:[vector([2*v[0]-v[1],3*v[0]-2*v[1]]) for v in lpm],
              keys['s1212']:[vector([v[0]-v[1],3*v[0]-2*v[1]]) for v in lpm],
              keys['s2121']:[vector([-2*v[0]+v[1],-3*v[0]+v[1]]) for v in lpm],
              keys['s12121']:[vector([-v[0],-3*v[0]+v[1]]) for v in lpm],
              keys['s21212']:[vector([v[0]-v[1],-v[1]]) for v in lpm],
              keys['s212121']:[vector([-v[0],-v[1]]) for v in lpm]}
    elif letter == 'A':
        rd = {keys['identity']:lpm,
              keys['s1']:[vector([v[1]-v[0],v[1]]) for v in lpm],
              keys['s2']:[vector([v[0],v[0]-v[1]]) for v in lpm],
              keys['s12']:[vector([-v[1],v[0]-v[1]]) for v in lpm],
              keys['s21']:[vector([v[1]-v[0],-v[0]]) for v in lpm],
              keys['s212']:[vector([-v[1],-v[0]]) for v in lpm]}
    elif letter == 'B':
        rd = {keys['identity']:lpm,
              keys['s1']:[vector([-v[0]+2*v[1],v[1]]) for v in lpm],
              keys['s2']:[vector([v[0],v[0]-v[1]]) for v in lpm],
              keys['s12']:[vector([v[0]-2*v[1],v[0]-v[1]]) for v in lpm],
              keys['s21']:[vector([-v[0]+2*v[1],-v[0]+v[1]]) for v in lpm],
              keys['s121']:[vector([-v[0],-v[0]+v[1]]) for v in lpm],
              keys['s212']:[vector([v[0]-2*v[1],-v[1]]) for v in lpm],
              keys['s2121']:[vector([-v[0],-v[1]]) for v in lpm]}
    elif letter == 'C':
        rd = {keys['identity']:lpm,
              keys['s1']:[vector([-v[0]+v[1],v[1]]) for v in lpm],
              keys['s2']:[vector([v[0],2*v[0]-v[1]]) for v in lpm],
              keys['s12']:[vector([v[0]-v[1],2*v[0]-v[1]]) for v in lpm],
              keys['s21']:[vector([-v[0]+v[1],-2*v[0]+v[1]]) for v in lpm],
              keys['s121']:[vector([-v[0],-2*v[0]+v[1]]) for v in lpm],
              keys['s212']:[vector([v[0]-v[1],-v[1]]) for v in lpm],
              keys['s2121']:[vector([-v[0],-v[1]]) for v in lpm]}
    for k in H0dict:
        if k in rd:
            H0dict[k]['reflect'] = rd[k]
    return H0dict


def _aff2_finite_cc(cp_tmp, H0dict):
    w = cp_tmp
    for g in H0dict.keys():
        c = g*w*g.inverse()
        H0dict[g]['finite_conjugate'] = tuple(c.reduced_word())
    return H0dict


def _aff2_eta(letter, crv_in_B, lpmul_i, cocon_elem, cp2, keys):
    l0 = lpmul_i[0]; l1 = lpmul_i[1]
    kn = _aff2_get_keyname(cp2)
    eta = []
    if kn == 'identity':
        eta = crv_in_B if (l0==0 and l1==0) else []
    elif letter == 'C':
        if kn=='s1':   eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1]
        elif kn=='s2': eta=[v for v in crv_in_B if 0==l0 and -2*v[0]+2*v[1]==l1]
        elif kn=='s12':eta=[v for v in crv_in_B if v[1]==l0 and -2*v[0]+2*v[1]==l1]
        elif kn=='s21':eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 2*v[0]==l1]
        elif kn=='s121':eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[0]==l1]
        elif kn=='s212':eta=[v for v in crv_in_B if v[1]==l0 and 2*v[1]==l1]
        elif kn=='s2121':eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[1]==l1]
    elif letter == 'A':
        if kn=='s1':   eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1]
        elif kn=='s2': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]==l1]
        elif kn=='s12':eta=[v for v in crv_in_B if v[0]+v[1]==l0 and -v[0]+2*v[1]==l1]
        elif kn=='s21':eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[1]==l1]
        elif kn=='s212':eta=[v for v in crv_in_B if v[0]+v[1]==l0 and v[0]+v[1]==l1]
    elif letter == 'B':
        if kn=='s1':   eta=[v for v in crv_in_B if 2*v[0]-2*v[1]==l0 and 0==l1]
        elif kn=='s2': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]==l1]
        elif kn=='s12':eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+3*v[1]==l1]
        elif kn=='s21':eta=[v for v in crv_in_B if 2*v[0]-2*v[1]==l0 and v[0]==l1]
        elif kn=='s121':eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]==l1]
        elif kn=='s212':eta=[v for v in crv_in_B if 2*v[1]==l0 and 2*v[1]==l1]
        elif kn=='s2121':eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[1]==l1]
    elif letter == 'G':
        if kn=='s1':   eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1]
        elif kn=='s2': eta=[v for v in crv_in_B if 0==l0 and -3*v[0]+2*v[1]==l1]
        elif kn=='s12':eta=[v for v in crv_in_B if -v[0]+v[1]==l0 and -3*v[0]+2*v[1]==l1]
        elif kn=='s21':eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 3*v[0]-v[1]==l1]
        elif kn=='s121':eta=[v for v in crv_in_B if 3*v[0]-v[1]==l0 and 3*v[0]-v[1]==l1]
        elif kn=='s212':eta=[v for v in crv_in_B if -v[0]+v[1]==l0 and -3*v[0]+3*v[1]==l1]
        elif kn=='s1212':eta=[v for v in crv_in_B if v[1]==l0 and -3*v[0]+3*v[1]==l1]
        elif kn=='s2121':eta=[v for v in crv_in_B if 3*v[0]-v[1]==l0 and 3*v[0]==l1]
        elif kn=='s12121':eta=[v for v in crv_in_B if 2*v[0]==l0 and 3*v[0]==l1]
        elif kn=='s21212':eta=[v for v in crv_in_B if v[1]==l0 and 2*v[1]==l1]
        elif kn=='s212121':eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[1]==l1]
    return eta, [[v, cocon_elem] for v in eta]


def _aff2_cocon_H0(cp1, cp2, H0dict):
    cocon = []
    for g in H0dict.keys():
        c = g*cp1*g.inverse()
        if c.reduced_word() == cp2.reduced_word():
            cocon.append(g)
    return cocon


def _aff2_cocon_H0_lambda(letter, tp1, tp2, cocon_H0, mod_set_h0, H0dict, keys):
    if letter == 'C':
        rd = {keys['identity']:tp1,
              keys['s1']:vector([-tp1[0]+tp1[1],tp1[1]]),
              keys['s2']:vector([tp1[0],2*tp1[0]-tp1[1]]),
              keys['s12']:vector([tp1[0]-tp1[1],2*tp1[0]-tp1[1]]),
              keys['s21']:vector([-tp1[0]+tp1[1],-2*tp1[0]+tp1[1]]),
              keys['s121']:vector([-tp1[0],-2*tp1[0]+tp1[1]]),
              keys['s212']:vector([tp1[0]-tp1[1],-tp1[1]]),
              keys['s2121']:vector([-tp1[0],-tp1[1]])}
    elif letter == 'A':
        rd = {keys['identity']:tp1,
              keys['s1']:vector([tp1[1]-tp1[0],tp1[1]]),
              keys['s2']:vector([tp1[0],tp1[0]-tp1[1]]),
              keys['s12']:vector([-tp1[1],tp1[0]-tp1[1]]),
              keys['s21']:vector([tp1[1]-tp1[0],-tp1[0]]),
              keys['s212']:vector([-tp1[1],-tp1[0]])}
    elif letter == 'B':
        rd = {keys['identity']:tp1,
              keys['s1']:vector([-tp1[0]+2*tp1[1],tp1[1]]),
              keys['s2']:vector([tp1[0],tp1[0]-tp1[1]]),
              keys['s12']:vector([tp1[0]-2*tp1[1],tp1[0]-tp1[1]]),
              keys['s21']:vector([-tp1[0]+2*tp1[1],-tp1[0]+tp1[1]]),
              keys['s121']:vector([-tp1[0],-tp1[0]+tp1[1]]),
              keys['s212']:vector([tp1[0]-2*tp1[1],-tp1[1]]),
              keys['s2121']:vector([-tp1[0],-tp1[1]])}
    elif letter == 'G':
        rd = {keys['identity']:tp1,
              keys['s1']:vector([-tp1[0]+tp1[1],tp1[1]]),
              keys['s2']:vector([tp1[0],3*tp1[0]-tp1[1]]),
              keys['s12']:vector([2*tp1[0]-tp1[1],3*tp1[0]-tp1[1]]),
              keys['s21']:vector([-tp1[0]+tp1[1],-3*tp1[0]+2*tp1[1]]),
              keys['s121']:vector([-2*tp1[0]+tp1[1],-3*tp1[0]+2*tp1[1]]),
              keys['s212']:vector([2*tp1[0]-tp1[1],3*tp1[0]-2*tp1[1]]),
              keys['s1212']:vector([tp1[0]-tp1[1],3*tp1[0]-2*tp1[1]]),
              keys['s2121']:vector([-2*tp1[0]+tp1[1],-3*tp1[0]+tp1[1]]),
              keys['s12121']:vector([-tp1[0],-3*tp1[0]+tp1[1]]),
              keys['s21212']:vector([tp1[0]-tp1[1],-tp1[1]]),
              keys['s212121']:vector([-tp1[0],-tp1[1]])}
    cocon_lambda = []; lpmul = []
    for g in cocon_H0:
        ulambda = rd[g]; diff = tp2 - ulambda
        if diff in mod_set_h0:
            cocon_lambda.append(g); lpmul.append(diff)
    return lpmul, cocon_lambda


# ============================================================
# Affine dimension 2  —  plot functions
# ============================================================

def _aff2_conjugacy_class_plot(letter, S_weyl, identity, G, B, H0dict, bounding_box,
                                plot_alcoves, plot_background):
    plotted = set()
    for k in H0dict.keys():
        tmp = []
        for vec in H0dict[k]['reflect']:
            check = B*vec
            if abs(check[0])<=bounding_box and abs(check[1])<=bounding_box:
                tmp.append(check)
        H0dict[k]['reflect_filtered'] = tmp
    h0_stripe = set()
    for k in H0dict.keys():
        fc = H0dict[k]['finite_conjugate']
        for vec in H0dict[k]['reflect_filtered']:
            vx = float(vec[0]); vy = float(vec[1])
            if abs(vx)<0.001 and abs(vy)<0.001:
                h0_stripe.add(fc); break
    for k in H0dict.keys():
        fc = H0dict[k]['finite_conjugate']
        if not fc: ak = identity
        else:
            ak = S_weyl[fc[0]]
            for idx in fc[1:]: ak = ak * S_weyl[idx]
        alcove = H0dict[ak]; verts = alcove['triangle_vertices']; color = alcove['color']
        rw_t = tuple(ak.reduced_word())
        for vec in H0dict[k]['reflect_filtered']:
            ok = all(abs(verts[i][j]+vec[j])<=bounding_box for i in range(3) for j in range(2))
            if ok:
                tv = [vector([verts[i][0]+vec[0],verts[i][1]+vec[1]]) for i in range(3)]
                G += polygon2d(tv, fill=True, color=color, alpha=1, zorder=1)
                vx = float(vec[0]); vy = float(vec[1])
                if abs(vx)<0.001 and abs(vy)<0.001:
                    G = _aff2_add_stripes_to_triangle(G, tv, zorder=2)
                plotted.add((round(vx,6), round(vy,6), rw_t))
    for key in H0dict:
        rw = tuple(key.reduced_word())
        if (0.0,0.0,rw) not in plotted:
            alcove = H0dict[key]; verts = alcove['triangle_vertices']; color = alcove['color']
            G += polygon2d(verts, fill=True, color=color, alpha=1, zorder=1)
    G += plot_alcoves + plot_background
    G.show(frame=False, axes=False, figsize=12)


def _aff2_plot_cocon(G, bounding_box, B, S_weyl, identity, H0dict, eta_cocon,
                     plot_alcoves, plot_background):
    def rw_tuple(item):
        return tuple(item.reduced_word()) if hasattr(item,'reduced_word') else tuple(item)
    plotted = set()
    filtered = []
    for elt in eta_cocon:
        vec = elt[0]; rw = rw_tuple(elt[1])
        check = B*vec
        if abs(check[0])<=bounding_box and abs(check[1])<=bounding_box:
            filtered.append([check, rw])
    for elt in filtered:
        rw = elt[1]
        if not rw: ak = identity
        else:
            ak = S_weyl[rw[0]]
            for idx in rw[1:]: ak = ak*S_weyl[idx]
        alcove = H0dict[ak]; verts = alcove['triangle_vertices']; color = alcove['color']
        vec = elt[0]
        ok = all(abs(verts[i][j]+vec[j])<=bounding_box for i in range(3) for j in range(2))
        if ok:
            tv = [vector([verts[i][0]+vec[0],verts[i][1]+vec[1]]) for i in range(3)]
            G += polygon2d(tv, fill=True, color=color, alpha=1, zorder=1)
            if vec[0]==0 and vec[1]==0:
                G = _aff2_add_stripes_to_triangle(G, tv, zorder=2)
            plotted.add((float(vec[0]),float(vec[1]),rw))
    for key in H0dict:
        rw = tuple(key.reduced_word())
        if (0.0,0.0,rw) not in plotted:
            alcove = H0dict[key]; verts = alcove['triangle_vertices']; color = alcove['color']
            G += polygon2d(verts, fill=True, color=color, alpha=1, zorder=1)
    G += plot_alcoves + plot_background
    G.show(frame=False, axes=False, figsize=12)


# ============================================================
# Affine dimension 2  —  main computation functions
# ============================================================

def _conjugacy_class_affine2(letter, elt, bounding_box):
    """Conjugacy class for affine dimension-2 groups (letter in A, B, C, G)."""
    E = ExtendedAffineWeylGroup([letter,2,1])
    PW0 = E.PW0(); FW = E.FW(); S_ext = FW.simple_reflections()
    L = RootSystem([letter,2,1]).ambient_space()
    W = L.weyl_group(); S_weyl = W.simple_reflections()
    identity = W.from_reduced_word([])
    # Check for t_(n1,n2)*s_XX format
    _parsed_aff2 = _parse_translation_input(elt)
    if _parsed_aff2 is not None:
        _coords2, _word2 = _parsed_aff2
        elt = _word2  # use just the word part for downstream
        # Override tp after geometry is set up — handled below
    ew = W.from_reduced_word([int(c) for c in elt])
    w = ew.reduced_word(); str_w = ''.join(str(x) for x in w)
    if str_w:
        ee = S_ext[int(str_w[0])]
        for i in range(1,len(str_w)): ee = ee * S_ext[int(str_w[i])]
    else:
        ee = FW.one()
    w_pw0 = PW0(ee); tmp_tp, cp = w_pw0.value
    cp_tmp = W.from_reduced_word(cp.reduced_word())
    cp_rw = tuple(cp_tmp.reduced_word())  # ensure plain tuple of ints
    vertex1, vertex2, vertexa, vertexb, keys, alcove_data, alcoves_list, lf, eps = \
        _aff2_get_geometry(letter, S_weyl, identity)
    B = matrix([vertex1,vertex2]).transpose()
    C = matrix([vertexa,vertexb]).transpose()
    coroot_vertices = [vector(i*vertex1+j*vertex2)
                       for i in range(-bounding_box,bounding_box)
                       for j in range(-bounding_box,bounding_box)
                       if abs(i*vertex1[0]+j*vertex2[0])<=bounding_box
                       and abs(i*vertex1[1]+j*vertex2[1])<=bounding_box]
    crv_in_B = [vector([ZZ(round((B.inverse()*v)[k].n(),0)) for k in range(2)])
                for v in coroot_vertices]
    H0dict = {}
    for key, verts, color in alcove_data:
        H0dict[key] = {'triangle_vertices': verts, 'color': color}
    H0_alcoves = [polygon2d(H0dict[k]['triangle_vertices'],fill=True,color=H0dict[k]['color'],alpha=1)
                  for k in H0dict]
    plot_alcoves = L.plot(fundamental_chamber=False,roots=False,coroots=False,
                          fundamental_weights=False,reflection_hyperplanes=False,
                          alcoves=alcoves_list,alcove_labels=False,bounding_box=5)
    plot_background = L.plot(fundamental_chamber=True,roots=False,coroots=False,
                             fundamental_weights=False,bounding_box=bounding_box,
                             reflection_hyperplanes=False)
    e1a,e1b,e2a,e2b,ha1,va1,ha2,va2 = eps
    G  = arrow((0,0),tuple(vertex1),color='red',width=3,zorder=2,arrowsize=4)
    G += arrow((0,0),tuple(vertex2),color='blue',width=3,zorder=2,arrowsize=4)
    G += text(r"$\alpha^\vee_1$",(vertex1[0]+e1a,vertex1[1]+e1b),fontsize=14,color='red',
              horizontal_alignment=ha1,vertical_alignment=va1)
    G += text(r"$\alpha^\vee_2$",(vertex2[0]+e2a,vertex2[1]+e2b),fontsize=14,color='blue',
              horizontal_alignment=ha2,vertical_alignment=va2)
    label_fontsize = lf(bounding_box)
    G += point2d(coroot_vertices,size=70,color='black')
    for h in H0_alcoves: G += h
    G = _aff2_add_H0_labels(G, H0dict, label_fontsize)
    tmp_le = tmp_tp.value
    tp_euclid = C*tmp_le.to_vector()
    tmp_v = B.inverse()*C*tmp_le.to_vector()
    tp = vector([ZZ(round(c.n(),0)) for c in tmp_v])
    # Override tp/tp_euclid if translation format was entered
    if _parsed_aff2 is not None:
        _n1, _n2 = _parsed_aff2[0][0], _parsed_aff2[0][1]
        tp = vector([ZZ(_n1), ZZ(_n2)])
        tp_euclid = B * tp
    if len(w)>0:
        indices = "".join(str(idx) for idx in w)
        lbl = f"$s_{{{indices}}}$"
        alc_v = H0dict[cp_tmp]['triangle_vertices']
        centroid = (alc_v[0]+alc_v[1]+alc_v[2])/3
        lp = tp_euclid + centroid
        G += text(lbl,(lp[0],lp[1]),fontsize=label_fontsize,color='black',zorder=45)
    G = _aff2_element_weyl_color(S_weyl, cp_rw, tp_euclid, identity, keys, H0dict, G, 'red')
    # If input is identity, add stripes to chosen element alcove
    if len(w)==0 and tp[0]==0 and tp[1]==0:
        _id_verts = H0dict[identity]['triangle_vertices']
        _tv_id = [vector([_id_verts[_i][0]+tp_euclid[0], _id_verts[_i][1]+tp_euclid[1]]) for _i in range(3)]
        G = _aff2_add_stripes_to_triangle(G, _tv_id, zorder=31)
    _aff2_mod_set(letter, vertex1, vertex2, B, H0dict, crv_in_B, cp_tmp, keys)
    _aff2_reflect_cc(letter, tp, cp_tmp, H0dict, keys)
    _aff2_finite_cc(cp_tmp, H0dict)
    # Print conjugacy class elements using PW0.from_translation (original notebook approach)
    _E_aff2 = ExtendedAffineWeylGroup([letter, 2, 1])
    _PW0_aff2 = _E_aff2.PW0(); _FW_aff2 = _E_aff2.FW(); _P_aff2 = _E_aff2.lattice()
    _S_FW_aff2 = _FW_aff2.simple_reflections()
    _W_aff2 = WeylGroup([letter, 2, 1]); _S_waff2 = _W_aff2.simple_reflections()
    # Print basis B first
    def _fmt_val(x):
        s = str(x)
        # If it already contains sqrt, pi, or a fraction, keep it exact
        if any(c in s for c in ['sqrt', 'pi', '/', 'e^']):
            return s
        try:
            f = float(x)
            return f"{f:.10g}"
        except Exception:
            return s
    def _fmt_vec(v):
        return [_fmt_val(x) for x in v]
    # Print input element info
    _rw_in = "s_" + "".join(str(x) for x in w) if w else "e"
    _fc_in = "s_" + "".join(str(x) for x in cp_rw) if cp_rw else "e"
    _tp_in = f"t_({tp[0]},{tp[1]})"
    print(f"Input element:  {_rw_in}  =  {_tp_in} · {_fc_in}")
    # Print basis B and elements as decomposition (no rw reconstruction)
    def _fmt_val_a2(x):
        s = str(x)
        if any(c in s for c in ["sqrt","pi","/"]): return s
        f = float(x)
        return str(int(f)) if f == int(f) else f"{f:.10g}"
    print(f"Coroot Basis (columns = simple coroots):")
    print(f"  b1 = ({', '.join([_fmt_val_a2(x) for x in vertex1])})")
    print(f"  b2 = ({', '.join([_fmt_val_a2(x) for x in vertex2])})")
    print(f"Translation vectors are in basis B coordinates (n1, n2) where t = n1*b1 + n2*b2")
    _E_a2 = ExtendedAffineWeylGroup([letter,2,1]); _PW0_a2 = _E_a2.PW0()
    _FW_a2 = _E_a2.FW(); _P_a2 = _E_a2.lattice(); _S_FW_a2 = _FW_a2.simple_reflections()
    _W_a2 = WeylGroup([letter,2,1]); _S_a2 = _W_a2.simple_reflections()
    def _build_rw_aff2(fc, a, b):
        if letter == "A": om1=2*a-b; om2=-a+2*b
        elif letter == "B": om1=2*a-b; om2=-2*a+2*b
        elif letter == "C": om1=2*a-2*b; om2=-a+2*b
        else:               om1=2*a-b; om2=-3*a+2*b  # G
        if fc:
            _ff = _S_FW_a2[fc[0]]
            for _x in fc[1:]: _ff = _ff * _S_FW_a2[_x]
        else: _ff = _FW_a2.one()
        _tl = _P_a2.zero()
        if om1: _tl = _tl + om1*_P_a2.fundamental_weight(1)
        if om2: _tl = _tl + om2*_P_a2.fundamental_weight(2)
        _ae_fw = list(_FW_a2(_PW0_a2.from_translation(_tl)*_PW0_a2(_ff)).to_affine_weyl_left().reduced_word())
        _ae_w = _W_a2.one()
        for _x in _ae_fw: _ae_w = _ae_w * _S_a2[_x]
        _rw = list(_ae_w.reduced_word())
        return "s_"+"".join(str(x) for x in _rw) if _rw else "e"
    _cc_aff2 = []
    _seen_aff2 = set()
    for _k in H0dict.keys():
        _fc = H0dict[_k].get("finite_conjugate", ())
        _fc_str = "s_" + "".join(str(x) for x in _fc) if _fc else "e"
        for _vec in H0dict[_k].get("reflect", []):
            _chk = B * _vec
            if abs(float(_chk[0])) > bounding_box or abs(float(_chk[1])) > bounding_box: continue
            _a, _b = int(round(float(_vec[0]),0)), int(round(float(_vec[1]),0))
            _key = (_a, _b, _fc)
            if _key in _seen_aff2: continue
            _seen_aff2.add(_key)
            _rw_str = _build_rw_aff2(_fc, _a, _b)
            _cc_aff2.append((_rw_str, f"t_({_a},{_b}) · {_fc_str}"))
    _cc_aff2.sort(key=lambda x:(len(x[0]),x[0]))
    _w2 = max(len(r) for r,_ in _cc_aff2) if _cc_aff2 else 4
    print(f"Conjugacy class contains {len(_cc_aff2)} element(s) in bounding box:")
    for _rw_str, _ts in _cc_aff2:
        print(f"  {_rw_str:<{_w2}}  =  {_ts}")
    _aff2_conjugacy_class_plot(letter, S_weyl, identity, G, B, H0dict, bounding_box,
                                plot_alcoves, plot_background)


def _Cocon_Elements_affine2(letter, elt1, elt2, bounding_box):
    """Coconjugation for affine dimension-2 groups (letter in A, B, C, G)."""
    E = ExtendedAffineWeylGroup([letter,2,1])
    PW0 = E.PW0(); FW = E.FW(); S_ext = FW.simple_reflections()
    L = RootSystem([letter,2,1]).ambient_space()
    W = L.weyl_group(); S_weyl = W.simple_reflections()
    identity = W.from_reduced_word([])
    vertex1, vertex2, vertexa, vertexb, keys, alcove_data, alcoves_list, lf, eps = \
        _aff2_get_geometry(letter, S_weyl, identity)
    B = matrix([vertex1,vertex2]).transpose()
    C = matrix([vertexa,vertexb]).transpose()
    coroot_vertices = [vector(i*vertex1+j*vertex2)
                       for i in range(-bounding_box,bounding_box)
                       for j in range(-bounding_box,bounding_box)
                       if abs(i*vertex1[0]+j*vertex2[0])<=bounding_box
                       and abs(i*vertex1[1]+j*vertex2[1])<=bounding_box]
    crv_in_B = [vector([ZZ(round((B.inverse()*v)[k].n(),0)) for k in range(2)])
                for v in coroot_vertices]
    H0dict = {}
    for key, verts, color in alcove_data:
        H0dict[key] = {'triangle_vertices': verts, 'color': color}
    H0_alcoves = [polygon2d(H0dict[k]['triangle_vertices'],fill=True,color=H0dict[k]['color'],alpha=1)
                  for k in H0dict]
    plot_alcoves = L.plot(fundamental_chamber=False,roots=False,coroots=False,
                          fundamental_weights=False,reflection_hyperplanes=False,
                          alcoves=alcoves_list,alcove_labels=False,bounding_box=5)
    plot_background = L.plot(fundamental_chamber=True,roots=False,coroots=False,
                             fundamental_weights=False,bounding_box=bounding_box,
                             reflection_hyperplanes=False)
    e1a,e1b,e2a,e2b,ha1,va1,ha2,va2 = eps
    label_fontsize = lf(bounding_box)
    G  = arrow((0,0),tuple(vertex1),color='red',width=3,zorder=2,arrowsize=4)
    G += arrow((0,0),tuple(vertex2),color='blue',width=3,zorder=2,arrowsize=4)
    G += text(r"$\alpha^\vee_1$",(vertex1[0]+e1a,vertex1[1]+e1b),fontsize=14,color='red',
              horizontal_alignment=ha1,vertical_alignment=va1)
    G += text(r"$\alpha^\vee_2$",(vertex2[0]+e2a,vertex2[1]+e2b),fontsize=14,color='blue',
              horizontal_alignment=ha2,vertical_alignment=va2)
    G += point2d(coroot_vertices,size=70,color='black')
    for h in H0_alcoves: G += h
    G = _aff2_add_H0_labels(G, H0dict, label_fontsize)

    def _parse_elt(e):
        # Check for t_(n1,n2)*s_XX format
        _parsed_e = _parse_translation_input(e) if e else None
        if _parsed_e is not None:
            _coords_e, _word_e = _parsed_e
            e = _word_e
        ew = W.from_reduced_word([Integer(c) for c in e])
        w = ew.reduced_word(); str_w = ''.join(str(x) for x in w)
        if str_w:
            ee = S_ext[int(str_w[0])]
            for i in range(1,len(str_w)): ee = ee * S_ext[int(str_w[i])]
        else:
            ee = FW.one()
        w_pw0 = PW0(ee); tmp_tp, cp = w_pw0.value
        cp_tmp = W.from_reduced_word(cp.reduced_word()); cp_rw = cp_tmp.reduced_word()
        tmp_le = tmp_tp.value
        tp_euclid = C*tmp_le.to_vector()
        tmp_v = B.inverse()*C*tmp_le.to_vector()
        tp = vector([ZZ(round(c.n(),0)) for c in tmp_v])
        # Override tp if translation format was entered
        if _parsed_e is not None:
            _n1_e, _n2_e = _parsed_e[0][0], _parsed_e[0][1]
            tp = vector([ZZ(_n1_e), ZZ(_n2_e)])
            tp_euclid = B * tp
        return ew, w, cp_tmp, cp_rw, tp_euclid, tp

    ew1, w1, cp1, cp_rw1, tpe1, tp1 = _parse_elt(elt1)
    ew2, w2, cp2, cp_rw2, tpe2, tp2 = _parse_elt(elt2)

    # Draw elements
    for ew, w, cp_tmp, cp_rw, tpe, tp, color in [
            (ew1,w1,cp1,cp_rw1,tpe1,tp1,'red'),
            (ew2,w2,cp2,cp_rw2,tpe2,tp2,'blue')]:
        if len(w)>0:
            idx_str = "".join(str(i) for i in w)
            lbl = f"$s_{{{idx_str}}}$"
            av = H0dict[cp_tmp]['triangle_vertices']
            c3 = (av[0]+av[1]+av[2])/3
            G += text(lbl,(tpe[0]+c3[0],tpe[1]+c3[1]),fontsize=label_fontsize,color='black',zorder=45)
        G = _aff2_element_weyl_color(S_weyl, cp_rw, tpe, identity, keys, H0dict, G, color)

    mod_set_h0p = _aff2_mod_set(letter, vertex1, vertex2, B, H0dict, crv_in_B, cp2, keys)
    cocon_H0 = _aff2_cocon_H0(cp1, cp2, H0dict)
    lpmul, cocon_H0_lambda = _aff2_cocon_H0_lambda(letter, tp1, tp2, cocon_H0, mod_set_h0p, H0dict, keys)
    if not cocon_H0_lambda:
        print("Coconjugacy set is empty.")
        G += plot_alcoves + plot_background
        G.show(frame=False, axes=False, figsize=12); return
    eta_cocon = []
    for i in range(len(lpmul)):
        _, tmp = _aff2_eta(letter, crv_in_B, lpmul[i], cocon_H0_lambda[i], cp_rw2, keys)
        eta_cocon.extend(tmp)
    # Print input elements and coconjugation set
    _E_cc2 = ExtendedAffineWeylGroup([letter,2,1]); _PW0_cc2 = _E_cc2.PW0()
    _FW_cc2 = _E_cc2.FW(); _P_cc2 = _E_cc2.lattice(); _S_FW_cc2 = _FW_cc2.simple_reflections()
    _W_cc2 = WeylGroup([letter,2,1]); _S_cc2 = _W_cc2.simple_reflections()
    def _build_rw_cc2(fc, a, b):
        if fc:
            _ff = _S_FW_cc2[fc[0]]
            for _x in fc[1:]: _ff = _ff * _S_FW_cc2[_x]
        else: _ff = _FW_cc2.one()
        _tl = _P_cc2.zero()
        if letter == "A": _om1 = 2*a-b; _om2 = -a+2*b
        elif letter == "B": _om1 = 2*a-b; _om2 = -2*a+2*b
        elif letter == "C": _om1 = 2*a-2*b; _om2 = -a+2*b
        else: _om1 = 2*a-b; _om2 = -3*a+2*b  # G
        if _om1: _tl = _tl + _om1*_P_cc2.fundamental_weight(1)
        if _om2: _tl = _tl + _om2*_P_cc2.fundamental_weight(2)
        _ae = list(_FW_cc2(_PW0_cc2.from_translation(_tl)*_PW0_cc2(_ff)).to_affine_weyl_left().reduced_word())
        _aw = _W_cc2.one()
        for _x in _ae: _aw = _aw * _S_cc2[_x]
        _r = list(_aw.reduced_word())
        return "s_"+"".join(str(x) for x in _r) if _r else "e"
    _rw1_str = "s_"+"".join(str(x) for x in w1) if w1 else "e"
    _rw2_str = "s_"+"".join(str(x) for x in w2) if w2 else "e"
    _fc1_str = "s_"+"".join(str(x) for x in cp_rw1) if cp_rw1 else "e"
    _fc2_str = "s_"+"".join(str(x) for x in cp_rw2) if cp_rw2 else "e"
    print(f"Input element 1: {_rw1_str}  =  t_({tp1[0]},{tp1[1]}) · {_fc1_str}")
    print(f"Input element 2: {_rw2_str}  =  t_({tp2[0]},{tp2[1]}) · {_fc2_str}")
    _cc_elems2 = []
    _seen_cc2 = set()
    for _eta in eta_cocon:
        _vec, _cp = _eta[0], _eta[1]
        _fc = tuple(_cp.reduced_word()) if hasattr(_cp,'reduced_word') else tuple(_cp)
        _chk = B * _vec
        if abs(float(_chk[0])) > bounding_box or abs(float(_chk[1])) > bounding_box: continue
        _a2, _b2 = int(round(float(_vec[0]),0)), int(round(float(_vec[1]),0))
        _key2 = (_a2, _b2, _fc)
        if _key2 in _seen_cc2: continue
        _seen_cc2.add(_key2)
        _fc_s = "s_"+"".join(str(x) for x in _fc) if _fc else "e"
        _rw_s = _build_rw_cc2(_fc, _a2, _b2)
        _cc_elems2.append((_rw_s, f"t_({_a2},{_b2}) · {_fc_s}"))
    _cc_elems2.sort(key=lambda x:(len(x[0]),x[0]))
    _w2p = max(len(r) for r,_ in _cc_elems2) if _cc_elems2 else 4
    print(f"Coconjugation set contains {len(_cc_elems2)} element(s) in bounding box:")
    for _rs,_ts in _cc_elems2:
        print(f"  {_rs:<{_w2p}}  =  {_ts}")
    _aff2_plot_cocon(G, bounding_box, B, S_weyl, identity, H0dict, eta_cocon,
                     plot_alcoves, plot_background)


# ============================================================
# Dimension 3  —  shared utilities (B3 and C3 share most code)
# ============================================================

def _r3_get_coroot_vertices_in_B(bounding_box):
    """Coroot lattice points in B-coordinates for dimension-3 groups."""
    verts = []
    for a in range(-bounding_box, bounding_box+1):
        for b in range(-bounding_box, bounding_box+1):
            for c in range(-bounding_box, bounding_box+1):
                verts.append(vector([a,b,c]))
    return verts


def _r3_format_label(reduced_word):
    subs = {'0':'₀','1':'₁','2':'₂','3':'₃'}
    if not reduced_word: return "e"
    return "s" + ''.join(subs[str(i)] for i in reduced_word)


def _r3_is_in_bounds(verts_3d, bbox):
    return all(abs(v[j]) <= bbox for v in verts_3d for j in range(3))


def _r3_add_stripes_to_tetrahedron(G, vertices, faces,
                                    stripe_color='#333333', num_stripes=4, stripe_width=3):
    for face in faces:
        v0,v1,v2 = [vertices[i] for i in face]
        for i in range(1, num_stripes+1):
            t = i/(num_stripes+1)
            p1 = v0 + t*(v2-v0); p2 = v1 + t*(v2-v1)
            G += line3d([p1,p2], color=stripe_color, thickness=stripe_width)
    return G


# ============================================================
# Ã₃  —  setup and helper functions
# ============================================================

def _a3_to_3d(vector_4_dim):
    v = vector(SR, vector_4_dim)
    v = v - sum(v)/4 * vector(SR, [1,1,1,1])
    e1 = vector(SR,[1,-1,0,0])/sqrt(2)
    e2 = vector(SR,[1,1,-2,0])/sqrt(6)
    e3 = vector(SR,[1,1,1,-3])/sqrt(12)
    return vector(SR,[v.dot_product(e1),v.dot_product(e2),v.dot_product(e3)])


def _a3_setup():
    R_fin = RootSystem(["A",3]); W_fin = R_fin.ambient_space().weyl_group()
    ambient = R_fin.ambient_space()
    E = ExtendedAffineWeylGroup(["A",3,1])
    PW0 = E.PW0(); FW = E.FW(); S_FW = FW.simple_reflections(); P = E.lattice()
    L = RootSystem(["A",3,1]).ambient_space()
    W_aff = WeylGroup(["A",3,1]); S_aff = W_aff.simple_reflections()
    simple_coroots = [ambient.simple_coroot(i) for i in [1,2,3]]
    coroots_3d = [vector(_a3_to_3d(list(cr.to_vector()))) for cr in simple_coroots]
    B = matrix(coroots_3d).transpose()
    fund_weights = [ambient.fundamental_weight(i) for i in [1,2,3]]
    fw_3d = [vector(_a3_to_3d(list(fw.to_vector()))) for fw in fund_weights]
    C = matrix(fw_3d).transpose()
    origin = ambient.zero()
    fw_ambient = [ambient.fundamental_weight(i) for i in [1,2,3]]
    alcove_verts = [origin] + fw_ambient
    S_weyl = W_fin.simple_reflections(); identity = W_fin.one()
    keys = {'identity': identity}
    for w in W_fin:
        rw = tuple(w.reduced_word())
        if len(rw)>0: keys['s'+''.join(str(x) for x in rw)] = w
    return {'R_finite':R_fin,'W_finite':W_fin,'ambient':ambient,'E':E,'PW0':PW0,'FW':FW,
            'S_FW':S_FW,'P':P,'L':L,'W_affine':W_aff,'S_affine':S_aff,
            'simple_coroots':simple_coroots,'B':B,'C':C,
            'alcove_verts_ambient':alcove_verts,'S_weyl':S_weyl,'identity':identity,'keys':keys}


def _a3_get_keyname(w):
    rw = tuple(w.reduced_word()) if hasattr(w,'reduced_word') else tuple(w)
    return 'identity' if not rw else 's'+''.join(str(x) for x in rw)


def _a3_H0_alcoves_dict(setup):
    W_fin = setup['W_finite']; alcove_verts = setup['alcove_verts_ambient']
    H0 = {}
    for w in W_fin:
        verts = [vector(_a3_to_3d(list(w.action(v).to_vector()))) for v in alcove_verts]
        H0[w] = {'vertices':verts,'color':(0.8,0.8,0.8),'reduced_word':tuple(w.reduced_word())}
    # Assign individual colors by index (matching original notebook)
    _colors = [
        (0.7,1.0,1.0),(1.0,0.7,1.0),(0.7,1.0,0.7),(1.0,1.0,0.7),
        (1.0,0.7,1.0),(0.7,1.0,0.7),(1.0,1.0,0.7),(0.7,0.8,1.0),
        (1.0,0.8,0.5),(0.8,0.7,1.0),(0.8,0.8,0.8),(1.0,0.7,0.7),
        (1.0,0.7,1.0),(0.7,1.0,0.7),(1.0,0.7,1.0),(0.7,1.0,0.7),
        (1.0,1.0,0.7),(1.0,1.0,0.7),(0.7,0.8,1.0),(0.8,0.7,1.0),
        (0.7,0.8,1.0),(1.0,0.8,0.5),(1.0,0.8,0.5),(0.8,0.7,1.0),
    ]
    for i,w in enumerate(H0.keys()):
        if i < len(_colors):
            H0[w]['color'] = _colors[i]
    return H0


def _parse_translation_input(elt_string):
    """Parse t_(n1,n2[,n3])*s_XXXX or t_(n1,n2[,n3])*e format.
    Returns (trans_coords, word_string) or None if not this format."""
    import re
    m = re.match(r't_\(([^)]+)\)\*(?:s_([0-9]*)|e)$', elt_string.strip())
    if not m: return None
    coords = [int(x.strip()) for x in m.group(1).split(',')]
    word = m.group(2) or ''
    return (coords, word)


def _a3_get_element_parts(elt_string, setup):
    W_fin = setup['W_finite']; E = setup['E']; PW0 = setup['PW0']; FW = setup['FW']
    ambient = setup['ambient']; W_aff = setup['W_affine']; S_aff = setup['S_affine']
    if not elt_string: return W_fin.one(), None, (), W_fin.one()
    # Check for t_(n1,n2,n3)*s_XXXX format
    _parsed = _parse_translation_input(elt_string)
    if _parsed is not None:
        _coords, _word = _parsed
        # Build fin_weyl from word
        if _word:
            _fin = W_fin.from_reduced_word([int(c) for c in _word])
            _S_aff2 = W_aff.simple_reflections()
            _ea = W_aff.one()
            for c in _word: _ea = _ea * _S_aff2[int(c)]
        else:
            _fin = W_fin.one()
            _ea = W_aff.one()
        # Build trans_vec from B-coords: trans_vec = sum(coords[i]*simple_coroot[i])
        _sc = setup['simple_coroots']
        _tv = ambient.zero()
        for i,n in enumerate(_coords):
            if n != 0: _tv = _tv + n*_sc[i]
        _trans_vec = None if _tv == ambient.zero() else _tv
        # Compute rw of full element: build in W_aff
        # t_(a,b,c)*w: need full affine element's reduced word
        # Use the PW0 approach with fin_weyl and translation
        _B = setup['B']
        if _trans_vec is not None:
            _t3 = vector(_a3_to_3d(list(_tv.to_vector())))
            _tB = _B.inverse()*vector(_t3)
            _tB = vector([ZZ(round(c.n(),0)) for c in _tB])
        else:
            _tB = vector([0,0,0])
        # For rw: build full element in W_aff via generators
        # This is the same as conjugacy class rw computation
        _S_FW2 = FW.simple_reflections(); _P2 = E.lattice()
        a2,b2,c2 = int(_tB[0]),int(_tB[1]),int(_tB[2])
        om1=2*a2-b2; om2=-a2+2*b2-c2; om3=-b2+2*c2
        if _word:
            _ff = _S_FW2[int(_word[0])]
            for ch in _word[1:]: _ff = _ff * _S_FW2[int(ch)]
        else: _ff = FW.one()
        _tl2 = _P2.zero()
        if om1: _tl2 = _tl2 + om1*_P2.fundamental_weight(1)
        if om2: _tl2 = _tl2 + om2*_P2.fundamental_weight(2)
        if om3: _tl2 = _tl2 + om3*_P2.fundamental_weight(3)
        _ae_fw = list(FW(PW0.from_translation(_tl2)*PW0(_ff)).to_affine_weyl_left().reduced_word())
        _ae_w = W_aff.one()
        for x in _ae_fw: _ae_w = _ae_w * S_aff[x]
        _rw = list(_ae_w.reduced_word())
        _fp = PW0(_S_FW2[int(_word[0])] if _word else FW.one())
        if _word:
            _fp_elt = _S_FW2[int(_word[0])]
            for ch in _word[1:]: _fp_elt = _fp_elt * _S_FW2[int(ch)]
        else: _fp_elt = FW.one()
        return _fin, _trans_vec, _rw, _fp_elt
    S_FW = FW.simple_reflections()
    ee = S_FW[int(elt_string[0])]
    for i in range(1,len(elt_string)): ee = ee*S_FW[int(elt_string[i])]
    w_pw0 = PW0(ee); tv, fp = w_pw0.value
    ea = W_aff.one()
    for ch in elt_string: ea = ea*S_aff[int(ch)]
    rw = ea.reduced_word()
    fin_weyl = W_fin.from_reduced_word(fp.reduced_word())
    tl = tv.value
    if tl.is_zero(): trans_vec = None
    else:
        trans_vec = ambient.zero()
        for i in tl.support():
            trans_vec = trans_vec + tl[i]*ambient.fundamental_weight(i)
    return fin_weyl, trans_vec, rw, fp


# ============================================================
# Ã₃  —  mod set, reflect, cocon, eta
# ============================================================

def _a3_mod_set(cp, crv_in_B, setup):
    VS = ZZ**3; kn = _a3_get_keyname(cp)
    mod_spans = {
        'identity':[vector([0,0,0])],
        's1':VS.span([vector([1,0,0])]), 's2':VS.span([vector([0,1,0])]),
        's3':VS.span([vector([0,0,1])]),
        's12':VS.span([vector([1,0,0]),vector([0,1,0])]),
        's21':VS.span([vector([1,0,0]),vector([0,1,0])]),
        's23':VS.span([vector([0,1,0]),vector([0,0,1])]),
        's32':VS.span([vector([0,1,0]),vector([0,0,1])]),
        's31':VS.span([vector([1,0,1]),vector([0,0,2])]),
        's121':VS.span([vector([1,1,0])]), 's232':VS.span([vector([0,1,1])]),
        's123':VS.span([vector([1,0,3]),vector([0,1,3]),vector([0,0,4])]),
        's321':VS.span([vector([1,0,3]),vector([0,1,3]),vector([0,0,4])]),
        's231':VS.span([vector([1,0,1]),vector([0,1,2]),vector([0,0,4])]),
        's312':VS.span([vector([1,0,1]),vector([0,1,2]),vector([0,0,4])]),
        's1231':VS.span([vector([1,1,0]),vector([0,0,1])]),
        's3121':VS.span([vector([1,1,0]),vector([0,0,1])]),
        's1232':VS.span([vector([1,0,0]),vector([0,1,1])]),
        's2321':VS.span([vector([1,0,0]),vector([0,1,1])]),
        's2312':VS.span([vector([1,0,-1]),vector([0,2,2])]),
        's12312':VS.span([vector([1,0,1]),vector([0,2,1]),vector([0,0,2])]),
        's23121':VS.span([vector([1,0,1]),vector([0,2,1]),vector([0,0,2])]),
        's12321':VS.span([vector([1,1,1])]),
        's123121':VS.span([vector([1,0,1]),vector([0,2,0])]),
    }
    if kn not in mod_spans: return [vector([0,0,0])],[vector([0,0,0])]
    sp = mod_spans[kn]
    mp = [vector([0,0,0])] if kn=='identity' else [v for v in crv_in_B if v in sp]
    return sp, mp


def _a3_reflect(v_B, g, setup):
    a,b,c = v_B[0],v_B[1],v_B[2]; kn = _a3_get_keyname(g)
    rd = {
        'identity':vector([a,b,c]),'s1':vector([-a+b,b,c]),'s2':vector([a,a-b+c,c]),
        's3':vector([a,b,b-c]),'s12':vector([-b+c,a-b+c,c]),'s21':vector([-a+b,-a+c,c]),
        's23':vector([a,a-c,b-c]),'s31':vector([-a+b,b,b-c]),'s32':vector([a,a-b+c,a-b]),
        's121':vector([-b+c,-a+c,c]),'s123':vector([-c,a-c,b-c]),
        's231':vector([-a+b,-a+b-c,b-c]),'s232':vector([a,a-c,a-b]),
        's312':vector([-b+c,a-b+c,a-b]),'s321':vector([-a+b,-a+c,-a]),
        's1231':vector([-c,-a+b-c,b-c]),'s1232':vector([-c,a-c,a-b]),
        's2312':vector([-b+c,-b,a-b]),'s2321':vector([-a+b,-a+b-c,-a]),
        's3121':vector([-b+c,-a+c,-a]),'s12312':vector([-c,-b,a-b]),
        's12321':vector([-c,-a+b-c,-a]),'s23121':vector([-b+c,-b,-a]),
        's123121':vector([-c,-b,-a]),
    }
    return rd.get(kn, vector([a,b,c]))


def _a3_cocon_H0(w1, w2, setup):
    return [g for g in setup['W_finite'] if
            tuple((g*w1*g.inverse()).reduced_word()) == tuple(w2.reduced_word())]


def _a3_cocon_H0_lambda(tb1, tb2, cocon, mod_sp, setup):
    lpmul=[]; cocon_lambda=[]
    for g in cocon:
        gl = _a3_reflect(tb1, g, setup); diff = tb2-gl
        inm = (diff==vector([0,0,0])) if isinstance(mod_sp,list) else (diff in mod_sp)
        if inm: lpmul.append(diff); cocon_lambda.append(g)
    return lpmul, cocon_lambda


def _a3_eta(lpmul_i, cocon_elem, cp2, crv_in_B, setup):
    l0,l1,l2 = lpmul_i[0],lpmul_i[1],lpmul_i[2]; kn = _a3_get_keyname(cp2)
    eta = []
    if kn=='identity': eta = crv_in_B if (l0==0 and l1==0 and l2==0) else []
    elif kn=='s1':  eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1 and 0==l2]
    elif kn=='s2':  eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]-v[2]==l1 and 0==l2]
    elif kn=='s3':  eta=[v for v in crv_in_B if 0==l0 and 0==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s12': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and -v[0]+2*v[1]-v[2]==l1 and 0==l2]
    elif kn=='s21': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[1]-v[2]==l1 and 0==l2]
    elif kn=='s23': eta=[v for v in crv_in_B if 0==l0 and -v[0]+v[1]+v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s32': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]-v[2]==l1 and -v[0]+v[1]+v[2]==l2]
    elif kn=='s31': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s121':eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and v[0]+v[1]-v[2]==l1 and 0==l2]
    elif kn=='s232':eta=[v for v in crv_in_B if 0==l0 and -v[0]+v[1]+v[2]==l1 and -v[0]+v[1]+v[2]==l2]
    elif kn=='s123':eta=[v for v in crv_in_B if v[0]+v[2]==l0 and -v[0]+v[1]+v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s321':eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[1]-v[2]==l1 and v[0]+v[2]==l2]
    elif kn=='s231':eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s312':eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and -v[0]+2*v[1]-v[2]==l1 and -v[0]+v[1]+v[2]==l2]
    elif kn=='s1231':eta=[v for v in crv_in_B if v[0]+v[2]==l0 and v[0]+v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s3121':eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and v[0]+v[1]-v[2]==l1 and v[0]+v[2]==l2]
    elif kn=='s1232':eta=[v for v in crv_in_B if v[0]+v[2]==l0 and -v[0]+v[1]+v[2]==l1 and -v[0]+v[1]+v[2]==l2]
    elif kn=='s2321':eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[2]==l1 and v[0]+v[2]==l2]
    elif kn=='s2312':eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and 2*v[1]==l1 and -v[0]+v[1]+v[2]==l2]
    elif kn=='s12312':eta=[v for v in crv_in_B if v[0]+v[2]==l0 and 2*v[1]==l1 and -v[0]+v[1]+v[2]==l2]
    elif kn=='s23121':eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and 2*v[1]==l1 and v[0]+v[2]==l2]
    elif kn=='s12321':eta=[v for v in crv_in_B if v[0]+v[2]==l0 and v[0]+v[2]==l1 and v[0]+v[2]==l2]
    elif kn=='s123121':eta=[v for v in crv_in_B if v[0]+v[2]==l0 and 2*v[1]==l1 and v[0]+v[2]==l2]
    return eta, [[v,cocon_elem] for v in eta]


# ============================================================
# Ã₃  —  plot and main functions
# ============================================================
TET_FACES_A3 = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]

def _a3_plot_element(G, fin_weyl, trans_vec, H0dict, setup, border_color, bbox, in_cocon=False, rw=None):
    B = setup['B']; alcove = H0dict[fin_weyl]; verts = alcove['vertices']; color = alcove['color']
    if rw is None: rw = tuple(fin_weyl.reduced_word())
    if trans_vec is not None:
        t3 = vector(_a3_to_3d(list(trans_vec.to_vector())))
        tv = [vector([float(v[0]+t3[0]),float(v[1]+t3[1]),float(v[2]+t3[2])]) for v in verts]
    else:
        tv = [vector([float(v[0]),float(v[1]),float(v[2])]) for v in verts]
    centroid = vector([sum(tv[i][j] for i in range(4))/4 for j in range(3)])
    if abs(centroid[0])>bbox or abs(centroid[1])>bbox or abs(centroid[2])>bbox: return G,False,None
    G += IndexFaceSet(TET_FACES_A3, tv, color=color, opacity=0.4)
    if in_cocon: G = _r3_add_stripes_to_tetrahedron(G, tv, TET_FACES_A3)
    for i in range(4):
        for j in range(i+1,4): G += line3d([tv[i],tv[j]],color=border_color,thickness=4)
    lbl = _r3_format_label(rw)
    G += text3d(lbl,tuple(centroid),color=border_color,fontsize=20)
    return G, True, tv


def _a3_plot_cocon(G, bbox, B, S_weyl, identity, H0dict, eta_cocon, setup, e1v=None, e2v=None):
    identity_color = (0.7,1.0,1.0); plotted = set()
    cocon_keys = set()
    for elt in eta_cocon:
        vB = elt[0]; rw = tuple(elt[1].reduced_word()) if hasattr(elt[1],'reduced_word') else tuple(elt[1])
        ve = B*vB
        ak = identity if not rw else S_weyl[rw[0]]
        if rw:
            for idx in rw[1:]: ak = ak*S_weyl[idx]
        vs = H0dict[ak]['vertices']
        tv = [vector([float(v[0]+ve[0]),float(v[1]+ve[1]),float(v[2]+ve[2])]) for v in vs]
        cocon_keys.add(tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in tv])))
    def in_cocon_check(vl):
        k = tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in vl]))
        return k in cocon_keys
    e1_in = in_cocon_check(e1v) if e1v else False
    e2_in = in_cocon_check(e2v) if e2v else False
    id_in = any(len(tuple(elt[1].reduced_word() if hasattr(elt[1],'reduced_word') else elt[1]))==0
                and all(abs(elt[0][j])<0.001 for j in range(3)) for elt in eta_cocon)
    id_v = H0dict[identity]['vertices']
    id_vf = [vector([float(v[0]),float(v[1]),float(v[2])]) for v in id_v]
    G += IndexFaceSet(TET_FACES_A3, id_vf, color=identity_color, opacity=0.4)
    if id_in: G = _r3_add_stripes_to_tetrahedron(G, id_vf, TET_FACES_A3)
    for i in range(4):
        for j in range(i+1,4): G += line3d([id_vf[i],id_vf[j]],color='black',thickness=3)
    ic = vector([sum(id_vf[i][j] for i in range(4))/4 for j in range(3)])
    G += text3d("e",tuple(ic),color='black',fontsize=20)
    plotted.add(tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in id_vf])))
    for elt in eta_cocon:
        vB = elt[0]; rw = tuple(elt[1].reduced_word()) if hasattr(elt[1],'reduced_word') else tuple(elt[1])
        ve = B*vB
        if abs(ve[0])>bbox or abs(ve[1])>bbox or abs(ve[2])>bbox: continue
        ak = identity if not rw else S_weyl[rw[0]]
        if rw:
            for idx in rw[1:]: ak = ak*S_weyl[idx]
        alcove = H0dict[ak]; vs = alcove['vertices']; color = alcove['color']
        tv = [vector([float(v[0]+ve[0]),float(v[1]+ve[1]),float(v[2]+ve[2])]) for v in vs]
        ck = tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in tv]))
        if ck in plotted: continue
        plotted.add(ck)
        centroid = vector([sum(tv[i][j] for i in range(4))/4 for j in range(3)])
        if all(abs(centroid[j])<=bbox for j in range(3)):
            G += IndexFaceSet(TET_FACES_A3, tv, color=color, opacity=0.3)
            for i in range(4):
                for j in range(i+1,4): G += line3d([tv[i],tv[j]],color=color,thickness=1)
    return G, e1_in, e2_in


def _a3_conjugacy_class_plot(W_fin, B, H0dict, bbox, alcove_verts, trans_B, fin_weyl):
    E = ExtendedAffineWeylGroup(["A",3,1]); PW0 = E.PW0(); FW = E.FW()
    S_FW = FW.simple_reflections(); P = E.lattice()
    G = Graphics(); S = W_fin.simple_reflections(); identity = W_fin.one()
    drawn = set()
    identity_verts = H0dict[identity]['vertices']
    ik = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in identity_verts]))
    drawn.add(ik)
    colored_by_finite = {}
    # If input is identity, include it in the element list
    if not fin_weyl.reduced_word() and all(c==0 for c in trans_B):
        colored_by_finite[()] = {(0,0,0): True}
    # Compute lpm once from fin_weyl's mod_plot (not each k's), then reflect by each k
    _fw_mp = H0dict[fin_weyl].get('mod_plot',[vector([0,0,0])])
    _base_lpm = [trans_B + vector(_fw_mp[i]) for i in range(len(_fw_mp))]
    for k in H0dict.keys():
        cp_tmp = k; lpm = _base_lpm
        kn = _a3_get_keyname(k)
        a,b,c = trans_B[0],trans_B[1],trans_B[2]
        reflect_map = {
            'identity':lpm,
            's1':[vector([-v[0]+v[1],v[1],v[2]]) for v in lpm],
            's2':[vector([v[0],v[0]-v[1]+v[2],v[2]]) for v in lpm],
            's3':[vector([v[0],v[1],v[1]-v[2]]) for v in lpm],
            's12':[vector([-v[1]+v[2],v[0]-v[1]+v[2],v[2]]) for v in lpm],
            's21':[vector([-v[0]+v[1],-v[0]+v[2],v[2]]) for v in lpm],
            's23':[vector([v[0],v[0]-v[2],v[1]-v[2]]) for v in lpm],
            's31':[vector([-v[0]+v[1],v[1],v[1]-v[2]]) for v in lpm],
            's32':[vector([v[0],v[0]-v[1]+v[2],v[0]-v[1]]) for v in lpm],
            's121':[vector([-v[1]+v[2],-v[0]+v[2],v[2]]) for v in lpm],
            's123':[vector([-v[2],v[0]-v[2],v[1]-v[2]]) for v in lpm],
            's231':[vector([-v[0]+v[1],-v[0]+v[1]-v[2],v[1]-v[2]]) for v in lpm],
            's232':[vector([v[0],v[0]-v[2],v[0]-v[1]]) for v in lpm],
            's312':[vector([-v[1]+v[2],v[0]-v[1]+v[2],v[0]-v[1]]) for v in lpm],
            's321':[vector([-v[0]+v[1],-v[0]+v[2],-v[0]]) for v in lpm],
            's1231':[vector([-v[2],-v[0]+v[1]-v[2],v[1]-v[2]]) for v in lpm],
            's1232':[vector([-v[2],v[0]-v[2],v[0]-v[1]]) for v in lpm],
            's2312':[vector([-v[1]+v[2],-v[1],v[0]-v[1]]) for v in lpm],
            's2321':[vector([-v[0]+v[1],-v[0]+v[1]-v[2],-v[0]]) for v in lpm],
            's3121':[vector([-v[1]+v[2],-v[0]+v[2],-v[0]]) for v in lpm],
            's12312':[vector([-v[2],-v[1],v[0]-v[1]]) for v in lpm],
            's12321':[vector([-v[2],-v[0]+v[1]-v[2],-v[0]]) for v in lpm],
            's23121':[vector([-v[1]+v[2],-v[1],-v[0]]) for v in lpm],
            's123121':[vector([-v[2],-v[1],-v[0]]) for v in lpm],
        }
        H0dict[k]['reflect'] = reflect_map.get(kn, lpm)
        conj = k * fin_weyl * k.inverse()
        H0dict[k]['finite_conjugate'] = tuple(conj.reduced_word())
    # Filter and plot
    for k in H0dict.keys():
        fc = H0dict[k]['finite_conjugate']
        if not fc: ak = identity
        else:
            ak = S[fc[0]]
            for idx in fc[1:]: ak = ak*S[idx]
        verts = H0dict[ak]['vertices']; color = H0dict[ak]['color']
        for vec in H0dict[k]['reflect']:
            check = B*vec
            if abs(check[0])>bbox or abs(check[1])>bbox or abs(check[2])>bbox: continue
            tv = [vector([float(verts[i][0]+check[0]),float(verts[i][1]+check[1]),
                          float(verts[i][2]+check[2])]) for i in range(4)]
            centroid = vector([sum(tv[i][j] for i in range(4))/4 for j in range(3)])
            if all(abs(centroid[j])<=bbox for j in range(3)):
                ck = tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in tv]))
                if ck in drawn: continue
                drawn.add(ck)
                G += IndexFaceSet(TET_FACES_A3, tv, color=color, opacity=0.2)
                for i in range(4):
                    for j in range(i+1,4): G += line3d([tv[i],tv[j]],color=color,thickness=1)
                # Track B-coords and full element for printing
                _vb = (int(round(float(vec[0]))),int(round(float(vec[1]))),int(round(float(vec[2]))))
                colored_by_finite.setdefault(fc, {})[_vb] = True
    # Print elements using only the drawn (filtered) alcoves
    _W_p = WeylGroup(["A",3,1]); _S_p = _W_p.simple_reflections()
    _E_p = ExtendedAffineWeylGroup(["A",3,1]); _PW0_p = _E_p.PW0()
    _FW_p = _E_p.FW(); _P_p = _E_p.lattice(); _S_FW_p = _FW_p.simple_reflections()
    def _fmt_b(x):
        s = str(x)
        if any(c in s for c in ["sqrt","pi","/"]): return s
        f = float(x)
        return str(int(f)) if f == int(f) else f"{f:.10g}"
    print("Coroot Basis (columns = simple coroots of A3):")
    for _i in range(3):
        print(f"  b{_i+1} = ({', '.join([_fmt_b(B[_j,_i]) for _j in range(3)])})")  
    print("Translation vectors are in basis B coordinates (n1, n2, n3)")
    def _build_rw(fc, a, b, c):
        om1=2*a-b; om2=-a+2*b-c; om3=-b+2*c
        if fc:
            _ff = _S_FW_p[fc[0]]
            for _x in fc[1:]: _ff = _ff * _S_FW_p[_x]
        else: _ff = _FW_p.one()
        _tl = _P_p.zero()
        if om1: _tl = _tl + om1*_P_p.fundamental_weight(1)
        if om2: _tl = _tl + om2*_P_p.fundamental_weight(2)
        if om3: _tl = _tl + om3*_P_p.fundamental_weight(3)
        _ae_fw = list(_FW_p(_PW0_p.from_translation(_tl)*_PW0_p(_ff)).to_affine_weyl_left().reduced_word())
        _ae_w = _W_p.one()
        for _x in _ae_fw: _ae_w = _ae_w * _S_p[_x]
        _rw = list(_ae_w.reduced_word())
        return "s_"+"".join(str(x) for x in _rw) if _rw else "e"
    _all_elems = []
    for fc, vecs in colored_by_finite.items():
        if not vecs: continue
        _fc_str = "s_"+"".join(str(x) for x in fc) if fc else "e"
        for (a,b,c) in sorted(vecs.keys()):
            _rw_str = _build_rw(fc, a, b, c)
            _all_elems.append((_rw_str, f"t_({a},{b},{c}) · {_fc_str}"))
    _all_elems.sort(key=lambda x:(len(x[0]),x[0]))
    _w = max(len(r) for r,_ in _all_elems) if _all_elems else 4
    print(f"Conjugacy class contains {len(_all_elems)} element(s) in bounding box:")
    for _rw_str, _ts in _all_elems:
        print(f"  {_rw_str:<{_w}}  =  {_ts}")
    return G


def _conjugacy_class_A3(elt, bbox):
    """Conjugacy class in Ã₃."""
    setup = _a3_setup(); B = setup['B']; W_fin = setup['W_finite']
    alcove_verts = setup['alcove_verts_ambient']
    H0dict = _a3_H0_alcoves_dict(setup)
    crv_in_B = _r3_get_coroot_vertices_in_B(bbox+1)
    fin_weyl, trans_vec, rw, fp = _a3_get_element_parts(elt, setup)
    # Set mod_plot for each element
    for k in H0dict:
        _, mp = _a3_mod_set(k, crv_in_B, setup)
        H0dict[k]['mod_plot'] = mp
    if trans_vec is not None:
        t3 = vector(_a3_to_3d(list(trans_vec.to_vector())))
        trans_B = B.inverse()*vector(t3)
        trans_B = vector([ZZ(round(c.n(),0)) for c in trans_B])
    else:
        trans_B = vector([0,0,0])
    # Print input element info
    _fc_str = "s_" + "".join(str(x) for x in fin_weyl.reduced_word()) if fin_weyl.reduced_word() else "e"
    _rw_str = "s_" + "".join(str(x) for x in rw) if rw else "e"
    _tb = trans_B
    _t_str = f"t_({_tb[0]},{_tb[1]},{_tb[2]})"
    print(f"Input element:  {_rw_str}  =  {_t_str} · {_fc_str}")
    # Wireframe + identity
    p_wire = Graphics(); p_col = Graphics(); p_lab = Graphics(); drawn_w = set()
    p_col, p_lab, drawn_w, _ = _r3_plot_colored_alcove_a3(
        p_col, p_lab, drawn_w, W_fin.one(), None, (0.7,1.0,1.0), 'black', bbox, alcove_verts)
    # If input is identity, add stripes to identity alcove (it is the chosen element)
    if not fin_weyl.reduced_word() and trans_vec is None:
        _id_verts_3d = [vector(_a3_to_3d(list(v.to_vector()))) for v in alcove_verts]
        p_col = _r3_add_stripes_to_tetrahedron(p_col, _id_verts_3d, TET_FACES_A3)
    elif elt:
        p_col, p_lab, drawn_w, _ = _r3_plot_colored_alcove_a3(
            p_col, p_lab, drawn_w, fin_weyl, trans_vec, (0.9,0.9,0.9), 'red', bbox, alcove_verts, rw)
    # Draw wireframe
    crv_euclid = []
    R_fin = setup['R_finite']; ambient = setup['ambient']
    simple_coroots = setup['simple_coroots']
    for vB in crv_in_B:
        ta = ambient.zero()
        for i,coef in enumerate(vB):
            if coef!=0: ta = ta + coef*simple_coroots[i]
        crv_euclid.append(ta)
    for t in crv_euclid:
        for w in W_fin:
            tv = [vector(_a3_to_3d(list((w.action(v)+t).to_vector()))) for v in alcove_verts]
            if not _r3_is_in_bounds(tv, bbox): continue
            ck = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in tv]))
            if ck in drawn_w: continue
            drawn_w.add(ck)
            for i in range(4):
                for j in range(i+1,4): p_wire += line3d([tv[i],tv[j]],color='#888888',thickness=1)
    G_conj = _a3_conjugacy_class_plot(W_fin, B, H0dict, bbox, alcove_verts, trans_B, fin_weyl)
    # Print conjugacy class elements using PW0.from_translation (original notebook approach)
    _E_a3 = ExtendedAffineWeylGroup(["A", 3, 1])
    _PW0_a3 = _E_a3.PW0(); _FW_a3 = _E_a3.FW(); _P_a3 = _E_a3.lattice()
    _S_FW_a3 = _FW_a3.simple_reflections()
    combined = p_wire + p_col + p_lab + G_conj
    combined += point3d((0,0,0),size=20,color='red')
    try:
        _html_dest = os.path.join(_WEYL_OUT_DIR, "result.html")
        combined.save(_html_dest, viewer="threejs", frame=True, online=True)
        _results.append(("html", _html_dest))
    except Exception as _e3d:
        print(f"3D save failed: {_e3d}")


def _r3_plot_colored_alcove_a3(p_col, p_lab, drawn, ew, trans_vec, color, border_color,
                                 bbox, alcove_verts, rw=None, label_fs=24):
    tv = []
    for va in alcove_verts:
        wv = ew.action(va)
        if trans_vec is not None: wv = wv + trans_vec
        tv.append(vector(_a3_to_3d(list(wv.to_vector()))))
    if rw is None: rw = ew.reduced_word()
    if not _r3_is_in_bounds(tv, bbox):
        print(f"Element {_r3_format_label(rw)} outside bounding box"); return p_col,p_lab,drawn,False
    ck = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in tv]))
    drawn.add(ck)
    p_col += IndexFaceSet(TET_FACES_A3, tv, color=color, opacity=0.3)
    for i in range(4):
        for j in range(i+1,4): p_col += line3d([tv[i],tv[j]],color=border_color,thickness=3)
    centroid = np.mean([list(v) for v in tv],axis=0)
    p_lab += text3d(_r3_format_label(rw),tuple(centroid),color='black',fontsize=label_fs)
    return p_col, p_lab, drawn, True


def _Cocon_Elements_A3(elt1, elt2, bbox):
    """Coconjugation in Ã₃."""
    setup = _a3_setup(); B = setup['B']; S_weyl = setup['S_weyl']; identity = setup['identity']
    H0dict = _a3_H0_alcoves_dict(setup)
    crv_in_B = _r3_get_coroot_vertices_in_B(bbox)
    fw1,tv1,rw1,fp1 = _a3_get_element_parts(elt1,setup)
    fw2,tv2,rw2,fp2 = _a3_get_element_parts(elt2,setup)
    def to_B(tv):
        if tv is None: return vector([0,0,0])
        t3 = vector(_a3_to_3d(list(tv.to_vector())))
        tb = B.inverse()*vector(t3)
        return vector([ZZ(round(c.n(),0)) for c in tb])
    tb1 = to_B(tv1); tb2 = to_B(tv2)
    mod_sp, _ = _a3_mod_set(fw2, crv_in_B, setup)
    cocon = _a3_cocon_H0(fw1, fw2, setup)
    lpmul, cocon_lambda = _a3_cocon_H0_lambda(tb1, tb2, cocon, mod_sp, setup)
    if not cocon_lambda: print("Coconjugacy set is empty!"); return
    eta_cocon = []
    for i in range(len(lpmul)):
        _,tmp = _a3_eta(lpmul[i], cocon_lambda[i], fw2, crv_in_B, setup)
        eta_cocon.extend(tmp)
    if not eta_cocon: print("No eta solutions found!"); return
    G = Graphics(); drawn_cc = set()
    # Draw wireframe for all alcoves in bounding box
    R_fin = setup['R_finite']; ambient = setup['ambient']
    simple_coroots = setup['simple_coroots']
    W_fin = setup['W_finite']
    alcove_verts = setup['alcove_verts_ambient']
    crv_euclid = []
    for vB in _r3_get_coroot_vertices_in_B(bbox+1):
        ta = ambient.zero()
        for i,coef in enumerate(vB):
            if coef!=0: ta = ta + coef*simple_coroots[i]
        crv_euclid.append(ta)
    for t in crv_euclid:
        for w in W_fin:
            tv_w = [vector(_a3_to_3d(list((w.action(v)+t).to_vector()))) for v in alcove_verts]
            if not _r3_is_in_bounds(tv_w, bbox): continue
            ck = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in tv_w]))
            if ck in drawn_cc: continue
            drawn_cc.add(ck)
            for i in range(4):
                for j in range(i+1,4): G += line3d([tv_w[i],tv_w[j]],color='#888888',thickness=1)
    def get_verts_3d(fw, tv):
        alcove = H0dict[fw]; vs = alcove['vertices']
        if tv is None: return [vector([float(v[0]),float(v[1]),float(v[2])]) for v in vs]
        t3 = vector(_a3_to_3d(list(tv.to_vector())))
        return [vector([float(v[0]+t3[0]),float(v[1]+t3[1]),float(v[2]+t3[2])]) for v in vs]
    e1v = get_verts_3d(fw1,tv1); e2v = get_verts_3d(fw2,tv2)
    G, e1_in, e2_in = _a3_plot_cocon(G,bbox,B,S_weyl,identity,H0dict,eta_cocon,setup,e1v,e2v)
    if elt1 and not (not fw1.reduced_word() and tv1 is None):
        G,_,_ = _a3_plot_element(G,fw1,tv1,H0dict,setup,'red',bbox,e1_in,rw=rw1)
    if elt2 and not (not fw2.reduced_word() and tv2 is None):
        G,_,_ = _a3_plot_element(G,fw2,tv2,H0dict,setup,'blue',bbox,e2_in,rw=rw2)
    G += point3d((0,0,0),size=20,color='black')
    # Print input elements and coconjugation set
    _W_cc = WeylGroup(["A",3,1]); _S_cc = _W_cc.simple_reflections()
    _E_cc = ExtendedAffineWeylGroup(["A",3,1]); _PW0_cc = _E_cc.PW0()
    _FW_cc = _E_cc.FW(); _P_cc = _E_cc.lattice(); _S_FW_cc = _FW_cc.simple_reflections()
    def _cc_rw(fc, tb):
        a,b,c = int(tb[0]),int(tb[1]),int(tb[2])
        om1=2*a-b; om2=-a+2*b-c; om3=-b+2*c
        if tuple(fc.reduced_word()) if hasattr(fc,"reduced_word") else fc:
            _rw_fc = tuple(fc.reduced_word()) if hasattr(fc,"reduced_word") else tuple(fc)
            _ff = _S_FW_cc[_rw_fc[0]]
            for _x in _rw_fc[1:]: _ff = _ff * _S_FW_cc[_x]
        else: _ff = _FW_cc.one()
        _tl = _P_cc.zero()
        if om1: _tl = _tl + om1*_P_cc.fundamental_weight(1)
        if om2: _tl = _tl + om2*_P_cc.fundamental_weight(2)
        if om3: _tl = _tl + om3*_P_cc.fundamental_weight(3)
        _ae = list(_FW_cc(_PW0_cc.from_translation(_tl)*_PW0_cc(_ff)).to_affine_weyl_left().reduced_word())
        _aw = _W_cc.one()
        for _x in _ae: _aw = _aw * _S_cc[_x]
        _r = list(_aw.reduced_word())
        return "s_"+"".join(str(x) for x in _r) if _r else "e"
    # Input element info
    _rw1 = "s_"+"".join(str(x) for x in rw1) if rw1 else "e"
    _rw2 = "s_"+"".join(str(x) for x in rw2) if rw2 else "e"
    _fc1 = "s_"+"".join(str(x) for x in fw1.reduced_word()) if fw1.reduced_word() else "e"
    _fc2 = "s_"+"".join(str(x) for x in fw2.reduced_word()) if fw2.reduced_word() else "e"
    print(f"Input element 1: {_rw1}  =  t_({tb1[0]},{tb1[1]},{tb1[2]}) · {_fc1}")
    print(f"Input element 2: {_rw2}  =  t_({tb2[0]},{tb2[1]},{tb2[2]}) · {_fc2}")
    # Coconjugation set elements
    _cc_elems = []
    for elt in eta_cocon:
        vB = elt[0]; _fw = elt[1]
        _rw_fc = tuple(_fw.reduced_word()) if hasattr(_fw,"reduced_word") else tuple(_fw)
        _fc_s = "s_"+"".join(str(x) for x in _rw_fc) if _rw_fc else "e"
        _rw_s = _cc_rw(_fw, vB)
        _cc_elems.append((_rw_s, f"t_({int(vB[0])},{int(vB[1])},{int(vB[2])}) · {_fc_s}"))
    _cc_elems.sort(key=lambda x:(len(x[0]),x[0]))
    _cc_elems = list(dict.fromkeys(_cc_elems))  # deduplicate
    _cw = max(len(r) for r,_ in _cc_elems) if _cc_elems else 4
    print(f"Coconjugation set contains {len(_cc_elems)} element(s) in bounding box:")
    for _rs,_ts in _cc_elems:
        print(f"  {_rs:<{_cw}}  =  {_ts}")
    try:
        _html_dest = os.path.join(_WEYL_OUT_DIR, "result.html")
        G.save(_html_dest, viewer="threejs", frame=True, axes=True, aspect_ratio=1, online=True)
        _results.append(("html", _html_dest))
    except Exception as _e3d:
        print(f"3D save error: {_e3d}")


# ============================================================
# B̃₃ and C̃₃  —  shared framework (parameterised by letter)
# ============================================================

def _bc3_setup(letter):
    """Setup for letter = 'B' or 'C'."""
    R_fin = RootSystem([letter,3]); ambient = R_fin.ambient_space()
    W_fin = ambient.weyl_group(); S_weyl = W_fin.simple_reflections(); identity = W_fin.one()
    E = ExtendedAffineWeylGroup([letter,3,1])
    PW0 = E.PW0(); WF = E.WF(); P = E.lattice()
    W_aff = WeylGroup([letter,3,1]); S_aff = W_aff.simple_reflections()
    simple_coroots = [ambient.simple_coroot(i) for i in [1,2,3]]
    coroots_3d = [vector(list(cr.to_vector())[:3]) for cr in simple_coroots]
    B = matrix(QQ, coroots_3d).transpose()
    return {'R_finite':R_fin,'ambient':ambient,'W_finite':W_fin,'S_weyl':S_weyl,
            'identity':identity,'E':E,'PW0':PW0,'WF':WF,'P':P,
            'W_affine':W_aff,'S_affine':S_aff,'B':B,'letter':letter}


def _bc3_get_keyname(w):
    rw = tuple(w.reduced_word()) if hasattr(w,'reduced_word') else tuple(w)
    return 'identity' if not rw else 's'+''.join(str(x) for x in rw)


def _bc3_get_element_parts(elt_string, setup):
    W_fin = setup['W_finite']; PW0 = setup['PW0']; WF = setup['WF']
    ambient = setup['ambient']; W_aff = setup['W_affine']; S_aff = setup['S_affine']
    if not elt_string: return W_fin.one(), None, (), W_fin.one()
    # Check for t_(n1,n2,n3)*s_XXXX format
    _parsed = _parse_translation_input(elt_string)
    if _parsed is not None:
        _coords, _word = _parsed
        _fin = W_fin.from_reduced_word([int(c) for c in _word]) if _word else W_fin.one()
        letter = setup['letter']
        _sc = [ambient.simple_coroot(i) for i in [1,2,3]]
        _tv = ambient.zero()
        for i,n in enumerate(_coords):
            if n != 0: _tv = _tv + n*_sc[i]
        _trans_vec = None if _tv == ambient.zero() else _tv
        _B = setup['B']
        if _trans_vec is not None:
            _t3 = vector(list(_tv.to_vector())[:3])
            _tB = _B.inverse()*_t3
            _tB = vector([ZZ(round(c,0)) for c in _tB])
        else:
            _tB = vector([0,0,0])
        _S_WF2 = WF.simple_reflections()
        a2,b2,c2 = int(_tB[0]),int(_tB[1]),int(_tB[2])
        if letter == 'B': om1=2*a2-b2; om2=-a2+2*b2-2*c2; om3=-b2+2*c2
        else:             om1=2*a2-b2; om2=-a2+2*b2-c2;   om3=-2*b2+2*c2
        _ff = WF.one()
        if _word:
            _ff = _S_WF2[int(_word[0])]
            for ch in _word[1:]: _ff = _ff * _S_WF2[int(ch)]
        E2 = setup['E']; P2 = E2.lattice(); PW02 = E2.PW0()
        _tl2 = P2.zero()
        if om1: _tl2 = _tl2 + om1*P2.fundamental_weight(1)
        if om2: _tl2 = _tl2 + om2*P2.fundamental_weight(2)
        if om3: _tl2 = _tl2 + om3*P2.fundamental_weight(3)
        _ae_fw = list(WF(PW02.from_translation(_tl2)*PW02(_ff)).to_affine_weyl_left().reduced_word())
        _ae_w = W_aff.one()
        for x in _ae_fw: _ae_w = _ae_w * S_aff[x]
        _rw = list(_ae_w.reduced_word())
        return _fin, _trans_vec, _rw, _ff
    S_WF = WF.simple_reflections()
    ee = S_WF[int(elt_string[0])]
    for i in range(1,len(elt_string)): ee = ee*S_WF[int(elt_string[i])]
    w_pw0 = PW0(ee); tv, fp = w_pw0.value
    ea = W_aff.one()
    for ch in elt_string: ea = ea*S_aff[int(ch)]
    rw = ea.reduced_word(); fin_weyl = W_fin.from_reduced_word(fp.reduced_word())
    tl = tv.value
    if tl.is_zero(): trans_vec = None
    else:
        cw1 = ambient.from_vector(vector([1,0,0]))
        cw2 = ambient.from_vector(vector([1,1,0]))
        cw3 = ambient.from_vector(vector([1,1,1]))
        cws = {1:cw1,2:cw2,3:cw3}
        trans_vec = ambient.zero()
        for i in tl.support(): trans_vec = trans_vec + tl[i]*cws[i]
    return fin_weyl, trans_vec, rw, fp


def _bc3_H0_alcoves_dict(setup):
    letter = setup['letter']; W_fin = setup['W_finite']; ambient = setup['ambient']
    if letter == 'B':
        origin = vector(QQ,[0,0,0]); omega3 = vector(QQ,[1/2,1/2,1/2])
        omega2 = vector(QQ,[1,1,0]); omega1 = vector(QQ,[1,0,0])
        fund_verts = [origin, omega3, omega2/2, omega1/2]
    else:  # C
        origin = vector(QQ,[0,0,0]); omega1 = vector(QQ,[1,0,0])
        omega2 = vector(QQ,[1,1,0]); omega3 = vector(QQ,[1,1,1])
        fund_verts = [origin, omega1/2, omega2/2, omega3/2]
    identity_color = (0.7,1.0,1.0)
    light_magenta = (1.0,0.7,1.0); light_lime = (0.7,1.0,0.7); light_yellow = (1.0,1.0,0.7)
    light_cb = (0.7,0.8,1.0); light_do = (1.0,0.8,0.5); light_dor = (0.8,0.7,1.0)
    light_gray = (0.8,0.8,0.8); light_red = (1.0,0.7,0.7); light_pink = (1.0,0.85,0.9)
    # color by conjugacy class (same for B3 and C3 since W(B3)=W(C3))
    cbw = {
        # Class 1: Identity
        ():(0.7,1.0,1.0),
        # Class 2: size 3
        (3,):light_magenta,(2,3,2):light_lime,(1,2,3,2,1):light_yellow,
        # Class 3: size 6
        (2,):light_magenta,(3,1,2,3,1):light_lime,(1,):light_yellow,
        (1,2,1):light_cb,(2,3,1,2,3,1,2):light_do,(3,2,3):light_dor,
        # Class 4: size 8
        (2,1):light_magenta,(2,3,1,2,3,2):light_lime,(3,1,2,3,1,2):light_yellow,
        (3,1,2,3):light_cb,(1,2):light_do,(2,3,1,2,3,1):light_dor,
        (3,2,3,1,2,1):light_gray,(3,2,3,1):light_red,
        # Class 5: size 6
        (3,2):light_magenta,(2,3):light_lime,(2,3,2,1):light_yellow,
        (3,1,2,1):light_cb,(1,2,3,2):light_do,(1,2,3,1):light_dor,
        # Class 6: size 8
        (3,2,1):light_magenta,(2,3,1):light_lime,(2,3,1,2,1):light_yellow,
        (3,2,3,1,2,3,2):light_cb,(3,1,2):light_do,(1,2,3):light_dor,
        (1,2,3,1,2):light_gray,(3,2,3,1,2,3,1):light_red,
        # Class 7: size 6
        (2,3,1,2):light_magenta,(3,2,3,1,2,3,2,1):light_lime,(3,1):light_yellow,
        (1,2,3,1,2,1):light_cb,(3,2,3,1,2,3,1,2):light_do,(3,2,3,1,2,3):light_dor,
        # Class 8: size 6
        (2,3,1,2,3,2,1):light_magenta,(3,1,2,3,1,2,1):light_lime,
        (3,1,2,3,2):light_yellow,(2,3,1,2,3):light_cb,
        (3,2,3,2,1):light_do,(3,2,3,1,2):light_dor,
        # Class 9: size 3
        (2,3,1,2,3,1,2,1):light_magenta,(3,1,2,3,2,1):light_lime,(3,2,3,2):light_yellow,
        # Class 10: central element
        (3,2,3,1,2,3,1,2,1):light_magenta,
    }
    H0 = {}
    for w in W_fin:
        rw = tuple(w.reduced_word())
        tv = []
        for v in fund_verts:
            av = ambient.from_vector(vector(QQ,list(v)))
            tv.append(vector(QQ,list(w.action(av).to_vector())[:3]))
        color = cbw.get(rw,(0.5,0.5,0.5))
        H0[w] = {'vertices':tv,'color':color,'reduced_word':rw}
    return H0


# ============================================================
# B̃₃ / C̃₃  —  mod set, reflect, cocon, eta (parameterised)
# ============================================================

def _bc3_mod_set(cp, crv_in_B, setup):
    letter = setup['letter']; VS = ZZ**3; kn = _bc3_get_keyname(cp)
    # B3 and C3 have distinct mod_spans but same structure
    if letter == 'B':
        mod_spans = {
            'identity':[vector([0,0,0])],
            's1':VS.span([vector([1,0,0])]),'s2':VS.span([vector([0,1,0])]),
            's3':VS.span([vector([0,0,1])]),'s12':VS.span([vector([1,0,0]),vector([0,1,0])]),
            's21':VS.span([vector([1,0,0]),vector([0,1,0])]),'s23':VS.span([vector([0,1,0]),vector([0,0,1])]),
            's32':VS.span([vector([0,1,0]),vector([0,0,1])]),'s31':VS.span([vector([1,0,1]),vector([0,0,2])]),
            's121':VS.span([vector([1,1,0])]),'s123':VS.span([vector([1,0,1]),vector([0,1,1]),vector([0,0,2])]),
            's232':VS.span([vector([0,2,1])]),'s231':VS.span([vector([1,0,1]),vector([0,1,0]),vector([0,0,2])]),
            's321':VS.span([vector([1,0,1]),vector([0,1,1]),vector([0,0,2])]),'s312':VS.span([vector([1,0,1]),vector([0,1,0]),vector([0,0,2])]),
            's323':VS.span([vector([0,1,1])]),'s1231':VS.span([vector([1,1,0]),vector([0,0,1])]),
            's1232':VS.span([vector([1,0,0]),vector([0,2,1])]),'s2312':VS.span([vector([1,3,1]),vector([0,4,2])]),
            's2321':VS.span([vector([1,0,0]),vector([0,2,1])]),'s3121':VS.span([vector([1,1,0]),vector([0,0,1])]),
            's3123':VS.span([vector([1,0,0]),vector([0,1,1])]),'s3231':VS.span([vector([1,0,0]),vector([0,1,1])]),
            's3232':VS.span([vector([0,2,0]),vector([0,0,1])]),'s12312':VS.span([vector([1,0,0]),vector([0,1,1]),vector([0,0,2])]),
            's12321':VS.span([vector([2,2,1])]),'s23121':VS.span([vector([1,0,0]),vector([0,1,1]),vector([0,0,2])]),
            's23123':VS.span([vector([1,3,0]),vector([0,4,0]),vector([0,0,1])]),'s31231':VS.span([vector([1,1,1])]),
            's31232':VS.span([vector([1,0,1]),vector([0,2,1]),vector([0,0,2])]),'s32312':VS.span([vector([1,3,0]),vector([0,4,0]),vector([0,0,1])]),
            's32321':VS.span([vector([1,0,1]),vector([0,2,1]),vector([0,0,2])]),'s123121':VS.span([vector([2,1,1]),vector([0,2,0])]),
            's231231':VS.span([vector([1,0,1]),vector([0,1,0])]),'s231232':VS.span([vector([1,0,-1]),vector([0,1,1])]),
            's312312':VS.span([vector([1,0,1]),vector([0,1,0])]),'s312321':VS.span([vector([2,2,0]),vector([0,0,1])]),
            's323121':VS.span([vector([1,0,-1]),vector([0,1,1])]),'s323123':VS.span([vector([1,3,2]),vector([0,4,2])]),
            's3123121':VS.span([vector([2,1,0]),vector([0,2,0]),vector([0,0,1])]),'s3231231':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),
            's3231232':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),'s2312312':VS.span([vector([1,2,1])]),
            's2312321':VS.span([vector([2,1,0]),vector([0,2,0]),vector([0,0,1])]),'s23123121':VS.span([vector([2,0,0]),vector([0,2,1])]),
            's32312312':VS.span([vector([1,2,0]),vector([0,0,2])]),'s32312321':VS.span([vector([2,1,0]),vector([0,2,2])]),
            's323123121':VS.span([vector([2,0,0]),vector([0,2,0]),vector([0,0,2])]),
        }
    else:  # C
        mod_spans = {
            'identity':[vector([0,0,0])],
            's1':VS.span([vector([1,0,0])]),'s2':VS.span([vector([0,1,0])]),
            's3':VS.span([vector([0,0,2])]),'s12':VS.span([vector([1,0,0]),vector([0,1,0])]),
            's21':VS.span([vector([1,0,0]),vector([0,1,0])]),'s23':VS.span([vector([0,1,0]),vector([0,0,2])]),
            's32':VS.span([vector([0,1,0]),vector([0,0,2])]),'s31':VS.span([vector([1,0,0]),vector([0,0,2])]),
            's121':VS.span([vector([1,1,0])]),'s123':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),
            's232':VS.span([vector([0,2,2])]),'s231':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),
            's321':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),'s312':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),
            's323':VS.span([vector([0,1,2])]),'s1231':VS.span([vector([1,1,0]),vector([0,0,2])]),
            's1232':VS.span([vector([1,0,0]),vector([0,2,2])]),'s2312':VS.span([vector([1,1,0]),vector([0,2,2])]),
            's2321':VS.span([vector([1,0,0]),vector([0,2,2])]),'s3121':VS.span([vector([1,1,0]),vector([0,0,2])]),
            's3123':VS.span([vector([1,0,0]),vector([0,1,2])]),'s3231':VS.span([vector([1,0,0]),vector([0,1,2])]),
            's3232':VS.span([vector([0,2,0]),vector([0,0,2])]),'s12312':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),
            's12321':VS.span([vector([2,2,2])]),'s23121':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),
            's23123':VS.span([vector([1,1,0]),vector([0,2,0]),vector([0,0,2])]),'s31231':VS.span([vector([1,1,2])]),
            's31232':VS.span([vector([1,0,0]),vector([0,2,0]),vector([0,0,2])]),'s32312':VS.span([vector([1,1,0]),vector([0,2,0]),vector([0,0,2])]),
            's32321':VS.span([vector([1,0,0]),vector([0,2,0]),vector([0,0,2])]),'s123121':VS.span([vector([2,0,2]),vector([0,1,0])]),
            's231231':VS.span([vector([1,0,2]),vector([0,1,0])]),'s231232':VS.span([vector([1,0,-2]),vector([0,1,2])]),
            's312312':VS.span([vector([1,0,2]),vector([0,1,0])]),'s312321':VS.span([vector([2,2,0]),vector([0,0,2])]),
            's323121':VS.span([vector([1,0,-2]),vector([0,1,2])]),'s323123':VS.span([vector([1,1,2]),vector([0,2,2])]),
            's3123121':VS.span([vector([2,0,0]),vector([0,1,0]),vector([0,0,2])]),'s3231231':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),
            's3231232':VS.span([vector([1,0,0]),vector([0,1,0]),vector([0,0,2])]),'s2312312':VS.span([vector([1,2,2])]),
            's2312321':VS.span([vector([2,0,0]),vector([0,1,0]),vector([0,0,2])]),'s23123121':VS.span([vector([2,0,0]),vector([0,2,2])]),
            's32312312':VS.span([vector([1,2,0]),vector([0,0,2])]),'s32312321':VS.span([vector([2,0,-2]),vector([0,1,2])]),
            's323123121':VS.span([vector([2,0,0]),vector([0,2,0]),vector([0,0,2])]),
        }
    if kn not in mod_spans: return [vector([0,0,0])],[vector([0,0,0])]
    sp = mod_spans[kn]
    mp = [vector([0,0,0])] if kn=='identity' else [v for v in crv_in_B if v in sp]
    return sp, mp


# ── Reflection formulas for B3 and C3 ─────────────────────────────────────

def _b3_reflect(vB, g, setup):
    a,b,c = vB[0],vB[1],vB[2]; kn = _bc3_get_keyname(g)
    rd = {
        'identity':vector([a,b,c]),'s1':vector([-a+b,b,c]),'s2':vector([a,a-b+2*c,c]),
        's3':vector([a,b,b-c]),'s12':vector([-b+2*c,a-b+2*c,c]),'s21':vector([-a+b,-a+2*c,c]),
        's23':vector([a,a+b-2*c,b-c]),'s31':vector([-a+b,b,b-c]),'s32':vector([a,a-b+2*c,a-b+c]),
        's121':vector([-b+2*c,-a+2*c,c]),'s123':vector([b-2*c,a+b-2*c,b-c]),
        's231':vector([-a+b,-a+2*b-2*c,b-c]),'s232':vector([a,2*a-b,a-b+c]),
        's312':vector([-b+2*c,a-b+2*c,a-b+c]),'s321':vector([-a+b,-a+2*c,-a+c]),
        's323':vector([a,a+b-2*c,a-c]),'s1231':vector([b-2*c,-a+2*b-2*c,b-c]),
        's1232':vector([a-b,2*a-b,a-b+c]),'s2312':vector([-b+2*c,a-2*b+2*c,a-b+c]),
        's2321':vector([-a+b,-2*a+b,-a+c]),'s3121':vector([-b+2*c,-a+2*c,-a+c]),
        's3123':vector([b-2*c,a+b-2*c,a-c]),'s3231':vector([-a+b,-a+2*b-2*c,-a+b-c]),
        's3232':vector([a,2*a-b,a-c]),'s12312':vector([a-b,a-2*b+2*c,a-b+c]),
        's12321':vector([-a,-2*a+b,-a+c]),'s23121':vector([-b+2*c,-a-b+2*c,-a+c]),
        's23123':vector([b-2*c,a-2*c,a-c]),'s31231':vector([b-2*c,-a+2*b-2*c,-a+b-c]),
        's31232':vector([a-b,2*a-b,a-c]),'s32312':vector([-b+2*c,a-2*b+2*c,-b+c]),
        's32321':vector([-a+b,-2*a+b,-a+b-c]),'s123121':vector([-a,-a-b+2*c,-a+c]),
        's231231':vector([b-2*c,-a+b-2*c,-a+b-c]),'s231232':vector([a-b,a-2*c,a-c]),
        's312312':vector([a-b,a-2*b+2*c,-b+c]),'s312321':vector([-a,-2*a+b,-a+b-c]),
        's323121':vector([-b+2*c,-a-b+2*c,-b+c]),'s323123':vector([b-2*c,a-2*c,-c]),
        's2312312':vector([a-b,-b,-b+c]),'s2312321':vector([-a,-a+b-2*c,-a+b-c]),
        's3123121':vector([-a,-a-b+2*c,-b+c]),'s3231231':vector([b-2*c,-a+b-2*c,-c]),
        's3231232':vector([a-b,a-2*c,-c]),'s23123121':vector([-a,-b,-b+c]),
        's32312312':vector([a-b,-b,-c]),'s32312321':vector([-a,-a+b-2*c,-c]),
        's323123121':vector([-a,-b,-c]),
    }
    return rd.get(kn, vector([a,b,c]))


def _c3_reflect(vB, g, setup):
    a,b,c = vB[0],vB[1],vB[2]; kn = _bc3_get_keyname(g)
    rd = {
        'identity':vector([a,b,c]),'s1':vector([-a+b,b,c]),'s2':vector([a,a-b+c,c]),
        's3':vector([a,b,2*b-c]),'s12':vector([-b+c,a-b+c,c]),'s21':vector([-a+b,-a+c,c]),
        's23':vector([a,a+b-c,2*b-c]),'s31':vector([-a+b,b,2*b-c]),'s32':vector([a,a-b+c,2*a-2*b+c]),
        's121':vector([-b+c,-a+c,c]),'s123':vector([b-c,a+b-c,2*b-c]),
        's231':vector([-a+b,-a+2*b-c,2*b-c]),'s232':vector([a,2*a-b,2*a-2*b+c]),
        's312':vector([-b+c,a-b+c,2*a-2*b+c]),'s321':vector([-a+b,-a+c,-2*a+c]),
        's323':vector([a,a+b-c,2*a-c]),'s1231':vector([b-c,-a+2*b-c,2*b-c]),
        's1232':vector([a-b,2*a-b,2*a-2*b+c]),'s2312':vector([-b+c,a-2*b+c,2*a-2*b+c]),
        's2321':vector([-a+b,-2*a+b,-2*a+c]),'s3121':vector([-b+c,-a+c,-2*a+c]),
        's3123':vector([b-c,a+b-c,2*a-c]),'s3231':vector([-a+b,-a+2*b-c,-2*a+2*b-c]),
        's3232':vector([a,2*a-b,2*a-c]),'s12312':vector([a-b,a-2*b+c,2*a-2*b+c]),
        's12321':vector([-a,-2*a+b,-2*a+c]),'s23121':vector([-b+c,-a-b+c,-2*a+c]),
        's23123':vector([b-c,a-c,2*a-c]),'s31231':vector([b-c,-a+2*b-c,-2*a+2*b-c]),
        's31232':vector([a-b,2*a-b,2*a-c]),'s32312':vector([-b+c,a-2*b+c,-2*b+c]),
        's32321':vector([-a+b,-2*a+b,-2*a+2*b-c]),'s123121':vector([-a,-a-b+c,-2*a+c]),
        's231231':vector([b-c,-a+b-c,-2*a+2*b-c]),'s231232':vector([a-b,a-c,2*a-c]),
        's312312':vector([a-b,a-2*b+c,-2*b+c]),'s312321':vector([-a,-2*a+b,-2*a+2*b-c]),
        's323121':vector([-b+c,-a-b+c,-2*b+c]),'s323123':vector([b-c,a-c,-c]),
        's2312312':vector([a-b,-b,-2*b+c]),'s2312321':vector([-a,-a+b-c,-2*a+2*b-c]),
        's3123121':vector([-a,-a-b+c,-2*b+c]),'s3231231':vector([b-c,-a+b-c,-c]),
        's3231232':vector([a-b,a-c,-c]),'s23123121':vector([-a,-b,-2*b+c]),
        's32312312':vector([a-b,-b,-c]),'s32312321':vector([-a,-a+b-c,-c]),
        's323123121':vector([-a,-b,-c]),
    }
    return rd.get(kn, vector([a,b,c]))


def _bc3_reflect(vB, g, setup):
    if setup['letter']=='B': return _b3_reflect(vB,g,setup)
    else: return _c3_reflect(vB,g,setup)


def _bc3_cocon_H0(w1, w2, setup):
    return [g for g in setup['W_finite'] if g*w1*g.inverse() == w2]


def _bc3_cocon_H0_lambda(tb1, tb2, cocon, mod_sp, setup):
    lpmul=[]; cocon_lambda=[]
    for g in cocon:
        gl = _bc3_reflect(tb1,g,setup); diff = tb2-gl
        inm = (diff==vector([0,0,0])) if isinstance(mod_sp,list) else (diff in mod_sp)
        if inm: lpmul.append(diff); cocon_lambda.append(g)
    return lpmul, cocon_lambda


# ── Eta equations for B3 and C3 ────────────────────────────────────────────

def _b3_eta(lpmul_i, cocon_elem, cp2, crv_in_B, setup):
    l0,l1,l2 = lpmul_i; kn = _bc3_get_keyname(cp2); eta=[]
    if kn=='identity': eta=crv_in_B if l0==0 and l1==0 and l2==0 else []
    elif kn=='s1': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1 and 0==l2]
    elif kn=='s2': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]-2*v[2]==l1 and 0==l2]
    elif kn=='s3': eta=[v for v in crv_in_B if 0==l0 and 0==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s12': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and -v[0]+2*v[1]-2*v[2]==l1 and 0==l2]
    elif kn=='s21': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[1]-2*v[2]==l1 and 0==l2]
    elif kn=='s23': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s31': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s32': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]-2*v[2]==l1 and -v[0]+v[1]==l2]
    elif kn=='s121': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and v[0]+v[1]-2*v[2]==l1 and 0==l2]
    elif kn=='s123': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and -v[0]+2*v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s231': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]-v[1]+2*v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s232': eta=[v for v in crv_in_B if 0==l0 and -2*v[0]+2*v[1]==l1 and -v[0]+v[1]==l2]
    elif kn=='s312': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and -v[0]+2*v[1]-2*v[2]==l1 and -v[0]+v[1]==l2]
    elif kn=='s321': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[1]-2*v[2]==l1 and v[0]==l2]
    elif kn=='s323': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[2]==l1 and -v[0]+2*v[2]==l2]
    elif kn=='s1231': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and v[0]-v[1]+2*v[2]==l1 and -v[1]+2*v[2]==l2]
    elif kn=='s1232': eta=[v for v in crv_in_B if v[1]==l0 and -2*v[0]+2*v[1]==l1 and -v[0]+v[1]==l2]
    elif kn=='s2312': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and -v[0]+3*v[1]-2*v[2]==l1 and -v[0]+v[1]==l2]
    elif kn=='s2321': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 2*v[0]==l1 and v[0]==l2]
    elif kn=='s3121': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and v[0]+v[1]-2*v[2]==l1 and v[0]==l2]
    elif kn=='s3123': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and -v[0]+2*v[2]==l1 and -v[0]+2*v[2]==l2]
    elif kn=='s3231': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]-v[1]+2*v[2]==l1 and v[0]-v[1]+2*v[2]==l2]
    elif kn=='s3232': eta=[v for v in crv_in_B if 0==l0 and -2*v[0]+2*v[1]==l1 and -v[0]+2*v[2]==l2]
    elif kn=='s12312': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+3*v[1]-2*v[2]==l1 and -v[0]+v[1]==l2]
    elif kn=='s12321': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[0]==l1 and v[0]==l2]
    elif kn=='s23121': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and v[0]+2*v[1]-2*v[2]==l1 and v[0]==l2]
    elif kn=='s23123': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and -v[0]+v[1]+2*v[2]==l1 and -v[0]+2*v[2]==l2]
    elif kn=='s31231': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and v[0]-v[1]+2*v[2]==l1 and v[0]-v[1]+2*v[2]==l2]
    elif kn=='s31232': eta=[v for v in crv_in_B if v[1]==l0 and -2*v[0]+2*v[1]==l1 and -v[0]+2*v[2]==l2]
    elif kn=='s32312': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and -v[0]+3*v[1]-2*v[2]==l1 and v[1]==l2]
    elif kn=='s32321': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 2*v[0]==l1 and v[0]-v[1]+2*v[2]==l2]
    elif kn=='s123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+2*v[1]-2*v[2]==l1 and v[0]==l2]
    elif kn=='s231231': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and v[0]+2*v[2]==l1 and v[0]-v[1]+2*v[2]==l2]
    elif kn=='s231232': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+v[1]+2*v[2]==l1 and -v[0]+2*v[2]==l2]
    elif kn=='s312312': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+3*v[1]-2*v[2]==l1 and v[1]==l2]
    elif kn=='s312321': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[0]==l1 and v[0]-v[1]+2*v[2]==l2]
    elif kn=='s323121': eta=[v for v in crv_in_B if v[0]+v[1]-2*v[2]==l0 and v[0]+2*v[1]-2*v[2]==l1 and v[1]==l2]
    elif kn=='s323123': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and -v[0]+v[1]+2*v[2]==l1 and 2*v[2]==l2]
    elif kn=='s2312312': eta=[v for v in crv_in_B if v[1]==l0 and 2*v[1]==l1 and v[1]==l2]
    elif kn=='s2312321': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+2*v[2]==l1 and v[0]-v[1]+2*v[2]==l2]
    elif kn=='s3123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+2*v[1]-2*v[2]==l1 and v[1]==l2]
    elif kn=='s3231231': eta=[v for v in crv_in_B if v[0]-v[1]+2*v[2]==l0 and v[0]+2*v[2]==l1 and 2*v[2]==l2]
    elif kn=='s3231232': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+v[1]+2*v[2]==l1 and 2*v[2]==l2]
    elif kn=='s23123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[1]==l1 and v[1]==l2]
    elif kn=='s32312312': eta=[v for v in crv_in_B if v[1]==l0 and 2*v[1]==l1 and 2*v[2]==l2]
    elif kn=='s32312321': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+2*v[2]==l1 and 2*v[2]==l2]
    elif kn=='s323123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[1]==l1 and 2*v[2]==l2]
    return eta, [[v,cocon_elem] for v in eta]


def _c3_eta(lpmul_i, cocon_elem, cp2, crv_in_B, setup):
    l0,l1,l2 = lpmul_i; kn = _bc3_get_keyname(cp2); eta=[]
    if kn=='identity': eta=crv_in_B if l0==0 and l1==0 and l2==0 else []
    elif kn=='s1': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1 and 0==l2]
    elif kn=='s2': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]-v[2]==l1 and 0==l2]
    elif kn=='s3': eta=[v for v in crv_in_B if 0==l0 and 0==l1 and -2*v[1]+2*v[2]==l2]
    elif kn=='s12': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and -v[0]+2*v[1]-v[2]==l1 and 0==l2]
    elif kn=='s21': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[1]-v[2]==l1 and 0==l2]
    elif kn=='s23': eta=[v for v in crv_in_B if 0==l0 and -v[0]+v[2]==l1 and -2*v[1]+2*v[2]==l2]
    elif kn=='s31': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 0==l1 and -2*v[1]+2*v[2]==l2]
    elif kn=='s32': eta=[v for v in crv_in_B if 0==l0 and -v[0]+2*v[1]-v[2]==l1 and -2*v[0]+2*v[1]==l2]
    elif kn=='s121': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and v[0]+v[1]-v[2]==l1 and 0==l2]
    elif kn=='s123': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and -v[0]+v[2]==l1 and -2*v[1]+2*v[2]==l2]
    elif kn=='s231': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]-v[1]+v[2]==l1 and -2*v[1]+2*v[2]==l2]
    elif kn=='s232': eta=[v for v in crv_in_B if 0==l0 and -2*v[0]+2*v[1]==l1 and -2*v[0]+2*v[1]==l2]
    elif kn=='s312': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and -v[0]+2*v[1]-v[2]==l1 and -2*v[0]+2*v[1]==l2]
    elif kn=='s321': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]+v[1]-v[2]==l1 and 2*v[0]==l2]
    elif kn=='s323': eta=[v for v in crv_in_B if 0==l0 and -v[0]+v[2]==l1 and -2*v[0]+2*v[2]==l2]
    elif kn=='s1231': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and v[0]-v[1]+v[2]==l1 and -2*v[1]+2*v[2]==l2]
    elif kn=='s1232': eta=[v for v in crv_in_B if v[1]==l0 and -2*v[0]+2*v[1]==l1 and -2*v[0]+2*v[1]==l2]
    elif kn=='s2312': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and -v[0]+3*v[1]-v[2]==l1 and -2*v[0]+2*v[1]==l2]
    elif kn=='s2321': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 2*v[0]==l1 and 2*v[0]==l2]
    elif kn=='s3121': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and v[0]+v[1]-v[2]==l1 and 2*v[0]==l2]
    elif kn=='s3123': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and -v[0]+v[2]==l1 and -2*v[0]+2*v[2]==l2]
    elif kn=='s3231': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and v[0]-v[1]+v[2]==l1 and 2*v[0]-2*v[1]+2*v[2]==l2]
    elif kn=='s3232': eta=[v for v in crv_in_B if 0==l0 and -2*v[0]+2*v[1]==l1 and -2*v[0]+2*v[2]==l2]
    elif kn=='s12312': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+3*v[1]-v[2]==l1 and -2*v[0]+2*v[1]==l2]
    elif kn=='s12321': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[0]==l1 and 2*v[0]==l2]
    elif kn=='s23121': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and v[0]+2*v[1]-v[2]==l1 and 2*v[0]==l2]
    elif kn=='s23123': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and -v[0]+v[1]+v[2]==l1 and -2*v[0]+2*v[2]==l2]
    elif kn=='s31231': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and v[0]-v[1]+v[2]==l1 and 2*v[0]-2*v[1]+2*v[2]==l2]
    elif kn=='s31232': eta=[v for v in crv_in_B if v[1]==l0 and -2*v[0]+2*v[1]==l1 and -2*v[0]+2*v[2]==l2]
    elif kn=='s32312': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and -v[0]+3*v[1]-v[2]==l1 and 2*v[1]==l2]
    elif kn=='s32321': eta=[v for v in crv_in_B if 2*v[0]-v[1]==l0 and 2*v[0]==l1 and 2*v[0]-2*v[1]+2*v[2]==l2]
    elif kn=='s123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+2*v[1]-v[2]==l1 and 2*v[0]==l2]
    elif kn=='s231231': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and v[0]+v[2]==l1 and 2*v[0]-2*v[1]+2*v[2]==l2]
    elif kn=='s231232': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+v[1]+v[2]==l1 and -2*v[0]+2*v[2]==l2]
    elif kn=='s312312': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+3*v[1]-v[2]==l1 and 2*v[1]==l2]
    elif kn=='s312321': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[0]==l1 and 2*v[0]-2*v[1]+2*v[2]==l2]
    elif kn=='s323121': eta=[v for v in crv_in_B if v[0]+v[1]-v[2]==l0 and v[0]+2*v[1]-v[2]==l1 and 2*v[1]==l2]
    elif kn=='s323123': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and -v[0]+v[1]+v[2]==l1 and 2*v[2]==l2]
    elif kn=='s2312312': eta=[v for v in crv_in_B if v[1]==l0 and 2*v[1]==l1 and 2*v[1]==l2]
    elif kn=='s2312321': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+v[2]==l1 and 2*v[0]-2*v[1]+2*v[2]==l2]
    elif kn=='s3123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+2*v[1]-v[2]==l1 and 2*v[1]==l2]
    elif kn=='s3231231': eta=[v for v in crv_in_B if v[0]-v[1]+v[2]==l0 and v[0]+v[2]==l1 and 2*v[2]==l2]
    elif kn=='s3231232': eta=[v for v in crv_in_B if v[1]==l0 and -v[0]+v[1]+v[2]==l1 and 2*v[2]==l2]
    elif kn=='s23123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[1]==l1 and 2*v[1]==l2]
    elif kn=='s32312312': eta=[v for v in crv_in_B if v[1]==l0 and 2*v[1]==l1 and 2*v[2]==l2]
    elif kn=='s32312321': eta=[v for v in crv_in_B if 2*v[0]==l0 and v[0]+v[2]==l1 and 2*v[2]==l2]
    elif kn=='s323123121': eta=[v for v in crv_in_B if 2*v[0]==l0 and 2*v[1]==l1 and 2*v[2]==l2]
    return eta, [[v,cocon_elem] for v in eta]


def _bc3_eta(lpmul_i, cocon_elem, cp2, crv_in_B, setup):
    if setup['letter']=='B': return _b3_eta(lpmul_i,cocon_elem,cp2,crv_in_B,setup)
    else: return _c3_eta(lpmul_i,cocon_elem,cp2,crv_in_B,setup)


# ── B3/C3 plot and main functions ──────────────────────────────────────────

TET_FACES_BC3 = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]

def _bc3_plot_element(G, fin_weyl, trans_vec, H0dict, setup, border_color, bbox, in_cocon=False, rw=None):
    alcove = H0dict[fin_weyl]; verts = alcove['vertices']; color = alcove['color']
    if rw is None: rw = tuple(fin_weyl.reduced_word())
    if trans_vec is not None:
        t3 = vector(list(trans_vec.to_vector())[:3])
        tv = [vector([float(v[0]+t3[0]),float(v[1]+t3[1]),float(v[2]+t3[2])]) for v in verts]
    else:
        tv = [vector([float(v[0]),float(v[1]),float(v[2])]) for v in verts]
    centroid = vector([sum(tv[i][j] for i in range(4))/4 for j in range(3)])
    if abs(centroid[0])>bbox or abs(centroid[1])>bbox or abs(centroid[2])>bbox: return G,False,None
    G += IndexFaceSet(TET_FACES_BC3, tv, color=color, opacity=0.4)
    if in_cocon: G = _r3_add_stripes_to_tetrahedron(G,tv,TET_FACES_BC3)
    for i in range(4):
        for j in range(i+1,4): G += line3d([tv[i],tv[j]],color=border_color,thickness=4)
    lbl = _r3_format_label(rw)
    G += text3d(lbl,tuple(centroid),color=border_color,fontsize=20)
    return G,True,tv


def _bc3_plot_cocon(G, bbox, B, S_weyl, identity, H0dict, eta_cocon, setup, e1v=None, e2v=None):
    identity_color = (0.7,1.0,1.0); plotted = set()
    cocon_keys = set()
    for elt in eta_cocon:
        vB = elt[0]; rw = tuple(elt[1].reduced_word()) if hasattr(elt[1],'reduced_word') else tuple(elt[1])
        ve = B*vB
        ak = identity if not rw else S_weyl[rw[0]]
        if rw:
            for idx in rw[1:]: ak = ak*S_weyl[idx]
        vs = H0dict[ak]['vertices']
        tv = [vector([float(v[0]+ve[0]),float(v[1]+ve[1]),float(v[2]+ve[2])]) for v in vs]
        cocon_keys.add(tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in tv])))
    def in_cocon(vl):
        return tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in vl])) in cocon_keys
    e1_in = in_cocon(e1v) if e1v else False
    e2_in = in_cocon(e2v) if e2v else False
    id_in = any(len(tuple(elt[1].reduced_word() if hasattr(elt[1],'reduced_word') else elt[1]))==0
                and all(abs(elt[0][j])<0.001 for j in range(3)) for elt in eta_cocon)
    id_vf = [vector([float(v[0]),float(v[1]),float(v[2])]) for v in H0dict[identity]['vertices']]
    G += IndexFaceSet(TET_FACES_BC3, id_vf, color=identity_color, opacity=0.4)
    if id_in: G = _r3_add_stripes_to_tetrahedron(G, id_vf, TET_FACES_BC3)
    for i in range(4):
        for j in range(i+1,4): G += line3d([id_vf[i],id_vf[j]],color='black',thickness=3)
    ic = vector([sum(id_vf[i][j] for i in range(4))/4 for j in range(3)])
    G += text3d("e",tuple(ic),color='black',fontsize=20)
    plotted.add(tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in id_vf])))
    for elt in eta_cocon:
        vB = elt[0]; rw = tuple(elt[1].reduced_word() if hasattr(elt[1],'reduced_word') else elt[1])
        ve = B*vB
        if abs(ve[0])>bbox or abs(ve[1])>bbox or abs(ve[2])>bbox: continue
        ak = identity if not rw else S_weyl[rw[0]]
        if rw:
            for idx in rw[1:]: ak = ak*S_weyl[idx]
        alcove = H0dict[ak]; vs = alcove['vertices']; color = alcove['color']
        tv = [vector([float(v[0]+ve[0]),float(v[1]+ve[1]),float(v[2]+ve[2])]) for v in vs]
        ck = tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in tv]))
        if ck in plotted: continue
        plotted.add(ck)
        centroid = vector([sum(tv[i][j] for i in range(4))/4 for j in range(3)])
        if all(abs(centroid[j])<=bbox for j in range(3)):
            G += IndexFaceSet(TET_FACES_BC3, tv, color=color, opacity=0.3)
            for i in range(4):
                for j in range(i+1,4): G += line3d([tv[i],tv[j]],color=color,thickness=1)
    return G, e1_in, e2_in


def _bc3_conjugacy_class_plot(letter, W_fin, B, H0dict, bbox, fin_weyl=None):
    E = ExtendedAffineWeylGroup([letter,3,1]); PW0 = E.PW0(); FW = E.FW()
    S_FW = FW.simple_reflections(); P = E.lattice()
    G = Graphics(); S = W_fin.simple_reflections(); identity = W_fin.one()
    drawn = set()
    id_verts = H0dict[identity]['vertices']
    drawn.add(tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in id_verts])))
    colored_by_finite = {}
    _fw_ref = fin_weyl if fin_weyl is not None else identity
    # If input is identity, include it in the element list
    if not _fw_ref.reduced_word():
        colored_by_finite[()] = {(0,0,0): True}
    # Use fin_weyl's mod_plot for all reflections (matches original notebook)
    _fw_mp_bc = H0dict[_fw_ref].get('mod_plot',[vector([0,0,0])])
    _trans_B_bc = H0dict[_fw_ref].get('trans_B', vector([0,0,0]))
    _base_lpm_bc = [_trans_B_bc + vector(_fw_mp_bc[i]) for i in range(len(_fw_mp_bc))]
    for k in H0dict.keys():
        kn = _bc3_get_keyname(k)
        trans_B = _trans_B_bc
        lpm = _base_lpm_bc
        refl_fn = _b3_reflect if letter=='B' else _c3_reflect
        setup_tmp = {'letter':letter}
        H0dict[k]['reflect'] = [refl_fn(v, k, setup_tmp) for v in lpm]
        conj = k * H0dict[k].get('_fin_weyl', identity) * k.inverse()
        H0dict[k]['finite_conjugate'] = tuple(conj.reduced_word())
    for k in H0dict.keys():
        fc = H0dict[k]['finite_conjugate']
        if not fc: ak = identity
        else:
            ak = S[fc[0]]
            for idx in fc[1:]: ak = ak*S[idx]
        verts = H0dict[ak]['vertices']; color = H0dict[ak]['color']
        for vec in H0dict[k]['reflect']:
            check = B*vec
            if abs(check[0])>bbox or abs(check[1])>bbox or abs(check[2])>bbox: continue
            tv = [vector([float(verts[i][0]+check[0]),float(verts[i][1]+check[1]),
                          float(verts[i][2]+check[2])]) for i in range(4)]
            centroid = vector([sum(tv[i][j] for i in range(4))/4 for j in range(3)])
            if all(abs(centroid[j])<=bbox for j in range(3)):
                ck = tuple(sorted([tuple(round(v[j],6) for j in range(3)) for v in tv]))
                if ck in drawn: continue
                drawn.add(ck)
                G += IndexFaceSet(TET_FACES_BC3, tv, color=color, opacity=0.2)
                for i in range(4):
                    for j in range(i+1,4): G += line3d([tv[i],tv[j]],color=color,thickness=1)
                _vb = (int(round(float(vec[0]))),int(round(float(vec[1]))),int(round(float(vec[2]))))
                colored_by_finite.setdefault(fc, {})[_vb] = True
    # Print elements using only drawn alcoves
    _E_p = ExtendedAffineWeylGroup([letter,3,1]); _PW0_p = _E_p.PW0()
    _FW_p = _E_p.FW(); _P_p = _E_p.lattice(); _S_FW_p = _FW_p.simple_reflections()
    _W_aff_p = WeylGroup([letter,3,1]); _S_aff_p = _W_aff_p.simple_reflections()
    def _fmt_b(x):
        s = str(x)
        if any(c in s for c in ["sqrt","pi","/"]): return s
        f = float(x)
        return str(int(f)) if f == int(f) else f"{f:.10g}"
    print(f"Coroot Basis (columns = simple coroots of {letter}3):")
    for _i in range(3):
        print(f"  b{_i+1} = ({', '.join([_fmt_b(B[_j,_i]) for _j in range(3)])})")
    print("Translation vectors are in basis B coordinates (n1, n2, n3)")
    _W_p = WeylGroup(["A",3,1]); _S_p = _W_p.simple_reflections()
    _E_p = ExtendedAffineWeylGroup(["A",3,1]); _PW0_p = _E_p.PW0()
    _FW_p = _E_p.FW(); _P_p = _E_p.lattice(); _S_FW_p = _FW_p.simple_reflections()
    def _build_rw(fc, a, b, c):
        if letter == "B": om1=2*a-b; om2=-a+2*b-2*c; om3=-b+2*c
        else:             om1=2*a-b; om2=-a+2*b-c;   om3=-2*b+2*c
        if fc:
            _ff = _S_FW_p[fc[0]]
            for _x in fc[1:]: _ff = _ff * _S_FW_p[_x]
        else: _ff = _FW_p.one()
        _tl = _P_p.zero()
        if om1: _tl = _tl + om1*_P_p.fundamental_weight(1)
        if om2: _tl = _tl + om2*_P_p.fundamental_weight(2)
        if om3: _tl = _tl + om3*_P_p.fundamental_weight(3)
        _ae_fw = list(_FW_p(_PW0_p.from_translation(_tl)*_PW0_p(_ff)).to_affine_weyl_left().reduced_word())
        _ae_w = _W_p.one()
        for _x in _ae_fw: _ae_w = _ae_w * _S_p[_x]
        _rw = list(_ae_w.reduced_word())
        return "s_"+"".join(str(x) for x in _rw) if _rw else "e"
    _all_elems = []
    for fc, vecs in colored_by_finite.items():
        if not vecs: continue
        _fc_str = "s_"+"".join(str(x) for x in fc) if fc else "e"
        for (a,b,c) in sorted(vecs.keys()):
            _rw_str = _build_rw(fc, a, b, c)
            _all_elems.append((_rw_str, f"t_({a},{b},{c}) · {_fc_str}"))
    _all_elems.sort(key=lambda x:(len(x[0]),x[0]))
    _w = max(len(r) for r,_ in _all_elems) if _all_elems else 4
    print(f"Conjugacy class contains {len(_all_elems)} element(s) in bounding box:")
    for _rw_str, _ts in _all_elems:
        print(f"  {_rw_str:<{_w}}  =  {_ts}")
    return G


def _conjugacy_class_BC3(letter, elt, bbox):
    """Conjugacy class in B̃₃ or C̃₃."""
    setup = _bc3_setup(letter); B = setup['B']; W_fin = setup['W_finite']
    H0dict = _bc3_H0_alcoves_dict(setup)
    crv_in_B = _r3_get_coroot_vertices_in_B(bbox+1)
    fin_weyl, trans_vec, rw, fp = _bc3_get_element_parts(elt, setup)
    for k in H0dict:
        _, mp = _bc3_mod_set(k, crv_in_B, setup)
        H0dict[k]['mod_plot'] = mp
    if trans_vec is not None:
        t3 = vector(list(trans_vec.to_vector())[:3])
        trans_B = B.inverse()*t3; trans_B = vector([ZZ(round(c,0)) for c in trans_B])
    else:
        trans_B = vector([0,0,0])
    # Print input element info
    _fc_str_bc = "s_" + "".join(str(x) for x in fin_weyl.reduced_word()) if fin_weyl.reduced_word() else "e"
    _rw_str_bc = "s_" + "".join(str(x) for x in rw) if rw else "e"
    _t_str_bc = f"t_({trans_B[0]},{trans_B[1]},{trans_B[2]})"
    print(f"Input element:  {_rw_str_bc}  =  {_t_str_bc} · {_fc_str_bc}")
    for k in H0dict:
        H0dict[k]['trans_B'] = trans_B
        H0dict[k]['_fin_weyl'] = fin_weyl
    # Build wireframe
    ambient = setup['ambient']; simple_coroots_list = [ambient.simple_coroot(i) for i in [1,2,3]]
    fund_weights_list = [ambient.fundamental_weight(i) for i in [1,2,3]]
    origin_amb = ambient.zero()
    if letter=='B':
        alcove_verts_amb = [origin_amb, ambient.from_vector(vector([1/2,1/2,1/2])),
                            ambient.from_vector(vector([1/2,1/2,0])), ambient.from_vector(vector([1/2,0,0]))]
    else:
        alcove_verts_amb = [origin_amb, ambient.from_vector(vector([1/2,0,0])),
                            ambient.from_vector(vector([1/2,1/2,0])), ambient.from_vector(vector([1/2,1/2,1/2]))]
    p_wire = Graphics(); p_col = Graphics(); p_lab = Graphics(); drawn_w = set()
    id_tv = [vector(list(W_fin.one().action(v).to_vector())[:3]) for v in alcove_verts_amb]
    ik = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in id_tv]))
    drawn_w.add(ik)
    p_col += IndexFaceSet(TET_FACES_BC3, id_tv, color=(0.7,1.0,1.0), opacity=0.3)
    for i in range(4):
        for j in range(i+1,4): p_col += line3d([id_tv[i],id_tv[j]],color='black',thickness=2)
    centroid = np.mean([list(v) for v in id_tv],axis=0)
    p_lab += text3d("e",tuple(centroid),color='black',fontsize=24)
    # If input is identity, add stripes to identity alcove (it is the chosen element)
    if not fin_weyl.reduced_word() and trans_vec is None:
        p_col = _r3_add_stripes_to_tetrahedron(p_col, id_tv, TET_FACES_BC3)
    if elt and not (not fin_weyl.reduced_word() and trans_vec is None):
        chosen_tv = []
        for va in alcove_verts_amb:
            wv = fin_weyl.action(va)
            if trans_vec is not None: wv = wv + trans_vec
            chosen_tv.append(vector(list(wv.to_vector())[:3]))
        if _r3_is_in_bounds(chosen_tv, bbox):
            ck = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in chosen_tv]))
            drawn_w.add(ck)
            p_col += IndexFaceSet(TET_FACES_BC3, chosen_tv, color=(0.9,0.9,0.9), opacity=0.3)
            for i in range(4):
                for j in range(i+1,4): p_col += line3d([chosen_tv[i],chosen_tv[j]],color='red',thickness=4)
            centroid = np.mean([list(v) for v in chosen_tv],axis=0)
            p_lab += text3d(_r3_format_label(rw),tuple(centroid),color='red',fontsize=24)
    for tB in crv_in_B:
        ta = ambient.zero()
        for i,coef in enumerate(tB):
            if coef!=0: ta = ta + coef*simple_coroots_list[i]
        for w in W_fin:
            tv = [vector(list((w.action(v)+ta).to_vector())[:3]) for v in alcove_verts_amb]
            if not _r3_is_in_bounds(tv, bbox): continue
            ck = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in tv]))
            if ck in drawn_w: continue
            drawn_w.add(ck)
            for i in range(4):
                for j in range(i+1,4): p_wire += line3d([tv[i],tv[j]],color='#888888',thickness=1)
    G_conj = _bc3_conjugacy_class_plot(letter, W_fin, B, H0dict, bbox, fin_weyl)
    # Print conjugacy class elements using PW0.from_translation (original notebook approach)
    _E_bc3 = ExtendedAffineWeylGroup([letter, 3, 1])
    _PW0_bc3 = _E_bc3.PW0(); _FW_bc3 = _E_bc3.FW(); _P_bc3 = _E_bc3.lattice()
    _S_FW_bc3 = _FW_bc3.simple_reflections()
    combined = p_wire + p_col + p_lab + G_conj
    combined += point3d((0,0,0),size=20,color='red')
    try:
        _html_dest = os.path.join(_WEYL_OUT_DIR, "result.html")
        combined.save(_html_dest, viewer="threejs", frame=True, online=True)
        _results.append(("html", _html_dest))
    except Exception as _e3d:
        print(f"3D save failed: {_e3d}")


def _Cocon_Elements_BC3(letter, elt1, elt2, bbox):
    """Coconjugation in B̃₃ or C̃₃."""
    setup = _bc3_setup(letter); B = setup['B']
    S_weyl = setup['S_weyl']; identity = setup['identity']
    H0dict = _bc3_H0_alcoves_dict(setup)
    crv_in_B = _r3_get_coroot_vertices_in_B(bbox)
    fw1,tv1,rw1,fp1 = _bc3_get_element_parts(elt1,setup)
    fw2,tv2,rw2,fp2 = _bc3_get_element_parts(elt2,setup)
    def to_B(tv):
        if tv is None: return vector([0,0,0])
        t3 = vector(list(tv.to_vector())[:3])
        tb = B.inverse()*t3; return vector([ZZ(round(c,0)) for c in tb])
    tb1 = to_B(tv1); tb2 = to_B(tv2)
    mod_sp,_ = _bc3_mod_set(fw2,crv_in_B,setup)
    cocon = _bc3_cocon_H0(fw1,fw2,setup)
    lpmul,cocon_lambda = _bc3_cocon_H0_lambda(tb1,tb2,cocon,mod_sp,setup)
    if not cocon_lambda: print("Coconjugacy set is empty!"); return
    eta_cocon = []
    for i in range(len(lpmul)):
        _,tmp = _bc3_eta(lpmul[i],cocon_lambda[i],fw2,crv_in_B,setup)
        eta_cocon.extend(tmp)
    if not eta_cocon: print("No eta solutions found!"); return
    G = Graphics(); drawn_cc = set()
    # Draw wireframe for all alcoves in bounding box
    W_fin = setup['W_finite']; ambient = setup['ambient']
    simple_coroots_list = [ambient.simple_coroot(i) for i in [1,2,3]]
    alcove_verts_amb = [ambient.from_vector(vector(QQ, list(v) if len(list(v))==3 else list(v)[:3])) for v in H0dict[identity]['vertices']]
    for tB in _r3_get_coroot_vertices_in_B(bbox+1):
        ta = ambient.zero()
        for i,coef in enumerate(tB):
            if coef!=0: ta = ta + coef*simple_coroots_list[i]
        for w in W_fin:
            tv_w = [vector(list((w.action(v)+ta).to_vector())[:3]) for v in alcove_verts_amb]
            if not _r3_is_in_bounds(tv_w, bbox): continue
            ck = tuple(sorted([tuple(round(float(v[j]),6) for j in range(3)) for v in tv_w]))
            if ck in drawn_cc: continue
            drawn_cc.add(ck)
            for i in range(4):
                for j in range(i+1,4): G += line3d([tv_w[i],tv_w[j]],color='#888888',thickness=1)
    def get_verts(fw,tv):
        vs = H0dict[fw]['vertices']
        if tv is None: return [vector([float(v[0]),float(v[1]),float(v[2])]) for v in vs]
        t3 = vector(list(tv.to_vector())[:3])
        return [vector([float(v[0]+t3[0]),float(v[1]+t3[1]),float(v[2]+t3[2])]) for v in vs]
    e1v = get_verts(fw1,tv1); e2v = get_verts(fw2,tv2)
    G,e1_in,e2_in = _bc3_plot_cocon(G,bbox,B,S_weyl,identity,H0dict,eta_cocon,setup,e1v,e2v)
    if elt1 and not (not fw1.reduced_word() and tv1 is None):
        G,_,_ = _bc3_plot_element(G,fw1,tv1,H0dict,setup,'red',bbox,e1_in,rw=rw1)
    if elt2 and not (not fw2.reduced_word() and tv2 is None):
        G,_,_ = _bc3_plot_element(G,fw2,tv2,H0dict,setup,'blue',bbox,e2_in,rw=rw2)
    G += point3d((0,0,0),size=20,color='black')
    # Print input elements and coconjugation set
    _W_cc = WeylGroup([letter,3,1]); _S_cc = _W_cc.simple_reflections()
    _E_cc = ExtendedAffineWeylGroup([letter,3,1]); _PW0_cc = _E_cc.PW0()
    _FW_cc = _E_cc.FW(); _P_cc = _E_cc.lattice(); _S_FW_cc = _FW_cc.simple_reflections()
    def _cc_rw_bc(fc, tb):
        a,b,c = int(tb[0]),int(tb[1]),int(tb[2])
        if letter == "B": om1=2*a-b; om2=-a+2*b-2*c; om3=-b+2*c
        else:             om1=2*a-b; om2=-a+2*b-c;   om3=-2*b+2*c
        _rw_fc = tuple(fc.reduced_word()) if hasattr(fc,"reduced_word") else tuple(fc)
        if _rw_fc:
            _ff = _S_FW_cc[_rw_fc[0]]
            for _x in _rw_fc[1:]: _ff = _ff * _S_FW_cc[_x]
        else: _ff = _FW_cc.one()
        _tl = _P_cc.zero()
        if om1: _tl = _tl + om1*_P_cc.fundamental_weight(1)
        if om2: _tl = _tl + om2*_P_cc.fundamental_weight(2)
        if om3: _tl = _tl + om3*_P_cc.fundamental_weight(3)
        _ae = list(_FW_cc(_PW0_cc.from_translation(_tl)*_PW0_cc(_ff)).to_affine_weyl_left().reduced_word())
        _aw = _W_cc.one()
        for _x in _ae: _aw = _aw * _S_cc[_x]
        _r = list(_aw.reduced_word())
        return "s_"+"".join(str(x) for x in _r) if _r else "e"
    _rw1 = "s_"+"".join(str(x) for x in rw1) if rw1 else "e"
    _rw2 = "s_"+"".join(str(x) for x in rw2) if rw2 else "e"
    _fc1 = "s_"+"".join(str(x) for x in fw1.reduced_word()) if fw1.reduced_word() else "e"
    _fc2 = "s_"+"".join(str(x) for x in fw2.reduced_word()) if fw2.reduced_word() else "e"
    print(f"Input element 1: {_rw1}  =  t_({tb1[0]},{tb1[1]},{tb1[2]}) · {_fc1}")
    print(f"Input element 2: {_rw2}  =  t_({tb2[0]},{tb2[1]},{tb2[2]}) · {_fc2}")
    _cc_elems = []
    for elt in eta_cocon:
        vB = elt[0]; _fw = elt[1]
        _rw_fc = tuple(_fw.reduced_word()) if hasattr(_fw,"reduced_word") else tuple(_fw)
        _fc_s = "s_"+"".join(str(x) for x in _rw_fc) if _rw_fc else "e"
        _rw_s = _cc_rw_bc(_fw, vB)
        _cc_elems.append((_rw_s, f"t_({int(vB[0])},{int(vB[1])},{int(vB[2])}) · {_fc_s}"))
    _cc_elems.sort(key=lambda x:(len(x[0]),x[0]))
    _cc_elems = list(dict.fromkeys(_cc_elems))
    _cw = max(len(r) for r,_ in _cc_elems) if _cc_elems else 4
    print(f"Coconjugation set contains {len(_cc_elems)} element(s) in bounding box:")
    for _rs,_ts in _cc_elems:
        print(f"  {_rs:<{_cw}}  =  {_ts}")
    try:
        _html_dest = os.path.join(_WEYL_OUT_DIR, "result.html")
        G.save(_html_dest, viewer="threejs", frame=True, axes=True, aspect_ratio=1, online=True)
        _results.append(("html", _html_dest))
    except Exception as _e3d:
        print(f"3D save error: {_e3d}")


# ============================================================
# APP  —  unified entry functions
# ============================================================
# All functions below call the helpers defined above.

def A1xA1_conjugacy_class(elt_1, elt_2, bbox=10):
    """Conjugacy class in Ã₁ × Ã₁.
    elt_1, elt_2 : reduced-word strings for the two Ã₁ components, e.g. '01', '010', ''
    bbox         : bounding box (default 10)
    Example: A1xA1_conjugacy_class('01', '01')
    """
    print(f"Computing conjugacy class in Ã₁ × Ã₁  for  ({elt_1 or 'e'}, {elt_2 or 'e'})  bbox={bbox}")
    print(f"  Bounding box: {bbox}")
    _Cox_Element_A1xA1(elt_1, elt_2, bbox)
    print("  Plot complete.")


def A1xA1_coconjugation(e1, e2, e3, e4, bbox=10):
    """Coconjugation in Ã₁ × Ã₁.
    Element 1 = (e1, e2),  Element 2 = (e3, e4)
    Example: A1xA1_coconjugation('01','01','10','01')
    """
    print(f"Computing coconjugation in Ã₁ × Ã₁")
    print(f"  Element 1 = ({e1 or 'e'}, {e2 or 'e'})   Element 2 = ({e3 or 'e'}, {e4 or 'e'})   bbox={bbox}")
    _Cocon_Elements_A1xA1(e1, e2, e3, e4, bbox)


def affine2_conjugacy_class(letter, elt, bbox=10):
    """Conjugacy class in a dimension-2 affine Weyl group."""
    letter = letter.upper()
    print(f"Computing conjugacy class in Ã-{letter}₂  for  s_{elt or 'e'}  bbox={bbox}")
    _conjugacy_class_affine2(letter, elt, bbox)


def affine2_coconjugation(letter, elt1, elt2, bbox=7):
    """Coconjugation in a dimension-2 affine Weyl group."""
    letter = letter.upper()
    print(f"Computing coconjugation in Ã-{letter}₂")
    print(f"  Element 1 = s_{elt1 or 'e'}   Element 2 = s_{elt2 or 'e'}   bbox={bbox}")
    _Cocon_Elements_affine2(letter, elt1, elt2, bbox)
    print("  Plot complete.")


def _reduce_word(elt, group_type):
    """Reduce a word string in the given affine Weyl group. Returns reduced word string."""
    if not elt: return elt
    W = RootSystem(group_type).ambient_space().weyl_group()
    S = W.simple_reflections()
    e = W.one()
    for c in elt: e = e * S[int(c)]
    return ''.join(str(x) for x in e.reduced_word())

def A3_conjugacy_class(elt, bbox=3):
    """Conjugacy class in Ã₃."""
    print(f"Computing conjugacy class in Ã₃  for  s_{elt or 'e'}  bbox={bbox}")
    print(f"  |W(A₃)| = 24  |  Bounding box: {bbox}")
    _conjugacy_class_A3(elt, bbox)
    print("  Plot complete.")


def A3_coconjugation(elt1, elt2, bbox=5):
    """Coconjugation in Ã₃."""
    print(f"Computing coconjugation in Ã₃")
    print(f"  Element 1 = s_{elt1 or 'e'}   Element 2 = s_{elt2 or 'e'}   bbox={bbox}")
    _Cocon_Elements_A3(elt1, elt2, bbox)
    print("  Plot complete.")


def B3_conjugacy_class(elt, bbox=3):
    """Conjugacy class in B̃₃."""
    print(f"Computing conjugacy class in B̃₃  for  s_{elt or 'e'}  bbox={bbox}")
    print(f"  |W(B₃)| = 48  |  Bounding box: {bbox}")
    _conjugacy_class_BC3('B', elt, bbox)
    print("  Plot complete.")


def B3_coconjugation(elt1, elt2, bbox=3):
    """Coconjugation in B̃₃."""
    print(f"Computing coconjugation in B̃₃")
    print(f"  Element 1 = s_{elt1 or 'e'}   Element 2 = s_{elt2 or 'e'}   bbox={bbox}")
    _Cocon_Elements_BC3('B', elt1, elt2, bbox)


def C3_conjugacy_class(elt, bbox=3):
    """Conjugacy class in C̃₃."""
    print(f"Computing conjugacy class in C̃₃  for  s_{elt or 'e'}  bbox={bbox}")
    print(f"  |W(C₃)| = 48  |  Bounding box: {bbox}")
    _conjugacy_class_BC3('C', elt, bbox)
    print("  Plot complete.")


def C3_coconjugation(elt1, elt2, bbox=3):
    """Coconjugation in C̃₃."""
    print(f"Computing coconjugation in C̃₃")
    print(f"  Element 1 = s_{elt1 or 'e'}   Element 2 = s_{elt2 or 'e'}   bbox={bbox}")
    _Cocon_Elements_BC3('C', elt1, elt2, bbox)

# numpy is available in SageMath but needs explicit import in .sage
import numpy as np
from sage.plot.plot3d.index_face_set import IndexFaceSet

# ── Dispatch ────────────────────────────────────────────────
def _dispatch(group, comp, inputs, bbox):
    inputs = [('' if v.strip() == 'e' else v.strip()) for v in inputs]
    if group == "A1xA1":
        if comp == "conj":
            A1xA1_conjugacy_class(inputs[0], inputs[1], bbox)
        else:
            A1xA1_coconjugation(inputs[0], inputs[1], inputs[2], inputs[3], bbox)
    elif group.startswith("aff2_"):
        letter = group[5:]
        if comp == "conj":
            affine2_conjugacy_class(letter, inputs[0], bbox)
        else:
            affine2_coconjugation(letter, inputs[0], inputs[1], bbox)
    elif group == "A3":
        A3_conjugacy_class(inputs[0], bbox) if comp=="conj" else A3_coconjugation(inputs[0], inputs[1], bbox)
    elif group == "B3":
        B3_conjugacy_class(inputs[0], bbox) if comp=="conj" else B3_coconjugation(inputs[0], inputs[1], bbox)
    elif group == "C3":
        C3_conjugacy_class(inputs[0], bbox) if comp=="conj" else C3_coconjugation(inputs[0], inputs[1], bbox)
    else:
        raise ValueError(f"Unknown group: {group!r}")

# ── Plot capture ────────────────────────────────────────────
@contextlib.contextmanager
def _capture_plot(out_dir):
    from sage.plot.graphics import Graphics
    _results = []
    _orig_2d = Graphics.show

    def _show_2d(self, **kw):
        _dest = os.path.join(out_dir, "result.png")
        import traceback as _tb
        try:
            self.save(_dest, figsize=10, axes=False, frame=False)
        except Exception as _e1:
            print(f"save attempt 1 failed: {_e1}")
            try:
                self.save(_dest, figsize=10)
            except Exception as _e2:
                print(f"save attempt 2 failed: {_e2}")
                print(_tb.format_exc())
        if os.path.exists(_dest):
            _results.append(("png", _dest))
        else:
            print(f"WARNING: no PNG produced at {_dest}")

    import webbrowser as _wb
    _orig_open     = _wb.open
    _orig_open_new = _wb.open_new
    _orig_open_tab = _wb.open_new_tab

    def _wb_intercept(url, *a, **kw):
        path = url[7:] if url.startswith("file://") else url
        if path.endswith(".html") and os.path.exists(path):
            dest = os.path.join(out_dir, "result.html")
            shutil.copy(path, dest)
            _results.append(("html", dest))
        else:
            _orig_open(url, *a, **kw)
        return True

    Graphics.show      = _show_2d
    _wb.open           = _wb_intercept
    _wb.open_new       = _wb_intercept
    _wb.open_new_tab   = _wb_intercept
    try:
        yield _results
    finally:
        Graphics.show    = _orig_2d
        _wb.open         = _orig_open
        _wb.open_new     = _orig_open_new
        _wb.open_new_tab = _orig_open_tab

# ── Run ─────────────────────────────────────────────────────
_console_file = open(os.path.join(_out_dir, "console.txt"), "w")
_error_msg   = ""
_result_type = ""
_result_path = ""

# Override print globally so all computation output goes to console file
import builtins as _builtins
_orig_print = _builtins.print
def _capturing_print(*args, **kwargs):
    kwargs["file"] = _console_file
    _orig_print(*args, **kwargs)
    _console_file.flush()
_builtins.print = _capturing_print

# Set global output dir so dimension-3 functions can save 3D plots directly
_WEYL_OUT_DIR = _out_dir
_results = []

try:
    with _capture_plot(_out_dir) as _results:
        _dispatch(_group, _comp, _inputs, _bbox)
    # Always write a summary regardless of internal prints
    _elem_parts = [f"s_{v}" for v in _inputs if v]
    _elem_str = " / ".join(_elem_parts) if _elem_parts else "e (identity)"
    _group_names = {
        "A1xA1": "Ã₁ × Ã₁",
        "aff2_A": "Affine Ã₂", "aff2_B": "Affine B̃₂",
        "aff2_C": "Affine C̃₂", "aff2_G": "Affine G̃₂",
        "A3": "Affine Ã₃", "B3": "Affine B̃₃", "C3": "Affine C̃₃",
    }
    _comp_names = {"conj": "Conjugacy Class", "cocon": "Coconjugation Set"}
    print("=" * 50)
    print(f"Group:       {_group_names.get(_group, _group)}")
    print(f"Computation: {_comp_names.get(_comp, _comp)}")
    print(f"Element(s):  {_elem_str}")
    print(f"Bounding box: {_bbox}")
    print(f"Result type: {_result_type if _results else 'none'}")
    print("=" * 50)
    if _results:
        _result_type, _result_path = _results[0]
except Exception as _exc:
    import traceback as _tb
    _error_msg = _tb.format_exc()
    print(f"ERROR: {_error_msg}")

# Restore print and close console file
_builtins.print = _orig_print
_console_file.close()

_result_json = {
    "type":  _result_type,
    "path":  _result_path,
    "error": _error_msg
}
with open(os.path.join(_out_dir, "result.json"), "w") as _rf:
    _json.dump(_result_json, _rf)

print("WEYL_DONE")
sys.stdout.flush()
