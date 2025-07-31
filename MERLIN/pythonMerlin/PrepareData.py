def PrepareData(Node, Panel, Supp, Load, BarCM, RotSpring, kpf, kpb, Abar, FoldKe):
    Bend = findbend(Panel, Node)
    Fold, Bdry, Trigl = findfdbd(Panel, Bend)
    Bars = np.vstack([Bend[:, :2], Fold[:, :2], Bdry])
    B, L = dirc3d(Node, Bars)

    if Supp.shape[0] == 0:
        rs = np.array([], dtype=int)
    else:
        all_dofs = np.vstack([
            Supp[:, 0]*3 - 2, Supp[:, 0]*3 - 1, Supp[:, 0]*3
        ]).T.flatten()
        values = Supp[:, 1:4].flatten()
        rs = all_dofs[values != 0]

    if np.isscalar(Abar):
        Abar = Abar * np.ones(Bars.shape[0])

    pf0 = np.array([FoldKe(Node, Fold[i], kpf, 0) for i in range(Fold.shape[0])])
    pb0 = np.array([FoldKe(Node, Bend[i], kpb, 0) for i in range(Bend.shape[0])])

    m = Node.shape[0]
    F = np.zeros(3 * m)
    indp = Load[:, 0].astype(int)
    F[3 * indp - 3] = Load[:, 1]
    F[3 * indp - 2] = Load[:, 2]
    F[3 * indp - 1] = Load[:, 3]

    truss = {
        'CM': BarCM,
        'Node': Node,
        'Bars': Bars,
        'Trigl': Trigl,
        'B': B,
        'L': L,
        'FixedDofs': np.unique(rs),
        'A': Abar
    }

    angles = {
        'CM': RotSpring,
        'fold': Fold,
        'bend': Bend,
        'kpf': kpf * np.ones(Fold.shape[0]),
        'kpb': kpb * np.ones(Bend.shape[0]),
        'pf0': pf0,
        'pb0': np.ones_like(pb0) * np.pi,
        'Panel': Panel
    }

    return truss, angles, F