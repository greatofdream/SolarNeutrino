import matplotlib.colors as mcolors
spectras = ['N13', 'O15', 'F17', 'pp', 'B8', 'hep']
labels = ['$^{13}$N', '$^{15}$O', '$^{17}$F', 'pp', '$^{8}$B', 'hep']
xs = [0.1, 0.1, 0.1, 0.5, 1, 1]
ys = [2, 0.1, 0.1, 0.2, 1.8, 1.8]
colors = list(mcolors.TABLEAU_COLORS.keys())
spectraConfig = dict([(s, {'label': label, 'color': c, 'x': x, 'y':y}) for s, label, c, x, y in zip(spectras, labels, colors, xs, ys)])
# append discrete spectra
