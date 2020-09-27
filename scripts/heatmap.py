import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def columnwise_heatmap(array, ax=None, expression=False, cmap_main='RdBu_r', cmap_reg='YlOrBr', **kw):
    ax = ax or plt.gca()
    premask = np.tile(np.arange(array.shape[1]), array.shape[0]).reshape(array.shape)
    images = []
    if expression:
        col = np.ma.array(array, mask=premask != 0)
        im = ax.imshow(col, cmap=cmap_reg, **kw)
        images.append(im)
        col = np.ma.array(array, mask=premask == 0)
    else:
        col = np.ma.array(array, mask=premask == -1)
    im = ax.imshow(col, cmap=cmap_main, **kw)
    images.append(im)
    return images


def add_colorbar(fig, ims, n):
    fig.subplots_adjust(left=0.06, right=0.65)
    rows = n
    cols = 1
    gs = GridSpec(rows, cols)
    gs.update(left=0.7, right=0.75, wspace=1, hspace=0.3)
    if n == 0:
        colorbar_names = ['']
        colorbar_images = []
    elif n == 1:
        colorbar_names = ['genes']
        colorbar_images = [-1]
    elif n == 2:
        colorbar_names = ['Regulator', 'Genes']
        colorbar_images = [0, 1]
    for i in colorbar_images:
        cax = fig.add_subplot(gs[i // cols, i % cols])
        fig.colorbar(ims[i], cax=cax)
        cax.set_title(colorbar_names[i], fontsize=10)


def draw_heatmap(names, values, output_name='output_ipage', expression=None, cmap_main='RdBu_r', cmap_reg='YlOrBr'):

    if type(names[0]) != list:
        df = pd.DataFrame(values, index=names)
    else:
        df = pd.DataFrame(values, index=names[0], columns=names[1])

    if expression:
        df.insert(0, 'rbp', expression)
    plt.rcParams.update({'font.weight': 'roman'})
    plt.rcParams.update({'ytick.labelsize': 10})
    fontsize_pt = plt.rcParams['ytick.labelsize']
    dpi = 72.27
    matrix_height_pt = (fontsize_pt+30/2) * df.shape[0]
    matrix_height_in = matrix_height_pt / dpi
    matrix_width_pt = (fontsize_pt+42.5/2) * df.shape[1]
    matrix_width_in = matrix_width_pt / dpi
    top_margin = 0.04  # in percentage of the figure height
    bottom_margin = 0.04  # in percentage of the figure height / (1 - top_margin - bottom_margin)
    figure_height = matrix_height_in
    figure_width = matrix_width_in

    fig, ax = plt.subplots(figsize=(figure_width, figure_height))

    ims = columnwise_heatmap(df.values, ax=ax, aspect="auto", expression=bool(expression),
                             cmap_main=cmap_main, cmap_reg=cmap_reg)

    if type(names[0]) == str:
        ax.set(xticks=np.arange(len(df.columns)), yticks=np.arange(len(df)),
               xticklabels=[''] * len(df.columns), yticklabels=df.index)
    else:
        ax.set(xticks=np.arange(len(df.columns)), yticks=np.arange(len(df)),
               xticklabels=df.columns, yticklabels=df.index)
        plt.xticks(rotation=90)

    # ax.tick_params(bottom=False, top=False,
    #               labelbottom=False, labeltop=True, left=False)
    if expression:
        n = 2
    else:
        n = 1
    add_colorbar(fig, ims, n)
    if output_name == 'stdout':
        plt.show(block=False)
    else:
        plt.savefig('%s.jpg' % output_name, bbox_inches='tight')
