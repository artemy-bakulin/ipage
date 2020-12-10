import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib



def columnwise_heatmap(array, ax=None, expression=False, cmap_main='RdBu_r', cmap_reg='YlOrBr', **kw):
    #ax = ax or plt.gca()
    images = []
    if expression:
        current_cmap = matplotlib.cm.get_cmap(cmap_reg)
        current_cmap.set_bad(color='black')
        im = ax[0].imshow(np.atleast_2d(array[:, 0]).T, cmap=current_cmap, **kw)
        images.append(im)
        im = ax[1].imshow(np.atleast_2d(array[:, 1:]), cmap=cmap_main, **kw)
    else:
        im = ax.imshow(np.atleast_2d(array[:, :]), cmap=cmap_main, **kw)

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
        colorbar_names = ['Regulon\'s \n enrichment']
        colorbar_images = [-1]
    elif n == 2:
        colorbar_names = ['Regulator\'s \n expression', 'Regulon\'s \n enrichment']
        colorbar_images = [0, 1]
    for i in colorbar_images:
        cax = fig.add_subplot(gs[i // cols, i % cols])
        fig.colorbar(ims[i], cax=cax)
        cax.set_title(colorbar_names[i], fontsize=10)


def draw_heatmap(names, values, output_name='output_ipage', expression=None, cmap_main='RdBu_r', cmap_reg='RdBu_r'):

    if type(names[0]) != list:
        df = pd.DataFrame(values, index=names)
    else:
        df = pd.DataFrame(values, index=names[0], columns=names[1])

    if expression:
        df.insert(0, 'regulator', expression)
    plt.rcParams.update({'font.weight': 'roman'})
    plt.rcParams.update({'ytick.labelsize': 10})
    fontsize_pt = plt.rcParams['ytick.labelsize']
    dpi = 72.27
    matrix_height_pt = (fontsize_pt+30/2) * df.shape[0]
    matrix_height_in = matrix_height_pt / dpi
    matrix_width_pt = (fontsize_pt+50/2) * df.shape[1]
    matrix_width_in = matrix_width_pt / dpi
    top_margin = 0.04  # in percentage of the figure height
    bottom_margin = 0.04  # in percentage of the figure height / (1 - top_margin - bottom_margin)
    figure_height = matrix_height_in
    figure_width = matrix_width_in

    if expression:
        fig, ax = plt.subplots(1, 2, figsize=(figure_width, figure_height), gridspec_kw={'width_ratios': [1.8, len(df)-1]})
    else:
        fig, ax = plt.subplots(1, 1, figsize=(figure_width, figure_height))

    ims = columnwise_heatmap(df.values, ax=ax, aspect="auto", expression=bool(expression),
                             cmap_main=cmap_main, cmap_reg=cmap_reg)
    if expression:
        ax[0].set(xticks=[], yticks=np.arange(len(df)), yticklabels=df.index, xlabel='Regulator')
        ax[0].xaxis.set_label_position('top')
        ax[1].set(xticks=[], yticks=[], xlabel='Regulon')
        ax[1].xaxis.set_label_position('top')
    else:

        ax.set(xticks=[], yticks=np.arange(len(df)), yticklabels=df.index, xlabel='Regulon')
        ax.xaxis.set_label_position('top')
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
        plt.close()
