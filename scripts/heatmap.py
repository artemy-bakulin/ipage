import matplotlib.pyplot as plt
import numpy as np
def draw_heatmap(p_values):
    up_regulated=lambda x: sum(p_values[x][:len(p_values[x])//2])<=sum(p_values[x][len(p_values[x])//2:])
    down_regulated=lambda x: sum(p_values[x][:len(p_values[x])//2])>=sum(p_values[x][len(p_values[x])//2:])
    p_names=list(filter(up_regulated, p_values))+list(filter(down_regulated,p_values))
    p_values=[p_values[name] for name in p_names]
    plt.rcParams.update({'font.weight': 'roman'})
    plt.rcParams.update({'ytick.labelsize': 20})
    fontsize_pt = plt.rcParams['ytick.labelsize']
    dpi = 72.27
    matrix_height_pt = fontsize_pt * len(p_values)+300
    matrix_height_in = matrix_height_pt / dpi
    matrix_width_pt=fontsize_pt*len(p_values[0])+300
    matrix_width_in=matrix_width_pt/dpi

    top_margin = 0.04  # in percentage of the figure height
    bottom_margin = 0.04 # in percentage of the figure height
    figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
    figure_width=matrix_width_in
    fig, ax = plt.subplots(
        figsize=(figure_width,figure_height), 
        gridspec_kw=dict(top=1-top_margin, bottom=bottom_margin))
    
    ax.set_xticks(np.arange(0))
    ax.set_yticks(np.arange(len(p_names)))
    ax.set_yticklabels(p_names)
    plt.tick_params(axis='y', which='both', labelleft=False, labelright=True,length=0)
    fig.autofmt_xdate()
    plt.set_cmap('RdBu_r')
    #ax.set_title('')
    im = ax.imshow(p_values)
    #fig.tight_layout()
    #color_bar_
    cbaxes = fig.add_axes([0, 0.3 , 0.35/figure_width, 0.4])  # This is the position of the colorbar
    hm=plt.colorbar(im, cax = cbaxes)
    #ax.set_xlabel('Expression increases –>', size=25,weight='roman')
    plt.savefig('heatmap.jpg', bbox_inches='tight')
    #plt.show()
    plt.close()