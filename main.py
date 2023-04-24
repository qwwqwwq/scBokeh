import logging
from functools import lru_cache
from os.path import dirname, join

import anndata
import colorcet as cc
import numpy as np
import pandas as pd
import scanpy as sc
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import BoxSelectTool, LassoSelectTool
from bokeh.models import CategoricalColorMapper
from bokeh.models import ColorBar
from bokeh.models import ColumnDataSource, Select, Legend, LegendItem
from bokeh.models import FixedTicker
from bokeh.models.widgets import AutocompleteInput, Div
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from scipy.stats.kde import gaussian_kde

DATA_DIR = dirname(__file__)

logger = logging.getLogger(__file__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())


@lru_cache()
def load_h5ad():
    fname = join(DATA_DIR, "pbmc68k_reduced.h5ad")
    return sc.read_h5ad(fname)


def get_umi(anndat, features):
    umi = anndat.obs_vector(features)
    return umi


def get_categorical_variables_and_colormaps(ad: anndata.AnnData,
                                            color_palette):
    categorical_vars = {}
    for k in anndat.obs:
        if ad.obs[k].dtype.name == 'category':
            factors = ad.obs[k].unique().to_list()
            categorical_vars[k] = CategoricalColorMapper(
                factors=factors,
                palette=color_palette[:len(factors)])

    return categorical_vars


class SingleCellViz:
    def __init__(self, ad: anndata.AnnData,
                 color_palette=cc.b_glasbey_bw_minc_20,
                 obsm_keys=('X_umap', 'X_pca'),
                 ):
        self.ad = ad
        self.n_cells = ad.shape[0]
        self.n_genes = ad.shape[1]
        self.color_palette = color_palette
        self.gene_symbol_input = AutocompleteInput(
            completions=anndat.var_names.tolist(),
            title="Enter Gene Name (e.g. HES4 ):", value="HES4")

        self.category_to_colormap = get_categorical_variables_and_colormaps(
            self.ad,
            self.color_palette)

        if len(self.category_to_colormap) == 0:
            self.categorical_variable_select = Select(
                title="Legend:",
                value=list(self.category_to_colormap.keys())[0],
                options=list(self.category_to_colormap.keys()))
        else:
            self.categorical_variable_select = None

        data_source_dict = {
            "expression": [0] * self.n_cells,
        }

        for attr in obsm_keys:
            data_source_dict[attr + "_1"] = self.ad.obsm[attr][:, 0].tolist()
            data_source_dict[attr + "_1"] = self.ad.obsm[attr][:, 1].tolist()

        for k in self.category_to_colormap:
            data_source_dict[k] = anndat.obs[k].to_list()

        self.data_source = ColumnDataSource(
            data_source_dict
        )

        self.obsm_figures = {}

        self.mapper = linear_cmap(
            field_name='expression',
            palette="Viridis256",
            low=low,
            high=high)

        self.color_bar = ColorBar(
            color_mapper=self.mapper['transform'],
            width=10,
            major_label_text_font_size="10pt",
            location=(0, 0))

        self.categorical_legend = Legend(items=[], location='top_right')

    def create_categorical_obsm_figure(self, attr):
        fig = figure(title="UMAP", tools=tools)
        fig.toolbar.logo = None
        fig.xaxis.axis_label = f"{attr}1"
        fig.yaxis.axis_label = f"{attr}2"
        fig.select(BoxSelectTool).select_every_mousemove = False
        fig.select(LassoSelectTool).select_every_mousemove = False
        fig.add_layout(categorical_legend, "right")

        scatter = fig.circle(f"{attr}1", f"{attr}2",
                             size=3,
                             source=source,
                             line_color=None,
                             fill_color={'field': self.categorical_variable_select.value,
                                         'transform': self.category_to_colormap[dd]},
                             selection_color="orange",
                             alpha=0.6,
                             nonselection_alpha=0.1,
                             selection_alpha=0.4)

        return fig, scatter

    def create_scalar_obsm_figure(self, attr):
        fig = figure(title="UMAP", tools=tools)
        fig.toolbar.logo = None
        fig.xaxis.axis_label = f"{attr}1"
        fig.yaxis.axis_label = f"{attr}2"
        fig.select(BoxSelectTool).select_every_mousemove = False
        fig.select(LassoSelectTool).select_every_mousemove = False
        fig.add_layout(categorical_legend, "right")

        scatter = u2.scatter(f"{attr}1", f"{attr}2",
                             size=3,
                             source=source,
                             fill_color=mapper,
                             line_color=None,
                             selection_color="orange",
                             alpha=0.6,
                             nonselection_alpha=0.1,
                             selection_alpha=0.4)
        
        u2.add_layout(color_bar, 'right')
        scatter = fig.circle('UMAP1', 'UMAP2',
                             size=3,
                             source=source,
                             line_color=None,
                             fill_color={'field': self.categorical_variable_select.value,
                                         'transform': self.category_to_colormap[dd]},
                             selection_color="orange",
                             alpha=0.6,
                             nonselection_alpha=0.1,
                             selection_alpha=0.4)

        return fig, scatter

    def update_legend(self):
        self.categorical_legend.items = []
        attr = self.select.value
        unique_categories = anndat.obs[attr].unique().to_list()
        for idx, cat in enumerate(sorted(unique_categories)):
            self.categorical_legend.items.append(
                LegendItem(label=str(cat),
                           renderers=[u1_scatter],
                           index=idx)
            )


anndat = load_h5ad()
umap = anndat.obsm['X_umap']
color_palette = cc.b_glasbey_bw_minc_20
category_to_colormap = get_categorical_variables_and_colormaps(anndat, color_palette)
# set up widgets
# stats = PreText(text='', width=500)
symbol = AutocompleteInput(completions=anndat.var_names.tolist(),
                           title="Enter Gene Name (e.g. HES4 ):", value="HES4")

select = Select(title="Legend:", value="phase",
                options=list(category_to_colormap.keys()))
# message box
message = Div(text="""Input Gene Name:\n Legend Option: """, width=200, height=100)

## setup data
dd = select.value
umis = get_umi(anndat, symbol.value)

source_dict = dict(
    color=[0] * umap.shape[0],
    umis=[0] * umap.shape[0],
    UMAP1=umap[:, 0].tolist(),
    UMAP2=umap[:, 1].tolist(), )

for k in category_to_colormap:
    source_dict[k] = anndat.obs[k].to_list()

source = ColumnDataSource(data=source_dict)

# PC3=pca[:,2].tolist(),))
source_vln = ColumnDataSource(data=dict(xs=[], ys=[], xj=[], yj=[], color=[]))
## setup figures
tools = 'reset,pan,wheel_zoom,box_select,lasso_select,save'

###
catogories = sorted(anndat.obs[dd].unique().tolist())
palette = color_palette[:len(catogories)]
low, high = 0, 1
## transforms
mapper = linear_cmap(field_name='umis', palette="Viridis256",
                     low=low, high=high)
color_bar = ColorBar(color_mapper=mapper['transform'],
                     width=10,
                     major_label_text_font_size="10pt",
                     location=(0, 0))

categorical_legend = Legend(items=[], location='top_right')

u1 = figure(width=1000, height=500, title="UMAP", tools=tools)
u1.toolbar.logo = None
u1.xaxis.axis_label = "UMAP1"
u1.yaxis.axis_label = "UMAP2"
u1.select(BoxSelectTool).select_every_mousemove = False
u1.select(LassoSelectTool).select_every_mousemove = False
u1.add_layout(categorical_legend, "right")

u2 = figure(width=1000, height=500, tools=tools, x_range=u1.x_range, y_range=u1.y_range)
u2.toolbar.logo = None
u2.xaxis.axis_label = "UMAP1"
u2.yaxis.axis_label = "UMAP2"
u2.select(BoxSelectTool).select_every_mousemove = False
u2.select(LassoSelectTool).select_every_mousemove = False

########### clustering plots ##################
## UMAP
u1_scatter = u1.circle('UMAP1', 'UMAP2',
                       size=3,
                       source=source,
                       line_color=None,
                       fill_color={'field': dd, 'transform': category_to_colormap[dd]},
                       selection_color="orange",
                       alpha=0.6,
                       nonselection_alpha=0.1,
                       selection_alpha=0.4)

###### gene expression plots ##################
# tSNE

# UMAP
u2_scatter = u2.scatter('UMAP1', 'UMAP2', size=3, source=source,
                        fill_color=mapper,
                        line_color=None,
                        selection_color="orange",
                        alpha=0.6,
                        nonselection_alpha=0.1,
                        selection_alpha=0.4)
u2.add_layout(color_bar, 'right')
#


## voilinplot
volin = figure(width=1000,
               height=500,
               tools='reset,pan,wheel_zoom,save')
volin.patches(xs='xs',
              ys='ys',
              alpha=0.6,
              fill_color='color',
              line_color='black',
              source=source_vln)

volin.toolbar.logo = None
volin.yaxis.axis_label = "Expression"

volin.outline_line_color = None
volin.background_fill_color = "#efefef"
volin.xaxis.major_label_orientation = np.pi / 4
volin.xgrid.grid_line_color = None
volin.ygrid.grid_line_color = "#dddddd"
volin.axis.minor_tick_line_color = None
volin.axis.major_tick_line_color = None
volin.axis.axis_line_color = None


def update_legend():
    categorical_legend.items = []
    dd = select.value
    unique_categories = anndat.obs[dd].unique().to_list()
    for idx, cat in enumerate(sorted(unique_categories)):
        categorical_legend.items.append(
            LegendItem(label=str(cat),
                       renderers=[u1_scatter],
                       index=idx)
        )


def volin_change(gene, catogory, umis, bins=1000, cut=2):
    # update axis
    cats = sorted(catogory.unique())
    x_range = (np.arange(0, len(cats)) * 3).tolist()
    ## update data
    color = []
    xs = []
    ys = []
    xj = []
    yj = []

    for idx, (i, cat) in enumerate(zip(x_range, cats)):
        umi = umis[catogory == cat]
        try:
            kde = gaussian_kde(umi)
        except (np.linalg.LinAlgError, ValueError):
            continue
        color.append(color_palette[:len(cats)][idx])
        # same default paramter same as seaborn.violinplot
        bw_used = kde.factor * umi.std(ddof=1) * cut
        support_min = umi.min() - bw_used
        support_max = umi.max() + bw_used
        kde_support = np.linspace(support_min, support_max, bins)

        # generate y and x values
        yy = np.concatenate([kde_support, kde_support[::-1]])  # note: think about how bokeh patches work 
        x = kde(kde_support)  ### you may change this if the violin not fit in
        x /= x.max()  # scale the relative area under the kde curve, resulting max density will be 1
        x2 = -x
        xx = np.concatenate([x, x2[::-1]]) + i
        xs.append(xx)
        ys.append(yy)
        # xj.append([i]*len(umi))
        # yj.append(umi)

    source_vln.data = dict(xs=xs, ys=ys, color=color)
    # source_vln.data = dict(xs=xs, ys=ys, color=color, xj=xj, yj=yj)

    volin.xaxis.ticker = FixedTicker(ticks=x_range)
    volin.xaxis.major_label_overrides = {k: str(v) for k, v in zip(x_range, cats)}
    volin.title.text = gene


# set up callbacks
def gene_change():
    ## update text input and select attribute
    gene = symbol.value.strip()
    dd = select.value

    message.text = "Input Gene Name: {g} \nLegend Option: {cat}".format(g=gene, cat=dd)
    # select gene expression value
    umis = get_umi(anndat, gene)
    ## update existing tranform
    vmin, vmax = np.percentile(umis, [2, 98])
    mapper['transform'].low = vmin
    mapper['transform'].high = vmax
    ## update title
    u2.title.text = gene
    ## update source data
    source.data.update(umis=umis, )

    # update volin
    clusters = anndat.obs[dd]
    try:
        volin_change(gene, clusters, umis, bins=1000)
    except Exception as e:
        logger.exception(e)


# set up callbacks
def factor_change():
    ## update text input and select attribute
    gene = symbol.value.strip()
    dd = select.value

    message.text = "Input Gene Name: {g} \nLegend Option: {cat}".format(g=gene, cat=dd)
    # select gene expression value
    umis = get_umi(anndat, gene)
    # update factor color 

    u1_scatter.glyph.fill_color = {'field': dd, 'transform': category_to_colormap[dd]}
    update_legend()

    try:
        volin_change(gene, anndat.obs[dd], umis, bins=1000)
    except Exception as e:
        logger.exception(e)


def on_select_change(attr, old, new):
    if len(new) == 0:
        factor_change()
        return

    ## update text input and select attribute
    gene = symbol.value.strip()
    dd = select.value

    message.text = "Input Gene Name: {g} \nLegend Option: {cat}".format(g=gene, cat=dd)
    # select gene expression value
    umis = get_umi(anndat, gene)
    ## update existing tranform
    vmin, vmax = np.percentile(umis, [2, 98])
    mapper['transform'].low = vmin
    mapper['transform'].high = vmax
    ## update title
    u2.title.text = gene
    ## update source data
    source.data.update(umis=umis, )

    # update volin
    in_selection = pd.Series(
        [x in anndat.obs.index[new] for x in anndat.obs.index],
        index=anndat.obs.index
    ).map({True: 'Selected', False: 'Not Selected'})

    try:
        volin_change(gene, in_selection, umis, bins=1000)
    except Exception as e:
        logger.exception(e)


source.selected.on_change("indices", on_select_change)

### input on change
symbol.on_change('value', lambda attr, old, new: gene_change())
select.on_change('value', lambda attr, old, new: factor_change())

# set up layout
sb = column(symbol, select, message)

layout = row(column(u1, u2, volin), sb, sizing_mode="stretch_both")

# initialize
# update()
gene_change()
factor_change()

curdoc().add_root(layout)
# curdoc.theme = 'dark_minimal'
curdoc().title = "ARPKD"
