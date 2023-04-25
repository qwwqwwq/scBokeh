import logging
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
            completions=self.ad.var_names.tolist(),
            title="Enter Gene Name (e.g. HES4 ):", value="HES4")

        self.category_to_colormap = self.get_categorical_variables_and_colormaps()

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
            data_source_dict[f"{attr}1"] = self.ad.obsm[attr][:, 0].tolist()
            data_source_dict[f"{attr}2"] = self.ad.obsm[attr][:, 1].tolist()

        for k in self.category_to_colormap:
            data_source_dict[k] = self.ad.obs[k].to_list()

        self.selected_points = None

        self.data_source = ColumnDataSource(
            data_source_dict
        )

        self.scalar_scatters = {}
        self.categorical_scatters = {}

        for attr in obsm_keys:
            categorical_fig, categorical_scatter = self.create_categorical_obsm_figure(attr)
            scalar_fig, scalar_scatter = self.create_scalar_obsm_figure(attr)

            self.categorical_scatters[attr] = (categorical_fig, categorical_scatter)
            self.scalar_scatters[attr] = (scalar_fig, scalar_scatter)

        self.violin_plot_source = ColumnDataSource(
            data=dict(xs=[], ys=[], xj=[], yj=[], color=[])
        )

        self.obsm_figures = {}

        self.mapper = linear_cmap(
            field_name='expression',
            palette="Viridis256",
            low=0,
            high=1)

        self.color_bar = ColorBar(
            color_mapper=self.mapper['transform'],
            width=10,
            major_label_text_font_size="10pt",
            location=(0, 0))

        self.categorical_legend = Legend(items=[], location='top_right')

        self.tools = "pan,wheel_zoom,box_zoom,reset,lasso_select,box_select"

        self.selected_umis = []

        self.violin_figure, self.violin_patches = self.create_violin_plot()

        def on_select_change(attr, old, new):
            if len(new) == 0:
                self.selected_points = None
                return

            self.selected_points = pd.Series(
                [x in self.ad.obs.index[new] for x in self.ad.obs.index],
                index=self.ad.obs.index
            ).map({True: 'Selected', False: 'Not Selected'})
            self.update_violin()

        self.data_source.selected.on_change("indices", on_select_change)

        def on_gene_symbol_change(attr, old, new):
            self.update_gene()

        def on_categorical_variable_change(attr, old, new):
            self.update_factor()

        ### input on change
        self.gene_symbol_input.on_change('value', on_gene_symbol_change)
        self.categorical_variable_select.on_change('value', on_categorical_variable_change)

        user_input_block = column(self.gene_symbol_input, self.categorical_variable_select)

        graph_column = []

        for attr in obsm_keys:
            graph_column.append(
                row(self.categorical_scatters[attr][0], self.scalar_scatters[attr][0], sizing_mode="stretch_both")
            )
        graph_column.append(self.violin_figure)

        graph_layout = column(*graph_column, sizing_mode="stretch_both")

        self.layout = row(graph_layout, user_input_block, sizing_mode="stretch_both")

        self.update_factor()
        self.update_gene()

    def get_umi(self, feature):
        return self.ad.obs_vector(feature)

    def create_categorical_obsm_figure(self, attr):
        fig = figure(title="UMAP", tools=self.tools)
        fig.toolbar.logo = None
        fig.xaxis.axis_label = f"{attr}1"
        fig.yaxis.axis_label = f"{attr}2"
        fig.select(BoxSelectTool).select_every_mousemove = False
        fig.select(LassoSelectTool).select_every_mousemove = False
        fig.add_layout(self.categorical_legend, "right")

        scatter = fig.circle(f"{attr}1", f"{attr}2",
                             size=3,
                             source=self.data_source,
                             line_color=None,
                             fill_color={'field': self.categorical_variable_select.value,
                                         'transform': self.category_to_colormap[
                                             self.categorical_variable_select.value]},
                             selection_color="orange",
                             alpha=0.6,
                             nonselection_alpha=0.1,
                             selection_alpha=0.4)

        return fig, scatter

    def create_scalar_obsm_figure(self, attr):
        fig = figure(title="UMAP", tools=self.tools)
        fig.toolbar.logo = None
        fig.xaxis.axis_label = f"{attr}1"
        fig.yaxis.axis_label = f"{attr}2"
        fig.select(BoxSelectTool).select_every_mousemove = False
        fig.select(LassoSelectTool).select_every_mousemove = False

        scatter = fig.circle(f"{attr}1", f"{attr}2",
                             size=3,
                             source=self.data_source,
                             line_color=None,
                             fill_color={'field': self.categorical_variable_select.value,
                                         'transform': self.category_to_colormap[attr]},
                             selection_color="orange",
                             alpha=0.6,
                             nonselection_alpha=0.1,
                             selection_alpha=0.4)

        scatter.add_layout(self.color_bar, 'right')

        return fig, scatter

    def create_violin_plot(self):
        fig = figure()
        patches = fig.patches(xs='xs',
                              ys='ys',
                              alpha=0.6,
                              fill_color='color',
                              line_color='black',
                              source=self.violin_plot_source)

        fig.toolbar.logo = None
        fig.yaxis.axis_label = "Expression"

        fig.outline_line_color = None
        fig.background_fill_color = "#efefef"
        fig.xaxis.major_label_orientation = np.pi / 4
        fig.xgrid.grid_line_color = None
        fig.ygrid.grid_line_color = "#dddddd"
        fig.axis.minor_tick_line_color = None
        fig.axis.major_tick_line_color = None
        fig.axis.axis_line_color = None

        return fig, patches

    def update_legend(self):
        self.categorical_legend.items = []
        attr = self.categorical_variable_select.value
        unique_categories = self.ad.obs[attr].unique().to_list()
        renderers = [x[1] for x in self.categorical_scatters.values()]

        for idx, cat in enumerate(sorted(unique_categories)):
            self.categorical_legend.items.append(
                LegendItem(label=str(cat),
                           renderers=renderers,
                           index=idx)
            )

    def update_factor(self):
        attr = self.categorical_variable_select.value

        for _, (_, scatter) in self.categorical_scatters.items():
            scatter.glyph.fill_color = {'field': attr, 'transform': self.category_to_colormap[attr]}

        self.update_legend()
        self.update_violin()

    def update_gene(self):
        gene_symbol = self.gene_symbol_input.value.strip()

        # select gene expression value
        umis = self.get_umi(gene_symbol)
        ## update existing tranform
        vmin, vmax = np.percentile(umis, [2, 98])
        self.mapper['transform'].low = vmin
        self.mapper['transform'].high = vmax
        ## update title
        for _, (fig, _) in self.categorical_scatters.items():
            fig.title.text = gene_symbol
        for _, (fig, _) in self.scalar_scatters.items():
            fig.title.text = gene_symbol
        ## update source data
        self.data_source.data.update(expression=umis, )

        self.update_violin()

    def update_violin(self, cut=2, bins=1000):
        gene_symbol = self.gene_symbol_input.value.strip()
        umis = self.get_umi(gene_symbol)

        if self.selected_points is None:
            categorical_array = self.ad.obs[self.categorical_variable_select.value]
        else:
            categorical_array = self.selected_points
        # update axis
        unique_categories = sorted(categorical_array.unique())
        x_range = (np.arange(0, len(unique_categories)) * 3).tolist()
        ## update data
        color = []
        xs = []
        ys = []

        for idx, (i, category_value) in enumerate(zip(x_range, unique_categories)):
            umi = umis[categorical_array == category_value]
            try:
                kde = gaussian_kde(umi)
            except (np.linalg.LinAlgError, ValueError):
                continue
            color.append(self.color_palette[:len(unique_categories)][idx])
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

        self.violin_plot_source.data = dict(xs=xs, ys=ys, color=color)
        self.violin_figure.xaxis.ticker = FixedTicker(ticks=x_range)
        self.violin_figure.xaxis.major_label_overrides = {k: str(v) for k, v in zip(x_range, unique_categories)}
        self.violin_figure.title.text = gene_symbol

    def get_categorical_variables_and_colormaps(self):
        categorical_vars = {}
        for k in self.ad.obs:
            if self.ad.obs[k].dtype.name == 'category':
                factors = self.ad.obs[k].unique().to_list()
                categorical_vars[k] = CategoricalColorMapper(
                    factors=factors,
                    palette=self.color_palette[:len(factors)])

        return categorical_vars


def main():
    logger.info("Running scBokeh")
    fname = join(DATA_DIR, "pbmc68k_reduced.h5ad")

    logger.info("Reading {}".format(fname))
    ad = sc.read_h5ad(fname)
    viz = SingleCellViz(ad)
    curdoc().add_root(viz.layout)
    curdoc().title = "scBokeh"


if __name__ == '__main__':
    main()
