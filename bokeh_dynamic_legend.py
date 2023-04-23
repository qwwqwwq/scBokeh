import pandas as pd
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, Select, Legend, LegendItem, CategoricalColorMapper
from bokeh.plotting import figure
from bokeh.layouts import column

# Assuming you have a pandas DataFrame named 'df' with columns 'x', 'y', 'typeA', 'typeB'
# Here's an example DataFrame with categorical data
data = {'x': [1, 2, 3, 4, 5],
        'y': [2, 4, 1, 3, 6],
        'typeA': ['X', 'Y', 'X', 'Y', 'X'],
        'typeB': ['A', 'A', 'B', 'B', 'B']}

df = pd.DataFrame(data)

source = ColumnDataSource(df)

# Define color mappers for 'typeA' and 'typeB'
color_mapper_typeA = CategoricalColorMapper(factors=['X', 'Y'], palette=['red', 'blue'])
color_mapper_typeB = CategoricalColorMapper(factors=['A', 'B'], palette=['green', 'yellow'])

# Create a scatter plot with dots colored by 'typeA' initially
plot = figure(title='Scatterplot with Color Update')
scatter = plot.scatter(x='x', y='y', source=source, color={'field': 'typeA', 'transform': color_mapper_typeA}, size=10)

# Create a dropdown widget with options 'typeA' and 'typeB'
dropdown = Select(title='Select Color Attribute:', value='typeA', options=['typeA', 'typeB'])

# Define a Python callback function to update the colors and legend based on the dropdown selection
def update(attr, old, new):
    selected_attribute = dropdown.value
    color_mapper = color_mapper_typeA if selected_attribute == 'typeA' else color_mapper_typeB
    scatter.glyph.fill_color = {'field': selected_attribute, 'transform': color_mapper}
    scatter.glyph.line_color = {'field': selected_attribute, 'transform': color_mapper}
    update_legend(selected_attribute, color_mapper)

def update_legend(attr, color_mapper):
    legend.items = []
    unique_categories = df[attr].unique()
    for category in unique_categories:
        color = color_mapper.palette[color_mapper.factors.index(category)]
        legend.items.append(LegendItem(label=category, renderers=[scatter], index=list(df[attr]).index(category)))

# Attach the callback function to the dropdown widget
dropdown.on_change('value', update)

# Create an empty legend
legend = Legend(items=[], location='top_right')
plot.add_layout(legend)

# Update the legend with the initial attribute 'typeA'
update_legend('typeA', color_mapper_typeA)

# Combine the plot and dropdown widget in a layout
layout = column(dropdown, plot)

# Add the layout to the current document and run the Bokeh server
curdoc().add_root(layout)