from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, Div, Label, Select
from bokeh.plotting import figure
from bokeh.models.formatters import FuncTickFormatter
from bokeh.themes import built_in_themes

from spectrum_calculator import *
from data_points import *

import numpy as np

###############
# Set up data #
###############
mods = ['Envelope', 'Semi-analytic', 'Numerical']

log10_T = -2.5
log10_H_on_beta = -1.
log10_alpha = -0.5
v = 0.5
mod = 'Semi-analytic'

f = np.logspace(-11., np.log10(3e-8))

Om_phi = h2_omega(f, log10_T, log10_H_on_beta, log10_alpha, 1, contr='bubble', mod=mod)
Om_sw = h2_omega(f, log10_T, log10_H_on_beta, log10_alpha, v, contr='sound')

source_phi = ColumnDataSource(data=dict(x=f, y=Om_phi))
source_sw = ColumnDataSource(data=dict(x=f, y=Om_sw))

##############
# Info boxes #
##############                 
info_box_k_sw = Label(x=15, y=525, x_units='screen', y_units='screen',
                 text='sound wave eff. factor = ' + "%.3f" % round(k_sw(v, 10**log10_alpha), 2),
                 render_mode='css',
                 border_line_color='white', border_line_alpha=0.0,
                 text_font_size='18px',
                 text_color='snow',
                 background_fill_color='yellow', background_fill_alpha=0.0)

###############
# Set up plot #
###############
plot = figure(plot_height=679, plot_width=877,
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[f[0], f[-1]], y_range=[10**-14, 10**-6],
              x_axis_type="log", y_axis_type="log",
              x_axis_label='f [Hz]', y_axis_label= 'h² Ω(f)',
              outline_line_color='snow')

plot.line('x', 'y', source=source_phi, line_width=5, line_alpha=0.9, legend_label='Bubble only', line_color="crimson")
plot.line('x', 'y', source=source_sw, line_width=5, line_alpha=0.9, legend_label='Sound-Wave only', line_color="mediumseagreen")


# layout 
plot.xaxis.axis_label_text_font_size = '16pt'
plot.yaxis.axis_label_text_font_size = '16pt'
plot.xaxis.major_label_text_font_size = "16pt"
plot.yaxis.major_label_text_font_size = "16pt"
plot.legend.label_text_font_size = '15pt'
plot.legend.background_fill_color = "white"
plot.legend.background_fill_alpha = 0.2

plot.title.text = "Interactive GW Spectrum"
plot.title.align = "right"
plot.title.text_color = "white"
plot.title.text_font_size = "23px"
plot.xgrid.visible = True
plot.ygrid.visible = True
plot.outline_line_width = 3

plot.add_layout(info_box_k_sw)


plot.yaxis[0].formatter = FuncTickFormatter(code="""
return 10 + (Math.log10(tick).toString()
             .split('')
             .map(function (d) { return d === '-' ? '⁻' : '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d]; })
             .join(''));
""")
plot.xaxis[0].formatter = FuncTickFormatter(code="""
return 10 + (Math.log10(tick).toString()
             .split('')
             .map(function (d) { return d === '-' ? '⁻' : '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d]; })
             .join(''));
""")

# NANOGrav data
plot.circle(freq_bins[:5], h2omega_median[:5], color='snow', size=8, line_alpha=0)

err_xs = []
err_ys = []

for x, y, yerr_low, yerr_hig in zip(freq_bins[:5], h2omega_median[:5], yerr[0], yerr[1]):
    err_xs.append((x, x))
    err_ys.append((y - yerr_low, y + yerr_hig))

plot.multi_line(err_xs, err_ys, color='snow', line_width=4, line_alpha=0.8)

###################
# Set up slsiders #
###################
T = Slider(title="log(T/GeV)", value=-2.5, start=-7.0, end=2.0, step=0.05)
H_on_beta = Slider(title="log(H/β)", value=-1.0, start=-2.0, end=0., step=0.05)
alpha = Slider(title="log ɑ", value=-0.5, start=-1., end=1., step=0.02)
v_w = Slider(title="v", value=0.5, start=0, end=1.0, step=0.02)
mode = Select(title='Bubble spectrum',value='Semi-analytic', options=mods)

##########################
# define update function #
##########################
def update_data(attrname, old, new):

    # Get the current slider values
    log10_T = T.value
    log10_H_on_beta = H_on_beta.value
    log10_alpha = alpha.value
    v = v_w.value
    mod = mode.value

    # Generate the new curve
    Om_phi = h2_omega(f, log10_T, log10_H_on_beta, log10_alpha, 1, contr='bubble', mod=mod)
    Om_sw = h2_omega(f, log10_T, log10_H_on_beta, log10_alpha, v, contr='sound')

    source_phi.data = dict(x=f, y=Om_phi)
    source_sw.data = dict(x=f, y=Om_sw)

    info_box_k_sw.text = 'sound wave eff. factor = ' + "%.3f" % round(k_sw(v, 10**log10_alpha), 2)

for w in [T, H_on_beta, alpha, v_w, mode]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = column(T, H_on_beta, alpha, v_w, mode)

curdoc().add_root(row(inputs, plot, width=1400))
curdoc().info_box = "Spectrum"
curdoc().theme = 'dark_minimal'