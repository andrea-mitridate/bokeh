from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, Div, Label
from bokeh.plotting import figure
from bokeh.models.formatters import FuncTickFormatter

from spectrum_calculator import *


# Set up data
log10_T = -2.5
log10_H_star_on_beta = -1.
log10_alpha = -0.5
log10_eta = -1.
h2_omega_bub, h2_omega_sw, h2_omega_turb, h2_omega_sum, k_sw, k_bub, v_w = pt_h2_omega(f, log10_T, log10_H_star_on_beta, log10_alpha, log10_eta)

source_tot = ColumnDataSource(data=dict(x=f, y=h2_omega_sum))
source_bub = ColumnDataSource(data=dict(x=f, y=h2_omega_bub))
source_sw = ColumnDataSource(data=dict(x=f, y=h2_omega_sw))
source_turb = ColumnDataSource(data=dict(x=f, y=h2_omega_turb))


info_box = Label(x=15, y=450, x_units='screen', y_units='screen',
                 text="runaway = " + str(runaway_Q(log10_alpha, log10_eta)) + r'   with   𝜅_ϕ=' + "%.3f" % round(k_bub, 2) + ', 𝜅_sw=' + "%.3f" % round(k_sw, 2) + ', and v=' + "%.3f" % round(v_w, 2), render_mode='css',
                 border_line_color='black', border_line_alpha=0.0,
                 text_font_size='20px',
                 background_fill_color='yellow', background_fill_alpha=0.3)

# Set up plot
plot = figure(plot_height=600, plot_width=800,
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[f[0], f[-1]], y_range=[10**-14, 10**-6],
              x_axis_type="log", y_axis_type="log",
              x_axis_label='f [Hz]', y_axis_label= 'h^2 Ω(f)',
              outline_line_color='black')

plot.title.text = "Interactive GW Spectrum"
plot.title.align = "right"
plot.title.text_color = "black"
plot.title.text_font_size = "25px"


plot.line('x', 'y', source=source_bub, line_width=5, line_alpha=0.6, legend_label='Bubble', line_color="orange")
plot.line('x', 'y', source=source_sw, line_width=5, line_alpha=0.6, legend_label='Sw', line_color="green")
plot.line('x', 'y', source=source_turb, line_width=5, line_alpha=0.6, legend_label='Turb', line_color="blue")
plot.line('x', 'y', source=source_tot, line_width=2, line_alpha=0.6, legend_label='Sum', line_color="black")
plot.circle(freq_bins[:5], h2omega_median[:5], color='red', size=5, line_alpha=0)
plot.xaxis.axis_label_text_font_size = '15pt'
plot.yaxis.axis_label_text_font_size = '15pt'
plot.xaxis.major_label_text_font_size = "15pt"
plot.yaxis.major_label_text_font_size = "15pt"
plot.legend.label_text_font_size = '15pt'

plot.add_layout(info_box)


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

err_xs = []
err_ys = []

for x, y, yerr_low, yerr_hig in zip(freq_bins[:5], h2omega_median[:5], yerr[0], yerr[1]):
    err_xs.append((x, x))
    err_ys.append((y - yerr_low, y + yerr_hig))

# plot them
plot.multi_line(err_xs, err_ys, color='red')

# Set up widgets

T = Slider(title="log T", value=-2.5, start=-7.0, end=2.0, step=0.1)
H_star_on_beta = Slider(title="log(H/β)$", value=-1.0, start=-5.0, end=0., step=0.1)
alpha = Slider(title="log ɑ", value=-0.5, start=-1., end=1., step=0.1)
eta = Slider(title="log η", value=-1.0, start=-1.5, end=1.0, step=0.1)


def update_data(attrname, old, new):

    # Get the current slider values
    log10_T = T.value
    log10_H_star_on_beta = H_star_on_beta.value
    log10_alpha = alpha.value
    log10_eta = eta.value

    # Generate the new curve
    h2_omega_bub, h2_omega_sw, h2_omega_turb, h2_omega_sum, k_sw, k_bub, v_w = pt_h2_omega(f, log10_T, log10_H_star_on_beta, log10_alpha, log10_eta)

    source_tot.data = dict(x=f, y=h2_omega_sum)
    source_bub.data = dict(x=f, y=h2_omega_bub)
    source_sw.data = dict(x=f, y=h2_omega_sw)
    source_turb.data = dict(x=f, y=h2_omega_turb)

    info_box.text = "runaway = " + str(runaway_Q(log10_alpha, log10_eta)) + r'   with   𝜅_ϕ=' + "%.3f" % round(k_bub, 2) + ', 𝜅_sw=' + "%.3f" % round(k_sw, 2) + ', and v=' + "%.3f" % round(v_w, 2)

for w in [T, H_star_on_beta, alpha, eta]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = column(T, H_star_on_beta, alpha, eta)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().info_box = "Spectrum"