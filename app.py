from bokeh.models.glyphs import Circle
from bokeh.models.layouts import ColSizing
import numpy as np
from bokeh.plotting import figure, output_file, show
from bokeh.models import *
from bokeh.io import curdoc
from bokeh.events import *
from math import pi
import scipy.signal as sc
import cmath
from bokeh.io import *
from bokeh.io import show
from bokeh.models import CustomJS, Select,ColumnDataSource


pole_flag = False
zero_flag = False
conj_flag = False
filter1_flag = False
filter2_flag = False
filter3_flag = False
Added_filter_flag= False

x_axis = [-1.5, 1.5]
y_axis = [0, 0]
X_axis = [0, 0]
Y_axis = [-1.5, 1.5]

TOOLS = "tap"

#   Add Plot
P = figure(
    title='Z_Transform',
    x_axis_label='Real',
    y_axis_label='Imaginary',
    plot_width=600,
    plot_height=500,
    match_aspect=True,
    tools=TOOLS
)
P_filter = figure(
    title='Customize your filter !',
    x_axis_label='Real',
    y_axis_label='Imaginary',
    plot_width=600,
    plot_height=500,
    match_aspect=True,
    tools=TOOLS
)

T = figure(
    title='Frequency Response',
    x_axis_label='Frequency',
    y_axis_label='Magnitude',
    y_range=(-60, 60),
    plot_width=600,
    plot_height=500,
    match_aspect=True,
)

H = figure(
    x_axis_label='Frequency',
    y_axis_label='Phase',
    y_range=(-60, 60),
    plot_width=600,
    plot_height=500,
    match_aspect=True,
)

def zplane(b, a):

    # The coefficients are less than 1, normalize the coeficients
    if np.max(b) > 1:
        kn = np.max(b)
        b = b/float(kn)
    else:
        kn = 1

    if np.max(a) > 1:
        kd = np.max(a)
        a = a/float(kd)
    else:
        kd = 1

    # Get the poles and zeros
    p = np.roots(a)
    z = np.roots(b)
    k = kn/float(kd)

    return z, p, k


# add a line renderer with legend and line thickness
P.line(x_axis, y_axis, line_width=2, line_color='grey')
P.line(X_axis, Y_axis, line_width=2, line_color='grey')

P_filter.line(x_axis, y_axis, line_width=2, line_color='grey')
P_filter.line(X_axis, Y_axis, line_width=2, line_color='grey')
#   Render glyph
P.circle(0, 0, radius=1.0, fill_color=None, line_color='red')

P_filter.circle(0, 0, radius=1.0, fill_color=None, line_color='red')

# pole_source_no_conj = ColumnDataSource(data=dict(x=[], y=[]))
# zero_source_no_conj = ColumnDataSource(data=dict(x=[], y=[]))

pole_source = ColumnDataSource(data=dict(x=[], y=[], x_conj=[], y_conj=[]))
zero_source = ColumnDataSource(data=dict(x=[], y=[], x_conj=[], y_conj=[]))

pole_addFilter_source = ColumnDataSource(data=dict(x=[], y=[], x_conj=[], y_conj=[]))
zero_addFilter_source = ColumnDataSource(data=dict(x=[], y=[], x_conj=[], y_conj=[]))
conj_button = Button(label="Add Conj", button_type="primary", width=200)
remove_all_button = Button(label="Remove All", button_type="danger", width=200)
remove_poles_button = Button(label="Remove Poles", button_type="danger", width=200)
remove_zeros_button = Button(label="Remove Zeros", button_type="danger",width=200)
AddFilter_button = Button(label="Add Filter", button_type="primary",width=200)

add_fig_button = Button(label="Add Graph", button_type="success", width=140)

render_pole = P.x(source=pole_source, x='x', y='y', size=10)
render_pole_conj = P.x(source=pole_source, x='x_conj', y='y_conj', size=10)
render_zero = P.circle(source=zero_source, x='x', y='y', size=10)
render_zero_conj = P.circle(source=zero_source, x='x_conj', y='y_conj',
                            size=10)
render_pole_addFilter = P_filter.x(source=pole_addFilter_source, x='x', y='y', size=10)
render_pole_conj_addFilter = P_filter.x(source=pole_addFilter_source, x='x_conj', y='y_conj', size=10)
render_zero_addFilter = P_filter.circle(source=zero_addFilter_source, x='x', y='y', size=10)
render_zero_conj_addFilter = P_filter.circle(source=zero_addFilter_source, x='x_conj', y='y_conj',
                            size=10)

draw_tool1 = PointDrawTool(renderers=[render_pole, render_pole_conj])
draw_tool2 = PointDrawTool(renderers=[render_zero, render_zero_conj])


draw_tool3 = PointDrawTool(renderers=[render_pole_addFilter, render_pole_conj_addFilter])
draw_tool4 = PointDrawTool(renderers=[render_zero_addFilter, render_zero_conj_addFilter])


P.add_tools(draw_tool1, draw_tool2)
P.toolbar.active_tap = draw_tool1

P_filter.add_tools(draw_tool3, draw_tool4)
P_filter.toolbar.active_tap = draw_tool3
poles_ADDfilter_x=[]
poles_ADDfilter_y=[]


zeros_ADDfilter_x = []
zeros_ADDfilter_y = []

def add_conj():
    global conj_flag
    conj_flag = True


def remove_all():
    remove_poles()
    remove_zeros()
    remove_poles_filter()
    remove_zeros_filter()

def remove_poles():
    pole_coordList = []
    pole_source.data = dict(x=[i[0] for i in pole_coordList],
                            y=[i[1] for i in pole_coordList],
                            x_conj=[i[0] for i in pole_coordList],
                            y_conj=[i[1] for i in pole_coordList])


def remove_zeros():
    zero_coordList = []
    zero_source.data = dict(x=[i[0] for i in zero_coordList],
                            y=[i[1] for i in zero_coordList],
                            x_conj=[i[0] for i in zero_coordList],
                            y_conj=[i[1] for i in zero_coordList])
def remove_poles_filter():
    pole_coordList = []
    pole_addFilter_source.data = dict(x=[i[0] for i in pole_coordList],
                            y=[i[1] for i in pole_coordList],
                            x_conj=[i[0] for i in pole_coordList],
                            y_conj=[i[1] for i in pole_coordList])
def remove_zeros_filter():
    zero_coordList = []
    zero_addFilter_source.data = dict(x=[i[0] for i in zero_coordList],
                            y=[i[1] for i in zero_coordList],
                            x_conj=[i[0] for i in zero_coordList],
                            y_conj=[i[1] for i in zero_coordList])
def ADD_FIlter():
    global poles_ADDfilter_x
    global poles_ADDfilter_y


    global zeros_ADDfilter_x 
    global zeros_ADDfilter_y 
    poles_ADDfilter_x=pole_addFilter_source.data.get('x')
    poles_ADDfilter_y=pole_addFilter_source.data.get('y')
    zeros_ADDfilter_x =zero_addFilter_source.data.get('x')
    zeros_ADDfilter_y =zero_addFilter_source.data.get('y')


# pole_button.on_click(add_pole)
# zero_button.on_click(add_zero)
conj_button.on_click(add_conj)
remove_all_button.on_click(remove_all)
remove_poles_button.on_click(remove_poles)
remove_zeros_button.on_click(remove_zeros)

AddFilter_button.on_click(ADD_FIlter)

y_conj_pole = []
y_conj_zero = []


def callback(event):
    global y_conj_pole
    y_conj_pole = []
    
    pole_source.data.update([('x_conj', pole_source.data.get('x'))])
    
    for i in pole_source.data.get('y'):
        y_conj_pole.append(-i)
    pole_source.data.update([('y_conj', y_conj_pole)])

    global y_conj_zero
    y_conj_zero = []
    zero_source.data.update([('x_conj', zero_source.data.get('x'))])
    for i in zero_source.data.get('y'):
        y_conj_zero.append(-i)
    zero_source.data.update([('y_conj', y_conj_zero)])


y_conj_pole_filter = []
y_conj_zero_filter = []  
def callback_filter(event):
    global y_conj_pole_filter
    y_conj_pole_filter = []
    
    pole_addFilter_source.data.update([('x_conj', pole_addFilter_source.data.get('x'))])
    
    for i in pole_addFilter_source.data.get('y'):
        y_conj_pole_filter.append(-i)
    pole_addFilter_source.data.update([('y_conj', y_conj_pole_filter)])

    global y_conj_zero_filter
    y_conj_zero_filter = []
    zero_addFilter_source.data.update([('x_conj', zero_addFilter_source.data.get('x'))])
    for i in zero_addFilter_source.data.get('y'):
        y_conj_zero_filter.append(-i)
    zero_addFilter_source.data.update([('y_conj', y_conj_zero_filter)])

def plot_graph():
    H.renderers = []
    T.renderers = []
    zeros = []
    poles = []

    zero_x = []
    zero_y = []

    for i in zero_source.data.get('x'):
        zero_x.append(i)
    for i in zero_source.data.get('y'):
        zero_y.append(i)

    pole_x = []
    pole_y = []
    for i in pole_source.data.get('x'):
        pole_x.append(i)
    for i in pole_source.data.get('y'):
        pole_y.append(i)

    counter = 0
    for i in pole_source.data.get('x'):
        print(pole_x[counter])
        pole = complex(round((pole_x[counter]), 5), round((pole_y[counter]), 5))
        poles.append(pole)
        counter += 1

    counter = 0
    for i in zero_source.data.get('x'):
        zero = complex(round((zero_x[counter]), 5), round((zero_y[counter]), 5))
        zeros.append(zero)
        counter += 1

    print("zeros")
    print(zeros)
    print("poles")
    print(poles)
    Numerator, Denominator = sc.zpk2tf(zeros, poles, 1.0)
    print("Numerator: ", Numerator)
    print("Denominator: ", Denominator)
    sr = 100
    z, p, k = zplane(b=Numerator, a=Denominator)
    frequencies, frequency_response = sc.freqz_zpk(z, p, k, fs=sr)

    # magnitude
    magnitude = 20 * np.log10(abs(frequency_response)) 
    T.line(frequencies, magnitude, line_width=2, line_color='grey')

    # Phase
    Phase = np.angle(frequency_response)
    H.line(frequencies, Phase, line_width=2, line_color='grey')

#############################################
    
poles_filter1 = [-1.0 + 0.0j , 0.2 + 0.25j , 0.7 + 0.4j]
poles_filter1_x=[-1.0,0.2,0.7]
poles_filter1_y=[0.0,0.25,0.4]
poles_filter2 = [0.3 - 0.3j , -0.25 - 0.25j , 1.0 + 0.0j]
poles_filter2_x=[0.3,-0.25,1.0]
poles_filter2_y=[-0.3,-0.25,0.0]
poles_filter3 = [-0.2 - 0.5j , -0.5 + 0.75j , 0.2 + 0.25j]
poles_filter3_x=[-0.2,-0.5,0.2]
poles_filter3_y=[-0.5,0.75,0.25]


zeros_filter1_x = []
zeros_filter1_y = []
zeros_filter2_x = []
zeros_filter2_y = []
zeros_filter3_x = []
zeros_filter3_y = []



def compute_zerors(list_of_poles,poles_x,poles_y):
    for i in range(3):
        if list_of_poles == poles_filter1 :
            zeros_filter1_x.append(poles_x[i]/abs(list_of_poles[i]))
            zeros_filter1_y.append(poles_y[i]/abs(list_of_poles[i]))

        elif list_of_poles == poles_filter2 :
            zeros_filter2_x.append(poles_x[i]/abs(list_of_poles[i]))
            zeros_filter2_y.append(poles_y[i]/abs(list_of_poles[i]))
        else:
            zeros_filter3_x.append(poles_x[i]/abs(list_of_poles[i]))
            zeros_filter3_y.append(poles_y[i]/abs(list_of_poles[i]))


compute_zerors(poles_filter1,poles_filter1_x,poles_filter1_y)
compute_zerors(poles_filter2,poles_filter2_x,poles_filter2_y)
compute_zerors(poles_filter3,poles_filter3_x,poles_filter3_y)


#############################################

menu = Select(title="all pass filters:", value="no selection", options=["no selection yet","filter 1", "filter 2", "filter 3","Added filter.."],width=200)
def SelectFilter(attr, old, new):
    global filter1_flag
    global filter2_flag
    global filter3_flag
    global Added_filter_flag
    if menu.value == 'filter 1': 
        filter1_flag = True
        pole_source.data.update([('x', pole_source.data.get('x')+poles_filter1_x)])
        pole_source.data.update([('y', pole_source.data.get('y')+poles_filter1_y)])
        zero_source.data.update([('x', zero_source.data.get('x')+zeros_filter1_x)])
        zero_source.data.update([('y', zero_source.data.get('y')+zeros_filter1_y)])

        draw_tool = PointDrawTool(renderers=[render_pole, render_pole_conj])

    elif menu.value == 'filter 2': 
        filter1_flag = True
        pole_source.data.update([('x', pole_source.data.get('x')+poles_filter2_x)])
        pole_source.data.update([('y', pole_source.data.get('y')+poles_filter2_y)])
        zero_source.data.update([('x', zero_source.data.get('x')+zeros_filter2_x)])
        zero_source.data.update([('y', zero_source.data.get('y')+zeros_filter2_y)])

        draw_tool = PointDrawTool(renderers=[render_pole, render_pole_conj])
    elif menu.value == 'filter 3': 
        filter1_flag = True
        pole_source.data.update([('x', pole_source.data.get('x')+poles_filter3_x)])
        pole_source.data.update([('y', pole_source.data.get('y')+poles_filter3_y)])
        zero_source.data.update([('x', zero_source.data.get('x')+zeros_filter3_x)])
        zero_source.data.update([('y', zero_source.data.get('y')+zeros_filter3_y)])

        draw_tool = PointDrawTool(renderers=[render_pole, render_pole_conj])
    elif menu.value == 'Added filter..': 
        Added_filter_flag = True
        pole_source.data.update([('x', pole_source.data.get('x')+poles_ADDfilter_x)])
        pole_source.data.update([('y', pole_source.data.get('y')+poles_ADDfilter_y)])
        zero_source.data.update([('x', zero_source.data.get('x')+zeros_ADDfilter_x)])
        zero_source.data.update([('y', zero_source.data.get('y')+zeros_ADDfilter_y)])

        draw_tool = PointDrawTool(renderers=[render_pole, render_pole_conj])
    elif menu.value == 'no selection': 
        remove_all()
        
menu.on_change('value', SelectFilter)




add_fig_button.on_click(plot_graph)
P.on_event(MouseMove, callback)
P_filter.on_event(MouseMove, callback_filter)
P_filter.on_event(MouseMove, callback)

UnitCircels = Column (P,menu,P_filter,AddFilter_button)
mag_freq_Plots = Column(T, H)
Buttons = Row(conj_button,remove_poles_button, remove_zeros_button, remove_all_button ,add_fig_button)
row = Row(UnitCircels,mag_freq_Plots)
app = Column (Buttons,row)
curdoc().add_root(app)