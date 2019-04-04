
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from collections import OrderedDict
def load_data(infile):
    import json
    data = None
    infile = infile.replace(".json", '')
    infile = infile.replace(".p", '')
    infile = infile.replace(".pickle", '')
    for ext in ['.json', '.pickle']:
        try:
            with open(infile + ext, 'r') as fp:
                data = json.load(fp)
        except Exception,e:
            print("Exception."),str(e)
            continue
    return data


#sqdm_prop = load_data("../out/tmp/usa_sqdm_properties.json")
sqdm_prop = load_data("../out/tmp/2USA_sqdm_prop.json")
#sqdm_prop = load_data("../out/tmp/2Psqdm_prop.json")


print sqdm_prop.keys()
print max(sqdm_prop['xkey_max_seg_holding'].values())

a =sqdm_prop["count_xcolumns_nsplits_yblocks"]
dic_count_nsplit_nyblocks = OrderedDict(sorted(a.items())) #sort from lowkey to high key.

seg_delx = sqdm_prop['seg_delx']
seg_dely = sqdm_prop['seg_dely']
unsignedseg_dely = [abs(int(dely)) for dely in seg_dely]
seg_lengths = sqdm_prop['seg_lengths']
split_count = sqdm_prop['split_counts']
slopes = []
for dely in seg_dely:
    if dely <0:
        slopes +=[2]
    elif dely > 0:
        slopes +=[4]
    else:
        slopes += [16]
##
cnyb = []
cns = []
for xkey,ns_nyb in dic_count_nsplit_nyblocks.items():
    ns, nyb = ns_nyb
    cnyb +=[nyb]
    cns +=[ns]

print max(cns),min(cns), max(cnyb), min(cnyb)

x = cns

figtitle="Histogram-Segment's split counts"+"(N="+str(len(x))+")"
figtitle="Histogram-Segment's X-spans"+"(N="+str(len(x))+")"
figtitle="Histogram-Segment Lengths"+"(N="+str(len(x))+")"
figtitle="Histogram-Number Of Y-Blocks in Vertical Columns"+"(N="+str(len(x))+")"
figtitle="Histogram-Number Of Splits Segments in Vertical Columns"+"(N="+str(len(x))+")"
xvarname = "Segment's Length"
xvarname = "Number of Y-Blocks"
xvarname = "Number of Splits"
print("max,min,len"), len(x), max(x), min(x)
##
nbins=10
#x = [math.log((v+1),2) for v in x]

N, bins, patches = plt.hist(x, color='#0504aa', alpha=0.7, rwidth=1, bins=nbins)
fracs = N / N.max()
# we need to normalize the data to 0..1 for the full range of the colormap
norm = colors.Normalize(fracs.min(), fracs.max())
# Now, we'll loop through our objects and set the color of each accordingly
for thisfrac, thispatch in zip(fracs, patches):
    color = plt.cm.viridis(norm(thisfrac))
    thispatch.set_facecolor(color)
plt.ylabel('Frequency')
#plt.xlabel('$Log_2$'+'('+xvarname+')')
plt.xlabel(xvarname)
plt.title(figtitle)
#plt.savefig("../out/tmp/"+figtitle+".png")
#plt.show()

def test():
    n=100
    N = len(cns[0:n])

    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars

    fig = plt.figure()
    ax = fig.add_subplot(111)


    rects1 = ax.bar(ind, cns[0:n], width, color='royalblue')
    rects2 = ax.bar(ind+width, cnyb[0:n], width, color='seagreen')

    # add some
    ax.set_ylabel('Vertical Columns at X')
    ax.set_title('Number of Segments and Y-Blocks in each Vertical Columns')
    #ax.set_xticks(ind + width / 2)
    #ax.set_xticklabels( [str(k)[0:3] for k in dic_count_nsplit_nyblocks.keys()[0:n]] )

    ax.legend( (rects1[0], rects2[0]), ('#SegmentSplits', '#Y-Blocks') )

    plt.show()

def stacked_bar(data, series_labels, category_labels=None,
                show_values=False, value_format="{}", y_label=None,
                grid=True, reverse=False):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list
                       containing data for each series in rows
    series_labels   -- list of series labels (these appear in
                       the legend)
    category_labels -- list of category labels (these appear
                       on the x-axis)
    show_values     -- If True then numeric value labels will
                       be shown on each bar
    value_format    -- Format string for numeric value labels
                       (default is "{}")
    y_label         -- Label for y-axis (str)
    grid            -- If True display grid
    reverse         -- If True reverse the order that the
                       series are displayed (left-to-right
                       or right-to-left)
    """

    ny = len(data[0])
    ind = list(range(ny))

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    for i, row_data in enumerate(data):
        axes.append(plt.bar(ind, row_data, bottom=cum_size,
                            label=series_labels[i]))
        cum_size += row_data

    if category_labels:
        plt.xticks(ind, category_labels)

    if y_label:
        plt.ylabel(y_label)

    plt.legend()

    #if grid:
    #    plt.grid()

    if show_values:
        for axis in axes:
            for bar in axis:
                w, h = bar.get_width(), bar.get_height()
                plt.text(bar.get_x() + w/2, bar.get_y() + h/2,
                         value_format.format(h), ha="center",
                         va="center")
import numpy as np
import matplotlib.pyplot as plt
def plot_cns_cnyb(cns,cnyb):
    plt.figure(figsize=(20, 20))

    series_labels = ['# Splits', '# Y-Blocks']

    data = [
        cns[0:60000],
        cnyb[0:60000]
    ]

    category_labels = [] #['Cat A', 'Cat B', 'Cat C', 'Cat D']

    stacked_bar(
        data,
        series_labels,
        category_labels=category_labels,
        show_values=False,
        value_format="{:.1f}",
        y_label="Numbers/Counts"
    )
    plt.title("Bar-Diagram: # of Splits & # of Y-Blocks in Vertical Columns")
    plt.xlabel("Vertical Columns of X-Values")
    plt.savefig('../out/tmp/Bar-Diagram- # of Splits & # of Y-Blocks.png')
    plt.show()
