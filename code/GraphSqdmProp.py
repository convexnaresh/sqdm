
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

def sqdm_prop():
    #sqdm_prop = load_data("../out/tmp/usa_sqdm_properties.json")
    sqdm_prop = load_data("../out/tmp/2USA_sqdm_prop.json")



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
    print("len xkey:ns,nyb"), dic_count_nsplit_nyblocks["625476645.0"]
    print("cns,cnyb"), max(cns), min(cns), max(cnyb), min(cnyb)

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

def compare_segment_lengths():
    Nsegs = 'all'
    home = "../out/tmp/Calif/"
    fC = home + "C_lengths-" + str(Nsegs) + ".json"  # original segment
    fCx = home + "ftCatx_lengths-" + str(Nsegs) + ".json"  # segment after split @x.
    fCxy = home + "ftCatxy_lengths-" + str(Nsegs) + ".json"  # segment after split @x and @y.
    fdftC = home + "delegated_ftC_lengths-" + str(Nsegs) + ".json" #segments in ftC after delegation
    DC = load_data(fC)
    DCx = load_data(fCx)
    DCxy = load_data(fCxy)
    DdftC = load_data(fdftC)

    print("Reading lenths form:")
    print fC
    print fCx
    print fCxy
    print fdftC


    print("len(DC)"), len(DC)
    print("len(DCx)"),len(DCx)
    print("len(DCxy)"), len(DCxy)
    print("len(DdC)"), len(DdftC)

    print DdftC.keys()[0:10]
    #print set(DdftC.keys()).difference(set(DCx.keys()))

    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111)

    print DC.keys()[0:5]
    print DCx.keys()[0:5]
    print DCxy.keys()[0:5]
    err_dc =[]

    keys = DC.keys()

    for key in keys:
        try:
            A = (DC[key] - DCx[key])
            B = (DC[key] - DCxy[key])
            C = abs(A) - abs(B)
            D = (DC[key] - DdftC[key])
            err_dc += [(int(key),A,B,C,D)]
        except:
            print("Missing matching keys in DC and DCx"), key

    err_dc = sorted(err_dc, key=lambda x: x[1]) #sort by (Si-Six)
    X = [t[0] for t  in err_dc ][0:] #id
    Y1x = [t[1] for t  in err_dc ][0:] #er
    Y2xy =  [t[2] for t  in err_dc ][0:]
    print X[0:10] #smallest

    err_dc = sorted(err_dc, key=lambda x: x[2]) #sort by (Si-Sixy)
    X = [t[0] for t  in err_dc ][0:] #id
    Y1x = [t[1] for t  in err_dc ][0:] #er
    Y2xy =  [t[2] for t  in err_dc ][0:]
    print X[0:10] #largest


    err_dc = sorted(err_dc, key=lambda x: x[3]) #sort by (Si-Six-(Si-Sixy))
    X = [t[0] for t  in err_dc ][0:] #id
    Y1x = [t[1] for t  in err_dc ][0:] #er
    Y2xy =  [t[2] for t  in err_dc ][0:]
    print X[0:10] #largest


    ax.plot(Y1x, X, 'bo-')
    ax.plot(Y2xy, X, 'ro')
    plt.title("Segment's Manhattan Distance Errors")
    plt.xlabel(r'$d\left(S_i\right)-\sum_{i=0}^ d\left(s_i\right)$'+r'$:s_i=split-of- S_i$')
    plt.ylabel("Segments (Si)")
    plt.legend(loc='upper left')
    # math text
    plt.title("Distribution of Length Errors after Splitting Segments")
    plt.grid()
    #ax.plot([int(k) for k in DC.keys()][0:2000],err_dc_dcx[0:2000], 'bo')  # plot x and y using blue circle markers
    #ax.plot([int(k) for k in DC.keys()][0:2000],err_dc_dcxy[0:2000], 'ro')  # plot x and y using blue circle markers

    #plt.show()
compare_segment_lengths()


