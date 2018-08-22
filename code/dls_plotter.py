# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import plotly
%matplotlib inline
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
plotly.tools.set_credentials_file(username='na542@msstate.edu', api_key='fC8lnTDCUQ5G3zDxrqeW')


# <codecell>

dfmscounty = pd.read_csv('./results/dls_stats/mscounty.csv',delimiter=',')
dfuscounty = pd.read_csv('./results/dls_stats/uscounty.csv',delimiter=',')
dfuscong = pd.read_csv('./results/dls_stats/cong.csv',delimiter=',')
dfusstate = pd.read_csv('./results/dls_stats/usstate.csv',delimiter=',')

dfmscounty.columns = ['data','pid','inputlines','splits','cgrayrec','cgreenrec','cwhiterec','dlsbytes']
dfuscounty.columns = ['data','pid','inputlines','splits','cgrayrec','cgreenrec','cwhiterec','dlsbytes']
dfusstate.columns = ['data','pid','inputlines','splits','cgrayrec','cgreenrec','cwhiterec','dlsbytes']
dfuscong.columns = ['data','pid','inputlines','splits','cgrayrec','cgreenrec','cwhiterec','dlsbytes']

frames= [dfmscounty,dfuscounty,dfusstate,dfuscong]
combDf = pd.concat(frames)
combDf.shape
combDf.head(4)

# <codecell>

import pandas as pd
import random
#df = pd.read_csv('./results/dls_stats/raw_data_all_cat.csv',delimiter=',')

dfselected = combDf[combDf['cgrayrec'] > 1]
print("selected all the rows with gray> 5"), dfselected.shape

mscountyDf =dfselected[dfselected['data'] == 'mscounty']
uscountyDf =dfselected[dfselected['data'] == 'uscounty']
uscongDf =dfselected[dfselected['data'] == 'cong']
usstateDf =dfselected[dfselected['data'] == 'usstate']

print("Selected data from each category")
print mscountyDf.shape, uscountyDf.shape, uscongDf.shape, usstateDf.shape
#randomly select data.(5, 8) (362, 8) (293, 8) (21, 8)

rl = random.sample(xrange(len(mscountyDf)), int(len(mscountyDf) * 40/100))
rl1 = random.sample(xrange(len(uscountyDf)), int(len(uscountyDf) * 5/100))
rl2 = random.sample(xrange(len(uscongDf)), int(len(uscongDf) * 5/100))
rl3 = random.sample(xrange(len(usstateDf)), int(len(usstateDf) * 5/100))

pmscDf = mscountyDf.iloc[rl,:]
puscDf = uscountyDf.iloc[rl1,:]
puscongcDf = uscongDf.iloc[rl2,:]
pusstateDf = usstateDf.iloc[rl3,:]

combPlotDf = pd.concat([pmscDf,puscDf, puscongcDf, pusstateDf])
print combPlotDf.shape

# <codecell>


# <codecell>

combDf.sort(['inputlines']).head(10)

# <codecell>

print combDf[combDf['data']=='cong'].min()
combDf.groupby(combDf['data'],).describe()

# <codecell>


# <codecell>

combDf.groupby(['data']).mean()

# <codecell>

#Determine max and min.
print("maximum")
print mscountyDf.max(), 
print 
print
print("minimum")
print mscountyDf.min()

# <codecell>


# <codecell>

category_names = list(combPlotDf['pid'])
print category_names
# Create the |general blog and the "subplots" i.e. the bars
f, ax1 = plt.subplots(1, figsize=(16,5))

# Set the bar width
bar_width = 0.5

# positions of the left bar-boundaries
bar_l = [i+1 for i in range(len(combPlotDf['inputlines']))]

# positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+(bar_width/2) for i in bar_l]

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        # using the pre_score data
        combPlotDf['inputlines'],
        # set the width
        width=bar_width,
        # with the label pre score
        label='Input-Lines',
        # with alpha 0.5
        alpha=0.5,
        # with color
        color='blue')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        combPlotDf['splits'],
        width=bar_width,
        bottom=combPlotDf['inputlines'],
        label='Split-Lines',
        alpha=0.5,
        color='red')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        combPlotDf['cgrayrec'],
        width=bar_width,
        bottom=[i+j for i,j in zip(combPlotDf['inputlines'],combPlotDf['splits'])],
        label='Gray',
        alpha=0.5,
        color='gray')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        combPlotDf['cgreenrec'],
        width=bar_width,
        bottom= [i+j+k for i,j,k in zip(combPlotDf['inputlines'],combPlotDf['splits'],combPlotDf['cgrayrec'])] ,
        label='Green',
        alpha=0.5,
        color='green')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        combPlotDf['cwhiterec'],
        width=bar_width,
        bottom=[i+j+k+l for i,j,k,l in zip(combPlotDf['inputlines'],combPlotDf['splits'],combPlotDf['cgrayrec'],combPlotDf['cgreenrec'])],
        label='White/Clear',
        alpha=0.5,
        color='yellow')

# set the x ticks with names
plt.xticks(tick_pos, category_names,rotation='vertical') #df3T.index)

# Set the label and legends
ax1.set_ylabel("# of Lines/Recs")
ax1.set_xlabel("Data Set")
plt.legend(loc='upper left')
# Set a buffer around the edge
plt.minorticks_on()
ax1 = plt.gca()
ax1.grid(True)
plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])

plt.savefig('./results/dls_stats/multi_poly.png', bbox_inches='tight')

# <codecell>

#Plotting means

df = combDf

#find out means for each of the columns
ms = list(df[df["data"] == 'mscounty'].iloc[0:,2:].mean(axis=0)) #ignore first 2
usst = list(df[df["data"] == 'uscounty'].iloc[0:,2:].mean(axis=0))
uscong = list(df[df["data"] == 'usstate'].iloc[0:,2:].mean(axis=0))
uscounty = list(df[df["data"] == 'cong'].iloc[0:,2:].mean(axis=0))

raw_data = {'stat': ['inputlines', 'splits','cgrayrec','cgreenrec','cwhiterec','dlsbytes'],
          #'MS-Counties', 'US-Counties', 'US-State', 'US-Congress'
        'msc':ms[0:],
        'usc':usst[0:],
        'uss': usst[0:],
        'uscong': uscong[0:]}

df3 = pd.DataFrame(raw_data, columns = ['stat','msc', 'usc', 'uss', 'uscong'])
df3T   = df3.T.iloc[1:,:]
df3T.columns= df3['stat']
df3T

# <codecell>

category_names = ['MS County', 'US County','US State', 'US Cong. D']

category_names = category_names
# Create the |general blog and the "subplots" i.e. the bars
f, ax1 = plt.subplots(1, figsize=(10,5))

# Set the bar width
bar_width = 0.5

# positions of the left bar-boundaries
bar_l = [i+1 for i in range(len(df3T['inputlines']))]

# positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+(bar_width/2) for i in bar_l]

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        # using the pre_score data
        df3T['inputlines'],
        # set the width
        width=bar_width,
        # with the label pre score
        label='Input-Lines',
        # with alpha 0.5
        alpha=0.5,
        # with color
        color='blue')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        df3T['splits'],
        width=bar_width,
        bottom=df3T['inputlines'],
        label='Split-Lines',
        alpha=0.5,
        color='red')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        df3T['cgrayrec'],
        width=bar_width,
        bottom=[i+j for i,j in zip(df3T['inputlines'],df3T['splits'])],
        label='Gray',
        alpha=0.5,
        color='gray')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        df3T['cgreenrec'],
        width=bar_width,
        bottom= [i+j+k for i,j,k in zip(df3T['inputlines'],df3T['splits'],df3T['cgrayrec'])] ,
        label='Green',
        alpha=0.5,
        color='green')

# Create a bar plot, in position bar_1
ax1.bar(bar_l,
        df3T['cwhiterec'],
        width=bar_width,
        bottom=[i+j+k+l for i,j,k,l in zip(df3T['inputlines'],df3T['splits'],df3T['cgrayrec'],df3T['cgreenrec'])],
        label='White/Clear',
        alpha=0.5,
        color='yellow')

# set the x ticks with names
plt.xticks(tick_pos, category_names) #df3T.index)

# Set the label and legends
ax1.set_ylabel("# of Lines/Recs")
ax1.set_xlabel("Data Set")
plt.legend(loc='upper left')
# Set a buffer around the edge
plt.minorticks_on()
ax1 = plt.gca()
ax1.grid(True)
plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])

plt.savefig('./results/dls_stats/means_plot.png', bbox_inches='tight')

# <codecell>

df3T

# <codecell>


