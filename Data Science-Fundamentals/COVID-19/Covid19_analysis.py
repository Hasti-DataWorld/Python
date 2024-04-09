## Analyzing Covid-19 data and creating visualizations, by HNateghi

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# read data and have a general look at the columns
df = pd.read_csv('owid-covid-data.csv')
#print(df)
#df.info()

df['date']= pd.to_datetime(df['date'],format='%Y-%m-%d')
df = df.loc[:,['date','new_cases','new_deaths', 'population', 'continent','location','female_smokers', 'male_smokers']]

newcase_perday = df.groupby(by=['date'])['new_cases'].sum()
deathcase_perday = df.groupby(by=['date'])['new_deaths'].sum()

fig,ax = plt.subplots(figsize=(12,6))

#rolling_avg = 14
#ax.plot(df.date, newcase_perday,color='red',label='total new cases per day')
newcase_perday.rolling(center=True, win_type='triang', window=15).mean().plot(ax=ax, x_compat=True, rot=35, color='r', linestyle='--', label='Global daily new cases')
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.set_yscale('log')
ax.set_ylabel('Total daily new cases')
ax.set_xlabel('Date')
ax.grid()

ax2=ax.twinx()
deathcase_perday.rolling(center=True, win_type='triang', window=15).mean().plot(ax=ax2, x_compat=True, rot=35, color='gold',label='Global daily new deaths' )
leg = ax.legend( loc='lower right', frameon=False) #bbox_to_anchor=(0.5, -0.10)
leg2=ax2.legend( loc='upper left', frameon=False) #,bbox_to_anchor=(0.49, -0.15)
ax2.grid(False) # turn off grid for second Y axis
ax2.set_ylabel('Total daily new deaths')
ax2.set_yscale('log')

plt.show()

#//////////////////////////////////////////////////////////////////////////////////
# total number of new cases in each continent normalized by population. 
#since we need to consider the factor of population of continents

df['continent'] = df['continent'].astype(str)
df['case_pop'] = df['new_cases']/df['population']
df['death_pop'] = df['new_deaths']/df['population']

casepop_perday = df.groupby(by=['continent'], dropna=True)['case_pop'].sum().dropna()
deathpop_perday = df.groupby(by=['continent'], dropna=True)['death_pop'].sum().dropna()
casepop_perday = casepop_perday.drop('nan')
deathpop_perday = deathpop_perday.drop('nan')

fig,ax = plt.subplots(figsize=(6,4))

ax.scatter( casepop_perday.index , casepop_perday.values, color='m', alpha=1, marker='o', label='infected cases normalized by population')
#'total number of new cases normalized by population (2020/02 - 2021/03)'
ax.set_ylabel('Total new cases/population (2020/02-2021/03)', fontsize=8)
ax.set_xlabel('Date-continent')
ax.set_xticklabels(casepop_perday.index, rotation=25, fontsize=7)

ax2=ax.twinx()
ax2.scatter(deathpop_perday.index ,deathpop_perday.values, color='c', alpha=0.9, marker='x', label='total death normalized by population')
#leg = ax.legend( loc='upper left', frameon=False) #bbox_to_anchor=(0.5, -0.10)
#leg2=ax2.legend( loc='upper right', frameon=False) #,bbox_to_anchor=(0.49, -0.15)
fig.legend(loc='upper right', fontsize=7)
#'total number of new deaths normalized by population (2020/02 - 2021/03)'
ax2.set_ylabel('Total new deaths/population (2020/02-2021/03)', fontsize=8)
plt.show()

#//////////////////////////////////////////////////////////////////////////////////

df['tot_smoker'] = df['female_smokers'] + df['male_smokers']

#df['population'] = df['population'].dropna()
#df['tot_smoker_norm'] = df['tot_smoker']/df['population'] 
df = df.dropna(subset=['tot_smoker'])
df['tot_smoker'] = df['tot_smoker'].sort_values()
#smoker = df.groupby(by=['location'])['tot_smoker'].value_counts()  >> counts the country repeating
smoker = df.groupby(by=['location'])['tot_smoker'].mean().sort_values(ascending=True)

df[df['location'] == 'India'] #22.5
df[df['location'] == 'Brazil'] #28
df[df['location'] == 'Greece'] #87
df[df['location'] == 'Russia'] #81
df[df['location'] == 'Mexico'] #28.3
df[df['location'] == 'Serbia'] #77.9

list_cont = ['India','Brazil','Mexico', 'Serbia','Russia','Greece']
list_smoker = ['India: 22.5', 'Brazil:28', '* Mexico:28.3', 'Serbia:77.9','Russia:81','Greece:87']
col = ['c','g','b','r','k','y']

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for i,cont in enumerate(list_cont):
    df_n = df[df['location'] == cont]
    ax1.plot(df_n['date'], (df_n['new_cases']/df_n['population']).rolling(center=True,window=15).mean(), color=col[i], label=list_smoker[i])
    ax2.plot(df_n['date'], (df_n['new_deaths']/df_n['population']).rolling(center=True,window=15).mean(), color=col[i], label=list_smoker[i])
#'relation between share of people who smoke and daily infections, most recent year available'
ax1.set_title('Relation between share of people who smoke and daily infections')
ax2.set_xlabel('Date')
ax1.set_ylabel("Daily new cases / population")
ax2.set_ylabel("Daily new deaths / population")
ax1.legend()
ax2.legend()
plt.show()

#just an idea
death_loc = df.groupby(by=['location'])['new_deaths'].mean()
smok_loc = df.groupby(by=['location'])['tot_smoker'].mean()

plt.scatter(death_loc.index, death_loc.values, c=smok_loc.values, cmap='gist_rainbow')
plt.ylim(-2,90)
plt.xticks(rotation=90, fontsize=4, ha='left')
plt.colorbar()
plt.colo
plt.show()


#//////////////////////////////////////////////////////////////////////////////////
# total number of new cases in each continent normalized by population. 
#we need to consider the factor of population of continents

df['continent'] = df['continent'].astype(str)
df = df.dropna(subset=['continent', 'new_cases', 'population', 'new_deaths', 'date'])
df['continent'] = df['continent'].dropna()
df['new_cases'] = df['new_cases'].dropna()
df['population'] = df['population'].dropna()
df['new_deaths'] = df['new_deaths'].dropna()
df['date'] = df['date'].dropna()
df['case_pop'] = df['new_cases']/df['population']
df['death_pop'] = df['new_deaths']/df['population']

casepop_perday = df.groupby(by=['date', 'continent'], dropna=True)['case_pop'].sum()
deathpop_perday = df.groupby(by=['date', 'continent'], dropna=True)['death_pop'].sum()

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)    

casepop_perday.unstack().rolling(center=True, window=7).mean().plot(ax=ax1, rot=45, alpha=1)
ax1.xaxis.set_major_locator(mdates.MonthLocator())
ax1.set_ylabel('number of new cases/population')
ax2.set_xlabel('Date')

deathpop_perday.unstack().rolling(center=True, window=7).mean().plot(ax=ax2, rot=45, alpha=0.8, linestyle='--')
ax2.xaxis.set_major_locator(mdates.MonthLocator())
leg = ax1.legend( loc='upper left', frameon=True) #bbox_to_anchor=(0.5, -0.10)
leg2=ax2.legend( loc='upper left', frameon=True) #,bbox_to_anchor=(0.49, -0.15)
ax2.set_ylabel('number of new deaths/population')
plt.show()



