
import matplotlib.pyplot as plt
import pandas as pd

inp_dir  = './data/' #input file directory
out_dir =  './plots/'

basics = pd.read_csv(inp_dir + 'title.basics.tsv.gz', sep='\t', na_values='\\N')
ratings = pd.read_csv(inp_dir + 'title.ratings.tsv.gz', sep='\t', na_values='\\N')

basics.sample(frac=0.1)
#basics.info()
#ratings.info()
#basics.describe()
#ratings.describe()
print(basics.info())
print(basics.describe())

movies = basics[basics['titleType'] == 'movie']
movies['runtimeMinutes'] = movies['runtimeMinutes'].astype(float)
#movies['runtimeMinutes'] = pd.to_numeric(movies['runtimeMinutes'])

movies_rating = pd.merge(movies, ratings, how='left', on='tconst')

#basics[basics'startYear']
#movies_rating.info()
st_war = movies_rating[movies_rating['primaryTitle'].str.contains( 'Star Wars' )]
# movies_rating[movies_rating['primaryTitle'].str.lower().apply( lambda x: 'Star Wars' in x )]
sorting = movies_rating.sort_values(by=['averageRating'], ascending=False)
sorting = sorting[sorting['numVotes'] > 10000]

#movies_rating.hist()
#plt.show()

movies_per_year = movies_rating['startYear'].value_counts().sort_index()
movies_per_year.plot()
movies_per_year.iloc[:-1].plot()
plt.savefig(out_dir + 'movies_counts')

##### transforming to series
movies_rating['genres'] = movies_rating['genres'].astype(str) 
sr = pd.Series(movies_rating['genres'])

sr_genresCount = sr.value_counts().sort_values(ascending=False) #sorting by value
sr_genresCount_indx = sr.value_counts().sort_index() # sorting by index

plot_genres = sr_genresCount.iloc[0:10] #sorted by values
plot_genres.plot(kind='bar', rot=45, color='r')
plt.savefig(out_dir + 'genres_counts')

################################ over time
genres_overtime = movies_rating.groupby(by=['startYear'])['genres'].value_counts().sort_values(ascending=False)
genres_overtime_indx = movies_rating.groupby(by=['startYear'])['genres'].value_counts().sort_index()

fig, ax = plt.subplots(figsize=(15,7))
genres_overtime.iloc[0:10].plot(ax=ax)
genres_overtime.iloc[0:10].plot(ax=ax, kind='bar', rot=25, color='b', fontsize=7, title='top 10 genre over time sorted by value') # sorted by value, the top 10 sorted by value include only 2 genres
plt.savefig(out_dir + 'genres_overtime_stacked')
genres_overtime.iloc[0:10].unstack().plot(ax=ax)
plt.savefig(out_dir + 'genres_overtime')

genres_overtime_indx.iloc[0:10].plot(kind='bar', rot=25, color='g', fontsize=5, title='top 10 genre over time sorted by index')
plt.savefig(out_dir + 'genres_overtime_indexSorted')

########### SORRY I forgot this step :(
proportions = (movies_rating['genres'].value_counts() / sum(movies_rating['genres'].value_counts())).sort_index()
proportion_overtime = movies_rating.groupby(by=['startYear'])['genres'].value_counts() / sum(movies_rating.groupby(by=['startYear'])['genres'].value_counts())
proportion_overtime_indxSorted = proportion_overtime.sort_index()
proportion_overtime_indxSorted[0:10].plot(kind='bar', rot=25, color='k', fontsize=6, title='proportion of top 10 genres over time sorted by index' )
plt.savefig(out_dir + 'proportion_overtime')

###############################################
movies_rating.groupby(by=['genres']).describe()
movies_rating.groupby(by=['genres']).mean()
movies_rating.groupby(by=['genres']).min()

movies_rating.plot(movies_rating['runtimeMinutes'],movies_rating['averageRating'], kind='scatter')
#####!!!!! I get error for above command and couldn't solve the issue

print()


