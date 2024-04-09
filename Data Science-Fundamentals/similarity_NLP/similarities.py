### inaugural similarities code - Hasti Nateghi

import pandas as pd
from nltk.corpus import inaugural
from nltk.tokenize import word_tokenize
import matplotlib.pylab as plt
import seaborn as sns
from nltk.stem.wordnet import WordNetLemmatizer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
# ---------------------------------------------------------------------------------------------------------------------

def get_wordnet_pos(word):
    # Get Part-of-speech (POS) tag of input word, and return the first POS tag character (which is
    # the character that lemmatize() accepts as input)

    from nltk import pos_tag
    from nltk.corpus import wordnet

    tag = pos_tag([word])[0][1][0].upper()
    tag_dict = {'J': wordnet.ADJ,
                'N': wordnet.NOUN,
                'V': wordnet.VERB,
                'R': wordnet.ADV}

    return tag_dict.get(tag, wordnet.NOUN)  
# ---------------------------------------------------------------------------------------------------------------------

df = pd.DataFrame([(x[:-4],inaugural.raw(x)) for x in inaugural.fileids()],
                  columns=["filenames","talks"])

df["words"] = df["talks"].apply(word_tokenize) # tokenize each inaugural address

df["word_length"] = df["words"].apply(len)

plt.figure(figsize=(12,8))
ax = sns.barplot(x="filenames",y="word_length",data=df)
plt.xticks(rotation=90)
plt.subplots_adjust(bottom=0.2)
plt.savefig('Inaugural.png')

##########################################################

# Lemmatize the corpus 
lem = WordNetLemmatizer()
df["words_lem"] = df["words"].apply(lambda l: [lem.lemmatize(word,pos=get_wordnet_pos(word)) for word in l])

# Obtain the TF-IDF and save in a new DataFrame called tfidf
tfidf_vectorizer = TfidfVectorizer(stop_words='english')
tfidf_feature = tfidf_vectorizer.fit_transform(df["talks"])
tfidf = pd.DataFrame(data=tfidf_feature.todense(),
                     columns=tfidf_vectorizer.get_feature_names())

row_numbers = range(0,len(df))
top5 = []
for row in row_numbers:
    a = tfidf.iloc[row]
    b = a.sort_values(ascending=False)[0:5].index
    top5.append(b)
df["tfidf_top5"] = top5

# obtain the cosine similarity
cosine_sim = cosine_similarity(tfidf)

# Plot the cosine similarity matrix
plt.figure(figsize=(15,15))
sns.heatmap(cosine_sim, cmap='Reds', linewidth=0.5)
plt.yticks(rotation=0)
plt.title("Cosine similarities for inaugural speeches")
plt.savefig('similarities_naugurals_exersice3.png');

