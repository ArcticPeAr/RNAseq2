  ####Library imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns




def BarPlot(goDict, versus):
    #plot horizontal bar plot of the number of up and down genes for each GO term with wrapped labels with gene count on top of bar
    #make a list of the number of up and down genes for each GO term
    upList = []
    downList = []
    for key in goDict:
        upList.append(len(goDict[key][0]))
        downList.append(len(goDict[key][1]))
    #make a list of the GO terms
    goList = list(goDict)
    #make a dataframe of the lists
    df = pd.DataFrame(list(zip(goList, upList, downList)), columns =['GO term', 'Up', 'Down'])
    #make a bar plot
    plt.figure(figsize=(18,10))
    ax = sns.barplot(y="GO term", x="value", hue="variable", data=pd.melt(df, ['GO term']))
    ax.set_yticklabels(ax.get_yticklabels(), rotation=70, ha="right")
    plt.xlabel('Gene count')
    plt.ylabel('GO term')
    plt.title(f'Number of up and down genes for each GO term searched for in {versus}')
    plt.tight_layout()
    #add gene count on top of bar
    for p in ax.patches:
        ax.annotate(str(p.get_width()), (p.get_width(), p.get_y()+0.55*p.get_height()))
    plt.show()