import requests
import json
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

class DisGeNet(object): 
    URL = "http://disgenet.org/api/gda/disease/"

    def __init__(self): 
        pass

    def query_disease_genes(self, disease, namespace="do"):
        r = requests.get(f"{self.URL}/{namespace}/{disease}")
        result = json.loads(r.content.decode("utf-8"))
        df = pd.DataFrame.from_dict(result)
        self.result = df

        return self

    def get_top_genes(self, method="second_derivative"): 
        if not hasattr(self, "result"): 
            raise ValueError("No result found. Please run a query first!")

        if method == 'second_derivative':
            y = np.sort(self.result['score'].values)
            ypp = y[2:] + y[-2] - 2*y[1:-1]

            index = np.argmax(ypp)
            
            try: 
                if len(index):
                    index = index[0]
            except TypeError: 
                pass

            top_n = len(y) - index

            return list(set(self.result.sort_values(by='score', ascending=False)['gene_symbol'].values[:top_n]))

    def plot_score(self): 
        y = np.sort(self.result['score'].values)
        x = np.arange(len(y))

        plt.plot(x,y)
