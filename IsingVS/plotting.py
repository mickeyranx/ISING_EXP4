import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 

path="C:/Users/Miki/uni/CPPrakt/ISING_V4/ISING_EXP4/IsingVS/"
#path="C:/Users/Miki/work/uni/CPPrakt/ISING_V4/ISING_EXP4/IsingVS/"
#determine optimal multihit parameter and thermal steps
def exercise3a():
    
    df = pd.read_table(path+"testing_params.txt") 
    e = df.iloc[:, 1]
    m = df.iloc[:, 3]
    it = np.arange(0,len(m), step = 1)

    fig, ax = plt.subplots()
    ax.scatter(it, m, marker="s", facecolors = "none", edgecolors="b")
    ax.set_xlabel("no. sweeps")
    ax.set_ylabel("magnetisation per particle")
    ax.set_title(r'$\beta = 0.6, N = 5$')
    ax.grid(True)
    plt.show()

    fig, ax = plt.subplots()
    ax.scatter(it, e, marker="s", facecolors = "none", edgecolors="b")
    ax.set_xlabel("no. sweeps")
    ax.set_ylabel("energy density")
    ax.set_title(r'$\beta = 0.6, N = 5$')
    ax.grid(True)
    plt.show()





exercise3a()




