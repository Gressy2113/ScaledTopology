import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def plot_prof(prof, color, xlim = (0, 1), ylim = (-30, 15)):
    fig, ax = plt.subplots(figsize = (5.5, 3))

    reaction = np.linspace(0, 1, len(prof))
    prof -= prof[0] #np.mean(prof[reaction < 25])
    prof = prof[::-1]
    plt.plot(reaction, prof, '.-', color = color,
             markeredgecolor = 'k', 
             markeredgewidth = 0.5
             )
    
    # plt.scatter(reaction, prof, s=10, color = color, 
    #         edgecolors='k', linewidth = 0.5,
    #         zorder=4)

    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(0, 1.1, 0.1)))
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(0, 1.1, 0.05)))

    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(-50, 50, 10)))
    ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(-50, 50, 5)))


    plt.grid(True, which='major', linestyle=  '-')
    plt.grid(True, which='minor', linestyle=  '-', lw=0.2)
    plt.ylim(ylim)
    plt.xlim(xlim)
    ax.set_ylabel("Free energy, kJ/mol")
    ax.set_xlabel("Reaction progress")

