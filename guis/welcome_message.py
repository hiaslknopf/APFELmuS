import numpy as np
import matplotlib.pyplot as plt

""" A nice little welcome message for the GUI plot windows """

def welcome_message(ax, text):
    # Randomly distribute colorful mus :)
    mus = [dict(x=0.5, y=0.5, fontsize=200, ha='center', va='center')]
    for _ in range(50):
        mus.append(dict(x=np.random.rand(), y=np.random.rand(), fontsize=np.random.randint(10, 30),
                        ha='center', va='center'))

    cmap = plt.get_cmap('hsv')
    for i, mu in enumerate(mus):
        ax.text(mu['x'], mu['y'], '$\mu$', fontsize=mu['fontsize'], ha=mu['ha'], va=mu['va'],
                color=cmap(np.random.rand()), transform=ax.transAxes, zorder=-1)
    ax.axis('off')

    # Add text
    rect = plt.Rectangle((0, 0.25), 1, 0.35, color='orange', alpha=0.5, transform=ax.transAxes)  # zorder to put it behind the text and with some
    ax.add_patch(rect)
    ax.text(0.5, 0.55, 'Welcome to', fontsize=20, ha='center', va='center', transform=ax.transAxes)
    ax.text(0.5, 0.425, 'APFELmuS', fontsize=30, ha='center', va='center', transform=ax.transAxes, weight='bold')
    ax.text(0.5, 0.36, 'A Python Framework for the EvaLuation of microdosimetric Spectra', fontsize=8, ha='center', va='center', transform=ax.transAxes)
    ax.text(0.5, 0.3, text, fontsize=15, ha='center', va='center', transform=ax.transAxes)