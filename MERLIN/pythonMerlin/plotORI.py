import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_ori(Node, Panel, color='skyblue', alpha=1.0):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Plot nodes
    ax.scatter(Node[:, 0], Node[:, 1], Node[:, 2], c='red', s=5, label='Nodes')

    # Plot panels
    faces = [Node[np.array(p)] for p in Panel]  # Already zero-based indexing
    poly = Poly3DCollection(faces, facecolors=color, edgecolors='k', linewidths=0.5)
    poly.set_alpha(alpha)
    ax.add_collection3d(poly)

    # Axes labels and title
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")
    ax.set_title("Interactive 3D Origami Plot")
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), borderaxespad=0.)

    # Equal aspect ratio
    max_range = (Node.max(axis=0) - Node.min(axis=0)).max()
    mid = (Node.max(axis=0) + Node.min(axis=0)) / 2
    ax.set_xlim(mid[0] - max_range/2, mid[0] + max_range/2)
    ax.set_ylim(mid[1] - max_range/2, mid[1] + max_range/2)
    ax.set_zlim(mid[2] - max_range/2, mid[2] + max_range/2)

    # Optional: turn axis on or off
    # ax.axis('off')

    plt.show()