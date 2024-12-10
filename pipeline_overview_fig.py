#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 09:29:11 2024

@author: jiangyanyu
"""

import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import numpy as np

# Define nodes and their positions
nodes = {
    "Input": (1, 1),
    "QC": (3, 2),
    "Alignment": (5, 3),
    "Analysis": (7, 2),
    "Report": (9, 1),
}

# Define edges (connections)
edges = [
    ("Input", "QC"),
    ("QC", "Alignment"),
    ("Alignment", "Analysis"),
    ("Analysis", "Report"),
]

# Define curved paths for edges
def curved_path(start, end, curvature=0.3):
    x1, y1 = start
    x2, y2 = end
    mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2 + curvature
    vertices = [start, (mid_x, mid_y), end]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    return Path(vertices, codes)

# Initialize the plot
fig, ax = plt.subplots(figsize=(10, 5))

# Plot edges with curved lines
for start, end in edges:
    path = curved_path(nodes[start], nodes[end], curvature=0.5)
    patch = PathPatch(path, edgecolor="blue", linewidth=2, linestyle="--", facecolor="none")
    ax.add_patch(patch)

# Plot nodes
for label, (x, y) in nodes.items():
    ax.scatter(x, y, s=100, color="red", zorder=3)  # Node
    ax.text(x, y, label, fontsize=12, ha="center", va="center", color="white", zorder=4, 
            bbox=dict(boxstyle="circle", facecolor="red", edgecolor="none"))

# Customize plot appearance
ax.set_xlim(0, 10)
ax.set_ylim(0, 5)
ax.axis("off")  # Hide axes

# Save and show the figure
plt.savefig("subway_map.png", dpi=300)
plt.show()

