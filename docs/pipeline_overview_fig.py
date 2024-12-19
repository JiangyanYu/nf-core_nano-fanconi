#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 09:29:11 2024

@author: jiangyanyu
"""

import matplotlib.pyplot as plt
import itertools

# Define nodes and their positions
nodes = {
    "fast5": (1, 4),
    "pod5": (0.3,3.5),
    "bam": (1.7,3.5),
    "dorado": (1, 3),
    "pycoqc": (1, 2),
    "minimap2": (2, 2),
    "samtools": (3, 2),
    "sniffles": (4, 2),
    "AnnotSV": (5, 3),
    "whatshap": (5, 1),
    "DeepVariant": (6.3, 1),
}

# Define edges (connections) as provided
edges = [
    ("fast5", "dorado"),
    ("pod5", "dorado"),
    ("bam", "dorado"),
    ("dorado", "pycoqc"),
    ("pycoqc", "minimap2"),
    ("minimap2", "samtools"),
    ("samtools", "sniffles"),
    ("sniffles", "whatshap"),
    ("sniffles","AnnotSV"),
    ("whatshap","DeepVariant")
]

# Node background colors
node_colors = {
    "fast5": "lightgrey",
    "pod5": "lightgrey",
    "bam": "lightgrey",
    "AnnotSV": "orange",
    "whatshap":"lightgreen",
    "DeepVariant":"orange",
    "default": "pink",
}

# Define legend labels and their colors
legend_labels = {
    #"Raw data input": "lightgrey",
    "Finished": "pink",
    "Ongoing": "lightgreen",
    "Pending": "orange",
}

# Initialize the plot
fig, ax = plt.subplots(figsize=(9, 5))

# Plot edges with straight lines
for start, end in edges:
    x1, y1 = nodes[start]
    x2, y2 = nodes[end]
    ax.plot([x1, x2], [y1, y2], color="black", linewidth=2)

# Plot nodes
for label, (x, y) in nodes.items():
    # Use specific color if defined, otherwise use default
    background_color = node_colors.get(label, node_colors["default"])
    ax.scatter(x, y, s=100, color="red", zorder=3)  # Node
    ax.text(
        x, y, label, fontsize=12, ha="center", va="center", color="black", zorder=4,
        bbox=dict(boxstyle="circle", facecolor=background_color, edgecolor="none")
    )
    
# Add "nano-Fanconi" text at coordinates (3, 3)
ax.text(3, 5, "nano-Fanconi", fontsize=25, ha="center", va="center", color="green", zorder=4)

# Add a legend
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label)
           for label, color in legend_labels.items()]
ax.legend(handles=handles, loc="upper right", title="Developing status", fontsize=10, title_fontsize=12)



# Customize plot appearance
ax.set_xlim(0, 7)
ax.set_ylim(0, 5)
ax.axis("off")  # Hide axes

# Save and show the figure
plt.savefig("workflow_complete_graph.png", dpi=300)
plt.show()
