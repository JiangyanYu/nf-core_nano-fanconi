#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 09:29:11 2024

@author: jiangyanyu
"""

import matplotlib.pyplot as plt

# Define nodes and their positions
nodes = {
    "fast5": (1, 4),
    "pod5": (0.3, 3.5),
    "bam": (1.7, 3.5),
    "dorado": (1, 3),
    "pycoqc": (1, 2),
    "pbmm2": (2, 2),
    "sniffles": (3, 3),
    "AnnotSV": (6, 2),
    "whatshap": (4.5, 1),
    "DeepVariant": (3, 1),
}

# Define edges (connections)
edges = [
    ("fast5", "dorado"),
    ("pod5", "dorado"),
    ("bam", "dorado"),
    ("dorado", "pycoqc"),
    ("pycoqc", "pbmm2"),
    ("pbmm2", "sniffles"),
    ("pbmm2", "DeepVariant"),
    ("DeepVariant", "whatshap"),
]

# Define special edges (square stepwise connections)
square_edges = [
    ("sniffles", "AnnotSV"),
    ("whatshap", "AnnotSV"),
]

# Nodes with background color
filled_nodes = {"fast5": "lightgrey", "pod5": "lightgrey", "bam": "lightgrey"}

# Initialize the plot
fig, ax = plt.subplots(figsize=(9, 5))

# Plot straight-line edges
for start, end in edges:
    x1, y1 = nodes[start]
    x2, y2 = nodes[end]
    ax.plot([x1, x2], [y1, y2], color="black", linewidth=2)

# Plot stepwise square edges
for start, end in square_edges:
    x1, y1 = nodes[start]
    x2, y2 = nodes[end]
    mid_x = (x2 + x2) / 2  # Midpoint for stepwise effect
    ax.plot([x1, mid_x, mid_x, x2], [y1, y1, y2, y2], color="black", linewidth=2)

# Plot nodes
for label, (x, y) in nodes.items():
    if label in filled_nodes:
        ax.scatter(x, y, s=100, color="red", zorder=3)  # Node with background
        ax.text(
            x, y, label, fontsize=12, ha="center", va="center", color="black", zorder=4,
            bbox=dict(boxstyle="circle", facecolor=filled_nodes[label], edgecolor="black")
        )
    else:
        ax.text(
            x, y, label, fontsize=12, ha="center", va="center", color="black", zorder=4,
            bbox=dict(boxstyle="square", facecolor="lightblue", edgecolor="black", linewidth=2)
        )

# Add "nano-Fanconi" text at coordinates (3, 5)
ax.text(3, 5, "nano-Fanconi", fontsize=25, ha="center", va="center", color="green", zorder=4)

# Customize plot appearance
ax.set_xlim(0, 7)
ax.set_ylim(0, 5)
ax.axis("off")  # Hide axes

# Save and show the figure
plt.savefig("workflow_complete_graph.png", dpi=300)
plt.show()
