#!/usr/bin/env python3
import csv

# 1) Print the fixed iTOL Color Strip header
print("""DATASET_COLORSTRIP
SEPARATOR COMMA
DATASET_LABEL,Type Annotation
COLOR,#ff0000

LEGEND_TITLE,Types
LEGEND_SHAPES,1,1,1,1,1,1,1
LEGEND_COLORS,#FF0000,#0000FF,#00FF00,#FFA500,#800080,#808080,#00FFFF
LEGEND_LABELS,gii,crispr,g2l,retron,DGR,Unknown,abi

DATA""")

# 2) Define your annotation->color mapping exactly in the same order as LEGEND_LABELS
color_map = {
    "gii":     "#FF0000",  # matches LEGEND_LABELS "gii"
    "crispr":  "#0000FF",  # matches "crispr"
    "g2l":     "#00FF00",  # matches "g2l"
    "retron":  "#FFA500",  # matches "retron"
    "DGR":     "#800080",  # matches "DGR"
    "Unknown": "#808080",  # matches "Unknown"
    "abi":     "#00FFFF",  # matches "abi"
}

# 3) Read your CSV and print the DATA lines
with open("all_headers.csv", "r", newline="") as infile:
    reader = csv.reader(infile)
    # If your CSV has a header row (e.g., "sequence_id,annotation"), you can skip it:
    # next(reader)

    for row in reader:
        sequence_id, annotation = row[0].strip(), row[1].strip()
        
        # Look up the color for this annotation; if not found, default to gray
        color = color_map.get(annotation, "#808080")
        
        # Print a line in iTOL format: "sequence_id,color,annotation"
        print(f"{sequence_id},{color},{annotation}")
