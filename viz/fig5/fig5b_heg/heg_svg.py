#!/usr/bin/env python3
"""
Plot rounded-edge domain cartoons for every *.gb file in the current directory,
stacking them in a single figure called domain_figure.svg.
"""

import glob
import os

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def parse_domain_features(gbk_file):
    features = []
    seq_len = 0
    for record in SeqIO.parse(gbk_file, "genbank"):
        seq_len = len(record.seq)
        for feat in record.features:
            if feat.type.lower() != "region":
                continue  # skip CDS, gene, etc.
            start = int(feat.location.start)
            end = int(feat.location.end)
            label = feat.qualifiers.get("region_name", [""])[0]
            features.append({"start": start, "end": end, "label": label})
    return features, seq_len


def plot_all_domains(gb_files, output_path="domain_figure.svg"):
    parsed = []
    global_max_len = 0
    for gb in gb_files:
        feats, seq_len = parse_domain_features(gb)
        parsed.append((os.path.basename(gb), feats, seq_len))
        global_max_len = max(global_max_len, seq_len)

    if not parsed:
        print("No Region features found in any .gb file.")
        return

    n_tracks   = len(parsed)
    bar_height = 0.6
    y_gap      = 1.0
    left_margin = 60

    fig, ax = plt.subplots(figsize=(12, 1 + 1.2 * n_tracks))

    for idx, (name, feats, seq_len) in enumerate(parsed):
        y = idx * y_gap

        # baseline bar for the whole protein
        ax.add_line(
            plt.Line2D([0, seq_len],
                       [y + bar_height / 2]*2,
                       linewidth=5, color="lightgrey", zorder=0)
        )

        for feat in feats:
            x0     = feat["start"]
            width  = feat["end"] - feat["start"]
            radius = bar_height / 2

            # centre rectangle
            body = mpatches.Rectangle(
                (x0 + radius, y),
                width - 2*radius,               
                bar_height,
                facecolor="grey",
                edgecolor="black",
                linewidth=1,
            )

            # left & right caps
            cap_left  = mpatches.Ellipse(
                (x0 + radius, y + radius),
                width=bar_height, height=bar_height,
                facecolor="grey", edgecolor="black", linewidth=1,
            )
            cap_right = mpatches.Ellipse(
                (x0 + width - radius, y + radius),
                width=bar_height, height=bar_height,
                facecolor="grey", edgecolor="black", linewidth=1,
            )

            for patch in (body, cap_left, cap_right):
                ax.add_patch(patch)

        # track label
        ax.text(-left_margin/2, y + bar_height/2, name,
                ha="right", va="center", fontsize=9)

    ax.set_xlim(-left_margin, global_max_len)
    ax.set_ylim(-0.2, (n_tracks - 1) * y_gap + bar_height + 0.2)
    ax.axis("off")
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved combined figure as {output_path}.")

if __name__ == "__main__":
    gb_files = glob.glob("*.gb")
    if not gb_files:
        print("No .gb files found in the current directory.")
    else:
        plot_all_domains(gb_files)
