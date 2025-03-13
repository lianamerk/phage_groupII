import glob
import os

from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def parse_gbk_features(gbk_file):
    """Parse a GenBank file and return a list of grey, rectangular GraphicFeatures at level=0."""
    features_list = []
    # Parse the file using Biopython's SeqIO
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feat in record.features:
            # Convert to integer positions
            start = int(feat.location.start)
            end = int(feat.location.end)
            # Create a rectangular, grey feature
            gf = GraphicFeature(
                start=start,
                end=end,
                strand=0,   # 0 means rectangular (no arrow)
                color="grey",
                linewidth=0,
                level=0     # all on the same horizontal band
            )
            features_list.append(gf)
    return features_list


def plot_domain_features(features, max_len, output_path):
    """Given a list of features and a known max length, plot them with a background arrow, saving to output_path."""
    # 1) Determine the bounding box of the features
    starts = [f.start for f in features]
    ends = [f.end for f in features]
    if not starts or not ends:
        # If no features, default to 0..1 to avoid errors
        min_start, max_end = 0, 1
    else:
        min_start, max_end = min(starts), max(ends)

    # 2) Create a GraphicRecord with the desired total length
    record = GraphicRecord(
        sequence_length=max_len,  # accommodate the largest genome
        features=features
    )

    fig, ax = plt.subplots(figsize=(10, 3))

    # 3) Plot features first (no ruler, baseline, or sequence text)
    record.plot(
        ax=ax,
        with_ruler=True,
        draw_line=False,
        plot_sequence=False
        # annotation_height_ratio removed for compatibility
    )

    # 4) Force x-limits so all figures share the same scale from [0..max_len]
    ax.set_xlim(0, max_len)

    # 5) Check y-limits after plotting
    ylim = ax.get_ylim()
    arrow_center = (ylim[0] + ylim[1]) / 2
    arrow_thickness = (ylim[1] - ylim[0]) * 0.8  # how tall the arrow is

    # 6) Add a background arrow from min_start..max_end
    arrow_length = max_end - min_start if max_end > min_start else 1
    background_arrow = mpatches.FancyArrow(
        min_start,
        1,
        arrow_length,
        0,
        width=0.3,
        head_width=0.3,
        head_length=10,
        length_includes_head=True,
        color="lightgrey",
        shape='full',
        zorder=-1
    )
    ax.add_patch(background_arrow)

    # 7) Restore y-limits
    # ax.set_ylim(ylim)

    # 8) Remove ticks and labels
    # ax.set_xticks([])
    # ax.set_yticks([])
    # ax.set_xticklabels([])
    # ax.set_yticklabels([])

    # plt.show()
    # 9) Save the figure
    plt.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main():
    """Plot domain figures for each .gb file in the current directory, sharing the same scale."""
    gb_files = glob.glob("*.gb")
    if not gb_files:
        print("No .gb files found.")
        return

    # 1) Find the largest genome length among all .gb files
    global_max_len = 0
    for gbf in gb_files:
        for record in SeqIO.parse(gbf, "genbank"):
            seq_len = len(record.seq)
            if seq_len > global_max_len:
                global_max_len = seq_len

    if global_max_len == 0:
        # if no sequences or empty files
        global_max_len = 1

    # 2) For each GenBank, parse features, plot them, save the figure
    for gbf in gb_files:
        basename = os.path.splitext(os.path.basename(gbf))[0]
        out_png = f"{basename}_domain.svg"

        # Parse the features from the file
        feats = parse_gbk_features(gbf)

        # Plot and save
        plot_domain_features(feats, global_max_len, out_png)
        print(f"Saved figure for {gbf} as {out_png}.")


if __name__ == "__main__":
    main()
