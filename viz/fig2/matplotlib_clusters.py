import matplotlib.pyplot as plt

target_names = [
    "ON107264", "MN270259", "MK448731", "NC_007581", "MW248466", "LC680885",
    "MN091626", "NC_031039", "NC_043027", "MW749003", "MZ422438", "MZ422438",
    "MZ422438", "MZ422438", "EU982300", "ON470608", "MK448228", "HQ906664",
    "KY695241", "OQ555808"
]

clusters = [
    1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 3
]

genome_lengths = [
    133818, 106491, 35793, 185683, 266637, 267055, 268748, 251042, 252197,
    252262, 223580, 223580, 223580, 223580, 22689, 46651, 46342, 0, 0, 62663
]

hosts = [
    "Bacillus", "Streptococcus", "Streptococcus", "Clostridium", "Staphylococcus",
    "Staphylococcus", "Staphylococcus", "Bacillus", "Bacillus", "Bacillus",
    "Listeria", "Listeria", "Listeria", "Listeria", "Pseudomonas", "Escherichia",
    "Klebsiella", "Wolbachia", "Wolbachia", "Mycobacterium"
]


# Define colors for each unique host
color_map = {
    "Escherichia": "#17558F",
    "Klebsiella": "#0583D2",
    "Pseudomonas": "#61B0B7",
    "Wolbachia": "#B8E3FF",
    
    "Bacillus": "#6E7C39",
    "Clostridium": "#A4DE03",
    "Enterococcus": "#76BA1C",
    "Listeria": "#4C9A2A",
    "Staphylococcus": "#ACDF87",
    "Streptococcus": "#65B556",
    
    "Mycobacterium": "orange",

    "Unspecified": "white",
}

# Count the number of clusters
num_clusters = len(set(clusters)) - 1  # Exclude the "-" cluster

# Calculate the maximum genome length (excluding zeros)
max_genome_length = max(length for length in genome_lengths if length != 0) if any(length != 0 for length in genome_lengths) else 1

# Create a figure and axis object
fig, ax = plt.subplots(figsize=(10, 8))

# Iterate through each cluster in reverse order to plot the 1 cluster first
for i in range(num_clusters, 0, -1):
    # Find the indices of targets in the current cluster
    indices = [j for j, cluster in enumerate(clusters) if cluster == i]

    # Plot a triangle for the current cluster
    if i != "-":
        cluster_start = len(target_names) - max(indices) - 0.1  # Adjusted to shift the triangles down
        cluster_end = len(target_names) - min(indices) + 0.1   # Adjusted to shift the triangles down
        cluster_center = (cluster_start + cluster_end) / 2
        ax.fill([0, 0.2, 0], [cluster_start-1, cluster_center-1, cluster_end-1], color='black')

# Define the size ranges and corresponding legend labels
size_ranges = [(10000, 50000), (50000, 100000), (100000, 150000), (150000, 200000), (200000, float('inf'))]
legend_labels = ['10-50', '50-100', '100-150', '150-200', '>200']

# Plot circles and squares with colors based on host
for i, (length, host) in enumerate(zip(genome_lengths, hosts)):
    if length == 0:
        # Plot square for dashes
        ax.plot(0.3, len(target_names) - i - 1, marker='s', markersize=8, color=color_map[host], markeredgecolor='black')
    else:
        # Plot circle with sizes based on genome size ranges
        for idx, size_range in enumerate(size_ranges):
            if size_range[0] <= length < size_range[1]:
                # Calculate the circle radius based on the index of the size range
                radius = idx + 2 +0.2
                ax.plot(0.3, len(target_names) - i - 1, marker='o', markersize=radius*3, color=color_map[host], markeredgecolor='black')
                break  # Exit loop once the size range is found

# Create legend for the size of circles and genome size ranges
legend1 = ax.legend(handles=[
    plt.Line2D([0], [0], marker='o', markersize=(idx + 2) * 3, color='black', linestyle='None', markeredgewidth=1, markeredgecolor='black') 
    for idx, _ in enumerate(size_ranges)],
    labels=legend_labels, title='Genome Size (kb)', loc='upper right', labelspacing=0.8)

ax.add_artist(legend1)

# Set x-axis tick positions and labels
ax.set_xticks([0])
ax.set_xticklabels([''])
ax.set_xlim(0, 0.5)  # Adjusted for circles and squares

# Set y-axis limits and labels
ax.set_ylim(-0.5, len(target_names) - 0.5)
# Reverse the order of target names
ax.set_yticks(range(len(target_names)))
ax.set_yticklabels(target_names[::-1])

# Create legend for host color information
legend2 = ax.legend(handles=[plt.Line2D([0], [0], marker='o', markersize=10, color=color, markeredgecolor='black', linestyle='None') for host, color in color_map.items()],
                    labels=list(color_map.keys()), title='Host', loc='lower right')
ax.add_artist(legend2)

# Remove spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

plt.tight_layout()

plt.savefig("clusters_plot1.svg", format="svg", dpi=600)

# If working in a notebook, show plot:
# plt.show()

print('Saved clusters_plot1.svg!')