import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from Bio import SeqIO
from collections import Counter
import json
import numpy as np
import matplotlib.cm as cm
import arcadia_pycolor as apc
import matplotlib.pyplot as plt
import seaborn as sns
import arcadia_pycolor as apc
apc.mpl.setup()




def calc_branch_length_to_root_node(tree, node_num):

    # Find the target internal node by confidence label (actually the node number)
    target_node = None
    node_num = int(node_num)
    for clade in tree.find_clades():
        if hasattr(clade, 'confidence') and clade.confidence == node_num:
            target_node = clade
            break

    if target_node is None:
        print(f"Node with confidence {node_num} not found.")
        return None

    # Calculate distance from target node to root
    distance_to_root = tree.distance(tree.root, target_node)
    return distance_to_root

def calc_branch_length_to_root_leaf(tree, leaf_name):
    return tree.distance(tree.root, leaf_name)

def get_node_labels_leaf_to_root(tree, leaf_name):
    """
    Retrieves the node labels along the path from a given leaf to the root.

    Parameters:
        tree (Phylo.BaseTree.Tree): The phylogenetic tree.
        leaf_name (str): The name of the leaf node.

    Returns:
        list: A list of node labels from leaf to root.
    """
    # Find the leaf
    leaf = tree.find_any(name=leaf_name)
    if leaf is None:
        print(f"Leaf {leaf_name} not found.")
        return []

    # Traverse from leaf to root
    labels = []
    path = tree.get_path(leaf)  # Returns list of clades from root to the leaf

    # Go in reverse to get from leaf to root
    for node in reversed(path):
        label = getattr(node, "confidence", None) or getattr(node, "name", None)
        labels.append(label)

    return ['node' + str(x) for x in labels[1:]]

def plot_evo_path(scores_df, tree, leaf_name, labels = False):
    these_nodes = get_node_labels_leaf_to_root(tree, leaf_name)

    # Subset the relevant data
    for_plot = scores_df[scores_df['sequence_id'].isin(these_nodes)]

    plt.figure(figsize=(16, 4))

    # # Scatter plot with color representing ML probability
    # scatter = plt.scatter(
    #     for_plot['bl_to_root'],
    #     for_plot['pseudo_perplexity'],
    #     c=for_plot['ML prob'],
    #     cmap=apc.gradients.viridis.to_mpl_cmap().reversed(),
    #     marker='o',
    #     s = 100
    # )
    scatter = plt.scatter(
    for_plot['bl_to_root'],
    for_plot['pseudo_perplexity'],
    c=for_plot['ML prob'],
    cmap=apc.gradients.viridis.to_mpl_cmap().reversed(),
    marker='o',
    s=100,
    vmin=0.77,  # Set to the minimum value you want to standardize on
    vmax=1   # Set to the maximum value you want to standardize on
)

    #Add labels for each point using the sequence_id
    if labels == True: 
        for i, row in for_plot.iterrows():
            plt.text(
                row['bl_to_root'],
                row['pseudo_perplexity'],
                row['sequence_id'],
                fontsize=9,
                ha='right',
                va='bottom'
            )

    # Add the leaf point in pink
    x = calc_branch_length_to_root_leaf(tree, leaf_name)
    y = scores_df[scores_df['sequence_id'] == leaf_name]['pseudo_perplexity'].to_list()[0]
    plt.scatter(x, y, color='#F28360', label='Leaf', zorder=10, s = 100)

    # Add colorbar
    plt.colorbar(scatter, label='ASR Mean Posterior Probability')

    # Labels and layout
    plt.xlabel('Distance to Tree Root (subs/site)')
    plt.ylabel('ESM2 Pseudo Perplexity')
    plt.title(leaf_name)
    plt.tight_layout()
    plt.show()

def plot_image_with_arrows(
    img_path,
    x_starts,
    x_ends,
    arrow_ys,
    labels,
    text_offsets=None,
    arrow_color='#F28360',
    arrow_linewidth=2,
    figsize=(12, 9),
    fontsize=14
):
    """
    Plot image with multiple left-pointing arrows and labels.

    Parameters:
    - img_path: str, path to image file.
    - x_starts: list of float, tail x positions (fraction of width).
    - x_ends: list of float, head x positions (fraction of width).
    - arrow_ys: list of float, y positions (fraction of height).
    - labels: list of str, text labels for each arrow.
    - text_offsets: list of float or None, horizontal offset for each label (fraction of width).
                    If None, defaults to 0.02 for all.
    - arrow_color: str, color of arrows and text.
    - arrow_linewidth: int, thickness of arrow lines.
    - figsize: tuple, figure size.
    - fontsize: int, font size for labels.
    """

    if not (len(x_starts) == len(x_ends) == len(arrow_ys) == len(labels)):
        raise ValueError("All input lists must have the same length.")

    if text_offsets is None:
        text_offsets = [0.02] * len(labels)
    elif len(text_offsets) != len(labels):
        raise ValueError("text_offsets must be None or have the same length as labels.")

    # Load image
    img = mpimg.imread(img_path)
    height, width = img.shape[0], img.shape[1]

    plt.figure(figsize=figsize)
    plt.imshow(img)
    plt.axis('off')

    # Plot each arrow + label
    for x_start, x_end, y_frac, label, txt_off in zip(x_starts, x_ends, arrow_ys, labels, text_offsets):
        tail_x = width * x_start
        head_x = width * x_end
        y_pos = height * y_frac

        # Draw arrow
        plt.annotate(
            '',
            xy=(head_x, y_pos),
            xytext=(tail_x, y_pos),
            arrowprops=dict(facecolor=arrow_color, edgecolor=arrow_color, arrowstyle='->', linewidth=arrow_linewidth)
        )

        # Draw label
        plt.text(
            tail_x + width * txt_off,
            y_pos,
            label,
            color=arrow_color,
            fontsize=fontsize,
            verticalalignment='center',
            horizontalalignment='left'
        )

    plt.show()

def plot_image_with_arrow_and_circles(
    img_path,
    x_start=0.85,
    x_end=0.75,
    arrow_y=0.05,
    circle_positions=[(0.5, 0.5)],
    circle_radius=10,
    circle_color='#F28360',
    filled=True,
    arrow_color='#F28360',
    arrow_linewidth=2,
    figsize=(8, 6)
):
    """
    Plot image with one arrow and multiple circles.

    Parameters:
    - img_path: str, path to image file.
    - x_start, x_end: float, tail and head x position of arrow (fraction of width).
    - arrow_y: float, y position of arrow (fraction of height).
    - circle_positions: list of (x, y) tuples (fractions of width/height).
    - circle_radius: int, radius of circles in pixels.
    - circle_color: str, color of circles.
    - filled: bool, whether the circles are filled (True) or outlined only (False).
    - arrow_color: str, color of arrow.
    - arrow_linewidth: int, arrow thickness.
    - figsize: tuple, figure size.
    """
    img = mpimg.imread(img_path)
    height, width = img.shape[0], img.shape[1]

    plt.figure(figsize=figsize)
    plt.imshow(img)
    plt.axis('off')

    # Arrow
    tail_x = width * x_start
    head_x = width * x_end
    y_pos = height * arrow_y

    plt.annotate(
        '',
        xy=(head_x, y_pos),
        xytext=(tail_x, y_pos),
        arrowprops=dict(facecolor=arrow_color, edgecolor=arrow_color, arrowstyle='->', linewidth=arrow_linewidth)
    )

    # Circles
    for x_frac, y_frac in circle_positions:
        circ_x = width * x_frac
        circ_y = height * y_frac
        circle = plt.Circle(
            (circ_x, circ_y),
            radius=circle_radius,
            facecolor=circle_color if filled else 'none',
            edgecolor=circle_color,
            linewidth=2
        )
        plt.gca().add_patch(circle)

    plt.show()



def get_descendant_leaves(tree, node_num):
    node_num = int(node_num.replace('node',''))

    # Find the target internal node by confidence label
    target_node = None
    for clade in tree.find_clades():
        if getattr(clade, 'confidence', None) == node_num:
            target_node = clade
            break

    if target_node is None:
        print(f"Node with confidence {node_num} not found.")
        return []

    # Get all terminal (leaf) names under this node
    return [leaf.name for leaf in target_node.get_terminals()]

def load_name_conversion(name_conv_dict):
    name_dict = {}
    rev_dict = {}
    with open(name_conv_dict) as file:
        for line in file:
            key, val = line.rstrip().split('\t')
            name_dict[val] = key
            rev_dict[key] = val
    return name_dict, rev_dict

def extract_sequences(leaf_names, alignment_file, name_dict, rev_dict):
    recoded_names = [name_dict[name] for name in leaf_names]
    output_dict = {}
    for entry in SeqIO.parse(alignment_file, 'fasta'):
        if entry.id in recoded_names:
            output_dict[rev_dict[entry.id]] = str(entry.seq)
    return list(output_dict.values())

def make_consensus(sequences, gap_profile):
    consensus = ''
    for i, is_gap in enumerate(gap_profile):
        if not is_gap:
            residues = [seq[i] for seq in sequences]
            counts = Counter(residues).most_common()
            if counts[0][0] == '-' and len(counts) > 1:
                consensus += counts[1][0]
            else:
                consensus += counts[0][0]
    return consensus

def generate_node_consensus(tree, node_num, alignment_file, name_conv_dict, gap_dict_json):
    leaf_names = get_descendant_leaves(tree, node_num)
    name_dict, rev_dict = load_name_conversion(name_conv_dict)
    sequences = extract_sequences(leaf_names, alignment_file, name_dict, rev_dict)
    with open(gap_dict_json, 'r') as f:
        gap_profile = json.load(f)[str(node_num.replace('node',''))]
    return make_consensus(sequences, gap_profile)


def plot_evo_path_quiver(scores_df, tree, leaf_name):
    these_nodes = get_node_labels_leaf_to_root(tree, leaf_name)

    # Subset the relevant data
    for_plot = scores_df[scores_df['sequence_id'].isin(these_nodes)].copy()

    # Get required columns
    x_vals = for_plot['bl_to_root'].values
    y_old = for_plot['pseudo_perplexity'].values
    y_new = for_plot['consensus_pppl'].values
    dy = y_new - y_old
    dx = np.zeros_like(dy)

    # Split indices by direction of change
    up_indices = dy > 0
    down_indices = dy < 0

    plt.figure(figsize=(18, 4))

    # Plot original points
    scatter = plt.scatter(
        x_vals,
        y_old,
        color = 'black',
        marker='o',
        label='ML ancestor',
        s = 80
    )

    # Upward arrows in green
    plt.quiver(
        x_vals[up_indices], y_old[up_indices],
        dx[up_indices], dy[up_indices],
        angles='xy', scale_units='xy', scale=1,
        color='#73B5E3', width=0.003, headwidth=3, headlength=5
    )

    # Downward arrows in red
    plt.quiver(
        x_vals[down_indices], y_old[down_indices],
        dx[down_indices], dy[down_indices],
        angles='xy', scale_units='xy', scale=1,
        color='#F7B846', width=0.003, headwidth=3, headlength=5
    )

    # Plot new (consensus) values
    plt.scatter(
        x_vals,
        y_new,
        color='black',
        marker='X',
        label='Consensus ancestor',
        alpha=0.8, 
        s = 80
    )

    # Highlight the leaf
    x_leaf = calc_branch_length_to_root_leaf(tree, leaf_name)
    y_leaf = scores_df[scores_df['sequence_id'] == leaf_name]['pseudo_perplexity'].values[0]
    plt.scatter(x_leaf, y_leaf, color='#F28360', label='Leaf', zorder=10, s = 80)

    # Labels and layout
    plt.xlabel('Distance to Tree Root (subs/site)')
    plt.ylabel('ESM2 Pseudo Perplexity')
    plt.title(leaf_name)
    plt.scatter([], [], color='#73B5E3', label='ML ancestor lower pppl', marker='^')
    plt.scatter([], [], color='#F7B846', label='Consensus ancestor lower pppl', marker='v')
    plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    plt.show()

def plot_multiple_evo_lines(score_dfs_leaves_labels, tree, normalize=True):
    """
    Parameters:
        score_dfs_leaves_labels: list of tuples (scores_df, leaf_name, label_name)
        tree: Phylogenetic tree used to compute distances
        normalize: whether to normalize pseudo perplexity to the leaf value
    """
    plt.figure(figsize=(18, 3))
    
    n_lines = len(score_dfs_leaves_labels)
    colors = cm.rainbow_r(np.linspace(0, 1.0, n_lines))  # rainbow colormap
    colors = ['#C85152', '#FFB984', '#97CD78', '#73B5E3', '#7A77AB']

    for i, (scores_df, leaf_name, label_name) in enumerate(score_dfs_leaves_labels):
        color = colors[i]

        these_nodes = get_node_labels_leaf_to_root(tree, leaf_name)
        for_plot = scores_df[scores_df['sequence_id'].isin(these_nodes)].copy()
        for_plot.sort_values('bl_to_root', inplace=True)

        # Get leaf value
        leaf_value = scores_df[scores_df['sequence_id'] == leaf_name]['pseudo_perplexity'].to_list()[0]

        # Normalize if requested
        if normalize:
            for_plot['plot_value'] = for_plot['pseudo_perplexity'] / leaf_value
            y_leaf = 1.0
        else:
            for_plot['plot_value'] = for_plot['pseudo_perplexity']
            y_leaf = leaf_value

        x_leaf = calc_branch_length_to_root_leaf(tree, leaf_name)

        # Main line
        plt.plot(
            for_plot['bl_to_root'],
            for_plot['plot_value'],
            color=color,
            label=label_name
        )

        # Invisible marker at leaf
        plt.scatter(x_leaf, y_leaf, color=color, s=0.1, zorder=10)

        # Line to closest point
        closest_idx = (for_plot['bl_to_root'] - x_leaf).abs().idxmin()
        x_closest = for_plot.loc[closest_idx, 'bl_to_root']
        y_closest = for_plot.loc[closest_idx, 'plot_value']

        plt.plot(
            [x_leaf, x_closest],
            [y_leaf, y_closest],
            color=color
        )

    # Labels and layout
    plt.xlabel('Distance to Tree Root (subs/site)')
    y_label = 'Normalized\nPseudo Perplexity\n(relative to leaf)' if normalize else 'ESM2 Pseudo Perplexity'
    plt.ylabel(y_label)
    plt.title(leaf_name)
    plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    plt.tight_layout()
    plt.show()

def violin_plot(orig_input_df, bins, bin_labels):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    input_df = orig_input_df.copy()

    custom_colors = {
    '< 0.8': '#fffdbd',
    '0.80–0.85': '#9dd07c',
    '0.85–0.90': '#649bb0',
    '0.90–0.95': '#446b9f',
    '0.95–1.00': '#313f65',
    'extant': '#F28360'
}

    # Bin the ML prob values using custom labels
    input_df['ML_prob_bin'] = pd.cut(
        input_df['ML prob'],
        bins=bins,
        labels=bin_labels,
        include_lowest=True
    )

    # Drop rows with missing bin assignments or pseudo_perplexity
    input_df = input_df.dropna(subset=['ML_prob_bin', 'pseudo_perplexity'])

    # Count samples per bin
    bin_counts = input_df['ML_prob_bin'].value_counts().sort_index()

    # Calculate medians for each bin
    group_stats = input_df.groupby('ML_prob_bin', observed=False)['pseudo_perplexity'].agg('median').values

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.violinplot(
        x='ML_prob_bin',
        y='pseudo_perplexity',
        hue='ML_prob_bin',
        data=input_df,
        ax=ax,
        inner=None,
        cut=0,
        palette=custom_colors,
        legend=False
    )

    # Add horizontal line for medians
    for i, stat in enumerate(group_stats):
        ax.plot([i - 0.2, i + 0.2], [stat, stat], color='black', linewidth=2)

    # Labeling
    ax.set_xlabel('ASR mean posterior probability bin')
    ax.set_ylabel('ESM2 Pseudo Perplexity')
    plt.suptitle("")

    # Add group sizes 
    y_max = ax.get_ylim()[1]
    y_min = ax.get_ylim()[0]
    offset = 0.03 * (y_max - y_min)

    for i, (__, count) in enumerate(bin_counts.items()):
        ax.text(i, y_max - offset, f'n={count}', ha='center', va='top', fontsize=12)

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def kruskal_wallis_with_significant_posthoc(data_df, bins, bin_labels):
    from scipy.stats import kruskal
    import scikit_posthocs as sp
    import pandas as pd

    # Bin the ML posterior probabilities
    data_df['ML_prob_bin'] = pd.cut(
        data_df['ML prob'],
        bins=bins,
        labels=bin_labels,
        include_lowest=True
    )

    # Drop missing values
    data_df = data_df.dropna(subset=['ML_prob_bin', 'pseudo_perplexity'])

    # Group pseudo_perplexity values by bin
    groups = [group['pseudo_perplexity'].values for _, group in data_df.groupby('ML_prob_bin', observed=False)]

    # Kruskal-Wallis test
    H, p = kruskal(*groups)
    print(f"Kruskal-Wallis test result: H = {H:.3f}, p = {p:.3e}\n")

    # Dunn's post hoc test with Bonferroni correction
    posthoc = sp.posthoc_dunn(data_df, val_col='pseudo_perplexity', group_col='ML_prob_bin', p_adjust='bonferroni')

    # Find significant pairs (p < 0.05)
    sig_pairs = []
    bins_list = list(bin_labels)
    for i, g1 in enumerate(bins_list):
        for j, g2 in enumerate(bins_list):
            if i < j:  # upper triangle only
                pval = posthoc.loc[g1, g2]
                if pval < 0.05:
                    sig_pairs.append((g1, g2, pval))

    # Print significant comparisons in table format
    if sig_pairs:
        print("Significant Dunn post hoc comparisons (adjusted p < 0.05):")
        print(f"{'Group 1':<12} {'Group 2':<12} {'Adj. p-value':<12}")
        print("-" * 38)
        for g1, g2, pval in sig_pairs:
            print(f"{g1:<12} {g2:<12} {pval:<12.3e}")
    else:
        print("No significant pairwise differences found in post hoc test.")