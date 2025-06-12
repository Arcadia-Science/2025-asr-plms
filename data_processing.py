import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from Bio import SeqIO
from collections import Counter
import json
import numpy as np



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

    plt.figure(figsize=(12, 3))

    # Scatter plot with color representing ML probability
    scatter = plt.scatter(
        for_plot['bl_to_root'],
        for_plot['pseudo_perplexity'],
        c=for_plot['ML prob'],
        cmap='viridis_r',
        marker='o'
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
    plt.scatter(x, y, color='deeppink', label='Leaf', zorder=10)

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
    arrow_color='deeppink',
    arrow_linewidth=2,
    figsize=(8, 6),
    fontsize=10
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
    circle_color='deeppink',
    filled=True,
    arrow_color='deeppink',
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

    plt.figure(figsize=(12, 3))

    # Plot original points
    scatter = plt.scatter(
        x_vals,
        y_old,
        color = 'black',
        marker='o',
        label='ML'
    )

    # Upward arrows in green
    plt.quiver(
        x_vals[up_indices], y_old[up_indices],
        dx[up_indices], dy[up_indices],
        angles='xy', scale_units='xy', scale=1,
        color='mediumseagreen', width=0.003, headwidth=3, headlength=5
    )

    # Downward arrows in red
    plt.quiver(
        x_vals[down_indices], y_old[down_indices],
        dx[down_indices], dy[down_indices],
        angles='xy', scale_units='xy', scale=1,
        color='tomato', width=0.003, headwidth=3, headlength=5
    )

    # Plot new (consensus) values
    plt.scatter(
        x_vals,
        y_new,
        color='black',
        marker='X',
        label='Consensus',
        alpha=0.8
    )

    # Highlight the leaf
    x_leaf = calc_branch_length_to_root_leaf(tree, leaf_name)
    y_leaf = scores_df[scores_df['sequence_id'] == leaf_name]['pseudo_perplexity'].values[0]
    plt.scatter(x_leaf, y_leaf, color='deeppink', label='Leaf', zorder=10)

    # Labels and layout
    plt.xlabel('Distance to Tree Root (subs/site)')
    plt.ylabel('ESM2 Pseudo Perplexity')
    plt.title(leaf_name)
    plt.legend()
    plt.tight_layout()
    plt.show()

