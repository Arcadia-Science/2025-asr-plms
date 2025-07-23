#!/usr/bin/env python3
"""
Script to calculate both standard perplexity and pseudo perplexity scores
for protein sequences using ESM-2.
"""

import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor

import esm
import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate standard and pseudo perplexity scores for protein sequences using ESM-2"
        )
    )
    parser.add_argument("--input", required=True, help="Input FASTA file with protein sequences")
    parser.add_argument("--output", required=True, help="Output CSV file for results")
    parser.add_argument("--model", default="esm2_t33_650M_UR50D", help="ESM-2 model to use")
    parser.add_argument("--batch_size", type=int, default=4, help="Batch size for inference")
    parser.add_argument(
        "--mask_batch_size", type=int, default=32, help="Batch size for masked inference"
    )
    parser.add_argument(
        "--device",
        default="cuda:0" if torch.cuda.is_available() else "cpu",
        help="Device to run on",
    )
    parser.add_argument(
        "--no_pll", action="store_true", help="Skip pseudo perplexity calculation (much faster)"
    )
    parser.add_argument(
        "--max_seq_len",
        type=int,
        default=None,
        help="Maximum sequence length to process (for chunking)",
    )
    parser.add_argument(
        "--workers", type=int, default=1, help="Number of worker threads for data preprocessing"
    )
    return parser.parse_args()


def calculate_perplexity(model, batch_converter, sequences, device, batch_size):
    """Calculate standard perplexity for a list of protein sequences using ESM-2."""
    results = []

    # Process in batches
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i : i + batch_size]
        batch_labels, batch_strs, batch_tokens = batch_converter(batch)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            # Forward pass with no gradient computation
            results_batch = model(batch_tokens, repr_layers=[])
            logits = results_batch["logits"]

            # Calculate token likelihoods
            for seq_idx, (seq_label, seq) in enumerate(batch):
                seq_len = len(seq)
                token_probs = torch.log_softmax(logits[seq_idx, 1 : seq_len + 1], dim=-1)

                # Extract the probability of the target token at each position
                target_indices = batch_tokens[seq_idx, 1 : seq_len + 1]
                token_log_probs = torch.gather(token_probs, 1, target_indices.unsqueeze(1)).squeeze(
                    1
                )

                # Calculate mean log probability and perplexity
                mean_log_prob = token_log_probs.mean().item()
                seq_perplexity = np.exp(-mean_log_prob)

                # Store in results
                results.append(
                    {
                        "sequence_id": seq_label,
                        "sequence": seq,
                        "length": seq_len,
                        "mean_log_prob": mean_log_prob,
                        "perplexity": seq_perplexity,
                        "token_log_probs": token_log_probs.cpu().numpy().tolist(),
                    }
                )

    return results


def prepare_masked_batches(sequences, mask_token, max_batch_size=32):
    """Prepare batches of masked sequences efficiently."""
    all_batches = []

    for seq_label, seq in sequences:
        seq_len = len(seq)
        # Create masked sequences
        masked_batch = []

        for pos in range(seq_len):
            # Create a masked sequence
            masked_seq = list(seq)
            masked_seq[pos] = mask_token
            masked_seq = "".join(masked_seq)
            masked_batch.append((f"{seq_label}_pos{pos}", masked_seq, pos))

            # If we've reached the batch size, add to all_batches and reset
            if len(masked_batch) >= max_batch_size:
                all_batches.append(masked_batch)
                masked_batch = []

        # Add any remaining sequences
        if masked_batch:
            all_batches.append(masked_batch)

    return all_batches


def process_masked_batch(model, batch_converter, masked_batch, device, alphabet, original_seqs):
    """Process a batch of masked sequences and return position-wise log probabilities."""
    # Convert batch format for the model
    batch_for_model = [(label, seq) for label, seq, _ in masked_batch]
    _, _, tokens = batch_converter(batch_for_model)
    tokens = tokens.to(device)

    # Get the mask token index
    mask_idx = alphabet.tok_to_idx["<mask>"]

    results = {}

    with torch.no_grad():
        # Single forward pass for the entire batch
        output = model(tokens, repr_layers=[])
        logits = output["logits"]

        # Process each masked sequence
        for i, (label, _, pos) in enumerate(masked_batch):
            # Extract sequence ID from batch label
            seq_id = label.split("_pos")[0]

            # Find where the mask token is
            mask_positions = (tokens[i] == mask_idx).nonzero()

            if len(mask_positions) > 0:
                mask_pos = mask_positions[0].item()

                # Get the logits for the masked position
                pos_logits = logits[i, mask_pos]
                pos_log_probs = torch.log_softmax(pos_logits, dim=-1)

                # Get the original amino acid
                original_seq = next(seq for id, seq in original_seqs if id == seq_id)
                correct_aa = original_seq[pos]
                correct_aa_idx = alphabet.get_idx(correct_aa)

                # Get the log probability
                aa_log_prob = pos_log_probs[correct_aa_idx].item()

                # Store result
                if seq_id not in results:
                    results[seq_id] = {}
                results[seq_id][pos] = aa_log_prob

    return results


def calculate_pseudo_perplexity(
    model,
    batch_converter,
    alphabet,
    sequences,
    device,
    mask_batch_size,
    workers=1,
    max_seq_len=None,
):
    """
    Optimized version of pseudo perplexity calculation using batched processing and multi-threading.
    """
    # If max_seq_len is specified, filter or chunk sequences
    if max_seq_len is not None:
        filtered_sequences = []
        for seq_id, seq in sequences:
            if len(seq) <= max_seq_len:
                filtered_sequences.append((seq_id, seq))
            else:
                print(
                    f"Warning: Sequence {seq_id} exceeds max length ({len(seq)} > {max_seq_len}),\
                    processing first {max_seq_len} residues"
                )
                filtered_sequences.append((seq_id, seq[:max_seq_len]))
        sequences = filtered_sequences

    # Prepare all masked batches in advance using threading
    print("Preparing masked sequences...")
    start_time = time.time()

    # Function for ThreadPoolExecutor
    def chunk_prepare(seq_chunk):
        return prepare_masked_batches(seq_chunk, "<mask>", mask_batch_size)

    # Divide sequences among workers
    chunk_size = max(1, len(sequences) // workers)
    sequence_chunks = [sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)]

    all_masked_batches = []
    if workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            chunk_results = list(executor.map(chunk_prepare, sequence_chunks))
            for result in chunk_results:
                all_masked_batches.extend(result)
    else:
        all_masked_batches = prepare_masked_batches(sequences, "<mask>", mask_batch_size)

    print(f"Prepared {len(all_masked_batches)} batches in {time.time() - start_time:.2f} seconds")

    # Process all batches with tqdm progress bar
    position_scores = {}

    for batch in tqdm(all_masked_batches, desc="Calculating pseudo perplexity"):
        batch_results = process_masked_batch(
            model, batch_converter, batch, device, alphabet, sequences
        )

        # Update position_scores with batch results
        for seq_id, pos_scores in batch_results.items():
            if seq_id not in position_scores:
                position_scores[seq_id] = {}
            position_scores[seq_id].update(pos_scores)

    # Calculate the final metrics for each sequence
    results_with_pseudo = []

    for seq_id, seq in sequences:
        seq_len = len(seq)

        # Get all position scores for this sequence
        seq_position_scores = position_scores.get(seq_id, {})

        # Convert to array
        pll_token_probs = np.array([seq_position_scores.get(pos, 0.0) for pos in range(seq_len)])

        # Calculate the total pseudo log-likelihood and mean
        pll_score = pll_token_probs.sum()
        mean_pll = pll_score / seq_len

        # Calculate pseudo perplexity (exp(-mean_pll))
        pseudo_perplexity = np.exp(-mean_pll)

        # Find the standard perplexity result
        std_result = next((r for r in standard_results if r["sequence_id"] == seq_id), None)

        if std_result:
            # Add pseudo perplexity to standard result
            result = std_result.copy()
            result["pll"] = pll_score
            result["mean_pll"] = mean_pll
            result["pseudo_perplexity"] = pseudo_perplexity
            result["pll_token_probs"] = pll_token_probs.tolist()

            results_with_pseudo.append(result)

    return results_with_pseudo


def main(args):
    start_total = time.time()
    print(f"Loading ESM-2 model: {args.model}")
    model, alphabet = esm.pretrained.load_model_and_alphabet(args.model)

    model = model.half()  # ← ✅ Reduce memory with half precision
    model = model.to(args.device)
    model.eval()

    batch_converter = alphabet.get_batch_converter()

    print(f"Reading sequences from {args.input}")
    sequences = []
    for record in SeqIO.parse(args.input, "fasta"):
        sequences.append((record.id, str(record.seq)))

    print(f"Processing {len(sequences)} sequences")

    # First calculate standard perplexity
    print("Calculating standard perplexity...")
    start_time = time.time()
    global standard_results
    standard_results = calculate_perplexity(
        model, batch_converter, sequences, args.device, args.batch_size
    )
    print(f"Standard perplexity calculation completed in {time.time() - start_time:.2f} seconds")

    # Optionally calculate pseudo perplexity
    if not args.no_pll:
        print("Calculating pseudo perplexity (this will take longer)...")
        start_time = time.time()
        results_with_pseudo = calculate_pseudo_perplexity(
            model,
            batch_converter,
            alphabet,
            sequences,
            args.device,
            args.mask_batch_size,
            args.workers,
            args.max_seq_len,
        )
        print(f"Pseudo perplexity calculation completed in {time.time() - start_time:.2f} seconds")

        final_results = results_with_pseudo
    else:
        print("Skipping pseudo perplexity calculation")
        final_results = standard_results

    # Save results to CSV
    print("Saving results...")
    if args.no_pll:
        df = pd.DataFrame(
            [
                {
                    "sequence_id": r["sequence_id"],
                    "sequence": r["sequence"],
                    "length": r["length"],
                    "perplexity": r["perplexity"],
                    "mean_log_prob": r["mean_log_prob"],
                }
                for r in final_results
            ]
        )
    else:
        df = pd.DataFrame(
            [
                {
                    "sequence_id": r["sequence_id"],
                    "sequence": r["sequence"],
                    "length": r["length"],
                    "perplexity": r["perplexity"],
                    "mean_log_prob": r["mean_log_prob"],
                    "pseudo_perplexity": r.get("pseudo_perplexity"),
                    "mean_pll": r.get("mean_pll"),
                    "pll": r.get("pll"),
                }
                for r in final_results
            ]
        )

    df.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")

    # Save detailed token-level scores if we calculated them
    if not args.no_pll:
        detailed_dir = os.path.splitext(args.output)[0] + "_detailed"
        os.makedirs(detailed_dir, exist_ok=True)

        for result in final_results:
            seq_id = result["sequence_id"]
            seq = result["sequence"]

            df_detailed = pd.DataFrame(
                {
                    "position": list(range(1, len(seq) + 1)),
                    "amino_acid": list(seq),
                    "log_prob": result["token_log_probs"],
                    "pll_log_prob": result.get("pll_token_probs", [0] * len(seq)),
                }
            )

            df_detailed.to_csv(f"{detailed_dir}/{seq_id}.csv", index=False)

        print(f"Detailed token-level scores saved to {detailed_dir}/")

    # Show a summary comparison table
    print("\nSequence Comparison Summary:")
    print("-----------------------------")
    if args.no_pll:
        print(f"{'Sequence ID':<20} {'Length':<8} {'Perplexity':<12}")
        print(f"{'-' * 20} {'-' * 8} {'-' * 12}")
        for r in final_results:
            print(f"{r['sequence_id']:<20} {r['length']:<8} {r['perplexity']:<12.2f}")
    else:
        print(f"{'Sequence ID':<20} {'Length':<8} {'Perplexity':<12} {'Pseudo Perplexity':<18}")
        print(f"{'-' * 20} {'-' * 8} {'-' * 12} {'-' * 18}")
        for r in final_results:
            print(
                f"{r['sequence_id']:<20} {r['length']:<8} {r['perplexity']:<12.2f}\
                    {r.get('pseudo_perplexity', 0):<18.2f}"
            )

    print(f"\nTotal execution time: {time.time() - start_total:.2f} seconds")


if __name__ == "__main__":
    args = parse_args()
    main(args)
