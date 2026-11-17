# -*- coding: utf-8 -*-
"""
CODON USAGE-BASED EXPRESSION LEVEL ESTIMATOR


import math
import matplotlib.pyplot as plt


# Standard genetic code (DNA codons)
STANDARD_GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Built once: amino acid -> list of synonymous codons (stop codons excluded)
AA_TO_CODONS = {}
for _codon, _aa in STANDARD_GENETIC_CODE.items():
    if _aa != '*':
        AA_TO_CODONS.setdefault(_aa, []).append(_codon)

# Codons with no synonymous alternatives (excluded from CAI and RSCU plots)
SINGLE_CODON_AMINO_ACIDS = ('ATG', 'TGG')


def clean_and_validate_sequence(sequence, require_complete_codons=False):
    """Normalize a DNA/mRNA sequence and validate its nucleotide characters."""
    if not isinstance(sequence, str):
        raise TypeError("Sequence must be provided as a string.")

    # Whitespace removed, uppercased, and U converted to T so mRNA works too.
    clean_seq = ''.join(sequence.split()).upper().replace('U', 'T')

    if not clean_seq:
        raise ValueError("Empty sequence provided.")

    invalid_chars = set(clean_seq) - {'A', 'C', 'G', 'T'}
    if invalid_chars:
        raise ValueError(
            f"Sequence contains invalid nucleotides: {sorted(invalid_chars)}"
        )

    remainder = len(clean_seq) % 3
    if remainder:
        message = (
            f"Sequence length ({len(clean_seq)}) is not divisible by 3. "
            f"The final {remainder} nucleotide(s) do not form a complete codon."
        )
        if require_complete_codons:
            raise ValueError(message)
        print(f"[Warning] {message} They will be ignored.")

    return clean_seq


def split_into_codons(sequence):
    """
    Return in-frame codons, stopping at the first stop codon.

    Codons after a premature stop are dropped (translation would not
    reach them). Warnings are printed for a missing start codon, a
    premature stop, or a missing terminal stop.
    """
    usable_length = len(sequence) - (len(sequence) % 3)
    codons = [sequence[i:i + 3] for i in range(0, usable_length, 3)]

    if codons and codons[0] != 'ATG':
        print(f"[Warning] Sequence does not start with ATG (starts with {codons[0]}).")

    for index, codon in enumerate(codons):
        if STANDARD_GENETIC_CODE[codon] == '*':
            remaining = len(codons) - index - 1
            if remaining:
                print(
                    f"[Warning] Premature stop codon {codon} at codon "
                    f"{index + 1}; {remaining} codon(s) after it ignored."
                )
            return codons[:index + 1]

    print("[Warning] No terminal stop codon found.")
    return codons


class CodonUsageAnalyzer:
    """
    Performs codon counting and RSCU analysis on a coding DNA sequence.

    If reference_weights are supplied, relative adaptiveness (wi) is taken
    from the reference dataset (None where the reference has no data for
    that amino acid). Otherwise, it is calculated from the input sequence
    itself and should be treated as descriptive, not as a reference-based
    expression measure.
    """

    STANDARD_GENETIC_CODE = STANDARD_GENETIC_CODE
    aa_to_codons = AA_TO_CODONS

    def validate_sequence(self, sequence):
        return clean_and_validate_sequence(sequence)

    def analyze(self, sequence, reference_weights=None):
        validated_seq = self.validate_sequence(sequence)
        codons = split_into_codons(validated_seq)

        codon_counts = {
            codon: 0
            for codon, amino_acid in self.STANDARD_GENETIC_CODE.items()
            if amino_acid != '*'
        }
        stop_codon_count = 0
        stop_positions = []

        for position, codon in enumerate(codons, start=1):
            if self.STANDARD_GENETIC_CODE[codon] == '*':
                stop_codon_count += 1
                stop_positions.append(position)
            else:
                codon_counts[codon] += 1

        aa_counts = {
            amino_acid: sum(codon_counts[c] for c in synonymous_codons)
            for amino_acid, synonymous_codons in self.aa_to_codons.items()
        }

        rscu_scores = {}
        sequence_adaptiveness = {}

        for amino_acid, synonymous_codons in self.aa_to_codons.items():
            total = aa_counts[amino_acid]
            n_syn = len(synonymous_codons)
            maximum_count = max(codon_counts[c] for c in synonymous_codons)

            for codon in synonymous_codons:
                count = codon_counts[codon]
                rscu_scores[codon] = count / (total / n_syn) if total > 0 else 0.0
                sequence_adaptiveness[codon] = (
                    count / maximum_count if maximum_count > 0 else 0.0
                )

        if reference_weights is not None:
            relative_adaptiveness = {
                codon: reference_weights.get(codon) for codon in codon_counts
            }
            adaptiveness_source = "reference dataset"
        else:
            relative_adaptiveness = sequence_adaptiveness
            adaptiveness_source = "input sequence (descriptive only)"

        return {
            'sequence_length': len(validated_seq),
            'complete_nucleotides_analyzed': len(codons) * 3,
            'codons_analyzed': len(codons),
            'stop_codon_count': stop_codon_count,
            'stop_codon_positions': stop_positions,
            'codon_counts': codon_counts,
            'aa_counts': aa_counts,
            'rscu_scores': rscu_scores,
            'relative_adaptiveness': relative_adaptiveness,
            'adaptiveness_source': adaptiveness_source
        }

    def display_results(self, results):
        print("=" * 78)
        print("                 CODON USAGE ANALYSIS RESULTS")
        print("=" * 78)
        print(f"Sequence length (nt)       : {results['sequence_length']}")
        print(
            "Complete nucleotides used : "
            f"{results['complete_nucleotides_analyzed']}"
        )
        print(f"Complete codons analyzed   : {results['codons_analyzed']}")
        print(f"Stop codons identified     : {results['stop_codon_count']}")
        if results['stop_codon_positions']:
            print(
                "Stop codon positions       : "
                f"{results['stop_codon_positions']}"
            )
        print(f"Relative adaptiveness source: {results['adaptiveness_source']}")
        print("-" * 78)
        print(
            f"{'Codon':<8}{'AA':<6}{'Count':<9}"
            f"{'RSCU':<12}{'Relative adaptiveness (wi)':<27}"
        )
        print("-" * 78)

        sorted_codons = sorted(
            results['codon_counts'],
            key=lambda c: (self.STANDARD_GENETIC_CODE[c], c)
        )

        for codon in sorted_codons:
            count = results['codon_counts'][codon]
            if count == 0:
                continue
            rscu = results['rscu_scores'][codon]
            weight = results['relative_adaptiveness'][codon]
            weight_text = "n/a (no reference)" if weight is None else f"{weight:.3f}"
            amino_acid = self.STANDARD_GENETIC_CODE[codon]
            print(
                f"{codon:<8}{amino_acid:<6}{count:<9}"
                f"{rscu:<12.3f}{weight_text:<27}"
            )
        print("=" * 78)

    def plot_rscu(self, results, save_path=None, show=True):
        """
        Plot RSCU for observed codons with synonymous alternatives
        (ATG and TGG are skipped: their RSCU is always 1).
        """
        active_codons = {
            codon: value
            for codon, value in results['rscu_scores'].items()
            if results['codon_counts'][codon] > 0
            and codon not in SINGLE_CODON_AMINO_ACIDS
        }

        if not active_codons:
            print("No codons with synonymous alternatives are available to plot.")
            return

        sorted_items = sorted(
            active_codons.items(),
            key=lambda item: (self.STANDARD_GENETIC_CODE[item[0]], item[0])
        )
        codons = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]
        labels = [
            f"{codon}\n({self.STANDARD_GENETIC_CODE[codon]})"
            for codon in codons
        ]

        plt.figure(figsize=(12, 5))
        plt.bar(labels, values, edgecolor='black', alpha=0.8)
        plt.axhline(
            y=1.0,
            linestyle='--',
            label='RSCU = 1 (equal synonymous usage)'
        )
        plt.xlabel('Codon (amino acid)')
        plt.ylabel('RSCU score')
        plt.title('Relative Synonymous Codon Usage of Input Sequence')
        plt.xticks(rotation=45)
        plt.legend()
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=200)
        if show:
            plt.show()
        plt.close()


class ExpressionEstimator:
    """
    Calculates CAI using codon frequencies from a reference set of highly
    expressed genes.

    Reference frequencies should ideally come from the same organism or
    expression host as the candidate gene. Amino acids with no reference
    data at all get weight None and are skipped in CAI (rather than being
    treated as perfectly adapted).
    """

    STANDARD_GENETIC_CODE = STANDARD_GENETIC_CODE
    aa_to_codons = AA_TO_CODONS

    def __init__(self, reference_frequencies, pseudocount=1.0):
        if not isinstance(reference_frequencies, dict):
            raise TypeError("Reference frequencies must be provided as a dictionary.")

        if not isinstance(pseudocount, (int, float)) or pseudocount <= 0:
            raise ValueError("Pseudocount must be a positive number.")

        self.pseudocount = float(pseudocount)
        self.reference_frequencies = {}

        for codon, count in reference_frequencies.items():
            codon = codon.strip().upper().replace('U', 'T')

            if codon not in self.STANDARD_GENETIC_CODE:
                raise ValueError(f"Invalid reference codon: {codon}")

            if not isinstance(count, (int, float)) or count < 0:
                raise ValueError(
                    f"Reference count for {codon} must be a non-negative number."
                )

            if self.STANDARD_GENETIC_CODE[codon] == '*':
                continue  # stop codons are not used in CAI weights

            self.reference_frequencies[codon] = float(count)

        self.weights = self._calculate_reference_weights()

    def _calculate_reference_weights(self):
        """
        wi = frequency of a codon / maximum frequency among synonymous codons.

        A pseudocount is added to each synonymous codon so an unobserved
        codon does not receive a zero weight, but ONLY for amino acids that
        have at least one observed codon in the reference. Otherwise the
        weights are None (no reference information).
        """
        weights = {}

        for amino_acid, synonymous_codons in self.aa_to_codons.items():
            raw_counts = [
                self.reference_frequencies.get(c, 0.0) for c in synonymous_codons
            ]

            if sum(raw_counts) == 0:
                for codon in synonymous_codons:
                    weights[codon] = None
                continue

            adjusted = {
                codon: self.reference_frequencies.get(codon, 0.0) + self.pseudocount
                for codon in synonymous_codons
            }
            maximum_count = max(adjusted.values())

            for codon in synonymous_codons:
                weights[codon] = adjusted[codon] / maximum_count

        return weights

    def validate_sequence(self, sequence):
        return clean_and_validate_sequence(sequence)

    def calculate_cai(self, sequence, return_details=False):
        """
        CAI = exp(mean(ln(wi))).

        Excluded from the calculation: stop codons, ATG, TGG (no synonymous
        alternatives), and codons whose amino acid has no reference data.
        Returns None if no codons can be evaluated (CAI is undefined).
        """
        validated_seq = self.validate_sequence(sequence)
        codons = split_into_codons(validated_seq)

        log_sum = 0.0
        evaluated_codons = []
        excluded_codons = []
        skipped_no_reference = []

        for codon in codons:
            amino_acid = self.STANDARD_GENETIC_CODE[codon]

            if amino_acid == '*' or codon in SINGLE_CODON_AMINO_ACIDS:
                excluded_codons.append(codon)
                continue

            weight = self.weights[codon]

            if weight is None:
                skipped_no_reference.append(codon)
                continue

            if weight <= 0:
                raise ValueError(
                    f"Non-positive adaptiveness weight encountered for {codon}."
                )

            log_sum += math.log(weight)
            evaluated_codons.append(codon)

        if evaluated_codons:
            cai = math.exp(log_sum / len(evaluated_codons))
        else:
            cai = None

        if return_details:
            return {
                'cai': cai,
                'evaluated_codon_count': len(evaluated_codons),
                'excluded_codons': excluded_codons,
                'skipped_no_reference': skipped_no_reference,
                'evaluated_codons': evaluated_codons
            }

        return cai

    def classify_expression(self, cai, high_threshold=0.8, moderate_threshold=0.5):
        """
        Approximate project categories, not universal biological cut-offs.
        Adjust thresholds only with a suitable organism-specific dataset.
        """
        if cai is None:
            return "Not enough synonymous codons with reference data"

        if not (0 <= cai <= 1):
            raise ValueError("CAI must be between 0 and 1.")

        if not (0 <= moderate_threshold <= high_threshold <= 1):
            raise ValueError(
                "Thresholds must satisfy 0 <= moderate <= high <= 1."
            )

        if cai >= high_threshold:
            return "High codon adaptation"
        if cai >= moderate_threshold:
            return "Moderate codon adaptation"
        return "Low codon adaptation"

    def evaluate_gene(self, sequence):
        details = self.calculate_cai(sequence, return_details=True)
        return {
            'cai': details['cai'],
            'classification': self.classify_expression(details['cai']),
            'evaluated_codon_count': details['evaluated_codon_count'],
            'excluded_codons': details['excluded_codons'],
            'skipped_no_reference': details['skipped_no_reference']
        }

    def display_reference_weights(self):
        """Print all reference-derived relative adaptiveness weights."""
        print("\nREFERENCE-BASED RELATIVE ADAPTIVENESS WEIGHTS")
        print("-" * 45)
        for codon in sorted(self.weights):
            amino_acid = self.STANDARD_GENETIC_CODE[codon]
            weight = self.weights[codon]
            text = "n/a (no reference data)" if weight is None else f"{weight:.3f}"
            print(f"{codon} ({amino_acid}): {text}")


# ---------------------------------------------------------------------
# TESTS / EXAMPLE RUN
# ---------------------------------------------------------------------

if __name__ == "__main__":

    # Illustrative reference counts only.
    # Replace these with counts from highly expressed genes of your
    # organism or expression host for meaningful biological interpretation.
    reference_frequencies = {
        'ATG': 100, 'TGG': 80,
        'GCC': 120, 'GCT': 20, 'GCA': 10, 'GCG': 5,
        'TTC': 90, 'TTT': 10,
        'GAC': 85, 'GAT': 15,
        'GAA': 110, 'GAG': 15
    }

    # Candidate coding DNA sequences (illustrative examples)
    valid_sequence = "ATGGCCTTCGACGAAGCTGCCGCAATGGCCTTCGACGAATAA"
    premature_stop_sequence = "ATGGCCTTCGACGAATAAGGC"        # TAA then GGC
    incomplete_codon_sequence = valid_sequence + "GG"        # 2 stray nt
    high_adaptation_candidate = "ATGGCCGACGAAATGGCCTTCGACGAATAA"
    low_adaptation_candidate = "ATGGCGGCGGCGGCGGCGGCGTAA"
    no_reference_candidate = "ATGCTGCTGTCCGGCTAA"            # Leu/Ser/Gly: no reference
    no_synonymous_candidate = "ATGTGGTAA"                    # only Met, Trp, stop

    analyzer = CodonUsageAnalyzer()
    estimator = ExpressionEstimator(
        reference_frequencies=reference_frequencies,
        pseudocount=1
    )

    print("\n--- REFERENCE DATASET / WEIGHTS ---")
    estimator.display_reference_weights()

    print("\n--- CODON USAGE ANALYSIS ---")
    results = analyzer.analyze(
        valid_sequence,
        reference_weights=estimator.weights
    )
    analyzer.display_results(results)
    analyzer.plot_rscu(results, save_path="rscu_plot.png", show=False)

    print("\n--- EDGE CASE: PREMATURE STOP ---")
    analyzer.display_results(
        analyzer.analyze(premature_stop_sequence,
                         reference_weights=estimator.weights)
    )

    print("\n--- EDGE CASE: INCOMPLETE FINAL CODON ---")
    analyzer.display_results(
        analyzer.analyze(incomplete_codon_sequence,
                         reference_weights=estimator.weights)
    )

    print("\n--- EXPRESSION / CAI ESTIMATION ---")
    for label, candidate in [
        ("High-adaptation example", high_adaptation_candidate),
        ("Low-adaptation example", low_adaptation_candidate),
        ("No reference data for its amino acids", no_reference_candidate),
        ("Only Met/Trp/stop", no_synonymous_candidate),
    ]:
        evaluation = estimator.evaluate_gene(candidate)
        cai = evaluation['cai']
        print(f"\n{label}")
        print(f"Sequence: {candidate}")
        print(f"CAI: {'undefined' if cai is None else f'{cai:.4f}'}")
        print(f"Category: {evaluation['classification']}")
        print(f"Codons included in CAI: {evaluation['evaluated_codon_count']}")
        print(f"Excluded codons: {evaluation['excluded_codons']}")
        print(f"Skipped (no reference data): {evaluation['skipped_no_reference']}")

    print("\nNote: These categories describe codon adaptation,")
    print("not experimentally measured gene expression.")
