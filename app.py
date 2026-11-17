
import json
import streamlit as st
import pandas as pd

from codon_usage_based_expression_level_estimator import (
    CodonUsageAnalyzer,
    ExpressionEstimator,
    STANDARD_GENETIC_CODE
)

# --------------------------------------------------
# PAGE CONFIGURATION
# --------------------------------------------------

st.set_page_config(
    page_title="GeneExpress",
    page_icon="🧬",
    layout="wide"
)

# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.title("🧬 GeneExpress")
st.subheader("Codon Usage-Based Expression Level Estimator")

st.write(
    "Analyze coding DNA sequences using codon frequency, "
    "RSCU, relative adaptiveness and Codon Adaptation Index (CAI)."
)

st.divider()

# --------------------------------------------------
# REFERENCE FREQUENCIES
# --------------------------------------------------

default_reference = {
    "ATG": 100,
    "TGG": 80,
    "GCC": 120,
    "GCT": 20,
    "GCA": 10,
    "GCG": 5,
    "TTC": 90,
    "TTT": 10,
    "GAC": 85,
    "GAT": 15,
    "GAA": 110,
    "GAG": 15
}

st.sidebar.header("Reference Dataset")

st.sidebar.caption(
    "Enter codon counts from highly expressed genes "
    "of your selected organism or expression host."
)

reference_json = st.sidebar.text_area(
    "Reference codon frequencies (JSON)",
    value=json.dumps(default_reference, indent=4),
    height=350
)

st.sidebar.info(
    "The pre-filled values are illustrative examples, "
    "not validated organism-specific reference data."
)

# --------------------------------------------------
# DNA SEQUENCE INPUT
# --------------------------------------------------

st.header("1. Enter your DNA sequence")

sequence = st.text_area(
    "Paste your coding DNA sequence (CDS)",
    height=180,
    placeholder="Example: ATGGCCTTCGACGAAGCTGCCGCAATGGCCTAA"
)

st.caption(
    "You can enter DNA (A, T, G, C) or RNA (A, U, G, C). "
    "Spaces and line breaks are accepted."
)

analyze_button = st.button(
    "🔬 Analyze Sequence",
    type="primary",
    use_container_width=True
)

# --------------------------------------------------
# ANALYSIS
# --------------------------------------------------

if analyze_button:

    if not sequence.strip():
        st.error("Please enter a DNA sequence first.")

    else:

        try:
            # Load reference dataset
            reference_frequencies = json.loads(reference_json)

            if not isinstance(reference_frequencies, dict):
                raise ValueError(
                    "Reference data must be a JSON dictionary."
                )

            # Initialize analysis modules
            analyzer = CodonUsageAnalyzer()

            estimator = ExpressionEstimator(
                reference_frequencies=reference_frequencies,
                pseudocount=1
            )

            # Perform codon usage analysis
            results = analyzer.analyze(
                sequence,
                reference_weights=estimator.weights
            )

            # Perform CAI estimation
            evaluation = estimator.evaluate_gene(sequence)

            # --------------------------------------------------
            # SUMMARY
            # --------------------------------------------------

            st.divider()
            st.header("2. Analysis Summary")

            col1, col2, col3, col4 = st.columns(4)

            col1.metric(
                "Sequence Length",
                f"{results['sequence_length']} nt"
            )

            col2.metric(
                "Codons Analyzed",
                results["codons_analyzed"]
            )

            col3.metric(
                "Stop Codons",
                results["stop_codon_count"]
            )

            cai = evaluation["cai"]

            col4.metric(
                "CAI Score",
                "Undefined" if cai is None else f"{cai:.4f}"
            )

            st.write(
                "**Adaptiveness source:**",
                results["adaptiveness_source"]
            )

            if results["stop_codon_positions"]:
                st.warning(
                    "Stop codon positions: "
                    + str(results["stop_codon_positions"])
                )

            # --------------------------------------------------
            # CODON USAGE TABLE
            # --------------------------------------------------

            st.divider()
            st.header("3. Codon Usage Analysis")

            codon_rows = []

            for codon, count in results["codon_counts"].items():

                if count > 0:

                    codon_rows.append({
                        "Codon": codon,
                        "Amino Acid": STANDARD_GENETIC_CODE[codon],
                        "Frequency": count,
                        "RSCU": round(
                            results["rscu_scores"][codon], 4
                        ),
                        "Relative Adaptiveness (wi)": (
                            None
                            if results["relative_adaptiveness"][codon] is None
                            else round(
                                results["relative_adaptiveness"][codon],
                                4
                            )
                        )
                    })

            codon_df = pd.DataFrame(codon_rows)

            st.dataframe(
                codon_df,
                use_container_width=True,
                hide_index=True
            )

            # --------------------------------------------------
            # RSCU GRAPH
            # --------------------------------------------------

            st.divider()
            st.header("4. RSCU Visualization")

            plot_data = codon_df[
                ~codon_df["Codon"].isin(["ATG", "TGG"])
            ][["Codon", "RSCU"]].set_index("Codon")

            if not plot_data.empty:
                st.bar_chart(plot_data)
                st.caption(
                    "RSCU = 1 indicates equal synonymous codon usage. "
                    "Values above or below 1 indicate relative preference "
                    "or underuse."
                )
            else:
                st.info(
                    "No observed codons with synonymous alternatives "
                    "are available for plotting."
                )

            # --------------------------------------------------
            # REFERENCE ADAPTIVENESS
            # --------------------------------------------------

            st.divider()
            st.header("5. Relative Adaptiveness (wi)")

            st.write(
                "These weights are calculated using the supplied "
                "reference codon frequencies."
            )

            weight_rows = []

            for codon, weight in estimator.weights.items():

                weight_rows.append({
                    "Codon": codon,
                    "Amino Acid": STANDARD_GENETIC_CODE[codon],
                    "Reference Count": reference_frequencies.get(codon, 0),
                    "Relative Adaptiveness (wi)": (
                        None if weight is None
                        else round(weight, 4)
                    )
                })

            weight_df = pd.DataFrame(weight_rows)

            st.dataframe(
                weight_df,
                use_container_width=True,
                hide_index=True
            )

            # --------------------------------------------------
            # CAI AND EXPRESSION ESTIMATION
            # --------------------------------------------------

            st.divider()
            st.header("6. Expression Estimation")

            if cai is not None:

                st.metric(
                    "Codon Adaptation Index (CAI)",
                    f"{cai:.4f}"
                )

                classification = evaluation["classification"]

                if classification == "High codon adaptation":
                    st.success(classification)

                elif classification == "Moderate codon adaptation":
                    st.warning(classification)

                else:
                    st.info(classification)

            else:
                st.warning(
                    "CAI could not be calculated because there were "
                    "not enough synonymous codons with reference data."
                )

            st.write(
                "**Codons included in CAI:**",
                evaluation["evaluated_codon_count"]
            )

            st.write(
                "**Excluded codons:**",
                evaluation["excluded_codons"]
            )

            st.write(
                "**Skipped due to missing reference data:**",
                evaluation["skipped_no_reference"]
            )

            # --------------------------------------------------
            # SCIENTIFIC NOTE
            # --------------------------------------------------

            st.divider()

            st.info(
                "**Scientific interpretation:** CAI estimates how closely "
                "a gene's synonymous codon usage matches a reference set "
                "of highly expressed genes. It is an indicator of codon "
                "adaptation, not a direct measurement of mRNA or protein "
                "expression."
            )

            st.success("Sequence analysis completed successfully!")

        except (ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
            st.error(f"Analysis error: {error}")

        except Exception as error:
            st.error(
                f"Unexpected error: {error}"
            )

# --------------------------------------------------
# FOOTER
# --------------------------------------------------

st.divider()

st.caption(
    "GeneExpress | Codon Usage-Based Expression Level Estimator"
)
