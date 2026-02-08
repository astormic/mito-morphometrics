# mito-morphometrics

A Python toolkit for batch processing and statistical analysis of mitochondrial morphology data from [Nellie](https://github.com/aeleftherios/nellie) segmentation outputs.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

## Overview

This tool processes Nellie's `features_organelles.csv` files to compute comprehensive cell-level mitochondrial metrics including size distributions, shape descriptors, and network properties. It automates the analysis of multiple samples and provides publication-ready statistical comparisons between experimental groups.

## Features

- **Batch Processing**: Automatically processes all CSV files in a directory
- **Comprehensive Metrics**: Calculates 15+ morphological parameters including:
  - **Basic metrics**: mitochondrial count, total/mean/median area, length distributions
  - **Shape descriptors**: solidity, extent, aspect ratio
  - **Network properties**: fragmentation index, tubularity score, network dominance, shape heterogeneity
- **Statistical Analysis**: Automated replicate comparison with appropriate statistical tests
  - Shapiro-Wilk normality testing
  - Welch's t-test for normally distributed data
  - Mann-Whitney U test for non-parametric data
- **Visualisation**: Publication-ready bar plots with significance indicators
- **Multiple Output Formats**: PNG and SVG for figures, CSV for data

## Installation

### Requirements

- Python 3.7 or higher
- Required packages:
  ```bash
  pip install numpy pandas scipy matplotlib seaborn
  ```

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/mito-morphometrics.git
cd mito-morphometrics

# Install dependencies
pip install -r requirements.txt

# Make scripts executable (optional, Unix/Mac)
chmod +x scripts/*.py
```

## Usage

### 1. Batch Analysis

Process multiple Nellie output CSV files at once:

**Interactive GUI mode:**
```bash
python scripts/batch_analyse.py
```

**Command line mode:**
```bash
python scripts/batch_analyse.py --in /path/to/nellie/csvs --out /path/to/output
```

**Recursive search in subdirectories:**
```bash
python scripts/batch_analyse.py --in /path/to/data --out /path/to/output --recursive
```

**Output:**
- Individual `*_cell_metrics.csv` files for each input
- `combined_cell_metrics.csv` with all results

### 2. Statistical Comparison

Compare two experimental replicates or conditions:

```bash
python scripts/compare_replicates.py \
  --rep1 replicate1_combined.csv \
  --rep2 replicate2_combined.csv \
  --out results/comparison
```

**Output:**
- `comparison_statistics.csv` - detailed statistical results
- `comparison_all_metrics.png/svg` - comprehensive bar plots
- `comparison_key_metrics.png/svg` - focused visualisation of key parameters

## Input Format

Expects CSV files from Nellie with the following columns:
- `organelle_area_raw`
- `organelle_axis_length_maj_raw`
- `organelle_axis_length_min_raw`
- `organelle_extent_raw`
- `organelle_solidity_raw`

## Output Metrics

### Basic Morphology
- `n_mitochondria` - Total number of detected mitochondria
- `total_mito_area_raw` - Sum of all mitochondrial areas
- `mean_mito_area_raw` - Average mitochondrial area
- `median_mito_area_raw` - Median mitochondrial area
- `sd_mito_area_raw` - Standard deviation of area

### Length Metrics
- `mean_mito_length_raw` - Average major axis length
- `max_mito_length_raw` - Maximum length (longest mitochondrion)
- `sd_mito_length_raw` - Standard deviation of length

### Shape Descriptors
- `mean_solidity` - Compactness measure (area/convex hull area)
- `mean_extent` - Proportion of bounding box filled
- `mean_aspect_ratio` - Elongation measure (major/minor axis)

### Network Properties
- `fragmentation_index` - Mitochondrial count normalised to total area (higher = more fragmented)
- `network_dominance_index` - Maximum length relative to mean (higher = less uniform network)
- `shape_heterogeneity_index` - Coefficient of variation of length (higher = more variable shapes)
- `tubularity_score` - Composite measure: aspect_ratio × (1 - solidity) (higher = more tubular)

## Example Workflow

```bash
# 1. Process Nellie outputs for all your samples
python scripts/batch_analyse.py --in ./raw_data --out ./processed

# 2. Compare experimental groups
python scripts/compare_replicates.py \
  --rep1 ./processed/control_combined.csv \
  --rep2 ./processed/treatment_combined.csv \
  --out ./results/control_vs_treatment
```

## Citation

If you use this tool in your research, please cite:

**Nellie:**
> [Nellie citation - check their repository for the most current citation]

**This tool:**
> Rahmani. A., (2026). mito-morphometrics: Batch analysis toolkit for mitochondrial morphology. GitHub: https://github.com/astormic/mito-morphometrics

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This tool processes data generated by [Nellie](https://github.com/aeleftherios/nellie). We are grateful to the Nellie developers for their excellent open-source segmentation tool.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any issues or have questions:
1. Check the [documentation](docs/)
2. Open an [issue](https://github.com/yourusername/mito-morphometrics/issues)
3. Contact: [your email]

## Roadmap

- [ ] Add support for time-series analysis
- [ ] Implement additional network metrics
- [ ] Add GUI for easier use
- [ ] Integration with other segmentation tools
- [ ] Docker containerisation

## Version History

### v1.0.0 (2026)
- Initial release
- Batch processing of Nellie outputs
- Statistical comparison tools
- Publication-ready visualisations
