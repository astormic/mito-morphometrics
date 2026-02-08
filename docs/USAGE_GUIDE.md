# User Guide

## Quick Start Tutorial

### Step 1: Prepare Your Data

After running Nellie on your cell images, you should have `features_organelles.csv` files. Organise them like this:

```
project/
├── raw_data/
│   ├── control/
│   │   ├── cell001_features_organelles.csv
│   │   ├── cell002_features_organelles.csv
│   │   └── ...
│   └── treatment/
│       ├── cell001_features_organelles.csv
│       ├── cell002_features_organelles.csv
│       └── ...
└── analysis/
```

### Step 2: Batch Process Your Data

Process all control samples:
```bash
python scripts/batch_analyse.py \
  --in raw_data/control \
  --out analysis/control
```

Process all treatment samples:
```bash
python scripts/batch_analyse.py \
  --in raw_data/treatment \
  --out analysis/treatment
```

### Step 3: Statistical Comparison

Compare the two groups:
```bash
python scripts/compare_replicates.py \
  --rep1 analysis/control/combined_cell_metrics.csv \
  --rep2 analysis/treatment/combined_cell_metrics.csv \
  --out analysis/control_vs_treatment
```

### Step 4: Interpret Results

Open the generated files:
- `control_vs_treatment_statistics.csv` - numerical results
- `control_vs_treatment_key_metrics.svg` - main figure for publication
- `control_vs_treatment_all_metrics.svg` - comprehensive results

## Understanding the Metrics

### Fragmentation Index
**What it measures:** Number of mitochondria normalised to total area

**Interpretation:**
- **Higher values** = more fragmented network (many small mitochondria)
- **Lower values** = more fused network (fewer, larger mitochondria)

**Biological relevance:** Increases with mitochondrial fission, apoptosis, stress

### Tubularity Score
**What it measures:** Composite of elongation and compactness

**Calculation:** `mean_aspect_ratio × (1 - mean_solidity)`

**Interpretation:**
- **Higher values** = more elongated, tubular mitochondria
- **Lower values** = more rounded, punctate mitochondria

**Biological relevance:** Healthy cells typically have tubular mitochondria; rounded mitochondria may indicate dysfunction

### Network Dominance Index
**What it measures:** Maximum length relative to mean length

**Interpretation:**
- **Higher values** = network dominated by few long mitochondria
- **Lower values** = more uniform distribution of lengths

**Biological relevance:** Indicates heterogeneity in the mitochondrial network

### Shape Heterogeneity Index
**What it measures:** Coefficient of variation of mitochondrial length

**Interpretation:**
- **Higher values** = mixed population of sizes
- **Lower values** = uniform population

**Biological relevance:** May reflect dynamic remodeling or mixed populations

## Statistical Methods

### Test Selection
The tool automatically chooses the appropriate test:

1. **Shapiro-Wilk test** checks if data are normally distributed
2. If both groups are normal (p > 0.05):
   - **Welch's t-test** (doesn't assume equal variances)
3. If either group is non-normal:
   - **Mann-Whitney U test** (non-parametric alternative)

### Significance Levels
- `***` p < 0.001 (highly significant)
- `**` p < 0.01 (very significant)
- `*` p < 0.05 (significant)
- `ns` p ≥ 0.05 (not significant)

## Troubleshooting

### Common Issues

**1. "Missing required columns" error**
- Make sure you're using Nellie's `features_organelles.csv` output
- Check that column names match exactly (case-sensitive)

**2. No CSV files found**
- Check your input path
- Use `--recursive` flag if files are in subdirectories

**3. All metrics show "inf" or "nan"**
- Some cells may have unusual values (e.g., division by zero)
- Check your source images for quality
- These values are automatically excluded from statistical tests

**4. Import errors**
- Install required packages: `pip install -r requirements.txt`
- Make sure you're using Python 3.7+

### Getting Help

If you encounter issues:
1. Check this documentation
2. Look at example data in `examples/`
3. Open an issue on GitHub with:
   - Error message
   - Your Python version
   - First few rows of your CSV file

## Best Practices

### Sample Size
- Aim for **n ≥ 15 cells** per group for reliable statistics
- More samples = more statistical power
- Document any exclusion criteria

### Image Quality
- Use consistent imaging parameters across all samples
- Ensure proper background subtraction in Nellie
- Document microscope settings and acquisition parameters

### Data Organisation
- Use clear, descriptive filenames
- Keep raw data separate from processed data
- Document experimental conditions in a spreadsheet

### Reporting Results
Include in your methods:
1. Number of cells analyzed per group
2. Statistical test used (automatically chosen by tool)
3. Significance threshold (typically p < 0.05)
4. Which metrics showed significant differences

Example text:
> "Mitochondrial morphology was quantified using mito-morphometrics (v1.0). 
> Statistical comparisons between groups were performed using Welch's t-test 
> or Mann-Whitney U test as appropriate based on normality testing (n=50 
> control cells, n=16 treatment cells). Significance was set at p < 0.05."

## Advanced Usage

### Custom Metric Selection

Modify `compare_replicates.py` to focus on specific metrics:

```python
# Edit the metrics list (around line 150)
metrics = [
    'n_mitochondria',
    'fragmentation_index',
    'tubularity_score',
    # Add or remove metrics as needed
]
```

### Batch Comparison of Multiple Groups

Create a shell script to compare multiple pairs:

```bash
#!/bin/bash
# compare_all.sh

for condition in treatment1 treatment2 treatment3; do
  python scripts/compare_replicates.py \
    --rep1 analysis/control/combined_cell_metrics.csv \
    --rep2 analysis/${condition}/combined_cell_metrics.csv \
    --out results/control_vs_${condition}
done
```

### Integration with Other Tools

Export combined CSV files to:
- **GraphPad Prism** for additional statistical tests
- **R** for more complex analyses
- **Excel** for manual inspection

The CSV format is universal and works with any statistical software.
