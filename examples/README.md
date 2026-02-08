# Example Data

This directory contains example outputs and expected results for testing the toolkit.

## Files

- `example_features_organelles.csv` - Sample Nellie output for one cell
- `example_cell_metrics.csv` - Expected output from batch_analyse.py
- `example_comparison.png` - Sample visualisation

## Running the Examples

### Test Batch Analysis

```bash
python scripts/batch_analyse.py \
  --in examples/ \
  --out examples/test_output/
```

This should produce a `*_cell_metrics.csv` file that matches the expected output.

### Test Statistical Comparison

First, create two example replicate files by copying the example:

```bash
cp examples/example_cell_metrics.csv examples/replicate1.csv
cp examples/example_cell_metrics.csv examples/replicate2.csv
```

Then run the comparison:

```bash
python scripts/compare_replicates.py \
  --rep1 examples/replicate1.csv \
  --rep2 examples/replicate2.csv \
  --out examples/test_comparison
```

## Notes

- Example data represents typical Nellie output
- Metrics should fall within normal biological ranges
- Use these examples to verify installation
