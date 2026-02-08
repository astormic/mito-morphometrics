# Contributing to mito-morphometrics

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:
- Clear description of the problem
- Steps to reproduce
- Expected vs actual behavior
- Your Python version and OS
- Sample data (if possible)

### Suggesting Enhancements

Open an issue describing:
- The enhancement you'd like to see
- Why it would be useful
- How it might work

### Code Contributions

1. **Fork the repository**
2. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. **Make your changes**
   - Follow PEP 8 style guidelines
   - Add docstrings to functions
   - Comment complex logic
4. **Test your changes**
   - Run on example data
   - Verify outputs are correct
5. **Commit with clear messages**
   ```bash
   git commit -m "Add feature: brief description"
   ```
6. **Push and create a Pull Request**

## Code Style

- Follow PEP 8
- Use type hints where appropriate
- Write descriptive variable names
- Add docstrings to all functions
- Keep functions focused and small

### Example:
```python
def calculate_metric(data: np.ndarray) -> float:
    """
    Calculate a specific morphological metric.
    
    Args:
        data: Array of measurements
    
    Returns:
        Calculated metric value
    """
    return np.mean(data)
```

## Testing

Before submitting:
- Test on example data
- Verify CSV outputs are correct
- Check that plots generate properly
- Test with different input scenarios

## Documentation

If adding features:
- Update README.md
- Update USAGE_GUIDE.md
- Add examples if relevant
- Update docstrings

## Questions?

Feel free to open an issue for discussion before starting major work.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
