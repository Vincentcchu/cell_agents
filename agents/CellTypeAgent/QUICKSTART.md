# CellTypeAgent

Specialized agent for automated cell type identification.

## Overview

CellTypeAgent provides automated cell type annotation using machine learning and LLM-based approaches.

## Prerequisites

- Python 3.8+
- Required packages (see `requirements.txt`)
- Input data in CellTypeAgent format

## Setup

1. **Navigate to agent directory**:
   ```bash
   cd agents/celltypeagent
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Install CellTypeAgent**:
   ```bash
   cd CellTypeAgent
   # Follow installation instructions in CellTypeAgent/README.md
   ```

## Input Data

- **Expected format**: CellTypeAgent-specific format
- **Preprocessing required**: Yes - see `../../preprocessing/agent_specific/celltypeagent/`
- **Conversion script**: `../../preprocessing/agent_specific/celltypeagent/format_convert_celltypeagent.ipynb`

### Data Preparation:

1. Convert h5ad to CellTypeAgent format:
   ```bash
   cd ../../preprocessing/agent_specific/celltypeagent/
   jupyter notebook format_convert_celltypeagent.ipynb
   ```

2. Prepared data will be ready for CellTypeAgent input

## Running

Navigate to the CellTypeAgent package directory and follow the execution instructions:

```bash
cd agents/celltypeagent/CellTypeAgent
# Run according to CellTypeAgent documentation
```

## Output

Results are saved to:
- **Default location**: `../../outputs/celltypeagent/`
- **Output format**: Cell type predictions with confidence scores
- **Includes**: Annotation results, analysis logs

## Notes

- **Requires preprocessing**: Standard h5ad files must be converted to CellTypeAgent format
- This agent may be a git submodule pointing to an external repository
- See original CellTypeAgent documentation: [CellTypeAgent/README.md](CellTypeAgent/README.md)

## Troubleshooting

**Issue**: Format incompatibility
- **Solution**: Use the preprocessing notebook to convert data

**Issue**: Module import errors
- **Solution**: Ensure CellTypeAgent package is properly installed

**Issue**: Path errors
- **Solution**: Verify paths in preprocessing scripts point to correct data locations

## Additional Resources

- Main documentation: [README.md](README.md)
- CellTypeAgent docs: [CellTypeAgent/README.md](CellTypeAgent/README.md)
- Figures: [`Figs/`](Figs/)
