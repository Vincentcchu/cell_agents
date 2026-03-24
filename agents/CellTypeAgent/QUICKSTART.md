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

### Batch automation across all tissues

If your formatted files are organized as:

```text
data/<tissue>/celltypeagent_format/*_formatted.csv
```

you can run all datasets without manually editing `get_prediction.py`:

```bash
cd cell_agents/agents/CellTypeAgent
python run_all_tissues.py --model gpt-5.1 --species human
```

Useful options:

```bash
# Only selected tissues
python run_all_tissues.py --tissues brain breast

# Resume from a dataset
python run_all_tissues.py --start-from brain/Data_Choudhury2022_Brain_formatted

# Skip datasets already processed for this model
python run_all_tissues.py --skip-existing
```

## Output

Results are saved to:
- **Default location**: `batch_outputs/`
- **Organized by**: `batch_outputs/<tissue>/<dataset>/...`
- **Includes**:
   - prediction CSVs and logs from CellTypeAgent
   - `run_metrics.json` (timing/tokens/cost)
   - `top_<n>_max_<markers>.txt` per dataset
   - batch-level summaries: `BATCH_SUMMARY.json` and `BATCH_SUMMARY.csv`

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
