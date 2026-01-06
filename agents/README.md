# Cell Type Annotation Agents

This directory contains implementations of five different LLM-based agents for cell type annotation. Each agent is independently developed with its own approach, dependencies, and execution methods.

## Available Agents

### 1. BioMaster
**Location**: [`biomaster/`](biomaster/)  
**Entrypoint**: [`biomaster/run.py`](biomaster/run.py)  
**Type**: Python script  
**Documentation**: [biomaster/README.md](biomaster/README.md)

Multi-agent system with RAG capabilities for cell annotation.

### 2. Biomni
**Location**: [`Biomni/`](Biomni/)  
**Entrypoint**: [`Biomni/run.ipynb`](Biomni/run.ipynb)  
**Type**: Jupyter notebook  
**Documentation**: [Biomni/QUICKSTART.md](Biomni/QUICKSTART.md)

Interactive notebook-based agent with multiple experimental variations.

### 3. CASSIA
**Location**: [`CASSIA/`](CASSIA/)  
**Entrypoint**: [`CASSIA/CASSIA_python/`](CASSIA/CASSIA_python/)  
**Type**: Python package  
**Documentation**: [CASSIA/QUICKSTART.md](CASSIA/QUICKSTART.md)

Cell-type Annotation using Semi-supervised Iterative Algorithm.

### 4. CellTypeAgent
**Location**: [`CellTypeAgent/`](CellTypeAgent/)  
**Entrypoint**: [`CellTypeAgent/CellTypeAgent/`](CellTypeAgent/CellTypeAgent/)  
**Type**: Python package  
**Documentation**: [CellTypeAgent/QUICKSTART.md](CellTypeAgent/QUICKSTART.md)

Specialized agent for automated cell type identification.

### 5. mLLMCellType
**Location**: [`mLLMCellType/`](mLLMCellType/)  
**Entrypoint**: [`mLLMCellType/run.py`](mLLMCellType/run.py)  
**Type**: Python script  
**Documentation**: [mLLMCellType/QUICKSTART.md](mLLMCellType/QUICKSTART.md)

Multi-modal LLM approach for cell type classification.

## General Usage Pattern

Each agent follows a similar workflow:

1. **Prepare environment**: Install dependencies (see agent README)
2. **Preprocess data**: Some agents need specific input formats (see `../../preprocessing/`)
3. **Configure**: Set API keys, model parameters (agent-specific)
4. **Run**: Execute the entrypoint script/notebook
5. **Collect outputs**: Results saved to `../../outputs/{agent}/`

## Agent Independence

- Each agent has its own dependencies and requirements
- No standardized interface is enforced
- Preprocessing pipelines are agent-specific
- Output formats may vary

## Quick Reference

| Agent | Main File | Format Required | Outputs Location |
|-------|-----------|----------------|------------------|
| BioMaster | `run.py` | Standard h5ad | `outputs/biomaster/` |
| Biomni | `run.ipynb` | Standard h5ad | `outputs/biomni/` |
| CASSIA | Python package | CASSIA format | `outputs/cassia/` |
| CellTypeAgent | Python package | CellTypeAgent format | `outputs/celltypeagent/` |
| mLLMCellType | `run.py` | Standard h5ad | `outputs/mllmcelltype/` |

See individual agent READMEs for detailed instructions.
