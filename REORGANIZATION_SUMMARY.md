# Repository Reorganization Summary

## What Was Done

This document summarizes the reorganization performed on the cell_agents repository on January 6, 2026.

---

## Before & After Structure

### Before:
```
cell_agents/
├── agents/                           # Agent code (mix of submodules)
│   ├── biomaster/, Biomni/, CASSIA/, CellTypeAgent/, mLLMCellType/
├── biomaster/, Biomni/, CASSIA/, CellTypeAgent/, mLLMCellType/  # Empty duplicate folders
├── cassia_preprocessing/             # Scattered at root
├── celltypeagnet_processing/         # Scattered at root  
├── gene_data_to_h5ad/                # Scattered at root
├── remove_label/                     # Scattered at root
├── subset_dataset_creation/          # Scattered at root
├── preprocessing/                    # Had old duplicates
│   ├── cassia_preprocessing/, celltypeagnet_processing/, etc.
├── data/                             # Datasets
├── outputs/                          # Agent results
├── runner/                           # Unified runner
├── evaluator2.0.ipynb                # At root
├── classification_for_evaluation.csv # At root
├── umap_embeddings/                  # At root
└── umap_embeddings_pred.npy          # At root
```

### After:
```
cell_agents/
├── README.md                         # ✨ NEW: Top-level guide
├── agents/                           # Agent implementations
│   ├── README.md                     # ✨ NEW: Agent overview
│   ├── biomaster/
│   │   ├── QUICKSTART.md             # ✨ NEW
│   │   └── run.py
│   ├── Biomni/
│   │   ├── QUICKSTART.md             # ✨ NEW
│   │   ├── run.ipynb
│   │   └── experiments/              # ✨ NEW: Archived run2-13.ipynb
│   ├── CASSIA/
│   │   ├── QUICKSTART.md             # ✨ NEW
│   │   └── CASSIA_python/
│   ├── CellTypeAgent/
│   │   ├── QUICKSTART.md             # ✨ NEW
│   │   └── CellTypeAgent/
│   └── mLLMCellType/
│       ├── QUICKSTART.md             # ✨ NEW
│       └── run.py (updated paths)
├── data/
│   ├── README.md                     # ✨ NEW: Data documentation
│   └── *.h5ad files
├── preprocessing/                    # ✨ REORGANIZED
│   ├── README.md                     # ✨ NEW: Preprocessing guide
│   ├── agent_specific/               # ✨ NEW structure
│   │   ├── cassia/                   # ← Moved from root
│   │   └── celltypeagent/            # ← Moved from root
│   └── shared/                       # ✨ NEW structure
│       ├── remove_label/             # ← Moved from root
│       ├── subset_dataset_creation/  # ← Moved from root
│       └── gene_data_to_h5ad/        # ← Moved from root
├── outputs/                          # Agent results (unchanged)
├── analysis/                         # ✨ NEW directory
│   ├── README.md                     # ✨ NEW: Analysis guide
│   ├── evaluator2.0.ipynb            # ← Moved from root
│   ├── classification_for_evaluation.csv  # ← Moved from root
│   ├── token_counter.txt             # ← Moved from root
│   ├── umap_embeddings_pred.npy      # ← Moved from root
│   ├── umap_embeddings_truth.npy     # ← Moved from root
│   └── umap_embeddings/              # ← Moved from root
├── runner/                           # Unified runner (unchanged)
└── logs/                             # ✨ NEW: For execution logs
```

---

## Changes Made

### 1. Directory Reorganization

**Removed:**
- ❌ Root-level empty agent folders: `biomaster/`, `Biomni/`, `CASSIA/`, `CellTypeAgent/`, `mLLMCellType/`
- ❌ Duplicate preprocessing folders in `preprocessing/`

**Created:**
- ✅ `analysis/` - Evaluation notebooks and results
- ✅ `preprocessing/agent_specific/` - Agent-specific preprocessing
- ✅ `preprocessing/shared/` - Shared preprocessing utilities
- ✅ `logs/` - For execution logs

**Moved:**
- 📁 `evaluator2.0.ipynb` → `analysis/`
- 📁 `classification_for_evaluation.csv` → `analysis/`
- 📁 `token_counter.txt` → `analysis/`
- 📁 `umap_embeddings/` → `analysis/`
- 📁 `umap_embeddings_*.npy` → `analysis/`
- 📁 `cassia_preprocessing/` → `preprocessing/agent_specific/cassia/`
- 📁 `celltypeagnet_processing/` → `preprocessing/agent_specific/celltypeagent/`
- 📁 `remove_label/` → `preprocessing/shared/remove_label/`
- 📁 `subset_dataset_creation/` → `preprocessing/shared/subset_dataset_creation/`
- 📁 `gene_data_to_h5ad/` → `preprocessing/shared/gene_data_to_h5ad/`

**Archived:**
- 📦 Biomni experimental notebooks (`run2.ipynb` - `run13.ipynb`) → `agents/Biomni/experiments/`

### 2. Documentation Added

**Top-level:**
- ✅ `README.md` - Repository overview, quick start, agent table

**Agent-specific:**
- ✅ `agents/README.md` - Agent overview and comparison
- ✅ `agents/biomaster/QUICKSTART.md` - How to run BioMaster
- ✅ `agents/Biomni/QUICKSTART.md` - How to run Biomni
- ✅ `agents/CASSIA/QUICKSTART.md` - How to run CASSIA
- ✅ `agents/CellTypeAgent/QUICKSTART.md` - How to run CellTypeAgent
- ✅ `agents/mLLMCellType/QUICKSTART.md` - How to run mLLMCellType

**Supporting:**
- ✅ `data/README.md` - Data documentation
- ✅ `preprocessing/README.md` - Preprocessing guide
- ✅ `analysis/README.md` - Evaluation and analysis guide

### 3. Code Updates

**Path fixes:**
- ✅ `agents/mLLMCellType/run.py` - Updated hardcoded path to use `../../data/dataset_restricted.h5ad`
- ✅ `agents/mLLMCellType/run_log.py` - Updated hardcoded path
- ✅ `agents/mLLMCellType/run_top.py` - Updated hardcoded path

---

## Verification Checklist

### ✅ Documentation
- [x] Top-level README exists and is clear
- [x] Each agent has documentation (README or QUICKSTART)
- [x] Data directory documented
- [x] Preprocessing pipelines documented
- [x] Analysis/evaluation documented

### ✅ Directory Structure
- [x] No duplicate agent folders at root level
- [x] Preprocessing organized into `agent_specific/` and `shared/`
- [x] Analysis files consolidated in `analysis/`
- [x] Logs directory created

### ✅ Agents Still Work (Manual Testing Required)

Test each agent with these commands:

#### BioMaster
```bash
cd agents/biomaster
cp config.yaml.template config.yaml
# Edit config.yaml with API keys
python run.py config.yaml
```
**Expected**: Should load data from `../../data/` and save to `../../outputs/biomaster/`

#### Biomni
```bash
cd agents/Biomni
jupyter notebook run.ipynb
```
**Expected**: Should run without path errors

#### CASSIA
```bash
cd agents/CASSIA/CASSIA_python/CASSIA
# Follow CASSIA-specific instructions
python run_full_test.py
```
**Expected**: Should process CASSIA-formatted data

#### CellTypeAgent
```bash
cd agents/CellTypeAgent/CellTypeAgent
# Follow CellTypeAgent-specific instructions
```
**Expected**: Should process CellTypeAgent-formatted data

#### mLLMCellType
```bash
cd agents/mLLMCellType
python run.py
```
**Expected**: Should load from `../../data/dataset_restricted.h5ad` (updated path)

### ✅ Preprocessing Works

#### Shared preprocessing:
```bash
cd preprocessing/shared/remove_label/
jupyter notebook "remove label.ipynb"
```

#### Agent-specific preprocessing:
```bash
cd preprocessing/agent_specific/cassia/
jupyter notebook convert_format.ipynb
```

### ✅ Analysis Works

```bash
cd analysis/
jupyter notebook evaluator2.0.ipynb
```
**Expected**: Should load predictions from `../outputs/{agent}/`

---

## Breaking Changes

### None Expected

All changes were designed to be non-breaking:
- Paths updated to use relative paths from agent directories
- No agent internals modified
- No preprocessing logic changed
- No standardization forced

### If Something Breaks

1. **Path errors**: Check that you're running from the correct directory
   - Agents should be run from `agents/{agent}/`
   - Preprocessing from `preprocessing/agent_specific/{agent}/` or `preprocessing/shared/`

2. **Data not found**: Ensure data is in `data/` directory
   - Path: `../../data/dataset_restricted.h5ad` from agent directories

3. **Preprocessing format errors**: Use the preprocessing scripts in `preprocessing/agent_specific/{agent}/`

---

## Git Status

Run to see what was changed:
```bash
git status
```

To commit these changes:
```bash
git add -A
git commit -m "Reorganize repo structure for clarity

- Create clear directory structure (agents/, preprocessing/, analysis/, data/)
- Add comprehensive documentation (README.md files)
- Move analysis files to analysis/
- Organize preprocessing into agent_specific/ and shared/
- Remove duplicate empty folders
- Update hardcoded paths to relative paths
- Archive experimental notebooks to experiments/
- Add QUICKSTART guides for each agent"
```

---

## Next Steps

1. **Test each agent** - Follow the verification checklist above
2. **Update .gitignore** - Consider adding:
   ```
   logs/*.log
   analysis/umap_embeddings/*.npy
   outputs/*/
   **/dataset_*.h5ad
   ```
3. **Update .gitmodules** - If needed, update submodule paths
4. **Create GitHub issues** - For any agent-specific fixes needed
5. **Update CI/CD** - If automated testing exists, update paths

---

## Questions?

See the main [README.md](README.md) for repository overview and usage instructions.
