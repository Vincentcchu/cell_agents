#!/bin/bash
# Quick verification script for repository reorganization
# Run from the root of cell_agents repository

echo "=== Cell Agents Repository Verification ==="
echo ""

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

check_exists() {
    if [ -e "$1" ]; then
        echo -e "${GREEN}✓${NC} $1"
        return 0
    else
        echo -e "${RED}✗${NC} $1 (missing)"
        return 1
    fi
}

check_not_exists() {
    if [ ! -e "$1" ]; then
        echo -e "${GREEN}✓${NC} $1 (correctly removed)"
        return 0
    else
        echo -e "${YELLOW}⚠${NC} $1 (still exists - should be removed)"
        return 1
    fi
}

echo "1. Checking Documentation..."
check_exists "README.md"
check_exists "agents/README.md"
check_exists "agents/biomaster/QUICKSTART.md"
check_exists "agents/Biomni/QUICKSTART.md"
check_exists "agents/CASSIA/QUICKSTART.md"
check_exists "agents/CellTypeAgent/QUICKSTART.md"
check_exists "agents/mLLMCellType/QUICKSTART.md"
check_exists "data/README.md"
check_exists "preprocessing/README.md"
check_exists "analysis/README.md"
check_exists "REORGANIZATION_SUMMARY.md"
echo ""

echo "2. Checking Directory Structure..."
check_exists "agents"
check_exists "data"
check_exists "preprocessing"
check_exists "outputs"
check_exists "analysis"
check_exists "runner"
check_exists "logs"
echo ""

echo "3. Checking Agent Directories..."
check_exists "agents/biomaster"
check_exists "agents/Biomni"
check_exists "agents/CASSIA"
check_exists "agents/CellTypeAgent"
check_exists "agents/mLLMCellType"
echo ""

echo "4. Checking Preprocessing Organization..."
check_exists "preprocessing/agent_specific"
check_exists "preprocessing/agent_specific/cassia"
check_exists "preprocessing/agent_specific/celltypeagent"
check_exists "preprocessing/shared"
check_exists "preprocessing/shared/remove_label"
check_exists "preprocessing/shared/subset_dataset_creation"
check_exists "preprocessing/shared/gene_data_to_h5ad"
echo ""

echo "5. Checking Analysis Files..."
check_exists "analysis/evaluator2.0.ipynb"
check_exists "analysis/classification_for_evaluation.csv"
check_exists "analysis/token_counter.txt"
check_exists "analysis/umap_embeddings"
echo ""

echo "6. Checking for Removed Duplicates..."
check_not_exists "biomaster"
check_not_exists "Biomni"
check_not_exists "CASSIA"
check_not_exists "CellTypeAgent"
check_not_exists "mLLMCellType"
check_not_exists "cassia_preprocessing"
check_not_exists "celltypeagnet_processing"
check_not_exists "remove_label"
check_not_exists "subset_dataset_creation"
check_not_exists "gene_data_to_h5ad"
echo ""

echo "7. Checking Agent Entrypoints..."
check_exists "agents/biomaster/run.py"
check_exists "agents/Biomni/run.ipynb"
check_exists "agents/CASSIA/CASSIA_python"
check_exists "agents/CellTypeAgent/CellTypeAgent"
check_exists "agents/mLLMCellType/run.py"
echo ""

echo "8. Checking Biomni Experiments Archive..."
check_exists "agents/Biomni/experiments"
echo ""

echo "=== Verification Complete ==="
echo ""
echo "Next steps:"
echo "1. Test each agent manually (see REORGANIZATION_SUMMARY.md)"
echo "2. Run: cd agents/{agent}/ && {run command}"
echo "3. Check outputs in: outputs/{agent}/"
echo ""
