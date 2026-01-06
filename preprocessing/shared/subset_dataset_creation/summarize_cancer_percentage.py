#!/usr/bin/env python3
"""
Script to summarize the cancer cell percentages from the cancer_analysis.txt file
"""

import re

def extract_percentages_from_file(file_path):
    """Extract cancer cell percentages from the analysis file"""
    results = {}
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Extract malignancy status info
        malignancy_match = re.search(r'Malignancy Percentages:\s+malignant:\s+([\d,]+)\s+cells\s+\(([\d.]+)%\)', content)
        if malignancy_match:
            count = int(malignancy_match.group(1).replace(',', ''))
            percentage = float(malignancy_match.group(2))
            results['malignancy'] = {
                'count': count,
                'percentage': percentage
            }
        
        # Extract disease status info
        disease_match = re.search(r'Disease Percentages:\s+cancer:\s+([\d,]+)\s+cells\s+\(([\d.]+)%\)', content)
        if disease_match:
            count = int(disease_match.group(1).replace(',', ''))
            percentage = float(disease_match.group(2))
            results['disease'] = {
                'count': count,
                'percentage': percentage
            }
        
        # Extract malignant cell annotation info
        annotation_match = re.search(r'Malignant cells:\s+([\d,]+)\s+\(([\d.]+)%\)', content)
        if annotation_match:
            count = int(annotation_match.group(1).replace(',', ''))
            percentage = float(annotation_match.group(2))
            results['cell_annotation'] = {
                'count': count,
                'percentage': percentage
            }
        
        # Get total cells
        total_match = re.search(r'Total cells:\s+([\d,]+)', content)
        if total_match:
            results['total_cells'] = int(total_match.group(1).replace(',', ''))
        
        return results
    
    except Exception as e:
        print(f"Error extracting percentages: {e}")
        return None

def summarize_cancer_percentages(analysis_results):
    """Create a summary of the cancer cell percentages"""
    summary = "Cancer Cell Percentage Summary\n"
    summary += "============================\n\n"
    
    if 'total_cells' in analysis_results:
        summary += f"Total cells in dataset: {analysis_results['total_cells']:,}\n\n"
    
    summary += "Cancer cell percentages by different classifications:\n\n"
    
    if 'malignancy' in analysis_results:
        summary += f"1. Based on malignancy status: {analysis_results['malignancy']['percentage']:.2f}%\n"
        summary += f"   ({analysis_results['malignancy']['count']:,} cells marked as malignant)\n\n"
    
    if 'disease' in analysis_results:
        summary += f"2. Based on disease status: {analysis_results['disease']['percentage']:.2f}%\n"
        summary += f"   ({analysis_results['disease']['count']:,} cells marked with cancer)\n\n"
    
    if 'cell_annotation' in analysis_results:
        summary += f"3. Based on explicit cell annotation: {analysis_results['cell_annotation']['percentage']:.2f}%\n"
        summary += f"   ({analysis_results['cell_annotation']['count']:,} cells labeled as malignant)\n\n"
    
    # Determine the most reliable metric
    reliable_metric = None
    explanation = ""
    
    if 'disease' in analysis_results and analysis_results['disease']['count'] > 0:
        reliable_metric = 'disease'
        explanation = "The disease status directly classifies cells as cancer or normal."
    elif 'malignancy' in analysis_results and analysis_results['malignancy']['count'] > 0:
        reliable_metric = 'malignancy'
        explanation = "The malignancy status directly identifies malignant cells."
    elif 'cell_annotation' in analysis_results and analysis_results['cell_annotation']['count'] > 0:
        reliable_metric = 'cell_annotation'
        explanation = "Cell annotation explicitly labels some cells as malignant."
    
    if reliable_metric:
        summary += "CONCLUSION:\n"
        summary += f"The most reliable measure indicates that {analysis_results[reliable_metric]['percentage']:.2f}% "
        summary += f"of cells in the dataset are cancer cells.\n"
        summary += f"({analysis_results[reliable_metric]['count']:,} out of {analysis_results.get('total_cells', 'unknown')} cells)\n\n"
        summary += f"Reasoning: {explanation}\n\n"
    
    summary += "NOTE: Different classification methods yield different results due to:\n"
    summary += "1. Varying annotation completeness across cells\n"
    summary += "2. Different criteria used for classifying cells\n"
    summary += "3. Some cells may lack explicit cancer/normal labeling\n"
    
    return summary

if __name__ == "__main__":
    analysis_file = "cancer_analysis.txt"
    output_file = "cancer_percentage_summary.txt"
    
    results = extract_percentages_from_file(analysis_file)
    if results:
        summary = summarize_cancer_percentages(results)
        
        with open(output_file, 'w') as f:
            f.write(summary)
        
        print(f"Summary saved to {output_file}")
        print("\nSUMMARY:")
        print("--------")
        print(summary)
    else:
        print(f"Could not extract cancer percentages from {analysis_file}")
