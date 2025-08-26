#!/usr/bin/env python3
"""
Download basis sets from Basis Set Exchange for DFT calculations
"""

import json
import sys
import requests
import argparse

def download_basis_set(basis_name, elements, output_file=None):
    """
    Download basis set from BSE API
    
    Args:
        basis_name: Name of basis set (e.g., 'def2-svp', 'sto-3g', '6-31g')
        elements: List of element symbols (e.g., ['H', 'C', 'O'])
        output_file: Output filename (default: basis_name.json)
    """
    
    # BSE API endpoint
    base_url = "https://www.basissetexchange.org/api"
    
    # Convert element symbols to atomic numbers
    element_map = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
        'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
        'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
        # Add more as needed
    }
    
    element_nums = []
    for elem in elements:
        if elem in element_map:
            element_nums.append(str(element_map[elem]))
        else:
            print(f"Warning: Unknown element {elem}, skipping")
    
    if not element_nums:
        print("Error: No valid elements specified")
        return False
    
    # Build request URL
    elements_str = ','.join(element_nums)
    url = f"{base_url}/basis/{basis_name}/format/json"
    params = {
        'elements': elements_str,
        'uncontract_general': 'false',
        'uncontract_segmented': 'false',
        'uncontract_spdf': 'false',
        'optimize_general': 'true'
    }
    
    print(f"Downloading {basis_name} basis set for elements {elements}...")
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        
        # Parse JSON response
        basis_data = response.json()
        
        # Save to file
        if output_file is None:
            element_str = ''.join(elements)
            output_file = f"{basis_name}.{element_str}.json"
        
        with open(output_file, 'w') as f:
            json.dump(basis_data, f, indent=2)
        
        print(f"Successfully saved basis set to {output_file}")
        
        # Print summary
        if 'elements' in basis_data:
            n_shells = 0
            n_primitives = 0
            for elem_num, elem_data in basis_data['elements'].items():
                if 'electron_shells' in elem_data:
                    n_shells += len(elem_data['electron_shells'])
                    for shell in elem_data['electron_shells']:
                        if 'exponents' in shell:
                            n_primitives += len(shell['exponents'])
            
            print(f"Basis set contains {n_shells} shells with {n_primitives} primitives total")
        
        return True
        
    except requests.exceptions.RequestException as e:
        print(f"Error downloading basis set: {e}")
        return False
    except json.JSONDecodeError as e:
        print(f"Error parsing response: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Download basis sets from Basis Set Exchange')
    parser.add_argument('basis', help='Basis set name (e.g., def2-svp, sto-3g, 6-31g)')
    parser.add_argument('elements', nargs='+', help='Element symbols (e.g., H C O)')
    parser.add_argument('-o', '--output', help='Output filename (default: basis.elements.json)')
    
    args = parser.parse_args()
    
    download_basis_set(args.basis, args.elements, args.output)

if __name__ == '__main__':
    main()
