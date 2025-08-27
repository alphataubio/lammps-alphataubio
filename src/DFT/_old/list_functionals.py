#!/usr/bin/env python3
"""
List all available functionals in libxc
"""

import subprocess
import re

def list_libxc_functionals():
    """List all available functionals from libxc"""
    
    # Common functional categories
    categories = {
        'LDA': [],
        'GGA': [], 
        'Meta-GGA': [],
        'Hybrid': [],
        'Range-separated': []
    }
    
    # Try to get list from xc-info utility if available
    try:
        result = subprocess.run(['xc-info', '--list'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            lines = result.stdout.split('\n')
            for line in lines:
                if line.strip():
                    print(line)
    except FileNotFoundError:
        print("xc-info not found, showing common functionals instead:\n")
        
        # Common functionals by category
        common = {
            'LDA': [
                'LDA_X - Slater exchange',
                'LDA_C_PW - Perdew-Wang correlation',
                'LDA_C_VWN - VWN correlation',
                'LDA_C_PZ - Perdew-Zunger correlation'
            ],
            'GGA': [
                'GGA_X_PBE - PBE exchange',
                'GGA_C_PBE - PBE correlation', 
                'GGA_X_B88 - Becke 88 exchange',
                'GGA_C_LYP - Lee-Yang-Parr correlation',
                'GGA_C_P86 - Perdew 86 correlation',
                'GGA_X_PW91 - PW91 exchange',
                'GGA_C_PW91 - PW91 correlation'
            ],
            'Hybrid GGA': [
                'HYB_GGA_XC_B3LYP - B3LYP',
                'HYB_GGA_XC_PBEH - PBE0',
                'HYB_GGA_XC_HSE06 - HSE06'
            ],
            'Meta-GGA': [
                'MGGA_X_TPSS - TPSS exchange',
                'MGGA_C_TPSS - TPSS correlation',
                'MGGA_X_SCAN - SCAN exchange',
                'MGGA_C_SCAN - SCAN correlation',
                'MGGA_X_MS2 - MS2 exchange',
                'MGGA_C_MS2 - MS2 correlation'
            ],
            'Hybrid Meta-GGA': [
                'HYB_MGGA_XC_M06 - Minnesota M06',
                'HYB_MGGA_XC_M06_2X - M06-2X',
                'HYB_MGGA_XC_M06_HF - M06-HF'
            ]
        }
        
        for category, functionals in common.items():
            print(f"\n{category}:")
            print("-" * 40)
            for func in functionals:
                print(f"  {func}")
    
    print("\n" + "=" * 50)
    print("Common functional mappings in DFT package:")
    print("=" * 50)
    
    mappings = [
        ("PBE", "GGA_X_PBE+GGA_C_PBE"),
        ("B3LYP", "HYB_GGA_XC_B3LYP"),
        ("LDA", "LDA_X+LDA_C_PW"),
        ("BLYP", "GGA_X_B88+GGA_C_LYP"),
        ("BP86", "GGA_X_B88+GGA_C_P86"),
        ("PBE0", "HYB_GGA_XC_PBEH"),
        ("HSE06", "HYB_GGA_XC_HSE06"),
        ("TPSS", "MGGA_X_TPSS+MGGA_C_TPSS"),
        ("SCAN", "MGGA_X_SCAN+MGGA_C_SCAN"),
        ("M06", "HYB_MGGA_XC_M06"),
        ("M06-2X", "HYB_MGGA_XC_M06_2X"),
        ("wB97M-V", "(custom implementation)")
    ]
    
    for common, libxc in mappings:
        print(f"  {common:10s} -> {libxc}")

if __name__ == '__main__':
    list_libxc_functionals()
