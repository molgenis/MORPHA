#!/usr/bin/env python3
"""
mmCIF to FoldX5 Analysis Pipeline

This script performs the following operations:
1. Converts mmCIF file to PDB format
2. Repairs the PDB using FoldX5's repairpdb command
3. Analyzes the complex using FoldX5's analysecomplex command
4. Extracts and saves energy results to a file

Requirements:
- Biopython (for mmCIF to PDB conversion)
- FoldX5 executable in PATH or specified

Usage:
    python mmcif_to_foldx.py input.mmcif [--foldx-path /path/to/foldx] [--output-dir ./output]
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

try:
    from Bio.PDB import MMCIFParser, PDBIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False


def convert_mmcif_to_pdb(mmcif_path: str, pdb_path: str) -> bool:
    """
    Convert mmCIF file to PDB format using Biopython.
    
    Args:
        mmcif_path: Path to input mmCIF file
        pdb_path: Path to output PDB file
        
    Returns:
        True if conversion succeeded, False otherwise
    """
    if not HAS_BIOPYTHON:
        print("Error: Biopython is required for mmCIF to PDB conversion.")
        print("Install it with: pip install biopython")
        return False
    
    try:
        # Parse mmCIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("protein", mmcif_path)
        
        # Save as PDB
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        
        print(f"Successfully converted {mmcif_path} to {pdb_path}")
        return True
        
    except Exception as e:
        print(f"Error converting mmCIF to PDB: {e}")
        return False


def run_foldx_command(foldx_path: str, command: str, pdb_path: str, 
                      output_dir: str, work_dir: str) -> str:
    """
    Run a FoldX command and return the output file path.
    
    Args:
        foldx_path: Path to FoldX5 executable
        command: FoldX command (repairpdb, analysecomplex, etc.)
        pdb_path: Path to input PDB file
        output_dir: Directory to save output files (absolute path)
        work_dir: Working directory for FoldX (absolute path)
        
    Returns:
        Path to the output file (absolute path)
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)
    
    # Copy PDB to work directory
    pdb_name = os.path.basename(pdb_path)
    work_pdb = os.path.join(work_dir, pdb_name)
    
    try:
        import shutil
        shutil.copy2(pdb_path, work_pdb)
        print(f"Copied PDB to: {work_pdb}")
    except Exception as e:
        print(f"Error copying PDB to work directory: {e}")
        return ""
    
    # FoldX 5 syntax: use absolute paths for output-dir
    pdb_basename = os.path.splitext(pdb_name)[0]
    
    # Build FoldX command
    # For AnalyseComplex, we need to specify an output-file tag
    # FoldX will generate: Indiv_energies_TAG_AC.fxout, Interaction_TAG_AC.fxout, etc.
    if command == "repairpdb":
        cmd = [
            foldx_path,
            "--command=RepairPDB",
            f"--pdb={pdb_name}",
            f"--output-dir={output_dir}"
        ]
        output_file = os.path.join(output_dir, f"{pdb_basename}_Repair.pdb")
        
    elif command == "analysecomplex":
        # Don't specify --output-file; FoldX will use the PDB filename as the tag
        cmd = [
            foldx_path,
            "--command=AnalyseComplex",
            f"--pdb={pdb_name}",
            f"--output-dir={output_dir}"
        ]
        # FoldX generates files with tag based on PDB filename
        # We'll find them automatically in parse_energy_results
        output_file = os.path.join(output_dir, "placeholder.fxout")
        
    else:
        print(f"Unknown FoldX command: {command}")
        return ""
    
    # Run command from the work directory where the PDB file is located
    print(f"Running: {' '.join(cmd)}")
    print(f"Working directory: {work_dir}")
    print(f"Files in work directory: {os.listdir(work_dir)}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            check=True
        )
        print(f"FoldX {command} completed successfully")
        if result.stdout:
            print(f"Output: {result.stdout}")
        if result.stderr:
            print(f"Warnings: {result.stderr}")
        
        # Verify output file exists
        if not os.path.exists(output_file):
            print(f"Warning: Expected output file not found at {output_file}")
            # Check if it was created in work_dir instead
            work_output = os.path.join(work_dir, os.path.basename(output_file))
            if os.path.exists(work_output):
                print(f"Found output at: {work_output}")
                output_file = work_output
        
        return output_file
        
    except subprocess.CalledProcessError as e:
        print(f"Error running FoldX {command}: {e}")
        print(f"Return code: {e.returncode}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return ""


def find_foldx_output_files(output_dir: str) -> dict:
    """
    Find all FoldX output files in a directory.
    
    FoldX AnalyseComplex generates multiple .fxout files with patterns:
    - Summary_*.fxout
    - Indiv_energies_*_AC.fxout
    - Interaction_*_AC.fxout
    - Interface_Residues_*_AC.fxout
    
    Args:
        output_dir: Directory to search for .fxout files
        
    Returns:
        Dictionary with file types as keys and paths as values
    """
    files = {
        'summary': None,
        'indiv_energies': None,
        'interaction': None,
        'interface_residues': None
    }
    
    if not os.path.exists(output_dir):
        return files
    
    # Find all .fxout files in the directory
    for filename in os.listdir(output_dir):
        if filename.endswith('.fxout'):
            if filename.startswith('Summary_'):
                files['summary'] = os.path.join(output_dir, filename)
            elif filename.startswith('Indiv_energies_') and '_AC.fxout' in filename:
                files['indiv_energies'] = os.path.join(output_dir, filename)
            elif filename.startswith('Interaction_') and '_AC.fxout' in filename:
                files['interaction'] = os.path.join(output_dir, filename)
            elif filename.startswith('Interface_Residues_') and '_AC.fxout' in filename:
                files['interface_residues'] = os.path.join(output_dir, filename)
    
    return files


def parse_energy_results(output_dir: str, output_tag: str = None) -> dict:
    """
    Parse FoldX analysis output files to extract energy terms.
    
    Automatically finds all .fxout files in the output directory.
    
    Args:
        output_dir: Directory containing FoldX output files
        output_tag: Optional tag used for output files (for debugging)
        
    Returns:
        Dictionary containing energy terms and values
    """
    energies = {}
    
    # Find all FoldX output files
    fxout_files = find_foldx_output_files(output_dir)
    
    print(f"Looking for FoldX output files in: {output_dir}")
    for file_type, filepath in fxout_files.items():
        if filepath:
            print(f"  Found {file_type}: {filepath}")
        else:
            print(f"  Warning: {file_type} file not found")
    
    # Parse Summary file
    if fxout_files['summary']:
        print(f"\nParsing summary file: {fxout_files['summary']}")
        try:
            with open(fxout_files['summary'], 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('Number'):
                        continue
                    # Summary file has format: Description, Value
                    if ',' in line:
                        parts = [p.strip() for p in line.split(',')]
                        if len(parts) >= 2:
                            try:
                                energies[f"Summary_{parts[0]}"] = float(parts[1])
                            except ValueError:
                                pass
        except Exception as e:
            print(f"Error parsing summary file: {e}")
    
    # Parse Individual Energies file
    if fxout_files['indiv_energies']:
        print(f"\nParsing individual energies file: {fxout_files['indiv_energies']}")
        try:
            with open(fxout_files['indiv_energies'], 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('Pdb'):
                        continue
                    # Format: Chain, Residue, EnergyTerm, Value
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            chain = parts[0]
                            residue = parts[1]
                            term = parts[2]
                            value = float(parts[3])
                            key = f"Chain_{chain}_{residue}_{term}"
                            energies[key] = value
                        except ValueError:
                            pass
        except Exception as e:
            print(f"Error parsing individual energies file: {e}")
    
    # Parse Interaction file
    if fxout_files['interaction']:
        print(f"\nParsing interaction file: {fxout_files['interaction']}")
        try:
            with open(fxout_files['interaction'], 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('Pdb'):
                        continue
                    # Format: Chain1, Chain2, EnergyTerm, Value
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            chain1 = parts[0]
                            chain2 = parts[1]
                            term = parts[2]
                            value = float(parts[3])
                            key = f"Interaction_{chain1}-{chain2}_{term}"
                            energies[key] = value
                        except ValueError:
                            pass
        except Exception as e:
            print(f"Error parsing interaction file: {e}")
    
    print(f"\nExtracted {len(energies)} energy terms from FoldX output files")
    return energies


def save_energy_results(energies: dict, output_file: str, 
                        original_pdb: str, repaired_pdb: str) -> None:
    """
    Save energy results to a structured CSV file.
    
    Args:
        energies: Dictionary of energy terms and values
        output_file: Path to save the results
        original_pdb: Path to original PDB (for reference)
        repaired_pdb: Path to repaired PDB (for reference)
    """
    try:
        with open(output_file, 'w') as f:
            f.write("# FoldX5 Energy Analysis Results\n")
            f.write(f"# Original PDB: {original_pdb}\n")
            f.write(f"# Repaired PDB: {repaired_pdb}\n")
            f.write(f"# Generated: {__import__('datetime').datetime.now().isoformat()}\n")
            f.write("\n")
            f.write("Energy_Term,Value\n")
            
            for term, value in sorted(energies.items()):
                # Escape commas in term names
                term_escaped = f'"{term}"' if ',' in term else term
                f.write(f"{term_escaped},{value}\n")
        
        print(f"Energy results saved to {output_file}")
        
    except Exception as e:
        print(f"Error saving energy results: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert mmCIF to PDB, repair with FoldX5, and analyze energies"
    )
    parser.add_argument(
        'mmcif_file',
        type=str,
        help='Path to input mmCIF file'
    )
    parser.add_argument(
        '--foldx-path',
        type=str,
        default='foldx',
        help='Path to FoldX5 executable (default: foldx in PATH)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./foldx_output',
        help='Directory to save output files (default: ./foldx_output)'
    )
    parser.add_argument(
        '--work-dir',
        type=str,
        default='./foldx_work',
        help='Working directory for FoldX (default: ./foldx_work)'
    )
    parser.add_argument(
        '--keep-intermediate',
        action='store_true',
        help='Keep intermediate PDB files'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.mmcif_file):
        print(f"Error: mmCIF file not found: {args.mmcif_file}")
        sys.exit(1)
    
    if not os.path.exists(args.foldx_path):
        print(f"Error: FoldX executable not found: {args.foldx_path}")
        print("Please specify the correct path with --foldx-path")
        sys.exit(1)
    
    # Create output and work directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)
    
    # Step 1: Convert mmCIF to PDB
    pdb_name = os.path.splitext(os.path.basename(args.mmcif_file))[0] + ".pdb"
    pdb_path = os.path.join(args.output_dir, pdb_name)
    
    print("=" * 60)
    print("Step 1: Converting mmCIF to PDB")
    print("=" * 60)
    
    if not convert_mmcif_to_pdb(args.mmcif_file, pdb_path):
        print("Failed to convert mmCIF to PDB. Exiting.")
        sys.exit(1)
    
    # Step 2: Repair PDB with FoldX5
    print("\n" + "=" * 60)
    print("Step 2: Repairing PDB with FoldX5")
    print("=" * 60)
    
    # Use absolute paths
    abs_pdb_path = os.path.abspath(pdb_path)
    abs_output_dir = os.path.abspath(args.output_dir)
    abs_work_dir = os.path.abspath(args.work_dir)
    
    repaired_pdb = run_foldx_command(
        args.foldx_path,
        "repairpdb",
        abs_pdb_path,
        abs_output_dir,
        abs_work_dir
    )
    
    if not repaired_pdb:
        print("Failed to repair PDB. Exiting.")
        sys.exit(1)
    
    # Use the repaired PDB for analysis - construct absolute path
    pdb_basename = os.path.splitext(pdb_name)[0]
    repaired_pdb_path = os.path.join(abs_output_dir, f"{pdb_basename}_Repair.pdb")
    
    # Verify the repaired PDB exists before proceeding
    if not os.path.exists(repaired_pdb_path):
        print(f"Error: Repaired PDB not found at {repaired_pdb_path}")
        print(f"Directory contents:")
        for f in os.listdir(abs_output_dir):
            print(f"  {f}")
        sys.exit(1)
    
    print(f"Using repaired PDB: {repaired_pdb_path}")
    
    # Step 3: Analyze complex with FoldX5
    print("\n" + "=" * 60)
    print("Step 3: Analyzing complex with FoldX5")
    print("=" * 60)
    
    analysis_file = run_foldx_command(
        args.foldx_path,
        "analysecomplex",
        repaired_pdb_path,
        abs_output_dir,
        abs_work_dir
    )
    
    if not analysis_file:
        print("Failed to analyze complex. Exiting.")
        sys.exit(1)


if __name__ == "__main__":
    main()